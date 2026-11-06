#!/usr/bin/env python3
"""Batched, exact BigWig means for BED3/6; compute once, plot from cache.

The bin semantics target deepTools 3.5.1's mean mode (not its entire CLI).
Dependencies: numpy, pyBigWig; matplotlib only for plotting.
"""
import argparse
import csv
import copy
import gzip
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

VERSION = "3.0.1"

SIGNAL_BLUE = ["#F8FBFD", "#D4E6F1", "#96C3DB", "#5B9EC2", "#2C78A3", "#145474", "#10384D"]
GROUP_COLORS = ["#3B6FB6", "#D98058", "#4A927D", "#9A79AD", "#C6A04A", "#608B9D"]


def read_beds(paths, labels=None):
    labels = labels or [Path(str(p)[:-3] if str(p).endswith(".gz") else p).stem for p in paths]
    if len(labels) != len(paths) or len(set(labels)) != len(labels):
        raise ValueError("BED group labels must be unique and match the BED files")
    rows = []
    for group, path in zip(labels, paths):
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for lineno, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if not 3 <= len(fields) <= 6:
                    raise ValueError(f"{path}:{lineno}: only tab-separated BED3–6 supported")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                strand = fields[5] if len(fields) == 6 else "."
                if start < 0 or end <= start or strand not in ("+", "-", "."):
                    raise ValueError(f"{path}:{lineno}: invalid coordinates or strand")
                name = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"
                rows.append((chrom, start, end, name, strand, group))
    if not rows:
        raise ValueError("No BED regions")
    return rows, labels


def validate_config(cfg):
    bs = cfg["bin_size"]
    lengths = [cfg["upstream"], cfg["downstream"]]
    if cfg["mode"] == "scale-regions":
        lengths.append(cfg["body_length"])
        if cfg["body_length"] <= 0:
            raise ValueError("body_length must be positive")
    if bs <= 0 or any(x < 0 or x % bs for x in lengths) or sum(lengths) == 0:
        raise ValueError("Window lengths must be nonnegative multiples of bin_size")
    return sum(lengths) // bs


def bin_bounds(rows, cfg):
    """Genomic low-to-high bins; reverse negative-strand rows after extraction.

    Mirror linspace's float step and cumulative zone offsets, including its
    rounding behavior; extending repeated edges to 1 bp matches deepTools.
    """
    starts = np.array([r[1] for r in rows], dtype=np.int64)
    ends = np.array([r[2] for r in rows], dtype=np.int64)
    negative = np.array([r[4] == "-" for r in rows])
    before = np.where(negative, cfg["downstream"], cfg["upstream"])
    after = np.where(negative, cfg["upstream"], cfg["downstream"])
    nbins = validate_config(cfg)
    lo = np.empty((len(rows), nbins), dtype=np.int64)
    hi = np.empty_like(lo)
    short = np.zeros(len(rows), dtype=bool)
    if cfg["mode"] == "reference-point":
        point = cfg["reference_point"]
        if point == "TSS":
            anchor = np.where(negative, ends, starts)
        elif point == "TES":
            anchor = np.where(negative, starts, ends)
        else:
            anchor = starts + (ends - starts) // 2
        lo[:] = anchor[:, None] - before[:, None] + np.arange(nbins) * cfg["bin_size"]
        hi[:] = lo + cfg["bin_size"]
    else:
        short = ends - starts < cfg["bin_size"]
        for reverse in (False, True):
            idx = np.flatnonzero(negative == reverse)
            left = cfg["downstream"] if reverse else cfg["upstream"]
            right = cfg["upstream"] if reverse else cfg["downstream"]
            counts = [left // cfg["bin_size"], cfg["body_length"] // cfg["bin_size"],
                      right // cfg["bin_size"]]
            origin = starts[idx] - left
            lengths = [np.full(len(idx), left), ends[idx] - starts[idx], np.full(len(idx), right)]
            offset = np.zeros(len(idx), dtype=np.int64)
            col = 0
            for length, count in zip(lengths, counts):
                if count:
                    edges = (offset[:, None] + (length / float(count))[:, None] * np.arange(count)).astype(np.int64)
                    last = (offset + length)[:, None]
                    stop = np.maximum(np.concatenate((edges[:, 1:], last), axis=1), edges + 1)
                    lo[idx, col:col + count] = origin[:, None] + edges
                    hi[idx, col:col + count] = origin[:, None] + stop
                    col += count
                offset += length
    return lo, hi, negative, short


def batches(lo, hi, chunk_bp, batch_rows, merge_gap):
    """One BigWig read per nearby batch, with bounds on span and row count."""
    left, right = lo.min(axis=1), hi.max(axis=1)
    order = np.argsort(left, kind="stable")
    pos = 0
    while pos < len(order):
        stop = pos + 1
        begin = left[order[pos]]
        end = right[order[pos]]
        while stop < min(len(order), pos + batch_rows):
            candidate = max(end, right[order[stop]])
            if candidate - begin > chunk_bp or left[order[stop]] > end + merge_gap:
                break
            end = candidate
            stop += 1
        yield order[pos:stop]
        pos = stop


def fill_sample(matrix, sample, bw_path, rows, cfg, chunk_bp, batch_rows, merge_gap=10000, row_ids=None):
    import pyBigWig

    by_chrom = {}
    for i, row in enumerate(rows):
        by_chrom.setdefault(row[0], []).append(i)
    with pyBigWig.open(str(bw_path)) as bw:
        sizes = bw.chroms()
        missing = sorted(set(by_chrom) - set(sizes))
        if missing:
            raise ValueError(f"{bw_path}: missing chromosomes {missing}; use matching genome/names")
        for chrom, indexes in by_chrom.items():
            ids = np.asarray(indexes)
            # Bound coordinate arrays as well as BigWig buffers for large BEDs.
            ids = ids[np.argsort([rows[i][1] for i in ids], kind="stable")]
            for block in range(0, len(ids), batch_rows):
                block_ids = ids[block:block + batch_rows]
                region_rows = [rows[i] for i in block_ids]
                if any(r[2] > sizes[chrom] for r in region_rows):
                    raise ValueError(f"BED region exceeds {chrom} length in {bw_path}")
                lo, hi, reverse, short = bin_bounds(region_rows, cfg)
                for part in batches(lo, hi, chunk_bp, batch_rows, merge_gap):
                    lower, upper = lo[part], hi[part]
                    start = max(0, min(sizes[chrom], int(lower.min())))
                    end = max(0, min(sizes[chrom], int(upper.max())))
                    if end > start:
                        values = bw.values(chrom, start, end, numpy=True).astype(np.float64)
                    else:
                        values = np.empty(0, dtype=np.float64)
                    valid = np.isfinite(values)
                    values[~valid] = 0
                    sums = np.concatenate(([0.0], np.cumsum(values)))
                    lidx = np.clip(lower - start, 0, len(values))
                    ridx = np.clip(upper - start, 0, len(values))
                    total = sums[ridx] - sums[lidx]
                    if cfg["missing_as_zero"]:
                        denom = upper - lower
                    else:
                        counts = np.concatenate(([0], np.cumsum(valid, dtype=np.int64)))
                        denom = counts[ridx] - counts[lidx]
                    result = np.full(total.shape, np.nan)
                    np.divide(total, denom, out=result, where=denom > 0)
                    result[short[part], :] = 0 if cfg["missing_as_zero"] else np.nan
                    result[reverse[part], :] = result[reverse[part], ::-1]
                    target_ids = block_ids[part] if row_ids is None else row_ids[block_ids[part]]
                    matrix[target_ids, sample, :] = result


def compute_task(task):
    """Workers write disjoint rows/sample slices of the shared disk matrix."""
    matrix_path, sample, bw, ids, rows, cfg, chunk_bp, batch_rows, merge_gap = task
    matrix = np.load(matrix_path, mmap_mode="r+", allow_pickle=False)
    fill_sample(matrix, sample, bw, rows, cfg, chunk_bp, batch_rows, merge_gap, np.asarray(ids))
    matrix.flush()
    return len(rows)


def parallel_tasks(matrix_path, bws, rows, cfg, chunk_bp, batch_rows, merge_gap):
    by_chrom = {}
    for i, row in enumerate(rows):
        by_chrom.setdefault(row[0], []).append(i)
    for sample, bw in enumerate(bws):
        for ids in by_chrom.values():
            ids.sort(key=lambda i: rows[i][1])
            # Balance very large chromosomes without changing row identities.
            for pos in range(0, len(ids), 20000):
                part = ids[pos:pos + 20000]
                yield (str(matrix_path), sample, bw, part, [rows[i] for i in part],
                       cfg, chunk_bp, batch_rows, merge_gap)


def row_filter(block, skip_zeros, skip_nan=False):
    valid = np.isfinite(block)
    keep = np.ones(len(block), dtype=bool)
    if skip_zeros:
        # deepTools uses masked row mean != 0, including signed cancellation.
        keep &= valid.any(axis=(1, 2)) & (np.where(valid, block, 0).sum(axis=(1, 2)) != 0)
    if skip_nan:
        keep &= valid.all(axis=(1, 2))
    return keep


def profile_summary(matrix, rows, groups, keep, batch_rows=2048):
    result = {}
    for group in groups:
        ids = np.array([i for i, r in enumerate(rows) if r[5] == group and keep[i]], dtype=int)
        shape = matrix.shape[1:]
        n = np.zeros(shape, dtype=np.int64)
        mean = np.zeros(shape)
        m2 = np.zeros(shape)
        for pos in range(0, len(ids), batch_rows):
            x = np.asarray(matrix[ids[pos:pos + batch_rows]], dtype=np.float64)
            valid = np.isfinite(x)
            count = valid.sum(axis=0)
            average = np.divide(np.where(valid, x, 0).sum(axis=0), count,
                                out=np.zeros(shape), where=count > 0)
            variance = np.where(valid, x - average, 0)
            variance = (variance * variance).sum(axis=0)
            both = n + count
            delta = average - mean
            fraction = np.divide(count, both, out=np.zeros(shape), where=both > 0)
            m2 += variance + delta * delta * n * fraction
            mean += delta * fraction
            n = both
        mean[n == 0] = np.nan
        se = np.sqrt(np.divide(m2, n * np.maximum(n - 1, 0),
                              out=np.full(shape, np.nan), where=n > 1))
        result[group] = (mean, se, n, len(ids))
    return result


def write_summary(path, summary, samples):
    with open(path, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["group", "sample", "bin", "n_regions", "n_valid", "mean", "sem"])
        for group, (mean, se, n, total) in summary.items():
            for sample, label in enumerate(samples):
                for b in range(mean.shape[1]):
                    writer.writerow([group, label, b + 1, total, n[sample, b], mean[sample, b], se[sample, b]])


def compute(args):
    start_time = time.perf_counter()
    default_flank = 3000 if args.mode == "scale-regions" else 2000
    args.upstream = default_flank if args.upstream is None else args.upstream
    args.downstream = default_flank if args.downstream is None else args.downstream
    cfg = dict(mode=args.mode, reference_point=args.reference_point, upstream=args.upstream,
               downstream=args.downstream, body_length=args.body_length, bin_size=args.bin_size,
               missing_as_zero=args.missing_as_zero, skip_zeros=args.skip_zeros, skip_nan=args.skip_nan)
    nbins = validate_config(cfg)
    if args.batch_rows <= 0 or args.chunk_bp <= 0 or args.merge_gap < 0 or args.processors <= 0:
        raise ValueError("batch_rows/chunk_bp/processors must be positive and merge_gap nonnegative")
    rows, groups = read_beds(args.regions, args.group_labels)
    samples = args.sample_labels or [Path(p).stem for p in args.signal]
    if len(samples) != len(args.signal) or len(set(samples)) != len(samples):
        raise ValueError("Sample labels must be unique and match BigWig files")
    inputs = []
    for p in args.regions + args.signal:
        stat = Path(p).stat()
        inputs.append(dict(path=str(Path(p).resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns))
    out = Path(args.output)
    if out.exists():
        raise ValueError(f"Cache exists: {out}; choose a new output directory")
    out.mkdir(parents=True)
    shape = (len(rows), len(samples), nbins)
    matrix = np.lib.format.open_memmap(out / "matrix.npy", mode="w+", dtype=np.float64, shape=shape)
    if args.processors == 1:
        for sample, bw in enumerate(args.signal):
            fill_sample(matrix, sample, bw, rows, cfg, args.chunk_bp, args.batch_rows, args.merge_gap)
            print(f"Computed {samples[sample]} ({len(rows)} regions)", file=sys.stderr, flush=True)
    else:
        tasks = parallel_tasks(out / "matrix.npy", args.signal, rows, cfg,
                               args.chunk_bp, args.batch_rows, args.merge_gap)
        with ProcessPoolExecutor(max_workers=args.processors) as pool:
            for done, _ in enumerate(pool.map(compute_task, tasks), 1):
                print(f"Computed batch {done}", file=sys.stderr, flush=True)
    matrix.flush()
    keep = np.zeros(len(rows), dtype=bool)
    for pos in range(0, len(rows), args.batch_rows):
        block = matrix[pos:pos + args.batch_rows]
        keep[pos:pos + args.batch_rows] = row_filter(block, args.skip_zeros, args.skip_nan)
    np.save(out / "keep.npy", keep, allow_pickle=False)
    with open(out / "regions.tsv", "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "name", "strand", "group", "kept"])
        writer.writerows([*row, int(k)] for row, k in zip(rows, keep))
    meta = dict(version=VERSION, parameters=cfg, shape=shape, samples=samples, groups=groups,
                inputs=inputs, dtype="float64", n_kept=int(keep.sum()),
                elapsed_seconds=time.perf_counter() - start_time,
                io=dict(chunk_bp=args.chunk_bp, batch_rows=args.batch_rows,
                        merge_gap=args.merge_gap, processors=args.processors),
                status="complete", semantics="deepTools 3.5.1 BED3/6 mean subset")
    if getattr(args, "sample_table", None) is not None:
        meta["sample_table"] = args.sample_table
    if getattr(args, "region_kind", None) is not None:
        meta["region_kind"] = args.region_kind
    np.savez_compressed(out / "metadata.npz", metadata=json.dumps(meta, ensure_ascii=False))
    print(f"Cache saved: {out}")
    return out


def load_cache(path):
    path = Path(path)
    if (path / "metadata.npz").exists():
        with np.load(path / "metadata.npz", allow_pickle=False) as archive:
            meta = json.loads(archive["metadata"].item())
    else:
        # Read JSON metadata without modifying the cache.
        meta = json.loads((path / "metadata.json").read_text())
    if meta.get("status") != "complete":
        raise ValueError("Incomplete cache")
    with open(path / "regions.tsv") as handle:
        rows = [(r["chrom"], int(r["start"]), int(r["end"]), r["name"], r["strand"], r["group"])
                for r in csv.DictReader(handle, delimiter="\t")]
    matrix = np.load(path / "matrix.npy", mmap_mode="r", allow_pickle=False)
    keep = np.load(path / "keep.npy", allow_pickle=False)
    if tuple(meta["shape"]) != matrix.shape or len(rows) != len(matrix) or keep.shape != (len(rows),):
        raise ValueError("Inconsistent cache files")
    return matrix, rows, keep, meta


def axis_ticks(cfg, kind=None):
    before = cfg["upstream"] / cfg["bin_size"]
    total = validate_config(cfg)
    ticks, labels = [], []
    if cfg["upstream"]:
        ticks.append(0)
        labels.append(f'-{cfg["upstream"] / 1000:g} kb')
    ticks.append(before)
    labels.append(("\nStart" if cfg["upstream"] else "Start")
                  if cfg["mode"] == "scale-regions" else (
                      {"te": "TE center", "peak": "Peak center"}.get(kind, "Center")
                      if cfg["reference_point"] == "center" else cfg["reference_point"]))
    if cfg["mode"] == "scale-regions":
        ticks.append(before + cfg["body_length"] / cfg["bin_size"])
        labels.append("\nEnd" if cfg["downstream"] else "End")
    if cfg["downstream"]:
        ticks.append(total)
        labels.append(f'+{cfg["downstream"] / 1000:g} kb')
    return ticks, labels


def sample_records(samples):
    """Accept a records list or DataFrame; stable IDs are distinct from labels."""
    records = samples.to_dict("records") if hasattr(samples, "to_dict") else list(samples)
    if not records or any(not r.get("sample") or not r.get("group") for r in records):
        raise ValueError("Samples need nonempty sample and group columns")
    if len({r["sample"] for r in records}) != len(records):
        raise ValueError("Sample IDs must be unique")
    return [dict(r, path=str(r["path"])) if "path" in r else dict(r) for r in records]


def mode_configs(kind="te", modes=None, window=None):
    """Public modes describe the view; TSS center still anchors the gene's 5' end."""
    if kind not in ("te", "tss", "peak", "region"):
        raise ValueError("kind must be te/tss/peak/region")
    modes = (["center"] if kind == "tss" else ["center", "body"]) if modes is None else modes
    modes = [modes] if isinstance(modes, str) else list(modes)
    if not modes or len(set(modes)) != len(modes) or set(modes) - {"center", "body"}:
        raise ValueError("modes must contain center and/or body without duplicates")
    if kind == "tss" and modes != ["center"]:
        raise ValueError("TSS only supports center; use kind='region' for gene bodies")
    if set(window or {}) - {"upstream", "downstream", "body_length", "bin_size"}:
        raise ValueError("window supports upstream/downstream/body_length/bin_size")
    result = {}
    for mode in modes:
        cfg = dict(mode="scale-regions" if mode == "body" else "reference-point",
                   reference_point="TSS" if kind == "tss" else "center",
                   upstream=2000 if kind == "tss" else 3000,
                   downstream=2000 if kind == "tss" else 3000,
                   body_length=5000, bin_size=50,
                   missing_as_zero=True, skip_zeros=False, skip_nan=False)
        cfg.update(window or {})
        if mode == "center":
            cfg["body_length"] = 0
        validate_config(cfg)
        result[mode] = cfg
    return result


def compute_profiles(samples, regions, cache, *, kind="te", modes=None, workers=4, window=None):
    """Return {mode: cache_path}, preserving every BED row in each mode.

    Default TE/peak/region: center + body. TSS: center anchored to TSS.
    Filtering and common sorting happen in make_plots, without re-reading BigWigs.
    """
    records = sample_records(samples)
    if kind == "tss" and any(r[4] == "." for r in read_beds(list(regions.values()), list(regions))[0]):
        raise ValueError("TSS requires stranded gene BED")
    return {mode: _compute_profile_mode(records, regions, Path(cache) / mode, kind, cfg, workers)
            for mode, cfg in mode_configs(kind, modes, window).items()}


def _compute_profile_mode(records, regions, cache, kind, cfg, workers):
    labels, beds = list(regions), list(map(str, regions.values()))
    bws, names = [r["path"] for r in records], [r["sample"] for r in records]
    cache = Path(cache)
    if cache.exists():
        _, old_rows, _, meta = load_cache(cache)
        rows, _ = read_beds(beds, labels)
        fingerprints = [dict(path=str(Path(p).resolve()), size=Path(p).stat().st_size,
                             mtime_ns=Path(p).stat().st_mtime_ns) for p in beds + bws]
        if (meta["parameters"] != cfg or meta.get("region_kind") != kind
                or meta["samples"] != names or old_rows != rows
                or meta["inputs"] != fingerprints):
            raise ValueError("Cache inputs/parameters changed; choose a new cache directory")
        return cache
    args = argparse.Namespace(**cfg, signal=bws, regions=beds, group_labels=labels,
        sample_labels=names, sample_table=records, region_kind=kind, output=str(cache), processors=workers,
        chunk_bp=2000000, batch_rows=2048, merge_gap=10000)
    return compute(args)


def nice_color_scale(upper):
    """Round a positive upper limit upward; return integer ticks + one multiplier."""
    if not np.isfinite(upper) or upper <= 0:
        raise ValueError("Color upper limit must be finite and positive")
    exponent = int(np.floor(np.log10(upper)))
    unit = 10.0 ** exponent
    top = int(np.ceil(upper / unit - 1e-12))
    step = 1 if top <= 5 else 2
    top = int(np.ceil(top / step)) * step
    integers = np.arange(0, top + 1, step)
    return float(top * unit), integers * unit, [str(x) for x in integers], exponent


def prepare_heatmap(values, quantile=.99, normalization="shared_max", clim=None):
    """Normalize complete selected data before display aggregation; never clip it."""
    if not 0 < quantile <= 1:
        raise ValueError("quantile must be in (0, 1]")
    finite = values[np.isfinite(values)]
    if not len(finite):
        raise ValueError("No finite signal in selected regions")
    maximum = float(finite.max())
    if normalization == "shared_max":
        if finite.min() < 0:
            raise ValueError("shared_max requires nonnegative signal; use normalization='none' for signed tracks")
        divisor = maximum if maximum > 0 else 1.0
    elif normalization == "none":
        divisor = 1.0
    else:
        raise ValueError("normalization must be 'shared_max' or 'none'")
    normalized = values / divisor
    signed = finite.min() < 0
    threshold = float(np.quantile(np.abs(finite) if signed else finite, quantile)) / divisor
    fallback = float(np.abs(finite).max()) / divisor if signed else maximum / divisor
    upper, ticks, labels, exponent = nice_color_scale(threshold if threshold > 0 else (fallback or 1.0))
    lower = 0.0
    if signed:
        lower = -upper
        ticks = np.concatenate((-ticks[:0:-1], ticks))
        labels = [f"{x / (10.0 ** exponent):.0f}" for x in ticks]
    if clim is not None:
        lower, upper = map(float, clim)
        if not np.isfinite([lower, upper]).all() or lower >= upper:
            raise ValueError("clim requires two finite increasing bounds")
        ticks = np.linspace(lower, upper, 5)
        labels = [f"{x:g}" for x in ticks]
        exponent = 0
    info = dict(normalization=normalization, shared_max=maximum, divisor=divisor,
        normalization_scope="all selected loci, bins and samples within this figure",
        quantile=quantile, percentile_value=threshold, color_limits=[lower, upper],
        color_ticks=ticks.tolist(), color_tick_labels=labels, color_exponent=exponent,
        color_scale="linear", colorbar_extend="neither",
        saturated_fraction=float(np.mean((finite / divisor < lower) | (finite / divisor > upper))))
    return normalized, info


def paired_rows(caches, *, regions=None, sort_sample=None, filter_all_zero=False):
    """Validate paired raw caches and derive one keep mask/order for each view."""
    if not caches or set(caches) - {"center", "body"}:
        raise ValueError("Expected a center/body cache mapping")
    loaded = {mode: load_cache(path) for mode, path in caches.items()}
    reference = "center" if "center" in loaded else "body"
    matrix, rows, _, meta = loaded[reference]
    sort_sample = sort_sample or meta["samples"][0]
    if sort_sample not in meta["samples"]:
        raise ValueError("sort_sample must be a cached sample ID")
    nonzero = {}
    for mode, (values, mode_rows, keep, config) in loaded.items():
        if mode_rows != rows or config["samples"] != meta["samples"]:
            raise ValueError("Modes must have identical BED rows and sample order")
        if config["inputs"] != meta["inputs"]:
            raise ValueError("Modes must use the same BED and BigWig input fingerprints")
        if not config["parameters"]["missing_as_zero"] or not keep.all():
            raise ValueError("Paired modes require filled missing values and unfiltered caches")
        expected_mode = "reference-point" if mode == "center" else "scale-regions"
        if config["parameters"]["mode"] != expected_mode:
            raise ValueError("Cache mode label disagrees with computation parameters")
        flags = np.zeros(len(rows), dtype=bool)
        for start in range(0, len(rows), 2048):
            block = values[start:start + 2048]
            if not np.isfinite(block).all():
                raise ValueError("Paired matrices must be finite after missing-as-zero")
            flags[start:start + 2048] = np.any(block != 0, axis=(1, 2))
        nonzero[mode] = flags
    any_signal = np.logical_or.reduce(list(nonzero.values()))
    views = regions if regions is not None else {
        g: np.array([r[5] == g for r in rows]) for g in meta["groups"]}
    result = {}
    for name, mask in views.items():
        mask = np.asarray(mask)
        if mask.dtype != bool or mask.shape != (len(rows),):
            raise ValueError("Region masks must be Boolean, in original BED row order")
        keep = mask & any_signal if filter_all_zero else mask.copy()
        ids = np.flatnonzero(keep)
        scores = matrix[ids, meta["samples"].index(sort_sample), :].mean(axis=1)
        order = ids[np.argsort(-scores, kind="stable")]
        counts = dict(selected_count=int(mask.sum()), kept_count=int(keep.sum()),
                      removed_count=int(mask.sum() - keep.sum()),
                      all_modes_zero=int((mask & ~any_signal).sum()),
                      filter_all_zero=bool(filter_all_zero))
        if set(nonzero) == {"center", "body"}:
            c, b = nonzero["center"], nonzero["body"]
            counts.update(center_only=int((mask & c & ~b).sum()),
                          body_only=int((mask & b & ~c).sum()),
                          both_nonzero=int((mask & c & b).sum()))
        result[name] = dict(keep=keep, order=order, counts=counts,
                            sort_mode=reference, sort_sample=sort_sample)
    return result


def make_plots(cache, *, samples=None, regions=None, style=None, filter_all_zero=False):
    """A mode mapping returns {mode: {view: figures}} with shared rows and order.

    A single cache path retains the single-mode return layout for CLI consumers.
    Zero rows are retained by default; optional filtering removes only joint zeros.
    """
    if not isinstance(cache, dict):
        if filter_all_zero:
            raise ValueError("Use a mode mapping to request post-compute zero filtering")
        return _make_single_plots(cache, samples=samples, regions=regions, style=style)
    records = sample_records(samples) if samples is not None else None
    sort_sample = (style or {}).get("sort_sample") or (records[0]["sample"] if records else None)
    layout = paired_rows(cache, regions=regions, sort_sample=sort_sample, filter_all_zero=filter_all_zero)
    output = {}
    for mode, path in cache.items():
        output[mode] = _make_single_plots(path, samples=samples, regions=regions, style=style,
                                          layout=layout)
        for name, item in output[mode].items():
            item["info"].update(mode=mode, sort_mode=layout[name]["sort_mode"],
                                  **layout[name]["counts"])
    return output


def _make_single_plots(cache, *, samples=None, regions=None, style=None, layout=None):
    """Return {region_set: {profile, heatmap, info}} without writing any files.

    samples supplies sample/group and optional label (path optional for plotting).
    regions optionally maps view names to Boolean masks in original cache row order.
    Each view gets a shared normalization/scale across all samples, its own figure pair.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    defaults = dict(cmap=SIGNAL_BLUE, colors=None, profile_size=(60, 45), heatmap_size=None,
        quantile=.99, sort_sample=None, max_rows=2000, signal_label="BigWig signal",
        sem=False, normalization="shared_max", clim=None, title=None)
    if set(style or {}) - set(defaults):
        raise ValueError(f"Unknown style keys: {set(style) - set(defaults)}")
    settings = dict(defaults, **(style or {}))
    matrix, rows, keep, meta = load_cache(cache)
    records = sample_records(samples if samples is not None else meta.get("sample_table", [
        dict(sample=s, group=s) for s in meta["samples"]]))
    if set(r["sample"] for r in records) != set(meta["samples"]):
        raise ValueError("Sample table must match all cached sample IDs exactly")
    # Stable group order follows first occurrence, never averages biological replicates.
    group_order = list(dict.fromkeys(r["group"] for r in records))
    records = [r for g in group_order for r in records if r["group"] == g]
    sample_order = [meta["samples"].index(r["sample"]) for r in records]
    sort_sample = settings["sort_sample"] or records[0]["sample"]
    if sort_sample not in meta["samples"]:
        raise ValueError("sort_sample must be a cached sample ID")
    sort_index = [r["sample"] for r in records].index(sort_sample)
    max_rows = settings["max_rows"]
    if not isinstance(max_rows, int) or max_rows <= 0:
        raise ValueError("max_rows must be a positive integer")
    colors = dict(zip(group_order, [GROUP_COLORS[i % len(GROUP_COLORS)] for i in range(len(group_order))]))
    colors.update(settings["colors"] or {})
    palette = settings["cmap"]
    palette = (LinearSegmentedColormap.from_list("signal_blue", palette)
               if isinstance(palette, (list, tuple)) else copy.copy(plt.get_cmap(palette)))
    palette.set_bad("0.85")
    views = regions if regions is not None else {
        g: np.array([r[5] == g for r in rows]) for g in meta["groups"]}
    ticks, tick_labels = axis_ticks(meta["parameters"], meta.get("region_kind"))
    nbins = matrix.shape[2]
    width, height = settings["heatmap_size"] or (20 + 35 * len(records), 80)
    pw, ph = settings["profile_size"]
    if min(width, height, pw, ph) <= 0:
        raise ValueError("Figure dimensions must be positive")
    rc = {"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
          "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
          "axes.linewidth": .5, "pdf.fonttype": 42}
    result = {}
    with plt.rc_context(rc):
        for name, mask in views.items():
            mask = np.asarray(mask)
            if mask.dtype != bool or mask.shape != keep.shape:
                raise ValueError("Region masks must be Boolean, in original cache row order")
            ids = np.flatnonzero(layout[name]["keep"] if layout is not None else mask & keep)
            if not len(ids):
                continue
            values = np.asarray(matrix[ids][:, sample_order, :])
            valid = np.isfinite(values)
            n_valid = valid.sum(axis=0)
            mean = np.divide(np.where(valid, values, 0).sum(axis=0), n_valid,
                             out=np.full(values.shape[1:], np.nan), where=n_valid > 0)
            residual = np.where(valid, values - mean, 0)
            sem = np.sqrt(np.divide((residual * residual).sum(axis=0), n_valid * (n_valid - 1),
                                   out=np.zeros_like(mean), where=n_valid > 1))
            sort_values = values[:, sort_index]
            counts = np.isfinite(sort_values).sum(axis=1)
            scores = np.divide(np.where(np.isfinite(sort_values), sort_values, 0).sum(axis=1), counts,
                               out=np.full(len(ids), -np.inf), where=counts > 0)
            order = (np.searchsorted(ids, layout[name]["order"]) if layout is not None
                     else np.argsort(-scores, kind="stable"))
            normalized, info = prepare_heatmap(values, settings["quantile"], settings["normalization"], settings["clim"])
            normalized = normalized[order]
            edges = np.linspace(0, len(ids), min(len(ids), max_rows) + 1, dtype=int)
            shown = []
            for a, b in zip(edges[:-1], edges[1:]):
                block = normalized[a:b]
                count = np.isfinite(block).sum(axis=0)
                shown.append(np.divide(np.where(np.isfinite(block), block, 0).sum(axis=0), count,
                                       out=np.full(block.shape[1:], np.nan), where=count > 0))
            shown = np.stack(shown)
            title = f"{settings['title']} · {name}" if settings["title"] else str(name)
            profile, ax = plt.subplots(figsize=(pw / 25.4, ph / 25.4), constrained_layout=True)
            seen = {}
            for s, record in enumerate(records):
                g = record["group"]
                repeat = seen.get(g, 0)
                seen[g] = repeat + 1
                color = colors.get(record["sample"], colors[g])
                ax.plot(np.arange(nbins) + .5, mean[s], label=record.get("label", record["sample"]),
                        color=color, lw=.8, ls=["-", "--", ":", "-."][repeat % 4])
                if settings["sem"]:
                    ax.fill_between(np.arange(nbins) + .5, mean[s] - sem[s], mean[s] + sem[s],
                                    color=color, alpha=.13, linewidth=0)
            for boundary in ticks:
                if 0 < boundary < nbins:
                    ax.axvline(boundary, color=".75", lw=.5, ls="--", zorder=0)
            ax.set(title=title, ylabel=settings["signal_label"], xlim=(0, nbins),
                   xticks=ticks, xticklabels=tick_labels)
            if meta["parameters"]["mode"] == "scale-regions":
                ax.set_xlabel("Scaled body")
            if np.nanmin(values) >= 0:
                ax.set_ylim(bottom=0)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.tick_params(length=2.5, width=.5)
            ax.legend(frameon=False, loc="best")
            heatmap, axes = plt.subplots(1, len(records), squeeze=False,
                figsize=(width / 25.4, height / 25.4), constrained_layout=True)
            for s, ax in enumerate(axes[0]):
                image = ax.imshow(shown[:, s], aspect="auto", interpolation="nearest", cmap=palette,
                    vmin=info["color_limits"][0], vmax=info["color_limits"][1],
                    extent=(0, nbins, len(ids), 0), rasterized=True)
                ax.set(title=records[s].get("label", records[s]["sample"]),
                       xticks=ticks, xticklabels=tick_labels)
                if meta["parameters"]["mode"] == "scale-regions":
                    ax.set_xlabel("Scaled body")
                ax.set_yticks([0, len(ids)])
                ax.set_yticklabels(["0", f"{len(ids):,}"] if s == 0 else [])
                if s == 0:
                    ax.set_ylabel("Loci (display bins)" if len(shown) < len(ids) else "Loci")
                ax.tick_params(length=2.5, width=.5)
            bar = heatmap.colorbar(image, ax=axes.ravel().tolist(), shrink=.65, pad=.025,
                                   ticks=info["color_ticks"], extend="neither")
            bar.ax.set_yticklabels(info["color_tick_labels"])
            if info["color_exponent"]:
                bar.ax.set_title(r"$\times10^{%d}$" % info["color_exponent"], fontsize=6, pad=5)
            bar.set_label("Normalized signal" if settings["normalization"] == "shared_max" else settings["signal_label"])
            bar.ax.tick_params(labelsize=6, length=2.5, width=.5)
            bar.outline.set_linewidth(.5)
            heatmap.suptitle(title, fontsize=7)
            info.update(region_set=str(name), selected_count=int(mask.sum()), kept_count=len(ids),
                sample_ids=[r["sample"] for r in records], sample_labels=[r.get("label", r["sample"]) for r in records],
                sample_groups=[r["group"] for r in records], sort_sample=sort_sample,
                row_indices=ids[order].tolist(), display_rows=len(shown),
                display_edges=edges.tolist(),
                profile_size_mm=[pw, ph], heatmap_size_mm=[width, height],
                signal_label=settings["signal_label"], colors=colors,
                cmap=settings["cmap"], window=meta["parameters"],
                uncertainty="cross-region SEM" if settings["sem"] else None)
            plt.close(profile)
            plt.close(heatmap)
            result[name] = dict(profile=profile, heatmap=heatmap, info=info,
                                data=dict(mean=mean, sem=sem, n_valid=n_valid))
    if not result:
        raise ValueError("No regions survive selection/filtering")
    return result


def save_plots(plots, output, *, dpi=600):
    """Save previewed objects as PDF + PNG only; keep export dpi in info."""
    import matplotlib.pyplot as plt

    output = Path(output)
    if all("profile" not in item for item in plots.values()):
        for mode, views in plots.items():
            if mode not in ("center", "body"):
                raise ValueError("Expected center/body figure mapping")
            save_plots(views, output / mode, dpi=dpi)
        return
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"pdf.fonttype": 42}):
        for name, item in plots.items():
            if not str(name) or Path(str(name)).name != str(name) or str(name) in (".", ".."):
                raise ValueError("Region-set names must be simple output basenames")
            for feature in ("profile", "heatmap"):
                figure = item[feature]
                for suffix in ("pdf", "png"):
                    figure.savefig(output / f"{name}_{feature}.{suffix}", dpi=dpi)
            item["info"]["dpi"] = dpi


def save_profile_data(plots, caches, output, *, normalized=False):
    """Export plotted summaries and a shared row table; optionally full normalized arrays."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    common = {}
    for mode, views in plots.items():
        matrix, rows, _, meta = load_cache(caches[mode])
        directory = output / mode
        directory.mkdir(parents=True, exist_ok=True)
        for name, item in views.items():
            if Path(str(name)).name != str(name) or name in (".", ".."):
                raise ValueError("View names must be simple output basenames")
            info, data = item["info"], item["data"]
            ids = np.array(info["row_indices"], dtype=int)
            if name in common:
                if common[name]["row_indices"] != info["row_indices"]:
                    raise ValueError("Modes have inconsistent export row order")
            else:
                common[name] = info
                keep = np.zeros(len(rows), dtype=bool)
                keep[ids] = True
                np.save(output / f"{name}_keep.npy", keep, allow_pickle=False)
                with (output / f"{name}_rows.tsv").open("w") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["row_index", "chrom", "start", "end", "name", "strand", "group"])
                    writer.writerows([[int(i), *rows[i]] for i in ids])
            summary = {name: (data["mean"], data["sem"], data["n_valid"], len(ids))}
            write_summary(directory / f"{name}_profile.tsv", summary, info["sample_ids"])
            if normalized:
                sample_ids = [meta["samples"].index(s) for s in info["sample_ids"]]
                values = np.array(matrix[ids][:, sample_ids, :])
                if info["normalization"] == "shared_max" and info["shared_max"] > 0:
                    values /= info["shared_max"]
                np.save(directory / f"{name}_heatmap_01.npy", values, allow_pickle=False)


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    _, _, _, meta = load_cache(args.matrix)
    if not 0 <= args.sort_sample < len(meta["samples"]):
        raise ValueError("sort-sample is outside cached samples")
    style = dict(sort_sample=meta["samples"][args.sort_sample], max_rows=args.max_rows,
                 signal_label=args.signal_label, sem=args.sem, quantile=args.quantile,
                 normalization=args.normalization)
    if args.cmap is not None:
        style["cmap"] = args.cmap
    if args.width_mm is not None or args.height_mm is not None:
        style["heatmap_size"] = (args.width_mm or (20 + 35 * len(meta["samples"])), args.height_mm or 80)
    if args.vmin is not None or args.vmax is not None:
        if args.vmin is None or args.vmax is None:
            raise ValueError("Provide --vmin and --vmax together")
        style["clim"] = (args.vmin, args.vmax)
    figures = make_plots(args.matrix, style=style)
    save_plots(figures, args.output)
    print(f"Figures saved: {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    sub = parser.add_subparsers(dest="command", required=True)
    comp = sub.add_parser("compute", help="Build a reusable region × sample × bin matrix")
    comp.add_argument("-S", "--signal", nargs="+", required=True)
    comp.add_argument("-R", "--regions", nargs="+", required=True)
    comp.add_argument("-o", "--output", required=True, help="New cache directory")
    comp.add_argument("--sample-labels", nargs="+")
    comp.add_argument("--group-labels", nargs="+")
    comp.add_argument("--mode", choices=["scale-regions", "reference-point"], default="scale-regions")
    comp.add_argument("--reference-point", choices=["TSS", "TES", "center"], default="TSS")
    comp.add_argument("-b", "--upstream", type=int, default=None)
    comp.add_argument("-a", "--downstream", type=int, default=None)
    comp.add_argument("-m", "--body-length", type=int, default=5000)
    comp.add_argument("-bs", "--bin-size", type=int, default=50)
    comp.add_argument("--missingDataAsZero", dest="missing_as_zero", action="store_true", default=True)
    comp.add_argument("--keep-missing", dest="missing_as_zero", action="store_false")
    comp.add_argument("--skipZeros", dest="skip_zeros", action="store_true", default=False)
    comp.add_argument("--keep-zero-rows", dest="skip_zeros", action="store_false")
    comp.add_argument("--skip-nan", action="store_true", help="Additionally remove rows with any missing bin")
    comp.add_argument("--chunk-bp", type=int, default=2000000)
    comp.add_argument("--batch-rows", type=int, default=2048)
    comp.add_argument("--merge-gap", type=int, default=10000, help="Maximum gap between windows read together")
    comp.add_argument("-p", "--processors", type=int, default=1, help="Number of worker processes")
    comp.set_defaults(func=compute)
    pl = sub.add_parser("plot", help="Plot cached matrix without reading BigWigs")
    pl.add_argument("--matrix", required=True)
    pl.add_argument("-o", "--output", required=True, help="Figure directory")
    pl.add_argument("--sort-sample", type=int, default=0, help="Zero-based sample index; shared row order")
    pl.add_argument("--max-rows", type=int, default=2000, help="Max displayed rows per group; numeric matrix stays full")
    pl.add_argument("--cmap", help="Matplotlib palette; default is the shared signal-blue palette")
    pl.add_argument("--quantile", type=float, default=.99)
    pl.add_argument("--normalization", choices=["shared_max", "none"], default="shared_max")
    pl.add_argument("--vmin", type=float)
    pl.add_argument("--vmax", type=float)
    pl.add_argument("--width-mm", type=float)
    pl.add_argument("--height-mm", type=float)
    pl.add_argument("--signal-label", default="BigWig signal")
    pl.add_argument("--sem", action="store_true")
    pl.set_defaults(func=plot)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
