#!/usr/bin/env python3
"""Numerical oracle tests against installed deepTools; all output in --workdir.

Run in a Python environment with deepTools and the runtime dependencies installed.
Optional --benchmark-bed/--benchmark-bw compare actual computeMatrix wall time.
"""
import argparse
import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pyBigWig
from deeptools.heatmapper import heatmapper

import fast_chip_profiles as fast


def run_deeptools(work, bed, bws, cfg, name):
    target = work / (name + ".gz")
    command = [sys.executable, str(Path(sys.executable).parent / "computeMatrix"), cfg["mode"],
               "-R", str(bed), "-S", *map(str, bws), "-b", str(cfg["upstream"]),
               "-a", str(cfg["downstream"]), "-bs", str(cfg["bin_size"]),
               "--sortRegions", "keep", "-p", "1", "-o", str(target), "--quiet"]
    if cfg["mode"] == "scale-regions":
        command += ["-m", str(cfg["body_length"])]
    else:
        command += ["--referencePoint", cfg["reference_point"]]
    if cfg["missing_as_zero"]:
        command.append("--missingDataAsZero")
    if cfg["skip_zeros"]:
        command.append("--skipZeros")
    start = time.perf_counter()
    proc = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(proc.stderr)
    records = {}
    with gzip.open(target, "rt") as handle:
        header = json.loads(next(handle)[1:])
        for line in handle:
            fields = line.rstrip().split("\t")
            records[fields[3]] = np.asarray(fields[6:], dtype=float).reshape(len(bws), -1)
    return records, elapsed, header


def compare(work, bed, bws, cfg, name):
    rows, groups = fast.read_beds([bed])
    matrix = np.empty((len(rows), len(bws), fast.validate_config(cfg)))
    start = time.perf_counter()
    for sample, bw in enumerate(bws):
        fast.fill_sample(matrix, sample, bw, rows, cfg, 2000000, 2048)
    keep = fast.row_filter(matrix, cfg["skip_zeros"])
    fast_time = time.perf_counter() - start
    truth, dt_time, _ = run_deeptools(work, bed, bws, cfg, name)
    names = [r[3] for r, k in zip(rows, keep) if k]
    assert set(names) == set(truth), (name, "filtered row names differ", set(names) ^ set(truth))
    error = 0.0
    for row, values, included in zip(rows, matrix, keep):
        if not included:
            continue
        expected = truth[row[3]]
        np.testing.assert_array_equal(np.isnan(values), np.isnan(expected))
        np.testing.assert_allclose(values, expected, rtol=1e-7, atol=5.1e-7, equal_nan=True)
        valid = np.isfinite(expected)
        if valid.any():
            error = max(error, float(np.max(np.abs(values[valid] - expected[valid]))))
    result = dict(case=name, n_regions=len(rows), n_kept=int(keep.sum()),
                  fast_core_seconds=fast_time, deeptools_cli_seconds=dt_time,
                  max_abs_error_vs_6decimal_export=error)
    print(json.dumps(result), flush=True)
    return result


def synthetic(work):
    n = 20000
    base = np.arange(n)
    bw_paths = []
    for sample in range(2):
        path = work / f"sample{sample}.bw"
        values = np.sin(base / 31.0) + 2 + sample
        values[6000:8000] = 0
        if sample == 1:
            values[6000:7000] = 1
        valid = (base < 4000) | (base >= 4500)
        with pyBigWig.open(str(path), "w") as bw:
            bw.addHeader([("chrT", n)])
            idx = base[valid]
            bw.addEntries(["chrT"] * len(idx), idx.tolist(), ends=(idx + 1).tolist(), values=values[valid].tolist())
        bw_paths.append(path)
    bed = work / "synthetic.bed"
    rows = []
    for i, length in enumerate([1, 49, 50, 51, 70, 99, 100, 121, 200, 507]):
        for strand in ["+", "-"]:
            start = 1000 + i * 510 + (0 if strand == "+" else 11)
            rows.append((start, start + length, strand))
    rows += [(0, 77, "+"), (5, 116, "-"), (19800, 20000, "+"), (19810, 19999, "-"),
             (6400, 6600, "+"), (7400, 7600, "-"), (4100, 4250, "+"),
             (3990, 4503, "-"), (7990, 8015, "+")]
    with open(bed, "w") as handle:
        for i, (start, end, strand) in enumerate(rows):
            handle.write(f"chrT\t{start}\t{end}\tr{i}\t0\t{strand}\n")
    cfg = dict(mode="scale-regions", reference_point="TSS", upstream=100, downstream=150,
               body_length=5000, bin_size=50, missing_as_zero=True, skip_zeros=True)
    results = []
    for mode, point in [("scale-regions", "TSS"), ("reference-point", "TSS"),
                        ("reference-point", "TES"), ("reference-point", "center")]:
        for missing in [True, False]:
            for skip in [True, False]:
                test_cfg = dict(cfg, mode=mode, reference_point=point, missing_as_zero=missing, skip_zeros=skip)
                name = f"{mode}_{point}_zero{int(missing)}_skip{int(skip)}"
                results.append(compare(work, bed, bw_paths, test_cfg, name))
    # Signed tracks: match deepTools' mean-zero rule; do not remove zero bins.
    signed = np.array([[[1.0, -1.0]], [[0, 1]], [[np.nan, 2]], [[np.nan, np.nan]], [[0, 0]]])
    np.testing.assert_array_equal(fast.row_filter(signed, True), [False, True, True, False, False])
    np.testing.assert_array_equal(fast.row_filter(signed, True, True), [False, True, False, False, False])
    # Test raw bin means directly (no 6-decimal text rounding) against library.
    test_rows, _ = fast.read_beds([bed])
    valid_rows = [r for r in test_rows if r[2] - r[1] >= 50 and r[1] >= 150 and r[2] <= n - 150]
    matrix = np.empty((len(valid_rows), 1, fast.validate_config(cfg)))
    fast.fill_sample(matrix, 0, bw_paths[0], valid_rows, cfg, 2000000, 2048)
    with pyBigWig.open(str(bw_paths[0])) as bw:
        for i, r in enumerate(valid_rows):
            left, right = (150, 100) if r[4] == "-" else (100, 150)
            zones = [([(r[1] - left, r[1])], left // 50), ([(r[1], r[2])], 100),
                     ([(r[2], r[2] + right)], right // 50)]
            expected = heatmapper.coverage_from_big_wig(bw, "chrT", zones, 50, "mean", True, False)
            if r[4] == "-":
                expected = expected[::-1]
            np.testing.assert_allclose(matrix[i, 0], expected, atol=1e-10, rtol=1e-10)
    return results, bed, bw_paths


def validate_plot_api(work, bed, bws, workers):
    samples = [dict(path=str(bws[i]), sample=f"s{i}", group="WT", label=f"WT rep{i+1}") for i in range(2)]
    caches = fast.compute_profiles(samples, {"TE": bed}, work / "api_te", workers=workers)
    cache = caches["body"]
    matrix, rows, keep, meta = fast.load_cache(cache)
    assert {p.name for p in cache.iterdir()} == {"matrix.npy", "keep.npy", "regions.tsv", "metadata.npz"}
    assert matrix.shape[2] == 220 and meta["parameters"]["upstream"] == 3000
    stamp = (cache / "matrix.npy").stat().st_mtime_ns
    assert fast.compute_profiles(samples, {"TE": bed}, work / "api_te") == caches
    assert (cache / "matrix.npy").stat().st_mtime_ns == stamp
    try:
        fast.compute_profiles(samples, {"TE": bed}, work / "api_te", window={"upstream": 1000})
    except ValueError:
        pass
    else:
        raise AssertionError("Stale cache was reused")
    tss = fast.compute_profiles(samples, {"TSS": bed}, work / "api_tss", kind="tss", workers=1)
    tm, _, _, tmeta = fast.load_cache(tss["center"])
    assert tm.shape[2] == 80 and tmeta["parameters"]["upstream"] == tmeta["parameters"]["downstream"] == 2000

    values = np.ones((100, 2, 10))
    values[:, 1] = 2
    values[0, 1, 0] = 1000
    original = values.copy()
    normalized, info = fast.prepare_heatmap(values)
    np.testing.assert_array_equal(values, original)
    np.testing.assert_array_equal(normalized, values / 1000)
    assert normalized.max() == 1 and info["color_limits"] == [0, .002]
    assert info["color_tick_labels"] == ["0", "1", "2"] and info["color_exponent"] == -3
    assert info["saturated_fraction"] > 0
    normalized, _ = fast.prepare_heatmap(np.array([[[np.nan, 0, 2]]]))
    assert np.isnan(normalized[0, 0, 0]) and normalized[0, 0, 2] == 1
    normalized, info = fast.prepare_heatmap(np.zeros((2, 2, 3)))
    assert (normalized == 0).all() and info["color_limits"] == [0, 1]
    sparse = np.zeros((100, 2, 10)); sparse[0, 0, 0] = 5
    _, info = fast.prepare_heatmap(sparse)
    assert info["percentile_value"] == 0 and info["color_limits"] == [0, 1]
    try:
        fast.prepare_heatmap(np.array([[[-1., 1.]]]))
    except ValueError:
        pass
    else:
        raise AssertionError("Signed signal silently max-normalized")
    raw, _ = fast.prepare_heatmap(np.array([[[-1., 1.]]]), normalization="none")
    np.testing.assert_array_equal(raw, [[[-1., 1.]]])

    mask = np.arange(len(rows)) % 2 == 0
    figures = fast.make_plots(cache, samples=samples[::-1], regions={"subset": mask},
        style={"colors": {"WT": "#4A927D"}, "max_rows": 3,
               "cmap": ["white", "#4A927D"], "sem": True})
    item = figures["subset"]
    info = item["info"]
    assert info["sample_ids"] == ["s1", "s0"] and info["sample_labels"] == ["WT rep2", "WT rep1"]
    expected = matrix[keep & mask][:, [1, 0], :]
    profile_lines = item["profile"].axes[0].lines[:2]
    for sample, line in enumerate(profile_lines):
        np.testing.assert_allclose(line.get_ydata(), expected[:, sample].mean(axis=0))
        assert line.get_color() == "#4A927D"
    assert profile_lines[0].get_linestyle() != profile_lines[1].get_linestyle()
    assert info["shared_max"] == expected.max()
    np.testing.assert_allclose(info["percentile_value"], np.quantile(expected / expected.max(), .99))
    ids = np.flatnonzero(keep & mask)
    order = np.argsort(-expected[:, 0].mean(axis=1), kind="stable")
    np.testing.assert_array_equal(info["row_indices"], ids[order])
    assert info["display_rows"] == 3
    image = item["heatmap"].axes[-2].images[0]
    assert image.colorbar.extend == "neither"
    assert image.get_clim() == tuple(info["color_limits"])
    np.testing.assert_allclose(item["profile"].get_size_inches() * 25.4, [60, 45])
    np.testing.assert_allclose(item["heatmap"].get_size_inches() * 25.4, [90, 80])
    fast.save_plots(figures, work / "api_figs")
    settings = figures["subset"]["info"]
    assert settings["dpi"] == 600 and settings["colorbar_extend"] == "neither"
    assert {p.suffix for p in (work / "api_figs").iterdir()} == {".pdf", ".png"}
    assert all((work / "api_figs" / f"subset_{feature}.{suffix}").stat().st_size > 1000
               for feature in ["profile", "heatmap"] for suffix in ["pdf", "png"])
    try:
        fast.make_plots(cache, regions={"bad": np.ones(2, dtype=bool)})
    except ValueError:
        pass
    else:
        raise AssertionError("Misaligned mask accepted")
    assert keep.all() and not meta["parameters"]["skip_zeros"]
    assert set(caches) == {"center", "body"} and set(tss) == {"center"}
    assert tmeta["parameters"]["reference_point"] == "TSS"
    validate_paired(work, caches, samples)
    assert not list(work.rglob("*.json"))
    print("PASS: public API, defaults, reuse, normalization, colorbar, samples, sorting and exports")


def validate_paired(work, caches, samples):
    # Deliberately different center/body signals expose independent filtering/sorting.
    for mode, cache in caches.items():
        values = np.load(cache / "matrix.npy", mmap_mode="r+")
        values[:] = 0
        if mode == "center":
            values[0, 0, :] = 3
            values[2, 0, :] = 3
            values[4, 1, :] = 2
        else:
            values[1, 1, :] = 7
            values[2, 0, :] = 1
            values[4, 0, :] = 4
        values.flush()
    mask = np.arange(len(values)) < 5
    layout = fast.paired_rows(caches, regions={"view": mask}, sort_sample="s0")
    stats = layout["view"]["counts"]
    assert stats["center_only"] == 1 and stats["body_only"] == 1
    assert stats["both_nonzero"] == 2 and stats["all_modes_zero"] == 1
    assert stats["kept_count"] == 5 and stats["removed_count"] == 0
    np.testing.assert_array_equal(layout["view"]["order"], [0, 2, 1, 3, 4])
    filtered = fast.paired_rows(caches, regions={"view": mask}, filter_all_zero=True)
    np.testing.assert_array_equal(filtered["view"]["order"], [0, 2, 1, 4])
    figures = fast.make_plots(caches, samples=samples, regions={"view": mask},
                              style={"max_rows": 2, "sort_sample": "s0"})
    for mode, views in figures.items():
        item = views["view"]
        np.testing.assert_array_equal(item["info"]["row_indices"], [0, 2, 1, 3, 4])
        assert item["info"]["display_edges"] == [0, 2, 5]
        raw = np.load(caches[mode] / "matrix.npy")
        np.testing.assert_allclose(item["data"]["mean"], raw[:5].mean(axis=0))
        assert item["info"]["kept_count"] == 5 and item["info"]["sort_mode"] == "center"
    fast.save_plots(figures, work / "paired_figs")
    fast.save_profile_data(figures, caches, work / "paired_data", normalized=True)
    assert {p.name for p in (work / "paired_data").iterdir()} == {
        "center", "body", "view_keep.npy", "view_rows.tsv"}
    for mode in caches:
        assert (work / "paired_figs" / mode / "view_profile.pdf").exists()
        array = np.load(work / "paired_data" / mode / "view_heatmap_01.npy")
        assert array.shape[0] == 5 and np.all(array[3] == 0)
    # Signed cancellation is signal, not an all-zero row.
    values = np.load(caches["center"] / "matrix.npy", mmap_mode="r+")
    values[3, 0, :2] = [-1, 1]
    values.flush()
    assert fast.paired_rows(caches, regions={"view": mask}, filter_all_zero=True)["view"]["keep"].sum() == 5
    print("PASS: paired zero filtering, stable shared ordering, aggregation, profiles and exports")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--benchmark-bed")
    parser.add_argument("--benchmark-bw")
    parser.add_argument("--processors", type=int, default=1,
                        help="Also validate grouped cache with this worker count")
    args = parser.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    _, bed, bws = synthetic(work)
    # Exercise cache, profile summary and both rendering entry points.
    command = [sys.executable, str(Path(fast.__file__)), "compute", "-S", *map(str, bws),
               "-R", str(bed), "-o", str(work / "cache"), "-b", "100", "-a", "150"]
    subprocess.run(command, check=True)
    matrix, rows, keep, meta = fast.load_cache(work / "cache")
    means = fast.profile_summary(matrix, rows, meta["groups"], keep)
    expected = matrix[keep].mean(axis=0)
    np.testing.assert_allclose(means[meta["groups"][0]][0], expected, atol=1e-12)
    subprocess.run([sys.executable, str(Path(fast.__file__)), "plot", "--matrix", str(work / "cache"),
                    "-o", str(work / "figs"), "--max-rows", "10", "--sem"], check=True)
    # Multiple region groups, both samples and shared row identities.
    lines = bed.read_text().splitlines(True)
    bed1, bed2 = work / "group1.bed", work / "group2.bed"
    bed1.write_text("".join(lines[:14]))
    bed2.write_text("".join(lines[14:]))
    grouped = command[:command.index("-R")] + ["-R", str(bed1), str(bed2), "-o",
               str(work / "grouped"), "-b", "100", "-a", "150", "-p", str(args.processors)]
    subprocess.run(grouped, check=True)
    gm, gr, gk, gmeta = fast.load_cache(work / "grouped")
    np.testing.assert_allclose(matrix, gm, atol=1e-12)
    np.testing.assert_array_equal(keep, gk)
    grouped_means = fast.profile_summary(gm, gr, gmeta["groups"], gk)
    for group in gmeta["groups"]:
        mask = np.array([r[5] == group for r in gr]) & gk
        np.testing.assert_allclose(grouped_means[group][0], gm[mask].mean(axis=0), atol=1e-12)
    subprocess.run([sys.executable, str(Path(fast.__file__)), "plot", "--matrix", str(work / "grouped"),
                    "-o", str(work / "grouped_figs"), "--max-rows", "10", "--sem"], check=True)
    validate_plot_api(work, bed, bws, args.processors)
    if args.benchmark_bed and args.benchmark_bw:
        cfg = dict(mode="scale-regions", reference_point="TSS", upstream=1000, downstream=1000,
                   body_length=5000, bin_size=50, missing_as_zero=True, skip_zeros=True)
        compare(work, args.benchmark_bed, [args.benchmark_bw], cfg, "real_H3K9me3")
    assert not list(work.rglob("*.json"))
    print("PASS: numerical comparisons, filtering, cache, summary, profile and heatmap")


if __name__ == "__main__":
    main()
