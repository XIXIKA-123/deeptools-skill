#!/usr/bin/env python3
"""Reproduce README benchmarks and default-style examples using synthetic data."""
import argparse
import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'skills/deeptools/scripts'))
from fast_chip_profiles import compute_profiles, make_plots, save_plots


def inputs(work, count):
    rng = np.random.default_rng(20261106)
    length = 5_000_000
    x = np.arange(0, length, 50)
    phase = (x + 25) % 10_000 - 5_000
    paths = []
    for sample, amplitude, shift in [('A', 8., 0), ('B', 5., 250)]:
        signal = .2 + amplitude * np.exp(-.5 * ((phase - shift) / 600.) ** 2)
        signal *= rng.uniform(.8, 1.2, len(x))
        path = work / (sample + '.bw')
        with pyBigWig.open(str(path), 'w') as bw:
            bw.addHeader([('chrSynthetic', length)])
            bw.addEntries(['chrSynthetic'] * len(x), x.tolist(),
                          ends=(x + 50).tolist(), values=signal.tolist())
        paths.append(path)
    centers = rng.integers(1, 499, count) * 10_000 + 5_000 + rng.integers(-900, 901, count)
    widths = rng.integers(300, 2501, count)
    rows = [('chrSynthetic', int(c - w // 2), int(c + (w + 1) // 2),
             f'region_{i}', 0, '+' if i % 2 == 0 else '-')
            for i, (c, w) in enumerate(zip(centers, widths))]
    return paths, rows


def write_bed(path, rows):
    path.write_text(''.join('\t'.join(map(str, row)) + '\n' for row in rows))


def commands(bed, bws, work, workers):
    common = ['-R', str(bed), '-S', *map(str, bws), '-b', '3000', '-a', '3000',
              '-m', '5000', '-bs', '50', '--missingDataAsZero', '-p', str(workers)]
    return {
        'deepTools': [sys.executable, str(Path(sys.executable).parent / 'computeMatrix'),
                     'scale-regions', *common, '--sortRegions', 'keep', '--quiet',
                     '-o', str(work / 'deeptools.gz')],
        'deeptools-skill': [sys.executable, str(ROOT / 'skills/deeptools/scripts/fast_chip_profiles.py'),
                           'compute', '--mode', 'scale-regions', *common,
                           '--keep-zero-rows', '-o', str(work / 'skill')],
    }


def verify(work, rows):
    matrix = np.load(work / 'skill/matrix.npy', mmap_mode='r')
    expected_ids = {row[3]: i for i, row in enumerate(rows)}
    seen = set()
    max_error = 0.
    with gzip.open(work / 'deeptools.gz', 'rt') as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip().split('\t')
            row_id = expected_ids[fields[3]]
            actual = matrix[row_id].reshape(-1)
            expected = np.asarray(fields[6:], dtype=float)
            np.testing.assert_allclose(actual, expected, atol=5.1e-7, rtol=1e-7, equal_nan=True)
            max_error = max(max_error, float(np.nanmax(np.abs(actual - expected))))
            seen.add(row_id)
    assert len(seen) == len(rows)
    return max_error


def benchmark(work, output, sizes, repeats, workers):
    bws, rows = inputs(work, max(sizes))
    bed = work / 'regions.bed'
    results = []
    for size in sizes:
        selected = rows[:size]
        write_bed(bed, selected)
        cmds = commands(bed, bws, work, workers)
        # Warm both executables and input pages; do not flush operating-system caches.
        for tool in cmds:
            subprocess.run(cmds[tool], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        verify(work, selected)
        shutil.rmtree(work / 'skill')
        for repeat in range(1, repeats + 1):
            order = list(cmds) if repeat % 2 else list(reversed(cmds))
            for tool in order:
                start = time.perf_counter()
                subprocess.run(cmds[tool], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                seconds = time.perf_counter() - start
                results.append(dict(tool=tool, regions=size, repeat=repeat, seconds=seconds, workers=workers))
                print(f'{size:>6} regions | {tool:>15} | run {repeat} | {seconds:.3f} s', flush=True)
            error = verify(work, selected)
            print(f'  Matrix agreement: max absolute error {error:.3g}', flush=True)
            shutil.rmtree(work / 'skill')
        with (output / 'benchmark.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)
    return bws, rows, results


def speed_figure(output, records):
    sizes = sorted({r['regions'] for r in records})
    x = np.asarray(sizes) / 1000
    with plt.rc_context({'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6,
                         'ytick.labelsize': 6, 'legend.fontsize': 6,
                         'axes.linewidth': .5, 'pdf.fonttype': 42}):
        fig, ax = plt.subplots(figsize=(89 / 25.4, 60 / 25.4))
        fig.subplots_adjust(left=.16, right=.96, bottom=.20, top=.95)
        for tool, color, marker in [('deepTools', '#555555', 'o'),
                                    ('deeptools-skill', '#3B6FB6', 's')]:
            values = [np.asarray([r['seconds'] for r in records if r['tool'] == tool and r['regions'] == n]) for n in sizes]
            median = np.asarray([np.median(v) for v in values])
            lo = np.asarray([v.min() for v in values])
            hi = np.asarray([v.max() for v in values])
            ax.fill_between(x, lo, hi, color=color, alpha=.15, linewidth=0)
            ax.plot(x, median, marker=marker, color=color, linewidth=.8,
                    markersize=2.5, markeredgewidth=.4, label=tool, clip_on=False)
            for position, seconds in zip(x, median):
                first_skill = tool == 'deeptools-skill' and position == x[0]
                offset = (8, 2) if first_skill else (0, 6 if tool == 'deepTools' else 4)
                ax.annotate(f'{seconds:.2f}', (position, seconds), xytext=offset,
                            textcoords='offset points', ha='left' if first_skill else 'center',
                            va='bottom', fontsize=6, color='black')
        ax.set(xlabel='BED regions (×10³)', ylabel='Time (s)', ylim=(0, None),
               xticks=x, xticklabels=[f'{n / 1000:g}' for n in sizes])
        for side in ['top', 'right']: ax.spines[side].set_visible(False)
        ax.tick_params(direction='out', length=2.5, width=.5)
        ax.legend(loc='center right',
                  frameon=False, handlelength=1.5, handletextpad=.5, borderaxespad=0)
        for ext in ['png', 'pdf']: fig.savefig(output / f'benchmark.{ext}', dpi=600)
        plt.close(fig)


def style_figures(work, output, bws, rows):
    bed = work / 'demo.bed'
    write_bed(bed, rows[:500])
    samples = [dict(path=str(p), sample=p.stem, group=p.stem, label='Sample ' + p.stem) for p in bws]
    caches = compute_profiles(samples, {'Synthetic peaks': str(bed)}, work / 'demo-cache', kind='peak', workers=1)
    plots = make_plots(caches, samples=samples, style={'signal_label': 'Signal (a.u.)'})
    for views in plots.values():
        for figures in views.values():
            figures['profile'].axes[0].set_title('')
            figures['heatmap'].suptitle('')
    save_plots(plots, output / 'style')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'assets')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='deeptools-readme-') as temp:
        work = Path(temp)
        bws, rows, records = benchmark(work, args.output, [500, 2000, 5000, 10000], 3, 1)
        speed_figure(args.output, records)
        style_figures(work, args.output, bws, rows)
    print('README assets complete.', flush=True)


if __name__ == '__main__':
    main()
