#!/usr/bin/env python3
"""Plot everything bam2prof measures, in one figure by default:

  A  damage: the substitution profiles at both fragment ends (all 12 types, C>T and G>A highlighted)
  B  DNA composition around the breaks (needs a -comp run; with -fa it includes the reference flank)
  C  fragment length distribution (needs --isize: the file given to bam2prof -is)

    python3 plot_bam2prof.py OUT_DIR [OUTPUT_PREFIX] [--isize FILE] [--title TITLE]

writes OUTPUT_PREFIX_summary.png and .pdf (default prefix: OUT_DIR/bam2prof). Panels with no data are left out;
--only damage|composition|isize draws just one of them, as OUTPUT_PREFIX_<panel>.png/.pdf. OUT_DIR is the directory
given to bam2prof -o and must hold one profile (one -classic run, or one reference of a -meta run).
"""
import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

RED, BLUE, ORANGE, GREY = '#C0392B', '#2F6DB5', '#E8803A', '#8C8C88'
# colours for the 10 substitution types other than C>T / G>A
CYCLE = ['#2E9E5B', '#E8A33D', '#7A5BA6', '#17A2B8', '#D2691E', '#8D6E63', '#E377C2', '#7F7F7F', '#BCBD22', '#1B7F79']
BASES = ['A', 'C', 'G', 'T']
BASE_COLORS = {'A': '#2E9E5B', 'C': '#3E7CB1', 'G': '#E8A33D', 'T': '#C0392B'}
PAIRED_BLUE, MERGED_ORANGE, TOTAL_INK = '#2F6DB5', '#E8803A', '#222222'


def set_style():
    plt.rcParams.update({
        'font.size': 12,
        'axes.edgecolor': '#888888',
        'axes.labelcolor': '#222222',
        'text.color': '#222222',
        'xtick.color': '#444444',
        'ytick.color': '#444444',
        'pdf.fonttype': 42,  # editable text in the PDFs
    })


def clean_axes(ax):
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------- damage (substitution) profiles

def load_damage(path):
    return pd.read_csv(path, delimiter='\t')


def damage_ylim(fivep, threep, ylim=None):
    return tuple(ylim) if ylim else (0, max(fivep.max().max(), threep.max().max()) * 1.12)


def _colour_for(column, others, grey_others):
    if column == 'C>T':
        return RED
    if column == 'G>A':
        return BLUE
    return GREY if grey_others else CYCLE[others.index(column) % len(CYCLE)]


def draw_damage(ax, data, end, ylim, xlim=None, annotate=False, percent=False, grey_others=False):
    """One end of the substitution profile. x is the distance from the fragment end (1 = terminal base); the 3'
    axis runs right to left so that both fragment ends sit on the outside of the figure, as in the fragment."""
    x = np.arange(1, len(data) + 1)
    others = [c for c in data.columns if c not in ('C>T', 'G>A')]
    for c in others:  # the 10 other substitution types underneath
        ax.plot(x, data[c], color=_colour_for(c, others, grey_others), linewidth=1.1 if grey_others else 1.5,
                alpha=0.6 if grey_others else 0.9, zorder=2, solid_capstyle='round')
    for c, col in (('G>A', BLUE), ('C>T', RED)):
        if c in data.columns:
            ax.plot(x, data[c], color=col, linewidth=2.6, marker='o', markersize=4.5, zorder=3)
            if annotate:
                # value at the terminal base, labelled directly on the plot; the G>A label sits further in and up
                # (with a leader line) so it clears the steeper C>T curve
                dx, dy = (10, 4) if c == 'C>T' else (62, 26)
                value = f'{100 * data[c].iloc[0]:.1f}%' if percent else f'{data[c].iloc[0]:.3f}'
                ax.annotate(f'{c.replace(">", "→")}  {value}', xy=(1, data[c].iloc[0]),
                            xytext=(dx if end == '5' else -dx, dy), textcoords='offset points',
                            ha='left' if end == '5' else 'right', va='bottom', color=col, fontsize=11, fontweight='bold',
                            arrowprops=dict(arrowstyle='-', color=col, lw=0.8, shrinkA=2, shrinkB=4) if c == 'G>A' else None,
                            path_effects=[pe.withStroke(linewidth=3, foreground='white')], zorder=4)
    ax.set_xlim(xlim if xlim else (0.5, len(data) + 0.5))
    if end == '3':
        ax.invert_xaxis()
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(2 if len(data) > 12 else 1))
    if percent:
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=1 if ylim[1] < 0.1 else 0))
    ax.set_xlabel(f"Distance from the {end}' end (bp)")
    ax.grid(axis='y', color='#EDEDED', linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    clean_axes(ax)
    ax.set_title(f"{end}' end", fontsize=13, fontweight='bold', pad=10)


def damage_legend_handles(data, grey_others=False):
    others = [c for c in data.columns if c not in ('C>T', 'G>A')]
    h = [Line2D([], [], color=RED, lw=2.6, marker='o', ms=4.5, label='C→T'),
         Line2D([], [], color=BLUE, lw=2.6, marker='o', ms=4.5, label='G→A')]
    if grey_others:
        h.append(Line2D([], [], color=GREY, lw=1.1, alpha=0.6, label=f'other substitutions ({len(others)} types)'))
    else:
        h += [Line2D([], [], color=_colour_for(c, others, False), lw=1.5, alpha=0.9, label=c.replace('>', '→')) for c in others]
    return h


# ---------------------------------------------------------------- base composition around the fragment ends

def load_comp(path):
    data = pd.read_csv(path, delimiter='\t')
    counts = data[BASES]
    totals = counts.sum(axis=1)
    freq = counts.div(totals.replace(0, np.nan), axis=0)
    data = data.assign(**{b: freq[b] for b in BASES})
    return data.sort_values('pos').reset_index(drop=True)


def draw_comp(ax, data, end_label, inside_is_negative):
    """end_label: '5' or '3'. inside_is_negative: True for the 3' panel, whose in-fragment
    positions are negative (distance back from the last base) while any reference flank is positive."""
    pos = data['pos'].to_numpy()
    has_flank = (inside_is_negative and (pos > 0).any()) or (not inside_is_negative and (pos < 0).any())

    if has_flank:
        boundary = -0.5 if not inside_is_negative else 0.0
        xmin, xmax = pos.min() - 0.5, pos.max() + 0.5
        if inside_is_negative:
            ax.axvspan(0.5, xmax, color='#F2F2F2', zorder=0)
            ax.text(0.5 + (xmax - 0.5) / 2, 1.02, 'reference flank', transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=9, color='#999999')
        else:
            ax.axvspan(xmin, -0.5, color='#F2F2F2', zorder=0)
            ax.text(xmin + (-0.5 - xmin) / 2, 1.02, 'reference flank', transform=ax.get_xaxis_transform(),
                    ha='center', va='bottom', fontsize=9, color='#999999')
        ax.axvline(boundary, color='#999999', linestyle=':', linewidth=1.2, zorder=1)

    ax.axhline(0.25, color='#CCCCCC', linestyle='--', linewidth=1, zorder=1)

    for b in BASES:
        ax.plot(pos, data[b], label=b, color=BASE_COLORS[b], linewidth=2.2, marker='o', markersize=3.5, zorder=3)

    ax.set_xlabel(f"Position relative to the {end_label}' end (bp)")
    ax.set_ylabel('Base frequency')
    ax.set_ylim(0, max(0.55, np.nanmax(data[BASES].to_numpy()) * 1.15))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    ax.grid(axis='y', color='#EDEDED', linewidth=0.8, zorder=0)
    clean_axes(ax)
    ax.set_title(f"{end_label}' end", fontsize=13, fontweight='bold', pad=14)


def comp_legend_handles():
    return [Line2D([], [], color=BASE_COLORS[b], lw=2.2, marker='o', ms=3.5, label=b) for b in BASES]


# ---------------------------------------------------------------- fragment length distribution

def _load_counts(path):
    """-> (lengths, counts) from 'count<TAB>length' lines, or None if the file is missing / empty"""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    a = np.loadtxt(path, dtype=np.int64, ndmin=2)
    return a[:, 1], a[:, 0]


def _median(lengths, counts):
    order = np.argsort(lengths)
    cum = np.cumsum(counts[order])
    return lengths[order][np.searchsorted(cum, cum[-1] / 2)]


def _percentile(lengths, counts, q):
    order = np.argsort(lengths)
    cum = np.cumsum(counts[order])
    return lengths[order][np.searchsorted(cum, cum[-1] * q)]


def load_isize_series(isize_file):
    """Fragment length classes written by bam2prof -is: [(label, colour, (lengths, counts)), ...]. Properly paired
    fragments (<file>.properly_paired, or all pairs from <file>.paired with -is-allpaired) and merged / single-end
    molecules (<file>.merged); falls back to the combined <file> as a single series."""
    series = []
    paired, paired_label = _load_counts(isize_file + '.properly_paired'), 'properly paired'
    if paired is None:
        paired, paired_label = _load_counts(isize_file + '.paired'), 'paired'
    merged = _load_counts(isize_file + '.merged')
    if paired is not None:
        series.append((paired_label, PAIRED_BLUE, paired))
    if merged is not None:
        series.append(('merged / single-end', MERGED_ORANGE, merged))
    if not series:
        combined = _load_counts(isize_file)
        if combined is None:
            raise SystemExit(f'No fragment lengths found in {isize_file}(.properly_paired/.paired/.merged)')
        series.append(('all fragments', GREY, combined))
    return series


def draw_isize(ax, series, xlim=None, fraction=False, log=False, legend_fontsize=11):
    if xlim:
        xmin, xmax = xlim
    else:
        # from the lowest 0.1th to the highest 99.9th percentile of any series, so a small class is not clipped by a big one
        xmin = max(0, min(_percentile(l, c, 0.001) for _, _, (l, c) in series) - 5)
        xmax = max(_percentile(l, c, 0.999) for _, _, (l, c) in series) + 5

    for name, colour, (lengths, counts) in series:
        order = np.argsort(lengths)
        x, y = lengths[order], counts[order].astype(float)
        n = counts.sum()
        if fraction:
            y = y / n
        # a continuous curve over every length (lengths with no fragments are 0, not missing)
        grid = np.arange(int(x.min()), int(x.max()) + 1)
        yy = np.zeros(len(grid)); yy[x - grid[0]] = y
        if log:
            yy[yy == 0] = np.nan  # no fragments of this length: a gap, not a drop to zero, on a log axis
        label = f'{name}  (n = {n:,}, median {_median(lengths, counts)} bp)'
        ax.fill_between(grid, yy, step='mid', color=colour, alpha=0.22, linewidth=0)
        ax.step(grid, yy, where='mid', color=colour, linewidth=1.8, label=label)

    if len(series) > 1 and not fraction:
        # all molecules together: the two classes are one continuous distribution (merged below the merging limit, pairs above)
        tot = {}
        for _, _, (lengths, counts) in series:
            for l, c in zip(lengths, counts):
                tot[l] = tot.get(l, 0) + c
        tl = np.arange(min(tot), max(tot) + 1)
        ty = np.array([tot.get(l, 0) for l in tl], dtype=float)
        n_all = int(ty.sum())
        if log:
            ty[ty == 0] = np.nan
        ax.step(tl, ty, where='mid', color=TOTAL_INK, linewidth=1.0, linestyle='--', label=f'all fragments  (n = {n_all:,})', zorder=4)

    ax.set_xlim(xmin, xmax)
    if log:
        ax.set_yscale('log')
    else:
        ax.set_ylim(0, ax.get_ylim()[1] * 1.3)  # headroom for the legend
        if not fraction:
            ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
            if max(c.max() for *_, (_, c) in series) >= 1000:
                ax.yaxis.set_major_formatter(mticker.EngFormatter())
    ax.set_xlabel('Fragment length (bp)')
    ax.set_ylabel('Fraction of fragments in class' if fraction else 'Number of fragments')
    ax.grid(axis='y', color='#EDEDED', linewidth=0.8)
    ax.set_axisbelow(True)
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=legend_fontsize, loc='upper right')


# ---------------------------------------------------------------- command line

def _find(directory, suffix):
    """the one file ending in `suffix` in `directory`, or None; several (e.g. -meta runs) is an error"""
    hits = sorted(glob.glob(os.path.join(directory, '*' + suffix)))
    if len(hits) > 1:
        raise SystemExit(f'Several files match *{suffix} in {directory}: {", ".join(map(os.path.basename, hits))}\n'
                         'Put the profile you want to plot in its own directory, or name the files with --prof5/--prof3/--comp5/--comp3.')
    return hits[0] if hits else None


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('out_dir', help='The directory given to bam2prof -o (holds *_5p.prof, *_3p.prof and, with -comp, *_5p_comp.prof / *_3p_comp.prof)')
    p.add_argument('output_prefix', nargs='?', default=None, help='Prefix for the output files (default: <out_dir>/bam2prof)')
    p.add_argument('--isize', default=None, metavar='FILE', help='The file given to bam2prof -is; adds the fragment length panel '
                   '(properly paired and merged / single-end molecules, read from FILE.properly_paired / FILE.paired and FILE.merged)')
    p.add_argument('--only', choices=['damage', 'composition', 'isize'], default=None, help='Draw only this panel instead of the summary figure')
    p.add_argument('--title', default=None, help='Title of the figure')
    p.add_argument('--prof5', metavar='FILE', help="5' substitution profile, if not the one found in out_dir")
    p.add_argument('--prof3', metavar='FILE', help="3' substitution profile, if not the one found in out_dir")
    p.add_argument('--comp5', metavar='FILE', help="5' composition profile, if not the one found in out_dir")
    p.add_argument('--comp3', metavar='FILE', help="3' composition profile, if not the one found in out_dir")
    d = p.add_argument_group('damage panel')
    d.add_argument('--xlim', type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=None, help='Fixed x range in bp from the fragment end (1 = terminal base)')
    d.add_argument('--ylim', type=float, nargs=2, metavar=('YMIN', 'YMAX'), default=None, help='Fixed y range (frequencies, e.g. 0 0.05)')
    d.add_argument('--annotate', action='store_true', help='Write the terminal-base C>T and G>A values on the plot')
    d.add_argument('--percent', action='store_true', help='Show the y axis (and --annotate values) as percentages instead of frequencies')
    d.add_argument('--grey-others', action='store_true', help='Draw the 10 substitution types other than C>T and G>A as quiet grey lines')
    i = p.add_argument_group('fragment length panel')
    i.add_argument('--isize-xlim', type=float, nargs=2, metavar=('XMIN', 'XMAX'), default=None, help='Fixed x range in bp (default: 0.1th to 99.9th percentile)')
    i.add_argument('--isize-log', action='store_true', help='Log-scale y axis; keeps a small class (e.g. properly paired fragments) visible next to a much larger one')
    i.add_argument('--isize-fraction', action='store_true', help='Fraction of each class (each sums to 1) instead of counts, to compare shapes')
    return p.parse_args()


def make_figure(args, damage, comp, series):
    """damage: (5' data, 3' data) or None; comp: (5', 3') or None; series: fragment length series or None"""
    rows = [k for k, v in (('damage', damage), ('composition', comp)) if v is not None]
    ncols = (2 if rows else 0) + (1 if series else 0)
    width_ratios = ([1, 1] if rows else []) + ([1.2] if series else [])
    nrows = max(len(rows), 1)
    fig = plt.figure(figsize=(max(sum(width_ratios) * 6.2, 9), 5.0 * nrows + (0.8 if nrows == 1 else 1.4)))
    fig.patch.set_facecolor('white')
    gs = fig.add_gridspec(nrows, ncols, width_ratios=width_ratios, hspace=0.62, wspace=0.14 if ncols == 2 else 0.2)

    axes_by_row, legends, legend_ncol = [], [], []
    for r, kind in enumerate(rows):
        if kind == 'damage':
            ylim = damage_ylim(damage[0], damage[1], args.ylim)
            a5 = fig.add_subplot(gs[r, 0]); a3 = fig.add_subplot(gs[r, 1], sharey=a5)
            opts = dict(xlim=args.xlim, annotate=args.annotate, percent=args.percent, grey_others=args.grey_others)
            draw_damage(a5, damage[0], '5', ylim, **opts)
            draw_damage(a3, damage[1], '3', ylim, **opts)
            a5.set_ylabel('Substitution frequency')
            a3.tick_params(labelleft=False)
            legends.append(damage_legend_handles(damage[0], args.grey_others))
            legend_ncol.append(3 if args.grey_others else 6)
        else:
            a5 = fig.add_subplot(gs[r, 0]); a3 = fig.add_subplot(gs[r, 1])
            draw_comp(a5, comp[0], '5', inside_is_negative=False)
            draw_comp(a3, comp[1], '3', inside_is_negative=True)
            legends.append(comp_legend_handles())
            legend_ncol.append(4)
        axes_by_row.append((a5, a3))

    ix = None
    if series:
        ix = fig.add_subplot(gs[:, ncols - 1])
        draw_isize(ix, series, xlim=args.isize_xlim, fraction=args.isize_fraction, log=args.isize_log,
                   legend_fontsize=9.5 if rows else 11)
        if rows:
            ix.set_title('Fragment length', fontsize=13, fontweight='bold', pad=10)

    fig.subplots_adjust(left=0.06 if ncols == 3 else 0.07, right=0.985,
                        bottom=0.08 if nrows == 2 else 0.14,
                        top=0.84 if nrows == 2 else ((0.74 if rows[0] == 'damage' else 0.78) if rows else 0.9))
    for (left_ax, right_ax), handles, nc in zip(axes_by_row, legends, legend_ncol):
        l, r, top = left_ax.get_position().x0, right_ax.get_position().x1, left_ax.get_position().y1
        fig.legend(handles=handles, loc='lower center', ncol=nc, bbox_to_anchor=((l + r) / 2, top + 0.045), frameon=False, fontsize=11.5)
    blocks = [row[0] for row in axes_by_row] + ([ix] if ix is not None else [])
    if len(blocks) > 1:  # panel letters
        for k, ax in enumerate(blocks):
            pos = ax.get_position()
            x = 0.008 if k < len(axes_by_row) else pos.x0 - 0.055
            fig.text(x, pos.y1 + 0.075, 'ABC'[k], fontsize=18, fontweight='bold', va='bottom')
    if args.title:
        fig.suptitle(args.title, fontsize=16, y=0.985)
    return fig


def main():
    args = parse_args()
    set_style()
    want = lambda panel: args.only in (None, panel)

    damage = comp = series = None
    if want('damage'):
        f5 = args.prof5 or _find(args.out_dir, '_5p.prof')
        f3 = args.prof3 or _find(args.out_dir, '_3p.prof')
        if not (f5 and f3):
            raise SystemExit(f'No *_5p.prof / *_3p.prof files found in {args.out_dir}')
        damage = (load_damage(f5), load_damage(f3))
    if want('composition'):
        c5 = args.comp5 or _find(args.out_dir, '_5p_comp.prof')
        c3 = args.comp3 or _find(args.out_dir, '_3p_comp.prof')
        if c5 and c3:
            comp = (load_comp(c5), load_comp(c3))
        elif args.only == 'composition':
            raise SystemExit(f'No *_5p_comp.prof / *_3p_comp.prof files found in {args.out_dir} (run bam2prof with -comp)')
    if want('isize'):
        if args.isize:
            series = load_isize_series(args.isize)
        elif args.only == 'isize':
            raise SystemExit('--only isize needs --isize FILE (the file given to bam2prof -is)')

    prefix = args.output_prefix or os.path.join(args.out_dir, 'bam2prof')
    name = 'summary' if args.only is None else args.only
    fig = make_figure(args, damage, comp, series)
    fig.savefig(f'{prefix}_{name}.pdf')
    fig.savefig(f'{prefix}_{name}.png', dpi=200)
    print(f'Wrote {prefix}_{name}.pdf')
    print(f'Wrote {prefix}_{name}.png')


if __name__ == '__main__':
    main()
