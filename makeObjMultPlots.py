import argparse
import uproot
import numpy as np
import boost_histogram as bh
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.patches as mpatches
import mplhep as hep
import os

hep.style.use('CMS')

FONTSIZE = 8

triggers = [
    "L1_ZeroBias",
    "L1_DoubleEG_LooseIso20_LooseIso12_er1p5",
    "L1_DoubleMu_15_7",
    "L1_HTT280er",
    "L1_AXO_Nominal",
    "L1_CICADA_Medium",
]

TRIGGER_LABELS = {
    "L1_ZeroBias": "Zero Bias",
    "L1_DoubleEG_LooseIso20_LooseIso12_er1p5": "Double $e\gamma$",
    "L1_DoubleMu_15_7": "Double Muon",
    "L1_HTT280er": "Jet $H_T$",
    "L1_AXO_Nominal": "AXOL1TL",
    "L1_CICADA_Medium": "CICADA",
}

TRIGGER_COLORS = {
    'L1_ZeroBias' : '#1845fb',
    'L1_HTT280er': '#ff5e02',
    'L1_CICADA_Medium': '#c91f16',
    'L1_DoubleMu_15_7': '#578dff',
    'L1_DoubleEG_LooseIso20_LooseIso12_er1p5': '#adad7d',
    'L1_AXO_Nominal': '#86c8dd'
}

TRIGGER_LINEWIDTH = {
    'L1_ZeroBias' : 1,
    'L1_HTT280er': 1,
    'L1_CICADA_Medium': 3,
    'L1_DoubleMu_15_7': 1,
    'L1_DoubleEG_LooseIso20_LooseIso12_er1p5': 1,
    'L1_AXO_Nominal': 3
}

# Default axis limits and x-labels per object type
OBJ_DEFAULTS = {
    "L1Mu": {
    "hist_key": "L1Mu_mult",
    "x_label": r"L1 Muon Multiplicity",
    "x_min": -0.5, "x_max": 8.5,
    "y_min": 5e-10, "y_max": 5e6,
    "y_min_ratio": 1e-1, "y_max_ratio": 1e4,
    },
    "L1EG": {
    "hist_key": "L1EG_mult",
    "x_label": r"L1 EG Multiplicity",
    "x_min": -0.5, "x_max": 12.5,
    "y_min": 5e-5, "y_max": 5e1,
    "y_min_ratio": 1e-3, "y_max_ratio": 1e2,
    },
    "L1Jet": {
    "hist_key": "L1Jet_mult",
    "x_label": "L1 Jet multiplicity",
    "x_min": -0.5, "x_max": 10.5,
    "y_min": 5e-8, "y_max": 5e4,
    "y_min_ratio": 1e-5, "y_max_ratio": 1e5,
    },
    "L1HT": {
    "hist_key": "l1_ht",
    "x_label": r"L1 $H_T$ [GeV]",
    "x_min": 0, "x_max": 1000,
    "y_min": 5e-10, "y_max": 5e2,
    "y_min_ratio": 1e-3, "y_max_ratio": 1e5,
    },
    "L1MET": {
    "hist_key": "l1_met",
    "x_label": r"L1 $p_T^{\text{miss}}$ [GeV]",
    "x_min": 0, "x_max": 180,
    "y_min": 5e-10, "y_max": 5e2,
    "y_min_ratio": 1e-2, "y_max_ratio": 1e5,
    },
}


def load_root_hists(root_file, hist_key, triggers):
    """
    Load histograms from ROOT file into a dict keyed by trigger name.
    Returns dict: {trigger: (counts, bins)}
    """
    hists = {}
    with uproot.open(root_file) as f:
        for trigger in triggers:
            key = f"{trigger}_{hist_key}"
            if key in f:
                h = f[key]
                counts, bins = h.to_numpy()
                hists[trigger] = (counts, bins)
            else:
                print(f"  WARNING: key '{key}' not found in ROOT file, skipping.")
    return hists


def draw_hist1d(counts, bins, ax=None, label="",
        norm=False, linestyle='solid', color=None, linewidth=1):

    if norm:
        norm_factor = np.sum(counts) * np.diff(bins)
    else:
        norm_factor = 1

    _counts = counts / norm_factor if norm else counts
    errs = np.sqrt(counts) / norm_factor if norm else np.sqrt(counts)
    _errs = np.where(_counts == 0, 0, errs)

    bin_centres = 0.5 * (bins[1:] + bins[:-1])

    if color is not None:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="", color=color, linewidth=linewidth)
    else:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="", linewidth=linewidth)
    color = l[0].get_color()
    ax.errorbar(
        x=bins, y=np.append(_counts, _counts[-1]), drawstyle="steps-post", label=label,
        color=color, linestyle=linestyle, linewidth=linewidth
    )
    return l


def draw_ratio(counts_num, bins_num, counts_denom, bins_denom, ax=None, color=None, linewidth=1,
           label="", norm=False):

    norm_factor_denom = np.sum(counts_denom) * np.diff(bins_denom) if norm else 1
    counts_denom = counts_denom / norm_factor_denom if norm else counts_denom
    errs_denom = np.sqrt(counts_denom * (1 - counts_denom / norm_factor_denom)) / norm_factor_denom if norm else np.sqrt(counts_denom)

    norm_factor_num = np.sum(counts_num) * np.diff(bins_num) if norm else 1
    counts_num = counts_num / norm_factor_num if norm else counts_num
    errs_num = np.sqrt(counts_num * (1 - counts_num / norm_factor_num)) / norm_factor_num if norm else np.sqrt(counts_num)

    denom = np.where(counts_denom == 0, np.nan, counts_denom)
    ratio = counts_num / denom
    x = 0.5 * (bins_num[:-1] + bins_num[1:])

    error = ratio * np.sqrt((errs_num / np.where(counts_num == 0, np.nan, counts_num))**2 +
                (errs_denom / np.where(counts_denom == 0, np.nan, counts_denom))**2)

    if color is not None:
        l = ax.errorbar(x=x, y=ratio, yerr=error, linestyle="", color=color)
    else:
        l = ax.errorbar(x=x, y=ratio, yerr=error, linestyle="")
    color = l[0].get_color()
    ax.errorbar(
        x=bins_num, y=np.append(ratio, ratio[-1]), drawstyle="steps-post", label=label,
        color=color, linestyle='solid', linewidth=linewidth,
    )
    return l


def make_plot(hists, triggers, x_label,
          x_min, x_max, y_min, y_max, y_min_ratio, y_max_ratio, output,
          log_scale=True, norm=False, leg_loc='upper right'):

    fig, ax = plt.subplots(2, figsize=(2.1, 2.1), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.12)
    ax[1].plot(np.linspace(x_min, x_max, 10), np.ones(10), '--', color='darkgray')

    counts_denom, bins_denom = hists["L1_ZeroBias"]

    # Compute all ratio values to determine y range for ratio panel
    all_ratios = []
    for trigger in triggers:
        if trigger == "L1_ZeroBias" or trigger not in hists:
            continue
        counts, bins = hists[trigger]
        denom = np.where(counts_denom == 0, np.nan, counts_denom)
        ratio = counts / denom
        finite = ratio[np.isfinite(ratio) & (ratio > 0)]
        if len(finite):
            all_ratios.extend(finite)

    for trigger in triggers:
        if trigger not in hists:
            continue
        counts, bins = hists[trigger]
        color = TRIGGER_COLORS[trigger]
        trigger_label = TRIGGER_LABELS[trigger]

        l = draw_hist1d(
            counts, bins,
            ax=ax[0],
            label=trigger_label,
            norm=norm,
            color=color,
            linewidth=TRIGGER_LINEWIDTH[trigger],
        )
        if trigger != "L1_ZeroBias":
            color = l[0].get_color()
            draw_ratio(
                counts, bins,
                counts_denom, bins_denom,
                ax=ax[1],
                label=trigger_label,
                norm=norm,
                color=color,
                linewidth=TRIGGER_LINEWIDTH[trigger]
            )

    ax[0].set_xlim([x_min, x_max])
    ax[0].set_ylim([y_min, y_max])
    ax[1].set_ylim([y_min_ratio, y_max_ratio])

    if log_scale:
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
    ax[0].set_ylabel(f"Events{' [A.U.]' if norm else ''}", loc="top", fontsize=FONTSIZE, labelpad=0)
    ax[1].set_ylabel("Ratio to Zero Bias", loc="top", fontsize=FONTSIZE, labelpad=0)
    ax[1].set_xlabel(x_label, loc="right", fontsize=FONTSIZE, labelpad=2)
    ax[1].yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=[1.0], numticks=10))
    ax[1].yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    exp_min = int(np.floor(np.log10(y_min_ratio)))
    exp_max = int(np.ceil(np.log10(y_max_ratio)))
    labeled_exps = set(range(exp_min, exp_max + 1, 2))
    ax[1].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda x, pos: f'$10^{{{int(round(np.log10(x)))}}}$' if int(round(np.log10(x))) in labeled_exps else ''
    ))
    ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    legend_handles = [
        mpatches.Rectangle(
            (0, 0), 1, 1,
            fill=False,
            edgecolor=TRIGGER_COLORS[t],
            linewidth=TRIGGER_LINEWIDTH[t],
            label=TRIGGER_LABELS[t],
        )
        for t in triggers if t in hists
    ]
    ax[0].legend(handles=legend_handles, loc=leg_loc, frameon=False, fontsize=FONTSIZE, ncols=2, columnspacing=0.7)

    hep.cms.label(
        "Preliminary",
        data=True,
        rlabel="2024 (13.6 TeV)",
        fontsize=FONTSIZE,
        ax=ax[0],
    )

    ax[0].tick_params(axis='both', which='major', labelsize=FONTSIZE, length=4, pad=2)
    ax[1].tick_params(axis='both', which='major', labelsize=FONTSIZE, length=4, pad=2)
    ax[0].minorticks_off()
    for ax_ in ax:
        for spine in ax_.spines.values():
            spine.set_linewidth(0.8)

    fig.subplots_adjust(left=0.18, right=0.98, top=0.93, bottom=0.12, hspace=0)

    out_dir = os.path.dirname(output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(f"{output}.pdf", format="pdf")
    fig.savefig(f"{output}.png", format="png")
    print(f"Saved {output}.pdf and {output}.png")


def main(args):

    defaults = OBJ_DEFAULTS[args.object]

    hists = load_root_hists(args.input, defaults["hist_key"], triggers)

    x_min = args.x_min if args.x_min is not None else defaults["x_min"]
    x_max = args.x_max if args.x_max is not None else defaults["x_max"]
    y_min = args.y_min if args.y_min is not None else defaults["y_min"]
    y_max = args.y_max if args.y_max is not None else defaults["y_max"]
    y_min_ratio = args.y_min_ratio if args.y_min_ratio is not None else defaults["y_min_ratio"]
    y_max_ratio = args.y_max_ratio if args.y_max_ratio is not None else defaults["y_max_ratio"]

    make_plot(
        hists,
        triggers,
        defaults["x_label"],
        x_min, x_max, y_min, y_max, y_min_ratio, y_max_ratio,
        args.output,
        log_scale=True,
        norm=True,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Make a multiplicity plot for a single L1 object type (L1Mu, L1EG, or L1Jet)."
    )
    parser.add_argument(
        "--object",
        required=True,
        choices=["L1Mu", "L1EG", "L1Jet", "L1HT", "L1MET"],
        help="Which object type to plot"
    )
    parser.add_argument(
        "--input",
        default="histograms/hist.root",
        help="Input .root histogram file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Full output path prefix, e.g. plots/mult_L1Mu (extensions .pdf/.png added automatically)"
    )
    parser.add_argument("--x-min", type=float, default=None, help="x-axis minimum")
    parser.add_argument("--x-max", type=float, default=None, help="x-axis maximum")
    parser.add_argument("--y-min", type=float, default=None, help="y-axis minimum")
    parser.add_argument("--y-max", type=float, default=None, help="y-axis maximum")
    parser.add_argument("--y-min-ratio", type=float, default=None, help="ratio y-axis minimum")
    parser.add_argument("--y-max-ratio", type=float, default=None, help="ratio y-axis maximum")

    args = parser.parse_args()
    main(args)
