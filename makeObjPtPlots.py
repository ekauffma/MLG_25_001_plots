import argparse
import uproot
import numpy as np
import boost_histogram as bh
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mplhep as hep
import os

hep.style.use('CMS')

FONTSIZE = 8

triggers = [
    "DST_PFScouting_ZeroBias",
    "DST_PFScouting_JetHT",
    "DST_PFScouting_CICADAMedium",
    "DST_PFScouting_DoubleMuon",
    "DST_PFScouting_DoubleEG",
    "DST_PFScouting_AXONominal",
]

TRIGGER_LABELS = {
    'DST_PFScouting_ZeroBias' : 'Zero Bias',
    'DST_PFScouting_JetHT': 'Jet $H_T$',
    'DST_PFScouting_CICADAMedium': 'CICADA',
    'DST_PFScouting_DoubleMuon': 'Double Muon',
    'DST_PFScouting_DoubleEG': 'Double $e\gamma$',
    'DST_PFScouting_AXONominal': 'AXO'
}

TRIGGER_COLORS = {
    'DST_PFScouting_ZeroBias' : '#1845fb',
    'DST_PFScouting_JetHT': '#ff5e02',
    'DST_PFScouting_CICADAMedium': '#c91f16',
    'DST_PFScouting_DoubleMuon': '#578dff',
    'DST_PFScouting_DoubleEG': '#adad7d',
    'DST_PFScouting_AXONominal': '#86c8dd'
}

# Default axis limits and x-labels per object type
OBJ_DEFAULTS = {
    "L1Mu": {
        "hist_key": "L1Mu_0_pt",
        "x_label": r"L1 Muon $p_T$ [GeV]",
        "x_min": 0.0, "x_max": 100,
        "y_min": 5e-7, "y_max": 5e5,
    },
    "ScoutingMuonVtx": {
        "hist_key": "ScoutingMuonVtx_0_pt",
        "x_label": r"Scouting Muon $p_T$ [GeV]",
        "x_min": 0.0, "x_max": 100,
        "y_min": 5e-7, "y_max": 5e5,
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
                norm=False, linestyle='solid', color=None):

    if norm:
        norm_factor = np.sum(counts) * np.diff(bins)
    else:
        norm_factor = 1

    _counts = counts / norm_factor if norm else counts
    errs = np.sqrt(counts) / norm_factor if norm else np.sqrt(counts)
    _errs = np.where(_counts == 0, 0, errs)

    bin_centres = 0.5 * (bins[1:] + bins[:-1])

    if color is not None:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="", color=color)
    else:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="")
    color = l[0].get_color()
    ax.errorbar(
        x=bins, y=np.append(_counts, _counts[-1]), drawstyle="steps-post", label=label,
        color=color, linestyle=linestyle, linewidth=1,
    )
    return l


def draw_ratio(counts_num, bins_num, counts_denom, bins_denom, ax=None, color=None,
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
        color=color, linestyle='solid', linewidth=1,
    )
    return l


def make_plot(hists, triggers, x_label,
              x_min, x_max, y_min, y_max, output,
              log_scale=True, norm=False, leg_loc='upper right'):

    fig, ax = plt.subplots(2, figsize=(2.1, 2.1), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.12)
    ax[1].plot(np.linspace(x_min, x_max, 10), np.ones(10), '--', color='darkgray')

    counts_denom, bins_denom = hists["DST_PFScouting_ZeroBias"]

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
            color=color
        )
        if trigger != "DST_PFScouting_ZeroBias":
            color = l[0].get_color()
            draw_ratio(
                counts, bins,
                counts_denom, bins_denom,
                ax=ax[1],
                label=trigger_label,
                norm=norm,
                color=color
            )

    ax[0].set_xlim([x_min, x_max])
    ax[0].set_ylim([y_min, y_max])

    if log_scale:
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
    ax[0].set_ylabel(f"Events{' [A.U.]' if norm else ''}", loc="top", fontsize=FONTSIZE, labelpad=0)
    ax[1].set_ylabel("Ratio to Zero Bias", loc="top", fontsize=FONTSIZE, labelpad=0)
    ax[1].set_xlabel(x_label, loc="right", fontsize=FONTSIZE, labelpad=2)
    legend_handles = [
        mpatches.Rectangle(
            (0, 0), 1, 1,
            fill=False,
            edgecolor=TRIGGER_COLORS[t],
            linewidth=1,
            label=TRIGGER_LABELS[t],
        )
        for t in triggers if t in hists
    ]
    ax[0].legend(handles=legend_handles, loc=leg_loc, frameon=False, fontsize=FONTSIZE, ncols=2, columnspacing=0.7)

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi=None,
        year="2024",
        com=13.6,
        fontsize=FONTSIZE,
        ax=ax[0],
    )

    ax[0].tick_params(axis='both', which='major', labelsize=FONTSIZE, length=4, pad=2)
    ax[1].tick_params(axis='both', which='major', labelsize=FONTSIZE, length=4, pad=2)
    ax[0].minorticks_off()
    ax[1].minorticks_off()
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

    make_plot(
        hists,
        triggers,
        defaults["x_label"],
        x_min, x_max, y_min, y_max,
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
        choices=["L1Mu", "ScoutingMuonVtx"],
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

    args = parser.parse_args()
    main(args)
