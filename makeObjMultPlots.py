import argparse
import uproot
import numpy as np
import boost_histogram as bh
import matplotlib.pyplot as plt
import mplhep as hep
import os

hep.style.use('CMS')

triggers = [
    "DST_PFScouting_ZeroBias",
    "DST_PFScouting_JetHT",
    "DST_PFScouting_CICADAMedium",
    "DST_PFScouting_DoubleMuon",
    "DST_PFScouting_AXONominal",
]

TRIGGER_LABELS = {
    'DST_PFScouting_ZeroBias' : 'Zero Bias',
    'DST_PFScouting_JetHT': 'Jet HT',
    'DST_PFScouting_CICADAMedium': 'CICADA Medium',
    'DST_PFScouting_DoubleMuon': 'Double Muon',
    'DST_PFScouting_AXONominal': 'AXO Medium'
}

# Default axis limits and x-labels per object type
OBJ_DEFAULTS = {
    "L1Mu": {
        "hist_key": "L1Mu_mult",
        "x_label": r"$N_{\text{L1Mu}}$",
        "x_min": -0.5, "x_max": 8.5,
        "y_min": 5e-10, "y_max": 5e2,
    },
    "L1EG": {
        "hist_key": "L1EG_mult",
        "x_label": r"$N_{\text{L1EG}}$",
        "x_min": -0.5, "x_max": 12.5,
        "y_min": 5e-5, "y_max": 5e0,
    },
    "L1Jet": {
        "hist_key": "L1Jet_mult",
        "x_label": r"$N_{\text{L1Jet}}$",
        "x_min": -0.5, "x_max": 12.5,
        "y_min": 5e-8, "y_max": 5e1,
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
        color=color, linestyle=linestyle
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
        color=color, linestyle='solid'
    )
    return l


def make_plot(hists, triggers, x_label,
              x_min, x_max, y_min, y_max, output,
              log_scale=True, norm=False, leg_loc='upper right'):

    fig, ax = plt.subplots(2, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    ax[1].plot(np.linspace(x_min, x_max, 10), np.ones(10), '--', color='darkgray')

    counts_denom, bins_denom = hists["DST_PFScouting_ZeroBias"]

    for trigger in triggers:
        if trigger not in hists:
            continue
        counts, bins = hists[trigger]
        color = 'k' if trigger == "DST_PFScouting_ZeroBias" else None
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
    ax[0].set_ylabel(f"Events{' [A.U.]' if norm else ''}", loc="top", fontsize=20)
    ax[1].set_ylabel("Ratio to Zero Bias", loc="top", fontsize=20)
    ax[1].set_xlabel(x_label, loc="right", fontsize=20)
    ax[0].legend(loc=leg_loc, frameon=False, fontsize=16, ncols=3, columnspacing=0.7)

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi=11.45,
        year="2024",
        com=13.6,
        fontsize=20,
        ax=ax[0],
    )

    fig.subplots_adjust(hspace=0)

    out_dir = os.path.dirname(output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(f"{output}.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(f"{output}.png", format="png", bbox_inches="tight")
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
        choices=["L1Mu", "L1EG", "L1Jet"],
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
