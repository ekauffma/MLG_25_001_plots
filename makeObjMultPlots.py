import argparse
from coffea.util import load
from hist.tag import rebin as reb
import matplotlib.pyplot as plt
import numpy as np
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
        "y_min": 5e-10, "y_max": 1e2,
    },
    "L1EG": {
        "hist_key": "L1EG_mult",
        "x_label": r"$N_{\text{L1EG}}$",
        "x_min": -0.5, "x_max": 12.5,
        "y_min": 5e-5, "y_max": 5e1,
    },
    "L1Jet": {
        "hist_key": "L1Jet_mult",
        "x_label": r"$N_{\text{L1Jet}}$",
        "x_min": -0.5, "x_max": 12.5,
        "y_min": 5e-8, "y_max": 5e2,
    },
}


def draw_hist1d(hist_in=None, ax=None, label="", rebin=1,
                obj=None, norm=False, gg=None, linestyle='solid', color=None):

    hist_in = hist_in[reb(rebin)]

    hist_data = hist_in.to_numpy()
    ndim = len(hist_data)

    if ndim == 2:
        counts, bins = hist_data
    else:
        counts, _, bins = hist_data
        counts = np.sum(counts, axis=0)

    if len(counts) > 0:
        norm_factor = np.sum(counts) * np.diff(bins) if norm else 1
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
    else:
        if color is not None:
            l = ax.errorbar(x=[], y=[], yerr=[], drawstyle="steps-post", color=color)
        else:
            l = ax.errorbar(x=[], y=[], yerr=[], drawstyle="steps-post")
        color = l[0].get_color()
        ax.errorbar(x=[], y=[], drawstyle="steps-post", label=label, color=color, linestyle=linestyle)

    return l


def draw_ratio(hist_num, hist_denom, ax=None, color=None,
               label="", rebin=1, norm=False, gg=None):

    hist_num = hist_num[reb(rebin)]
    hist_denom = hist_denom[reb(rebin)]

    data_num = hist_num.to_numpy()
    data_denom = hist_denom.to_numpy()

    if len(data_denom) == 2:
        counts_denom, bins = data_denom
        counts_num, _ = data_num
    else:
        counts_denom, _, bins = data_denom
        counts_num, _, _ = data_num
        counts_denom = np.sum(counts_denom, axis=0)
        counts_num = np.sum(counts_num, axis=0)

    norm_factor_denom = np.sum(counts_denom) * np.diff(bins) if norm else 1
    counts_denom = counts_denom / norm_factor_denom if norm else counts_denom
    errs_denom = np.sqrt(counts_denom * (1 - counts_denom / norm_factor_denom)) / norm_factor_denom if norm else np.sqrt(counts_denom)

    norm_factor_num = np.sum(counts_num) * np.diff(bins) if norm else 1
    counts_num = counts_num / norm_factor_num if norm else counts_num
    errs_num = np.sqrt(counts_num * (1 - counts_num / norm_factor_num)) / norm_factor_num if norm else np.sqrt(counts_num)

    denom = np.where(counts_denom == 0, np.nan, counts_denom)
    ratio = counts_num / denom
    x = 0.5 * (bins[:-1] + bins[1:])

    error = ratio * np.sqrt((errs_num / counts_num)**2 + (errs_denom / counts_denom)**2)

    if color is not None:
        l = ax.errorbar(x=x, y=ratio, yerr=error, linestyle="", color=color)
    else:
        l = ax.errorbar(x=x, y=ratio, yerr=error, linestyle="")
    color = l[0].get_color()
    ax.errorbar(
        x=bins, y=np.append(ratio, ratio[-1]), drawstyle="steps-post", label=label,
        color=color, linestyle='solid'
    )
    return l

def make_plot(hist_result, thing_to_plot, triggers, x_label,
              x_min, x_max, y_min, y_max, output,
              rebin=1, log_scale=True, norm=False, leg_loc='upper right'):
    print(thing_to_plot)

    fig, ax = plt.subplots(2, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    ax[1].plot(np.linspace(x_min, x_max, 10), np.ones(10), '--', color='darkgray')
    hist_denom = hist_result['hists'][thing_to_plot][:, 'DST_PFScouting_ZeroBias', :].integrate("dataset")

    for trigger in triggers:
        print(trigger)
        color = 'k' if trigger == "DST_PFScouting_ZeroBias" else None
        trigger_label = TRIGGER_LABELS[trigger]
        l = draw_hist1d(
            hist_in=hist_result['hists'][thing_to_plot][:, trigger, :].integrate("dataset"),
            ax=ax[0],
            label=trigger_label,
            rebin=rebin,
            norm=norm,
            color=color
        )
        if trigger != "DST_PFScouting_ZeroBias":
            color = l[0].get_color()
            hist_num = hist_result['hists'][thing_to_plot][:, trigger, :].integrate("dataset")
            draw_ratio(
                hist_num,
                hist_denom,
                ax=ax[1],
                label=trigger_label,
                rebin=rebin,
                norm=norm,
                color=color
            )

    ax[0].set_xlim([x_min, x_max])
    ax[0].set_ylim([y_min, y_max])

    if log_scale:
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
    ax[0].set_ylabel(f"Events{' [A.U.]' if norm else ''}", loc="top", fontsize=14)
    ax[1].set_ylabel("Ratio to Zero Bias", loc="top", fontsize=14)
    ax[1].set_xlabel(x_label, fontsize=14, loc="right")
    ax[0].legend(loc=leg_loc, frameon=False, fontsize=12)

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi=11.45,
        year="2024",
        com=13.6,
        fontsize=12,
        ax=ax[0],
    )

    fig.subplots_adjust(hspace=0)

    # output is used directly as the full path prefix
    out_dir = os.path.dirname(output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(f"{output}.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(f"{output}.png", format="png", bbox_inches="tight")
    print(f"Saved {output}.pdf and {output}.png")


def main(args):

    hist_result = load(args.input)

    defaults = OBJ_DEFAULTS[args.object]

    x_min = args.x_min if args.x_min is not None else defaults["x_min"]
    x_max = args.x_max if args.x_max is not None else defaults["x_max"]
    y_min = args.y_min if args.y_min is not None else defaults["y_min"]
    y_max = args.y_max if args.y_max is not None else defaults["y_max"]

    make_plot(
        hist_result,
        defaults["hist_key"],
        triggers,
        defaults["x_label"],
        x_min, x_max, y_min, y_max,
        args.output,
        rebin=1,
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
        default="histograms/hist_result_plotA_plotB_plotC.pkl",
        help="Input .pkl histogram file"
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