import argparse
from coffea.util import load
from hist.tag import rebin as reb
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os

hep.style.use('CMS')

# Default axis limits and x-labels per object type
DEFAULTS = {
    "x_min": 0, "x_max": 2000,
    "y_min": 5e0, "y_max": 4e7,
}

TRIGGER_LABELS = {
    "DST_PFScouting_AXONominal": 'AXO Medium',
    "DST_PFScouting_CICADAMedium": 'CICADA Medium',
    "pure_L1_DST_PFScouting_AXONominal": 'AXO Medium pure w.r.t. L1',
    "pure_L1_DST_PFScouting_CICADAMedium": 'CICADA Medium pure w.r.t. L1',
}

NORM = False

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

def main(args):

    hist_result = load(args.input)

    defaults = DEFAULTS

    x_min = args.x_min if args.x_min is not None else defaults["x_min"]
    x_max = args.x_max if args.x_max is not None else defaults["x_max"]
    y_min = args.y_min if args.y_min is not None else defaults["y_min"]
    y_max = args.y_max if args.y_max is not None else defaults["y_max"]

    fig, ax = plt.subplots(figsize=(9, 7))

    triggers = [
        "DST_PFScouting_AXONominal",
        "DST_PFScouting_CICADAMedium",
        "pure_L1_DST_PFScouting_AXONominal",
        "pure_L1_DST_PFScouting_CICADAMedium"
    ]

    for trigger in triggers:

        draw_hist1d(
            hist_in=hist_result["hists"]["l1_ht"]["2024I_paper", trigger, :],
            ax=ax,
            label=TRIGGER_LABELS[trigger],
            rebin=5,
            norm=NORM,
        )

    ax.set_yscale("log")
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.legend(loc="upper right", frameon=False, fontsize=16)
    ax.set_ylabel(f"Events{' [A.U.]' if NORM else ''}", loc="top", fontsize=25)
    ax.set_xlabel(r"L1 $H_T$ [GeV]", fontsize=25)
    
    hep.cms.label(
        "Preliminary",
        data=True,
        lumi=0.83,
        year="2024 Run 386924",
        com=13.6,
        fontsize=18,
    )

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{args.output}.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(f"{args.output}.png", format="png", bbox_inches="tight")
    print(f"Saved {args.output}.pdf and {args.output}.png")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Make a plot of L1 HT for AXO and CICADA including pure events"
    )
    parser.add_argument(
        "--input",
        default="histograms/hist_result_plotF.pkl",
        help="Input .pkl histogram file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Full output path prefix, e.g. plots/l1_ht_purity (extensions .pdf/.png added automatically)"
    )
    parser.add_argument("--x-min", type=float, default=None, help="x-axis minimum")
    parser.add_argument("--x-max", type=float, default=None, help="x-axis maximum")
    parser.add_argument("--y-min", type=float, default=None, help="y-axis minimum")
    parser.add_argument("--y-max", type=float, default=None, help="y-axis maximum")

    args = parser.parse_args()
    main(args)