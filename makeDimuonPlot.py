import argparse
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os

hep.style.use('CMS')

DEFAULTS = {
    "x_min": 5e-2, "x_max": 1e3,
    "y_min": 5e-1, "y_max": 1e6,
}

TRIGGER_LABELS = {
    'DST_PFScouting_ZeroBias': 'Zero Bias',
    'DST_PFScouting_AXONominal': 'AXOL1TL V4 Medium',
    'DST_PFScouting_AXOVTight': 'AXOL1TL V4 Very Tight',
}

NORM = False


def load_root_hists(root_file, hist_key, triggers):
    hists = {}
    with uproot.open(root_file) as f:
        for trigger in triggers:
            key = f"{trigger}_{hist_key}"
            if key in f:
                counts, bins = f[key].to_numpy()
                hists[trigger] = (counts, bins)
            else:
                print(f"  WARNING: key '{key}' not found in ROOT file, skipping.")
    return hists


def draw_hist1d(counts, bins, ax=None, label="", rebin=1,
                norm=False, linestyle='solid', color=None):

    if rebin > 1:
        counts = counts[:len(counts) - len(counts) % rebin].reshape(-1, rebin).sum(axis=1)
        bins = bins[::rebin]
        if len(bins) != len(counts) + 1:
            bins = np.append(bins[:len(counts)], bins[len(counts)])

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
    return l


def main(args):

    x_min = args.x_min if args.x_min is not None else DEFAULTS["x_min"]
    x_max = args.x_max if args.x_max is not None else DEFAULTS["x_max"]
    y_min = args.y_min if args.y_min is not None else DEFAULTS["y_min"]
    y_max = args.y_max if args.y_max is not None else DEFAULTS["y_max"]

    triggers = [
        "DST_PFScouting_AXONominal",
        "DST_PFScouting_AXOVTight",
        "DST_PFScouting_ZeroBias",
    ]

    hists = load_root_hists(args.input, "ScoutingMuonVtx_ScoutingMuonVtx_mass", triggers)

    fig, ax = plt.subplots(figsize=(14, 8))

    for trigger in triggers:
        if trigger not in hists:
            continue
        counts, bins = hists[trigger]
        draw_hist1d(counts, bins, ax=ax, label=TRIGGER_LABELS[trigger], rebin=3, norm=NORM)

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.legend(loc="upper right", frameon=False, fontsize=16)
    ax.set_ylabel(f"Events{' [A.U.]' if NORM else ''}", loc="top", fontsize=25)
    ax.set_xlabel(r"HLT Scouting $m_{\mu\mu}$ [GeV]", fontsize=25)
    ax.text(60, 4e4, r"$p_T^\mu>3$ GeV, $|\eta|<2.4$", fontsize=16)

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi=1.62,
        year="2024",
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
        description="Make a plot of diobject invariant mass"
    )
    parser.add_argument(
        "--input",
        default="histograms/hist.root",
        help="Input .root histogram file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Full output path prefix, e.g. plots/dimuon_mass (extensions .pdf/.png added automatically)"
    )
    parser.add_argument("--x-min", type=float, default=None)
    parser.add_argument("--x-max", type=float, default=None)
    parser.add_argument("--y-min", type=float, default=None)
    parser.add_argument("--y-max", type=float, default=None)

    args = parser.parse_args()
    main(args)