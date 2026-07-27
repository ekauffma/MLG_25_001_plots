import argparse
import uproot
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.special import erf

hep.style.use('CMS')

DEFAULTS = {
    "x_min": 1e0, "x_max": 1e3,
    "y_min": 1e-8, "y_max": 1e1,
}

TRIGGER_LABELS = {
    'DST_PFScouting_AXONominal': 'AXOL1TL',
    'DST_PFScouting_ZeroBias': 'Zero Bias',
}

TRIGGER_COLORS = {
    'DST_PFScouting_ZeroBias' : '#1845fb',
    'DST_PFScouting_AXONominal': '#86c8dd',
}

TRIGGER_SCALING = {
    'DST_PFScouting_ZeroBias' : 1,
    'DST_PFScouting_AXONominal': 1,
}

RESONANCE_FITS = {
    "J/psi": {
        "fit_range": (2.8, 3.4),
        "p0": (500, 3.1, 0.05, 200, 0.5),
        "bounds": (
            [0,      2.9,   0.001, 0,      0.0001],   # lower bounds
            [50000,  3.3,   0.5,   50000,  10]        # upper bounds
        ),
    },
    "Z": {
        "fit_range": (70.0, 110.0),
        "p0": (100, 91.2, 2.5, 50, 0.02),
        "bounds": (
            [0,     85,   0.1,   0,      0.0001],
            [50000, 105,  10,    50000,  10]
        )
    },
}

TRIGGER_ZORDER = {
    'DST_PFScouting_ZeroBias': 2,
    'DST_PFScouting_AXONominal': 3,
}

NORM = True

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
                norm=False, linestyle='solid', color=None, scale=1.0, zorder=2):

    if rebin > 1:
        counts = counts[:len(counts) - len(counts) % rebin].reshape(-1, rebin).sum(axis=1)
        bins = bins[::rebin]
        if len(bins) != len(counts) + 1:
            bins = np.append(bins[:len(counts)], bins[len(counts)])

    # Compute Poisson errors on raw counts first, then apply scale factor
    norm_factor = np.sum(counts) * np.diff(bins) if norm else 1
    _counts = counts * scale / norm_factor if norm else counts * scale
    errs = np.sqrt(counts) * scale / norm_factor if norm else np.sqrt(counts) * scale
    _errs = np.where(_counts == 0, 0, errs)

    bin_centres = 0.5 * (bins[1:] + bins[:-1])

    if color is not None:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="", color=color, zorder=zorder+1)
    else:
        l = ax.errorbar(x=bin_centres, y=_counts, yerr=_errs, linestyle="", zorder=zorder+1)
    color = l[0].get_color()
    ax.step(
        x=bins[:-1], y=_counts, where="post", label=label,
        color=color, linestyle=linestyle, zorder=zorder, linewidth=2
    )
    return l

def gauss_plus_exp(x, N_sig, mu, sigma, N_bkg, lam):
    """Gaussian signal + falling exponential background."""
    signal = N_sig * np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    background = N_bkg * lam * np.exp(-lam * x)
    return signal + background

def fit_resonance(counts, bins, fit_range, p0, bounds, rebin=1, label="", trigger="", out_prefix=""):
    """
    Fit a Gaussian + exponential model to a histogram in fit_range.
    Returns (signal_efficiency, fit_params, fit_errors) and saves a debug plot.

    signal_efficiency = integral of Gaussian within ±2sigma / total counts in fit_range
    """
    # Apply same rebinning as main plot
    if rebin > 1:
        counts = counts[:len(counts) - len(counts) % rebin].reshape(-1, rebin).sum(axis=1)
        bins = bins[::rebin]
        if len(bins) != len(counts) + 1:
            bins = np.append(bins[:len(counts)], bins[len(counts)])

    bin_centres = 0.5 * (bins[1:] + bins[:-1])

    # Mask to fit range
    mask = (bin_centres >= fit_range[0]) & (bin_centres <= fit_range[1])
    x = bin_centres[mask]
    y = counts[mask]
    yerr = np.where(y > 0, np.sqrt(y), 1.0)  # Poisson errors, floor at 1

    try:
        popt, pcov = curve_fit(
            gauss_plus_exp, x, y,
            p0=p0,
            sigma=yerr,
            absolute_sigma=True,
            maxfev=10000,
            bounds=bounds,
        )
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError as e:
        print(f"  [WARNING] Fit failed for {trigger} / {label}: {e}")
        return None, None, None

    N_sig, mu, sigma, N_bkg, lam = popt

    # Signal efficiency: fraction of total counts in fit window that are signal
    # Signal integral within ±2sigma window
    total_in_window = np.sum(y)
    # Gaussian is normalised in gauss_plus_exp, so integral over all x = N_sig
    # Within ±2sigma: erf(2/sqrt(2)) ~ 0.9545
    sig_fraction_in_2sigma = erf(2.0 / np.sqrt(2))
    signal_in_2sigma = N_sig * sig_fraction_in_2sigma

    # Background integral within ±2sigma window
    lo, hi = mu - 2 * sigma, mu + 2 * sigma
    bkg_in_2sigma = N_bkg * (np.exp(-lam * lo) - np.exp(-lam * hi))

    efficiency = signal_in_2sigma / (signal_in_2sigma + bkg_in_2sigma) if (signal_in_2sigma + bkg_in_2sigma) > 0 else 0.0

    # --- Uncertainty propagation via numerical gradients ---
    def efficiency_from_params(params):
        N_sig, mu, sigma, N_bkg, lam = params
    
        sig_fraction_in_2sigma = erf(2.0 / np.sqrt(2))
        S = N_sig * sig_fraction_in_2sigma
    
        lo, hi = mu - 2 * sigma, mu + 2 * sigma
        B = N_bkg * (np.exp(-lam * lo) - np.exp(-lam * hi))
    
        denom = S + B
        return S / denom if denom > 0 else 0.0
    
    
    # Numerical gradient
    eps = 1e-6
    params = np.array(popt)
    grad = np.zeros_like(params)
    
    for i in range(len(params)):
        dp = np.zeros_like(params)
        dp[i] = eps * (abs(params[i]) + 1e-6)
    
        f_plus = efficiency_from_params(params + dp)
        f_minus = efficiency_from_params(params - dp)
    
        grad[i] = (f_plus - f_minus) / (2 * dp[i])
    
    # Variance: g^T Cov g
    eff_variance = grad @ pcov @ grad
    eff_error = np.sqrt(eff_variance) if eff_variance > 0 else 0.0

    # --- Debug plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    hep.style.use('CMS')
    x_fine = np.linspace(fit_range[0], fit_range[1], 500)
    ax.errorbar(x, y, yerr=yerr, fmt='o', color='black', markersize=4, label='Data')
    ax.plot(x_fine, gauss_plus_exp(x_fine, *popt), 'r-', label='Signal + Bkg fit')
    ax.plot(
        x_fine,
        popt[3] * popt[4] * np.exp(-popt[4] * x_fine),
        'b--', label='Background only'
    )
    ax.axvspan(lo, hi, alpha=0.15, color='green', label=r'$\pm2\sigma$ window')
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    ax.set_ylabel("Counts")
    ax.set_title(f"{TRIGGER_LABELS.get(trigger, trigger)} — {label} fit\n"
                 f"$\\mu={mu:.3f}$, $\\sigma={sigma:.3f}$, S/(S+B)|_{{2\\sigma}}={efficiency:.3f}")
    ax.legend(fontsize=10)
    if out_prefix:
        os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
        safe_trigger = trigger.replace("/", "_")
        safe_label = label.replace("/", "_").replace("$", "").replace("\\", "").replace("{", "").replace("}", "")
        debug_path = f"{out_prefix}_debug_{safe_trigger}_{safe_label}"
        fig.savefig(f"{debug_path}.png", bbox_inches="tight")
        fig.savefig(f"{debug_path}.pdf", bbox_inches="tight")
        print(f"  Saved debug plot: {debug_path}.png")
    plt.close(fig)

    return efficiency, eff_error, popt, perr

def run_signal_efficiency_fits(hists, triggers, out_prefix, rebin=3):
    print("\n" + "=" * 60)
    print("Signal efficiency fits (Gaussian + exponential)")
    print("Efficiency = S / (S+B) within ±2sigma of fitted mean")
    print("=" * 60)
 
    for res_name, cfg in RESONANCE_FITS.items():
        print(f"\n  Resonance: {res_name}")
        for trigger in triggers:
            if trigger not in hists:
                continue
            counts, bins = hists[trigger]
            eff, eff_err, popt, perr  = fit_resonance(
                counts.copy(), bins.copy(),
                fit_range=cfg["fit_range"],
                p0=cfg["p0"],
                bounds=cfg["bounds"],
                rebin=rebin,
                label=res_name,
                trigger=trigger,
                out_prefix=out_prefix,
            )
            if eff is None:
                print(f"    {TRIGGER_LABELS.get(trigger, trigger):20s}: fit failed")
            else:
                N_sig, mu, sigma, N_bkg, lam = popt
                N_sig_err, mu_err, sigma_err, _, _ = perr
                print(f"    {TRIGGER_LABELS.get(trigger, trigger):20s}: "
                      f"eff = {eff:.4f} +/- {eff_err:.4f} |  "
                      f"mu = {mu:.3f} ± {mu_err:.3f} GeV  |  "
                      f"sigma = {sigma:.3f} ± {sigma_err:.3f} GeV  |  "
                      f"N_sig = {N_sig:.1f} ± {N_sig_err:.1f}")
    print("=" * 60 + "\n")

def main(args):

    x_min = args.x_min if args.x_min is not None else DEFAULTS["x_min"]
    x_max = args.x_max if args.x_max is not None else DEFAULTS["x_max"]
    y_min = args.y_min if args.y_min is not None else DEFAULTS["y_min"]
    y_max = args.y_max if args.y_max is not None else DEFAULTS["y_max"]

    triggers = [
        "DST_PFScouting_ZeroBias",
        "DST_PFScouting_AXONominal",
    ]

    hists = load_root_hists(args.input, "ScoutingMuonVtx_ScoutingMuonVtx_mass", triggers)

    run_signal_efficiency_fits(hists, triggers, "outputs/dimuon_efficiency_", rebin=1)

    fig, ax = plt.subplots(figsize=(14, 6))

    for trigger in triggers:
        if trigger not in hists:
            continue
        color = TRIGGER_COLORS[trigger]
        counts, bins = hists[trigger]
        scaled_counts = counts * TRIGGER_SCALING.get(trigger, 1)
        draw_hist1d(scaled_counts, bins, ax=ax, label=TRIGGER_LABELS[trigger], rebin=6, norm=NORM, color=color, zorder=TRIGGER_ZORDER[trigger])

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])

    resonances = [
        (3.1, r"$J/\psi$"),
        (91.2, r"$Z$")
    ]
    
    for x, label in resonances:
        ax.axvline(x=x, color='black', linestyle='dotted', linewidth=1.5)
        ax.text(
            x * 0.95, 1e0, label,
            ha='right',
            rotation=0,
            verticalalignment='bottom',
            fontsize=18
        )


    legend_order = [
        "DST_PFScouting_ZeroBias",
        "DST_PFScouting_AXONominal",
    ]
    legend_handles = [
        mpatches.Rectangle(
            (0, 0), 1, 1,
            fill=False,
            edgecolor=TRIGGER_COLORS[t],
            linewidth=2,
            label=TRIGGER_LABELS[t],
        )
        for t in legend_order if t in hists
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=16)
    ax.set_ylabel(f"Events{' [A.U.]' if NORM else ''}", loc="top", fontsize=25)
    ax.set_xlabel(r"Reconstructed Muon $m_{\mu\mu}$ [GeV]", fontsize=25)
    ax.text(250, 5e-2, r"$p_T^\mu>3$ GeV, $|\eta|<2.4$", fontsize=16)

    hep.cms.label(
        "Preliminary",
        data=True,
        rlabel="2024 (13.6 TeV)",
        fontsize=22,
        ax=ax,
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
