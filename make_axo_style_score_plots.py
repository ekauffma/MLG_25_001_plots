import argparse
import mplhep as hep
import numpy as np
import ROOT
import json
import pickle as pkl
import matplotlib.pyplot as plt

from rich.console import Console

console = Console()

FONTSIZE=8

def draw_axo_style_score_plot(
        hist_dict,
        output_path,
        score_name,
        x_axis_bounds=(0., 180.0),
        y_axis_bounds=(1e-1, 1e8),
        x_axis_label="Emulated CICADA Score",
        working_point_label = "CICADA Nominal",
        pure_label = "CICADA Unique",
):
    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(2.1, 2.1))
    fig.subplots_adjust(left=0.15, right=0.95, top=0.92, bottom=0.12)
    hep.cms.label(
        "Preliminary",
        data=True,
        rlabel="2024 (13.6 TeV)",
        fontsize=FONTSIZE,
        ax=ax,
    )

    overall_hist = hist_dict["overall"]
    working_point_hist = hist_dict['working']
    pure_hist = hist_dict['pure']

    overall_fig = hep.histplot(
        overall_hist,
        label='All Zero Bias',
        color="#5790FC",
        linewidth=1,
    )

    working_point_fig = hep.histplot(
        working_point_hist,
        label=working_point_label,
        linestyle='--',
        color="#F89C20",
        linewidth=1,
    )
    pure_score_fig = hep.histplot(
        pure_hist,
        label=pure_label,
        linestyle='--',
        color="#E42536",
        linewidth=1,
    )

    plt.subplots_adjust(left=0.18, right=0.98, top=0.93, bottom=0.12)
    plt.legend(loc='upper right', title='Zero Bias Triggered Events', frameon=False, fontsize=FONTSIZE, title_fontsize=FONTSIZE)
    plt.xlabel(x_axis_label, loc="right", fontsize=FONTSIZE, labelpad=2)
    plt.ylabel('Events', loc="top", fontsize=FONTSIZE, labelpad=0)
    plt.tick_params(axis='both', which='major', labelsize=FONTSIZE, length=4, pad=2)
    plt.minorticks_off()
    plt.yscale('log')
    plt.ylim(y_axis_bounds)
    for spine in plt.gca().spines.values():
        spine.set_linewidth(0.8)

    hist_name = f'{score_name}_axo_style_score_plot'

    plt.savefig(
        f'{output_path}/{hist_name}.png'
    )
    plt.savefig(
        f'{output_path}/{hist_name}.pdf'
    )
    plt.close()

def main(args):
    # Get the input file information we need
    console.log("Making AXO style score plots")
    # cicada_df = ROOT.RDataFrame("CICADA2024_CICADAScore_plot_info", args.input)
    # axo_df = ROOT.RDataFrame("axol1tl_v3_AXOScore_plot_info", args.input)
    with open("inputs/CICADA2024_CICADAScore_plot_info.pkl", 'rb') as theFile:
        cicada_plot_dict = pkl.load(theFile)
    with open("inputs/axol1tl_v4_AXOScore_plot_info.pkl", "rb") as theFile:
        axo_plot_dict = pkl.load(theFile)


    # Hand each off to the drawing function
    draw_axo_style_score_plot(
        cicada_plot_dict,
        args.output,
        score_name="CICADA_2024",
        x_axis_bounds=(0., 180.0),
        y_axis_bounds=(1e-1, 1e12),
        x_axis_label="Emulated CICADA Score",
        working_point_label = "CICADA Medium",
        pure_label = "CICADA Unique",
    )

    draw_axo_style_score_plot(
        axo_plot_dict,
        args.output,
        score_name = "AXOL1TL_v4",
        x_axis_bounds=(0., 2000.0),
        y_axis_bounds=(1e-1, 1e9),
        x_axis_label="Emulated AXOL1TL Score",
        working_point_label = "AXOL1TL Medium",
        pure_label = "AXOL1TL Unique",
    )


    # Done!
    console.log("Done with AXO style score plots")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument(
    #     '--input',
    #     required=True,
    #     nargs='?',
    #     help='Input file to run drawing from'
    # )

    parser.add_argument(
        '--output',
        required=True,
        nargs='?',
        help='output directory to store output image files to'
    )

    args = parser.parse_args()

    main(args)
