import argparse
import mplhep as hep
import numpy as np
import ROOT
import json
import pickle as pkl
import matplotlib.pyplot as plt

from rich.console import Console

console = Console()

def draw_axo_style_score_plot(
        hist_dict,
        output_path,
        score_name,
        x_axis_bounds=(0., 180.0),
        x_axis_label="Emulated CICADA Score",
        working_point_label = "CICADA Nominal",
        pure_label = "CICADA Unique",
):
    hep.style.use("CMS")
    hep.cms.text("Preliminary", loc=2)

    overall_hist = hist_dict["overall"]
    working_point_hist = hist_dict['working']
    pure_hist = hist_dict['pure']

    overall_fig = hep.histplot(
        overall_hist,
        label='All Zero Bias',
        color="#5790FC"
    )
    working_point_fig = hep.histplot(
        working_point_hist,
        label=working_point_label,
        color="#F89C20"
    )
    pure_score_fig = hep.histplot(
        pure_hist,
        label=pure_label,
        linestyle='--',
        color="#E42536",
    )

    plt.legend(loc='upper right', title='Zero Bias Triggered Events')
    plt.xlabel(x_axis_label)
    plt.ylabel('Events')
    plt.yscale('log')
    plt.ylim(1.0, np.max(overall_hist[0])*100.0)

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
    with open("inputs/axol1tl_v3_AXOScore_plot_info.pkl", "rb") as theFile:
        axo_plot_dict = pkl.load(theFile)
    
    
    # Hand each off to the drawing function
    draw_axo_style_score_plot(
        cicada_plot_dict,
        args.output,
        score_name="CICADA_2024",
        x_axis_bounds=(0., 180.0),
        x_axis_label="Emulated CICADA Score",
        working_point_label = "CICADA Nominal",
        pure_label = "CICADA Unique",
    )

    draw_axo_style_score_plot(
        axo_plot_dict,
        args.output,
        score_name = "AXOL1TL_v3",
        x_axis_bounds=(0., 2000.0),
        x_axis_label="Emulated AXOL1TL Score",
        working_point_label = "AXOL1TL Nominal",
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
