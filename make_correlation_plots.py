import argparse
import mplhep as hep
import numpy as np
import ROOT
import json
import pickle as pkl
import matplotlib.pyplot as plt

from rich.console import Console

console = Console()

def make_1D_correlation_plot(
        snapshot_dict,
        output_path
):
    x_labels = list(snapshot_dict.keys())
    correlations = list(snapshot_dict.values())
    x_points = np.arange(len(x_labels))

    plt.plot(
        x_points,
        correlations,
        linestyle='None',
        marker='o',
        color='#5790FC',
        markersize=10,
    )

    vertical_line = plt.axvline(
        x=2.0 - 0.5,
        color='grey',
        linestyle='--',
    )

    plt.xticks(x_points, x_labels, rotation=90)
    plt.ylabel('Correlation Coeff.')
    plt.xlabel('Signals')

    plt.subplots_adjust(bottom=0.35)

    hist_name = f'1D_correlation_plot'

    plt.savefig(
        f'{output_path}/{hist_name}.png'
    )
    plt.savefig(
        f'{output_path}/{hist_name}.pdf'
    )
    plt.close()
    
    

def main(args):
    # Get the input file information we need
    console.log("Making 1D correlation plots")
    with open(args.input, "rb") as theFile:
        snapshot_dict = pkl.load(theFile)

    #Pass the histogram snapshot off to the drawing function
    make_1D_correlation_plot(
        snapshot_dict,
        args.output
    )

    console.log("Done making 1D correlation plots")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
        nargs="?",
        help="Input file to make correlation plots from"
    )

    parser.add_argument(
        "--output",
        required=True,
        nargs="?",
        help="Output directory to store output image files to"
    )

    args = parser.parse_args()

    main(args)
