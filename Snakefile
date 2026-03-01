rule all:
   input:
      "outputs/AXOL1TL_v3_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v3_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",

rule axo_style_score_plots:
   input:
      "inputs/CICADA2024_CICADAScore_plot_info.pkl",
      "inputs/axol1tl_v3_AXOScore_plot_info.pkl",
   output:
      "outputs/AXOL1TL_v3_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v3_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
   shell:
      "python3 make_axo_style_score_plots.py --output outputs/"

rule correlation_plots:
   input:
      "inputs/correlation_dict.pkl",
   output:
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
   shell:
      "python3 make_correlation_plots.py --input inputs/correlation_dict.pkl --output outputs/"