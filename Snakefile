rule all:
   input:
      "outputs/AXOL1TL_v4_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v4_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
      "outputs/L1Jet_mult_nPV10.pdf",
      "outputs/L1Jet_mult_nPV10.png",
      "outputs/L1EG_mult_nPV10.pdf",
      "outputs/L1EG_mult_nPV10.png",
      "outputs/L1Mu_mult_nPV10.pdf",
      "outputs/L1Mu_mult_nPV10.png",
      "outputs/l1_ht_dist_nPV10.pdf",
      "outputs/l1_ht_dist_nPV10.png",
      "outputs/l1_met_dist_nPV10.pdf",
      "outputs/l1_met_dist_nPV10.png",
      "outputs/l1_ht_purity_nPV10.pdf",
      "outputs/l1_ht_purity_nPV10.png",
      "outputs/dimuon_mass_nPV10.pdf",
      "outputs/dimuon_mass_nPV10.png",

rule axo_style_score_plots:
   input:
      "inputs/CICADA2024_CICADAScore_plot_info.pkl",
      "inputs/axol1tl_v4_AXOScore_plot_info.pkl",
      "make_axo_style_score_plots.py",
   output:
      "outputs/AXOL1TL_v4_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v4_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
   shell:
      "python3 make_axo_style_score_plots.py --output outputs/"

rule correlation_plots:
   input:
      "inputs/correlation_dict.pkl",
      "make_correlation_plots.py",
   output:
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
   shell:
      "python3 make_correlation_plots.py --input inputs/correlation_dict.pkl --output outputs/"

rule obj_mult_plots:
   input:
      "inputs/hists_plotA_plotB_plotC.root",
      "makeObjMultPlots.py",
   output:
      "outputs/L1Jet_mult_nPV10.pdf",
      "outputs/L1Jet_mult_nPV10.png",
      "outputs/L1EG_mult_nPV10.pdf",
      "outputs/L1EG_mult_nPV10.png",
      "outputs/L1Mu_mult_nPV10.pdf",
      "outputs/L1Mu_mult_nPV10.png",
   shell:
      "python3 makeObjMultPlots.py --object L1Jet --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Jet_mult_nPV10 && "
      "python3 makeObjMultPlots.py --object L1EG  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1EG_mult_nPV10 && "
      "python3 makeObjMultPlots.py --object L1Mu  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Mu_mult_nPV10"

rule l1_dist_plots:
   input:
      "inputs/hists_plotD_plotE.root",
      "makeObjMultPlots.py",
   output:
      "outputs/l1_ht_dist_nPV10.pdf",
      "outputs/l1_ht_dist_nPV10.png",
      "outputs/l1_met_dist_nPV10.pdf",
      "outputs/l1_met_dist_nPV10.png",
   shell:
      "python3 makeObjMultPlots.py --input inputs/hists_plotD_plotE.root --output outputs/l1_ht_dist_nPV10 --object L1HT && "
      "python3 makeObjMultPlots.py --input inputs/hists_plotD_plotE.root --output outputs/l1_met_dist_nPV10 --object L1MET"

rule ht_purity_plot:
   input:
      "inputs/hists_plotF.root",
      "makeHTPurityPlot.py",
   output:
      "outputs/l1_ht_purity_nPV10.pdf",
      "outputs/l1_ht_purity_nPV10.png",
   shell:
      "python3 makeHTPurityPlot.py --input inputs/hists_plotF.root --output outputs/l1_ht_purity_nPV10"

rule ht_purity_plot_norm:
   input:
      "inputs/hists_plotF.root",
      "makeHTPurityPlot_norm.py",
   output:
      "outputs/l1_ht_purity_nPV10_norm.pdf",
      "outputs/l1_ht_purity_nPV10_norm.png",
   shell:
      "python3 makeHTPurityPlot_norm.py --input inputs/hists_plotF.root --output outputs/l1_ht_purity_nPV10_norm"

rule dimuon_mass_plot:
   input:
      "inputs/hists_plotG.root",
      "makeDimuonPlot.py",
   output:
      "outputs/dimuon_mass_nPV10.pdf",
      "outputs/dimuon_mass_nPV10.png",
   shell:
      "python3 makeDimuonPlot.py --input inputs/hists_plotG.root --output outputs/dimuon_mass_nPV10"
