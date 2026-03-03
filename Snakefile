rule all:
   input:
      "outputs/AXOL1TL_v3_axo_style_score_plot.pdf",
      "outputs/AXOL1TL_v3_axo_style_score_plot.png",
      "outputs/CICADA_2024_axo_style_score_plot.pdf",
      "outputs/CICADA_2024_axo_style_score_plot.png",
      "outputs/1D_correlation_plot.pdf",
      "outputs/1D_correlation_plot.png",
      "outputs/L1Jet_mult.pdf",
      "outputs/L1Jet_mult.png",
      "outputs/L1EG_mult.pdf",
      "outputs/L1EG_mult.png",
      "outputs/L1Mu_mult.pdf",
      "outputs/L1Mu_mult.png",
      "outputs/l1_ht_dist.pdf",
      "outputs/l1_ht_dist.png",
      "outputs/l1_met_dist.pdf",
      "outputs/l1_met_dist.png",
      "outputs/l1_ht_purity.pdf",
      "outputs/l1_ht_purity.png",
      "outputs/dimuon_mass.pdf",
      "outputs/dimuon_mass.png",

rule axo_style_score_plots:
   input:
      "inputs/CICADA2024_CICADAScore_plot_info.pkl",
      "inputs/axol1tl_v3_AXOScore_plot_info.pkl",
      "make_axo_style_score_plots.py",
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
      "outputs/L1Jet_mult.pdf",
      "outputs/L1Jet_mult.png",
      "outputs/L1EG_mult.pdf",
      "outputs/L1EG_mult.png",
      "outputs/L1Mu_mult.pdf",
      "outputs/L1Mu_mult.png",
   shell:
      "python3 makeObjMultPlots.py --object L1Jet --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Jet_mult && "
      "python3 makeObjMultPlots.py --object L1EG  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1EG_mult && "
      "python3 makeObjMultPlots.py --object L1Mu  --input inputs/hists_plotA_plotB_plotC.root --output outputs/L1Mu_mult"

rule l1_dist_plots:
   input:
      "inputs/hists_plotD_plotE.root",
      "makeL1DistPlot.py",
   output:
      "outputs/l1_ht_dist.pdf",
      "outputs/l1_ht_dist.png",
      "outputs/l1_met_dist.pdf",
      "outputs/l1_met_dist.png",
   shell:
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE.root --output outputs/l1_ht_dist --observable ht && "
      "python3 makeL1DistPlot.py --input inputs/hists_plotD_plotE.root --output outputs/l1_met_dist --observable met"

rule ht_purity_plot:
   input:
      "inputs/hists_plotF.root",
      "makeHTPurityPlot.py",
   output:
      "outputs/l1_ht_purity.pdf",
      "outputs/l1_ht_purity.png",
   shell:
      "python3 makeHTPurityPlot.py --input inputs/hists_plotF.root --output outputs/l1_ht_purity"

rule dimuon_mass_plot:
   input:
      "inputs/hists_plotG.root",
      "makeDimuonPlot.py",
   output:
      "outputs/dimuon_mass.pdf",
      "outputs/dimuon_mass.png",
   shell:
      "python3 makeDimuonPlot.py --input inputs/hists_plotG.root --output outputs/dimuon_mass"