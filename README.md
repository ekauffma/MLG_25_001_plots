# Environment

```
cmsrel CMSSW_15_1_0_pre1
cd CMSSW_15_1_0_pre1/src/
git clone https://github.com/aloeliger/MLG_25_001_plots.git
cd MLG_15_001_plots/
python3 -m venv plot_env
source plot_env/bin/activate
python3 -m pip install -r requirements.txt
```

# Running all plots
## Score distributions, rate plots, and ROC Curves (Lino)
Skimmed ntuples, containing only the cariables needed for the plots are stored here
```
/eos/user/l/ligerlac/cicada_data/skimmed-2025-07-09/
```
If needed, these can be (re-)created like so (takes roughly 1 hour):
```python
python skim-inputs-mp.py --input "/eos/cms/store/group/phys_exotica/axol1tl/CICADANtuples/ZeroBias/*01May*/*/*/*.root" --output ZB.root
python skim-inputs-mp.py --input "/eos/cms/store/group/phys_exotica/axol1tl/CICADANtuples/TT_TuneCP5_13p6TeV_powheg-pythia8/*01May*/*/*/*.root" --output TT.root
...
```
The actual plotting happens in the notebook `ad-paper-plots.ipynb`

## Other plots (Andrew & Elliott)
Changed plots can be run via snakemake commands. In particular:

```
snakemake all -c<number of threads>
```

Should remake any plots for which inputs or code has changed.
