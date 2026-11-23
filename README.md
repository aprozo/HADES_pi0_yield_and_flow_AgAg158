# π⁰ → γγ in Ag+Ag 1.58A GeV (HADES mar19) with ECAL

Neutral pion yields and flow (v1, v2) from diphotons, compared to charged pions and models.

## Layout
| dir | what |
|---|---|
| `ecal_calibration/` | ECAL energy + time calibration → `mar19EmcCalParNew{Energy,Time}.txt` |
| `efficiency/` | ECAL efficiency and acceptance maps |
| `pi0_peak_calibration/` | energy non-linearity correction from π⁰ peak |
| `photon_beta/` | photon β cut per cell |
| `dst_analysis/` | DST tasks: embedding eff×acc, lepton β cut, event plane |
| `correction_data/` | inputs for the analysis (`PATH_TO_CORRECTION_DATA`) |
| `analysis/` | main pipeline + systematics |

## Pipeline (`analysis/`)
```
exp_sim_to_photons → photons_to_histograms(_flow) → histograms_to_pt_y → ana_results → compareToChargedPions
urqmd_to_pt_y, embedding_to_effacc: reference and eff×acc
```
```bash
bash run.sh parameters.conf            # one cut set (set RUN_* flags inside)
bash run_all_parameters_sets.sh        # all sets in parameters_files/ + ana_results
```
