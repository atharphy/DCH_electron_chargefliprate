# Flip Rate Framework

Author: Athar Ahmad  
Contact: athar.ahmad@cern.ch

## Overview

This repository contains a small ROOT-based framework for measuring and validating the electron charge-flip rate in Run 2 skims.

Electron charge mismeasurement can convert a true opposite-sign `Z → ee` event into a reconstructed same-sign `ee` event. This effect is small, but it matters for analyses that use same-sign electrons or categories where charge misidentification can create background contamination.

The framework has two main parts:

1. Measure a global electron charge-flip probability in Data and MC.
2. Validate the correction by making Data/MC stacked plots of the dielectron mass before and after applying the charge-flip correction.

## Repository structure

```text
flip_rate_framework/
├── runChargeFlipSF_unbinned_allYears.C
├── ChargeFlipRateTest.C
├── filemap/
│   └── FileMap.h
└── include/
    ├── MyBranch.C
    ├── Kinematics.C
    ├── Kinematics.h
    ├── Xsections.C
    ├── MET_split.C
    └── other helper files
```

## Input skims

The charge-flip rate measurement reads skims from:

```text
/eos/uscms/store/user/aahmad2/run2_files/fliprate_skims
```

The script uses the FNAL EOS redirector:

```text
root://cmseos.fnal.gov/
```

so the files can be accessed from outside the LPC interactive nodes if XRootD access works.

The expected input layout is:

```text
Data_2016preVFP/
Data_2016postVFP/
Data_2017/
Data_2018/
MC_2016preVFP/
MC_2016postVFP/
MC_2017/
MC_2018/
```

The framework treats the Run 2 eras separately:

```text
2016preVFP
2016postVFP
2017
2018
```

The combined `Run2` result is made by adding the already-normalized era histograms. No separate Run2 skim is read.

## Why charge-flip rates are needed

A true `Z → ee` event has opposite-sign electrons. If one electron has its charge reconstructed incorrectly, the event becomes same-sign at reconstruction level.

For a per-electron charge-flip probability `f`:

```math
P(\mathrm{SS}) = 2f(1-f)
```

because either electron can flip, but not both.

The reconstructed opposite-sign probability is:

```math
P(\mathrm{OS}) = (1-f)^2 + f^2
```

because either neither electron flips, or both electrons flip.

Therefore, the observed ratio

```math
R = \frac{N_{\mathrm{SS}}}{N_{\mathrm{OS}}}
```

is related to `f` by

```math
R = \frac{2f(1-f)}{(1-f)^2 + f^2}.
```

Solving this equation gives the charge-flip probability. The script uses the smaller physical solution, since charge-flip probabilities are expected to be much smaller than 50%.

The Data/MC scale factor is then:

```math
SF = \frac{f_{\mathrm{Data}}}{f_{\mathrm{MC}}}.
```

## Script 1: measuring the charge-flip rate

The measurement script is:

```text
runChargeFlipSF_unbinned_allYears.C
```

It measures one global charge-flip probability per era using clean `Z → ee` events.

Selection:

- dielectron channel only,
- both selected objects are already tight electrons in the skim,
- MC requires both electrons to be prompt using:
  ```cpp
  genPartFlav_1 == 1 && genPartFlav_2 == 1
  ```
- invariant mass requirement:
  ```text
  |m_ee - 91.2| < 10 GeV
  ```

The script counts:

```text
N_SS = same-sign ee events in the Z window
N_OS = opposite-sign ee events in the Z window
```

and converts `N_SS / N_OS` into a per-electron flip probability.

### Run the charge-flip measurement

From the repository directory:

```bash
root -l -b -q 'runChargeFlipSF_unbinned_allYears.C'
```

or inside ROOT:

```cpp
.L runChargeFlipSF_unbinned_allYears.C
runChargeFlipSF_unbinned_allYears()
```

### Output

The output files are written to:

```text
/eos/uscms/store/user/aahmad2/flip_rate_results
```

One file is produced per era:

```text
ChargeFlipSF_SIMPLE_2016preVFP.root
ChargeFlipSF_SIMPLE_2016postVFP.root
ChargeFlipSF_SIMPLE_2017.root
ChargeFlipSF_SIMPLE_2018.root
```

Each file contains:

```text
MC_rate
Data_rate
ChargeFlipSF
```

where:

```text
MC_rate      = global MC charge-flip probability
Data_rate    = global Data charge-flip probability
ChargeFlipSF = Data_rate / MC_rate
```

## Script 2: Data/MC stack validation

The validation script is:

```text
ChargeFlipRateTest.C
```

It makes Data/MC stacked plots of `m_ee` in the `ee` channel for three regions:

```text
ss_Zwin          same-sign ee in the Z window
os_Zwin          opposite-sign ee in the Z window
inclusive_Zwin   sign-inclusive ee in the Z window
```

The histogram binning is:

```text
40 bins from 70 to 110 GeV
```

so each bin is 1 GeV wide.

The script produces both nominal and charge-flip-corrected histograms in the same event loop if requested.

## Charge-flip correction method

The correction is applied only to MC. Data is left unchanged.

The script does not simply scale the same-sign yield. Instead, it performs stochastic charge migration.

If the measured Data charge-flip probability is larger than MC, some MC opposite-sign events are randomly migrated into the same-sign category. If the measured Data charge-flip probability is smaller than MC, some same-sign MC events are migrated back into the opposite-sign category.

This is important because increasing the same-sign yield should also remove a corresponding small number of events from the opposite-sign yield.

The migration uses:

```math
P_{\mathrm{SS}}(f) = 2f(1-f)
```

and

```math
P_{\mathrm{OS}}(f) = (1-f)^2 + f^2.
```

If Data has a larger same-sign probability than MC, the OS to SS migration probability is:

```math
P(\mathrm{OS}\rightarrow\mathrm{SS}) =
\frac{P_{\mathrm{SS}}(f_{\mathrm{Data}}) - P_{\mathrm{SS}}(f_{\mathrm{MC}})}
     {P_{\mathrm{OS}}(f_{\mathrm{MC}})}.
```

If Data has a smaller same-sign probability than MC, the SS to OS migration probability is:

```math
P(\mathrm{SS}\rightarrow\mathrm{OS}) =
\frac{P_{\mathrm{SS}}(f_{\mathrm{MC}}) - P_{\mathrm{SS}}(f_{\mathrm{Data}})}
     {P_{\mathrm{SS}}(f_{\mathrm{MC}})}.
```

Only the reconstructed charges are changed. The event kinematics and event weight are unchanged.

## Data/MC normalization

For MC, the event weight is:

```math
w =
w_{\mathrm{gen}}
\times w_{\mathrm{br}}
\times \frac{\sigma}{N_{\mathrm{gen}}}
\times \mathcal{L}
\times w_{\mathrm{PU}}
\times w_{\mathrm{prefiring}}
\times SF_{\mathrm{ID}}
\times SF_{\mathrm{ISO}}
\times SF_{\mathrm{trigger}}.
```

In the code, this is built from:

```cpp
Generator_weight
brWeight
xsecWeight
lumiPB
weightPUtruejson
L1PreFiringWeight_Nom
IDSF_1 * IDSF_2
ISOSF_1 * ISOSF_2
TrigSF_1 or TrigSF_2
```

Data events are assigned weight 1.

## FileMap

The validation script reads samples from:

```text
filemap/FileMap.h
```

The process names in `FileMap.h` should match the process labels used in the stack order:

```cpp
DY
DY10_50
VV
ZZ
TTbar
ttV
WJ
ST
VVV
QCD
other
data
```

For Data, `FileMap.h` may include both electron and muon datasets. The validation script keeps only electron-triggered Data files:

```cpp
SingleElectron
EGamma
```

This avoids mixing in `SingleMuon` data for the `ee` validation region.

## Run the validation plots

From the repository directory:

```bash
root -l -b -q 'ChargeFlipRateTest.C("true")'
```

The safer way is to open ROOT and run:

```cpp
.L ChargeFlipRateTest.C
ChargeFlipRateTest(true)
```

This produces both:

```text
nominal
withChargeFlip
```

outputs in one pass over the skims.

For nominal only:

```cpp
.L ChargeFlipRateTest.C
ChargeFlipRateTest(false)
```

The script will print a warning at the beginning because it loops over all Run 2 Data and MC files and may take a while.

## Validation output

Outputs are written to:

```text
/eos/uscms/store/user/aahmad2/flip_rate_stacks
```

ROOT files:

```text
Zmass_DataMCStacks_nominal.root
Zmass_DataMCStacks_withChargeFlip.root
```

PNG plots:

```text
/eos/uscms/store/user/aahmad2/flip_rate_stacks/nominal/<year>/
/eos/uscms/store/user/aahmad2/flip_rate_stacks/withChargeFlip/<year>/
```

where `<year>` is one of:

```text
2016preVFP
2016postVFP
2017
2018
Run2
```

Each year directory contains plots for:

```text
mEE_<year>_ss_Zwin.png
mEE_<year>_os_Zwin.png
mEE_<year>_inclusive_Zwin.png
```

## Recommended workflow

First measure the charge-flip rates:

```cpp
.L runChargeFlipSF_unbinned_allYears.C
runChargeFlipSF_unbinned_allYears()
```

Then run the Data/MC validation:

```cpp
.L ChargeFlipRateTest.C
ChargeFlipRateTest(true)
```

Then compare:

```text
nominal/Run2/
withChargeFlip/Run2/
```

especially:

```text
ss_Zwin
os_Zwin
inclusive_Zwin
```

The expected behavior is that the charge-flip correction mainly changes the same-sign and opposite-sign distributions. The inclusive sign distribution should change very little, since the method mostly migrates events between OS and SS.

## Notes

- The skims are assumed to already contain tight selected objects.
- No extra electron WP90 requirement is applied in these scripts.
- The charge-flip measurement uses only prompt DY MC electrons.
- The stack validation uses all MC backgrounds from `FileMap.h`.
- Run2 is obtained by summing the four era histograms after era-specific normalization.
- The stochastic charge migration uses a fixed random seed:
  ```cpp
  TRandom3 randCF(12345);
  ```
  so repeated runs should be reproducible.
