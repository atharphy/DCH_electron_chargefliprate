#include "include/MyBranch.C"
#include "include/Kinematics.C"

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

using namespace std;

// ============================================================================
// Global charge-flip scale factor measurement
//
// This macro measures a single unbinned electron charge-flip probability in
// Data and MC using Z->ee events. The output scale factor is
//
//        SF = f(Data) / f(MC)
//
// where f is the measured per-electron charge-flip probability.
// ============================================================================

// Files are accessed through the FNAL EOS XRootD redirector so this macro can
// run from outside the LPC interactive nodes, as long as XRootD access works.
TString redirector = "root://cmseos.fnal.gov/";
TString baseInDir  = "/eos/uscms/store/user/aahmad2/run2_files/fliprate_skims";
TString baseOutDir = "/eos/uscms/store/user/aahmad2/flip_rate_results";

// ----------------------------------------------------------------------------
// Return the list of input ROOT files for a given data-taking period.
//
// For Data, the measurement uses the SingleElectron datasets, or EGamma in
// 2018, because the charge-flip rate is measured from reconstructed electron
// pairs in Z->ee events.
//
// For MC, only the inclusive DY->ll sample is needed. DY provides a clean
// prompt Z->ee sample, so the same-sign pairs mostly come from detector charge
// mismeasurement rather than from physics backgrounds.
// ----------------------------------------------------------------------------
vector<TString> getSamples(string year, bool isData) {

    if (!isData) {
        return {Form("DYJetsToLLM50_%s.root", year.c_str())};
    }

    if (year == "2016preVFP") {
        return {
            "SingleElectronB_2016preVFP.root",
            "SingleElectronC_2016preVFP.root",
            "SingleElectronD_2016preVFP.root",
            "SingleElectronE_2016preVFP.root",
            "SingleElectronF_2016preVFP.root"
        };
    }

    if (year == "2016postVFP") {
        return {
            "SingleElectronF_2016postVFP.root",
            "SingleElectronG_2016postVFP.root",
            "SingleElectronH_2016postVFP.root"
        };
    }

    if (year == "2017") {
        return {
            "SingleElectronB_2017.root",
            "SingleElectronC_2017.root",
            "SingleElectronD_2017.root",
            "SingleElectronE_2017.root",
            "SingleElectronF_2017.root"
        };
    }

    if (year == "2018") {
        return {
            "EGammaA_2018.root",
            "EGammaB_2018.root",
            "EGammaC_2018.root",
            "EGammaD_2018.root"
        };
    }

    cout << "Unknown year: " << year << endl;
    return {};
}

// ----------------------------------------------------------------------------
// Measure the global electron charge-flip probability.
//
// The basic idea is:
//   1. Select clean Z->ee events.
//   2. Count opposite-sign electron pairs.
//   3. Count same-sign electron pairs.
//
// A true Z->ee decay is opposite-sign. Same-sign reconstructed pairs mostly
// appear when one electron charge is reconstructed incorrectly. The observed
// SS/OS ratio is then converted into a per-electron charge-flip probability.
// ----------------------------------------------------------------------------
double computeSimpleCF(string year, bool isData) {

    TString skimDir = Form("%s_%s", isData ? "Data" : "MC", year.c_str());
    vector<TString> samples = getSamples(year, isData);

    // One-bin histograms are used as event counters. We only need the total
    // number of same-sign and opposite-sign Z candidates for this unbinned rate.
    TH1D *h_ss = new TH1D("h_ss", "", 1, 0, 1);
    TH1D *h_os = new TH1D("h_os", "", 1, 0, 1);
    h_ss->Sumw2();
    h_os->Sumw2();

    // Loop over all samples for the requested year. Data has several run-period
    // files, while MC has one inclusive DY sample.
    for (auto &s : samples) {

        cout << " -> Processing sample: " << s;
        cout << (isData ? "  [DATA]" : "  [MC]") << endl;

        TString fileName = Form(
            "%s/%s/%s/%s",
            redirector.Data(),
            baseInDir.Data(),
            skimDir.Data(),
            s.Data()
        );

        TFile *fin = TFile::Open(fileName);
        if (!fin || fin->IsZombie()) {
            cout << "   !! Could not open file, skipping: " << fileName << endl;
            continue;
        }

        TTree *t = (TTree*)fin->Get("Events");
        if (!t) {
            cout << "   !! Events tree missing, skipping." << endl;
            fin->Close();
            continue;
        }

        MyBranch(t);
        Long64_t n = t->GetEntries();

        cout << "   Entries: " << n << endl;

        // These skims already contain events passing the tight object selection
        // used in the analysis. For this study, the _1 and _2 branches already
        // correspond to the two tight selected electrons, so no extra WP90 cut
        // is applied here.
        for (Long64_t i = 0; i < n; i++) {
            t->GetEntry(i);

            if (i % 1000000 == 0 && i > 0)
                cout << "     ... processed " << i << " events" << endl;

            // Keep only dielectron events. In this analysis framework, the ee
            // category peaks at cat = 40, and numberToCat(cat) maps it to "ee".
            if (std::string(numberToCat(cat)) != "ee") continue;

            // For MC, require both reconstructed electrons to be prompt
            // generator-level electrons. In the new skims this information is
            // stored in genPartFlav_*; value 1 corresponds to prompt electron.
            if (!isData && !(genPartFlav_1 == 1 && genPartFlav_2 == 1)) continue;

            // Select a clean Z->ee sample by requiring the dilepton invariant
            // mass to be within +/-10 GeV of the Z boson mass.
            double mll = (LepV(1) + LepV(2)).M();
            if (fabs(mll - 91.2) > 10) continue;

            // Opposite-sign pairs are the expected Z->ee topology. Same-sign
            // pairs are charge-flip candidates.
            if (q_1 * q_2 > 0) h_ss->Fill(0.5);
            else               h_os->Fill(0.5);
        }

        fin->Close();
        cout << "   Finished " << s << endl;
    }

    double nSS = h_ss->GetBinContent(1);
    double nOS = h_os->GetBinContent(1);

    // ------------------------------------------------------------------------
    // Convert the measured SS/OS ratio into the per-electron charge-flip rate.
    //
    // Let f be the probability that one electron has its charge mismeasured.
    // For a true opposite-sign Z->ee event:
    //
    //   P(reco SS) = 2f(1-f)
    //
    // because either electron can flip, but not both.
    //
    //   P(reco OS) = (1-f)^2 + f^2
    //
    // because either neither electron flips, or both electrons flip.
    //
    // Therefore:
    //
    //   R = N_SS / N_OS = 2f(1-f) / [(1-f)^2 + f^2]
    //
    // Solving this equation for f gives the expression below. The smaller
    // physical solution is chosen because charge-flip rates are expected to be
    // much smaller than 50%.
    // ------------------------------------------------------------------------
    double flipRate = 0.0;

    if (nOS > 0) {
        double R = nSS / nOS;
        double disc = 2 * R + 1 - R * R;

        if (disc >= 0 && R > 0)
            flipRate = ((R + 1) - sqrt(disc)) / (2 * R);
    }

    delete h_ss;
    delete h_os;

    return flipRate;
}

// ----------------------------------------------------------------------------
// Measure the charge-flip probability in MC and Data, then compute the global
// correction factor:
//
//        ChargeFlipSF = f(Data) / f(MC)
//
// The result is saved as one-bin histograms so downstream code can read the
// Data rate, MC rate, and scale factor directly from the output ROOT file.
// ----------------------------------------------------------------------------
void unbinnedChargeFlipSF(string year) {

    double fMC   = computeSimpleCF(year, false);
    double fData = computeSimpleCF(year, true);

    double SF = (fMC > 0 ? fData / fMC : 1.0);

    TString outName = Form(
        "%s/%s/ChargeFlipSF_SIMPLE_%s.root",
        redirector.Data(),
        baseOutDir.Data(),
        year.c_str()
    );

    TFile fout(outName, "RECREATE");

    TH1D *hMC   = new TH1D("MC_rate",      "MC global charge-flip rate", 1, 0, 1);
    TH1D *hData = new TH1D("Data_rate",    "Data global charge-flip rate", 1, 0, 1);
    TH1D *hSF   = new TH1D("ChargeFlipSF", "Data/MC scale factor", 1, 0, 1);

    hMC->SetBinContent(1, fMC);
    hData->SetBinContent(1, fData);
    hSF->SetBinContent(1, SF);

    hMC->Write();
    hData->Write();
    hSF->Write();

    fout.Close();

    cout << "Year " << year
         << "  f_MC = " << fMC
         << "  f_Data = " << fData
         << "  SF = " << SF
         << "\nSaved: " << outName << endl;
}

// ----------------------------------------------------------------------------
// Run the full measurement for all Run-2 eras.
//
// The 2016 sample is split into preVFP and postVFP because detector conditions
// changed between these periods, so separate charge-flip corrections are made.
// ----------------------------------------------------------------------------
void runChargeFlipSF_unbinned_allYears() {

    vector<string> years = {"2016preVFP", "2016postVFP", "2017", "2018"};

    for (auto &y : years) {
        unbinnedChargeFlipSF(y);
    }
}