// Author: Athar Ahmad (athar.ahmad@cern.ch)

#include "include/MyBranch.C"
#include "include/Kinematics.C"
#include "include/Xsections.C"
#include "filemap/FileMap.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLorentzVector.h"

// For random charge flipping
#include "TRandom3.h"
TRandom3 randCF(12345);

using namespace std;

// ============================================================================
// Data/MC stack for Z->ee charge-flip validation
//
// This macro is used to validate the electron charge-flip correction in the
// clean Z->ee control region. It makes stacked Data/MC plots of m_ee for:
//
//   1. same-sign ee events in the Z window,
//   2. opposite-sign ee events in the Z window,
//   3. sign-inclusive ee events in the Z window.
//
// The nominal histograms use the reconstructed charges directly from the skim.
// If fillChargeFlip=true, a second set of histograms is also filled in the same
// event loop after applying a stochastic MC-only charge migration.
//
// Run examples:
//
//   ChargeFlipRateTest(false);  // nominal only
//   ChargeFlipRateTest(true);   // nominal and charge-flip corrected
//
// Run2 is not read as a separate sample. Instead, the four era histograms are
// filled independently with their own luminosities and then added together.
// ============================================================================

static TString outTag = "nominal";

// The analysis is performed separately for each Run-2 era. The 2016 sample is
// split into preVFP and postVFP because the detector conditions changed.

static const vector<string> YEARS = {
    "2016preVFP",
    "2016postVFP",
    "2017",
    "2018"
};

// These are the three validation regions that will be plotted and written out.
// All three are restricted to |m_ee - 91.2| < 10 GeV (Z mass window) during the event loop.

static const vector<string> SCENARIOS = {
    "ss_Zwin",
    "os_Zwin",
    "inclusive_Zwin"
};

// Order in which MC backgrounds are stacked. This must match the process names
// used as keys in FileMap.h.

static const vector<string> BKG_ORDER = {
    "DY",
    "DY10_50",
    "VV",
    "ZZ",
    "TTbar",
    "ttV",
    "WJ",
    "ST",
    "VVV",
    "QCD",
    "other"
};

static map<string, int> fillColor = {
    {"DY",7},
    {"DY10_50",20},
    {"VV",8},
    {"ZZ",5},
    {"TTbar",46},
    {"ttV",4},
    {"WJ",9},
    {"ST",38},
    {"VVV",6},
    {"QCD", kOrange-3},
    {"other",3},
    {"ttH",28},
    {"signal",2}
};

// Change this if you want the output somewhere else.
static TString outDir =
    "/eos/uscms/store/user/aahmad2/flip_rate_stacks";

// Integrated luminosities in /pb.
// These are used only to normalize MC samples.
double getLumiPB(const string& year)
{
    if (year == "2016preVFP")  return 19500.0;
    if (year == "2016postVFP") return 16800.0;
    if (year == "2017")        return 41500.0;
    if (year == "2018")        return 59700.0;
    if (year == "Run2")        return 19500.0 + 16800.0 + 41500.0 + 59700.0;

    return 1.0;
}

string getLumiLabel(const string& year)
{
    if (year == "2016preVFP")  return "2016 preVFP, 19.5 fb^{-1}";
    if (year == "2016postVFP") return "2016 postVFP, 16.8 fb^{-1}";
    if (year == "2017")        return "2017, 41.5 fb^{-1}";
    if (year == "2018")        return "2018, 59.7 fb^{-1}";
    if (year == "Run2")        return "Run 2, 137.5 fb^{-1}";

    return year;
}

bool isElectronDataFile(const string& filename)
{
    return filename.find("SingleElectron") != string::npos ||
           filename.find("EGamma")         != string::npos;
}

TH1D* makeHist(const string& name)
{
    TH1D* h = new TH1D(
        name.c_str(),
        ";m_{ee} [GeV];Events",
        40,
        70.0,
        110.0
    );

    h->Sumw2();
    h->SetDirectory(nullptr);
    return h;
}

void addHist(TH1D*& target, TH1D* source, const string& newName)
{
    if (!source) return;

    if (!target) {
        target = (TH1D*)source->Clone(newName.c_str());
        target->SetDirectory(nullptr);
    } else {
        target->Add(source);
    }
}

// Convert each MC file into a luminosity-normalized sample weight.
// Equation for normalization is: (MC process cross section * Data Era Luminosity) / sum_of_generated_event_weights
//
// In this function normalization used is:
//
//      xsec / sum_of_generated_event_weights
//
// The luminosity factor is applied later in getEventWeight(). For Data, the
// weight is always 1.

double getMCWeight(const string& filename, TFile* f, bool isData)
{
    if (isData) return 1.0;

    TH1D* hnevts = (TH1D*)f->Get("hNWEvts");
    if (!hnevts) hnevts = (TH1D*)f->Get("hNEvts");

    double denom = hnevts ? hnevts->Integral() : 0.0;
    if (denom <= 0.0) {
        cout << "  [warning] Missing or empty hNWEvts/hNEvts for "
             << filename << ". Setting xsec weight to 0." << endl;
        return 0.0;
    }

    string base = gSystem->BaseName(filename.c_str());
    double xs = XSec(base);

    if (xs == 1.0) return 1.0;

    return xs / denom;
}

// Build the full event weight for MC events.
//
// The event weight includes:
//   - generator weight,
//   - branching/event weight from the skim,
//   - cross section normalization,
//   - luminosity,
//   - pileup weight,
//   - L1 prefiring weight,
//   - electron ID and isolation scale factors,
//   - trigger scale factor.
//
// Data events are not reweighted.

double getEventWeight(double xsecWeight, double lumiPB, bool isData)
{
    if (isData) return 1.0;

    double w = Generator_weight * brWeight * xsecWeight * lumiPB;

    w *= weightPUtruejson;
    w *= L1PreFiringWeight_Nom;

    // For the ee validation region, only the first two selected leptons are
    // used. The new skims already contain tight selected electrons in _1 and _2.
    w *= IDSF_1  * IDSF_2;
    w *= ISOSF_1 * ISOSF_2;

    // Trigger scale factors are applied in the same simple way as the main
    // histogram-making script for dilepton events.
    double trigSF = 1.0;

    if (isTrig_1 >= 1)       trigSF = TrigSF_1;
    else if (isTrig_1 == -1) trigSF = TrigSF_2;

    w *= trigSF;

    return w;
}

void drawCMSLabel(const string& year)
{
    TLatex cms;
    cms.SetNDC();
    cms.SetTextFont(61);
    cms.SetTextSize(0.050);
    cms.DrawLatex(0.12, 0.925, "CMS");

    TLatex prelim;
    prelim.SetNDC();
    prelim.SetTextFont(52);
    prelim.SetTextSize(0.045);
    prelim.DrawLatex(0.205, 0.925, "Preliminary");

    TLatex lumi;
    lumi.SetNDC();
    lumi.SetTextFont(42);
    lumi.SetTextSize(0.043);
    lumi.SetTextAlign(31);
    lumi.DrawLatex(0.68, 0.925, getLumiLabel(year).c_str());
}


// ============================================================================
// Charge-flip correction
//
// The charge-flip measurement produced one global flip probability in Data and
// one in MC for each era:
//
//      fData, fMC
//
// These are read from:
//      /eos/uscms/store/user/aahmad2/flip_rate_results/
//
// For the corrected MC histograms, the script does not simply scale the SS
// yield. Instead, it randomly migrates MC events between OS and SS according to
// the difference between fData and fMC. This is important because increasing the
// SS yield should also remove a small number of events from the OS yield.
// ============================================================================

struct ChargeFlipInfo {
    double fMC   = 0.0;
    double fData = 0.0;
    double SF    = 1.0;
};

ChargeFlipInfo loadChargeFlipInfo(const string& year)
{
    ChargeFlipInfo info;

    TString fileName = Form(
        "root://cmseos.fnal.gov//eos/uscms/store/user/aahmad2/flip_rate_results/ChargeFlipSF_SIMPLE_%s.root",
        year.c_str()
    );

    TFile* f = TFile::Open(fileName, "READ");

    if (!f || f->IsZombie()) {
        cout << "[warning] Could not open charge-flip file: "
             << fileName << endl;
        if (f) f->Close();
        return info;
    }

    TH1D* hMC   = (TH1D*)f->Get("MC_rate");
    TH1D* hData = (TH1D*)f->Get("Data_rate");
    TH1D* hSF   = (TH1D*)f->Get("ChargeFlipSF");

    if (hMC)   info.fMC   = hMC->GetBinContent(1);
    if (hData) info.fData = hData->GetBinContent(1);
    if (hSF)   info.SF    = hSF->GetBinContent(1);

    f->Close();

    cout << "[ChargeFlip] " << year
         << " fMC=" << info.fMC
         << " fData=" << info.fData
         << " SF=" << info.SF
         << endl;

    return info;
}

// Apply stochastic OS <-> SS migration to one MC ee event.
//
// For a true Z->ee event with per-electron flip probability f:
//
//      P(SS) = 2 f (1 - f)
//      P(OS) = (1 - f)^2 + f^2
//
// If Data has a larger SS probability than MC, some OS MC events are migrated
// to SS. If Data has a smaller SS probability than MC, some SS MC events are
// migrated back to OS.
//
// Only the charges are modified here. The event kinematics and weight are left
// unchanged.

void applyChargeFlipMigration(
    int& q1,
    int& q2,
    const ChargeFlipInfo& cf
) {
    double fMC   = cf.fMC;
    double fData = cf.fData;

    if (fMC <= 0.0 || fData <= 0.0) return;

    double PssMC   = 2.0 * fMC   * (1.0 - fMC);
    double PosMC   = (1.0 - fMC) * (1.0 - fMC) + fMC * fMC;

    double PssData = 2.0 * fData * (1.0 - fData);

    bool isOS = (q1 * q2 < 0);
    bool isSS = (q1 * q2 > 0);

    if (isOS && PssData > PssMC && PosMC > 0.0) {
        double pMove = (PssData - PssMC) / PosMC;

        if (randCF.Uniform() < pMove) {
            if (randCF.Uniform() < 0.5) q1 *= -1;
            else                        q2 *= -1;
        }
    }

    else if (isSS && PssData < PssMC && PssMC > 0.0) {
        double pMove = (PssMC - PssData) / PssMC;

        if (randCF.Uniform() < pMove) {
            if (randCF.Uniform() < 0.5) q1 *= -1;
            else                        q2 *= -1;
        }
    }
}

// Fill histograms for one data-taking era.
//
// This function reads all Data and MC files for a given year from FileMap.h.
// It fills the nominal histograms for every selected event. If fillChargeFlip
// is true, it also fills a second set of histograms after applying the MC-only
// charge-flip migration.
//
// Doing both nominal and corrected filling in the same event loop avoids reading
// the large skim files twice.

void fillYearHists(
    const string& year,
    map<string, map<string, TH1D*>>& HNominal,
    map<string, map<string, TH1D*>>& HChargeFlip,
    bool fillChargeFlip = true
) {
    cout << "\n==================================================" << endl;
    cout << "Filling histograms for year: " << year << endl;
    cout << "==================================================" << endl;

    map<string, vector<string>> fileMap = getFileMap(year);

    if (fileMap.empty()) {
        cout << "[error] FileMap is empty for year " << year << endl;
        return;
    }

    double lumiPB = getLumiPB(year);

    ChargeFlipInfo cfInfo;

    if (fillChargeFlip) {
        cfInfo = loadChargeFlipInfo(year);
    }

    for (const auto& procEntry : fileMap) {

        string proc = procEntry.first;
        bool isData = (proc == "data");

        for (const string& fileName : procEntry.second) {

            // The Data entry in FileMap contains both electron and muon datasets.
            // For this ee validation, use only SingleElectron/EGamma data to
            // avoid mixing in SingleMuon-triggered events.
            if (isData && !isElectronDataFile(fileName)) continue;

            cout << "\n[" << year << "] " << proc << " : " << fileName << endl;

            TFile* f = TFile::Open(fileName.c_str(), "READ");

            if (!f || f->IsZombie()) {
                cout << "  [skip] could not open file" << endl;
                if (f) f->Close();
                continue;
            }

            TTree* t = (TTree*)f->Get("Events");

            if (!t) {
                cout << "  [skip] missing Events tree" << endl;
                f->Close();
                continue;
            }

            // Speed optimization:
            // Only the branches needed for this validation plot are enabled. This greatly
            // reduces I/O compared to reading the full skim.

            t->SetBranchStatus("*", 0);

            const char* brs[] = {
                "cat",
                "pt_1", "pt_2",
                "eta_1", "eta_2",
                "phi_1", "phi_2",
                "m_1", "m_2",
                "q_1", "q_2",

                "Generator_weight",
                "brWeight",
                "weightPUtruejson",
                "L1PreFiringWeight_Nom",

                "IDSF_1", "IDSF_2",
                "ISOSF_1", "ISOSF_2",
                "TrigSF_1", "TrigSF_2",
                "isTrig_1"
            };

            for (auto br : brs) {
                if (t->GetBranch(br)) t->SetBranchStatus(br, 1);
            }

            MyBranch(t);

            double xsecWeight = getMCWeight(fileName, f, isData);

            Long64_t nEntries = t->GetEntries();
            Long64_t nPass = 0;

            for (Long64_t i = 0; i < nEntries; ++i) {
                t->GetEntry(i);

                if (i > 0 && i % 1000000 == 0) {
                    cout << "  processed " << i << " / " << nEntries << endl;
                }

                // Keep only the ee category. In this analysis framework cat == 40 corresponds
                // to selected dielectron events.
                if (cat != 40) continue;

                // Reconstruct m_ee from the two selected electrons. The skims already store
                // tight selected objects in the _1 and _2 branches, so no extra electron ID
                // requirement is applied here.
                TLorentzVector e1, e2;
                e1.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_1);
                e2.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_2);

                double mee = (e1 + e2).M();

                // The validation is focused on the Z peak region.
                if (fabs(mee - 91.2) > 10.0) continue;

                double w = getEventWeight(xsecWeight, lumiPB, isData);

                // Nominal filling:
                // use the reconstructed charges directly from the skim.
                bool isSS_nom = (q_1 * q_2 > 0);
                bool isOS_nom = (q_1 * q_2 < 0);

                if (isSS_nom) HNominal[proc]["ss_Zwin"]->Fill(mee, w);
                if (isOS_nom) HNominal[proc]["os_Zwin"]->Fill(mee, w);

                HNominal[proc]["inclusive_Zwin"]->Fill(mee, w);


                // Charge-flip corrected filling:
                // Data is copied unchanged. For MC, the reconstructed charges may be randomly
                // migrated between OS and SS according to the measured Data/MC charge-flip
                // difference.
                if (fillChargeFlip) {

                    int q1_corr = q_1;
                    int q2_corr = q_2;

                    if (!isData) {
                        applyChargeFlipMigration(q1_corr, q2_corr, cfInfo);
                    }

                    bool isSS_cf = (q1_corr * q2_corr > 0);
                    bool isOS_cf = (q1_corr * q2_corr < 0);

                    if (isSS_cf) HChargeFlip[proc]["ss_Zwin"]->Fill(mee, w);
                    if (isOS_cf) HChargeFlip[proc]["os_Zwin"]->Fill(mee, w);

                    HChargeFlip[proc]["inclusive_Zwin"]->Fill(mee, w);
                }

                ++nPass;
            }

            cout << "  passed ee Z-window selection: " << nPass << endl;

            f->Close();
        }
    }
}

// Write the per-process histograms, total MC histogram, Data histogram, and
// Data/MC ratio histogram to the output ROOT file.

void writeOutputHists(
    TFile* fout,
    const string& year,
    const string& scenario,
    map<string, map<string, TH1D*>>& H
) {
    fout->cd();

    TH1D* hData = nullptr;
    TH1D* hMC   = nullptr;

    if (H.count("data")) {
        hData = H["data"][scenario];
    }

    for (const string& proc : BKG_ORDER) {
        if (!H.count(proc)) continue;

        TH1D* h = H[proc][scenario];
        if (!h) continue;

        h->Write(Form("%s_%s_%s", year.c_str(), scenario.c_str(), proc.c_str()));
        addHist(hMC, h, Form("%s_%s_MCtotal", year.c_str(), scenario.c_str()));
    }

    if (hData) {
        hData->Write(Form("%s_%s_data", year.c_str(), scenario.c_str()));
    }

    if (hMC) {
        hMC->Write(Form("%s_%s_MCtotal", year.c_str(), scenario.c_str()));
    }

    if (hData && hMC) {
        TH1D* ratio = (TH1D*)hData->Clone(
            Form("%s_%s_ratio_data_over_mc", year.c_str(), scenario.c_str())
        );
        ratio->SetDirectory(nullptr);
        ratio->Divide(hMC);
        ratio->Write();
        delete ratio;
    }

    delete hMC;
}

// Draw one stacked Data/MC comparison plot.
//
// The canvas is split into three pads:
//   - upper-left: stacked MC and Data,
//   - lower-left: Data/Simulation ratio,
//   - right: legend with process yields.

void drawStack(
    const string& year,
    const string& scenario,
    map<string, map<string, TH1D*>>& H,
    TFile* fout
) {
    TH1D* hData = nullptr;

    if (H.count("data")) {
        hData = H["data"][scenario];
    }

    THStack* stack = new THStack(
        Form("stack_%s_%s", year.c_str(), scenario.c_str()),
        ""
    );

    TH1D* hMC = nullptr;

    for (const string& proc : BKG_ORDER) {
        if (!H.count(proc)) continue;

        TH1D* h = H[proc][scenario];
        if (!h || h->Integral() == 0.0) continue;

        h->SetFillColor(fillColor[proc]);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);

        stack->Add(h, "hist");
        addHist(hMC, h, Form("MCtotal_%s_%s", year.c_str(), scenario.c_str()));
    }

    if (!hMC) {
        cout << "[skip] no MC for " << year << " " << scenario << endl;
        delete stack;
        return;
    }

    double ymax = hMC->GetMaximum();
    if (hData) ymax = max(ymax, hData->GetMaximum());

    if (ymax <= 0.0) {
        cout << "[skip] empty plot for " << year << " " << scenario << endl;
        delete hMC;
        delete stack;
        return;
    }

    TCanvas* c = new TCanvas(
        Form("c_%s_%s", year.c_str(), scenario.c_str()),
        "",
        1000,
        800
    );

    TPad* p1 = new TPad("p1", "", 0.00, 0.32, 0.72, 1.00);
    TPad* p2 = new TPad("p2", "", 0.00, 0.00, 0.72, 0.32);
    TPad* p3 = new TPad("p3", "", 0.72, 0.00, 1.00, 1.00);

    p1->SetTopMargin(0.12);
    p1->SetBottomMargin(0.03);
    p1->SetLeftMargin(0.16);
    p1->SetRightMargin(0.04);
    p1->SetTicks(1, 1);
    p1->SetGridx(1);
    p1->SetGridy(1);

    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.35);
    p2->SetLeftMargin(0.16);
    p2->SetRightMargin(0.04);
    p2->SetTicks(1, 1);
    p2->SetGridx(1);
    p2->SetGridy(1);

    p3->SetLeftMargin(0.0);
    p3->SetRightMargin(0.05);
    p3->SetTopMargin(0.05);
    p3->SetBottomMargin(0.05);

    p1->Draw();
    p2->Draw();
    p3->Draw();

    p1->cd();

    stack->SetMaximum(1.45 * ymax);
    stack->SetMinimum(0.0);
    stack->Draw("hist");

    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleSize(0.052);
    stack->GetYaxis()->SetTitleOffset(1.5);
    stack->GetYaxis()->SetLabelSize(0.045);

    stack->GetXaxis()->SetTitle("");
    stack->GetXaxis()->SetLabelSize(0.0);
    stack->GetXaxis()->SetTitleSize(0.0);

    hMC->SetFillColor(kGray + 2);
    hMC->SetFillStyle(3004);
    hMC->SetMarkerSize(0);
    hMC->Draw("E2 SAME");

    if (hData) {
        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(1.0);
        hData->SetMarkerColor(kBlack);
        hData->SetLineColor(kBlack);
        hData->Draw("E SAME");
    }

    TLegend* leg = new TLegend(0.05, 0.05, 0.95, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.080);

    for (const string& proc : BKG_ORDER) {
        if (!H.count(proc)) continue;
        TH1D* h = H[proc][scenario];
        if (!h || h->Integral() == 0.0) continue;

        leg->AddEntry(
            h,
            Form("%s (%.1f)", proc.c_str(), h->Integral()),
            "f"
        );
    }

    if (hData) {
        leg->AddEntry(
            hData,
            Form("Data (%.1f)", hData->Integral()),
            "lep"
        );
    }

    leg->AddEntry(hMC, "MC stat. unc.", "f");

    p3->cd();
    leg->Draw();
    p1->cd();

    drawCMSLabel(year);

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.038);

    if (scenario == "ss_Zwin")
        label.DrawLatex(0.17, 0.82, "ee same-sign, |m_{ee}-91.2| < 10 GeV");
    else if (scenario == "os_Zwin")
        label.DrawLatex(0.17, 0.82, "ee opposite-sign, |m_{ee}-91.2| < 10 GeV");
    else
        label.DrawLatex(0.17, 0.82, "ee inclusive sign, |m_{ee}-91.2| < 10 GeV");

    p1->RedrawAxis();

    p2->cd();

    TH1D* ratio = nullptr;

    if (hData) {
        ratio = (TH1D*)hData->Clone(
            Form("ratio_%s_%s", year.c_str(), scenario.c_str())
        );
        ratio->SetDirectory(nullptr);
        ratio->Divide(hMC);
    } else {
        ratio = (TH1D*)hMC->Clone(
            Form("ratio_%s_%s", year.c_str(), scenario.c_str())
        );
        ratio->SetDirectory(nullptr);
        ratio->Reset();
    }

    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("Data/Simulation");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.100);
    ratio->GetYaxis()->SetTitleOffset(0.50);
    ratio->GetYaxis()->SetLabelSize(0.080);

    ratio->GetXaxis()->SetTitle("m_{ee} [GeV]");
    ratio->GetXaxis()->SetTitleSize(0.120);
    ratio->GetXaxis()->SetTitleOffset(1.05);
    ratio->GetXaxis()->SetLabelSize(0.090);

    ratio->SetMinimum(0.0);
    ratio->SetMaximum(2.0);
    ratio->SetMarkerStyle(20);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);

    ratio->Draw("E");

    TLine* line1 = new TLine(70.0, 1.0, 110.0, 1.0);
    line1->SetLineStyle(2);
    line1->Draw("SAME");

    TLine* line05 = new TLine(70.0, 0.5, 110.0, 0.5);
    line05->SetLineStyle(3);
    line05->Draw("SAME");

    TLine* line15 = new TLine(70.0, 1.5, 110.0, 1.5);
    line15->SetLineStyle(3);
    line15->Draw("SAME");

    p2->SetGridy();
    p2->RedrawAxis();

    TString plotDir = Form("%s/%s/%s", outDir.Data(), outTag.Data(), year.c_str());
    gSystem->mkdir(plotDir, kTRUE);

    TString pngName = Form(
        "%s/mEE_%s_%s.png",
        plotDir.Data(),
        year.c_str(),
        scenario.c_str()
    );

    c->SaveAs(pngName);

    writeOutputHists(fout, year, scenario, H);

    delete line1;
    delete line05;
    delete line15;
    delete ratio;
    delete leg;
    delete p1;
    delete p2;
    delete p3;
    delete c;
    delete hMC;
    delete stack;
}

void initializeHists(map<string, map<string, TH1D*>>& H,
                     const vector<string>& processes,
                     const string& year)
{
    for (const string& proc : processes) {
        for (const string& scenario : SCENARIOS) {
            H[proc][scenario] = makeHist(
                Form("h_%s_%s_%s", year.c_str(), scenario.c_str(), proc.c_str())
            );
        }
    }
}

// Add one era into the Run2 histogram collection.
//
// Run2 is produced by summing already-normalized era histograms. No Run2 files
// are read directly, which avoids double counting.

void addYearToRun2(
    map<string, map<string, TH1D*>>& run2,
    map<string, map<string, TH1D*>>& yearHists
) {
    for (auto& procEntry : yearHists) {
        const string& proc = procEntry.first;

        for (const string& scenario : SCENARIOS) {
            if (!procEntry.second.count(scenario)) continue;

            addHist(
                run2[proc][scenario],
                procEntry.second[scenario],
                Form("h_Run2_%s_%s", scenario.c_str(), proc.c_str())
            );
        }
    }
}

void deleteHists(map<string, map<string, TH1D*>>& H)
{
    for (auto& procEntry : H) {
        for (auto& histEntry : procEntry.second) {
            delete histEntry.second;
        }
    }

    H.clear();
}


// Main function.
//
// If fillChargeFlip=false:
//   - only nominal histograms and plots are produced.
//
// If fillChargeFlip=true:
//   - nominal and charge-flip-corrected histograms are both produced in the
//     same pass over the skims.
//
// Outputs:
//   - Zmass_DataMCStacks_nominal.root
//   - Zmass_DataMCStacks_withChargeFlip.root, if requested (if bool is true)
//   - PNG plots under outDir/nominal and outDir/withChargeFlip

void ChargeFlipRateTest(bool fillChargeFlip = true)
{

    cout << "\n==================================================" << endl;
    cout << "Starting ChargeFlipRateTest" << endl;
    cout << "This script will take quite a while to run." << endl;
    cout << "It loops over all Run-2 Data and MC files and writes both ROOT histograms and PNG plots." << endl;
    cout << "fillChargeFlip = " << (fillChargeFlip ? "true" : "false") << endl;
    cout << "Output directory: " << outDir << endl;
    cout << "==================================================\n" << endl;

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    gSystem->mkdir(outDir, kTRUE);

    TString outRootNominal = Form(
        "%s/Zmass_DataMCStacks_nominal.root",
        outDir.Data()
    );

    TString outRootCF = Form(
        "%s/Zmass_DataMCStacks_withChargeFlip.root",
        outDir.Data()
    );

    TFile* foutNominal = new TFile(outRootNominal, "RECREATE");
    TFile* foutCF = nullptr;

    if (fillChargeFlip) {
        foutCF = new TFile(outRootCF, "RECREATE");
    }

    vector<string> processes = BKG_ORDER;
    processes.push_back("data");

    map<string, map<string, TH1D*>> HRun2Nominal;
    map<string, map<string, TH1D*>> HRun2CF;

    for (const string& year : YEARS) {

        map<string, map<string, TH1D*>> HYearNominal;
        map<string, map<string, TH1D*>> HYearCF;

        initializeHists(HYearNominal, processes, year);

        if (fillChargeFlip) {
            initializeHists(HYearCF, processes, year);
        }

        fillYearHists(year, HYearNominal, HYearCF, fillChargeFlip);

        outTag = "nominal";
        for (const string& scenario : SCENARIOS) {
            drawStack(year, scenario, HYearNominal, foutNominal);
        }

        addYearToRun2(HRun2Nominal, HYearNominal);

        if (fillChargeFlip) {
            outTag = "withChargeFlip";

            for (const string& scenario : SCENARIOS) {
                drawStack(year, scenario, HYearCF, foutCF);
            }

            addYearToRun2(HRun2CF, HYearCF);
        }

        deleteHists(HYearNominal);

        if (fillChargeFlip) {
            deleteHists(HYearCF);
        }
    }

    outTag = "nominal";
    for (const string& scenario : SCENARIOS) {
        drawStack("Run2", scenario, HRun2Nominal, foutNominal);
    }

    if (fillChargeFlip) {
        outTag = "withChargeFlip";

        for (const string& scenario : SCENARIOS) {
            drawStack("Run2", scenario, HRun2CF, foutCF);
        }
    }

    deleteHists(HRun2Nominal);

    if (fillChargeFlip) {
        deleteHists(HRun2CF);
    }

    foutNominal->Close();

    if (foutCF) {
        foutCF->Close();
    }

    cout << "\nSaved nominal ROOT output: " << outRootNominal << endl;

    if (fillChargeFlip) {
        cout << "Saved charge-flip ROOT output: " << outRootCF << endl;
    }

    cout << "Saved PNG plots under: " << outDir << endl;
}