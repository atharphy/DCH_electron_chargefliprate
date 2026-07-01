#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
TRandom3 randCF(12345);


#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "include/MyBranch.C"
#include "include/Kinematics.C"
#include "include/MET_split.C"
#include "include/Xsections.C"

using std::string;
using std::vector;
using std::unordered_map;
using FRMap = std::unordered_map<std::string, std::vector<TH1D*>>;

const vector<string> channels = {
    "ee","em","et","mm","mt","tt",
    "eee","eem","eet","eme","emm","emt","ete","etm","ett",
    "mme","mmm","mmt","mte","mtm","mtt","tte","ttm","ttt",
    "eeee","eeem","eeet","eemm","eemt","eett",
    "emem","emet","emmm","emmt","emtt",
    "etet","etmm","etmt","ettt",
    "mmmm","mmmt","mmtt","mtmt","mttt","tttt"
};

// const vector<string> channels = {
//     "eet","mmt","eett","mmtt","ett","mtt","ettt","mttt"
// };

const int NlepMax = 4;

string cuts = "tight";
string year = "2016";
// bool useTightTauForMC = false;
enum TauWP { TIGHT, LOOSE, LOOSE_NOT_TIGHT };
bool doExtra2D = false; 

bool passLeptonCuts(char flavor, int idx) {
    if (cuts == "tight") {
        if (flavor == 'e') {
            if ((idx==1 && !EleID_WP90_1) ||
                (idx==2 && !EleID_WP90_2) ||
                (idx==3 && !EleID_WP90_3) ||
                (idx==4 && !EleID_WP90_4)) return false;
            if ((idx==1 && iso_1 >= 0.1) ||
                (idx==2 && iso_2 >= 0.1) ||
                (idx==3 && iso_3 >= 0.1) ||
                (idx==4 && iso_4 >= 0.1)) return false;
        } else if (flavor == 'm') {
            if ((idx==1 && MuID_1 != 3) ||
                (idx==2 && MuID_2 != 3) ||
                (idx==3 && MuID_3 != 3) ||
                (idx==4 && MuID_4 != 3)) return false;
            if ((idx==1 && iso_1 >= 0.15) ||
                (idx==2 && iso_2 >= 0.15) ||
                (idx==3 && iso_3 >= 0.15) ||
                (idx==4 && iso_4 >= 0.15)) return false;
        } else if (flavor == 't') {
            if ((idx==1 && !(TauIDe_1 > 8 && TauIDm_1 > 4 && TauIDj_1 > 16)) ||
                (idx==2 && !(TauIDe_2 > 8 && TauIDm_2 > 4 && TauIDj_2 > 16)) ||
                (idx==3 && !(TauIDe_3 > 8 && TauIDm_3 > 4 && TauIDj_3 > 16)) ||
                (idx==4 && !(TauIDe_4 > 8 && TauIDm_4 > 4 && TauIDj_4 > 16))) return false;
        }
    } else if (cuts == "loose") {
        if (flavor == 'e') {
            if ((idx==1 && iso_1 >= 0.1) ||
                (idx==2 && iso_2 >= 0.1) ||
                (idx==3 && iso_3 >= 0.1) ||
                (idx==4 && iso_4 >= 0.1)) return false;
        } else if (flavor == 'm') {
            if ((idx==1 && iso_1 >= 0.5) ||
                (idx==2 && iso_2 >= 0.5) ||
                (idx==3 && iso_3 >= 0.5) ||
                (idx==4 && iso_4 >= 0.5)) return false;
        } else if (flavor == 't') {
            if ((idx==1 && !(TauIDe_1 > 8 && TauIDm_1 > 4 && TauIDj_1 > 2)) ||
                (idx==2 && !(TauIDe_2 > 8 && TauIDm_2 > 4 && TauIDj_2 > 2)) ||
                (idx==3 && !(TauIDe_3 > 8 && TauIDm_3 > 4 && TauIDj_3 > 2)) ||
                (idx==4 && !(TauIDe_4 > 8 && TauIDm_4 > 4 && TauIDj_4 > 2))) return false;
        }
    }
    return true;
}
bool passLooseTau(int idx) {
    if (idx == 1) return (TauIDj_1 > 2 && TauIDm_1 > 4 && TauIDe_1 > 8);
    if (idx == 2) return (TauIDj_2 > 2 && TauIDm_2 > 4 && TauIDe_2 > 8);
    if (idx == 3) return (TauIDj_3 > 2 && TauIDm_3 > 4 && TauIDe_3 > 8);
    if (idx == 4) return (TauIDj_4 > 2 && TauIDm_4 > 4 && TauIDe_4 > 8);
    return false;
}

bool passTightTau(int idx) {
    if (idx == 1) return (TauIDj_1 > 16 && TauIDm_1 > 4 && TauIDe_1 > 8);
    if (idx == 2) return (TauIDj_2 > 16 && TauIDm_2 > 4 && TauIDe_2 > 8);
    if (idx == 3) return (TauIDj_3 > 16 && TauIDm_3 > 4 && TauIDe_3 > 8);
    if (idx == 4) return (TauIDj_4 > 16 && TauIDm_4 > 4 && TauIDe_4 > 8);
    return false;
}

bool passLooseElectron(int idx) {
    if (idx == 1) return (iso_1 < 0.10);
    if (idx == 2) return (iso_2 < 0.10);
    if (idx == 3) return (iso_3 < 0.10);
    if (idx == 4) return (iso_4 < 0.10);
    return false;
}

bool passTightElectron(int idx) {
    if (idx == 1) return (EleID_WP90_1 == 1 && iso_1 < 0.10);
    if (idx == 2) return (EleID_WP90_2 == 1 && iso_2 < 0.10);
    if (idx == 3) return (EleID_WP90_3 == 1 && iso_3 < 0.10);
    if (idx == 4) return (EleID_WP90_4 == 1 && iso_4 < 0.10);
    return false;
}

bool passLooseMuon(int idx) {
    if (idx == 1) return (iso_1 < 0.5);
    if (idx == 2) return (iso_2 < 0.5);
    if (idx == 3) return (iso_3 < 0.5);
    if (idx == 4) return (iso_4 < 0.5);
    return false;
}

bool passTightMuon(int idx) {
    if (idx == 1) return (MuID_1 == 3 && iso_1 < 0.15);
    if (idx == 2) return (MuID_2 == 3 && iso_2 < 0.15);
    if (idx == 3) return (MuID_3 == 3 && iso_3 < 0.15);
    if (idx == 4) return (MuID_4 == 3 && iso_4 < 0.15);
    return false;
}

bool passLooseNotTightTau(int idx) {
    return (passLooseTau(idx) && !passTightTau(idx));
}

bool passLooseNotTightElectron(int idx) {
    return passLooseElectron(idx) && !passTightElectron(idx);
}

bool passLooseNotTightMuon(int idx) {
    return passLooseMuon(idx) && !passTightMuon(idx);
}

bool passEventCuts(const string& catstr,
                   bool isData,
                   bool isDY,
                   bool isTT,
                   bool isWJ,
                   bool isQCD,
                   TauWP tauWP_MC = LOOSE) {

    auto passTau = [&](int idx){
        switch (tauWP_MC) {
            case TIGHT: return passTightTau(idx);
            case LOOSE: return passLooseTau(idx);
            case LOOSE_NOT_TIGHT: return passLooseNotTightTau(idx);
        }
        return false;
    };

    if (isData && 
        (catstr == "eet"  || catstr == "mmt" ||
         catstr == "eett" || catstr == "mmtt")) {

        bool firstTight  = passLeptonCuts(catstr[0], 1);
        bool secondTight = passLeptonCuts(catstr[1], 2);

        if (!firstTight || !secondTight)
            return false;

        if (q_1 * q_2 >= 0)
            return false;

        if (catstr == "eet" || catstr == "mmt") {

            if (!passTightTau(3))
                return false;

            return true;
        }

        if (catstr == "eett" || catstr == "mmtt") {

            if (!passTightTau(3))
                return false;

            if (!passTightTau(4))
                return false;

            return true;
        }
    }

    if (!isData &&
        (catstr == "eet"  || catstr == "mmt" ||
         catstr == "eett" || catstr == "mmtt")) {

        bool firstTight  = passLeptonCuts(catstr[0], 1);
        bool secondTight = passLeptonCuts(catstr[1], 2);

        if (!firstTight || !secondTight)
            return false;

        if (q_1 * q_2 >= 0)
            return false;

        if (catstr == "eet" || catstr == "mmt") {

            if (!passTau(3))
                return false;

            return true;
        }

        if (catstr == "eett" || catstr == "mmtt") {

            if (!passTau(3) || !passTau(4))
                return false;

            return true;
        }
    }

    if (!isData &&
        (catstr == "ett"  || catstr == "mtt" ||
         catstr == "ettt" || catstr == "mttt")) {

        bool firstTight = passLeptonCuts(catstr[0], 1);
        if (!firstTight)
            return false;

        for (size_t i = 1; i < catstr.size(); ++i) {

            if (catstr[i] != 't')
                return false;

            if (!passTau(i + 1))
                return false;
        }

        return true;
    }

    for (size_t i = 0; i < catstr.size(); ++i)
        if (!passLeptonCuts(catstr[i], i + 1))
            return false;

    return true;
}


void createHists(unordered_map<string,TH1D*>& h, const string& name,
                 const string& title, int n, double xmin, double xmax) {
    h.reserve(channels.size() * 3);
    for (auto& ch : channels) {
        string base = "h_" + name + "_" + ch;
        h[ch] = new TH1D(base.c_str(), (title + " " + ch).c_str(), n, xmin, xmax);
        h[ch]->Sumw2();

        if (ch=="ee"||ch=="em"||ch=="et"||ch=="mm"||ch=="mt"||ch=="tt") {
            h[ch+"_os"] = new TH1D((base+"_os").c_str(), (title+" "+ch+" OS").c_str(), n, xmin, xmax);
            h[ch+"_os"]->Sumw2();
            h[ch+"_ss"] = new TH1D((base+"_ss").c_str(), (title+" "+ch+" SS").c_str(), n, xmin, xmax);
            h[ch+"_ss"]->Sumw2();
            h[ch+"_os_Zwin"] = new TH1D((base+"_os_Zwin").c_str(), (title+" "+ch+" OS Zwin").c_str(), n, xmin, xmax);
            h[ch+"_os_Zwin"]->Sumw2();
            h[ch+"_os_Zveto"] = new TH1D((base+"_os_Zveto").c_str(), (title+" "+ch+" OS Zveto").c_str(), n, xmin, xmax);
            h[ch+"_os_Zveto"]->Sumw2();
            h[ch+"_ss_Zwin"] = new TH1D((base+"_ss_Zwin").c_str(), (title+" "+ch+" SS Zwin").c_str(), n, xmin, xmax);
            h[ch+"_ss_Zwin"]->Sumw2();
            h[ch+"_ss_Zveto"] = new TH1D((base+"_ss_Zveto").c_str(), (title+" "+ch+" SS Zveto").c_str(), n, xmin, xmax);
            h[ch+"_ss_Zveto"]->Sumw2();
        } else {
            h[ch+"_Zwin"] = new TH1D((base+"_Zwin").c_str(), (title+" "+ch+" Zwin").c_str(), n, xmin, xmax);
            h[ch+"_Zwin"]->Sumw2();
            h[ch+"_Zveto"] = new TH1D((base+"_Zveto").c_str(), (title+" "+ch+" Zveto").c_str(), n, xmin, xmax);
            h[ch+"_Zveto"]->Sumw2();
        }
    }
}

void createHists2D(unordered_map<string,TH2D*>& h,
                   const string& name,
                   const string& title,
                   int nx,double xmin,double xmax,
                   int ny,double ymin,double ymax)
{
    h.reserve(channels.size() * 6);

    for (auto& ch : channels) {

        string base = "h_" + name + "_" + ch;

        h[ch] = new TH2D(base.c_str(),
                         (title + " " + ch).c_str(),
                         nx,xmin,xmax,
                         ny,ymin,ymax);
        h[ch]->Sumw2();

        if (ch=="ee"||ch=="em"||ch=="et"||
            ch=="mm"||ch=="mt"||ch=="tt") {

            h[ch+"_os"] = new TH2D((base+"_os").c_str(),
                                   (title+" "+ch+" OS").c_str(),
                                   nx,xmin,xmax,
                                   ny,ymin,ymax);
            h[ch+"_os"]->Sumw2();

            h[ch+"_ss"] = new TH2D((base+"_ss").c_str(),
                                   (title+" "+ch+" SS").c_str(),
                                   nx,xmin,xmax,
                                   ny,ymin,ymax);
            h[ch+"_ss"]->Sumw2();

            h[ch+"_os_Zwin"] = new TH2D((base+"_os_Zwin").c_str(),
                                        (title+" "+ch+" OS Zwin").c_str(),
                                        nx,xmin,xmax,
                                        ny,ymin,ymax);
            h[ch+"_os_Zwin"]->Sumw2();

            h[ch+"_os_Zveto"] = new TH2D((base+"_os_Zveto").c_str(),
                                         (title+" "+ch+" OS Zveto").c_str(),
                                         nx,xmin,xmax,
                                         ny,ymin,ymax);
            h[ch+"_os_Zveto"]->Sumw2();

            h[ch+"_ss_Zwin"] = new TH2D((base+"_ss_Zwin").c_str(),
                                        (title+" "+ch+" SS Zwin").c_str(),
                                        nx,xmin,xmax,
                                        ny,ymin,ymax);
            h[ch+"_ss_Zwin"]->Sumw2();

            h[ch+"_ss_Zveto"] = new TH2D((base+"_ss_Zveto").c_str(),
                                         (title+" "+ch+" SS Zveto").c_str(),
                                         nx,xmin,xmax,
                                         ny,ymin,ymax);
            h[ch+"_ss_Zveto"]->Sumw2();
        }
        else {

            h[ch+"_Zwin"] = new TH2D((base+"_Zwin").c_str(),
                                     (title+" "+ch+" Zwin").c_str(),
                                     nx,xmin,xmax,
                                     ny,ymin,ymax);
            h[ch+"_Zwin"]->Sumw2();

            h[ch+"_Zveto"] = new TH2D((base+"_Zveto").c_str(),
                                      (title+" "+ch+" Zveto").c_str(),
                                      nx,xmin,xmax,
                                      ny,ymin,ymax);
            h[ch+"_Zveto"]->Sumw2();
        }
    }
}

inline void fillOne(unordered_map<string,TH1D*>& h, const string& key, double val, double w) {
    auto it = h.find(key);
    if (it != h.end()) it->second->Fill(val, w);
}

inline void fillDilep(unordered_map<string,TH1D*>& h, 
                      const string& cat, double val, double w, 
                      int q1, int q2, const std::string& Zflag = "")
{
    fillOne(h, cat, val, w);

    bool isOS = (q1 * q2 < 0);
    bool isSS = (q1 * q2 > 0);
    if (isOS) fillOne(h, cat + "_os", val, w);
    if (isSS) fillOne(h, cat + "_ss", val, w);
    if (!Zflag.empty()) {
        if (isOS) fillOne(h, cat + (Zflag == "Zwindow" ? "_os_Zwin" : "_os_Zveto"), val, w);
        if (isSS) fillOne(h, cat + (Zflag == "Zwindow" ? "_ss_Zwin" : "_ss_Zveto"), val, w);
    }
}

inline void fillMultilep(unordered_map<string,TH1D*>& h, const string& cat, double val, double w, const string& flag) {
    fillOne(h, cat, val, w);
    if (flag == "Zwindow") fillOne(h, cat + "_Zwin", val, w);
    else if (flag == "Zveto" || flag == "Zv") fillOne(h, cat + "_Zveto", val, w);
}

inline void fillDilep2D(unordered_map<string,TH2D*>& h,
                        const string& cat,
                        double x,double y,double w,
                        int q1,int q2,
                        const string& Zflag="")
{
    auto fillSafe = [&](const string& key){
        auto it = h.find(key);
        if(it != h.end()) it->second->Fill(x,y,w);
    };

    fillSafe(cat);

    bool isOS = (q1*q2 < 0);
    bool isSS = (q1*q2 > 0);

    if(isOS) fillSafe(cat+"_os");
    if(isSS) fillSafe(cat+"_ss");

    if(!Zflag.empty()) {
        if(isOS)
            fillSafe(cat + (Zflag=="Zwindow" ? "_os_Zwin" : "_os_Zveto"));
        if(isSS)
            fillSafe(cat + (Zflag=="Zwindow" ? "_ss_Zwin" : "_ss_Zveto"));
    }
}

inline void fillMultilep2D(unordered_map<string,TH2D*>& h,
                           const string& cat,
                           double x,double y,double w,
                           const string& flag)
{
    auto fillSafe = [&](const string& key){
        auto it = h.find(key);
        if(it != h.end()) it->second->Fill(x,y,w);
    };

    fillSafe(cat);

    if(flag=="Zwindow")
        fillSafe(cat+"_Zwin");
    else if(flag=="Zveto" || flag=="Zv")
        fillSafe(cat+"_Zveto");
}

void scaleAndWrite(unordered_map<string,TH1D*>& h) {
    for (auto& [_, hist] : h) hist->Write();
}

void write2D(unordered_map<string,TH2D*>& h) {
    for (auto& [_, hist] : h)
        hist->Write();
}

double clampPt(double pt) {
    if (pt < 20) return 20.01;
    if (pt > 200) return 199.9;
    return pt;
}

bool applyChargeFlipStochastic(vector<Lepton>& L, const string& catstr,
                               TH2D* hProbMC, TH2D* hSF)
{
    if (catstr != "ee" || L.size() < 2) return false;

    auto f_lookup = [&](double pt, double eta) {
        double ptc  = std::min(pt, hProbMC->GetXaxis()->GetXmax() - 1e-3);
        double aeta = fabs(eta);
        int binMC   = hProbMC->FindBin(ptc, aeta);
        int binSF   = hSF->FindBin(ptc, aeta);
        double fMC  = std::clamp(hProbMC->GetBinContent(binMC), 1e-6, 0.5);
        double SF   = std::clamp(hSF->GetBinContent(binSF), 0.5, 2.0);
        double fData = std::clamp(fMC * SF, 1e-6, 0.5);
        return std::make_pair(fMC, fData);
    };

    auto [fMC1, fData1] = f_lookup(L[0].pt, L[0].eta);
    auto [fMC2, fData2] = f_lookup(L[1].pt, L[1].eta);

    double Padd1    = std::max(0.0, (fData1 - fMC1) / (1.0 - fMC1));
    double Padd2    = std::max(0.0, (fData2 - fMC2) / (1.0 - fMC2));
    double Punflip1 = (fData1 < fMC1) ? (fMC1 - fData1) / fMC1 : 0.0;
    double Punflip2 = (fData2 < fMC2) ? (fMC2 - fData2) / fMC2 : 0.0;

    bool isOS = (L[0].q * L[1].q < 0);
    bool isSS = (L[0].q * L[1].q > 0);
    bool flipped = false;

    if (isOS) {
        double Ppair = Padd1*(1 - Padd2) + Padd2*(1 - Padd1);
        if (randCF.Uniform() < Ppair) {
            double w1 = Padd1*(1 - Padd2);
            double w2 = Padd2*(1 - Padd1);
            double r2 = randCF.Uniform();
            if (r2 < w1 / (w1 + w2)) L[0].q *= -1;
            else                     L[1].q *= -1;
            flipped = true;
        }
    }
    else if (isSS) {
        double Ppair = Punflip1*(1 - Punflip2) + Punflip2*(1 - Punflip1);
        if (randCF.Uniform() < Ppair) {
            double w1 = Punflip1*(1 - Punflip2);
            double w2 = Punflip2*(1 - Punflip1);
            double r2 = randCF.Uniform();
            if (r2 < w1 / (w1 + w2)) L[0].q *= -1;
            else                     L[1].q *= -1;
            flipped = true;
        }
    }

    return flipped;
}
double getTauFakeRate(double pt, double eta,
                      const std::string& year,
                      const std::string& channel,
                      FRMap& FR_MC)
{
    std::string key = year + "_" + channel;

    auto it = FR_MC.find(key);
    if (it == FR_MC.end()) {
        std::cerr << "ERROR: FR not found for " << key << std::endl;
        return 0.0;
    }

    int idx = (std::fabs(eta) < 1.5) ? 0 : 1;
    double ptc = clampPt(pt);

    TH1D* h = it->second[idx];
    int bin = std::clamp(h->FindBin(ptc), 1, h->GetNbinsX());

    return std::clamp(h->GetBinContent(bin), 0.0, 0.99);
}

void loadTauFR(FRMap& FR_MC,
               const std::string& year,
               const std::string& channel)
{
    TString fileName = Form(
        "/Volumes/Lexar/DCH_analysis/rederived_fakes/"
        "Tau_DY_Fake_rates_loose_not_tight/DY_tau_fake_rates_%s.root",
        channel.c_str()
    );

    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open FR file "
                  << fileName << std::endl;
        exit(1);
    }

    TString hBarrelName = Form(
        "DY_FR_%s_%s_eta0to1p5",
        channel.c_str(), year.c_str()
    );
    TString hEndcapName = Form(
        "DY_FR_%s_%s_eta1p5to2p5",
        channel.c_str(), year.c_str()
    );

    TH1D* hBarrel = (TH1D*)f->Get(hBarrelName);
    TH1D* hEndcap = (TH1D*)f->Get(hEndcapName);

    if (!hBarrel || !hEndcap) {
        std::cerr << "ERROR: missing FR histograms:\n"
                  << "  " << hBarrelName << "\n"
                  << "  " << hEndcapName << std::endl;
        exit(1);
    }

    std::string key = year + "_" + channel;
    FR_MC[key] = {
        (TH1D*)hBarrel->Clone(),
        (TH1D*)hEndcap->Clone()
    };

    for (auto* h : FR_MC[key])
        h->SetDirectory(0);

    f->Close();

    std::cout << "[FR][MC] loaded " << key << std::endl;
}

void loadTauFR_Data(FRMap& FR_Data,
                    const std::string& year,
                    const std::string& channel)
{
    TString fileName = Form(
        "/Volumes/Lexar/DCH_analysis/rederived_fakes/Final_tau_fake_rates/Data_tau_fake_rates_%s_only.root",
        channel.c_str()
    );

    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open DATA FR file "
                  << fileName << std::endl;
        exit(1);
    }

    TString hBarrelName = Form(
        "Data_FR_%s_only_%s_eta0to1p5",
        channel.c_str(), year.c_str()
    );
    TString hEndcapName = Form(
        "Data_FR_%s_only_%s_eta1p5to2p5",
        channel.c_str(), year.c_str()
    );

    TH1D* hBarrel = (TH1D*)f->Get(hBarrelName);
    TH1D* hEndcap = (TH1D*)f->Get(hEndcapName);

    if (!hBarrel || !hEndcap) {
        std::cerr << "ERROR: missing DATA FR histograms:\n"
                  << "  " << hBarrelName << "\n"
                  << "  " << hEndcapName << std::endl;
        exit(1);
    }

    std::string key = year + "_" + channel;
    FR_Data[key] = {
        (TH1D*)hBarrel->Clone(),
        (TH1D*)hEndcap->Clone()
    };

    for (auto* h : FR_Data[key])
        h->SetDirectory(0);

    f->Close();

    std::cout << "[FR][DATA] loaded " << key << std::endl;
}

void DCH_run2_tauFR(string inYear="2016", string inCuts="tight", const char* ext="root") {
    (void)ext;
    
    year = inYear;
    cuts = inCuts;

    FRMap FR_MC;
    FRMap FR_Data;

    loadTauFR(FR_MC, year, "eet");
    loadTauFR(FR_MC, year, "mmt");

    loadTauFR_Data(FR_Data, year, "eet");
    loadTauFR_Data(FR_Data, year, "mmt");

    TFile *fMC = TFile::Open(Form("/Volumes/Lexar/DCH_analysis/guru_offline/run2_chargeflip_updated/ChargeFlipRate_MC_%s.root", year.c_str()));
    TH2D *hProbMC = (TH2D*)fMC->Get("h_rate");

    TFile *fSF = TFile::Open(Form("/Volumes/Lexar/DCH_analysis/guru_offline/run2_chargeflip_updated/ChargeFlipSF_%s.root", year.c_str()));
    TH2D *hSF = (TH2D*)fSF->Get("ChargeFlipSF");

    TH2D *hMatrix = (TH2D*)fMC->Get("MatrixCorrection");
    if (!hMatrix) std::cout << "[Warning] MatrixCorrection missing for " << year << std::endl;

    const char* inDir = "/Volumes/Lexar/DCH_analysis/guru_offline/inputs/run2_skims";
    TString outDirBase = "/Volumes/Lexar/DCH_analysis/guru_offline/run2_hists_Tau_Fake_weights";
    TString outDir = Form("%s/%s/%s", outDirBase.Data(), year.c_str(), cuts.c_str());
    gSystem->mkdir(outDir, kTRUE);

    TSystemDirectory dir("indir", inDir);
    TList* fileList = dir.GetListOfFiles();
    vector<string> files;
    if (fileList) {
        TIter next(fileList);
        TObject* obj;
        while ((obj = next())) {
            string fname = obj->GetName();
            if (fname.rfind(".root") == string::npos) continue;
            if (fname.rfind("._", 0) == 0) continue; 
            if (year != "run2" && fname.find(year) == string::npos) continue;
            files.emplace_back(Form("%s/%s", inDir, fname.c_str()));
        }
    }

    double lumi = (year=="2016preVFP")?19500:(year=="2016postVFP")?16400:(year=="2017")?41500:(year=="2018")?58900:(19500+16400+41500+58900);

    for (size_t j=0; j<files.size(); ++j) {
        std::cout << "\n[file " << j+1 << "/" << files.size() << "] " << files[j] << std::endl;

        TFile* f = TFile::Open(files[j].c_str());
        if (!f || f->IsZombie()) { 
            std::cout << "  [skip] cannot open\n"; 
            if (f) f->Close();
            continue; 
        }

        TTree* t = (TTree*)f->Get("Events");
        if (!t) { 
            std::cout << "  [skip] no Events tree\n"; 
            f->Close(); 
            continue; 
        }

        t->SetBranchStatus("*", 0);

        const char* brs[] = {
            "pt_1","pt_2","pt_3","pt_4",
            "eta_1","eta_2","eta_3","eta_4",
            "phi_1","phi_2","phi_3","phi_4",
            "m_1","m_2","m_3","m_4",
            "q_1","q_2","q_3","q_4",
            "d0_1","d0_2","d0_3","d0_4",
            "dZ_1","dZ_2","dZ_3","dZ_4",
            "iso_1","iso_2","iso_3","iso_4",
            "gen_match_1","gen_match_2","gen_match_3","gen_match_4",
            "EleID_WP90_1","EleID_WP90_2","EleID_WP90_3","EleID_WP90_4",
            "MuID_1","MuID_2","MuID_3","MuID_4",
            "TauIDe_1","TauIDe_2","TauIDe_3","TauIDe_4",
            "TauIDm_1","TauIDm_2","TauIDm_3","TauIDm_4",
            "TauIDj_1","TauIDj_2","TauIDj_3","TauIDj_4",
            "IDSF_1","IDSF_2","IDSF_3","IDSF_4",
            "ISOSF_1","ISOSF_2","ISOSF_3","ISOSF_4",
            "TrigSF_1","TrigSF_2","TrigSF_3","TrigSF_4",
            "isTrig_1","isTrig_2",
            "TauVsEleIDSF_1","TauVsEleIDSF_2","TauVsEleIDSF_3","TauVsEleIDSF_4",
            "TauVsMuIDSF_1","TauVsMuIDSF_2","TauVsMuIDSF_3","TauVsMuIDSF_4",
            "TauVsJetIDSF_1","TauVsJetIDSF_2","TauVsJetIDSF_3","TauVsJetIDSF_4",
            "cat","met","run",
            "Generator_weight","brWeight","L1PreFiringWeight_Nom","weightPUtruejson"
        };

        for (auto br : brs) t->SetBranchStatus(br, 1);

        MyBranch(t);  

        TH1D* hnevts = (TH1D*)f->Get("hNWEvts");
        if (!hnevts) hnevts = (TH1D*)f->Get("hNEvts");
        std::string filename = files[j];

        double denom = hnevts ? hnevts->Integral() : 0;
        bool isData = (XSec(files[j].c_str()) == 1);
        bool isDY = (filename.find("DYJets") != std::string::npos);
        bool isWJ = (filename.find("WJets") != std::string::npos);
        bool isTT = (filename.find("TTTo") != std::string::npos);
        bool isQCD = (filename.find("QCD") != std::string::npos);
        bool isDY50 = (filename.find("DYJetsToLLM50") != std::string::npos);
        double xsw = (XSec(files[j].c_str())!=1 && denom>0)? lumi*XSec(files[j].c_str())/denom : 1.0;
        TString outName = Form("%s/hist_%s", outDir.Data(), gSystem->BaseName(files[j].c_str()));
        TFile* fout = new TFile(outName, "RECREATE");
        fout->SetCompressionLevel(0);

        unordered_map<string,TH1D*> h_mZ,h_mH,h_met;
        unordered_map<string,TH1D*> h_ST;
        unordered_map<string,TH1D*> h_pt[NlepMax],h_eta[NlepMax],h_phi[NlepMax],h_d0[NlepMax],h_dZ[NlepMax],h_iso[NlepMax];
        unordered_map<string,TH2D*> h_ST_mZ;
        unordered_map<string,TH2D*> h_pt_eta[NlepMax];
        unordered_map<string,TH1D*> h_genmatch[NlepMax];
        unordered_map<string,TH2D*> h_pt_phi[NlepMax];
        unordered_map<string,TH1D*> h_pt_etaBin1[NlepMax];  // |eta| < 0.5
        unordered_map<string,TH1D*> h_pt_etaBin2[NlepMax];  // 0.5–1.5
        unordered_map<string,TH1D*> h_pt_etaBin3[NlepMax];  // 1.5–2.5
        createHists(h_mZ,"mZ","mZ",1000,0,1000);
        createHists(h_mH,"mH","mll",1000,0,1000);
        createHists(h_met,"met","MET",1000,0,1000);
        createHists(h_ST,"ST","ST",200,0,2000);

        createHists2D(h_ST_mZ,"ST_mZ","ST vs mZ",20,0,2000,1000,0,1000);  

        for(int k=0;k<NlepMax;k++){
            string idx = std::to_string(k+1);
            createHists(h_pt[k],"pt"+idx,"pt"+idx,1000,0,1000);
            createHists(h_eta[k],"eta"+idx,"eta"+idx,30,-3,3);
            createHists(h_phi[k],"phi"+idx,"phi"+idx,30,-3.5,3.5);
            createHists(h_d0[k],"d0_"+idx,"d0_"+idx,30,-0.045,0.045);
            createHists(h_dZ[k],"dZ"+idx,"dZ"+idx,30,-0.1,0.1); 
            createHists(h_iso[k],"iso"+idx,"iso"+idx,30,0,0.30);
            createHists(h_genmatch[k],"genmatch"+idx,"gen_match_"+idx,25, 0, 25);
        }
        for(int k=0; k<NlepMax; k++) {
            string idx = std::to_string(k+1);

            if (doExtra2D) {
                createHists2D(h_pt_eta[k],
                            "pt_eta"+idx,
                            "pt vs eta "+idx,
                            100, 0, 500,
                            30, -3, 3);

                createHists2D(h_pt_phi[k],
                            "pt_phi"+idx,
                            "pt vs phi "+idx,
                            100,0,500,
                            30,-3.5,3.5);
            }

            createHists(h_pt_etaBin1[k],
                        "pt_eta0to0p5_"+idx,
                        "pt |eta|<0.5 "+idx,
                        100,0,500);

            createHists(h_pt_etaBin2[k],
                        "pt_eta0p5to1p5_"+idx,
                        "pt 0.5<|eta|<1.5 "+idx,
                        100,0,500);

            createHists(h_pt_etaBin3[k],
                        "pt_eta1p5to2p5_"+idx,
                        "pt 1.5<|eta|<2.5 "+idx,
                        100,0,500);  
        }
        Long64_t nEnt = t->GetEntriesFast();
        Long64_t pass = 0;

        for(Long64_t i=0;i<nEnt;i++){
            t->GetEntry(i);

            if (i % 1000000 == 0 && i > 0)
                std::cout << "    processed " << i << " / " << nEnt << std::endl;

            string catstr = numberToCat(cat);
            if (std::find(channels.begin(), channels.end(), catstr) == channels.end()) continue;
            if (!passEventCuts(catstr, isData, isDY, isTT, isWJ, isQCD, LOOSE_NOT_TIGHT)) continue;
            vector<Lepton> L;

            if(pt_1>0) L.emplace_back(pt_1,eta_1,phi_1,m_1,q_1,d0_1,dZ_1,iso_1);
            if(pt_2>0) L.emplace_back(pt_2,eta_2,phi_2,m_2,q_2,d0_2,dZ_2,iso_2);
            if(pt_3>0) L.emplace_back(pt_3,eta_3,phi_3,m_3,q_3,d0_3,dZ_3,iso_3);
            if(pt_4>0) L.emplace_back(pt_4,eta_4,phi_4,m_4,q_4,d0_4,dZ_4,iso_4);

            double ST = 0.0;
            for(const auto& lep : L)
                ST += lep.pt;

            bool isFakeTau = (!isData &&
                            (catstr == "eet" || catstr == "mmt") &&
                            gen_match_3 != 5);

            bool isRealTau = (!isData &&
                            (catstr == "eet" || catstr == "mmt") &&
                            gen_match_3 == 5);

            if (!isData && (catstr == "eet" || catstr == "mmt")) {

                if (isRealTau) {
                    if (!passTightTau(3))
                        continue;
                }

                if (isFakeTau) {
                    if (!passLooseNotTightTau(3))
                        continue;
                }
            }

            bool dup=false;
            for(size_t a=0;a<L.size() && !dup;a++){
                for(size_t b=a+1;b<L.size();b++){
                    if(isDuplicate(L[a],L[b])) { dup=true; break; }
                }
            }
            if(dup) continue;

            if (year=="2018" && run>=319077) {
                if(applyHEMveto(catstr)=="yes") continue;
            }

            applyTauES(catstr);

            double evtwt_nom = 1.0;
            double w_FR = 1.0;


            double xsw_corrected = xsw;
            // if (!isData && isDY50) xsw_corrected *= 0.9936;

            if (!isData) {
                            evtwt_nom = Generator_weight * brWeight * xsw_corrected;
                            evtwt_nom *= L1PreFiringWeight_Nom * weightPUtruejson;

                            double id_sf = 1.0, iso_sf = 1.0;
                            double tauEle_sf = 1.0, tauMu_sf = 1.0, tauJet_sf = 1.0;

                            int nLep = catstr.size();
                            
                            if (nLep >= 1) {
                                id_sf *= IDSF_1; iso_sf *= ISOSF_1; 
                                tauEle_sf *= TauVsEleIDSF_1; tauMu_sf *= TauVsMuIDSF_1; tauJet_sf *= TauVsJetIDSF_1;
                            }
                            if (nLep >= 2) {
                                id_sf *= IDSF_2; iso_sf *= ISOSF_2;
                                tauEle_sf *= TauVsEleIDSF_2; tauMu_sf *= TauVsMuIDSF_2; tauJet_sf *= TauVsJetIDSF_2;
                            }
                            if (nLep >= 3) {
                                id_sf *= IDSF_3; iso_sf *= ISOSF_3;
                                tauEle_sf *= TauVsEleIDSF_3; tauMu_sf *= TauVsMuIDSF_3; tauJet_sf *= TauVsJetIDSF_3;
                            }
                            if (nLep >= 4) {
                                id_sf *= IDSF_4; iso_sf *= ISOSF_4;
                                tauEle_sf *= TauVsEleIDSF_4; tauMu_sf *= TauVsMuIDSF_4; tauJet_sf *= TauVsJetIDSF_4;
                            }

                            evtwt_nom *= id_sf * iso_sf * tauEle_sf * tauMu_sf * tauJet_sf;

                            double trig_sf = 1.0;

                            if (nLep < 3) {
                                if (isTrig_1 >= 1)      trig_sf = TrigSF_1;
                                else if (isTrig_1 == -1) trig_sf = TrigSF_2;
                            } 
                            else {
                                if (isTrig_1 >= 1 && isTrig_2 == 0)      trig_sf = TrigSF_1;
                                else if (isTrig_1 == -1 && isTrig_2 == 0) trig_sf = TrigSF_2;
                                else if (isTrig_2 >= 1 && isTrig_1 == 0)  trig_sf = TrigSF_3;
                                else if (isTrig_2 == -1 && isTrig_1 == 0) trig_sf = TrigSF_4;
                                else if (isTrig_1 == 2 && isTrig_2 == 2)  trig_sf = TrigSF_1;
                            }

                            evtwt_nom *= trig_sf;
                        }

            if (isFakeTau) {
                double fData = getTauFakeRate(L[2].pt, L[2].eta, year, catstr, FR_Data);
                if (fData > 0 && fData < 1)
                    w_FR *= (fData / (1 - fData)); 
            }

            if (!isData && (catstr == "eett" || catstr == "mmtt") && passLooseNotTightTau(3) && passLooseNotTightTau(4)) {
                bool tau3_fake  = (gen_match_3 != 5);
                bool tau4_fake  = (gen_match_4 != 5);
                double f1 = 0.0, f2 = 0.0;
                std::string FR_channel = (catstr == "eett") ? "eet" : "mmt";

                if (tau3_fake) {
                    double r = getTauFakeRate(L[2].pt, L[2].eta, year, FR_channel, FR_Data);
                    f1 = r / (1 - r);
                }
                if (tau4_fake) {
                    double r = getTauFakeRate(L[3].pt, L[3].eta, year, FR_channel, FR_Data);
                    f2 = r / (1 - r);
                }

                double w_single = f1 * (1.0 - f2) + f2 * (1.0 - f1);
                double w_double = f1 * f2;

                if (tau3_fake && tau4_fake) w_FR *= w_double;
                else if (tau3_fake || tau4_fake) w_FR *= w_single; 
            }

            if (!isData && (catstr == "ett" || catstr == "mtt") && passLooseNotTightTau(2) && passLooseNotTightTau(3)) {
                bool tau2_fake  = (gen_match_2 != 5);
                bool tau3_fake  = (gen_match_3 != 5);
                double f1 = 0.0, f2 = 0.0;
                std::string FR_channel = (catstr == "ett") ? "eet" : "mmt";

                if (tau2_fake) {
                    double r = getTauFakeRate(L[1].pt, L[1].eta, year, FR_channel, FR_Data);
                    f1 = r / (1 - r);
                }
                if (tau3_fake) {
                    double r = getTauFakeRate(L[2].pt, L[2].eta, year, FR_channel, FR_Data);
                    f2 = r / (1 - r);
                }

                double w_single = f1 * (1.0 - f2) + f2 * (1.0 - f1);
                double w_double = f1 * f2;

                if (tau2_fake && tau3_fake) w_FR *= w_double;
                else if (tau2_fake || tau3_fake) w_FR *= w_single;
            }

            if (!isData && (catstr == "ettt" || catstr == "mttt") && passLooseNotTightTau(2) && passLooseNotTightTau(3) && passLooseNotTightTau(4)) {
                bool tau2_f = (gen_match_2 != 5);
                bool tau3_f = (gen_match_3 != 5);
                bool tau4_f = (gen_match_4 != 5);
                double f1 = 0.0, f2 = 0.0, f3 = 0.0;
                std::string FR_channel = (catstr == "ettt") ? "eet" : "mmt";

                if (tau2_f) { double r = getTauFakeRate(L[1].pt, L[1].eta, year, FR_channel, FR_Data); f1 = r / (1 - r); }
                if (tau3_f) { double r = getTauFakeRate(L[2].pt, L[2].eta, year, FR_channel, FR_Data); f2 = r / (1 - r); }
                if (tau4_f) { double r = getTauFakeRate(L[3].pt, L[3].eta, year, FR_channel, FR_Data); f3 = r / (1 - r); } 

                double w_single = f1*(1-f2)*(1-f3) + f2*(1-f1)*(1-f3) + f3*(1-f1)*(1-f2);
                double w_double = f1*f2*(1-f3) + f1*f3*(1-f2) + f2*f3*(1-f1);
                double w_triple = f1 * f2 * f3;

                int nFake = (int)tau2_f + (int)tau3_f + (int)tau4_f;
                if (nFake == 3)      w_FR *= w_triple;
                else if (nFake == 2) w_FR *= w_double;
                else if (nFake == 1) w_FR *= w_single;
            }

            if (!isData && isDY50 && catstr == "ee") {
                bool flipped = applyChargeFlipStochastic(L, catstr, hProbMC, hSF);
            }
            // w_FR = 1.0;
            double w = evtwt_nom * w_FR;

            if (catstr=="ee"||catstr=="em"||catstr=="et"||
                catstr=="mm"||catstr=="mt"||catstr=="tt") {

                double mll = (LepV(1) + LepV(2)).M();
                string flag = pairFunc_dilep(1, 2, catstr, 15.0);

                ULong64_t gen_matches[4] = {
                    gen_match_1,
                    gen_match_2,
                    gen_match_3,
                    gen_match_4
                };

                std::string Zflag = "";
                if (flag.find("Zwindow") != string::npos)      Zflag = "Zwindow";
                else if (flag.find("Zveto")  != string::npos)  Zflag = "Zveto";

                fillDilep(h_mZ,  catstr, mll, w, L[0].q, L[1].q, Zflag);
                fillDilep(h_mH,  catstr, mll, w, L[0].q, L[1].q, Zflag);
                fillDilep(h_met, catstr, met, w, L[0].q, L[1].q, Zflag);
                fillDilep(h_ST,  catstr, ST,  w, L[0].q, L[1].q, Zflag);
                fillDilep2D(h_ST_mZ, catstr, ST, mll, w, L[0].q, L[1].q, Zflag);

                for (size_t i = 0; i < L.size() && i < 2; ++i) {

                    fillDilep(h_pt[i],  catstr, L[i].pt,  w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_eta[i], catstr, L[i].eta, w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_phi[i], catstr, L[i].phi, w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_d0[i],  catstr, L[i].d0,  w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_dZ[i],  catstr, L[i].dZ,  w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_iso[i], catstr, L[i].iso, w, L[0].q, L[1].q, Zflag);
                    fillDilep(h_genmatch[i], catstr, gen_matches[i], w, L[0].q, L[1].q, Zflag);

                    if (doExtra2D) {
                        fillDilep2D(h_pt_eta[i], catstr,
                                    L[i].pt, L[i].eta,
                                    w, L[0].q, L[1].q, Zflag);

                        fillDilep2D(h_pt_phi[i], catstr,
                                    L[i].pt, L[i].phi,
                                    w, L[0].q, L[1].q, Zflag);
                    }

                    double abseta = std::fabs(L[i].eta);

                    if (abseta < 0.5)
                        fillDilep(h_pt_etaBin1[i], catstr, L[i].pt, w, L[0].q, L[1].q, Zflag);
                    else if (abseta < 1.5)
                        fillDilep(h_pt_etaBin2[i], catstr, L[i].pt, w, L[0].q, L[1].q, Zflag);
                    else if (abseta < 2.5)
                        fillDilep(h_pt_etaBin3[i], catstr, L[i].pt, w, L[0].q, L[1].q, Zflag);
                }
            }
            else {

                string eventFlag = ""; 

                const double mZ = 91.2;
                const double Zwindow = 20.0;

                if (catstr.size() == 3) {

                    double mass12 = (LepV(1) + LepV(2)).M();
                    bool isOS = (q_1 * q_2 < 0);
                    bool inWindow = (fabs(mass12 - mZ) < Zwindow);

                    if (isOS && inWindow)
                        eventFlag = "Zwindow";
                    else if (isOS && !inWindow)
                        eventFlag = "Zveto";

                    fillMultilep(h_mZ, catstr, mass12, w, eventFlag);
                    fillMultilep(h_mH, catstr, mass12, w, eventFlag);
                    fillMultilep2D(h_ST_mZ, catstr, ST, mass12, w, eventFlag);
                }

                else if (catstr.size() == 4) {

                    double mass12 = (LepV(1) + LepV(2)).M();
                    double mass34 = (LepV(3) + LepV(4)).M();

                    bool isOS12 = (q_1 * q_2 < 0);
                    bool isOS34 = (q_3 * q_4 < 0);

                    bool inWin12 = (fabs(mass12 - mZ) < Zwindow);
                    bool inWin34 = (fabs(mass34 - mZ) < Zwindow);

                    bool anyZwindow = (isOS12 && inWin12) || (isOS34 && inWin34);
                    bool anyZveto   = (isOS12 && !inWin12) || (isOS34 && !inWin34);

                    if (anyZwindow)
                        eventFlag = "Zwindow";
                    else if (!anyZwindow && anyZveto)
                        eventFlag = "Zveto";
                    else
                        eventFlag = "";  

                    fillMultilep(h_mZ, catstr, mass12, w, eventFlag);
                    fillMultilep(h_mZ, catstr, mass34, w, eventFlag);

                    fillMultilep(h_mH, catstr, mass12, w, eventFlag);
                    fillMultilep(h_mH, catstr, mass34, w, eventFlag);
                    fillMultilep2D(h_ST_mZ, catstr, ST, mass12, w, eventFlag);
                    fillMultilep2D(h_ST_mZ, catstr, ST, mass34, w, eventFlag);
                }

                fillMultilep(h_met, catstr, met, w, eventFlag);
                fillMultilep(h_ST,  catstr, ST,  w, eventFlag);

                ULong64_t gen_matches[4] = {
                    gen_match_1,
                    gen_match_2,
                    gen_match_3,
                    gen_match_4
                };

                for (size_t idx = 0; idx < L.size() && idx < NlepMax; ++idx) {

                    fillMultilep(h_pt[idx],  catstr, L[idx].pt,  w, eventFlag);
                    fillMultilep(h_eta[idx], catstr, L[idx].eta, w, eventFlag);
                    fillMultilep(h_phi[idx], catstr, L[idx].phi, w, eventFlag);
                    fillMultilep(h_d0[idx],  catstr, L[idx].d0,  w, eventFlag);
                    fillMultilep(h_dZ[idx],  catstr, L[idx].dZ,  w, eventFlag);
                    fillMultilep(h_iso[idx], catstr, L[idx].iso, w, eventFlag);
                    fillMultilep(h_genmatch[idx], catstr, gen_matches[idx], w, eventFlag);

                    if (doExtra2D) {
                        fillMultilep2D(h_pt_eta[idx], catstr,
                                    L[idx].pt, L[idx].eta,
                                    w, eventFlag);

                        fillMultilep2D(h_pt_phi[idx], catstr,
                                    L[idx].pt, L[idx].phi,
                                    w, eventFlag);
                    }

                    double abseta = std::fabs(L[idx].eta);

                    if (abseta < 0.5)
                        fillMultilep(h_pt_etaBin1[idx], catstr, L[idx].pt, w, eventFlag);
                    else if (abseta < 1.5)
                        fillMultilep(h_pt_etaBin2[idx], catstr, L[idx].pt, w, eventFlag);
                    else if (abseta < 2.5)
                        fillMultilep(h_pt_etaBin3[idx], catstr, L[idx].pt, w, eventFlag);
                }
            }

            pass++;
        } 

        scaleAndWrite(h_mZ);
        scaleAndWrite(h_mH);
        scaleAndWrite(h_met);
        scaleAndWrite(h_ST);

        for (int k = 0; k < NlepMax; ++k) {

            scaleAndWrite(h_pt[k]);
            scaleAndWrite(h_eta[k]);
            scaleAndWrite(h_phi[k]);
            scaleAndWrite(h_d0[k]);
            scaleAndWrite(h_dZ[k]);
            scaleAndWrite(h_iso[k]);
            scaleAndWrite(h_genmatch[k]);
            scaleAndWrite(h_pt_etaBin1[k]);
            scaleAndWrite(h_pt_etaBin2[k]);
            scaleAndWrite(h_pt_etaBin3[k]);

            if (doExtra2D) {
                write2D(h_pt_eta[k]);
                write2D(h_pt_phi[k]);
            }
        }
        write2D(h_ST_mZ);
        fout->Close();
        f->Close();

        std::cout << "  [done] processed=" << pass
                << " written -> " << outName << std::endl;
}
}
