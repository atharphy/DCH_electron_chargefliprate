#include "TMath.h"
#include <cmath>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrix.h"

auto MET_split(TLorentzVector Lep1, TLorentzVector Lep2, TLorentzVector Lep3, TLorentzVector Lep4, double met, double metphi){
	TLorentzVector MET, Vis1, Vis2;
	MET.SetPtEtaPhiM(met, 0, metphi, 0);
	Vis1 = Lep1+Lep2;
	Vis2 = Lep3+Lep4;
	TLorentzVector vis_total = Vis1+Vis2; 
	//cout<<MET.Pt()<<endl;
	TLorentzVector neutr_leg1, neutr_leg2;

	if (mll2 < 0){//if only one vis pair is found
		neutr_leg1 = MET;
	}
	else{
		double ptleg1 = MET.Pt()*0.5;//assumption 1
		//constraint 1: ptleg1+ptleg2 = MET by vector addition
		double A=1, B = 2*ptleg1*TMath::Cos(Vis1.Phi()-Vis2.Phi()), C = ptleg1*ptleg1 - MET.Pt()*MET.Pt(); 
		double ptleg2 = (-B + TMath::Sqrt( B*B-4*A*C ))/(2*A);
		if (ptleg2 < 0){
			ptleg2 = (-B - TMath::Sqrt( B*B-4*A*C ))/(2*A);
		}
		
		float MDCH = (mll+mll2)/2;
		double mleg2 = abs(MDCH-mll2);//assumption 2 
		//constraint 2: mDCH1 = mDCH2 
		double Eleg2 = TMath::Sqrt(mleg2*mleg2 + ptleg2*ptleg2*TMath::CosH(Vis2.Eta())*TMath::CosH(Vis2.Eta()));//E = sqrt(m^2+Pt^2*Cosh^2(eta))
		double J = (Eleg2+Vis2.E()), K = (ptleg2+Vis2.Pt())*TMath::CosH(Vis2.Eta()), L = (ptleg1+Vis1.Pt())*TMath::CosH(Vis1.Eta()); 
		double Eleg1 = -Vis1.E()+TMath::Sqrt(J*J - K*K + L*L);
		if (Eleg1 < 0){
			Eleg1 = -Vis1.E()-TMath::Sqrt(J*J - K*K + L*L);
		}
		
		neutr_leg1.SetPtEtaPhiE(ptleg1, Vis1.Eta(), Vis1.Phi(), Eleg1);
		neutr_leg2.SetPtEtaPhiE(ptleg2, Vis2.Eta(), Vis2.Phi(), Eleg2);	
			
		if(neutr_leg1.M() < 0 || neutr_leg1.E() < 0){//for only one leg 
			ptleg1 = 0; ptleg2 = MET.Pt(); 
			Eleg1 = 0;
			mleg2 = abs(mll-mll2);
			Eleg2 = TMath::Sqrt(mleg2*mleg2 + ptleg2*ptleg2*TMath::CosH(Vis2.Eta())*TMath::CosH(Vis2.Eta()));
			neutr_leg1.SetPtEtaPhiE(ptleg1, Vis1.Eta(), Vis1.Phi(), Eleg1);
			neutr_leg2.SetPtEtaPhiE(ptleg2, Vis2.Eta(), Vis2.Phi(), Eleg2);	
			//printf("Ahhh -ve energies! \n");
		}
	}
	
	struct nuleg_pair {TLorentzVector first; TLorentzVector second;};
	return nuleg_pair {neutr_leg1, neutr_leg2};
}
