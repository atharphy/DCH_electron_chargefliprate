#include "Kinematics.h"
#include "cat.h" //cat sort funs

float* SortPt(string cat, string tau){
	float pt[] = {pt_1,pt_2,pt_3,pt_4};
	int n = cat_lepCount(cat, 'e','m');
	std::vector<float> sort_pt;
	if (tau=="t"){
		for (int j = 0; j < cat.length(); j++){
			if (cat[j] == 't')
				sort_pt.push_back(pt[j]);
		}
	}	
	else{
		for (int j = 0; j < cat.length(); j++){
			if (cat[j] != 't')
				sort_pt.push_back(pt[j]);
		}
	}	
	float *pt_arr = new float[0]; //any number to initialize it.
	std::sort(sort_pt.begin(), sort_pt.end(), greater());//descending order
	std::copy(sort_pt.begin(), sort_pt.begin()+n, pt_arr);//vector->array	
	return pt_arr;
}
	
float ST(string cat){//for only leptons
	float st = 0;
	if (cat.find("e")+1==1 or cat.find("m")+1==1)
		st += pt_1;
	if (cat.find("e")+1==2 or cat.find("m")+1==2)
		st += pt_2; 
	if (cat.find("e")+1==3 or cat.find("m")+1==3)
		st += pt_3;
	if (cat.find("e")+1==4 or cat.find("m")+1==4)
		st += pt_4; 
	return st;
}
float getDR(float eta1, float phi1, float eta2, float phi2) {
    float pi = TMath::Pi();
    float dPhi = fmin(fabs(phi2 - phi1), 2.0 * pi - fabs(phi2 - phi1));
    float DR = sqrt(pow(dPhi, 2) + pow(eta2 - eta1, 2));
    return DR;
}

TLorentzVector LepV(int n){
	TLorentzVector lepV;
	if (n==1)
		lepV.SetPtEtaPhiM(pt_1, eta_1, phi_1, m_1);
	if (n==2)
		lepV.SetPtEtaPhiM(pt_2, eta_2, phi_2, m_2);
	if (n==3)
		lepV.SetPtEtaPhiM(pt_3, eta_3, phi_3, m_3);
	if (n==4)
		lepV.SetPtEtaPhiM(pt_4, eta_4, phi_4, m_4);
	return lepV;
}	

TLorentzVector *ZCandMaker(string cat, TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4, float Zwindow){
	float mZ = 91.2;
	TLorentzVector *Zpaired_lep = new TLorentzVector[4] ;
	TLorentzVector L1, L2, L3, L4;
	
	if(cat.length() == 3){
		if( cat[0]==cat[2] and abs((l1+l3).M()-mZ) < Zwindow){
			if(l1.Pt() >= l3.Pt()){
				L1 = l1; L2 = l3;
				L3 = l2; L4 = l4;
			}
			else{
				L1 = l3; L2 = l1;
				L3 = l2; L4 = l4;
			}
		}
		else if( cat[1]==cat[2] and abs((l2+l3).M()-mZ) < Zwindow){
			if(l2.Pt() >= l3.Pt()){
				L1 = l2; L2 = l3;
				L3 = l1; L4 = l4;
			}
			else{
				L1 = l3; L2 = l2;
				L3 = l1; L4 = l4;	
			}
		}
		if( cat[0]==cat[2] and cat[1]==cat[2] and abs((l1+l3).M()-mZ) < Zwindow and abs((l2+l3).M()-mZ) < Zwindow){
			if( l1.Pt() >= l3.Pt() and l3.Pt() >= l2.Pt()){
				L1 = l1; L2 = l3;
				L3 = l2; L4 = l4;
			}
			else if( l2.Pt() >= l3.Pt() and l3.Pt() >= l1.Pt()){
				L1 = l2; L2 = l3;
				L3 = l1; L4 = l4;	
			}
			else if( l3.Pt() >= l1.Pt() and l1.Pt() >= l1.Pt()){
				L1 = l3; L2 = l1;
				L3 = l2; L4 = l4;	
			}
			else if( l3.Pt() >= l2.Pt() and l2.Pt() >= l1.Pt()){
				L1 = l3; L2 = l2;
				L3 = l1; L4 = l4;	
			}
		}
	}
	else if(cat.length() == 4){//needs fixing!!
		if( cat[0]==cat[2] and abs((l1+l3).M()-mZ) < Zwindow){
			if(l1.Pt() >= l3.Pt()){
				L1 = l1; L2 = l3;
				L3 = l2; L4 = l4;
			}
			else{
				L1 = l3; L2 = l1;
				L3 = l2; L4 = l4;
			}
		}
		else if( cat[0]==cat[3] and abs((l1+l4).M()-mZ) < Zwindow){
			if(l1.Pt() >= l4.Pt()){
				L1 = l1; L2 = l4;
				L3 = l2; L4 = l3;
			}
			else{
				L1 = l4; L2 = l1;
				L3 = l2; L4 = l3;
			}
		}
		else if( cat[1]==cat[2] and abs((l2+l3).M()-mZ) < Zwindow){
			if(l2.Pt() >= l3.Pt()){
				L1 = l2; L2 = l3;
				L3 = l1; L4 = l4;
			}
			else{
				L1 = l3; L2 = l2;
				L3 = l1; L4 = l4;	
			}
		}
		else if( cat[1]==cat[3] and abs((l2+l4).M()-mZ) < Zwindow){
			if(l2.Pt() >= l4.Pt()){
				L1 = l2; L2 = l4;
				L3 = l1; L4 = l3;
			}
			else{
				L1 = l4; L2 = l2;
				L3 = l1; L4 = l3;	
			}
		}
	}
	Zpaired_lep[0] = L1;
	Zpaired_lep[1] = L2;
	Zpaired_lep[2] = L3;
	Zpaired_lep[3] = L4;
	return Zpaired_lep;
}		

TLorentzVector *ZVetoMaker(string cat, TLorentzVector l1, TLorentzVector l2, TLorentzVector l3, TLorentzVector l4, float Zwindow){
	float mZ = 91.2;
	TLorentzVector *Zveto_lep = new TLorentzVector[4] ;
	TLorentzVector L1, L2, L3, L4;
	
	if(cat.length() == 3){
		if( cat[0]==cat[2] and abs((l1+l3).M()-mZ) > Zwindow){
			if(l1.Pt() >= l3.Pt()){
				L1 = l1; L2 = l3;
				L3 = l2; L4 = l4;
			}
			else{
				L1 = l3; L2 = l1;
				L3 = l2; L4 = l4;
			}
		}
		else if( cat[1]==cat[2] and abs((l2+l3).M()-mZ) > Zwindow){
			if(l2.Pt() >= l3.Pt()){
				L1 = l2; L2 = l3;
				L3 = l1; L4 = l4;
			}
			else{
				L1 = l3; L2 = l2;
				L3 = l1; L4 = l4;	
			}
		}
		if( cat[0]==cat[2] and cat[1]==cat[2] and abs((l1+l3).M()-mZ) > Zwindow and abs((l2+l3).M()-mZ) > Zwindow){
			if( l1.Pt() >= l3.Pt() and l3.Pt() >= l2.Pt()){
				L1 = l1; L2 = l3;
				L3 = l2; L4 = l4;
			}
			else if( l2.Pt() >= l3.Pt() and l3.Pt() >= l1.Pt()){
				L1 = l2; L2 = l3;
				L3 = l1; L4 = l4;	
			}
			else if( l3.Pt() >= l1.Pt() and l1.Pt() >= l1.Pt()){
				L1 = l3; L2 = l1;
				L3 = l2; L4 = l4;	
			}
			else if( l3.Pt() >= l2.Pt() and l2.Pt() >= l1.Pt()){
				L1 = l3; L2 = l2;
				L3 = l1; L4 = l4;	
			}
		}
	}
	//else if(cat.length() == 4){
	
	Zveto_lep[0] = L1;
	Zveto_lep[1] = L2;
	Zveto_lep[2] = L3;
	Zveto_lep[3] = L4;
	return Zveto_lep;
}