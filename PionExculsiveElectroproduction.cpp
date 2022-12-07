#include"PionExculsiveElectroproduction.h"

#include<iostream>
#include<string.h>
#include<cmath>
#include"TMath.h"

using namespace std;


PionExculsiveElectroproduction::PionExculsiveElectroproduction(){
	cout<<"****EicC Meson Form Factor Project"<<endl;
	cout<<"****Coding issues, contact rwang@impcas.ac.cn"<<endl;
	cout<<endl<<endl;
	cout<<"    Simulation starting..."<<endl;
	cout<<"    Exclusive process: e p --> e n pi^+"<<endl;

	beam_cross_angle = 0.05; //// 50 mrad = 0.05rad
	///// the kinematical ranges for MC sampling
	xBmin = 0.0001;
	xBmax = 0.95;
	Q2min = 1;
	Q2max = 50;
	t0 = -0.001;
	t1 = -100;
	ymin = 0.01;
	ymax = 0.99;
	Tmin = 0.01;
	Tmax = 50;
	//// nucleon mass and electron mass
	mN = 0.938272;
	mpi = 0.139;
	me = 0.000511;
	cutoff_pi = 0.71; //GeV

	PI = TMath::Pi();  //// 3.141592653;
	alpha = 1.0/137.0;

	strFileName = new char[100];
	strcpy(strFileName, "DEMP_pion.root");

	//// EicC optimal collision energy
	eBeamE = 3.5;
	pBeamE = 20;
	eBeam = new TLorentzVector(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam = new TLorentzVector(p_px, 0, p_pz, pBeamE);
	BoostToEIC = new TVector3(0,0,0);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
	
	//// initialize the final-state particles
	elec_out = new TLorentzVector(0, 0, 0, 0);
	neut_out = new TLorentzVector(0, 0, 0, 0);
	pip_out = new TLorentzVector(0, 0, 0, 0);

	/// random seed
	random = new TRandom3(0);
	/// kinematic calculator
	kine = new KineCal();
	/// phase space generator
	eventGenerator = new TGenPhaseSpace();
}
PionExculsiveElectroproduction::~PionExculsiveElectroproduction(){
	tree->Write();
	fout->Write();
	fout->Close();
	cout<<"    Data file saved and closed~"<<endl<<endl;
	//delete random;
	//delete tree;
	//delete fout;
	//delete kine;
	//delete eBeam;
	//delete pBeam;
	//delete elec_out;
	//delete neut_out;
}

int PionExculsiveElectroproduction::Generate(int N = 20000){

	MakeROOTFile(strFileName);

	cout<<"    To generate "<<N<<" events..."<<endl;

	eBeam->Boost(-(*BoostToEIC));   /// boost to the proton-rest frame
	TVector3 zdirection = eBeam->Vect();
	TVector3 ydirection(0, 1, 0);
	TVector3 xdirection(zdirection.Z(), 0, -zdirection.X());

	for(int i=0; i<N; ){
		xB = random->Uniform(xBmin, xBmax);
		Q2 = random->Uniform(Q2min, Q2max);
		if(Q2>(xB*ymax*(s-mN*mN)))continue;
		W2 = (1.0/xB-1)*Q2 + mN*mN;
		//if(W2<2.6)continue;
		if(W2<3.5)continue;
		y = Q2 / xB / (s-mN*mN);
		epsilon = kine->calEpsilon(y, Q2, s);
		// get the physical range of t and compare to [-Tmax, -Tmin]
		t0 = kine->calTMin(W2, -Q2, mN*mN, mpi*mpi, mN*mN);
		if(!(t0==t0))continue;
		if(t0<(-Tmax))continue;
		if(t0>(-Tmin))t0 = -Tmin;
		t1 = kine->calTMax(W2, -Q2, mN*mN, mpi*mpi, mN*mN);
		if(!(t1==t1))continue;
		if(t1>(-Tmin))continue;
		if(t1<(-Tmax))t1 = -Tmax;
		//cout<<"t0 "<<t0<<endl<<t1<<endl;  //test code.
		t = random->Uniform(t1, t0);		

		// calculate the kinematics of the final-state electron in the proton-rest frame
		double nv = Q2/2.0/mN/xB;
		double eOutE = eBeamENRest - nv;
		//cout<<eBeamENRest<<"\t"<<nv<<endl;   // test code.
		double costheta_e = 1 - Q2/2.0/eBeamENRest/eOutE;
		double sintheta_e = sqrt(1 - costheta_e*costheta_e);
		double ephi = random->Uniform(0, 2*PI);
		double emom = sqrt(eOutE*eOutE - me*me);
		TVector3 emom_v3 = emom*sintheta_e*cos(ephi) * xdirection.Unit();
		emom_v3 += emom*sintheta_e*sin(ephi) * ydirection.Unit();
		emom_v3 += emom*costheta_e * zdirection.Unit();
		elec_out->SetXYZT(emom_v3.X(), emom_v3.Y(), emom_v3.Z(), eOutE);
		// calculate the kinematics of the final-state neutron and pion^+ in the proton-rest frame
		double En = mN - t/2.0/mN;
		double Epip = nv + mN - En;
		double Pn = sqrt(En*En - mN*mN);
		double Ppip = sqrt(Epip*Epip - mpi*mpi);
		TLorentzVector virtualPhoton = (*eBeam) - (*elec_out);
		TVector3 elec_out_v3 = elec_out->Vect(); 
		TVector3 virtualPhoton_v3 = virtualPhoton.Vect();
		TVector3 normal_v3 = virtualPhoton_v3.Cross(elec_out_v3);
		TVector3 leptonplane_transverse_v3 = normal_v3.Cross(virtualPhoton_v3);
		//// check the kinematics
		if(virtualPhoton_v3.Mag()>=(Ppip+Pn))continue;
		if(virtualPhoton_v3.Mag()<=fabs(Ppip-Pn))continue;
		double costheta_pip = (virtualPhoton_v3.Mag2()+Ppip*Ppip-Pn*Pn) /2.0 /virtualPhoton_v3.Mag() /Ppip;
		double costheta_n = (virtualPhoton_v3.Mag2()+Pn*Pn-Ppip*Ppip) /2.0 /virtualPhoton_v3.Mag() /Pn;
		double sintheta_pip = sqrt(1 - costheta_pip*costheta_pip);
		double sintheta_n = sqrt(1 - costheta_n*costheta_n);
		double Phi = random->Uniform(0, 2*PI); //the angle between the hadronic plane and the leptonic plane
		TVector3 pip_out_v3 = (Ppip*costheta_pip)*(virtualPhoton_v3.Unit());
		pip_out_v3 = pip_out_v3 + (Ppip*sintheta_pip*cos(Phi))*(leptonplane_transverse_v3.Unit());
		pip_out_v3 = pip_out_v3 + (Ppip*sintheta_pip*sin(Phi))*(normal_v3.Unit());
		TVector3 neut_out_v3 = (Pn*costheta_n)*(virtualPhoton_v3.Unit());
		neut_out_v3 = neut_out_v3 + (Pn*sintheta_n*cos(Phi+PI))*(leptonplane_transverse_v3.Unit());
		neut_out_v3 = neut_out_v3 + (Pn*sintheta_n*sin(Phi+PI))*(normal_v3.Unit());
		pip_out->SetXYZT(pip_out_v3.Px(), pip_out_v3.Py(), pip_out_v3.Pz(), Epip);
		neut_out->SetXYZT(neut_out_v3.Px(), neut_out_v3.Py(), neut_out_v3.Pz(), En);
		//cout<<costheta_pip<<"  "<<costheta_n<<"  "<<pip_out_v3.Pz()<<"  "<<neut_out_v3.Pz()<<endl;  //test code

		//// Boost from proton rest frame to collider frame!
		elec_out->Boost(*BoostToEIC);
		pip_out->Boost(*BoostToEIC);
		neut_out->Boost(*BoostToEIC);

		//// calculate the differential cross sections
		d4sigma = d4sigma_dQ2dxBdtdPhi(Q2, xB, t, 0);
		d3sigma = d3sigma_dQ2dxBdt(Q2, xB, t);

		tree->Fill();
		i++;
	}

	eBeam->Boost(*BoostToEIC);   /// the elec. beam boost back to the collider frame!!!

	cout<<"    Event generation done! "<<endl;
	return N;
}

//// Create a ROOT file and a TTree.
void PionExculsiveElectroproduction::MakeROOTFile(char *filename){
	//// create the output file and the output TTree
	cout<<"    Creating the output file: "<<filename<<endl;
	fout = new TFile(filename,"recreate");
	tree = new TTree("tree","Exclusive pion electroproduction");
	tree->Branch("xB", &xB, "xB/D");
	tree->Branch("Q2", &Q2, "Q2/D");
	tree->Branch("t", &t, "t/D");
	tree->Branch("y", &y, "y/D");
	tree->Branch("W2", &W2, "W2/D");
	tree->Branch("epsilon", &epsilon, "epsilon/D");
	tree->Branch("s", &s, "s/D");
	tree->Branch("d4sigma", &d4sigma, "d4sigma/D");
	tree->Branch("d3sigma", &d3sigma, "d3sigma/D");
	tree->Branch("elec_out", "TLorentzVector", elec_out);
	tree->Branch("neut_out", "TLorentzVector", neut_out);
	tree->Branch("pip_out" , "TLorentzVector", pip_out);
}
void PionExculsiveElectroproduction::SetOutputFileName(char *filename){
	strcpy(strFileName, filename);
}

void PionExculsiveElectroproduction::SetElecBeamEnergy(double ebeamenergy){
	if(ebeamenergy<0.001){cout<<"Error: electron beam energy is too small!!!"<<endl; return;}
	if(ebeamenergy>1e6){cout<<"Error: electron beam energy is too high!!!"<<endl; return;}
	eBeamE = ebeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double PionExculsiveElectroproduction::GetElecBeamEnergy(){return eBeamE;}
void PionExculsiveElectroproduction::SetProtBeamEnergy(double pbeamenergy){
	if(pbeamenergy<1){cout<<"Error: proton beam energy is too small!!!"<<endl; return;}
	if(pbeamenergy>1e6){cout<<"Error: proton beam energy is too high!!!"<<endl; return;}
	pBeamE = pbeamenergy;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double PionExculsiveElectroproduction::GetProtBeamEnergy(){return pBeamE;}
//// set beam crossing angle
void PionExculsiveElectroproduction::SetBeamCrossAngle(double _angle){
	beam_cross_angle = _angle;
	eBeam->SetXYZT(0, 0, sqrt(eBeamE*eBeamE-me*me), eBeamE);
	double p_pz = cos(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	double p_px = sin(PI-beam_cross_angle) * sqrt(pBeamE*pBeamE-mN*mN);
	pBeam->SetXYZT(p_px, 0, p_pz, pBeamE);
	*BoostToEIC = pBeam->BoostVector();
	s = (*eBeam + *pBeam).M2();
	eBeam->Boost(-(*BoostToEIC));
	eBeamENRest = eBeam->E();
	eBeam->Boost(*BoostToEIC);
}
double PionExculsiveElectroproduction::GetBeamCrossAngle(){return beam_cross_angle;}





//// the model of pion form factor
double PionExculsiveElectroproduction::FF_pion(double _Q2){
	return 1.0 / (1.0 + _Q2/cutoff_pi/cutoff_pi);
}
//// return cross-section in the unit of nb/GeV^4.
double PionExculsiveElectroproduction::d4sigma_dQ2dxBdtdPhi(double _Q2, double _xB, double _t, double _Phi)
{
	//// No implementation of the Phi-dependence yet
	return d3sigma_dQ2dxBdt(_Q2, _xB, _t) / 2.0 / PI;
}
//// return cross-section in the unit of nb/GeV^4.
double PionExculsiveElectroproduction::d3sigma_dQ2dxBdt(double _Q2, double _xB, double _t)
{
	double flux = PhotonFlux(y, _xB, epsilon, _Q2); 
	/// in the natural unit, i.e., GeV^-6
	double sigma = flux * (dsigmaT() + epsilon*dsigmaL());
	/// return in the unit of nb/GeV^4
	return 3.8809e5 * sigma;
}
double PionExculsiveElectroproduction::N_factor(double _W2, double _Q2){
	double W2_mN2 = _W2 - mN*mN;
	return 32*PI*W2_mN2 * sqrt(W2_mN2*W2_mN2 + _Q2*_Q2 + 2*_Q2*(W2+mN*mN) );
}
double PionExculsiveElectroproduction::g_piNN(double _t){
	return 13.4 * (cutoff_pi*cutoff_pi - mpi*mpi) / (cutoff_pi*cutoff_pi - _t);
}
double PionExculsiveElectroproduction::PhotonFlux(double _y, double _xB, double _epsilon, double _Q2){
	return alpha*_y*_y*(1-_xB)/2.0/PI/_xB/(1-_epsilon)/_Q2; 
}
double PionExculsiveElectroproduction::dsigmaT(){
	return 0; //// the transverse component can be ignored at very small |t| and high Q2
}
double PionExculsiveElectroproduction::dsigmaL(){
	double nfactor = N_factor(W2, Q2);
	double ffpion = FF_pion(Q2);
	double gpinn = g_piNN(t);
	double pipole = -t / (t-mpi*mpi) / (t-mpi*mpi);
	//// this is the Born-term contribution of pion pole, valid at small |t|.
	return 16*PI*alpha * gpinn*gpinn * pipole * Q2 * ffpion*ffpion / nfactor;
}




//// set sampling ranges
void PionExculsiveElectroproduction::SetxBmin(double min){xBmin = min;}
void PionExculsiveElectroproduction::SetxBmax(double max){xBmax = max;}
void PionExculsiveElectroproduction::SetQ2min(double min){Q2min = min;}
void PionExculsiveElectroproduction::SetQ2max(double max){Q2max = max;}
void PionExculsiveElectroproduction::SetTmin(double min){Tmin = min;}
void PionExculsiveElectroproduction::SetTmax(double max){Tmax = max;}
void PionExculsiveElectroproduction::Setymin(double min){ymin = min;}
void PionExculsiveElectroproduction::Setymax(double max){ymax = max;}


