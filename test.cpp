#include"KineCal.h"
#include"PionExclusiveElectroproduction.h"

#include"KineCal.cpp"
#include"PionExclusiveElectroproduction.cpp"


int test(){

	PionExclusiveElectroproduction demp_pion;

	demp_pion.SetTmax(0.7);
	demp_pion.SetTmin(0.01);
	demp_pion.SetQ2max(50);
	demp_pion.SetQ2min(1);
	demp_pion.SetxBmax(0.8);
	demp_pion.SetxBmin(0.001);

	//char filename[50] = "DEMP-pion-pole-at-EicC.root";
	//demp_pion.SetOutputFileName(filename);
	TString filename = "DEMP-pion-pole-at-EicC.root";
	demp_pion.SetOutputFileName(filename); 
	demp_pion.SetEvtFileOutput(1); 


	demp_pion.SetElecBeamEnergy(3.5);
	demp_pion.SetProtBeamEnergy(20);
	//demp_pion.SetBeamCrossAngle(0.0);
	demp_pion.SetBeamCrossAngle(0.05);

	//demp_pion.SetQuiet(1);
	demp_pion.SetSamplingMode(1);
	demp_pion.Generate(20000);

	//cout<<demp_pion.GetElecBeamEnergy()<<"  "<<demp_pion.GetProtBeamEnergy()<<"  "<<demp_pion.GetBeamCrossAngle()<<endl;

	return 0;
}


