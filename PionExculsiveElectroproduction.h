#ifndef _PionExculsiveElectroproduction_H
#define _PionExculsiveElectroproduction_H 1

#include"TRandom3.h"
#include"TFile.h"
#include"TTree.h"
#include"TLorentzVector.h"
#include"TVector3.h"
#include"TGenPhaseSpace.h"

#include"KineCal.h"


class PionExculsiveElectroproduction{

	public:
		PionExculsiveElectroproduction();
		~PionExculsiveElectroproduction();

		/// the function to generate and dump N events into root file
		int Generate(int N);

		void SetOutputFileName(char *filename);

		//// set beam energies and crossing angle
		void SetElecBeamEnergy(double ebeamenergy);
		void SetProtBeamEnergy(double pbeamenergy);
		void SetBeamCrossAngle(double _angle);
		double GetBeamCrossAngle();
		double GetElecBeamEnergy();
		double GetProtBeamEnergy();

		//// pion form factor parametrization and the cross-section model
		double FF_pion(double Q2);
		double d4sigma_dQ2dxBdtdPhi(double Q2, double xB, double t, double Phi);
		double d3sigma_dQ2dxBdt(double Q2, double xB, double t);


		//// set sampling ranges
		void SetxBmin(double min);
		void SetxBmax(double max);
		void SetQ2min(double min);
		void SetQ2max(double max);
		void SetTmin(double min);
		void SetTmax(double max);
		void Setymin(double min);
		void Setymax(double max);



	private:
		double me;
		double mpi;
		double mN;
		double cutoff_pi;

		double PI;
		double alpha;

		double xB;
		double Q2;
		double t;
		double y;
		double W2;
		double epsilon;
		double s;

		double d4sigma;
		double d3sigma;

		double beam_cross_angle;
		double xBmin;
		double xBmax;
		double Q2min;
		double Q2max;
		double t0;
		double t1;
		double ymin;
		double ymax;
		double Tmin;
		double Tmax;

		// for final-state electron
		TLorentzVector *elec_out;
		// for final-state neutron
		TLorentzVector *neut_out;
		// for final-state pion plus
		TLorentzVector *pip_out;

		double eBeamE;
		double pBeamE;
		double eBeamENRest;
		TLorentzVector *eBeam;
		TLorentzVector *pBeam;

		TVector3 *BoostToEIC;

		TRandom3 *random;
		TFile *fout;
		TTree *tree;
		KineCal *kine;

		TGenPhaseSpace *eventGenerator;

		char *strFileName;

		void MakeROOTFile(char *filename);

		double N_factor(double _W2, double _Q2);
		double g_piNN(double _t);
		double PhotonFlux(double _y, double _xB, double _epsilon, double _Q2);
		double dsigmaT();
		double dsigmaL();
};


#endif
