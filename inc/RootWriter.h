#ifndef ROOT_WRITER_H
#define ROOT_WRITER_H
#include <map>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TCutG.h"
#include "Particle.h"

const float pion_mass = 0.13957018; //GeV/c^2
const float proton_mass = 0.938272013; //GeV/c^2
const float nucleon_mass = 0.9389186795; //GeV/c^2

enum charge
{
	Neg = 0,
	All = 1,
	Pos = 2
};

struct Histos
{
	TH1I	*histCharged;
	TH1I	*histChargedNeg;
	TH1I	*histChargedPos;

	TH1F	*histMeanCharge;

	TH1F	*histAngle;
	TH1F	*histAngle1;
	TH1F	*histAngle2;
	TH1F	*histAngle3;
	TH1F	*histAngle4;
	TH1F	*histAngleNegNotrot;
	TH1F	*histAngleNeg;
	TH1F	*histAnglePos;
	TH1F	*histTheta;
	TH1F	*histThetaNeg;
	TH1F	*histThetaPos;
	TH1F	*histThetacms;
	TH1F	*histThetacmsNeg;
	TH1F	*histThetacmsPos;

	TH1F	*histYpi;
	TH1F	*histYpiNeg;
	TH1F	*histYpiPos;
	TH1F	*histYcms;
	TH1F	*histYcmsNeg;
	TH1F	*histYcmsPos;
	TH1F	*histYprotcms;
	TH1F	*histYprotcmsNeg;
	TH1F	*histYprotcmsPos;
	TH1F	*histEta;
	TH1F	*histEtaNeg;
	TH1F	*histEtaPos;
	TH1F	*histEtacms;
	TH1F	*histEtacmsNeg;
	TH1F	*histEtacmsPos;

	TH1F	*histPtWide;
	TH1F	*histPtAll;
	TH1F	*histPtNeg;
	TH1F	*histPtPos;
	TH1F	*histPzAll;
	TH1F	*histPzNeg;
	TH1F	*histPzPos;
	TH1F	*histPzcmsAll;
	TH1F	*histPzcmsNeg;
	TH1F	*histPzcmsPos;
	TH1F	*histMeanPt;
	TH1F	*histMeanPtNeg;
	TH1F	*histMeanPtPos;

	TH2F	*histPtVsYAll;
	TH2F	*histPtVsYNeg;
	TH2F	*histPtVsYPos;
	TH2F	*histPtVsYprotAll;
	TH2F	*histPtVsYprotNeg;
	TH2F	*histPtVsYprotPos;

	TH2F	*histPhiVsPtAll;
	TH2F	*histPhiVsPtPos;
	TH2F	*histPhiVsPtNeg;

	TH2F	*histDyDphiAll;
	TH2F	*histDyDphiPos;
	TH2F	*histDyDphiNeg;
	TH2F	*histDyDphiUnlike;
	TH2F	*histDetaDphiAll;
	TH2F	*histDetaDphiAll1;
	TH2F	*histDetaDphiAll2;
	TH2F	*histDetaDphiAll3;
	TH2F	*histDetaDphiAll4;
	TH2F	*histDetaDphiPos;
	TH2F	*histDetaDphiNeg;
	TH2F	*histDetaDphiUnlike;

	TH1D	*histInvMass;

	TH2F	*histDedx_DetaDphiUnlike_05;
	TH2F	*histDedx_DyDphiUnlike_05;
	TH2F	*histDedx_DetaDphiUnlike_025;
	TH2F	*histDedx_DyDphiUnlike_025;

	TH2F	*histDedx;
	TH2F	*histDedxPos;
	TH2F	*histDedxNeg;

	//TH2F	*histDedxSelected;
	//TH2F	*histDedxPosSelected;
	//TH2F	*histDedxNegSelected;

	TH2F	*histDedxVtpc1;
	TH2F	*histDedxVtpc1Pos;
	TH2F	*histDedxVtpc1Neg;

	TH2F	*histDedxVtpc2;
	TH2F	*histDedxVtpc2Pos;
	TH2F	*histDedxVtpc2Neg;

	TH2F	*histDedxMtpc;
	TH2F	*histDedxMtpcPos;
	TH2F	*histDedxMtpcNeg;

	TH1I	*histnDedx;
	TH1I	*histnDedxPos;
	TH1I	*histnDedxNeg;

	TH1I	*histnDedxVtpc1;
	TH1I	*histnDedxVtpc1Pos;
	TH1I	*histnDedxVtpc1Neg;

	TH1I	*histnDedxVtpc2;
	TH1I	*histnDedxVtpc2Pos;
	TH1I	*histnDedxVtpc2Neg;

	TH1I	*histnDedxMtpc;
	TH1I	*histnDedxMtpcPos;
	TH1I	*histnDedxMtpcNeg;

	TH3I	*histPartPopMatrixPos;
	TH3I	*histPartPopMatrixNeg;

	void init();
	void write();
	void clear();
	void LogBinning(TH2F*);
};

class Particles //klasa liczaca zmienne wszystkich trackow w JEDNYM evencie
{
	Histos *histos;
	Float_t	angle,
			theta, theta_cms,
			y_pi, y_pi_cms,
			y_proton_cms,
			eta, eta_cms,
			E_pi, E_proton,
			px, py, pz, p, pt, pz_cms,
			mean_pt[3];

	charge		particle_charge;
	UInt_t n[3];

public:
	static Float_t y_cms;
	static Float_t beta;
	static Float_t gamma;
	static Float_t gamma_beta_e;

	Particles() {}

	void init(Histos *histograms, const float ener);
	void newEvent(bool first = false);
	void analyze(Particle*, const int);
	static Float_t choose_dedx(Particle* particle)
	{
		static Int_t vtpc1_part;
		static Int_t vtpc2_part;
		static Int_t mtpc_part;

		vtpc1_part = 0;
		vtpc2_part = 0;
		mtpc_part = 0;

		vtpc1_part = particle->GetNdEdxVtpc1();
		vtpc2_part = particle->GetNdEdxVtpc2();
		mtpc_part = particle->GetNdEdxMtpc();

		//std::cout << "dE/dx: VTPC1 part: " << vtpc1_part << "\tVTPC2 part: " << vtpc2_part << "\tMTPC part: " << mtpc_part << std::endl;
		if((vtpc1_part == 0) && (vtpc2_part == 0) && (mtpc_part == 0))
		{
			std::cout << "WTF? Particle with no dE/dx information!" << std::endl;
			return 0;
		}
		else
		{
			if(mtpc_part > 0)
				return (particle->GetdEdxMtpc());
			else if(vtpc2_part >= vtpc1_part)
				return (particle->GetdEdxVtpc2());
			else
				return (particle->GetdEdxVtpc1());
		}
	}
	
	static float inline calc_beta(float ener) { return (ener/(ener+nucleon_mass));}
	static float inline calc_gamma(float ener) { return (1/(TMath::Sqrt(1-TMath::Power(calc_beta(ener),2))));}
	float inline calc_gbE(float ener) { return (gamma_beta_e = beta*gamma*ener);}
};

std::map<int,int> getDistro(TString);
float mk_angle3(float);
#endif
