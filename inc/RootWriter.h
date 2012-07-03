#ifndef ROOT_WRITER_H
#define ROOT_WRITER_H
#include <map>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
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
	TH1F	*histEta;
	TH1F	*histEtaNeg;
	TH1F	*histEtaPos;
	TH1F	*histEtacms;
	TH1F	*histEtacmsNeg;
	TH1F	*histEtacmsPos;

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
	TH2F	*histPtVsYAll_cut;

	TH2F	*histPhiVsPtAll;
	TH2F	*histPhiVsPtPos;
	TH2F	*histPhiVsPtNeg;

	TH2F	*histDyDphiAll;
	TH2F	*histDyDphiPos;
	TH2F	*histDyDphiNeg;
	TH2F	*histDyDphiUnlike;
	TH2F	*histDetaDphiAll;
	TH2F	*histDetaDphiPos;
	TH2F	*histDetaDphiNeg;
	TH2F	*histDetaDphiUnlike;

	TCutG	*cutg;

	void init();
	void write();
	void clear();
	void init_dEdxCut();
};

class Particles //klasa liczaca zmienne wszystkich trackow w JEDNYM evencie
{
	Histos *histos;
	Float_t	angle,
			theta, theta_cms,
			y_pi, y_pi_cms,
			eta, eta_cms,
			px, py, pz, p, pt, E, pz_cms,
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
	static float inline calc_beta(float ener) { return (ener/(ener+nucleon_mass));}
	static float inline calc_gamma(float ener) { return (1/(TMath::Sqrt(1-TMath::Power(calc_beta(ener),2))));}
	float inline calc_gbE(float ener) { return (gamma_beta_e = beta*gamma*ener);}
};

std::map<int,int> getDistro(TString);
float mk_angle3(float);
#endif
