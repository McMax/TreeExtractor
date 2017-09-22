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
const float neutron_mass = 0.939565346; //GeV/c^2

enum charge
{
	Neg = 0,
	All = 1,
	Pos = 2
};

//All histograms
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
	TH1F	*histPtot;

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
	TH2F	*histDetaDphiPos;
	TH2F	*histDetaDphiNeg;
	TH2F	*histDetaDphiUnlike;

	TH2F	*histDetaDphiAllReflected;
	TH2F	*histDetaDphiPosReflected;
	TH2F	*histDetaDphiNegReflected;
	TH2F	*histDetaDphiUnlikeReflected;

	TH1D	*histInvMass;

	TH2F	*histDedx;
	TH2F	*histDedxPos;
	TH2F	*histDedxNeg;

	TH2F	*histDedxVtpc1;
	TH2F	*histDedxVtpc1Pos;
	TH2F	*histDedxVtpc1Neg;

	TH2F	*histDedxVtpc2;
	TH2F	*histDedxVtpc2Pos;
	TH2F	*histDedxVtpc2Neg;

	TH2F	*histDedxMtpc;
	TH2F	*histDedxMtpcPos;
	TH2F	*histDedxMtpcNeg;

	TH3I	*histPartPopMatrixPos;
	TH3I	*histPartPopMatrixNeg;

	TH1D	*histTTAverageDistance;

	void init(const float momentum);
	void write();
	void clear();
	void LogBinning(TH2F*);
};

//All kinematic variables needed to describe the particle
class Particles
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

	int init(Histos *histograms, const TString system, const float momentum);
	void newEvent(bool first = false);
	void analyze(Particle*, const int);

	static float calc_beta_sym(float b_momentum)	//Calculation of c.m.s. beta with assumption of symmetric system
	{
		return (b_momentum/(TMath::Sqrt(b_momentum*b_momentum+nucleon_mass*nucleon_mass)+nucleon_mass));
	}
	
	//calculation of c.m.s. beta for assymetric system
	static float calc_beta_asym(float b_momentum, int beam_protons, int beam_neutrons, int target_protons, int target_neutrons)
	{ 
		double total_beam_momentum = b_momentum*(beam_protons+beam_neutrons);
		std::cout << "Total beam momentum: " << total_beam_momentum << std::endl;
		double beam_mass = beam_protons*proton_mass + beam_neutrons*neutron_mass;
		double target_mass = target_protons*proton_mass + target_neutrons*neutron_mass;
		double beam_energy = TMath::Sqrt(total_beam_momentum*total_beam_momentum + beam_mass*beam_mass);
		std::cout << "Beam energy: " << beam_energy << std::endl;

		return (total_beam_momentum/(beam_energy+target_mass));
	}
	static float inline calc_gamma(float momentum) { return (1/(TMath::Sqrt(1-TMath::Power(beta,2))));}
	float inline calc_gbE(float ener) { return (gamma_beta_e = beta*gamma*ener);}
};

float mk_angle3(float);
#endif
