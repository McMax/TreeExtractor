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

	TH1F	*histYreal;
	TH1F	*histYrealNeg;
	TH1F	*histYrealPos;
	TH1F	*histYpi;
	TH1F	*histYpiNeg;
	TH1F	*histYpiPos;
	TH1F	*histEta;
	TH1F	*histEtaNeg;
	TH1F	*histEtaPos;

	TH1F	*histPtWide;
	TH1F	*histPtAll;
	TH1F	*histPtNeg;
	TH1F	*histPtPos;
	TH1F	*histPzAll;
	TH1F	*histPzNeg;
	TH1F	*histPzPos;
	TH1F	*histMeanPt;
	TH1F	*histMeanPtNeg;
	TH1F	*histMeanPtPos;
	TH1F	*histPtot;

	TH1F	*histEtot;
	TH1F	*histEtotCALM;

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

	TH3I	*histPartPopMatrixPos;
	TH3I	*histPartPopMatrixNeg;

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
			theta,
			mass, y,
			y_pi,
			eta,
			E_real, E_pi, E_proton,
			px, py, pz, p, pt,
			mean_pt[3];

	charge		particle_charge;
	UInt_t n[3];
	Int_t gpid;

public:

	Particles() {}

	void init(Histos *histograms, const float momentum);
	void newEvent(bool first = false);
	void analyze(Particle*, const int);
};

float mk_angle3(float);
#endif
