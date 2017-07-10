#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "RootWriter.h"
#include "Particle.h"
#include "TFile.h"

void Histos::init(const float momentum)
{
	int detadphibins[2];
//Due to low statistics, pp20 and pp31 have wider bins. This is applied only in ALICE style distributions.
	if(momentum <= 31)
	{
		detadphibins[0] = 26;
		detadphibins[1] = 13;
	}
	else
	{
		detadphibins[0] = 50;
		detadphibins[1] = 25;
	}

	histCharged = new TH1I("histCharged","Multiplicity of charged;N",50,0,50);
	histChargedNeg = new TH1I("histChargedNeg","Multiplicity of negatively charged;N",40,0,40);
	histChargedPos = new TH1I("histChargedPos","Multiplicity of positively charged;N",40,0,40);
	histMeanCharge = new TH1F("histMeanCharge","Mean charge;charge",23,-1.15,1.15);
	histAngle = new TH1F("histAngle","Azimuthal angle;#phi [rad]",50,-TMath::Pi(), TMath::Pi());
	histAngleNegNotrot = new TH1F("histAngleNegNotrot","Azimuthal angle (notrot), neg.;#phi [rad]",50,-TMath::Pi(), TMath::Pi());
	histAngleNeg = new TH1F("histAngleNeg","Azimuthal angle, neg.;#phi [rad]",50,-TMath::Pi(),TMath::Pi());
	histAnglePos = new TH1F("histAnglePos","Azimuthal angle, pos.;#phi [rad]",50,-TMath::Pi(),TMath::Pi());
	histTheta = new TH1F("histTheta","Polar angle;#theta [rad]", 50, -1, 1);
	histThetaNeg = new TH1F("histThetaNeg","Polar angle, neg.;#theta [rad]", 50, -1, 1);
	histThetaPos = new TH1F("histThetaPos","Polar angle, pos.;#theta [rad]", 50, -1, 1);
	histThetacms = new TH1F("histThetacms","Polar angle, CMS;#theta [rad]", 100, -TMath::Pi(), 2*TMath::Pi());
	histThetacmsNeg = new TH1F("histThetacmsNeg","Polar angle, CMS, neg.;#theta [rad]", 50, -TMath::Pi(), 2*TMath::Pi());
	histThetacmsPos = new TH1F("histThetacmsPos","Polar angle, CMS, pos.;#theta [rad]", 50, -TMath::Pi(), 2*TMath::Pi());
	histYpi = new TH1F("histYpi","Rapidity with #pi mass;y_{#pi}",100,-2,8);
	histYpiNeg = new TH1F("histYpiNeg","Rapidity with #pi mass, neg.;y_{#pi}",100,-2,8);
	histYpiPos = new TH1F("histYpiPos","Rapidity with #pi mass, pos.;y_{#pi}",100,-2,8);
	histYcms = new TH1F("histYcms","CMS rapidity with #pi mass;y_{#pi}",100,-5,5);
	histYcmsNeg = new TH1F("histYcmsNeg","CMS rapidity with #pi mass, neg.;y_{#pi}",100,-5,5);
	histYcmsPos = new TH1F("histYcmsPos","CMS rapidity with #pi mass, pos.;y_{#pi}",100,-5,5);
	histYprotcms = new TH1F("histYprotcms","CMS rapidity with proton mass;y_{p}",100,-5,5);
	histYprotcmsNeg = new TH1F("histYprotcmsNeg","CMS rapidity with proton mass, neg.;y_{p}",100,-5,5);
	histYprotcmsPos = new TH1F("histYprotcmsPos","CMS rapidity with proton mass, pos.;y_{p}",100,-5,5);
	histEta = new TH1F("histEta","Pseudorapidity;#eta",100,-2,8);
	histEtaNeg = new TH1F("histEtaNeg","Pseudorapidity, negatively charged;#eta",100,-2,8);
	histEtaPos = new TH1F("histEtaPos","Pseudorapidity, positively charged;#eta",100,-2,8);
	histEtacms = new TH1F("histEtacms","Pseudorapidity;#eta",100,-4,6);
	histEtacmsNeg = new TH1F("histEtacmsNeg","Pseudorapidity, negatively charged;#eta",100,-4,6);
	histEtacmsPos = new TH1F("histEtacmsPos","Pseudorapidity, positively charged;#eta",100,-4,6);
	histPtWide = new TH1F("histPtWide","Transverse momentum (wide);p_{T} [GeV/c]",750,0,150);
	histPtAll = new TH1F("histPtAll","Transverse momentum;p_{T} [GeV/c]",150,0,3);
	histPtNeg = new TH1F("histPtNeg","Transverse momentum, neg.;p_{T} [GeV/c]",150,0,3);
	histPtPos = new TH1F("histPtPos","Transverse momentum, pos.;p_{T} [GeV/c]",150,0,3);
	histPzAll = new TH1F("histPzAll","Longitudinal momentum;p_{z} [GeV/c]",150,0,50);
	histPzNeg = new TH1F("histPzNeg","Longitudinal momentum, neg.;p_{z} [GeV/c]",150,0,50);
	histPzPos = new TH1F("histPzPos","Longitudinal momentum, pos.;p_{z} [GeV/c]",150,0,50);
	histPzcmsAll = new TH1F("histPzcmsAll","Longitudinal momentum, CMS;p_{z} [GeV/c]",200,-5,10);
	histPzcmsNeg = new TH1F("histPzcmsNeg","Longitudinal momentum, CMS, neg.;p_{z} [GeV/c]",200,-5,10);
	histPzcmsPos = new TH1F("histPzcmsPos","Longitudinal momentum, CMS, pos.;p_{z} [GeV/c]",200,-5,10);
	histPtot = new TH1F("histPtot","Total momentum;p_{tot} [GeV/c]",300, 0, 150);
	histMeanPt = new TH1F("histMeanPt","Mean transverse momentum (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtNeg = new TH1F("histMeanPtNeg","Mean transverse momentum, neg. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtPos = new TH1F("histMeanPtPos","Mean transverse momentum, pos. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histPtVsYAll = new TH2F("histPtVsYAll","Trans. momentum vs. rapidity; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYNeg = new TH2F("histPtVsYNeg","Trans. momentum vs. rapidity, neg.; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYPos = new TH2F("histPtVsYPos","Trans. momentum vs. rapidity, pos.; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYprotAll = new TH2F("histPtVsYprotAll","Trans. momentum vs. rapidity; y^{*}_{prot}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYprotNeg = new TH2F("histPtVsYprotNeg","Trans. momentum vs. rapidity, neg.; y^{*}_{prot}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYprotPos = new TH2F("histPtVsYprotPos","Trans. momentum vs. rapidity, pos.; y^{*}_{prot}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPhiVsPtAll = new TH2F("histPhiVsPtAll","Az. angle vs. transverse momentum; #phi [rad]; p_{T} [GeV/c]",50,-TMath::Pi(),TMath::Pi(),150,0,1.5);
	histPhiVsPtPos = new TH2F("histPhiVsPtPos","Az. angle vs. transverse momentum, pos.; #phi [rad]; p_{T} [GeV/c]",50,-TMath::Pi(),TMath::Pi(),150,0,1.5);
	histPhiVsPtNeg = new TH2F("histPhiVsPtNeg","Az. angle vs. transverse momentum, neg.; #phi [rad]; p_{T} [GeV/c]",50,-TMath::Pi(),TMath::Pi(),150,0,1.5);
	histDyDphiAll = new TH2F("histDyDphiAll","#Deltay versus #Delta#phi;#Delta#phi [rad]; #Deltay_{#pi}",200,0,TMath::Pi(),200,0,6.5);
	histDyDphiPos = new TH2F("histDyDphiPos","#Deltay versus #Delta#phi, pos.;#Delta#phi [rad]; #Deltay_{#pi}",200,0,TMath::Pi(),200,0,6.5);
	histDyDphiNeg = new TH2F("histDyDphiNeg","#Deltay versus #Delta#phi, neg.;#Delta#phi [rad]; #Deltay_{#pi}",200,0,TMath::Pi(),200,0,6.5);
	histDyDphiUnlike = new TH2F("histDyDphiUnlike","#Deltay versus #Delta#phi, unlike-sign;#Delta#phi [rad]; #Deltay_{#pi}",200,0,TMath::Pi(),200,0,6.5);

	histDetaDphiAll = new TH2F("histDetaDphiAll","#Delta#eta#Delta#phi;#Delta#phi [rad];#Delta#eta",192,0,TMath::Pi(),192,0,6);
	histDetaDphiPos = new TH2F("histDetaDphiPos","#Delta#eta#Delta#phi, pos.;#Delta#phi [rad];#Delta#eta",192,0,TMath::Pi(),192,0,6);
	histDetaDphiNeg = new TH2F("histDetaDphiNeg","#Delta#eta#Delta#phi, neg.;#Delta#phi [rad];#Delta#eta",192,0,TMath::Pi(),192,0,6);
	histDetaDphiUnlike = new TH2F("histDetaDphiUnlike","#Delta#eta#Delta#phi, unlike-sign;#Delta#phi [rad];#Delta#eta",192,0,TMath::Pi(),192,0,6);

	histDetaDphiAllReflected = new TH2F("histDetaDphiAllReflected","#Delta#eta#Delta#phi (Reflected);#Delta#phi [rad];#Delta#eta",detadphibins[0],-(TMath::Pi()/2),1.5*TMath::Pi(),detadphibins[1],-3,3);
	histDetaDphiPosReflected = new TH2F("histDetaDphiPosReflected","#Delta#eta#Delta#phi (Reflected), pos.;#Delta#phi [rad];#Delta#eta",detadphibins[0],-(TMath::Pi()/2),1.5*TMath::Pi(),detadphibins[1],-3,3);
	histDetaDphiNegReflected = new TH2F("histDetaDphiNegReflected","#Delta#eta#Delta#phi (Reflected), neg.;#Delta#phi [rad];#Delta#eta",detadphibins[0],-(TMath::Pi()/2),1.5*TMath::Pi(),detadphibins[1],-3,3);
	histDetaDphiUnlikeReflected = new TH2F("histDetaDphiUnlikeReflected","#Delta#eta#Delta#phi (Reflected), unlike-sign;#Delta#phi [rad];#Delta#eta",detadphibins[0],-(TMath::Pi()/2),1.5*TMath::Pi(),detadphibins[1],-3,3);

	histInvMass = new TH1D("histInvMass","Invariant mass (assumed #pi mass);m_{inv} [GeV/c^{2}]",5000,0,5);

	histDedx = new TH2F("histDedx","dE/dx (all charged)",400,-3,3,400,0,3);
	histDedxPos = new TH2F("histDedxPos","dE/dx (pos. charged)",400,-3,3,400,0,3);
	histDedxNeg = new TH2F("histDedxNeg","dE/dx (neg. charged)",400,-3,3,400,0,3);
	histDedxVtpc1 = new TH2F("histDedxVtpc1","dE/dx (VTPC1, all charged)",400,-3,3,400,0,3);
	histDedxVtpc1Pos = new TH2F("histDedxVtpc1Pos","dE/dx (VTPC1, pos. charged)",400,-3,3,400,0,3);
	histDedxVtpc1Neg = new TH2F("histDedxVtpc1Neg","dE/dx (VTPC1, neg. charged)",400,-3,3,400,0,3);
	histDedxVtpc2 = new TH2F("histDedxVtpc2","dE/dx (VTPC2, all charged)",400,-3,3,400,0,3);
	histDedxVtpc2Pos = new TH2F("histDedxVtpc2Pos","dE/dx (VTPC2, pos. charged)",400,-3,3,400,0,3);
	histDedxVtpc2Neg = new TH2F("histDedxVtpc2Neg","dE/dx (VTPC2, neg. charged)",400,-3,3,400,0,3);
	histDedxMtpc = new TH2F("histDedxMtpc","dE/dx (MTPC, all charged)",400,-3,3,400,0,3);
	histDedxMtpcPos = new TH2F("histDedxMtpcPos","dE/dx (MTPC, pos. charged)",400,-3,3,400,0,3);
	histDedxMtpcNeg = new TH2F("histDedxMtpcNeg","dE/dx (MTPC, neg. charged)",400,-3,3,400,0,3);

	histTTAverageDistance = new TH1D("histTTAverageDistance","Two track average distance;distance [cm]",500,0,25);

	histPartPopMatrixPos = new TH3I("histPartPopMatrixPos","Particle population matrix, pos. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,150,40,0,2,36,-TMath::Pi(),TMath::Pi());
	histPartPopMatrixNeg = new TH3I("histPartPopMatrixNeg","Particle population matrix, neg. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,150,40,0,2,36,-TMath::Pi(),TMath::Pi());

	LogBinning(histDedx);
	LogBinning(histDedxPos);
	LogBinning(histDedxNeg);
	LogBinning(histDedxVtpc1);
	LogBinning(histDedxVtpc1Pos);
	LogBinning(histDedxVtpc1Neg);
	LogBinning(histDedxVtpc2);
	LogBinning(histDedxVtpc2Pos);
	LogBinning(histDedxVtpc2Neg);
	LogBinning(histDedxMtpc);
	LogBinning(histDedxMtpcPos);
	LogBinning(histDedxMtpcNeg);
}

//Taken from Maja
void Histos::LogBinning(TH2F *hist)
{
	TAxis *axis = hist->GetXaxis();
	int bins = axis->GetNbins();

	Axis_t from = axis->GetXmin();
	Axis_t to = axis->GetXmax();
	Axis_t width = (to - from) / bins;
	Axis_t *new_bins = new Axis_t[bins + 1];

	for (int i = 0; i <= bins; i++) {
		new_bins[i] = TMath::Power(10, from + i * width);

	}
	axis->Set(bins, new_bins);
	delete new_bins;
}

void Histos::write()
{
	histCharged->Write();
	histChargedNeg->Write();
	histChargedPos->Write();
	histMeanCharge->Write();
	histAngle->Write();
	histAnglePos->Write();
	histAngleNeg->Write();
	histAngleNegNotrot->Write();
	histTheta->Write();
	histThetaPos->Write();
	histThetaNeg->Write();
	histThetacms->Write();
	histThetacmsPos->Write();
	histThetacmsNeg->Write();
	histYpi->Write();
	histYpiPos->Write();
	histYpiNeg->Write();
	histYcms->Write();
	histYcmsPos->Write();
	histYcmsNeg->Write();
	histYprotcms->Write();
	histYprotcmsPos->Write();
	histYprotcmsNeg->Write();
	histEta->Write();
	histEtaNeg->Write();
	histEtaPos->Write();
	histEtacms->Write();
	histEtacmsNeg->Write();
	histEtacmsPos->Write();
	histPtWide->Write();
	histPtAll->Write();
	histPtPos->Write();
	histPtNeg->Write();
	histPzAll->Write();
	histPzPos->Write();
	histPzNeg->Write();
	histPzcmsAll->Write();
	histPzcmsPos->Write();
	histPzcmsNeg->Write();
	histMeanPt->Write();
	histMeanPtNeg->Write();
	histMeanPtPos->Write();
	histPtot->Write();
	histPtVsYAll->Write();
	histPtVsYPos->Write();
	histPtVsYNeg->Write();
	histPtVsYprotAll->Write();
	histPtVsYprotPos->Write();
	histPtVsYprotNeg->Write();
	histPhiVsPtAll->Write();
	histPhiVsPtPos->Write();
	histPhiVsPtNeg->Write();
	histDyDphiAll->Write();
	histDyDphiPos->Write();
	histDyDphiNeg->Write();
	histDyDphiUnlike->Write();
	histDetaDphiAll->Write();
	histDetaDphiPos->Write();
	histDetaDphiNeg->Write();
	histDetaDphiUnlike->Write();
	histDetaDphiAllReflected->Write();
	histDetaDphiPosReflected->Write();
	histDetaDphiNegReflected->Write();
	histDetaDphiUnlikeReflected->Write();
	histInvMass->Write();
	histDedx->Write();
	histDedxPos->Write();
	histDedxNeg->Write();
	histDedxVtpc1->Write();
	histDedxVtpc1Pos->Write();
	histDedxVtpc1Neg->Write();
	histDedxVtpc2->Write();
	histDedxVtpc2Pos->Write();
	histDedxVtpc2Neg->Write();
	histDedxMtpc->Write();
	histDedxMtpcPos->Write();
	histDedxMtpcNeg->Write();
	histPartPopMatrixPos->Write();
	histPartPopMatrixNeg->Write();
	histTTAverageDistance->Write();
}

void Histos::clear()
{
	delete	histCharged;
	delete	histChargedNeg;
	delete	histChargedPos;
	delete	histMeanCharge;
	delete 	histAngle;
	delete 	histAnglePos;
	delete 	histAngleNeg;
	delete 	histAngleNegNotrot;
	delete 	histTheta;
	delete 	histThetaPos;
	delete 	histThetaNeg;
	delete 	histThetacms;
	delete 	histThetacmsPos;
	delete 	histThetacmsNeg;
	delete 	histYpi;
	delete 	histYpiPos;
	delete 	histYpiNeg;
	delete 	histYcms;
	delete 	histYcmsPos;
	delete 	histYcmsNeg;
	delete 	histYprotcms;
	delete 	histYprotcmsPos;
	delete 	histYprotcmsNeg;
	delete	histEta;
	delete	histEtaNeg;
	delete	histEtaPos;
	delete	histEtacms;
	delete	histEtacmsNeg;
	delete	histEtacmsPos;
	delete 	histPtWide;
	delete 	histPtAll;
	delete 	histPtPos;
	delete 	histPtNeg;
	delete 	histPzAll;
	delete 	histPzPos;
	delete 	histPzNeg;
	delete 	histPzcmsAll;
	delete 	histPzcmsPos;
	delete 	histPzcmsNeg;
	delete	histMeanPt;
	delete	histMeanPtNeg;
	delete	histMeanPtPos;
	delete	histPtot;
	delete	histPtVsYAll;
	delete	histPtVsYPos;
	delete	histPtVsYNeg;
	delete	histPtVsYprotAll;
	delete	histPtVsYprotPos;
	delete	histPtVsYprotNeg;
	delete histPhiVsPtAll;
	delete histPhiVsPtPos;
	delete histPhiVsPtNeg;
	delete histDyDphiAll;
	delete histDyDphiPos;
	delete histDyDphiNeg;
	delete histDyDphiUnlike;
	delete histDetaDphiAll;
	delete histDetaDphiPos;
	delete histDetaDphiNeg;
	delete histDetaDphiUnlike;
	delete histDetaDphiAllReflected;
	delete histDetaDphiPosReflected;
	delete histDetaDphiNegReflected;
	delete histDetaDphiUnlikeReflected;
	delete histInvMass;
	delete histDedx;
	delete histDedxPos;
	delete histDedxNeg;
	delete histDedxVtpc1;
	delete histDedxVtpc1Pos;
	delete histDedxVtpc1Neg;
	delete histDedxVtpc2;
	delete histDedxVtpc2Pos;
	delete histDedxVtpc2Neg;
	delete histDedxMtpc;
	delete histDedxMtpcPos;
	delete histDedxMtpcNeg;
	delete histPartPopMatrixPos;	
	delete histPartPopMatrixNeg;	
	delete histTTAverageDistance;
}

float mk_angle3(float  angle1)
{
	if(angle1 < 0.)		//Konwersja angle3
		angle1 += TMath::Pi();
	else
		angle1 -= TMath::Pi();

	return angle1;
}

Float_t Particles::y_cms = 0;
Float_t Particles::beta = 0;
Float_t Particles::gamma = 0;
Float_t Particles::gamma_beta_e = 0;

void Particles::init(Histos *histograms, const TString system, const float momentum)
{
	histos = histograms;
	angle = 0.;
	theta = theta_cms = 0.;
	y_pi = y_cms = y_pi_cms = y_proton_cms = 0.;
	eta = eta_cms = 0.;
	pz_cms = pt = 0.;
	particle_charge = All;
	for(int ch=0; ch<3; ch++)
	{
		n[ch] = 0;
		mean_pt[ch] = 0.;
	}

	if(system.EqualTo("pp"))
		//beta = calc_beta_asym(momentum, 1, 0, 1, 0);
		beta = calc_beta_sym(momentum);
	else if(system.EqualTo("BeBe"))
		//beta = calc_beta_asym(momentum, 4, 3, 4, 5);
		beta = calc_beta_sym(momentum);
	else
	{
		std::cout << "System not recognized. Exiting." << std::endl;
		return;
	}
	gamma = calc_gamma(momentum);
	y_cms = 0.5*TMath::Log((1+beta)/(1-beta));
	std::cout << "Beta factor: " << beta << std::endl;
	std::cout << "Gamma factor: " << gamma << std::endl;
	std::cout << "y_cms = " << y_cms << std::endl;
}

void Particles::newEvent(bool first)
{
	if(first)
	{
		n[All] = n[Neg] = n[Pos] = 0;
		mean_pt[All] = mean_pt[Neg] = mean_pt[Pos] = 0.;
	}
	else
	{
		if(n[All] != 0)
		{
			mean_pt[All] /= n[All];
			histos->histMeanPt->Fill(mean_pt[All]);
			histos->histMeanCharge->Fill((((double)n[Pos]-n[Neg]))/n[All]);

			if(n[Neg] != 0)
			{
				mean_pt[Neg] /= n[Neg];
				histos->histMeanPtNeg->Fill(mean_pt[Neg]);
			}
			else
				mean_pt[Neg] = 0.;

			if(n[Pos] != 0)
			{
				mean_pt[Pos] /= n[Pos];
				histos->histMeanPtPos->Fill(mean_pt[Pos]);
			}
			else
				mean_pt[Pos] = 0.;
		}
		else
			mean_pt[All] = mean_pt[Neg] = mean_pt[Pos] = 0.;

		histos->histCharged->Fill(n[All]);
		histos->histChargedNeg->Fill(n[Neg]);
		histos->histChargedPos->Fill(n[Pos]);

		//std::cout << "Mult. in Particles class: " << (n[All]) << std::endl;
		n[All] = n[Neg] = n[Pos] = 0;
		mean_pt[All] = mean_pt[Neg] = mean_pt[Pos] = 0.;
	}
}

//Calculation of kinematic variables and filling kinematic distributions.
//Some of them are calculed second time. Quite inefficient but had no time and self-denial to optimize it.
void Particles::analyze(Particle *particle, const int ener)
{
	if(particle->isPositive())
		particle_charge = Pos;
	else
		particle_charge = Neg;

	px = particle->GetPx();
	py = particle->GetPy();
	pz = particle->GetPz();
	//std::cout << "gamma = " << gamma << " | gamma_beta_e = " << gamma_beta_e << std::endl;
	pt = TMath::Sqrt(py*py+px*px);
	p = TMath::Sqrt(px*px+py*py+pz*pz);
	E_pi = TMath::Sqrt(pion_mass*pion_mass + p*p);
	E_proton = TMath::Sqrt(proton_mass*proton_mass + p*p);
	pz_cms = gamma*pz - calc_gbE(E_pi); 

	angle = TMath::ATan2(py,px);
	theta = TMath::ATan2(pt,pz);
	y_pi = 0.5*TMath::Log((E_pi+pz)/(E_pi-pz));
	y_pi_cms = y_pi - y_cms;
	y_proton_cms = 0.5*TMath::Log((E_proton+pz)/(E_proton-pz)) - y_cms;
	theta_cms = TMath::ATan2(pt,pz_cms);

	eta  = -TMath::Log(TMath::Tan(theta/2.));

	eta_cms = -TMath::Log(TMath::Tan(theta_cms/2.));

	n[All]++;

	mean_pt[All] += pt;

	histos->histYpi->Fill(y_pi);
	histos->histYcms->Fill(y_pi_cms);
	histos->histYprotcms->Fill(y_proton_cms);
	histos->histEta->Fill(eta);
	histos->histEtacms->Fill(eta_cms);
	histos->histAngle->Fill(angle);
	histos->histTheta->Fill(theta);
	histos->histThetacms->Fill(theta_cms);
	histos->histPtWide->Fill(pt);
	histos->histPtAll->Fill(pt);
	histos->histPzAll->Fill(pz);
	histos->histPtot->Fill(p);
	histos->histPzcmsAll->Fill(pz_cms);
	histos->histPtVsYAll->Fill(y_pi_cms, pt);
	histos->histPtVsYprotAll->Fill(y_proton_cms, pt);
	histos->histPhiVsPtAll->Fill(angle, pt);
	histos->histDedx->Fill(p,particle->GetdEdx());
	histos->histDedxVtpc1->Fill(p,particle->GetdEdxVtpc1());
	histos->histDedxVtpc2->Fill(p,particle->GetdEdxVtpc2());
	histos->histDedxMtpc->Fill(p,particle->GetdEdxMtpc());

	if(particle->isPositive())
	{
		n[Pos]++;
		mean_pt[Pos] += pt;
		histos->histPtPos->Fill(pt);
		histos->histPzPos->Fill(pz);
		histos->histPzcmsPos->Fill(pz_cms);
		histos->histYpiPos->Fill(y_pi);
		histos->histYcmsPos->Fill(y_pi_cms);
		histos->histYprotcmsPos->Fill(y_proton_cms);
		histos->histEtaPos->Fill(eta);
		histos->histEtacmsPos->Fill(eta_cms);
		histos->histPtVsYPos->Fill(y_pi_cms, pt);
		histos->histPtVsYprotPos->Fill(y_proton_cms, pt);
		histos->histAnglePos->Fill(angle);
		histos->histThetaPos->Fill(theta);
		histos->histThetacmsPos->Fill(theta_cms);
		histos->histPhiVsPtPos->Fill(angle, pt);
		histos->histPartPopMatrixPos->Fill(p,pt,angle);
		histos->histDedxPos->Fill(p,particle->GetdEdx());
		histos->histDedxVtpc1Pos->Fill(p,particle->GetdEdxVtpc1());
		histos->histDedxVtpc2Pos->Fill(p,particle->GetdEdxVtpc2());
		histos->histDedxMtpcPos->Fill(p,particle->GetdEdxMtpc());
	}
	else
	{
		n[Neg]++;
		histos->histPartPopMatrixNeg->Fill(p,pt,angle);
		histos->histAngleNegNotrot->Fill(angle);	//bez obrotu (angle2)
		angle = mk_angle3(angle);			//obrot (angle2 -> angle3)
		histos->histAngleNeg->Fill(angle);		//z obrotem (angle3)
		histos->histThetaNeg->Fill(theta);
		histos->histThetacmsNeg->Fill(theta_cms);
		histos->histPtNeg->Fill(pt);
		histos->histPzNeg->Fill(pz);
		histos->histPzcmsNeg->Fill(pz_cms);
		mean_pt[Neg] += pt;
		histos->histYpiNeg->Fill(y_pi);
		histos->histYcmsNeg->Fill(y_pi_cms);
		histos->histYprotcmsNeg->Fill(y_proton_cms);
		histos->histEtaNeg->Fill(eta);
		histos->histEtacmsNeg->Fill(eta_cms);
		histos->histPtVsYNeg->Fill(y_pi_cms,pt);
		histos->histPtVsYprotNeg->Fill(y_proton_cms,pt);
		histos->histPhiVsPtNeg->Fill(angle, pt);
		histos->histDedxNeg->Fill(p,particle->GetdEdx());
		histos->histDedxVtpc1Neg->Fill(p,particle->GetdEdxVtpc1());
		histos->histDedxVtpc2Neg->Fill(p,particle->GetdEdxVtpc2());
		histos->histDedxMtpcNeg->Fill(p,particle->GetdEdxMtpc());
	}
}
