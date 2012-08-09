#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "RootWriter.h"
#include "Particle.h"
#include "TFile.h"

void Histos::init()
{
	histCharged = new TH1I("histCharged","Multiplicity of charged;N",25,0,25);
	histChargedNeg = new TH1I("histChargedNeg","Multiplicity of negatively charged;N",25,0,25);
	histChargedPos = new TH1I("histChargedPos","Multiplicity of positively charged;N",25,0,25);
	histMeanCharge = new TH1F("histMeanCharge","Mean charge;charge",23,-1.15,1.15);
	histAngle = new TH1F("histAngle","Azimuthal angle;#phi [rad]",50,-3.25,3.25);
	histAngleNegNotrot = new TH1F("histAngleNegNotrot","Azimuthal angle (notrot), neg.;#phi [rad]",50,-3.25,3.25);
	histAngleNeg = new TH1F("histAngleNeg","Azimuthal angle, neg.;#phi [rad]",50,-3.25,3.25);
	histAnglePos = new TH1F("histAnglePos","Azimuthal angle, pos.;#phi [rad]",50,-3.25,3.25);
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
	histEta = new TH1F("histEta","Pseudorapidity;#eta",100,-2,8);
	histEtaNeg = new TH1F("histEtaNeg","Pseudorapidity, negatively charged;#eta",100,-2,8);
	histEtaPos = new TH1F("histEtaPos","Pseudorapidity, positively charged;#eta",100,-2,8);
	histEtacms = new TH1F("histEtacms","Pseudorapidity;#eta",100,-4,6);
	histEtacmsNeg = new TH1F("histEtacmsNeg","Pseudorapidity, negatively charged;#eta",100,-4,6);
	histEtacmsPos = new TH1F("histEtacmsPos","Pseudorapidity, positively charged;#eta",100,-4,6);
	histPtAll = new TH1F("histPtAll","Transverse momentum;p_{T} [GeV/c]",150,0,3);
	histPtNeg = new TH1F("histPtNeg","Transverse momentum, neg.;p_{T} [GeV/c]",150,0,3);
	histPtPos = new TH1F("histPtPos","Transverse momentum, pos.;p_{T} [GeV/c]",150,0,3);
	histPzAll = new TH1F("histPzAll","Longitudinal momentum;p_{z} [GeV/c]",150,0,50);
	histPzNeg = new TH1F("histPzNeg","Longitudinal momentum, neg.;p_{z} [GeV/c]",150,0,50);
	histPzPos = new TH1F("histPzPos","Longitudinal momentum, pos.;p_{z} [GeV/c]",150,0,50);
	histPzcmsAll = new TH1F("histPzcmsAll","Longitudinal momentum, CMS;p_{z} [GeV/c]",200,-5,10);
	histPzcmsNeg = new TH1F("histPzcmsNeg","Longitudinal momentum, CMS, neg.;p_{z} [GeV/c]",200,-5,10);
	histPzcmsPos = new TH1F("histPzcmsPos","Longitudinal momentum, CMS, pos.;p_{z} [GeV/c]",200,-5,10);
	histMeanPt = new TH1F("histMeanPt","Mean transverse momentum (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtNeg = new TH1F("histMeanPtNeg","Mean transverse momentum, neg. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtPos = new TH1F("histMeanPtPos","Mean transverse momentum, pos. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histPtVsYAll = new TH2F("histPtVsYAll","Trans. momentum vs. rapidity; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYNeg = new TH2F("histPtVsYNeg","Trans. momentum vs. rapidity, neg.; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPtVsYPos = new TH2F("histPtVsYPos","Trans. momentum vs. rapidity, pos.; y^{*}_{#pi}; p_{T} [GeV/c]",100,-2,8,150,0,1.5);
	histPhiVsPtAll = new TH2F("histPhiVsPtAll","Az. angle vs. transverse momentum; #phi [rad]; p_{T} [GeV/c]",50,-3.2,3.2,150,0,1.5);
	histPhiVsPtPos = new TH2F("histPhiVsPtPos","Az. angle vs. transverse momentum, pos.; #phi [rad]; p_{T} [GeV/c]",50,-3.2,3.2,150,0,1.5);
	histPhiVsPtNeg = new TH2F("histPhiVsPtNeg","Az. angle vs. transverse momentum, neg.; #phi [rad]; p_{T} [GeV/c]",50,-3.2,3.2,150,0,1.5);
	histDyDphiAll = new TH2F("histDyDphiAll","#Deltay versus #Delta#phi;#Delta#phi [rad]; #Deltay_{#pi}",200,0,3.2,200,0,6.5);
	histDyDphiPos = new TH2F("histDyDphiPos","#Deltay versus #Delta#phi, pos.;#Delta#phi [rad]; #Deltay_{#pi}",200,0,3.2,200,0,6.5);
	histDyDphiNeg = new TH2F("histDyDphiNeg","#Deltay versus #Delta#phi, neg.;#Delta#phi [rad]; #Deltay_{#pi}",200,0,3.2,200,0,6.5);
	histDyDphiUnlike = new TH2F("histDyDphiUnlike","#Deltay versus #Delta#phi, unlike-sign;#Delta#phi [rad]; #Deltay_{#pi}",200,0,3.2,200,0,6.5);
	histDetaDphiAll = new TH2F("histDetaDphiAll","#Delta#eta versus #Delta#phi;#Delta#phi [rad];#Delta#eta",200,0,3.2,200,0,6);
	histDetaDphiPos = new TH2F("histDetaDphiPos","#Delta#eta versus #Delta#phi, pos.;#Delta#phi [rad];#Delta#eta",200,0,3.2,200,0,6);
	histDetaDphiNeg = new TH2F("histDetaDphiNeg","#Delta#eta versus #Delta#phi, neg.;#Delta#phi [rad];#Delta#eta",200,0,3.2,200,0,6);
	histDetaDphiUnlike = new TH2F("histDetaDphiUnlike","#Delta#eta versus #Delta#phi, unlike-sign;#Delta#phi [rad];#Delta#eta",200,0,3.2,200,0,6);
	histPartPopMatrixPos = new TH3I("histPartPopMatrixPos","Particle population matrix, pos. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,149,40,0,2,36,-TMath::Pi(),TMath::Pi());
	histPartPopMatrixNeg = new TH3I("histPartPopMatrixNeg","Particle population matrix, neg. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,149,40,0,2,36,-TMath::Pi(),TMath::Pi());
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
	histEta->Write();
	histEtaNeg->Write();
	histEtaPos->Write();
	histEtacms->Write();
	histEtacmsNeg->Write();
	histEtacmsPos->Write();
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
	histPtVsYAll->Write();
	histPtVsYPos->Write();
	histPtVsYNeg->Write();
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
	histPartPopMatrixPos->Write();
	histPartPopMatrixNeg->Write();
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
	delete	histEta;
	delete	histEtaNeg;
	delete	histEtaPos;
	delete	histEtacms;
	delete	histEtacmsNeg;
	delete	histEtacmsPos;
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
	delete	histPtVsYAll;
	delete	histPtVsYPos;
	delete	histPtVsYNeg;
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
	delete histPartPopMatrixPos;	
	delete histPartPopMatrixNeg;	
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

void Particles::init(Histos *histograms, const float ener)
{
	histos = histograms;
	angle = 0.;
	theta = theta_cms = 0.;
	y_pi = y_cms = y_pi_cms = 0.;
	eta = eta_cms = 0.;
	pz_cms = pt = 0.;
	particle_charge = All;
	for(int ch=0; ch<3; ch++)
	{
		n[ch] = 0;
		mean_pt[ch] = 0.;
	}

	beta = calc_beta(ener);
	gamma = calc_gamma(ener);
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

		n[All] = n[Neg] = n[Pos] = 0;
		mean_pt[All] = mean_pt[Neg] = mean_pt[Pos] = 0.;
	}
}

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
	//E = TMath::Sqrt(0.1396*0.1396 + p*p);
	E = TMath::Sqrt(pion_mass*pion_mass + p*p);
	pz_cms = gamma*pz - calc_gbE(E); 

	angle = TMath::ATan2(py,px);
	theta = TMath::ATan2(pt,pz);
	y_pi = 0.5*TMath::Log((E+pz)/(E-pz));
	y_pi_cms = y_pi - y_cms;
	theta_cms = TMath::ATan2(pt,pz_cms);

	//std::cout << "Theta = " << theta << " | Theta_cms = " << theta_cms << std::endl;

	eta  = -TMath::Log(TMath::Tan(theta/2.));

	eta_cms = -TMath::Log(TMath::Tan(theta_cms/2.));

	n[All]++;

	mean_pt[All] += pt;

	histos->histYpi->Fill(y_pi);
	histos->histYcms->Fill(y_pi_cms);
	histos->histEta->Fill(eta);
	histos->histEtacms->Fill(eta_cms);
	histos->histAngle->Fill(angle);
	histos->histTheta->Fill(theta);
	histos->histThetacms->Fill(theta_cms);
	histos->histPtAll->Fill(pt);
	histos->histPzAll->Fill(pz);
	histos->histPzcmsAll->Fill(pz_cms);
	histos->histPtVsYAll->Fill(y_pi_cms, pt);
	histos->histPhiVsPtAll->Fill(angle, pt);

	if(particle->isPositive())
	{
		n[Pos]++;
		mean_pt[Pos] += pt;
		histos->histPtPos->Fill(pt);
		histos->histPzPos->Fill(pz);
		histos->histPzcmsPos->Fill(pz_cms);
		histos->histYpiPos->Fill(y_pi);
		histos->histYcmsPos->Fill(y_pi_cms);
		histos->histEtaPos->Fill(eta);
		histos->histEtacmsPos->Fill(eta_cms);
		histos->histPtVsYPos->Fill(y_pi_cms, pt);
		histos->histAnglePos->Fill(angle);
		histos->histThetaPos->Fill(theta);
		histos->histThetacmsPos->Fill(theta_cms);
		histos->histPhiVsPtPos->Fill(angle, pt);
		histos->histPartPopMatrixPos->Fill(p,pt,angle);
	}
	else
	{
		n[Neg]++;
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
		histos->histEtaNeg->Fill(eta);
		histos->histEtacmsNeg->Fill(eta_cms);
		histos->histPtVsYNeg->Fill(y_pi_cms,pt);
		histos->histPhiVsPtNeg->Fill(angle, pt);
		histos->histPartPopMatrixNeg->Fill(p,pt,angle);
	}
}

using namespace std;

map<int,int> getDistro(TString root_filename)
{
	TFile *rootfile = new TFile(root_filename);
	map<int,int> distribution;

	TH1I* histMult = (TH1I*)rootfile->Get("histCharged");
	cout << "\t\tdone" << endl;

	for(int i=0; i<histMult->GetNbinsX(); i++)
		distribution.insert(pair<int,int> (i,histMult->GetBinContent(i+1)));

	cout << distribution[0] << " events with 0 particles" << endl;

	return distribution;
}
