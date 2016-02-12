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

	histCharged = new TH1I("histCharged","Multiplicity of charged;N",20,0,20);
	histChargedNeg = new TH1I("histChargedNeg","Multiplicity of negatively charged;N",20,0,20);
	histChargedPos = new TH1I("histChargedPos","Multiplicity of positively charged;N",20,0,20);
	histMeanCharge = new TH1F("histMeanCharge","Mean charge;charge",23,-1.15,1.15);
	histAngle = new TH1F("histAngle","Azimuthal angle;#phi [rad]",50,-TMath::Pi(), TMath::Pi());
	histAngleNegNotrot = new TH1F("histAngleNegNotrot","Azimuthal angle (notrot), neg.;#phi [rad]",50,-TMath::Pi(), TMath::Pi());
	histAngleNeg = new TH1F("histAngleNeg","Azimuthal angle, neg.;#phi [rad]",50,-TMath::Pi(),TMath::Pi());
	histAnglePos = new TH1F("histAnglePos","Azimuthal angle, pos.;#phi [rad]",50,-TMath::Pi(),TMath::Pi());
	histTheta = new TH1F("histTheta","Polar angle;#theta [rad]", 50, -1, 1);
	histThetaNeg = new TH1F("histThetaNeg","Polar angle, neg.;#theta [rad]", 50, -1, 1);
	histThetaPos = new TH1F("histThetaPos","Polar angle, pos.;#theta [rad]", 50, -1, 1);
	histYreal = new TH1F("histYreal","Rapidity with real mass;y",200,-8,8);
	histYrealNeg = new TH1F("histYrealNeg","Rapidity with real mass, neg.;y",200,-8,8);
	histYrealPos = new TH1F("histYrealPos","Rapidity with real mass, pos.;y",200,-8,8);
	histYpi = new TH1F("histYpi","Rapidity with #pi mass;y_{#pi}",200,-10,10);
	histYpiNeg = new TH1F("histYpiNeg","Rapidity with #pi mass, neg.;y_{#pi}",200,-10,10);
	histYpiPos = new TH1F("histYpiPos","Rapidity with #pi mass, pos.;y_{#pi}",200,-10,10);
	histEta = new TH1F("histEta","Pseudorapidity;#eta",200,-8,8);
	histEtaNeg = new TH1F("histEtaNeg","Pseudorapidity, negatively charged;#eta",200,-8,8);
	histEtaPos = new TH1F("histEtaPos","Pseudorapidity, positively charged;#eta",200,-8,8);
	histPtWide = new TH1F("histPtWide","Transverse momentum (wide);p_{T} [GeV/c]",750,0,150);
	histPtAll = new TH1F("histPtAll","Transverse momentum;p_{T} [GeV/c]",150,0,3);
	histPtNeg = new TH1F("histPtNeg","Transverse momentum, neg.;p_{T} [GeV/c]",150,0,3);
	histPtPos = new TH1F("histPtPos","Transverse momentum, pos.;p_{T} [GeV/c]",150,0,3);
	histPzAll = new TH1F("histPzAll","Longitudinal momentum;p_{z} [GeV/c]",150,0,50);
	histPzNeg = new TH1F("histPzNeg","Longitudinal momentum, neg.;p_{z} [GeV/c]",150,0,50);
	histPzPos = new TH1F("histPzPos","Longitudinal momentum, pos.;p_{z} [GeV/c]",150,0,50);
	histPtot = new TH1F("histPtot","Total momentum;p_{tot} [GeV/c]",300, 0, 150);
	histEtot = new TH1F("histEtot","Total energy;E_{tot} [GeV]", 300, 0, 150);
	histEtotCALM = new TH1F("histEtotCALM","Total energy of #pi, K, p, n, #Lambda;E_{tot} [GeV]", 400, 0, 200);

	histMeanPt = new TH1F("histMeanPt","Mean transverse momentum (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtNeg = new TH1F("histMeanPtNeg","Mean transverse momentum, neg. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histMeanPtPos = new TH1F("histMeanPtPos","Mean transverse momentum, pos. (ev. without 0 mult.);M(p_{T}) [GeV/c]",100,0,2);
	histPtVsYAll = new TH2F("histPtVsYAll","Trans. momentum vs. rapidity; y^{*}_{#pi}; p_{T} [GeV/c]",100,-8,8,200,0,1.5);
	histPtVsYNeg = new TH2F("histPtVsYNeg","Trans. momentum vs. rapidity, neg.; y^{*}_{#pi}; p_{T} [GeV/c]",200,-8,8,150,0,1.5);
	histPtVsYPos = new TH2F("histPtVsYPos","Trans. momentum vs. rapidity, pos.; y^{*}_{#pi}; p_{T} [GeV/c]",200,-8,8,150,0,1.5);
	histPtVsYprotAll = new TH2F("histPtVsYprotAll","Trans. momentum vs. rapidity; y^{*}_{prot}; p_{T} [GeV/c]",200,-8,8,150,0,1.5);
	histPtVsYprotNeg = new TH2F("histPtVsYprotNeg","Trans. momentum vs. rapidity, neg.; y^{*}_{prot}; p_{T} [GeV/c]",200,-8,8,150,0,1.5);
	histPtVsYprotPos = new TH2F("histPtVsYprotPos","Trans. momentum vs. rapidity, pos.; y^{*}_{prot}; p_{T} [GeV/c]",200,-8,8,150,0,1.5);
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

	histPartPopMatrixPos = new TH3I("histPartPopMatrixPos","Particle population matrix, pos. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,150,40,0,2,36,-TMath::Pi(),TMath::Pi());
	histPartPopMatrixNeg = new TH3I("histPartPopMatrixNeg","Particle population matrix, neg. charged; p_{tot} [GeV/c]; p_{T} [GeV/c]; #phi [rad]",150,0,150,40,0,2,36,-TMath::Pi(),TMath::Pi());
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
	histYreal->Write();
	histYrealPos->Write();
	histYrealNeg->Write();
	histYpi->Write();
	histYpiPos->Write();
	histYpiNeg->Write();
	histEta->Write();
	histEtaNeg->Write();
	histEtaPos->Write();
	histPtWide->Write();
	histPtAll->Write();
	histPtPos->Write();
	histPtNeg->Write();
	histPzAll->Write();
	histPzPos->Write();
	histPzNeg->Write();
	histMeanPt->Write();
	histMeanPtNeg->Write();
	histMeanPtPos->Write();
	histPtot->Write();
	histEtot->Write();
	histEtotCALM->Write();
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
	delete 	histYreal;
	delete 	histYrealPos;
	delete 	histYrealNeg;
	delete 	histYpi;
	delete 	histYpiPos;
	delete 	histYpiNeg;
	delete	histEta;
	delete	histEtaNeg;
	delete	histEtaPos;
	delete 	histPtWide;
	delete 	histPtAll;
	delete 	histPtPos;
	delete 	histPtNeg;
	delete 	histPzAll;
	delete 	histPzPos;
	delete 	histPzNeg;
	delete	histMeanPt;
	delete	histMeanPtNeg;
	delete	histMeanPtPos;
	delete	histPtot;
	delete	histEtot;
	delete	histEtotCALM;
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

void Particles::init(Histos *histograms, const float momentum)
{
	histos = histograms;
	angle = 0.;
	theta = 0.;
	y = 0.;
	E_real = 0.;
	y_pi = 0.;
	eta = 0.;
	pt = 0.;
	particle_charge = All;
	for(int ch=0; ch<3; ch++)
	{
		n[ch] = 0;
		mean_pt[ch] = 0.;
	}

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
	mass = particle->GetMass();
	pt = TMath::Sqrt(py*py+px*px);
	p = TMath::Sqrt(px*px+py*py+pz*pz);
	E_real = TMath::Sqrt(mass*mass + p*p);
	E_pi = TMath::Sqrt(pion_mass*pion_mass + p*p);
	E_proton = TMath::Sqrt(proton_mass*proton_mass + p*p);

	angle = TMath::ATan2(py,px);
	theta = TMath::ATan2(pt,pz);
	y = 0.5*TMath::Log((E_real+pz)/(E_real-pz));
	y_pi = 0.5*TMath::Log((E_pi+pz)/(E_pi-pz));

	eta  = -TMath::Log(TMath::Tan(theta/2.));

	n[All]++;

	mean_pt[All] += pt;

	histos->histYreal->Fill(y);
	histos->histYpi->Fill(y_pi);
	histos->histEta->Fill(eta);
	histos->histAngle->Fill(angle);
	histos->histTheta->Fill(theta);
	histos->histPtWide->Fill(pt);
	histos->histPtAll->Fill(pt);
	histos->histPzAll->Fill(pz);
	histos->histPtot->Fill(p);
	histos->histEtot->Fill(E_real);
	histos->histPhiVsPtAll->Fill(angle, pt);

	if(particle->isPositive())
	{
		n[Pos]++;
		mean_pt[Pos] += pt;
		histos->histPtPos->Fill(pt);
		histos->histPzPos->Fill(pz);
		histos->histYrealPos->Fill(y);
		histos->histYpiPos->Fill(y_pi);
		histos->histEtaPos->Fill(eta);
		histos->histAnglePos->Fill(angle);
		histos->histThetaPos->Fill(theta);
		histos->histPhiVsPtPos->Fill(angle, pt);
		histos->histPartPopMatrixPos->Fill(p,pt,angle);
	}
	else
	{
		n[Neg]++;
		histos->histPartPopMatrixNeg->Fill(p,pt,angle);
		histos->histAngleNegNotrot->Fill(angle);	//bez obrotu (angle2)
		angle = mk_angle3(angle);			//obrot (angle2 -> angle3)
		histos->histAngleNeg->Fill(angle);		//z obrotem (angle3)
		histos->histThetaNeg->Fill(theta);
		histos->histPtNeg->Fill(pt);
		histos->histPzNeg->Fill(pz);
		mean_pt[Neg] += pt;
		histos->histYrealNeg->Fill(y);
		histos->histYpiNeg->Fill(y_pi);
		histos->histEtaNeg->Fill(eta);
		histos->histPhiVsPtNeg->Fill(angle, pt);
	}
}
