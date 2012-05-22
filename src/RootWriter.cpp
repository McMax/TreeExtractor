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
	histYpi = new TH1F("histYpi","Rapidity with #pi mass;y_{#pi}",100,-2,8);
	histYpiNeg = new TH1F("histYpiNeg","Rapidity with #pi mass, neg.;y_{#pi}",100,-2,8);
	histYpiPos = new TH1F("histYpiPos","Rapidity with #pi mass, pos.;y_{#pi}",100,-2,8);
	histYcms = new TH1F("histYcms","CMS rapidity with #pi mass;y_{#pi}",100,-5,5);
	histYcmsNeg = new TH1F("histYcmsNeg","CMS rapidity with #pi mass, neg.;y_{#pi}",100,-5,5);
	histYcmsPos = new TH1F("histYcmsPos","CMS rapidity with #pi mass, pos.;y_{#pi}",100,-5,5);
	histEta = new TH1F("histEta","Pseudorapidity;#eta",100,-2,8);
	histEtaNeg = new TH1F("histEtaNeg","Pseudorapidity, negatively charged;#eta",100,-2,8);
	histEtaPos = new TH1F("histEtaPos","Pseudorapidity, positively charged;#eta",100,-2,8);
	histPtAll = new TH1F("histPtAll","Transverse momentum;p_{T} [GeV/c]",150,0,3);
	histPtNeg = new TH1F("histPtNeg","Transverse momentum, neg.;p_{T} [GeV/c]",150,0,3);
	histPtPos = new TH1F("histPtPos","Transverse momentum, pos.;p_{T} [GeV/c]",150,0,3);
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
	histYpi->Write();
	histYpiPos->Write();
	histYpiNeg->Write();
	histYcms->Write();
	histYcmsPos->Write();
	histYcmsNeg->Write();
	histEta->Write();
	histEtaNeg->Write();
	histEtaPos->Write();
	histPtAll->Write();
	histPtPos->Write();
	histPtNeg->Write();
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
	delete 	histYpi;
	delete 	histYpiPos;
	delete 	histYpiNeg;
	delete 	histYcms;
	delete 	histYcmsPos;
	delete 	histYcmsNeg;
	delete	histEta;
	delete	histEtaNeg;
	delete	histEtaPos;
	delete 	histPtAll;
	delete 	histPtPos;
	delete 	histPtNeg;
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
}	

float mk_angle3(float  angle1)
{
	if(angle1 < 0.)		//Konwersja angle3
		angle1 += TMath::Pi();
	else
		angle1 -= TMath::Pi();

	return angle1;
}

void Particles::init(Histos *histograms)
{
	histos = histograms;
	angle = 0.;
	theta = 0.;
	y_pi = y_cms = y_beam = 0.;
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
	pt = TMath::Sqrt(py*py+px*px);
	p = TMath::Sqrt(px*px+py*py+pz*pz);
	E = TMath::Sqrt(0.1396*0.1396 + p*p);

	angle = TMath::ATan2(py,px);
	theta = TMath::ATan(pt/pz);
	y_pi = 0.5*TMath::Log((E+pz)/(E-pz));
	y_cms = lab2cms(y_pi,ener);

	if(theta > 0)
		eta = -TMath::Log(0.5*theta);
	else
		eta = 0.;

	n[All]++;

	mean_pt[All] += pt;

	histos->histYpi->Fill(y_pi);
	histos->histYcms->Fill(y_cms);
	histos->histEta->Fill(eta);
	histos->histAngle->Fill(angle);
	histos->histPtAll->Fill(pt);
	histos->histPtVsYAll->Fill(y_cms, pt);
	histos->histPhiVsPtAll->Fill(angle, pt);

	if(particle->isPositive())
	{
		n[Pos]++;
		mean_pt[Pos] += pt;
		histos->histPtPos->Fill(pt);
		histos->histYpiPos->Fill(y_pi);
		histos->histYcmsPos->Fill(y_cms);
		histos->histEtaPos->Fill(eta);
		histos->histPtVsYPos->Fill(y_cms, pt);
		histos->histAnglePos->Fill(angle);
		histos->histPhiVsPtPos->Fill(angle, pt);
	}
	else
	{
		n[Neg]++;
		histos->histAngleNegNotrot->Fill(angle);	//bez obrotu (angle2)
		angle = mk_angle3(angle);			//obrot (angle2 -> angle3)
		histos->histAngleNeg->Fill(angle);		//z obrotem (angle3)
		histos->histPtNeg->Fill(pt);
		mean_pt[Neg] += pt;
		histos->histYpiNeg->Fill(y_pi);
		histos->histYcmsNeg->Fill(y_cms);
		histos->histEtaNeg->Fill(eta);
		histos->histPtVsYNeg->Fill(y_cms,pt);
		histos->histPhiVsPtNeg->Fill(angle, pt);
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

double lab2cms(const float rapidity, const float ener)
{
	static float rapidity_beam;
	rapidity_beam = TMath::ATanH(TMath::Sqrt(1-TMath::Power(0.9315/ener,2)));
	return (rapidity - rapidity_beam/2);
}
