#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cstdio>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "Prefifi.h"
#include "RootWriter.h"
#include "Event.h"
#include "Particle.h"

using namespace std;

void mainanalyze(TTree *particletree, const float beam_momentum, const TString output_filename="Extracted_distributions.root")
{
	//ofstream debugfile("Debug.txt");
	float	angle, angle3,
		p1, p2,
		pt1, pt2,
		pz1, pz2,
		E1, E2,
		E_prot,
		inv_mass,
		y1, y2,
		eta1, eta2,
		theta1, theta2,
		angle_j,
		angle_diff,
		y_diff,
		eta_diff,
		mass,
		Etot;
	
	bool	positive,
		positive_j;

	int	n[3];
	unsigned int all_particles=0;
	UInt_t	i,j;

	TLorentzVector v1, v2, v;

	unsigned correlations = 0, pos_correlations = 0, neg_correlations = 0, all_correlations = 0, unlike_correlations = 0;

//Preparation of ParticleTree output file
	Event *event = new Event();
	Particle *particleA, *particleB;
	particletree->SetBranchAddress("event",&event);
	Long64_t treeNentries = particletree->GetEntries();
	cout << "Number of events: " << treeNentries << endl;
	Long64_t ev;

	Particles particles;
	Histos histos;
	TFile *root_output_file;

	histos.init(beam_momentum);
	particles.init(&histos, beam_momentum);
	particles.newEvent(true);
	root_output_file = new TFile(output_filename,"recreate");
//End of preparation

	cout << "Writing events" << endl;

//Loop over events
	for(ev=0; ev<treeNentries; ++ev)
	{
		particletree->GetEntry(ev);
		
		n[Neg] = n[All] = n[Pos] = 0;
		Etot = 0.;

		//debugfile << ev << "\t" << event->GetNpa() << endl;

//Loop over the first particle in two-particle pair
		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);

//Calculating kinematics for Delta-eta Delta-phi
			pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
			p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
			E1 = TMath::Sqrt(TMath::Power(pion_mass,2)+p1*p1);

			mass = particleA->GetMass();
			pz1 = particleA->GetPz();

//Adding the energy to total energy of pi. K, n, p, lambdas
			Etot += TMath::Sqrt(mass*mass + p1*p1);

			E_prot = TMath::Sqrt(proton_mass*proton_mass+p1*p1);
			v1.SetPxPyPzE(particleA->GetPx(),particleA->GetPy(),pz1,E1);

			y1 = 0.5*TMath::Log((E1+pz1)/(E1-pz1));
			angle = TMath::ATan2(particleA->GetPy(), particleA->GetPx());

//Calculating the rest of variables (RootWriter.cpp)
			particles.analyze(particleA,beam_momentum);

			//debugfile << i << ": " <<  particleA->GetPx() << "\t" << particleA->GetPy() << "\t" << particleA->GetPz() << endl;

			positive = particleA->isPositive();

			//debugfile << (positive ? "1 " : "-1 ") << angle << endl;
			if(!positive)
				angle3 = mk_angle3(angle);

			if(event->GetNpa() > 1)
			{
//Loop over the second particle in two-particle pair
				for(j=i+1; j<event->GetNpa(); ++j)
				{
					particleB = event->GetParticle(j);

					//cout << "Particle A [" << i << "]:" << particleA << " particleB[" << j << "]:" << particleB << endl;
					//cout << "Particle A: px=" << particleA->GetPx() << " py=" << particleA->GetPy() << " pz=" << particleA->GetPz() << endl;
					//cout << "Particle B: px=" << particleB->GetPx() << " py=" << particleB->GetPy() << " pz=" << particleB->GetPz() << endl;

					pt2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2));
					pz2 = particleB->GetPz();
					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));
					E2 = TMath::Sqrt(TMath::Power(pion_mass,2)+p2*p2);
					E_prot = TMath::Sqrt(proton_mass*proton_mass+p2*p2);
					v2.SetPxPyPzE(particleB->GetPx(),particleB->GetPy(),pz2,E2);

					v = v1 + v2;
					inv_mass = v.M();
					
					/*
					if(inv_mass < 0.285) //GeV dipion (280 MeV) + Coulomb interactions (5 MeV)
						continue;
						*/

					histos.histInvMass->Fill(inv_mass);

					//cout << "E1 = " << E1 << " | E2 = " << E2 << endl;


					y2 = 0.5*TMath::Log((E2+pz2)/(E2-pz2));

					angle_j = TMath::ATan2(particleB->GetPy(), particleB->GetPx());

					//cout << "y1 = " << y1 << " | y2 = " << y2 << endl;

					theta1 = TMath::Abs(TMath::ATan2(pt1,pz1));
					theta2 = TMath::Abs(TMath::ATan2(pt2,pz2));

					//debugfile << "theta1 = " << theta1 << " | theta2 = " << theta2 << endl;

					eta1 = -TMath::Log(TMath::Tan(0.5*theta1));
					eta2 = -TMath::Log(TMath::Tan(0.5*theta2));

					//debugfile << "eta1 = " << eta1 << " | eta2 = " << eta2 << endl;
					//cout << "angle1 = " << angle << " | angle2 = " << angle_j << endl;

					positive_j = particleB->isPositive();

//Calculating differences in azimuthal angle and rapidity
					if((angle_diff = TMath::Abs(angle-angle_j)) > TMath::Pi())
						angle_diff = 2*TMath::Pi()-angle_diff;
					y_diff = TMath::Abs(y1-y2);

//Filling dydphi and detadphi histograms
					histos.histDyDphiAll->Fill(angle_diff, (y_diff = TMath::Abs(y1-y2)));
					histos.histDetaDphiAll->Fill(angle_diff, (eta_diff = TMath::Abs(eta1-eta2)));

//Filling with automatic reflection (ALICE style)
					Fill4Times(histos.histDetaDphiAllReflected, eta_diff, angle_diff);

					//debugfile << "deta=" << eta_diff << endl;
					
					++all_correlations;

//The same for different charge combinations

					if((positive_j == true) && (positive == true))
					{
						++correlations;
						++pos_correlations;

						histos.histDyDphiPos->Fill(angle_diff, y_diff);
						histos.histDetaDphiPos->Fill(angle_diff, eta_diff);
						Fill4Times(histos.histDetaDphiPosReflected, eta_diff, angle_diff);
					}
					else if((positive_j == false) && (positive == false))
					{
						++correlations;
						++neg_correlations;
						histos.histDyDphiNeg->Fill(angle_diff, y_diff);
						histos.histDetaDphiNeg->Fill(angle_diff, eta_diff);
						Fill4Times(histos.histDetaDphiNegReflected, eta_diff, angle_diff);
					}
					else
					{
						++unlike_correlations;
						histos.histDyDphiUnlike->Fill(angle_diff, y_diff);
						histos.histDetaDphiUnlike->Fill(angle_diff, eta_diff);
						Fill4Times(histos.histDetaDphiUnlikeReflected, eta_diff, angle_diff);
					}
				}
			}

			all_particles++;
			n[All]++;

			if(positive)
				n[Pos]++;
			else
				n[Neg]++;
		}	

		//debugfile << "----------" << endl;

		//cout << "\rEvent " << ev;
		if(!(ev%5000))
			cout << "Event " << ev << endl;

		histos.histEtotCALM->Fill(Etot);

		particles.newEvent();
	}

	//event--;

	cout << "All correlations: " << all_correlations << endl;
	cout << "Like-sign correlations: " << correlations << endl;
	cout << "Positive correlations: " << pos_correlations << endl;
	cout << "Negative correlations: " << neg_correlations << endl;
	cout << "=======================" << endl << "All particles: " << all_particles << ", all events: " << ev << endl;
	cout << "Mean multiplicity: " << (((double)all_particles)/ev) << endl;
	//debugfile << "All correlations: " << all_correlations << endl;
	//debugfile << "\nLike-sign correlations: " << correlations << endl;
	//debugfile << "Positive correlations: " << pos_correlations << endl;
	//debugfile << "Negative correlations: " << neg_correlations << endl;

	//debugfile.close();

	//	histos.histCharged->ResetStats();
	//	histos.histChargedNeg->ResetStats();
	//	histos.histChargedPos->ResetStats();
	root_output_file->cd();
	histos.write();
	histos.clear();
	root_output_file->Close();
}

void Fill4Times(TH2F* hist, const float deta, const float dphi)
{
	hist->Fill(dphi,deta);
	hist->Fill(dphi,-deta);

	if((-1*dphi) < (-0.5*TMath::Pi()))	//If phi's reflection is smaller than -Pi/2, reflect it
	{
		hist->Fill((2*TMath::Pi()-dphi),deta);
		hist->Fill((2*TMath::Pi()-dphi),-deta);
	}
	else
	{

		hist->Fill(-dphi,deta);
		hist->Fill(-dphi,-deta);
	}
}

