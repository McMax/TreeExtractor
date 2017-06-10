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
#include "ClusterGraphs.h"
#include "Event.h"
#include "Particle.h"

using namespace std;

void mainanalyze(TTree *particletree, const TString system, const float beam_momentum, const TString output_filename="Extracted_distributions.root")
{
	//ofstream debugfile("Debug.txt");
	float	angle, angle3,
		p1, p2,
		pt1, pt2,
		pz_cms1, pz_cms2,
		E1, E2,
		E_prot,
		inv_mass,
		gbE1, gbE2,
		theta1, theta2,
		y1, y2,
		y_prot_cms,
		eta1, eta2,
		angle_j,
		angle_diff,
		y_diff,
		eta_diff;
	
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
	ClusterGraphs clustergraphs;
	TFile *root_output_file;

	histos.init(beam_momentum);
	particles.init(&histos, system, beam_momentum);
	particles.newEvent(true);
	root_output_file = new TFile(output_filename,"recreate");
	clustergraphs.setOtherHistFile(root_output_file);
	AdditionalInfo ai;
//End of preparation

	cout << "Writing events" << endl;

//Loop over events
	//for(ev=0; ev<treeNentries; ++ev)
	for(ev=0; ev<10000; ++ev)
	{
		particletree->GetEntry(ev);
		
		n[Neg] = n[All] = n[Pos] = 0;

		//debugfile << ev << "\t" << event->GetNpa() << endl;

//Loop over the first particle in two-particle pair
		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);

//Calculating kinematics for Delta-eta Delta-phi
			pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
			p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
			E1 = TMath::Sqrt(TMath::Power(pion_mass,2)+p1*p1);
			E_prot = TMath::Sqrt(proton_mass*proton_mass+p1*p1);
			y_prot_cms = 0.5*TMath::Log((E_prot+particleA->GetPz())/(E_prot-particleA->GetPz())) - particles.y_cms;
			v1.SetPxPyPzE(particleA->GetPx(),particleA->GetPy(),particleA->GetPz(),E1);

			y1 = 0.5*TMath::Log((E1+particleA->GetPz())/(E1-particleA->GetPz())) - particles.y_cms;
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

					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));
					E2 = TMath::Sqrt(TMath::Power(pion_mass,2)+p2*p2);
					E_prot = TMath::Sqrt(proton_mass*proton_mass+p2*p2);
					y_prot_cms = 0.5*TMath::Log((E_prot+particleB->GetPz())/(E_prot-particleB->GetPz())) - particles.y_cms;
					v2.SetPxPyPzE(particleB->GetPx(),particleB->GetPy(),particleB->GetPz(),E2);

					v = v1 + v2;
					inv_mass = v.M();
					
					/*
					if(inv_mass < 0.285) //GeV dipion (280 MeV) + Coulomb interactions (5 MeV)
						continue;
						*/

					histos.histInvMass->Fill(inv_mass);

					//cout << "E1 = " << E1 << " | E2 = " << E2 << endl;

//Gamma, beta, energy. All needed to shift particles rapidity to CMS. Calculation done in RootWriter.h. I agree, this is strange place of putting such code in header.
					gbE1 = particles.calc_gbE(E1);
					gbE2 = particles.calc_gbE(E2);

					//cout << "Beta factor: " << particles.beta << endl;
					//cout << "Gamma factor: " << particles.gamma << endl;
					//cout << "gamma*beta*E1: " << gbE1 << " | gamma*beta*E2: " << gbE2 << endl;

					pz_cms1 = particles.gamma*particleA->GetPz() - gbE1;
					pz_cms2 = particles.gamma*particleB->GetPz() - gbE2;

					//debugfile << "pz_cms1 = " << pz_cms1 << " | pz_cms2 = " << pz_cms2 << endl;

					y2 = 0.5*TMath::Log((E2+particleB->GetPz())/(E2-particleB->GetPz())) - particles.y_cms;

					angle_j = TMath::ATan2(particleB->GetPy(), particleB->GetPx());

					//cout << "y1 = " << y1 << " | y2 = " << y2 << endl;

					theta1 = TMath::Abs(TMath::ATan2(pt1,pz_cms1));
					theta2 = TMath::Abs(TMath::ATan2(pt2,pz_cms2));

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
					eta_diff = TMath::Abs(eta1-eta2);

					if((angle_diff < 0.00552) && (eta_diff < 0.00792))
					{
						ai.pz_cms1 = pz_cms1;
						ai.pz_cms2 = pz_cms2;
						ai.eta_cms1 = eta1;
						ai.eta_cms2 = eta2;
						ai.phi1 = angle;
						ai.phi2 = angle_j;
						ai.deta = eta_diff;
						ai.dphi = angle_diff;
						clustergraphs.addGraph(ev, particleA, particleB, ai);
					}

//Filling dydphi and detadphi histograms
					histos.histDyDphiAll->Fill(angle_diff, y_diff);
					histos.histDetaDphiAll->Fill(angle_diff, eta_diff);

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

		if((event->GetNpa()!=n[All]))
				cerr << "Event: " << ev << " ParticleTree: " << (event->GetNpa()) << " Extraction: " << (n[All]) << endl;

		//debugfile << "----------" << endl;

		//cout << "\rEvent " << ev;
		if(!(ev%500))
			cout << "Event " << ev << endl;

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
	clustergraphs.closeFile();
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

