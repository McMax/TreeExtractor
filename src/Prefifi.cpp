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
	cout << "Beta calculated for nucleon mass: " << nucleon_mass << " GeV/c^2" << endl;

	float 	phi[3],
		phiSq[3],
		angle,
		angle3,

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
	UInt_t	i,j;

	TLorentzVector v1, v2, v;

	unsigned correlations = 0, pos_correlations = 0, neg_correlations = 0, all_correlations = 0, unlike_correlations = 0;

	Event *event = new Event();
	Particle *particleA, *particleB;
	particletree->SetBranchAddress("event",&event);
	Long64_t treeNentries = particletree->GetEntries();
	cout << "Number of events: " << treeNentries << endl;
	Long64_t ev;

	Particles particles;
	Histos histos;
	TFile *root_output_file;

//	to dE/dx of particles in (deta,dphi) < (0.5,0.5)
	std::set<UInt_t> unique_particles_y;
	std::set<UInt_t> unique_particles_eta;
	std::set<UInt_t> unique_particles_y_025;
	std::set<UInt_t> unique_particles_eta_025;

	histos.init(beam_momentum);
	particles.init(&histos, beam_momentum);
	particles.newEvent(true);
	root_output_file = new TFile(output_filename,"recreate");

	cout << "Writing events" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		particletree->GetEntry(ev);
		
		phi[Neg] = phi[All] = phi[Pos]= 0.;
		phiSq[Neg] = phiSq[All] = phiSq[Pos] = 0.;
		n[Neg] = n[All] = n[Pos] = 0;


		//debugfile << ev << "\t" << event->GetNpa() << endl;

		unique_particles_y.clear();
		unique_particles_eta.clear();
		unique_particles_y_025.clear();
		unique_particles_eta_025.clear();

		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);

			//if((TMath::Abs(particleA->GetBx()) > 4) || (TMath::Abs(particleA->GetBy()) > 2))
			//	continue;
			pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
			p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
			E1 = TMath::Sqrt(TMath::Power(pion_mass,2)+p1*p1);
			E_prot = TMath::Sqrt(proton_mass*proton_mass+p1*p1);
			y_prot_cms = 0.5*TMath::Log((E_prot+particleA->GetPz())/(E_prot-particleA->GetPz())) - particles.y_cms;
			v1.SetPxPyPzE(particleA->GetPx(),particleA->GetPy(),particleA->GetPz(),E1);

			//if(y_prot_cms > (particles.y_cms - 0.5))		//Quick cross-check
			//	continue;

			y1 = 0.5*TMath::Log((E1+particleA->GetPz())/(E1-particleA->GetPz())) - particles.y_cms;
			angle = TMath::ATan2(particleA->GetPy(), particleA->GetPx());

			particles.analyze(particleA,beam_momentum);

			//debugfile << i << ": " <<  particleA->GetPx() << "\t" << particleA->GetPy() << "\t" << particleA->GetPz() << endl;

			//angle = TMath::Sqrt(TMath::Power(particleA->GetPx(),2) + TMath::Power(particleA->GetPy(),2));
			positive = particleA->isPositive();

			//debugfile << (positive ? "1 " : "-1 ") << angle << endl;
			if(!positive)
				angle3 = mk_angle3(angle);

			if(event->GetNpa() > 1)
			{
				for(j=i+1; j<event->GetNpa(); ++j)
				{
					particleB = event->GetParticle(j);

					//cout << "Particle A [" << i << "]:" << particleA << " particleB[" << j << "]:" << particleB << endl;
					//cout << "Particle A: px=" << particleA->GetPx() << " py=" << particleA->GetPy() << " pz=" << particleA->GetPz() << endl;
					//cout << "Particle B: px=" << particleB->GetPx() << " py=" << particleB->GetPy() << " pz=" << particleB->GetPz() << endl;

			//		if((TMath::Abs(particleB->GetBx()) > 4) || (TMath::Abs(particleB->GetBy()) > 2))
			//			continue;
					pt2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2));

					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));

					//cout << "p1 = " << p1 << " | p2 = " << p2 << endl;

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

					/*
					if(y_prot_cms > (particles.y_cms - 0.5))		//Quick cross-check
						continue;
						*/

					//cout << "E1 = " << E1 << " | E2 = " << E2 << endl;

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
					if((angle_diff = TMath::Abs(angle-angle_j)) > TMath::Pi())
						angle_diff = 2*TMath::Pi()-angle_diff;
					y_diff = TMath::Abs(y1-y2);

					histos.histDyDphiAll->Fill(angle_diff, (y_diff = TMath::Abs(y1-y2)));
					histos.histDetaDphiAll->Fill(angle_diff, (eta_diff = TMath::Abs(eta1-eta2)));
					Fill4Times(histos.histDetaDphiAllReflected, eta_diff, angle_diff);

					//debugfile << "deta=" << eta_diff << endl;
					
					++all_correlations;

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

			n[All]++;
			phi[All] += angle;
			phiSq[All] += angle*angle;

			if(positive)
			{
				n[Pos]++;
				phi[Pos] += angle;
				phiSq[Pos] += angle*angle;
			}
			else
			{
				n[Neg]++;
				phi[Neg] += angle3;
				phiSq[Neg] += angle3*angle3;
			}
		}	

		if((event->GetNpa()!=n[All]))
				cerr << "Event: " << ev << " ParticleTree: " << (event->GetNpa()) << " Extraction: " << (n[All]) << endl;

		//debugfile << "----------" << endl;

		//cout << "\rEvent " << ev;
		if(!(ev%5000))
			cout << "Event " << ev << endl;

		particles.newEvent();
	}

	//event--;

	cout << "All correlations: " << all_correlations << endl;
	cout << "Like-sign correlations: " << correlations << endl;
	cout << "Positive correlations: " << pos_correlations << endl;
	cout << "Negative correlations: " << neg_correlations << endl;
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

