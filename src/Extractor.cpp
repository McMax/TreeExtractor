#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"

#include "RootWriter.h"
#include "Particle.h"
#include "Event.h"
//#include "Prefifi.h"

using namespace std;

const float  bx[] = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0},
		by[] = {0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.25, 1.5, 2.0};

void mainanalyze(TTree *particletree, const unsigned zeros, bool write_to_root, const float energy, string fifivsbpar)
{
	ofstream prefifi_file("Pre_fifi.txt");
	ofstream debugfile("Debug.txt");
	ofstream prefifi_b_file[12];
	float 	phi_b[3][12],
		phiSq_b[3][12];
	int	n_b[3][12];
	string targettype;

	bool with_fifivsbpar;
	if(fifivsbpar.compare("NONE"))
		with_fifivsbpar = false;
	else if(fifivsbpar.compare("EMPTY"))
	{
		with_fifivsbpar = true;
		targettype = "EMPTY";
	}
	else if(fifivsbpar.compare("FULL"))
	{
		with_fifivsbpar = true;
		targettype = "FULL";
	}
	else
	{
		with_fifivsbpar = false;
		cout << "Target type unknown. Skipping fifivsbpar extraction." << endl;
	}

	if(with_fifivsbpar)
	{
		char prefifi_b_filename[50];
		char ch_targettype[targettype.size()];

		targettype.copy(ch_targettype, targettype.size());

		for(int j=0; j<12; j++)
		{
			sprintf(prefifi_b_filename,"PhiphiVsBpar_%s/B%d.txt",ch_targettype, j);
			prefifi_b_file[j].open(prefifi_b_filename);             
		}
		
	}

	float 	phi[3],
		phiSq[3],
		angle,
		angle3,

		p1, p2,
		pt1, pt2,
		E1, E2,
		theta1, theta2,
		y1, y2,
		eta1, eta2,
		angle_j,
		angle_diff,
		y_diff,
		eta_diff;
	
	bool	positive,
		positive_j;

	int	n[3];
	UInt_t	i,j;

	unsigned correlations = 0, pos_correlations = 0, neg_correlations = 0, all_correlations = 0;;

	Event *event = new Event();
	Particle *particleA, *particleB;
	particletree->SetBranchAddress("event",&event);
	Long64_t treeNentries = particletree->GetEntries();
	Long64_t ev;

	Particles particles;
	Histos histos;
	TFile *root_output_file;

	if(write_to_root)
	{
		histos.init();
		particles.init(&histos);
		particles.newEvent(true);
		root_output_file = new TFile("Extracted_distributions.root","recreate");
	}

	cout << "Writing events" << endl;

	for(ev=0; ev<treeNentries; ++ev)
	{
		particletree->GetEntry(ev);

		phi[Neg] = phi[All] = phi[Pos]= 0.;
		phiSq[Neg] = phiSq[All] = phiSq[Pos] = 0.;
		n[Neg] = n[All] = n[Pos] = 0;

		if(with_fifivsbpar)
		{
			for(j=0; j<12; ++j)
			{
				phi_b[Neg][j] = phi_b[All][j] = phi_b[Pos][j]= 0.;
				phiSq_b[Neg][j] = phiSq_b[All][j] = phiSq_b[Pos][j] = 0.;
				n_b[Neg][j] = n_b[All][j] = n_b[Pos][j] = 0;
			}
		}

		debugfile << ev << endl;

		for(i=0; i<event->GetNpa(); ++i)
		{
			particleA = event->GetParticle(i);
			if(write_to_root)
				particles.analyze(particleA,energy);

			angle = TMath::ATan2(particleA->GetPy(), particleA->GetPx());
			//angle = TMath::Sqrt(TMath::Power(particleA->GetPx(),2) + TMath::Power(particleA->GetPy(),2));
			positive = particleA->isPositive();

			debugfile << (positive ? "1 " : "-1 ") << angle << endl;
			if(!positive)
				angle3 = mk_angle3(angle);

			if(write_to_root && (event->GetNpa() > 1))
			{
				for(j=i+1; j<event->GetNpa(); ++j)
				{
					particleB = event->GetParticle(j);
					pt1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2));
					pt2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2));

					p1 = TMath::Sqrt(TMath::Power(particleA->GetPx(),2)+TMath::Power(particleA->GetPy(),2)+TMath::Power(particleA->GetPz(),2));
					p2 = TMath::Sqrt(TMath::Power(particleB->GetPx(),2)+TMath::Power(particleB->GetPy(),2)+TMath::Power(particleB->GetPz(),2));

					E1 = TMath::Sqrt(0.1396*0.1396+p1*p1);
					E2 = TMath::Sqrt(0.1396*0.1396+p2*p2);

					y1 = lab2cms(0.5*TMath::Log((E1+particleA->GetPz())/(E1-particleA->GetPz())),energy);
					y2 = lab2cms(0.5*TMath::Log((E2+particleB->GetPz())/(E2-particleB->GetPz())),energy);

					theta1 = TMath::ATan(pt1/(particleA->GetPz()));
					theta2 = TMath::ATan(pt2/(particleB->GetPz()));

					if(theta1 > 0)
						eta1 = -TMath::Log(0.5*theta1);
					else
						eta1 = 0.;

					if(theta2 > 0)
						eta2 = -TMath::Log(0.5*theta2);
					else
						eta2 = 0.;

					positive_j = particleB->isPositive();
					angle_j = TMath::ATan2(particleB->GetPy(), particleB->GetPx());
					if((angle_diff = TMath::Abs(angle-angle_j)) > TMath::Pi())
						angle_diff = 2*TMath::Pi()-angle_diff;
					y_diff = TMath::Abs(y1-y2);

					histos.histDyDphiAll->Fill(angle_diff, (y_diff = TMath::Abs(y1-y2)));
					histos.histDetaDphiAll->Fill(angle_diff, (eta_diff = TMath::Abs(eta1-eta2)));
					
					++all_correlations;

					if((positive_j == true) && (positive == true))
					{
						//debugfile << "E: " << ev << " pid1: " << particleA->GetPid() << " pid2: " << particleB->GetPid() << endl;
						++correlations;
						++pos_correlations;

						histos.histDyDphiPos->Fill(angle_diff, y_diff);
						histos.histDetaDphiPos->Fill(angle_diff, eta_diff);
					}
					else if((positive_j == false) && (positive == false))
					{
							++correlations;
							++neg_correlations;
							histos.histDyDphiNeg->Fill(angle_diff, y_diff);
							histos.histDetaDphiNeg->Fill(angle_diff, eta_diff);
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

			

			if(with_fifivsbpar)
			{
				for(j=0; j<12; ++j)
				{
					if((TMath::Abs(particleA->GetBx()) < bx[j]) && (TMath::Abs(particleA->GetBy()) < by[j]))
					{
						n_b[All][j]++;
						phi_b[All][j] += angle;
						phiSq_b[All][j] += angle*angle;

						if(positive)
						{
							n_b[Pos][j]++;
							phi_b[Pos][j] += angle;
							phiSq_b[Pos][j] += angle*angle;
						}
						else
						{
							n_b[Neg][j]++;
							phi_b[Neg][j] += angle3;
							phiSq_b[Neg][j] += angle3*angle3;
						}	
					}
				}
			}
		}	

		debugfile << "----------" << endl;

		prefifi_file << 10000 << "\t\t" << n[All] << "\t" << phi[All] << "\t" << phiSq[All] << "\t\t" << 
				n[Neg] << "\t" << phi[Neg] << "\t" << phiSq[Neg] << "\t\t" << 
				n[Pos] << "\t" << phi[Pos] << "\t" << phiSq[Pos] << endl;

		if(with_fifivsbpar)
		{
			for(j=0; j<12; ++j)
			{
				prefifi_b_file[j] << 10000 << "\t\t" << n_b[All][j] << "\t" << phi_b[All][j] << "\t" << phiSq_b[All][j] << "\t\t" << 
					n_b[Neg][j] << "\t" << phi_b[Neg][j] << "\t" << phiSq_b[Neg][j] << "\t\t" << 
					n_b[Pos][j] << "\t" << phi_b[Pos][j] << "\t" << phiSq_b[Pos][j] << endl;
			}
		}

		cout << "\rEvent " << ev;

		if(write_to_root)
			particles.newEvent();
	}

	//event--;

	cout << endl << "Filling with zeros" << endl;
	unsigned zero_event = 0;
	
	if(with_fifivsbpar)
	{
		while(zero_event+1 <= zeros)
		{
			zero_event++;
			prefifi_file << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;

			for(j=0; j<12; ++j)
				prefifi_b_file[j] << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;
			
			cout << "\rEvent " << zero_event;
		}
		cout << endl << (ev+zero_event) << " lines written to Pre_fifi" << endl;
		prefifi_file.close();
		for(j=0; j<12; ++j)
			prefifi_b_file[j].close();
	}
	else
	{
		while(zero_event+1 <= zeros)
		{
			zero_event++;
			prefifi_file << "10000\t\t0\t0\t0\t\t0\t0\t0\t\t0\t0\t0" << endl;
			cout << "\rEvent " << zero_event;
		}
		cout << endl << (ev+zero_event) << " lines written to Pre_fifi" << endl;
		prefifi_file.close();
	}

	cout << "All correlations: " << all_correlations << endl;
	cout << "Like-sign correlations: " << correlations << endl;
	cout << "Positive correlations: " << pos_correlations << endl;
	cout << "Negative correlations: " << neg_correlations << endl;
	debugfile << "All correlations: " << all_correlations << endl;
	debugfile << "\nLike-sign correlations: " << correlations << endl;
	debugfile << "Positive correlations: " << pos_correlations << endl;
	debugfile << "Negative correlations: " << neg_correlations << endl;

	debugfile.close();

	if(write_to_root)
	{
		//histos.histCharged->AddBinContent(1,zeros);
		//histos.histChargedNeg->AddBinContent(1,zeros);
		//histos.histChargedPos->AddBinContent(1,zeros);
		histos.histCharged->ResetStats();
		histos.histChargedNeg->ResetStats();
		histos.histChargedPos->ResetStats();
		histos.write();
		histos.clear();
		root_output_file->Close();
	}
}

int main(int argc, char** argv)
{
	if(argc < 5)
	{
		cout << "USAGE: extractor <path_to_particle_base> <energy> <root_file> <fifivsbpar=NONE>" << endl;
		return -1;
	}

	TFile *root_tree_file = new TFile(argv[1]);
	TTree *particletree = (TTree*)root_tree_file->Get("events");
	TString root_filename = argv[3];
	float energy = atof(argv[2]);
	string fifivsbpar = argv[4];

	cout << "Reading file" << endl;

	//Biore rozklad ile eventow mialo dany rozmiar
	map<int,int> distribution = getDistro(root_filename);
	mainanalyze(particletree, distribution[0], true, energy, fifivsbpar);

	root_tree_file->Close();
}
