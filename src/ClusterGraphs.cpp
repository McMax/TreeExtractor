#include "ClusterGraphs.h"
#include "Particle.h"
#include "TFile.h"


ClusterGraphs::ClusterGraphs()
{
	clusters_graphs_file = new TFile("Graphs_with_clusters.root","recreate");
}

void ClusterGraphs::closeFile()
{
	clusters_graphs_file->Close();
}

void ClusterGraphs::setOtherHistFile(TFile* hist_file)
{
	other_hist_file = hist_file;
}

void ClusterGraphs::addGraph(Int_t event_id, Particle* particle1, Particle* particle2)
{
	//First graph
	cluster_positions = particle1->GetClustersPositions();
	cluster_positions->SetBranchAddress("x",&clus_x);
	cluster_positions->SetBranchAddress("y",&clus_y);
	cluster_positions->SetBranchAddress("z",&clus_z);

	x1 = new Float_t[cluster_positions->GetEntries()];
	y1 = new Float_t[cluster_positions->GetEntries()];
	z1 = new Float_t[cluster_positions->GetEntries()];

	for(unsigned int clus=0; clus<cluster_positions->GetEntries(); clus++)
	{
		cluster_positions->GetEntry(clus);
		x1[clus] = clus_x;
		y1[clus] = clus_y;
		z1[clus] = clus_z;
	}

	TGraph *cluster_graph1 = new TGraph(cluster_positions->GetEntries(), z1, x1); //Z is drawn as horizontal axis and X as vertical
	cluster_graph1->SetName(TString::Format("%d", particle1->GetPid()));

	//Second graph
	cluster_positions = particle2->GetClustersPositions();
	cluster_positions->SetBranchAddress("x",&clus_x);
	cluster_positions->SetBranchAddress("y",&clus_y);
	cluster_positions->SetBranchAddress("z",&clus_z);

	x2 = new Float_t[cluster_positions->GetEntries()];
	y2 = new Float_t[cluster_positions->GetEntries()];
	z2 = new Float_t[cluster_positions->GetEntries()];

	for(unsigned int clus=0; clus<cluster_positions->GetEntries(); clus++)
	{
		cluster_positions->GetEntry(clus);
		x2[clus] = clus_x;
		y2[clus] = clus_y;
		z2[clus] = clus_z;
	}

	TGraph *cluster_graph2 = new TGraph(cluster_positions->GetEntries(), z2, x2); //Z is drawn as horizontal axis and X as vertical
	cluster_graph2->SetName(TString::Format("%d", particle2->GetPid()));

	//Saving to file
	path.Form("e%d_p%d_p%d", event_id, particle1->GetPid(), particle2->GetPid());
	clusters_graphs_file->mkdir(path);
	clusters_graphs_file->cd(path);
	cluster_graph1->Write();
	cluster_graph2->Write();

	other_hist_file->cd();
}
