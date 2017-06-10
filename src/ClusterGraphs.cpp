#include "ClusterGraphs.h"
#include "Particle.h"
#include "TFile.h"
#include "TPaveText.h"


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

void ClusterGraphs::addGraph(Int_t event_id, Particle* particle1, Particle* particle2, AdditionalInfo ai)
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

	TPaveText *text = new TPaveText(-700,150,-200,450,"NB");
	text->SetTextColor(kRed);
	text->AddText(TString::Format("ch=%d, p=(%.3f,%.3f,%.3f)",(particle1->isPositive()) ? 1 : -1,particle1->GetPx(),particle1->GetPy(),particle1->GetPz()));
	text->AddText(TString::Format("p_{z}^{cms}=%.3f, #eta^{cms}=%.3f, #phi=%.3f",ai.pz_cms1, ai.eta_cms1, ai.phi1));
	text->AddText("");
	text->SetTextColor(kBlue);
	text->AddText(TString::Format("ch=%d, p=(%.3f,%.3f,%.3f)",(particle2->isPositive()) ? 1 : -1,particle2->GetPx(),particle2->GetPy(),particle2->GetPz()));
	text->AddText(TString::Format("p_{z}^{cms}=%.3f, #eta^{cms}=%.3f, #phi=%.3f",ai.pz_cms2, ai.eta_cms2, ai.phi2));
	text->AddText("");
	text->SetTextColor(kBlack);
	text->AddText(TString::Format("#Delta#eta=%.4f, #Delta#phi=%.4f",ai.deta,ai.dphi));

	//Saving to file
	path.Form("e%d_p%d_p%d", event_id, particle1->GetPid(), particle2->GetPid());
	clusters_graphs_file->mkdir(path);
	clusters_graphs_file->cd(path);
	cluster_graph1->Write();
	cluster_graph2->Write();
	text->Write();

	other_hist_file->cd();
}
