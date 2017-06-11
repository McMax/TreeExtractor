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
	//First particle
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

	TGraph *cluster_graph_xz1 = new TGraph(cluster_positions->GetEntries(), z1, x1); //Z is drawn as horizontal axis and X as vertical
	TGraph *cluster_graph_yz1 = new TGraph(cluster_positions->GetEntries(), z1, y1); //Z is drawn as horizontal axis and Y as vertical
	TGraph *cluster_graph_xy1 = new TGraph(cluster_positions->GetEntries(), x1, y1);

	//Second particle
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

	TGraph *cluster_graph_xz2 = new TGraph(cluster_positions->GetEntries(), z2, x2);
	TGraph *cluster_graph_yz2 = new TGraph(cluster_positions->GetEntries(), z2, y2);
	TGraph *cluster_graph_xy2 = new TGraph(cluster_positions->GetEntries(), x2, y2);

	//-----------------------------------------------------
	//Drawing clusters

	//XZ plane
	cluster_graph_xz1->SetName(TString::Format("%d", particle1->GetPid()));
	cluster_graph_xz2->SetName(TString::Format("%d", particle2->GetPid()));

	//YZ plane
	cluster_graph_yz1->SetName(TString::Format("%d", particle1->GetPid()));
	cluster_graph_yz2->SetName(TString::Format("%d", particle2->GetPid()));

	//XY plane
	cluster_graph_xy1->SetName(TString::Format("%d", particle1->GetPid()));
	cluster_graph_xy2->SetName(TString::Format("%d", particle2->GetPid()));


	TPaveText *text = new TPaveText();
	text->SetName("info");
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

	localpath.Form("%s/xz", path.Data());
	clusters_graphs_file->mkdir(localpath);
	clusters_graphs_file->cd(localpath);
	cluster_graph_xz1->Write();
	cluster_graph_xz2->Write();

	localpath.Form("%s/yz", path.Data());
	clusters_graphs_file->mkdir(localpath);
	clusters_graphs_file->cd(localpath);
	cluster_graph_yz1->Write();
	cluster_graph_yz2->Write();

	localpath.Form("%s/xy", path.Data());
	clusters_graphs_file->mkdir(localpath);
	clusters_graphs_file->cd(localpath);
	cluster_graph_xy1->Write();
	cluster_graph_xy2->Write();

	clusters_graphs_file->cd(path);
	text->Write();

	other_hist_file->cd();
}
