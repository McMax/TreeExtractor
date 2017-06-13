#include "ClusterGraphs.h"
#include "Particle.h"
#include "TFile.h"
#include "TPaveText.h"
#include <iostream>

using namespace std;

Float_t VTPC1_X[] = {-100.12,99.88};
Float_t VTPC1_Y[] = {-48.99,49.01};
Float_t VTPC1_Z[] = {-504.77, -254.77};

Float_t GTPC_X[] = {-40.36, 41.14};
Float_t GTPC_Y[] = {-41.68, 28.32};
Float_t GTPC_Z[] = {-203.28,-173.28};

Float_t VTPC2_X[] = {-99.65, 100.35};
Float_t VTPC2_Y[] = {-49.01, 48.99};
Float_t VTPC2_Z[] = {-125.20, 124.80};

Float_t MTPCR_X[] = {-393.78, -3.78};
Float_t MTPCR_Y[] = {-89.97, 90.03};
Float_t MTPCR_Z[] = {352.61, 742.61};

Float_t MTPCL_X[] = {3.80, 393.80};
Float_t MTPCL_Y[] = {-89.99, 90.01};
Float_t MTPCL_Z[] = {352.29, 742.29};

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

void ClusterGraphs::addGraph(Int_t event_id, UInt_t na61run, UInt_t na61event, Particle* particle1, Particle* particle2, AdditionalInfo ai)
{
	//First particle
	cluster_positions = particle1->GetClustersPositions();
	cluster_positions->SetBranchAddress("x",&clus_x);
	cluster_positions->SetBranchAddress("y",&clus_y);
	cluster_positions->SetBranchAddress("z",&clus_z);

	x1 = new Float_t[cluster_positions->GetEntries()];
	y1 = new Float_t[cluster_positions->GetEntries()];
	z1 = new Float_t[cluster_positions->GetEntries()];

	for(unsigned int tpc=0; tpc<5; tpc++)
	{
		nClusters1[tpc] = 0;
		nClusters2[tpc] = 0;
	}


	for(unsigned int clus=0; clus<cluster_positions->GetEntries(); clus++)
	{
		cluster_positions->GetEntry(clus);
		x1[clus] = clus_x;
		y1[clus] = clus_y;
		z1[clus] = clus_z;
		
		//Counting points
		countClusters(clus_x, clus_y, clus_z, nClusters1);
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
		
		//Counting points
		countClusters(clus_x, clus_y, clus_z, nClusters2);
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
	text->AddText(TString::Format("run=%d, NA61 event: %d", na61run, na61event));
	text->AddText(TString::Format("ch=%d, p=(%.3f,%.3f,%.3f)",(particle1->isPositive()) ? 1 : -1,particle1->GetPx(),particle1->GetPy(),particle1->GetPz()));
	text->AddText(TString::Format("p_{z}^{cms}=%.3f, #eta^{cms}=%.3f, #phi=%.3f",ai.pz_cms1, ai.eta_cms1, ai.phi1));
	text->AddText(TString::Format("nVTPC1=%d, nGTPC=%d, nVTPC2=%d, nMTPCR=%d, nMTPCL=%d", nClusters1[eVTPC1], nClusters1[eGTPC], nClusters1[eVTPC2], nClusters1[eMTPCR], nClusters1[eMTPCL]));
	text->AddText(TString::Format("ppVTPC1=%d, ppGTPC=%d, ppVTPC2=%d, ppMTPC=%d", particle1->GetPPvtpc1(), particle1->GetPPgtpc(), particle1->GetPPvtpc2(), particle1->GetPPmtpc()));
	text->AddText("");
	text->AddText(TString::Format("ch=%d, p=(%.3f,%.3f,%.3f)",(particle2->isPositive()) ? 1 : -1,particle2->GetPx(),particle2->GetPy(),particle2->GetPz()));
	text->AddText(TString::Format("p_{z}^{cms}=%.3f, #eta^{cms}=%.3f, #phi=%.3f",ai.pz_cms2, ai.eta_cms2, ai.phi2));
	text->AddText(TString::Format("nVTPC1=%d, nGTPC=%d, nVTPC2=%d, nMTPCR=%d, nMTPCL=%d", nClusters2[eVTPC1], nClusters2[eGTPC], nClusters2[eVTPC2], nClusters2[eMTPCR], nClusters2[eMTPCL]));
	text->AddText(TString::Format("ppVTPC1=%d, ppGTPC=%d, ppVTPC2=%d, ppMTPC=%d", particle2->GetPPvtpc1(), particle2->GetPPgtpc(), particle2->GetPPvtpc2(), particle2->GetPPmtpc()));
	text->AddText("");
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

void ClusterGraphs::countClusters(Float_t x_pos, Float_t y_pos, Float_t z_pos, Int_t* clusters)
{
	bool pointx_in_vtpc1 = false, pointx_in_gtpc = false,  pointx_in_vtpc2 = false,  pointx_in_mtpcr = false, pointx_in_mtpcl = false;
	bool pointy_in_vtpc1 = false, pointy_in_gtpc = false,  pointy_in_vtpc2 = false,  pointy_in_mtpcr = false, pointy_in_mtpcl = false;
	bool pointz_in_vtpc1 = false, pointz_in_gtpc = false,  pointz_in_vtpc2 = false,  pointz_in_mtpcr = false, pointz_in_mtpcl = false;

	//X axis
	if((x_pos >= VTPC1_X[0]) && (x_pos <= VTPC1_X[1]))
		pointx_in_vtpc1 = true;
	if((x_pos >= GTPC_X[0]) && (x_pos <= GTPC_X[1]))
		pointx_in_gtpc = true;
	if((x_pos >= VTPC2_X[0]) && (x_pos <= VTPC2_X[1]))
		pointx_in_vtpc2 = true;
	if((x_pos >= MTPCR_X[0]) && (x_pos <= MTPCR_X[1]))
		pointx_in_mtpcr = true;
	else if((x_pos >= MTPCL_X[0]) && (x_pos <= MTPCL_X[1]))
		pointx_in_mtpcl = true;
	//Y ayis
	if((y_pos >= VTPC1_Y[0]) && (y_pos <= VTPC1_Y[1]))
		pointy_in_vtpc1 = true;
	if((y_pos >= GTPC_Y[0]) && (y_pos <= GTPC_Y[1]))
		pointy_in_gtpc = true;
	if((y_pos >= VTPC2_Y[0]) && (y_pos <= VTPC2_Y[1]))
		pointy_in_vtpc2 = true;
	if((y_pos >= MTPCR_Y[0]) && (y_pos <= MTPCR_Y[1]))
		pointy_in_mtpcr = true;
	if((y_pos >= MTPCL_Y[0]) && (y_pos <= MTPCL_Y[1]))
		pointy_in_mtpcl = true;
	//Z ayis
	if((z_pos >= VTPC1_Z[0]) && (z_pos <= VTPC1_Z[1]))
		pointz_in_vtpc1 = true;
	else if((z_pos >= GTPC_Z[0]) && (z_pos <= GTPC_Z[1]))
		pointz_in_gtpc = true;
	else if((z_pos >= VTPC2_Z[0]) && (z_pos <= VTPC2_Z[1]))
		pointz_in_vtpc2 = true;
	else
	{
		if((z_pos >= MTPCR_Z[0]) && (z_pos <= MTPCR_Z[1]))
			pointz_in_mtpcr = true;
		if((z_pos >= MTPCL_Z[0]) && (z_pos <= MTPCL_Z[1]))
			pointz_in_mtpcl = true;
	}

	//Adding points
	if(pointx_in_vtpc1 && pointy_in_vtpc1 && pointz_in_vtpc1)
		clusters[eVTPC1]++;
	else if(pointx_in_gtpc && pointy_in_gtpc && pointz_in_gtpc)
		clusters[eGTPC]++;
	else if(pointx_in_vtpc2 && pointy_in_vtpc2 && pointz_in_vtpc2)
		clusters[eVTPC2]++;
	else if(pointx_in_mtpcr && pointy_in_mtpcr && pointz_in_mtpcr)
		clusters[eMTPCR]++;
	else if(pointx_in_mtpcl && pointy_in_mtpcl && pointz_in_mtpcl)
		clusters[eMTPCL]++;
}
