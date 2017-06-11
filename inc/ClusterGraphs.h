#include "TFile.h"
#include "TNtuple.h"
#include "RootWriter.h"


struct AdditionalInfo
{
	Float_t pz_cms1;
	Float_t pz_cms2;
	Float_t eta_cms1;
	Float_t eta_cms2;
	Float_t phi1;
	Float_t phi2;
	Float_t deta;
	Float_t dphi;
};

class ClusterGraphs
{
	TFile *clusters_graphs_file;			//File for storing graphs
	TNtuple	*cluster_positions;				//TNtuple pointer to open tuple with cluster positions
	Float_t  clus_x, clus_y, clus_z;		//float variables to connect their addressess with tuple branches
	Float_t *x1, *y1, *z1,
			*x2, *y2, *z2;

	TString path;
	TString localpath;
	TFile *other_hist_file;

	enum TPCs
	{
		eVTPC1 = 0,
		eGTPC  = 1,
		eVTPC2 = 2,
		eMTPCR = 3,
		eMTPCL = 4
	};

	Int_t nClusters1[5];
	Int_t nClusters2[5];

public:
	
	ClusterGraphs();

	void closeFile();
	void setOtherHistFile(TFile*);
	void addGraph(Int_t, Particle*, Particle*, AdditionalInfo);
	void countClusters(Float_t, Float_t, Float_t, Int_t*);
};
