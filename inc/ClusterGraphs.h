#include "TFile.h"
#include "TNtuple.h"
#include "RootWriter.h"

class ClusterGraphs
{
	TFile *clusters_graphs_file;			//File for storing graphs
	TNtuple	*cluster_positions;				//TNtuple pointer to open tuple with cluster positions
	Float_t  clus_x, clus_y, clus_z;		//float variables to connect their addressess with tuple branches
	Float_t *x1, *y1, *z1,
			*x2, *y2, *z2;

	TString path;
	TFile *other_hist_file;

public:
	
	ClusterGraphs();

	void closeFile();
	void setOtherHistFile(TFile*);
	void addGraph(Int_t, Particle*, Particle*);
};
