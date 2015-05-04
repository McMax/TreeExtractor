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
#include "Prefifi.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc == 4)
	{
		TFile *root_tree_file = new TFile(argv[1]);
		TTree *particletree = (TTree*)root_tree_file->Get("events");
		float momentum = atof(argv[2]);
		TString root_output_filename = argv[3];

		cout << "Reading file" << endl;

		//Biore rozklad ile eventow mialo dany rozmiar
		mainanalyze(particletree, momentum, root_output_filename);

		root_tree_file->Close();
	}
	else
	{
		cout << "USAGE: extractor <path_to_particle_tree> <beam_momentum> <outputfile>" << endl;
		return -1;
	}
}
