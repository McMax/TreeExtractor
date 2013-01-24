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
	if(argc == 6)
	{
		TFile *root_tree_file = new TFile(argv[1]);
		TTree *particletree = (TTree*)root_tree_file->Get("events");
		TString root_filename = argv[3];
		float energy = atof(argv[2]);
		string fifivsbpar = argv[4];
		TString root_output_filename = argv[5];

		cout << "Reading file" << endl;

		//Biore rozklad ile eventow mialo dany rozmiar
		map<int, int> distribution = getDistro(root_filename);
		mainanalyze(particletree, distribution[0], true, energy, fifivsbpar, root_output_filename);

		root_tree_file->Close();
	}
	else if(argc == 4) //Pb-Pb mode
	{
		TFile *root_tree_file = new TFile(argv[1]);
		TTree *particletree = (TTree*)root_tree_file->Get("events");
		float energy = atof(argv[2]);
		TString root_output_filename = argv[3];

		cout << "Reading file" << endl;

		//Biore rozklad ile eventow mialo dany rozmiar
		mainanalyzePb(particletree, energy, root_output_filename);

		root_tree_file->Close();
	}
	else
	{
		cout << "USAGE: extractor <path_to_particle_tree> <energy> <root_file> <fifivsbpar> <outputfile>" << endl;
		cout << "OR: extractor <path_to_particle_tree> <energy> <outputfile>" << endl;
		return -1;
	}
}
