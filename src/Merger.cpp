#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TApplication.h"

#include "ParticleTree.h"

using namespace std;

int merge(TString output_filename, const vector<string> datafiles)
{
	ParticleTree output_tree(output_filename);	

	TString inputfile_path;
	TFile *input_file;
	TTree *input_tree;
	Event *event = new Event();
	Particle *particle;
	Int_t nentries;
	UInt_t global_event = 0;

	for(unsigned int i=0; i<datafiles.size(); i++)
	{
		inputfile_path = datafiles[i];
		cout << "Opening file: " << inputfile_path << endl;
		input_file = new TFile(inputfile_path);
		if(input_file->IsZombie())
		{
			cout << "Input file is zombie, skipping" << endl;
			continue;
		}

		input_tree = (TTree*)input_file->Get("events");
		input_tree->SetBranchAddress("event",&event);
		nentries = input_tree->GetEntries();
		if(nentries == 0)
		{
			cout << "Input file has no events, skipping" << endl;
			continue;
		}

		for(Int_t ev = 0; ev < nentries; ++ev)
		{
			input_tree->GetEntry(ev);
			++global_event;

			output_tree.BeginEvent();	//EID will be incremented internally by ParticleTree
			
			for(UInt_t part = 0; part < event->GetNpa(); part++)
			{
				particle = event->GetParticle(part);

				output_tree.AddParticle(
						particle->GetPDGpid(), particle->GetCharge(),
						particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetMass());

				//cout << "Pid: " << particle->GetPid() << ", Eid: " << particle->GetEid() << ", px= " << particle->GetPx() << ", py= " << particle->GetPy() << endl;
				//PID will be incremented internally by ParticleTree
			}

			if(!(global_event%10000))
				cout << "Global event: " << global_event << ", local event: " << ev << endl;

			output_tree.EndEvent();
		}

		input_file->Close();
	}

	cout << "Closing" << endl;
	output_tree.Close();

	return 0;
}

int main(int argc, char* argv[])
{
	int signal = -1;
	TString output_filename = argv[1];
	const vector<string> datafiles(argv + 2, argv + argc);

	cout << "Processing files: " << endl;
	for(unsigned int i=0;i<datafiles.size();i++)
		cout << datafiles[i] << endl;

	if(argc >= 3)
	{
		TApplication app("app",&argc,argv);
		signal = merge(output_filename, datafiles);
	}
	else
	{
		cout << "Usage: merger <output_file> <file1> [<file2> [<file3> ...]]" << endl;
		return 0;
	}

	return signal;
}
