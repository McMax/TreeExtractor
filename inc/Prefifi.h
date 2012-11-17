#ifndef PREFIFI_H
#define PREFIFI_H
#include <string>
#include "TTree.h"

void mainanalyzePb(TTree*, const float, const TString);

void mainanalyze(TTree*, const int, bool, const float, std::string, const TString);

const float  bx[] = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0},
		by[] = {0.05, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.25, 1.5, 2.0};

#endif
