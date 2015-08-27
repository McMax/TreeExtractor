#ifndef PREFIFI_H
#define PREFIFI_H
#include <string>
#include "TTree.h"
#include "TH2.h"

void Fill4Times(TH2F*, const float, const float);

void mainanalyze(TTree*, const float, const TString);

#endif
