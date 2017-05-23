#ifndef PREFIFI_H
#define PREFIFI_H
#include <string>
#include "TTree.h"
#include "TH2.h"

#include "Particle.h"

void Fill4Times(TH2F*, const float, const float);

void mainanalyze(TTree*, const TString, const float, const TString);

Double_t calculate_distance(Particle*,Particle*);

#endif
