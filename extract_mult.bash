#!/bin/bash
if [ $# -lt 2 ]; then
	echo "Provide two arguments:"
	echo "1) Charge"
	echo "2) Multiplicity"
else
	CHARGE=$1
	MULT=$2

	./extractor ../Testy_z_danymi/corrs2mult/Dane/mult/$MULT/$CHARGE/ParticleTree_mult$(expr $MULT)_$(expr $CHARGE).root 158 ../Testy_z_danymi/corrs2mult/Dane/pp-158-FULL-2009.root NONE && mv Extracted_distributions.root ../Testy_z_danymi/corrs2mult/Dane/mult/$MULT/$CHARGE/pp-158-nonmixed-distr.root

	./extractor ../Testy_z_danymi/corrs2mult/Dane/mult/$MULT/$CHARGE/MixedParticleTree_mult$(expr $MULT)_$(expr $CHARGE).root 158 ../Testy_z_danymi/corrs2mult/Dane/pp-158-FULL-2009.root NONE && mv Extracted_distributions.root ../Testy_z_danymi/corrs2mult/Dane/mult/$MULT/$CHARGE/pp-158-mixed-distr.root
fi
