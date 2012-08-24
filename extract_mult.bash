#!/bin/bash
if [ $# -lt 3 ]; then
	echo "Provide three arguments:"
	echo "1) Charge"
	echo "2) Multiplicity"
	echo "3) Dataset (directory name)"
else
	CHARGE=$1
	MULT=$2
	DATASET=$3

	./extractor ../Testy_z_danymi/corrs2mult/$DATASET/mult/$MULT/$CHARGE/ParticleTree_mult$(expr $MULT)_$(expr $CHARGE).root 158 ../Testy_z_danymi/corrs2mult/$DATASET/pp-158-FULL-2009.root NONE && mv Extracted_distributions.root ../Testy_z_danymi/corrs2mult/$DATASET/mult/$MULT/$CHARGE/pp-158-nonmixed-distr.root

	./extractor ../Testy_z_danymi/corrs2mult/$DATASET/mult/$MULT/$CHARGE/MixedParticleTree_mult$(expr $MULT)_$(expr $CHARGE).root 158 ../Testy_z_danymi/corrs2mult/$DATASET/pp-158-FULL-2009.root NONE && mv Extracted_distributions.root ../Testy_z_danymi/corrs2mult/$DATASET/mult/$MULT/$CHARGE/pp-158-mixed-distr.root
fi
