#!/bin/bash
FIFI_PATH=$ANALYSIS_DIR/Dane/Skrypty

for i in 20 31 40 80 158
do
	echo "./extractor files/particlebasefiles/VGCALOR/pp-$i-nonmixed-base $i files/rootfiles/VGCALOR/pp-$i-MC-2009.root NONE"
	#echo "./extractor files/particlebasefiles/pp-$i-nonmixed-base $i files/rootfiles/pp-$i-FULL-2009.root NONE"
	sleep 1

#BEFORE MIX

	./extractor files/particlebasefiles/VGCALOR/pp-$i-nonmixed-base $i files/rootfiles/VGCALOR/pp-$i-MC-2009.root NONE
	#./extractor files/particlebasefiles/pp-$i-nonmixed-base $i files/rootfiles/pp-$i-FULL-2009.root NONE
	mv Extracted_distributions.root Extracted_distributions/pp-$i-nonmixed-distr.root

	echo "Shuffling"
	shuf Pre_fifi.txt > Pre_fifi.shuf && $FIFI_PATH/./fifi < Pre_fifi.shuf > Fifi_$i-nonmixed.txt && rm Pre_fifi.shuf
	tail -n9 Fifi_$i-nonmixed.txt

#AFTER MIX

	./extractor files/particlebasefiles/VGCALOR/pp-$i-mixed-base $i files/rootfiles/VGCALOR/pp-$i-MC-2009.root NONE
	#./extractor files/particlebasefiles/pp-$i-mixed-base $i files/rootfiles/pp-$i-FULL-2009.root NONE
	mv Extracted_distributions.root Extracted_distributions/pp-$i-mixed-distr.root

	echo "Shuffling"
	shuf Pre_fifi.txt > Pre_fifi.shuf && $FIFI_PATH/./fifi < Pre_fifi.shuf > Fifi_$i-mixed.txt && rm Pre_fifi.shuf
	tail -n9 Fifi_$i-mixed.txt
done
