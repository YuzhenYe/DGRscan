#!/bin/bash

echo ">>>Example 1: De novo prediction"
#Example 1: de novo prediction
#input small.seq
#command:
../DGRscan.py -inseq small.seq -summary small-DGRscan.summary

if !((command -v blastall >/dev/null 2>&1) & (command -v blastn >/dev/null 2>&1)); then
	echo "Stop: blastall/blastn needs to be set up for Example 2 & 3"
	exit 1
fi

echo
echo ">>>Example 2: De novo prediction constrained by putative RTs"
#Example 2: de novo prediction constrained by RTs (using -rev_hom option)
#step 1: run blast to search for potential RTs
#make sure that you can use blastall or blastx
if command -v blastall >/dev/null 2>&1; then
	blastall -F F -p blastx -e 1e-10 -m 8 -d ../data/RT-155 -i FQ312004.fna -o FQ312004-vs-RT155.m8
else
	blastx -evalue 1e-3 -outfmt 6 -db ../data/RT-155 -query FQ312004.fna -out FQ312004-vs-RT155.m8
fi
#step 2: run DGRscan
../DGRscan.py -inseq FQ312004.fna -summary FQ312004-DGRscan.summary -rev_hom FQ312004-vs-RT155.m8

echo
echo ">>>Example 3: Homology based prediction"
#Example 3: homology based (using -tr_hom option)
#step 1: run blast to search for potential TRs and VRs (that are similar to reference TRs)
#make sure that you can use blastall
if command -v blastall >/dev/null 2>&1; then
	blastall -F F -p blastn -e 1e-3 -m 8 -d ../data/ref-hv-hmg-TR-nr0.95 -i FQ312004.fna -o FQ312004-vs-TR.m8
else
	#use this command instead if you have blast+ instead of the legacy blast suite
	blastn -evalue 1e-3 -outfmt 6 -db ../data/ref-hv-hmg-TR-nr0.95 -query FQ312004.fna -out FQ312004-vs-TR.m8
fi
#step 2: run DGRscan (using -tr_hom option)
../DGRscan.py -inseq FQ312004.fna -summary FQ312004-DGRscan-hom.summary -tr_hom FQ312004-vs-TR.m8
