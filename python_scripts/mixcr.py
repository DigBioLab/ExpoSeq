#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 18:36:36 2021

@author: chvis
"""

# mixcr batch alignment on the unpacked fastq files - with help from habib
for i in *2.fastq
do
input=$i
output=${i##NG-28545_Library_}.vdjca
report=${i##NG-28545_Library_}_AlignmentReport.txt
mixcr align -p kaligner2 --species hsa --threads 7 --report ./$report $input ./$output
done

# Chris test with different parameters - only using V library due to the iontas library being contrstructed from these
for i in *2.fastq
do
input=$i
output=${i##NG-28545_}.vdjca
report=${i##NG-28545_}_AlignmentReport.txt
mixcr align -p default --species hsa --report ./$report $input ./$output
done

# mixcr assemble command -
# This function does not cluster sequences with identical CDR3s if they are unique elsewhere
for i in *2.fastq.vdjca
do
input=$i
output=${i%%.fastq.vdjca}_UniqueSeq.clns
report=${i%%.fastq.vdjca}_AssemRep_UniqueSeq.txt
mixcr assemble -OseparateByC=true -OseparateByV=true -OseparateByJ=true --report ./$report $input ./$output
done

# mixcr assemble command - Same as above with different clustering/assembling
# This function clusters sequences with identical CDR3s even though they have unique regions elsewhere
for i in *2.fastq.vdjca
do
input=$i
output=${i%%.fastq.vdjca}_UniqueCDR3.clns
report=${i%%.fastq.vdjca}_AssemRep_UniqueCDR3.txt
mixcr assemble --report ./$report $input ./$output
done
#with first assemble i had 5% of the seq lost to failed mapping and 2% in PCR error corection

# Export 1 - exportclones with unique CDR3s only
for i in *2_UniqueCDR3.clns
do
input=$i
output=${i%%.clns}_Exp_UniqueCDR3.txt
mixcr exportClones -cloneId -count -fraction -lengthOf CDR3 -nFeature CDR3 -aaFeature CDR3 -avrgFeatureQuality CDR3 ./$input ./$output
done

# Export 2 - exportclones with unique seqs
for i in *2_UniqueSeq.clns
do
input=$i
output=${i%%.clns}_Exp_UniqueSeq.txt
mixcr exportClones -cloneId -count -fraction -lengthOf CDR3 -nFeature CDR3 -aaFeature CDR3 -avrgFeatureQuality CDR3 ./$input ./$output
done

# I assume the next parts is to start visualzing some of the data, but i am unsure exactly how
# Tried this: brew install deeptools - didn't work