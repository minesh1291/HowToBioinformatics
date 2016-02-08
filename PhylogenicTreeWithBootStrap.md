
###
#  How Could We Generate A Tree With Bootstrap Value
###

# Requirements

#    (1) Muscle Tool
#            http://www.drive5.com/muscle/, http://www.ncbi.nlm.nih.gov/pubmed/15034147
#    (2) embassy-phylip
#            http://emboss.sourceforge.net/apps/cvs/embassy/phylipnew/
#    (3) FastTreeMP64
#            http://www.microbesonline.org/fasttree/#Install, http://mbe.oxfordjournals.org/content/26/7/1641.full
#    (4) formatTreeFile.pl
#            https://github.com/minesh1291/Phylogeny-Utilities/blob/master/formatTreeFile.pl
#    (5) ETE Python
#            http://etetoolkit.org/
#    (6) DrawTree.py
#            https://github.com/minesh1291/Phylogeny-Utilities/blob/master/DrawTree.py

#1 Fasta File with Multiple Sequences
# eg. Seqs.fa
#2 Multiple Sequence Alignment with Muscle 
muscle -in Seqs.fa -out Seqs.aln
#3 Bootstraping Multiple Sequence Alignment
fseqboot -sequence Seqs.aln -outfile Seqs.boot
#4 Calculating Phylogenic Tree (with multiple threading)
FastTreeMP64 -n 100 <  Seqs.aln.boot > Seqs.trees
#5 Calculating Consensus Tree
fconsense -intreefile Seqs.trees -outfile Seqs.treefile  -method ml -trout -treeprint Y -outtreefile Seqs.tree
#6 Formatting Tree File
perl formatTreeFile.pl Seqs.tree Seqs.BS.tree
#7 Draw Tree PNG image
python DrawTree.py Seqs.BS.tree
