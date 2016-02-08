
#  How Can We Generate A Tree With Bootstrap Value

### Requirements

1. [Muscle Tool](http://www.drive5.com/muscle/), [reference](http://www.ncbi.nlm.nih.gov/pubmed/15034147)
2. [embassy-phylip](http://emboss.sourceforge.net/apps/cvs/embassy/phylipnew/)
3. [FastTreeMP64](http://www.microbesonline.org/fasttree/#Install), [reference](http://mbe.oxfordjournals.org/content/26/7/1641.full)
4. [formatTreeFile.pl](https://github.com/minesh1291/Phylogeny-Utilities/blob/master/formatTreeFile.pl)
5. [ETE Python](http://etetoolkit.org/)
6. [DrawTree.py](https://github.com/minesh1291/Phylogeny-Utilities/blob/master/DrawTree.py)

### Steps

1. Fasta File with Multiple Sequences. eg. Seqs.fasta

2. Multiple Sequence Alignment with Muscle 
  ```bash
  muscle -in Seqs.fa -out Seqs.aln
  ```
3. Bootstraping Multiple Sequence Alignment
  ```bash
  fseqboot -sequence Seqs.aln -outfile Seqs.boot
  ```
4. Calculating Phylogenic Tree (with multiple threading)
  ```bash
  FastTreeMP64 -n 100 <  Seqs.aln.boot > Seqs.trees
  ```
5. Calculating Consensus Tree
  ```bash
  fconsense -intreefile Seqs.trees -outfile Seqs.treefile  -method ml -trout -treeprint Y -outtreefile Seqs.tree
  ```
6. Formatting Tree File
  ```bash
  perl formatTreeFile.pl Seqs.tree Seqs.BS.tree
  ```
7. Draw Tree PNG image
  ```bash
  python DrawTree.py Seqs.BS.tree
  ```
