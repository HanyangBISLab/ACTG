# ACTG

##ver1.09
- attribute in flat file was modified (easier to understand).

##ver1.08
- mapping peptides using threads.

##ver1.07
- catch wrong file formats in peptide list and VCF.

##ver1.06
- Add frame number in a GFF result file.
- Modify genomic coordinates representation as straight forward way.

##ver1.05
- Fix that overlapping peptide sequence could not be annotated correctly (contributed by KIST Jounghun Yeom).

##ver1.04
- Fix deletion issue.
- Largescale deletion issue has not been solved yet.

##ver1.03
- Fix infinite loop issue when taken vcf file as an input.

##ver1.02
- Fix annotating regional information.

##ver1.01
- Added frame shift annotation for regional information.

##ver1.00
- This version supports three methods to map peptides (listed below).
  protein database mapping, variant splice graph mapping, six-frame transloation mapping

- Structural variations can be detected by automatical way if a user wants (listed below).
  single exon skipping, splice junction, exon extension, frame shift

- Mutations can be detected if a user gives VCF as an input.

- Proposing a regional information which is the best fix on genome annotation (listed below).
  CDS, utr, pseudo, intron, and unknown
