/****************************************************************/<br>
<br>
ACTG v1.20 : Novel peptide mapping onto gene models<br>
Release: Nov 17, 2020<br>
Hanyang University, Seoul, Korea<br>
<br>
/****************************************************************/<br>

##ver1.20
- Recursive routine was changed to Iterative routine to avoid "STACK OVERFLOW."

- Nucleotide sequence correspoding to a peptide was displayed in FLAT result file.

##ver1.11
- Fixed wrong of frame shift annotation.

##ver1.10
- Added threads mapping.

##ver1.09
- attribute in flat file was modified (easier to understand).

##ver1.08
- mapping peptides using threads.

##ver1.07
- catched wrong file formats in peptide list and VCF.

##ver1.06
- Added frame number in a GFF result file.
- Modified genomic coordinates representation as straight forward way.

##ver1.05
- Fixed that overlapping peptide sequence could not be annotated correctly (contributed by Jounghun Yeom, KIST).

##ver1.04
- Fixed deletion issue.
- Largescale deletion issue has not been solved yet.

##ver1.03
- Fixed infinite loop issue when taken vcf file as an input.

##ver1.02
- Fixed annotating regional information.

##ver1.01
- Added frame shift annotation for regional information.

##ver1.00
- This version supported three methods to map peptides (list below).
  protein database mapping, variant splice graph mapping, six-frame transloation mapping

- Structural variations can be detected by automatical way if a user wants (list below).
  single exon skipping, splice junction, exon extension, frame shift

- Mutations could be detected if a user gives VCF as an input.

