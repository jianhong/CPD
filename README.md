# CPD
CRISPR/Cas9 gRNA designer

## Introduction

CPD is a R/Bioconductor package for CRISPR/Cas9 gRNA design. 
The input is GenomicRanges object. 
Output will be multiple excel sheet of primers and off targets.
This is a parse of CRISPRseek. By input a bed file, users can
search gRNA in upstream or downstream of the coordinates.

## Installation

```
library(devtools)
install_github("jianhong/CPD")
```

## Usage

### Install dependence packages

Now it is support mm10, mm9, hg19, hg38, danRer10. Here is the packages should be installed for each assembly.

|          | genome                          | txdb                                 | orgAnn      
| -------- | ------------------------------- | ------------------------------------ | ------------
| mm10     | BSgenome.Mmusculus.UCSC.mm10    | TxDb.Mmusculus.UCSC.mm10.knownGene   | org.Mm.eg.db
| mm9      | BSgenome.Mmusculus.UCSC.mm9     | TxDb.Mmusculus.UCSC.mm9.knownGene    | org.Mm.eg.db
| hg19     | BSgenome.Hsapiens.UCSC.hg19     | TxDb.Hsapiens.UCSC.hg19.knownGene    | org.Hs.eg.db
| hg38     | BSgenome.Hsapiens.UCSC.hg38     | TxDb.Hsapiens.UCSC.hg38.knownGene    | org.Hs.eg.db
| danRer10 | BSgenome.Drerio.UCSC.danRer10   | TxDb.Drerio.UCSC.danRer10.refGene    | org.Dr.eg.db

For example: for mm10
```
library(BiocInstaller)
biocLite(c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))
```

### Default

By default, the gRNA search is setting as upstream = 2000, downstream = 2000, 
scoring.method = "CFDscore", annotatePaired = FALSE, max.mismatch=1.

 scoring.method = "CFDscore" means not only the mismatch position but also the mismatch type is considered for scoring.
 annotatePaired = FALSE means do NOT search potential paired gRNAs.
 max.mismatch =1 means max mismatch of off target position is 1.
 
 
```
## load the library
library(CPD)
## set the bed file path
gr <- system.file('extdata', 'Lamp3.bed', package = 'CPD')
## search gRNA
## anchor: Anchor point of upstream and downstream. 
##         set to TSS: search gRNA for promoter of given coordinates.
##         set to TES: search gRNA for upstream and downstream of transcript end site.
##         set to ALL: search gRNA for all region including upstream, body and downstream.
##         set to onlyUpstream&Downstream: search gRNA for upstream and downstream of given region (without the given region).
gRNA(gr, anchor="onlyUpstream&Downstream", species="mm10", outputDir="test")
```

### Finding gRNAs with given restriction enzyme cut site(s)

The restriction enzyme cut sites could be set by REpatternFile parameter. It is the file path of given enzymes in fasta file.
The content of restriction enzyme cut site pattern file will be like this:
<pre>
>BamHI
GGATCC
>EcoRI
GAATTC
</pre>

If the file name of the pattern file path is `path/to/enzyme.fa`, set REpatternFile="path/to/enzyme.fa".

```
gRNA(gr, anchor="onlyUpstream&Downstream", species="mm10", outputDir="test", REpatternFile="path/to/enzyme.fa")
```

###  Quick gRNA finding without off-target analysis

Set max.mismatch = 0 to find gRNAs without off-target analysis.

```
gRNA(gr, anchor="onlyUpstream&Downstream", species="mm10", outputDir="test", max.mismatch=0)
```

### More

Get more helps, please refer the vignette and help files of [https://bioconductor.org/packages/release/bioc/html/CRISPRseek.html](CRISPRseek).


