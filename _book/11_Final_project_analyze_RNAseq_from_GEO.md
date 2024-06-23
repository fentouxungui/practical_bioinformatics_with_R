
--- 
always_allow_html: true
editor_options: 
  chunk_output_type: console
---


# Final Project: Analyzing RNAseq Data from GEO

## Final Project Overview

In this final project, we will work with RNAseq data obtained from the `GEO` database, specifically the dataset with ITPR3 and RELB knockout in the SW480 cell line under varying oxygen conditions. Our primary objectives are to clean and prepare the metadata, identify differentially expressed genes, explore data distribution, and visualize gene expression patterns.

The steps we'll take:

1. To start, we will download the RNAseq count table and the associated metadata. With the help of the `dplyr` package in R, we will clean and organize the metadata. Our focus will be on selecting and subsetting the samples that are most relevant to our analysis - specifically, two wild-type samples under normoxia and two under hypoxia. For the initial phase of the project, we will disregard knockout samples to simplify our analysis.

2. Next, we will delve into exploring the RNAseq count matrix. Without considering knockout samples, we will calculate essential summary statistics to gain insights into gene expression levels' variability and distribution between the two conditions - normoxia and hypoxia.

3. After understanding the dataset's characteristics, we will proceed to identify differentially expressed genes. This step is crucial for uncovering genes that exhibit significant expression changes in response to varying oxygen levels. We will utilize well-established tool, `DESeq2`, to perform differential gene expression analysis.

4. Following the identification of differentially expressed genes, we will continue exploring the count matrix by calculating summary statistics specifically for these genes. This allows us to gain a deeper understanding of how their expression patterns differ between normoxia and hypoxia conditions.

5. To visualize data p-value distribution more effectively, we will create histograms. These histograms will offer a visual representation of how p values are distributed.

6. Additionally, we will use boxplots to compare gene expression levels between normoxia and hypoxia. Boxplots provide a concise summary of the data's central tendency, spread, and potential outliers, aiding in the identification of expression differences.

7. Principal Component Analysis (PCA) will be employed to obtain insights into the overall structure of the data and any potential clustering patterns. This dimensionality reduction technique will help us visualize how samples group based on their gene expression profiles.

8. Finally, we will create a heatmap using R. This heatmap will visualize the expression patterns of the identified differentially expressed genes across samples, providing a comprehensive view of how these genes respond to changes in oxygen levels in the SW480 cell line.

In summary, this project will take us through the full process of RNAseq data analysis, focusing on the hypoxia vs normoxia comparison in the SW480 cell line. We'll clean the data, pinpoint differentially expressed genes, explore their distributions, and visualize gene expression patterns.

Are you ready? Let's go!

## How to pre-process RNAseq data

This is a bonus section on how to pre-process RNAseq data. In this course, we mainly focus on how to analyze RNAseq data for downstream analysis. We will start with a count matrix (next lesson) downloaded from `GEO`.

However, in real-world data analysis, sequencing data comes as a `FASTQ` file. FASTQ files are just normal text files with 4 lines for each read. Go to the link to understand the format.

Watch this video:


```{=html}
<div class="vembedr">
<div>
<iframe src="https://www.youtube.com/embed/Gc-Hzvt7KVQ" width="533" height="300" frameborder="0" allowfullscreen="" data-external="1"></iframe>
</div>
</div>
```


### RNAseq pre-processing steps

1. The first step is to do Quality control of the [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

2. Trim adaptors and low-quality bases using tools such as [trimmomatic](https://github.com/timflutre/trimmomatic) or [fastp](https://github.com/OpenGene/fastp). Trimming of the reads is optional.

3. Align the reads to transcriptome using [`STAR`](https://github.com/alexdobin/STAR). The single-cell RNAseq version is called [STAR-solo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) from the same lab.

4. Quantify the number of reads fall into each gene using [FeatureCounts](https://subread.sourceforge.net/).

I have written a Snakemake pipeline to pre-process RNAseq fastq file to get a count matrix at https://github.com/crazyhottommy/pyflow-RNAseq

The bash script to align the fastq files to transcriptome using `STAR`:


```bash
STAR --runMode alignReads \
		--runThreadN 5 \
		--genomeDir /path/to/the/STAR/index \
		--genomeLoad NoSharedMemory \
		--readFilesIn mysample_R1.fastq.gz mysample_R2.fastq.gz \
		--readFilesCommand zcat \
		--twopassMode Basic \
		--runRNGseed 777 \
		--outFilterType Normal \
		--outFilterMultimapNmax 20 \
		--outFilterMismatchNmax 10 \
		--outFilterMultimapScoreRange 1 \
		--outFilterMatchNminOverLread 0.33 \
		--outFilterScoreMinOverLread 0.33 \
		--outReadsUnmapped None \
		--alignIntronMin 20 \
		--alignIntronMax 500000 \
		--alignMatesGapMax 1000000 \
		--alignSJoverhangMin 8 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--sjdbScore 2 \
		--alignSJDBoverhangMin 1 \
		--sjdbOverhang 100 \
		--chimSegmentMin 20 \
		--chimJunctionOverhangMin 20 \
		--chimSegmentReadGapMax 3 \
		--quantMode GeneCounts \
		--outMultimapperOrder Random \
		--outSAMstrandField intronMotif \
		--outSAMattributes All \
		--outSAMunmapped Within KeepPairs \
		--outSAMtype BAM Unsorted \
		--limitBAMsortRAM 30000000000 \
		--outSAMmode Full \
		--outSAMheaderHD @HD VN:1.4 \
		--outFileNamePrefix mysample
```


Then quantifying using `FeatureCounts`:


```bash
featureCounts -T 5 -p -t exon -g gene_id -a gene.gtf -o mysample_featureCount.txt mysampleAligned.out.bam

```

`mysample_featureCount.txt` will be a count table for one sample.

Alternative Alignment-free RNAseq quantification tools such as [salmon](https://combine-lab.github.io/salmon/getting_started/) and [kallisto](https://pachterlab.github.io/kallisto/) are also very popular. I recommend you to read the tutorial of `STAR`, `FeatureCounts`, `salmon` and `kallisto` to learn how to use those command line tools.

### How to use salmon to preprocess GEO fastq to counts

Please refer to this [blog post](https://divingintogeneticsandgenomics.com/post/how-to-preprocess-geo-bulk-rnaseq-data-with-salmon/) and this youtube video if you want to learn more:


```{=html}
<div class="vembedr">
<div>
<iframe src="https://www.youtube.com/embed/_Q8fYokTCTs" width="533" height="300" frameborder="0" allowfullscreen="" data-external="1"></iframe>
</div>
</div>
```

## Download and subset Count Matrix

In this lesson, you will learn how to explore a count matrix in R, a common task in data analysis. We'll cover downloading the data, reading it into R, examining the data's dimensions and column names, subsetting the data, converting it into a matrix, and adding row names.

### Downloading the Count Matrix

To begin, we need to download the count matrix, which is a tab-separated values (TSV) file containing gene expression data. You can obtain it from the following FTP address: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197576/suppl/GSE197576_raw_gene_counts_matrix.tsv.gz

You can use the wget command in Unix to download the processed count matrix file as follows:


```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197576/suppl/GSE197576_raw_gene_counts_matrix.tsv.gz
```

>If you don't have access to a Unix-like command-line environment, you can download the file manually through your web browser. Simply open your web browser and go to the following URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197576/suppl/GSE197576_raw_gene_counts_matrix.tsv.gz

### Reading the Count Matrix in R
Now that you have the data, let's read it into R using the readr package. The first step is to load the required libraries and read the TSV file:


```r
library(dplyr)
library(readr)

raw_counts <- read_tsv("~/Downloads/GSE197576_raw_gene_counts_matrix.tsv.gz")
```

### Examining the Data
It's crucial to understand the data structure. We'll start by examining the dimensions of the data and the column names:


```r
# Check the dimensions of the data
dim(raw_counts)  # This will show the number of rows and columns
```

```
## [1] 43809    13
```


```r
# List the column names
colnames(raw_counts)  # This will display the names of all columns
```

```
##  [1] "gene"                 "01_SW_sgCTRL_Norm"    "02_SW_sgCTRL_Norm"   
##  [4] "03_SW_sgITPR3_1_Norm" "04_SW_sgITPR3_1_Norm" "07_SW_sgRELB_3_Norm" 
##  [7] "08_SW_sgRELB_3_Norm"  "11_SW_sgCTRL_Hyp"     "12_SW_sgCTRL_Hyp"    
## [10] "13_SW_sgITPR3_1_Hyp"  "14_SW_sgITPR3_1_Hyp"  "17_SW_sgRELB_3_Hyp"  
## [13] "18_SW_sgRELB_3_Hyp"
```

The first six rows:


```r
head(raw_counts)
```

```
## # A tibble: 6 x 13
##   gene        `01_SW_sgCTRL_Norm` `02_SW_sgCTRL_Norm` `03_SW_sgITPR3_1_Norm`
##   <chr>                     <dbl>               <dbl>                  <dbl>
## 1 DDX11L1                       0                   0                      0
## 2 WASH7P                       18                  11                     28
## 3 MIR6859-1                     5                   1                      6
## 4 MIR1302-2HG                   0                   0                      0
## 5 MIR1302-2                     0                   0                      0
## 6 FAM138A                       0                   0                      0
## # i 9 more variables: `04_SW_sgITPR3_1_Norm` <dbl>,
## #   `07_SW_sgRELB_3_Norm` <dbl>, `08_SW_sgRELB_3_Norm` <dbl>,
## #   `11_SW_sgCTRL_Hyp` <dbl>, `12_SW_sgCTRL_Hyp` <dbl>,
## #   `13_SW_sgITPR3_1_Hyp` <dbl>, `14_SW_sgITPR3_1_Hyp` <dbl>,
## #   `17_SW_sgRELB_3_Hyp` <dbl>, `18_SW_sgRELB_3_Hyp` <dbl>
```

>Notice that the first column name is the gene name, and the other 12 columns are the sample names.

### Subsetting the Data

Next, let's narrow down our data to the specific samples we need for comparison. To do this, we'll create a logical vector by matching column names that contain "sgCTRL" or "gene":


```r
columns_to_select <- colnames(raw_counts) %>%
  stringr::str_detect("sgCTRL|gene")

columns_to_select
```

```
##  [1]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
## [13] FALSE
```

The `stringr::str_detect` function searches for patterns in the column names. The resulting columns_to_select variable is a logical vector that helps us select the relevant columns.

Now, let's use this logical vector to subset the data frame:


```r
counts_sub <- raw_counts[, columns_to_select]

head(counts_sub)
```

```
## # A tibble: 6 x 5
##   gene        `01_SW_sgCTRL_Norm` `02_SW_sgCTRL_Norm` `11_SW_sgCTRL_Hyp`
##   <chr>                     <dbl>               <dbl>              <dbl>
## 1 DDX11L1                       0                   0                  0
## 2 WASH7P                       18                  11                 23
## 3 MIR6859-1                     5                   1                  8
## 4 MIR1302-2HG                   0                   0                  1
## 5 MIR1302-2                     0                   0                  0
## 6 FAM138A                       0                   0                  0
## # i 1 more variable: `12_SW_sgCTRL_Hyp` <dbl>
```

### Converting to a Matrix
To perform various analyses, it's often more convenient to work with a matrix. We'll remove the first column (gene names) and convert the data frame to a matrix:


```r
#subset the dataframe by removing the first column using negative index
# and then use as.matrix to convert it to a matrix
raw_counts_mat<- counts_sub[, -1] %>% 
              as.matrix()

head(raw_counts_mat)
```

```
##      01_SW_sgCTRL_Norm 02_SW_sgCTRL_Norm 11_SW_sgCTRL_Hyp 12_SW_sgCTRL_Hyp
## [1,]                 0                 0                0                0
## [2,]                18                11               23               45
## [3,]                 5                 1                8               12
## [4,]                 0                 0                1                0
## [5,]                 0                 0                0                0
## [6,]                 0                 0                0                0
```

Here, we use the `%>%` (pipe) operator to perform multiple operations sequentially. The `as.matrix()` function converts the data frame to a matrix.

###  Adding Row Names

The matrix lacks row names, which can be crucial for identifying genes. We can add the gene names as row names:

```r
rownames(raw_counts_mat) <- raw_counts$gene

head(raw_counts_mat)
```

```
##             01_SW_sgCTRL_Norm 02_SW_sgCTRL_Norm 11_SW_sgCTRL_Hyp
## DDX11L1                     0                 0                0
## WASH7P                     18                11               23
## MIR6859-1                   5                 1                8
## MIR1302-2HG                 0                 0                1
## MIR1302-2                   0                 0                0
## FAM138A                     0                 0                0
##             12_SW_sgCTRL_Hyp
## DDX11L1                    0
## WASH7P                    45
## MIR6859-1                 12
## MIR1302-2HG                0
## MIR1302-2                  0
## FAM138A                    0
```

Now, our matrix has gene names associated with each row.

## Calculate the total exon length per gene

We want to find the differentially expressed genes between hypoxia and normoxia. The ideal workflow is to use the `DESeq2` R package, which models the count data with the negative binomial distribution. For now, let’s use a t-test to compare and we will use DESeq2 later.

However, because different samples have different sequencing depths, and different genes have different lengths, we must first normalize the counts to transcript per million (TPM).

### Why Normalize Gene Counts to TPM?

When working with gene expression data, it's essential to account for variations in sequencing depths (the number of reads obtained for each sample) and gene lengths (some genes are longer than others). Normalization helps ensure that our gene expression values are comparable across different samples and genes.

Transcript Per Million (TPM) is a commonly used normalization method in RNA-seq analysis. It scales the gene counts to a common unit (per million) based on gene length and sequencing depth.

**The general process**:

1. For each gene, we divide its raw count values by its total exon length. This step essentially converts the raw counts into counts per unit exon length.

2. We sum up the values for each column (sample). This calculation gives us the total counts for each sample.

3. For each gene in each sample, we divide the value obtained in Step 1 by the total count for that sample calculated in Step 2. This step scales the counts relative to the sequencing depth of each sample.

4. To bring the values to a common scale and make them more interpretable, we multiply the result from Step 3 by 1,000,000 (1e6). This step converts the values to TPM, where the final unit is "Transcripts Per Million."

The normalization process ensures that the gene expression values are now comparable across different samples, making it easier to identify genes that are differentially expressed between conditions (e.g., hypoxia and normoxia).

### Creating a function for normalization

Let’s write a function. Before that, we need to know the gene length of all the genes in the data frame.

We will use the Bioconductor package `TxDb.Hsapiens.UCSC.hg19.knownGene` to get the gene length, or more exactly, the total exon lengths for each gene (most of the RNAseq reads are from the exons).


```r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

TxDb.Hsapiens.UCSC.hg19.knownGene
```

```
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: UCSC
## # Genome: hg19
## # Organism: Homo sapiens
## # Taxonomy ID: 9606
## # UCSC Table: knownGene
## # Resource URL: http://genome.ucsc.edu/
## # Type of Gene ID: Entrez Gene ID
## # Full dataset: yes
## # miRBase build ID: GRCh37
## # transcript_nrow: 82960
## # exon_nrow: 289969
## # cds_nrow: 237533
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2015-10-07 18:11:28 +0000 (Wed, 07 Oct 2015)
## # GenomicFeatures version at creation time: 1.21.30
## # RSQLite version at creation time: 1.0.0
## # DBSCHEMAVERSION: 1.1
```

It is a `TxDb` object. We can use functions such as genes and exons to get the genes or exons.


```r
# make it a shorter name
txdb<- TxDb.Hsapiens.UCSC.hg19.knownGene

genes(txdb)
```

```
## GRanges object with 23056 ranges and 1 metadata column:
##         seqnames              ranges strand |     gene_id
##            <Rle>           <IRanges>  <Rle> | <character>
##       1    chr19   58858172-58874214      - |           1
##      10     chr8   18248755-18258723      + |          10
##     100    chr20   43248163-43280376      - |         100
##    1000    chr18   25530930-25757445      - |        1000
##   10000     chr1 243651535-244006886      - |       10000
##     ...      ...                 ...    ... .         ...
##    9991     chr9 114979995-115095944      - |        9991
##    9992    chr21   35736323-35743440      + |        9992
##    9993    chr22   19023795-19109967      - |        9993
##    9994     chr6   90539619-90584155      + |        9994
##    9997    chr22   50961997-50964905      - |        9997
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
```

It returns a `GRanges` object with the chromosome name, start, end, strand and the gene_id. Note that `gene_id` here is the ENTREZ ID.

However, to be accurate, we want the exons, not the whole genes which contain introns. Let's get the exons:


```r
exons<- exonsBy(txdb, by = "gene")
exons
```

```
## GRangesList object of length 23459:
## $`1`
## GRanges object with 15 ranges and 2 metadata columns:
##        seqnames            ranges strand |   exon_id   exon_name
##           <Rle>         <IRanges>  <Rle> | <integer> <character>
##    [1]    chr19 58858172-58858395      - |    250809        <NA>
##    [2]    chr19 58858719-58859006      - |    250810        <NA>
##    [3]    chr19 58859832-58860494      - |    250811        <NA>
##    [4]    chr19 58860934-58862017      - |    250812        <NA>
##    [5]    chr19 58861736-58862017      - |    250813        <NA>
##    ...      ...               ...    ... .       ...         ...
##   [11]    chr19 58868951-58869015      - |    250821        <NA>
##   [12]    chr19 58869318-58869652      - |    250822        <NA>
##   [13]    chr19 58869855-58869951      - |    250823        <NA>
##   [14]    chr19 58870563-58870689      - |    250824        <NA>
##   [15]    chr19 58874043-58874214      - |    250825        <NA>
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
## 
## $`10`
## GRanges object with 2 ranges and 2 metadata columns:
##       seqnames            ranges strand |   exon_id   exon_name
##          <Rle>         <IRanges>  <Rle> | <integer> <character>
##   [1]     chr8 18248755-18248855      + |    113603        <NA>
##   [2]     chr8 18257508-18258723      + |    113604        <NA>
##   -------
##   seqinfo: 93 sequences (1 circular) from hg19 genome
## 
## ...
## <23457 more elements>
```

This returns a GRangesList object and each element of the list is a GRanges containing all the exons for that gene.

Let’s calculate the total exon lengths for each gene by the `width` function:


```r
# width of exons per gene
width(exons)
```

```
## IntegerList of length 23459
## [["1"]] 224 288 663 1084 282 297 273 270 36 96 65 335 97 127 172
## [["10"]] 101 1216
## [["100"]] 326 103 130 65 102 72 128 116 144 123 62 161
## [["1000"]] 1394 165 140 234 234 143 254 186 138 173 145 156 147 227 106 112 519
## [["10000"]] 218 68 5616 50 103 88 215 129 123 69 66 132 145 112 126 158 41
## [["100008586"]] 92 121 126 127
## [["100009676"]] 2784
## [["10001"]] 704 28 116 478 801 109 130 83 92 160 52
## [["10002"]] 308 127 104 222 176 201 46 106 918 705
## [["10003"]] 191 112 187 102 126 187 94 1193 99 ... 64 68 92 91 265 82 93 1054
## ...
## <23449 more elements>
```


```r
# sum the exons width up per gene
head(sum(width(exons)))
```

```
##         1        10       100      1000     10000 100008586 
##      4309      1317      1532      4473      7459       466
```

Note that the `width` and `su`m functions are vectorized. It will calculate across all genes.

Let’s turn it into a tibble using the `enframe` function:


```r
exon_len<- sum(width(exons)) %>%
      tibble::enframe(name = "ENTREZID", value = "exon_length")

head(exon_len)
```

```
## # A tibble: 6 x 2
##   ENTREZID  exon_length
##   <chr>           <int>
## 1 1                4309
## 2 10               1317
## 3 100              1532
## 4 1000             4473
## 5 10000            7459
## 6 100008586         466
```

Next, let’s map the ENTREZID to the official gene symbol so we can match the rownames of the RNAseq count matrix. We will need the `org.Hs.eg.db` Bioconductor package (install it if you
do not have it).


```r
library(org.Hs.eg.db)

map<- AnnotationDbi::select(org.Hs.eg.db, 
                            keys = exon_len$ENTREZID, 
                            columns= "SYMBOL",
                            keytype = "ENTREZID")

head(map)
```

```
##    ENTREZID  SYMBOL
## 1         1    A1BG
## 2        10    NAT2
## 3       100     ADA
## 4      1000    CDH2
## 5     10000    AKT3
## 6 100008586 GAGE12F
```

### join the exon length table

Read this [article](https://dplyr.tidyverse.org/reference/mutate-joins.html) to understand different join functions in `dplyr`.


```r
map<- left_join(exon_len, map)
head(map)
```

```
## # A tibble: 6 x 3
##   ENTREZID  exon_length SYMBOL 
##   <chr>           <int> <chr>  
## 1 1                4309 A1BG   
## 2 10               1317 NAT2   
## 3 100              1532 ADA    
## 4 1000             4473 CDH2   
## 5 10000            7459 AKT3   
## 6 100008586         466 GAGE12F
```

One of the key problems with genomics is that gene IDs are not always 1:1 mappable. Different versions of the genome (hg19 vs hg38 for humans) may have slightly different gene symbols.


```r
table(rownames(raw_counts_mat) %in% map$SYMBOL)
```

```
## 
## FALSE  TRUE 
## 20559 23250
```

### what genes are not in the mapping table?


```r
base::setdiff(rownames(raw_counts_mat), map$SYMBOL) %>%
  head(n = 20)
```

```
##  [1] "MIR6859-1"    "MIR1302-2HG"  "MIR1302-2"    "FAM138A"      "LOC100996442"
##  [6] "DDX11L17"     "WASH9P"       "MIR6859-2"    "LOC107985721" "LOC112268260"
## [11] "LOC100132287" "LOC105378947" "LOC101928626" "MIR12136"     "LINC01409"   
## [16] "FAM87B"       "LOC107984850" "LOC284600"    "LOC107985728" "LOC100288175"
```

Most of the differences are from non-coding RNA (LOC genes) or microRNAs. Many of those genes have a limited number of counts, we can ignore them for the moment.


```r
not_in_map<- setdiff(rownames(raw_counts_mat), map$SYMBOL)

raw_counts_mat[not_in_map, ] %>%
  head(n = 15)
```

```
##              01_SW_sgCTRL_Norm 02_SW_sgCTRL_Norm 11_SW_sgCTRL_Hyp
## MIR6859-1                    5                 1                8
## MIR1302-2HG                  0                 0                1
## MIR1302-2                    0                 0                0
## FAM138A                      0                 0                0
## LOC100996442                 9                 3               17
## DDX11L17                     0                 0                0
## WASH9P                      52                32               68
## MIR6859-2                    0                 0                0
## LOC107985721                 0                 0                0
## LOC112268260                 0                 0                0
## LOC100132287                 0                 0                0
## LOC105378947                 0                 0                0
## LOC101928626                 0                 0                0
## MIR12136                     3                 1                0
## LINC01409                    7                 8               19
##              12_SW_sgCTRL_Hyp
## MIR6859-1                  12
## MIR1302-2HG                 0
## MIR1302-2                   0
## FAM138A                     0
## LOC100996442               15
## DDX11L17                    0
## WASH9P                     76
## MIR6859-2                   0
## LOC107985721                0
## LOC112268260                0
## LOC100132287                0
## LOC105378947                0
## LOC101928626                0
## MIR12136                    0
## LINC01409                  24
```

subset only the common genes for the map file and the count matrix. Make sure the order of the genes is the same for both data.


```r
common_genes<- intersect(rownames(raw_counts_mat), map$SYMBOL)

## select only the common genes and re-order them by common_genes
map<- map %>%
  dplyr::slice(match(common_genes, SYMBOL))

# subset the common genes and re-order them by common_genes
raw_counts_mat<- raw_counts_mat[common_genes, ]

head(map)
```

```
## # A tibble: 6 x 3
##   ENTREZID  exon_length SYMBOL   
##   <chr>           <int> <chr>    
## 1 100287102        2838 DDX11L1  
## 2 653635           8050 WASH7P   
## 3 79501             918 OR4F5    
## 4 729737           5474 LOC729737
## 5 729759           1878 OR4F29   
## 6 81399             939 OR4F16
```

The order of the genes is the same for `map` and `raw_counts_mat`.


```r
head(raw_counts_mat)
```

```
##           01_SW_sgCTRL_Norm 02_SW_sgCTRL_Norm 11_SW_sgCTRL_Hyp 12_SW_sgCTRL_Hyp
## DDX11L1                   0                 0                0                0
## WASH7P                   18                11               23               45
## OR4F5                     0                 0                0                0
## LOC729737                 3                 3               16               16
## OR4F29                    0                 0                0                0
## OR4F16                    1                 0                0                3
```

## Normalizing Raw Counts to Transcripts per Million (TPM)

TPM normalization is essential for comparing gene expression levels across different samples and genes. We will write a function in the R programming language to perform this conversion and explain the steps involved.

### The `count2tpm` Function:

We will create a function called count2tpm in R to perform TPM normalization. This function takes two arguments: a count matrix and a vector of exon lengths. Let's break down the code step by step.


```r
count2tpm <- function(count_matrix, exon_length) {
  # Calculate reads per base pair per gene
  reads_per_bp_gene <- count_matrix / exon_length
  
  # Calculate the total reads per base pair for each sample
  reads_per_bp_sample <- colSums(reads_per_bp_gene)
  
  # Normalize to the library size and calculate TPM
  tpm_matrix <- t(t(reads_per_bp_gene) / reads_per_bp_sample) * 1000000
  return(tpm_matrix)
}
```

1. We start by defining the `count2tpm` function, which takes two arguments: `count_matrix` (raw gene expression counts) and `exon_length` (a vector of gene exon lengths).

2. We calculate the number of reads per base pair for each gene by dividing the count matrix by the exon length vector. This step helps us account for gene length differences.

3. We sum the reads per base pair values for each sample (column-wise) to calculate the total reads per base pair for each sample. This is crucial for library size normalization.

4. To normalize the data to library size, we divide the transposed `reads_per_bp_gene` matrix by the `reads_per_bp_sample` vector. The transposition allows us to perform element-wise division efficiently. Finally, we multiply the result by 1,000,000 to obtain TPM values.

### Applying the Function

Now, let's apply the `count2tpm` function to our raw count matrix and exon length vector. Here's how you can do it:


```r
tpm <- count2tpm(raw_counts_mat, map$exon_length)

head(tpm)
```

```
##           01_SW_sgCTRL_Norm 02_SW_sgCTRL_Norm 11_SW_sgCTRL_Hyp 12_SW_sgCTRL_Hyp
## DDX11L1          0.00000000        0.00000000        0.0000000        0.0000000
## WASH7P           0.31352121        0.24557646        0.4016983        0.7781915
## OR4F5            0.00000000        0.00000000        0.0000000        0.0000000
## LOC729737        0.07684343        0.09849323        0.4109445        0.4068975
## OR4F29           0.00000000        0.00000000        0.0000000        0.0000000
## OR4F16           0.14932231        0.00000000        0.0000000        0.4447598
```

These values represent the TPM-normalized gene expression levels for each gene in different samples.

### Conclusions

In this lesson, we have learned how to normalize raw gene expression counts to Transcripts per Million (TPM) using the `count2tpm` function in R. This normalization is crucial for comparing gene expression levels accurately across samples and genes, taking into account library size and gene length.

## Analyzing Gene Expression Data Using t-Tests

Gene expression data provides valuable insights into how genes are activated or deactivated under different conditions, such as in response to diseases or environmental changes.

A t-test is a statistical test that helps us determine whether there is a significant difference between the means of two groups. In the context of gene expression analysis, we can use t-tests to identify genes that are differentially expressed between two experimental conditions. For example, we might want to know which genes are upregulated or downregulated in response to hypoxia (low oxygen levels) compared to normoxia (normal oxygen levels).

### Hypothesis Testing

Before we dive into the code, let's understand the key components of hypothesis testing:

1. **Null Hypothesis (H0)**: This is the default assumption that there is no significant difference between the groups. In gene expression analysis, it means that the gene is not differentially expressed.

2. **Alternative Hypothesis (Ha)**: This is the hypothesis we want to test. It suggests that there is a significant difference between the groups. In gene expression analysis, it implies that the gene is differentially expressed.

3. **p-value**: The p-value represents the probability of observing the data, assuming that the null hypothesis is true. A small p-value (typically less than 0.05) suggests that we can reject the null hypothesis and accept the alternative hypothesis.

### Analyzing Specific Genes

Now, let's use t-tests to analyze the expression of specific genes and understand how they respond to hypoxia.

We'll start by examining the gene `WASH7P`. We perform a t-test to compare its expression levels between normoxia and hypoxia samples.


```r
t.test(tpm["WASH7P", c(1,2)], tpm["WASH7P", c(3,4)])
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  tpm["WASH7P", c(1, 2)] and tpm["WASH7P", c(3, 4)]
## t = -1.6227, df = 1.0651, p-value = 0.3404
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -2.416652  1.795860
## sample estimates:
## mean of x mean of y 
## 0.2795488 0.5899449
```

In this case, the p-value is 0.3404, suggesting that there is no significant difference in the expression of the WASH7P gene between normoxia and hypoxia.

Next, we examine the `VEGFA` gene, which is known to be a key regulator of angiogenesis in response to hypoxia.


```r
t.test(tpm["VEGFA", c(1,2)], tpm["VEGFA", c(3,4)])
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  tpm["VEGFA", c(1, 2)] and tpm["VEGFA", c(3, 4)]
## t = -31.953, df = 1.8939, p-value = 0.00132
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -56.58778 -42.50235
## sample estimates:
## mean of x mean of y 
##  13.58050  63.12557
```

Here, the p-value is 0.00132, indicating a significant difference in VEGFA expression between normoxia and hypoxia. This suggests that VEGFA is likely upregulated under hypoxic conditions.

Now, let's analyze the `SLC2A1` (GLUT1) gene, which plays a role in glucose transport during anaerobic glycolysis.


```r
t.test(tpm["SLC2A1", c(1,2)], tpm["SLC2A1", c(3,4)])
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  tpm["SLC2A1", c(1, 2)] and tpm["SLC2A1", c(3, 4)]
## t = -31.938, df = 1.1196, p-value = 0.01354
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -2460.877 -1295.615
## sample estimates:
## mean of x mean of y 
##  1360.405  3238.651
```

The p-value is 0.01354, indicating a significant difference in `SLC2A1` expression between the two conditions. This suggests that `SLC2A1` may be upregulated under hypoxic conditions as well.

### Analyzing All Genes
To analyze all genes in our dataset, we can create a custom function called `mytest` that performs t-tests for each gene pair and extracts the p-values.


```r
mytest <- function(x) t.test(x[c(1,2)], x[c(3,4)], var.equal = TRUE)$p.value
pvals <- apply(tpm, 1, mytest)

head(pvals)
```

```
##     DDX11L1      WASH7P       OR4F5   LOC729737      OR4F29      OR4F16 
##         NaN 0.246130682         NaN 0.001173022         NaN 0.593224964
```

Here, we apply the `mytest` function to each row (gene) in our gene expression dataset (`tpm`) to calculate p-values.

Finally, we count how many genes have p-values smaller than 0.01 to identify differentially expressed genes:


```r
sum(pvals < 0.01, na.rm = TRUE)
```

```
## [1] 3378
```

`pvals< 0.01` returns a logical vector of TRUE and FALSE. TRUE is 1 and FALSE is 0 under the hood in R. If you sum them up `sum(pvals< 0.01, na.rm = TRUE)` will tell you how many TRUEs are in the vector.

There are 3378 genes with p-values smaller than 0.01!

### Conclusion

In this lesson, we learned how to use t-tests to analyze gene expression data and identify differentially expressed genes. We examined specific genes and performed t-tests, understanding the significance of p-values and hypothesis testing. Additionally, we applied t-tests to all genes in the dataset to identify potential candidates for further investigation.

## Analyzing Gene Expression Data with ggplot2

In this lesson, we will explore how to analyze gene expression data using the powerful ggplot2 library. We will focus on visualizing p-value distributions, comparing differentially expressed genes, and creating boxplots to gain insights into gene expression changes under different conditions.

### Import Libraries and Load Data

First, let's import the necessary libraries and load your gene expression data. In this lesson, we assume you have a dataset containing gene names and p-values representing their significance.


```r
# Import required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load your p-values data (replace 'pvals' with your actual dataset)
pval_df <- pvals %>%
  tibble::enframe(name = "gene", value = "pvalue")

# Display the first few rows of the p-value data
head(pval_df)
```

```
## # A tibble: 6 x 2
##   gene         pvalue
##   <chr>         <dbl>
## 1 DDX11L1   NaN      
## 2 WASH7P      0.246  
## 3 OR4F5     NaN      
## 4 LOC729737   0.00117
## 5 OR4F29    NaN      
## 6 OR4F16      0.593
```

This code converts your gene names and corresponding p-values into a data frame, making it easier to work with in ggplot2.

### Visualizing P-Value Distribution
Next, let's create a histogram to visualize the distribution of p-values:


```r
# Create a histogram of p-values
ggplot(pval_df, aes(x = pvalue)) +
  geom_histogram(color = "white") +
  theme_bw(base_size = 14) +
  ggtitle("p-value distribution")
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-35-1.pdf)<!-- --> 

This code uses ggplot2 to create a histogram, providing insights into the distribution of p-values across your genes. Understanding p-value distribution can help assess the significance of your results.

### Identifying Differentially Expressed Genes
Now, we'll focus on comparing gene expression between hypoxia and normoxia conditions, specifically looking at up-regulated genes. We'll start by calculating the average expression levels for both conditions and identifying the up-regulated genes.


```r
# Calculate average expression for normoxia and hypoxia conditions
avg_normoxia <- rowMeans(tpm[, c(1, 2)])
avg_hypoxia <- rowMeans(tpm[, c(3, 4)])

# Identify up-regulated genes
up_genes <- (avg_hypoxia - avg_normoxia) > 0

# Get the names of up-regulated genes
up_gene_names <- rownames(tpm)[up_genes]

head(up_gene_names)
```

```
## [1] "WASH7P"       "LOC729737"    "OR4F16"       "LOC100288069" "LINC02593"   
## [6] "SAMD11"
```

Here, we calculate the average expression for normoxia and hypoxia conditions and then identify up-regulated genes by comparing the averages. `up_gene_names` contains the names of these genes.

### Selecting Differentially Expressed Genes
Not all up-regulated genes may be significantly different. Let's find the intersection of up-regulated genes with those that have a p-value less than 0.01 (significant up-regulation):


```r
# Select differentially expressed genes (intersection of up-regulated genes and significant p-values)
differential_genes <- pvals[pvals < 0.01 & !is.na(pvals)] %>%
  names()

# Find the intersection
differential_up_genes <- intersect(differential_genes, up_gene_names)

length(differential_up_genes)
```

```
## [1] 1990
```

`differential_up_genes` now contains the names of genes that are both up-regulated and significantly different under hypoxia.

### boxplot to visualize gene expression changes between the two conditions.

#### preparing the data 

Now, let's prepare our data for creating a boxplot to visualize gene expression changes between the two conditions. We will convert the gene expression data into a long format suitable for `ggplot2`:


```r
# Convert gene expression data to long format
tpm[differential_genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  tidyr::pivot_longer(-1, names_to = "sample", values_to = "tpm")
```

```
## # A tibble: 13,512 x 3
##    gene      sample                 tpm
##    <chr>     <chr>                <dbl>
##  1 LOC729737 01_SW_sgCTRL_Norm   0.0768
##  2 LOC729737 02_SW_sgCTRL_Norm   0.0985
##  3 LOC729737 11_SW_sgCTRL_Hyp    0.411 
##  4 LOC729737 12_SW_sgCTRL_Hyp    0.407 
##  5 NOC2L     01_SW_sgCTRL_Norm 128.    
##  6 NOC2L     02_SW_sgCTRL_Norm 124.    
##  7 NOC2L     11_SW_sgCTRL_Hyp   73.4   
##  8 NOC2L     12_SW_sgCTRL_Hyp   70.9   
##  9 PERM1     01_SW_sgCTRL_Norm   0.176 
## 10 PERM1     02_SW_sgCTRL_Norm   0.113 
## # i 13,502 more rows
```

This code converts the gene expression data into a long format with columns for gene names, sample names, and expression values.

Add another column to denote the condition by separating the sample column to two columns: sample and condition


```r
tpm_df<- tpm[differential_genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  tidyr::pivot_longer(-1, names_to = "sample", values_to = "tpm") %>%
  tidyr::separate(sample, into = c("sample", "condition"), sep = "_sgCTRL_")

head(tpm_df)
```

```
## # A tibble: 6 x 4
##   gene      sample condition      tpm
##   <chr>     <chr>  <chr>        <dbl>
## 1 LOC729737 01_SW  Norm        0.0768
## 2 LOC729737 02_SW  Norm        0.0985
## 3 LOC729737 11_SW  Hyp         0.411 
## 4 LOC729737 12_SW  Hyp         0.407 
## 5 NOC2L     01_SW  Norm      128.    
## 6 NOC2L     02_SW  Norm      124.
```

#### Customizing Boxplot Order
By default, ggplot2 orders boxplots alphabetically. To change the order, convert the condition column to a factor with the desired order:


```r
# Define the order for boxplot
tpm_df$condition <- factor(tpm_df$condition, levels = c("Norm", "Hyp"))
```
This code ensures that the boxplot orders the conditions as "Norm" followed by "Hyp."

### Creating the Boxplot

Now, we can create the boxplot to visualize gene expression changes between the two conditions:


```r
# Create the boxplot
ggplot(tpm_df, aes(x = condition, y = log2(tpm + 1))) +
  geom_boxplot() +
  theme_bw(base_size = 14)
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-41-1.pdf)<!-- --> 

This code uses `ggplot2` to create a boxplot, showing the distribution of gene expression values between the "Norm" and "Hyp" conditions. The `log2(tpm + 1)` transformation is often used to visualize RNA-Seq data.

### Visualizing Raw Expression Values
Sometimes, it's essential to examine the raw expression values to identify outliers. Here, we use a boxplot to visualize the raw values of selected genes:


```r
# Select specific genes (e.g., "VEGFA" and "SLC2A1") for visualization
ggplot(tpm_df %>%
         filter(gene %in% c("VEGFA", "SLC2A1")), 
       aes(x = condition, y = log2(tpm + 1))) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(~ gene) + 
  theme_bw(base_size = 14)
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-42-1.pdf)<!-- --> 

This code creates scatter plots for selected genes and overlays boxplots to help visualize
the distribution of the gene expression levels.

### Conclusion
In this lesson, we covered the entire process of analyzing gene expression data using ggplot2, from loading data to visualizing differential expression. Understanding the steps involved and customizing plots can provide valuable insights into your gene expression analysis.

## Correcting for Multiple Comparisons in Statistical Analysis

This is a critical step when conducting hypothesis tests on a large number of data points, such as in genomics research. We will cover the need for correction, different methods to control errors, and demonstrate how to implement one of the widely-used methods, the False Discovery Rate (FDR) correction.

### Why Correct for Multiple Comparisons?

Imagine you are a scientist studying the gene expression levels of thousands of genes in response to a treatment. You perform statistical tests to identify which genes are significantly differentially expressed. If you run these tests without correction, you are likely to encounter a problem known as the "multiple comparisons problem." In essence, the more tests you perform, the higher the chance of obtaining false positives (i.e., incorrectly identifying genes as significant).

To address this issue, we need correction methods that control the probability of making at least one false discovery while testing multiple hypotheses. In this project, we will focus on the `False Discovery Rate (FDR)` correction method, specifically using the `Benjamini & Hochberg (BH)` procedure.

### The Multiple Comparisons Example
Let's illustrate the concept with a practical example. Suppose you have gene expression data from 23,250 genes and you want to identify those that are differentially expressed between two conditions (e.g., control and treatment). You perform a statistical test for each gene and obtain p-values.


```r
# Sample code to generate random p-values for demonstration purposes
m <- 23250  # Number of genes
n <- 100    # Number of comparisons
randomData <- matrix(rnorm(n * m), m, n)
nullpvalues <- apply(randomData, 1, mytest)  # Simulated p-values
hist(nullpvalues)
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-43-1.pdf)<!-- --> 

If you were to plot the histogram of these p-values, you might expect them to follow a uniform distribution (a flat line) under the null hypothesis (no differential expression). However, due to the nature of p-values as random variables, you will still observe some p-values below the commonly used significance level of 0.05, even when no genes are differentially expressed.

Compare this histogram with the histogram for the real data. what do you see? Even if we randomly generated the data, you still see some p values are smaller than 0.05!! We randomly generated data, there should be No genes that deferentially expressed. However, we see a flat line across different p values.

p values are random variables. Mathematically, one can demonstrate that under the null hypothesis (and some assumptions are met, in this case, the test statistic T follows standard normal distribution), p-values follow a uniform (0,1) distribution, which means that P(p < p1) = p1.

This means that the probability we see a p value smaller than p1 is equal to p1. That being said, with a 100 t-tests, under the null (no difference between control and treatment), we will see 1 test with a p value smaller than 0.01. And we will see 2 tests with a p value smaller than 0.02 etc.

We have 23250 genes in the matrix, and we did 23250 comparisons at one time. This explains why we see `23250 * 0.05 = 1162` p-values are smaller than 0.05 in our randomly generated numbers. That’s exactly what we see in the null distribution of the p-values.

>In fact, checking the p-value distribution by histogram is a very important step during data analysis. You may want to read a blog post by David Robinson: [How to interpret a p-value histogram](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/).

### Correcting for Multiple Comparisons: False Discovery Rate (FDR)
How do we control the false positives for multiple comparisons? One way is to use the Bonferroni correction to correct the familywise error rate (FWER): define a particular comparison as statistically significant only when the P value is less than alpha(often 0.05) divided by the number (m) of comparisons (p < alpha/m).

Say we computed 100 t-tests, and got 100 p values, we only consider the genes with a p value smaller than 0.05/100 as significant. This approach is very conservative and is used in Genome-wide association studies (GWAS). Since we often compare millions of genetic variations between (tens of thousands) cases and controls, this threshold will be very small!

Alternatively, we can use False Discovery Rate (FDR) to report the gene list. FDR = #false positives/# called significant. This approach does not use the term statistically significant but instead use the term discovery. Let’s control FDR for a gene list with FDR = 0.05. It means that of all the discoveries, 5% of them is expected to be false positives.

Benjamini & Hochberg (BH method) in 1995 proposed a way to control FDR: Let k be the largest i such that `p(i) <= (i/m) * alpha`, (m is the number of comparisons) then reject H(i) for i =1, 2, …k

This process controls the FDR at level alpha. The method sets a different threshold p value for each comparison. Say we computed 100 t-tests, and got 100 p values, and we want to control the FDR =0.05. We then rank the p values from small to big. if `p(1) <= 1/100 * 0.05`, we then reject null hypothesis and accept the alternative. if `p(2) < = 2/100 * 0.05`, we then reject the null and accept the alternative.


```r
#remove the NAs
pvals<- pvals[!is.na(pvals)]

## order the pvals computed above and plot it.
alpha<- 0.05

#m is the number of comparisons 
m<- length(pvals)

# let's arrange the p-value from small to big and get only the first 5000  
top_5000_pvalue<- pval_df %>%
  dplyr::arrange(pvalue) %>%
  mutate(rank = row_number()) %>%
  dplyr::slice(1:5000)

head(top_5000_pvalue)  
```

```
## # A tibble: 6 x 3
##   gene         pvalue  rank
##   <chr>         <dbl> <int>
## 1 PPFIA4  0.000000907     1
## 2 GINS2   0.00000173      2
## 3 ZFP36L2 0.00000197      3
## 4 ECHS1   0.00000297      4
## 5 MOB3A   0.00000392      5
## 6 FLYWCH2 0.00000661      6
```

let's plot


```r
ggplot(top_5000_pvalue, aes(x= rank, y = pvalue))+
  geom_point() + 
  geom_abline(slope = alpha/m, intercept = 0, color = "red", linetype = 2) 
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-45-1.pdf)<!-- --> 

p values that are below the red dotted line are controlled at FDR of 0.05.

We will use p.adjust function and the method “fdr” or “BH” to correct the p value, what the p.adjust function does is to recalculate the p-values.

With the FDR definition, p value is only significant if `p(i)<= (i/m) * alpha` We can rewrite it to `p(i) * m/i <= alpha`. The p.adjust function returns `p(i) * m/i` the adjusted p-value. We can then only accept the returned the p values if `p.adjust(pvals) <= alpha`.


```r
top_5000_pvalue %>%
  mutate(padj = pvalue * m/rank) %>%
  head()
```

```
## # A tibble: 6 x 4
##   gene         pvalue  rank   padj
##   <chr>         <dbl> <int>  <dbl>
## 1 PPFIA4  0.000000907     1 0.0161
## 2 GINS2   0.00000173      2 0.0154
## 3 ZFP36L2 0.00000197      3 0.0116
## 4 ECHS1   0.00000297      4 0.0132
## 5 MOB3A   0.00000392      5 0.0139
## 6 FLYWCH2 0.00000661      6 0.0196
```

How many of those p-values are below the dotted red line?


```r
top_5000_pvalue %>%
  mutate(padj = pvalue * m/rank) %>%
  filter(padj <= alpha) %>%
  filter(rank == which.max(rank))
```

```
## # A tibble: 1 x 4
##   gene      pvalue  rank   padj
##   <chr>      <dbl> <int>  <dbl>
## 1 RNASEH2A 0.00893  3173 0.0500
```

There are total 3173 p-values that are significant after FDR correction.

We can verify it using the p.adjust function in R:


```r
adjusted_pvalues<- p.adjust(pvals, method="fdr")

sum(adjusted_pvalues < 0.05)
```

```
## [1] 3173
```

This is the same as what we calculated manually!

We can plot a vertical line on the p-value ranking plot:


```r
ggplot(top_5000_pvalue, aes(x= rank, y = pvalue))+
  geom_point() + 
  geom_abline(slope = alpha/m, intercept = 0, color = "red", linetype = 2) +
  geom_vline(xintercept = 3173, linetype = 2, color = "red")
```

![](11_Final_project_analyze_RNAseq_from_GEO_files/figure-latex/unnamed-chunk-49-1.pdf)<!-- --> 

### Conclusion
Correcting for multiple comparisons is essential when conducting statistical tests on a large number of hypotheses. The False Discovery Rate (FDR) correction, such as the Benjamini & Hochberg method, allows us to control the rate of false discoveries while identifying significant results. This approach is valuable in various scientific disciplines to ensure the reliability of statistical findings.
