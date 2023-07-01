
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LIRA

<!-- badges: start -->
<!-- badges: end -->

LIRA is an R package to predict the effectiveness of immunotherapy in
patients with NSCLC.

## 1.Introduction

1.LIRA was designed to predict the reponse to immunotherapy in patients
with NSCLC. 2.This package consists a random survival model based on
expression profiles of 50 genes, which were selected through univariate
Cox regression and subgroup analysis. 3.This package providers functions
for different scenarios and allows visual assessment to determine
whether patients are suitable for immunotherapy.

## 2.Installation

It is essential that you have R 3.6.3 or above already installed on your
computer or server. Before installing LIRA, please install all
dependencies by executing the following command in R console:

The dependencies
includes`<IOBR>`,`<crayon>`,`<ggplot2>`,`<randomForest>`,`<randomForestSRC>`.

#### Graphical abstract for construction and clinical application of LIRA

![LIRA logo](./man/figures/LIRA-workflow.png)

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c("IOBR", "crayon", "ggplot2", "randomForest", "randomForestSRC")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")
```

Then, you can start to install IOBR from github by typing the following
code into your R session:

``` r
if (!requireNamespace("LIRA", quietly = TRUE))
  remotes::install_github("LiaoWJLab/LIRA")
```

Load the IOBR package in your R session after the installation is
complete:

``` r
library(LIRA)
library(IOBR)
```

## 3.Usage

### Obtain dataset from TCGA using UCSCXenaTools R package.

For transcriptomic data of TCGA data sets, we strongly recommend user to
use UCSCXenaTools R package. Here, we download counts data of TCGA-LUAD
from UCSC using UCSCXenaTools R package

``` r
if (!requireNamespace("UCSCXenaTools", quietly = TRUE))
    BiocManager::install("UCSCXenaTools")
library(UCSCXenaTools)
#> Warning: 程辑包'UCSCXenaTools'是用R版本4.2.1 来建造的
#> =========================================================================================
#> UCSCXenaTools version 1.4.8
#> Project URL: https://github.com/ropensci/UCSCXenaTools
#> Usages: https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html
#> 
#> If you use it in published research, please cite:
#> Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
#>   from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq.
#>   Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627
#> =========================================================================================
#>                               --Enjoy it--

eset_LUAD<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Lung Adenocarcinoma (LUAD)") %>% 
  XenaFilter(filterDatasets    = "TCGA-LUAD.htseq_counts.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()
#> This will check url status, please be patient.
#> All downloaded files will under directory C:\Users\inter\AppData\Local\Temp\RtmpEt6GED.
#> The 'trans_slash' option is FALSE, keep same directory structure as Xena.
#> Creating directories for datasets...
#> Downloading TCGA-LUAD.htseq_counts.tsv.gz
#> Warning in download.file(url, destfile, ...): downloaded length 2384423 !=
#> reported length 66869176
#> Warning in download.file(url, destfile, ...): URL 'https://gdc-hub.s3.us-
#> east-1.amazonaws.com:443/download/TCGA-LUAD.htseq_counts.tsv.gz': Timeout of 60
#> seconds was reached
#> ==> Trying #2
#> Warning in download.file(url, destfile, ...): downloaded length 887335 !=
#> reported length 66869176

#> Warning in download.file(url, destfile, ...): URL 'https://gdc-hub.s3.us-
#> east-1.amazonaws.com:443/download/TCGA-LUAD.htseq_counts.tsv.gz': Timeout of 60
#> seconds was reached
#> ==> Trying #3
#> Warning in download.file(url, destfile, ...): downloaded length 725434 !=
#> reported length 66869176

#> Warning in download.file(url, destfile, ...): URL 'https://gdc-hub.s3.us-
#> east-1.amazonaws.com:443/download/TCGA-LUAD.htseq_counts.tsv.gz': Timeout of 60
#> seconds was reached
#> Tried 3 times but failed, please check your internet connection!

# Remove the version numbers in Ensembl ID.
eset_LUAD$Ensembl_ID<-substring(eset_LUAD$Ensembl_ID, 1, 15)
eset_LUAD<-column_to_rownames(eset_LUAD, var = "Ensembl_ID")
# Revert back to original format because the data from UCSC was log2(x+1)transformed.
eset_LUAD<-(2^eset_LUAD)-1
head(eset_LUAD[1:5,1:5])
#>                 TCGA-97-7938-01A TCGA-55-7574-01A TCGA-05-4250-01A
#> ENSG00000000003             2032             1000             5355
#> ENSG00000000005               15                0                5
#> ENSG00000000419             1220              744             2898
#> ENSG00000000457              876              560              734
#> ENSG00000000460              250              325              785
#>                 TCGA-55-6979-11A TCGA-95-A4VK-01A
#> ENSG00000000003              516             2269
#> ENSG00000000005                0                1
#> ENSG00000000419              589              819
#> ENSG00000000457              361             1477
#> ENSG00000000460               88              327
```

### Annotate genes in expression matrix and remove duplicate genes using IOBR package.

``` r
eset_LUAD <- anno_eset(eset = eset_LUAD, annotation = anno_grch38, probe = "id", symbol = "symbol")
#> Row number of original eset:
#> >>>>  244
#> >>> 100.00% of probe in expression set was annotated
#> Row number after filtering duplicated gene symbol:
#> >>>>  244
```

## Calculate riskscore of TCGA

We offer two methods to calculate the riskscore for samples in different
scenarios. One method involves calculating the score for individual
samples by inputting them into the LIRA algorithm, while the other
method involves calculating the score for multiple samples by inputting
them together. Users can choose the appropriate method based on your
needs, using the `<loop>` parameter. We utilize data from the OAK and
POPLAR studies as reference data to help users understand the
distribution of scores for their samples.

``` r
# reference
head(ref_lira_score)
#>              ID riskscore riskscore_all     cohort score2               id
#> 1 EA-002d8e81a8  26.25412     215.66614 validation   high PAT-a11bdfce3174
#> 2 EA-00b68d2219  15.12424      82.07215   training    low PAT-fdf7aaeefdc7
#> 3 EA-01813c5ad1  15.11904     135.49451   training    low PAT-17d3715fdbda
#> 4 EA-022952077d  15.14425     160.59087   training    low PAT-d0e1de13f4e8
#> 5 EA-02352606b6  29.34257     222.25082   training   high PAT-df63e2224d1d
#> 6 EA-033a6dffad  17.04630     119.55705   training    low PAT-5020088d94e4
#>     UNI_TRNC    ACTARM         HIST OS_months OS_status PFS_months PFS_status
#> 1 a11bdfce31 MPDL3280A NON-SQUAMOUS 21.585216         0  20.073922          0
#> 2 fdf7aaeefd MPDL3280A NON-SQUAMOUS 21.347952         1   8.499647          1
#> 3 17d3715fdb MPDL3280A NON-SQUAMOUS  5.026694         1   1.182752          1
#> 4 d0e1de13f4 MPDL3280A NON-SQUAMOUS  7.063655         1   1.412731          1
#> 5 df63e2224d MPDL3280A NON-SQUAMOUS  6.028820         1   2.800271          1
#> 6 5020088d94 MPDL3280A NON-SQUAMOUS  3.942505         0   3.942505          1
#>   BOR STUDYNAME PDL1_TC_22C3 tTMB STK11_status KEAP1_status EGFR_status
#> 1  PR       OAK               <16           WT           WT          WT
#> 2  SD    POPLAR         <NA> <NA>         <NA>         <NA>        <NA>
#> 3  PD       OAK               <16           WT           WT          WT
#> 4  PD       OAK       [1,50)  <16           WT           WT          WT
#> 5  SD    POPLAR         <NA> <NA>         <NA>         <NA>        <NA>
#> 6  SD       OAK       [1,50)  <16           WT          MUT          WT
#>   trunc_anonymized_patient_id title description bioSampleId caseOrControl
#> 1            PAT-a11bdfce3174    NA          NA          NA            NA
#> 2                        <NA>    NA          NA          NA            NA
#> 3            PAT-17d3715fdbda    NA          NA          NA            NA
#> 4            PAT-d0e1de13f4e8    NA          NA          NA            NA
#> 5                        <NA>    NA          NA          NA            NA
#> 6            PAT-5020088d94e4    NA          NA          NA            NA
#>   gender organismPart cellLine region phenotype          id2 ARM  group4
#> 1   male           NA       NA     NA     NSCLC a11bdfce3174  IO high+IO
#> 2   male           NA       NA     NA     NSCLC fdf7aaeefdc7  IO  low+IO
#> 3   male           NA       NA     NA     NSCLC 17d3715fdbda  IO  low+IO
#> 4 female           NA       NA     NA     NSCLC d0e1de13f4e8  IO  low+IO
#> 5   male           NA       NA     NA     NSCLC df63e2224d1d  IO high+IO
#> 6 female           NA       NA     NA     NSCLC 5020088d94e4  IO  low+IO
#>   treat_score
#> 1     high+IO
#> 2      low+IO
#> 3      low+IO
#> 4      low+IO
#> 5     high+IO
#> 6      low+IO

# 选择前5个样本，并命名为validation，输入计算
eset <- eset_LUAD[, 1:5]
colnames(eset)[1:5]<- paste0("validation", 1:5)

# 单个样本计算score
res_loop <- lira_model2(eset       = eset[,1:5],
                        ref        = T, 
                        loop       = FALSE,
                        check_eset = T, 
                        pdata      = NULL,
                        id_pdata   = "ID",
                        scale      = TRUE,
                        plot       = FALSE,
                        return_all = T)
#> >>>-- Predicting new data with combined gene expression data...
#> >>>-- Scaling data...
#> >>>-- Removing outlier genes...
#> [1] ">>> Is NA exist:  0"
#> [1] ">>>> Is nonnumeric variables exist ? >>>>"
#>    Mode   FALSE 
#> logical     240 
#> [1] ">>>> Is -Inf variables exist ? >>>>"
#>    Mode   FALSE 
#> logical     240 
#> [1] ">>>> Is Inf variables exist ? >>>>"
#>    Mode   FALSE 
#> logical     240 
#> [1] ">>> Variables with sd = 0 :  "
#>    Mode   FALSE 
#> logical     240 
#> >>>-- Predicting new data with LIRA model...
#> >>>-- Only 0.00% of model genes appear on gene matrix,
#>  interpret results with caution 
#> >>>-- DONE!
head(res_loop)
#>              ID riskscore     cohort
#> 1   validation1  35.92336 validation
#> 2   validation2  22.42322 validation
#> 3   validation3  28.35451 validation
#> 4   validation4  35.47925 validation
#> 5   validation5  21.19544 validation
#> 6 EA-d9d7aa94b3  22.53870  reference
```

#### Check the distribution of the samples in the reference data

``` r
# OS
lira_score_location(res_loop[1, ], 
                    pat_id       = "ID",
                    col_score    = "riskscore",
                    ref_score    = ref_lira_score,
                    palette      = "nrc",
                    cols         = NULL, 
                    palette_line = "jama",
                    panel        = "OS",
                    path         = "/man/figures")
```

#### Check the relative expression levels of each sample on the model genes via heatmap.

We used a subset of samples from OAK and POPLAR as references

``` r
lira_heatmap(eset            =  eset[,1:2],
             pdata           = NULL,
             id_pdata        = "ID",
             BOR             = "BOR", 
             palette_heatmap = 6,
             path            = "/man/figures")
```

## References

1.Zeng D, Fang Y, Chen G, …, Liao W; Predicting NSCLC immunotherapy
resistance based on Random Forest Model. (2023) *Under Review*.

## Reporting bugs

Please report bugs to the [Github issues
page](https://github.com/LiaoLab/LIRA/issues)

E-mail any questions to <dongqiangzeng0808@gmail.com>
