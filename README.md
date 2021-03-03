# UKBB_GWAS-pipeline
Python Pipeline for analysing UKBB data


## Introduction

The UK-Biobank genome wide association study (GWAS) pipeline has been optimized to perform GWAS on the imputed genetic dataset of the full 500,000 from UK Biobank quickly, efficiently and in a standardized manner.

We developed a flexible in-house software pipeline for automated analysis of UKBB genome-wide data set of quantitative/qualitative phenotypes. The software pipeline is implemented in Python and R in an object-oriented style, makes extensively use of the open source PLINK/BOLT libraries and runs on our Linux computer cluster with a multi-cores processing system. Parameters, options and data sets can be combined in a flexible way and are currently read from shell. An attached plotting routine facilitates a fast evaluation of summary statistics (e.g. Manhattan and quantile-quantile (Q-Q) plots). 


## Design principles


The goal of the UKB-GWAS pipeline is to provide a portable and robust pipeline for reproducble genome-wide association studies.

A GWAS requires a complex set of analyses with complex dependancies between the analyses. We want to support GWAS work by supporting

- reproducibility: we can rerun the entire analysis from start to finish;
- reusability:  we can run the entire analysis with different parameters in an efficient and consistent way;

We took the following into account:

- Users are heterogeneous both in terms of their bioinformatics needs and the computing environments they will use.
- Each GWAS is different -- it must be customisable to allow the bioinformaticists to set different parameters.


Users can run the pipeline from *~Programs/ukb_gwas.py*


## Phenotype and covariate files 


### IDs

The pipeline uses a filtered version of the UK Biobank genetic data and therefore both the phenotype and covariate files need to have the genetic ids of UK Biobank. 


### File Format

Users doesn't specify the location of Imputed genetic data, samples file, harcall set , ect... . Pipeline can load automatically BGEN/PGEN data and a set of default covariates.
Users must prepare phenotype and covariate files ( or a combine pheno/covar file)

- These files should be space delimited TAB files.
- The first two columns must be FID and IID (the PLINK identifiers of an individual); any number of columns may follow.
- Values in the column should be numeric.
- Case/control phenotypes should be encoded as 1=unaffected (control), 2=affected (case) [PLINK2] or 0=unaffected (control), 1=affected (case) [BOLT-LMM]
- Missing values in pheno/covariate files should be encoded as NA or -9. 


An example of a combined phenotype and covariate file:


```{}
FID IID phenotype1 covariate1 covariate2
123456 123456 100 1 1
234567 234567 20 2 1
345678 345678 50 3 NA

```




### Default Covariates

The pipeline load some covariates directly from ukb_sqc_v2.txt file. By default sex, genotyping array and first 40 principal component are loaded and ready to use.
Default covariates, if none are specified, for BOLT-LMM are **sex and genotyping array**; for PLINK, **sex, genotyping array and the first 10 principal components**. The pipeline automatically load default covariates, you don't need to include them into pheno/covar file but only flag *--cname* option ( see some examples below)



## UK Biobank Genetic Data for GWAS pipeline - Scenario 


The UK Biobank cohort includes a large number of participants from a wide range of ethnic backgrounds. Most of participants (88.26%) are 'British'. In addition, a 30% of participants are inferred to be related. Broadly speaking, 4 scenario were implemented:

- Scenario 1. Analyse All samples: a non homogeneus population( White British + others) with a certain degree of kinship. (**--model1**). GWAS perfomed by BOLT
- Scenario 2. Analyse All unrelated samples. (**--model2**). GWAS performed by PLINK
- Scenario 3. Analyse All 'White british' samples. (**--model3**). GWAS perfomed by BOLT
- Scenario 4. Analyse unrelated 'White British' samples only.(**--model4**). GWAS performed by PLINK

PLINK allows for an analysis to be performed in a homogeneous and unrelated population.
BOLT-LMM uses a linear mixed model (LMM) to account for both relatedness and population stratification, therefore allowing a wider range of individuals to be included in terms of relatedness and ancestry [https://www.nature.com/articles/ng.3190]() .


### BOLT-LMM binary trait analysis
 

- BOLT-LMM performs a *linear regression* and therefore the output betas will need to be transformed to obtain log odds ratios using the following formula:

```
logOR = beta_bolt/(mu*(1-mu)) , where mu=case/(case+control)

```

The standard errors are adjusted using the same method: 

```
se_bolt/(mu*(1-mu))
```

- case-control balance: for traits with a case fraction of at least 10%, BOLT-LMM test statistics are well-calibrated for SNPs with MAF>0.1%. For highly unbalanced case-control the authors recommend using SAIGE.




### Restriction to unrelated population: PLINK

Within PLINK, the association between SNPs and a binary/quantitative outcome (value **1** = unaffected and value **2** = affected; '0' and '-9' represent missing) has been restricted to the 'All' and 'White British' subsets and related individuals have been excluded. To test All or White-British  samples only use option *--model2* or *--model4* respectively).


### Association analysis: statistical methods


**PLINK**:
Genome-wide association analysis (GWAS) is conducted using linear/logistic regression, implemented using the software PLINKv2.00. Genotype array, sex and the first 10 PCs (out of 40) supplied by UKBiobank were fitted as covariates in the model.

**BOLT-LMM**:
To model population structure in the sample we used 126,508 directly genotyped SNPs, obtained after filtering on MAF > 0.05; genotyping rate > 95%; Hardy-Weinberg equilibrium p-value < 0.0001 and LD pruning to an r2 threshold of 0.2 using PLINK. Genotype array and sex were adjusted for in the model. BOLT-LMM association statistics are on the linear scale.


## User options

To get a list of basic options, run:


python3 /GWASdir/ukb_gwas.py -h


```{}

usage: UKB-GWAS [-h] -p PHENO -ph PHNAME [-c COVAR] [-cname COVARNAME] -o
                OUTPUTFILE [-plot GWASPLOT] [-chr CHR] [-maf MAF] [-pca PCA]
                [-rsq {0.3,0.6,0.8}] [-cpu CORES] [--adj {yes,no}] -m
                [{model1,model2,model3,model4}] [-v]

A python script to run GWAS analysis of Binary/Quantitative phenotypes using
UKBB imputed data

optional arguments:
  -h, --help            show this help message and exit
  -p PHENO, --pheno PHENO
                        Specify phenotype/covariate file (default: None)
  -ph PHNAME, --pheno-name PHNAME
                        Only load the designated phenotype from the -p/--pheno
                        (default: None)
  -c COVAR, --covar COVAR
                        Specify covariate file (default: None)
  -cname COVARNAME, --cov-name COVARNAME
                        Covariates name to include into the model (default:
                        None)
  -o OUTPUTFILE, --out OUTPUTFILE
                        the name of the output file (default: None)
  -plot GWASPLOT, --plot GWASPLOT
                        Save Manhattan and Q-Q plots of your GWAS data
                        (default: yes)
  -chr CHR, --chr CHR   Include all variants on the given chromosome(s)
                        (default: ['1:22'])
  -maf MAF, --maf MAF   Exclude variants with nonmajor allele frequency lower
                        than a threshold (default: 1e-05)
  -pca PCA, --pca PCA   Set Number of PCA to include into the model (default:
                        10)
  -rsq {0.3,0.6,0.8}, --rsq {0.3,0.6,0.8}
                        Set Imputation Quality threshold (default: 0.3)
  -cpu CORES, --cpu CORES
                        Set maximum number of compute threads (default: 45)
  --adj {yes,no}, --adj {yes,no}
                        Flag to include sex, array and PCs{1:10} as covariates
                        (default: None)
  -m [{model1,model2,model3,model4}], --model [{model1,model2,model3,model4}]
                        See WIKI page on K:Drive (default: None)
  -v, --version         Shows the version number of the program
```




<br>

Users have to specify phenotype/covariate file and phenotype'name (mandatory). In addition, option -adj has to be set to *yes* or *no* evert time.

### Some examples


To run a gwas using all UKB samples, including sex+array+10PCs+age as covariates in BOLT, including only SNPs with MAF>=5% and using 60CPUs :
```
python3 ukb_gwas.py -p bmi.pheno -o test -m model1 -ph bmi -cname age -maf 0.05 -cpu 60 -adj yes
```

<br>


To run a Case-Control gwas on unrelated White British samples only, including sex+array+10PCs+age as covariates, using PLINK and  60 CPUs :

```
python3 ukb_gwas.py -p pheno.txt -o test -ph CAD -cname age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,array  -adj yes --model model4 -cpu 60
```

To tun GWAS on binary trait without including covariates into the model:


```
python3 ukb_gwas.py -p pheno.txt -o test -ph CAD --model model4 -adj no
```

<br>


## Output files


### BOLT-LMM

BOLT-LMM association statistics are output in a compressed tab-delimited file (with suffix '.BOLT.GWAS.gz') with the following fields, one line per SNP:

```
SNP:     rs number or ID string
CHR:     chromosome
BP:      physical (base pair) position
GENPOS:  genetic position either from bim file or interpolated from genetic map
ALLELE1: first allele in bim file (usually the minor allele), used as the effect allele
ALLELE0: second allele in bim file, used as the reference allele
A1FREQ:  frequency of first allele
F_MISS:  fraction of individuals with missing genotype at this SNP
BETA:    effect size from BOLT-LMM approximation to infinitesimal mixed model
SE:      standard error of effect size
P_BOLT_LMM_INF: infinitesimal mixed model association test p-value
P_BOLT_LMM:     non-infinitesimal mixed model association test p-value
logOR:          log odds ratio (only for binary gwas)
SE_logOR:       standard error of log odds ratio (only for binary gwas)
MAF:            Minor allele frequency in UKBB (all samples)
RSQ:            Imputation Info Score

```

### PLINK2

The Plink-output (with a case/control or quantitative phenotype) is a text file  with a header line, and then one line per variant with the following columns:

```
CHROM  : Chromosome code
POS    : Base-pair coordinate
ID     : Variant ID
REF    : Reference allele
ALT    : Alternate alleles
A1     : Counted allele in regression
TEST   : Test identifier (e.g., Additive...)
OBS_CT : Number of samples in the regression
BETA   : Regression coefficient (for A1 allele)
OR     : Odds ratio (for A1 allele)
SE     : Standard error of beta (or log-odds)
L95    : Bottom of symmetric approx. confidence interval
U95    : Top of symmetric approx. confidence interval
T_STAT : t-statistic for linear regression, Z-score for logistic
P      : Asymptotic p-value (or -log10(p)) for T/chisq-stat
MAF    : Minor allele frequency in UKBB (all samples)
RSQ    : Imputation Info Score

```



## Summary Statistics Plots


After a GWAS, the pipeline saves Manhattan and quantile-quantile (Q-Q) plots (without filtering for MAF or Imputation quality).



## Software version

That pipeline works with **Python 3**, PLINK2 (PLINK v2.00a2LM 64-bit Intel (2 Jan 2019)), BOLT-LMM v2.3.2 and R (>=3.2). If you don't have python 3 installed into your linux environment you can't use the pipeline. 


## Links

1. BOLT-LMM: https://data.broadinstitute.org/alkesgroup/BOLT-LMM/
2. PLINK2: https://www.cog-genomics.org/plink/2.0/
3. Python: https://www.python.org/
4. Anaconda Distribution :https://www.anaconda.com/distribution/







