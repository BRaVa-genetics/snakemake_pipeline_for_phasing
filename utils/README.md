## Assessing imputed homozygotes after phasing
The scripts in this directory aim to quantify how many missing genotypes in the raw/input VCF are imputed as homozygotes (especially non-reference) during phasing with Shapeit5. Imputing missing genotypes as homozygotes can introduce errors, so we flag these cases to assess the reliability of the phasing and imputation process.

Imputation of missing genotypes is an important step during phasing, but may introduce errors if the goal is to survey rare knockouts. Here we systematically compare raw and phased data to estimate the rate of potentially incorrect homozygote imputations and improve confidence in downstream analyses that rely on phased genotypes. Nonetheless, getting an accurate estimate of the fraction of wrongly imputed genotypes is difficult, as that depends on several factors, including overapping samples between target/reference, available variants, multiallelic sites, or QC guidelines.

### Requirements
* Phased and raw (in VCF?) genotypes, and the output from `called_chet` (for a given variant set).
* Python including the module `pysam` for reading VCF/BCF files.

### Workflow
#### 1. Prepare List of Genotypes to Assess
For each chromosome, extract a list of sample-variant pairs where the genotype is homozygous alternate (hom:alt) in the phased data, given a set of annotated variants, e.g. pLoF. The input here should be the output of `called_chets`, the software we used previously to generate all biallelic genotypes.

#### 2. Compare Raw and Phased Genotypes
For each sample-variant pair, compare the genotype in the raw (unphased) VCF and the phased BCF (or VCF). This is performed with `debug_homozygotes.py` which reports cases for further inspection and prints summaries with number of:
1. cases where both raw and phased genotypes are homozygous (mathed).
2. cases where the raw genotype is missing but the phased genotype is homozygous (imputed, potentially wrong). 
3. other cases, like when a missing getotype is imputed as hom:ref (which could also be wrong).

#### 3. Summarize Results

Aggregate results across chromosomes to assess the overall extent of imputation. This could be done with something like this:
```
PREFIX="assess_homz.report.gnh_39k"
echo "Working for $PREFIX"
grep Checked $PREFIX.{1..22}.* | cut -d: -f3 | awk '{s=s+$1} END {print s, "("NR" files)"}'
a=$(grep Correct $PREFIX.{1..22}.* | cut -d: -f3 | awk '{s=s+$1} END {print s}')
b=$(grep -E 'Flagged.*ALT' $PREFIX.{1..22}.* | cut -d: -f4 | awk '{s=s+$1} END {print s}')
c=$(grep Other $PREFIX.{1..22}.* | cut -d: -f3 | awk '{s=s+$1} END {print s}')
echo -e "Matched homozygotes: $a\nPotentially wrong: $b + $c\nTotal: $((a+b+c))"
```

Note that not all imputed genotypes will be wrong - hopefully a small set - and those imputed as hom:ref will be more alarming (i.e. "false positives") than those imputed as hom:alt ("false negatives"). Also note that "other" (=`$c`) will catch cases with multiple homozygotes or cases where the variant is not observed at all, so please inspect those manually (they're printed, if any).

For reference, I report the numbers in the two versions of G&H in the next table, while also giving the uppper bound of false positives, which is the fraction of "missing -> 1|1" out of all homozygotes. 

### Case study: Genes & Health
|  | v1 (39k) | v2 (49k) |
|--:|:-----:|:-----:|
Total | 13821 | 17126 |
Matched homozygotes| 13465 | 17056 |
Potentially wrong| 338 + 18 | 0 + 70 |
Upper bound for FP | 2.58% | 0.4% |

### Case study: 100,000 Genomes Project
This phenomenemon is also observed in the 100kGP, where Shi et al (2024) followed a totally different approach for phasing (though they worked with Shapeit v4, not v5). In particular, and after adapting the pipeline to work on 1,371 chunks instead of whole-chromosome VCFs, the above quantification resulted in the following numbers:

|  | EUR | SAS |
|--:|:-----:|:-----:|
Total | 8978 | 4837 |
Matched homozygotes| 8915 | 4822 |
Potentially wrong| 0 + 63 | 2 + 13 |
Upper bound for FP | 0.8% | 0.4% |
