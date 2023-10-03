# Genomic prediction in small population breeds with a joint reference population

RD WP6: Genomic selection strategies

- XY MC TM SW

## Abstract

## Introduction



### Objective

Validate genomic prediction based on a joint reference population of ERDB

## Materials and methods
### Genotypes
Description of the source data.
Imputation.			  <!-- used beagle, I organized the data -->

Genotype organization steps:
1. map preparation
2. merge pedigrees from the Netherlands and Germany
3. organize genotypes, take the autosome subset.
4. QC of the genotypes, and remove duplicates
5. Imputation, Norwegian data using Beagle, Dutch and German data using FImpute


- `FImpute` was used, which incorporated pedigree information. 

On Thursday afternoon I proposed to Xijiang that I could try to impute the MRY genotypes using FImpute, making explicit use of the pedigree information. Xijiang then made a genotype and pedigree file available including all Dutch and German genotyped animals and their ancestors. We have now combined the Dutch and German data (for the imputation), because these populations are very closely related.

In short: the imputed data from FImpute looks very good. The number of genotyped parent-offspring pairs with too many opposing homozygote SNPs was quite low, and for most of those, a cow involved had > 1% of her (genotyped) genotypes changed during imputation. I removed a total of 15 cows involved in these conflicts, and thereafter repeated the imputation. (I added some more details in the attached Word document.) Now there are no conflicts left (the maximum % of opposing homozygote SNPs between parents and offspring is now 0.27%). See attached image that shows the genomic versus pedigree relationships in this data. This looks very clean to me, and the earlier inflated genomic inbreeding coefficients of the Dutch MRY bulls disappeared.

So, I propose that we proceed with this imputed data. Using this list of 3321 genotyped/imputed Dutch and German animals, over the next few days I will do the deregressions separately for the set of genotyped Dutch and German animals, respectively.

@Xijiang: I uploaded the imputed data (MRY_genotypes_imp_final_calc_grm.txt.gz), and some additional files (including a shell script with all the steps taken), to the following directory on your server: /home/guest/two/FImputed . The IDs are the numerical IDs as in the file “all.ped”. From the file you can subset the genotypes for the SNPs on the list as agreed earlier. Note that the file contains only 1 space between ID and genotypes, meaning that with longer IDs the genotype for the first SNP is given on the next position. Use the following commands to visualize this:


G-matrix calculation. <!-- used Mario's program -->

- Mario's `calc_grm` was used.

A-matrix

### Phenotypes

#### DRP (deregressed EBV)
DRP were computed using so-called matrix deregression (Jairath et al., 1998; Calus et al., 2016). 
During the deregression, the computed mean EBV was subtracted from all EBV, such that the subsequent GEBV, computed based on the DRP, should have a value close to zero.
Before computing the DRP, the EBV in each country were divided by the genetic standard deviation (as shown in Table 1), such that all DRP across the three countries have a genetic standard deviation of 1, and thus are on the same scale. 
As a result, given that now $h^2 = \frac{\sigma_a^2}{\sigma_a^2+\sigma_e^2}=\frac{1}{1+\sigma_e^2}$, the (rescaled) residual variance for each country becomes: $\sigma_e^2 = \frac{1-h^2}{h^2}$.

#### ERC (ERC (Effective Record Contributions = weights of the DRP)

ERC were obtained from the reliabilities (REL; assuming a 0-1 scale) of the EBV, using: $\mathrm{ERC}=\lambda\frac{\mathrm{REL}}{1-\mathrm{REL}}$, with $\lambda=\frac{\sigma_e^2}{\sigma_a^2}=\frac{1-h^2}{h^2}$, where $h^2$ is the heritability for a specific country scale (as shown in Table 1) that are used together with the ERC values in the computation of the DRP.

Note that in any subsequent analyses based on the DRP, for an individual animal $i$ its ERC value is used to scale the residual variance used for its own information (i.e. its DRP), as: $\sigma_{e_i}^2 = \frac{\sigma_e^2}{\mathrm{ERC}_i}$. 
Due to the standardization of the original EBV, $\sigma_e^2 = \frac{1-h^2}{h^2}$. 
To standardize the $\sigma_e^2$ values across the three countries, within each country both $\sigma_e^2$ and the ERC values were multiplied by $\frac{h^2}{1-h^2}$. 
Note that these adjusted ERC are thus: $\mathrm{ERC} = \frac{\mathrm{REL}}{1-\mathrm{REL}}$. 
After this final adjustment, for each country applies that $\sigma_a^2=1$, $\sigma_e^2=1$, $h^2=\frac{1}{2}$, and the original $h^2$ (as in Table 1) applies to an individual with an ERC of $\frac{h^2}{1-h^2}$.

- A summary of data

| Country | Total ID | Males | Females |
| -- | --: | --: | --: |
| Dutch | 2445 | 2073 | 372 |
| Germany | 716 | 0 | 716 |
| Norway | 7159 | 6527 | 632 |

- Heritabilities of the EBV in each country

| Trait | Netherlands | Germany | Norway |
| -- | --: | --: | --: |
| Milk yield | 0.51 | 0.39 | 0.277 |
| Fat yield | 0.52 | 0.30 | 0.213 |
| Protein yield | 0.44 | 0.30 | 0.235 |
| Somatic cell score | 0.37 | 0.23 | 0.136 |
| Interval calving to first insemination | 0.176 | 0.039 | 0.142 |

- Genetic standard deviations of the EBV in each country

| Trait | Netherlands | Germany | Norway |
| -- | --: | --: | --: |
| Milk yield | 687 | 690 | 12 |
| Fat yield | 28.0 | 25.1 | 12 |
| Protein yield | 19.0 | 19.7 | 12 |
| Somatic cell score | 4 | 12 | 12 |
| Interval calving to first insemination | 4 | 12 | 4 |

#### Cross-validation setup
General strategy for crossvalidation

We use about 500 cows for crossvalidation. The cows used for crossvalidation should not be daughters of bulls that are in the training data, i.e. we try to generate genetic distance between validation cows and the training population. The remaining data is used for training in two way: within country and all data including foreign data. We compare the correlation between the validation cows’ DYDs and their predicted EBV. Specifically we compare how much this correlation increases by including the foreign data.
Special for Norwegian data

There are only ~600 cows. If the daughters of training bulls are excuded this results perhaps in <500 cows. If this is the case, perhaps some bulls with many daughters in the validation set could be excluded from the training set.
Special for German data

There are not many cows and no bulls. We try to go for 50% of the cows in the training and 50% in the validation group. Daughters of Dutch bulls that are in the training go into the training set. Hopefully we find a reasonable balance?
Special for the Dutch data

The general strategy above seems to fit best for the Dutch data. If there are problems in finding enough validation animals, we follow the same strategy as for the Norwegian data.

Xijiang: can you find out whether the above strategies for finding validation and training animals works out OK? If you have questions let me know.

After 

- filtering ID with no breeding values:
- remove ID with strange chips.

| Set | Dutch | Germany | Norway |
| -- | --: | --: | --: |
| Training | 1925 | 358 | 6021 |
| Validation | 515 | 358 | 587 |
| Sum | 2440 | 716 | 6608 |

### Models

Assuming genetic correlations are 1 between countries.

Single trait GBLUP model was used:

$$\mathbf{y} = \mathrm{country} + \mathrm{animal} + \mathbf{e}$$

where,
$\mathbf{y}$ are the deregressed EBV, weighted by EDC.
Country was included as a fixed effect.

## Results
### Milk yield
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| WT |
| D | 0.5 | 0.0017 | 0.16 | 0.49 | 0.51 | 0.1 | 0.5 |
| G | 0.1 | 0.31 | 0.032 | 0.27 | 0.11 | 0.3 | 0.27 |
| N | 0.15 | 0.049 | 0.65 | 0.15 | 0.65 | 0.65 | 0.65 |
| NW|
| D | 0.51 | 0.00076 | 0.056 | 0.5 | 0.5 | 0.062 | 0.5 |
| G | 0.12 | 0.31 | 0.11 | 0.31 | 0.14 | 0.31 | 0.32 |
| N | 0.17 | 0.053 | 0.65 | 0.17 | 0.66 | 0.66 | 0.66 |

### Protein yield
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| WT |
| D | 0.53 | 0.083 | 0.16 | 0.53 | 0.54 | 0.16 | 0.53 |
| G | 0.21 | 0.44 | -0.012 | 0.35 | 0.21 | 0.36 | 0.34 |
| N | 0.099 | 0.053 | 0.67 | 0.1 | 0.67 | 0.67 | 0.67 |
| NW|
| D | 0.55 | 0.088 | 0.042 | 0.55 | 0.51 | 0.095 | 0.51 |
| G | 0.22 | 0.44 | 0.032 | 0.39 | 0.2 | 0.29 | 0.35 |
| N | 0.14 | 0.056 | 0.66 | 0.15 | 0.66 | 0.66 | 0.66 |

### Fat yield
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| WT |
| D | 0.57 | 0.35 | 0.056 | 0.58 | 0.57 | 0.32 | 0.58 |
| G | 0.3 | 0.47 | -0.05 | 0.45 | 0.29 | 0.4 | 0.42 |
| N | 0.0019 | 0.11 | 0.69 | 0.053 | 0.68 | 0.69 | 0.68 |
| D | 0.6 | 0.35 | 0.0053 | 0.61 | 0.57 | 0.28 | 0.58 |
| G | 0.32 | 0.48 | -0.06 | 0.48 | 0.27 | 0.38 | 0.43 |
| N | 0.013 | 0.11 | 0.7 | 0.079 | 0.71 | 0.71 | 0.71 |

### Somatic cell score
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| WT |
| D | 0.56 | 0.41 | 0.11 | 0.57 | 0.57 | 0.38 | 0.57 |
| G | 0.55 | 0.71 | 0.012 | 0.67 | 0.53 | 0.68 | 0.66 |
| N | -0.0022 | -0.0097 | 0.64 | 0.01 | 0.64 | 0.64 | 0.64 |
| NW|
| D | 0.6 | 0.41 | -0.092 | 0.61 | 0.35 | 0.014 | 0.38 |
| G | 0.53 | 0.71 | 0.016 | 0.67 | 0.33 | 0.43 | 0.52 |
| N | -0.033 | -0.017 | 0.38 | -0.034 | 0.38 | 0.38 | 0.38 |

### Calving to 1st insemination
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| WT |
| D | 0.57 | 0.24 | -0.096 | 0.54 | 0.55 | 0.23 | 0.52 |
| G | 0.071 | 0.13 | -0.035 | 0.1 | 0.063 | 0.12 | 0.096 |
| N | 0.12 | -0.00018 | 0.51 | 0.094 | 0.51 | 0.51 | 0.51 |
| NW|
| D | 0.59 | 0.25 | -0.022 | 0.55 | 0.41 | 0.15 | 0.4 |
| G | 0.11 | 0.12 | -0.036 | 0.12 | 0.045 | 0.092 | 0.11 |
| N | 0.16 | 0.0063 | 0.32 | 0.12 | 0.32 | 0.32 | 0.32 |

## Discussion

### Summary

- Accuracies of across-breed predictions suggest that
  - RDN benefits from MRY, and vice versa
  - NR: no predictive ablity for RDN and MRY, and vice versa
  - all in line with the $M_e$ results.
- Accuracies never improved when adding another breed:
  - Hardly changed for MRY and NR
  - Decreased for RDN
  - Not in line with the $M_e$ results
  
### Why no benefit for combining breeds

- SNP density (50k and less?)
  - should be enough for MRY/RDN
  - Higher density could help NR versus MRY/RDN
  - but the breeds may simply be too divergent ($M_e$)
- Genetic correlations between countries <1?
- Base of different populations not (sufficiently) aligned?
- Small breed RDN overwhelmed by larger MRY/NR
- Fitting a multi-trait model may overcome all of those.

## Conclusion

- MRY and GRW show high relatedness
- More benefit expected for smaller MRY and GRW breed than NR

## References
