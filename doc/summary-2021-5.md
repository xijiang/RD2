# ReDiverse summary
<!-- Note: many descriptions are as function docs leading the header of a function -->

## Data Summary

### Pedigree

The Dutch and German population share a few ID.
They were merged first, and then merged with the Norwegian pedigree.

### Genotypes

#### Merge

There are 7 genotyping platforms:
Illumina cow 50k v1-3, 777k, and 3 platforms (10690, 10993, 11483) for Dutch population.
Illumina cow 50k-v3 was used as target platform. 
All other genotypes were imputed to this level.
Extra SNP are excluded.
Sex chromosome and SNP with unknown chromosomes were removed.

For Norwegian data, duplicate ID of lower density or genotyped earlier were removed.

Raw data were of various formats.
They were merged into platform specific plink and VCF files,
for ease of downstream manipulations.

#### Quality control

Visually inspect the plot of missing ratio and HWE of the genotypes.
SNP of 
- genotype missing ratio > 0.1
- MAF < 0.01
- P-value of HWE statistics <0.0001

ID of
- genotype missing ration >0.1

are removed.
Duplicate SNP (different name, same chomosome and base pair positon) were removed.

#### Imputation

Imputation was with `FImpute`, which utilized pedigree information.
Norwegian data were imputed independantly.
Dutch and German data were imputed together.

#### G and A calculation
These were done with Mario's `calc_grm` program.

### Phenotypes

To date, the Dutch population has 16 traits.
The German population has 18 traits.
The Norwegian population has 14 traits.

Four traits, milk, fat, production and SCS were analysed.
One or two traits may be added later.

The phenotypes were deregressed, and using their reliability as weights.

## Cross-validation setup
### Dutch data
  - exclude ID not in the G
    - 0 excluded
    - total 2445 ID
  - put parent/offspring pairs in the training set first
    - 556 ID are parents of ID in current dataset
    - training set initialized with 1819 ID
  - put daughters not in 1819 into validation set
    - 515 cows, which is close 500. So the validation set is done
  - put the rest into the training set
    - +111 ID -> 1930 ID for training

### German scenario
  - exclude ID not in the G matrix. 
    - 0 excluded.
    - total 716 ID
  - put mothers and daughters in the training set
    - 98 ID are dams of ID in the current dataset.
    - training set is updated with 204 initial ID
    - shuffle the rest $(716-204)= 512$ ID, to let the training and validation sets have 358 ID each.

### Norwegian scenario

Similar as above.
Only have had daughters' parents removed from the training set.

## Cross-validation results
Countries were fitted as a fixed effect for breeding value estimation of the validation set.
correlation between EBV and breeding values in raw data as a the CV indicator.

Note below, the 2nd half of a table are results with no weight.

### The milk trait
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| D | 0.5 | 0.0017 | 0.16 | 0.49 | 0.51 | 0.11 | 0.5 |
| G | 0.1 | 0.31 | 0.032 | 0.27 | 0.11 | 0.29 | 0.27 |
| N | 0.15 | 0.049 | 0.65 | 0.14 | 0.65 | 0.65 | 0.65 |
| |
| D | 0.51 | 0.00076 | 0.056 | 0.5 | 0.5 | 0.065 | 0.49 |
| G | 0.12 | 0.31 | 0.11 | 0.31 | 0.14 | 0.31 | 0.31 |
| N | 0.17 | 0.053 | 0.65 | 0.17 | 0.66 | 0.66 | 0.66 |

### The protein trait
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| D | 0.53 | 0.083 | 0.16 | 0.53 | 0.53 | 0.16 | 0.53 |
| G | 0.21 | 0.44 | -0.012 | 0.35 | 0.21 | 0.34 | 0.34 |
| N | 0.099 | 0.053 | 0.67 | 0.1 | 0.67 | 0.67 | 0.67 |
| |
| D | 0.55 | 0.088 | 0.042 | 0.55 | 0.5 | 0.092 | 0.5 |
| G | 0.22 | 0.44 | 0.032 | 0.39 | 0.19 | 0.27 | 0.34 |
| N | 0.14 | 0.056 | 0.66 | 0.15 | 0.66 | 0.66 | 0.66 |

### The fat trait
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| D | 0.57 | 0.35 | 0.056 | 0.58 | 0.57 | 0.31 | 0.58 |
| G | 0.3 | 0.47 | -0.05 | 0.45 | 0.29 | 0.38 | 0.42 |
| N | 0.0019 | 0.11 | 0.69 | 0.052 | 0.68 | 0.69 | 0.68 |
| |
| D | 0.6 | 0.35 | 0.0053 | 0.61 | 0.56 | 0.26 | 0.57 |
| G | 0.32 | 0.48 | -0.06 | 0.48 | 0.26 | 0.37 | 0.42 |
| N | 0.013 | 0.11 | 0.7 | 0.077 | 0.71 | 0.71 | 0.71 |

### The SCS trait
| V\\T | D | G | N | D + G | D + N | G + N | D + G + N |
| -- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| D | 0.56 | 0.41 | 0.11 | 0.57 | 0.57 | 0.38 | 0.57 |
| G | 0.55 | 0.71 | 0.012 | 0.67 | 0.53 | 0.67 | 0.66 |
| N | -0.0022 | -0.0097 | 0.64 | 0.012 | 0.64 | 0.64 | 0.63 |
| |
| D | 0.6 | 0.41 | -0.092 | 0.61 | 0.32 | 0.0015 | 0.35 |
| G | 0.53 | 0.71 | 0.016 | 0.67 | 0.31 | 0.4 | 0.49 |
| N | -0.033 | -0.017 | 0.38 | -0.035 | 0.38 | 0.38 | 0.38 |
