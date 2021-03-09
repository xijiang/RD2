# Data organization

Data were stored in four files of JLD2 (compressed) format.
- cv-setup.jld
- drp-training.jld
- ebv-all.jld
- GAID.jld


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
  - =>did.txt

### German scenario
  - exclude ID not in the G matrix. 
    - 0 excluded.
    - total 716 ID
  - put mothers and daughters in the training set
    - 98 ID are dams of ID in the current dataset.
    - training set is updated with 204 initial ID
    - shuffle the rest (716-204)= 512 ID, to let the training and validation sets have 358 ID each.
  - => gid.txt

### Norwegian scenario
Similar as above.  Only have had daughters' parents removed from the training set.

## DRP training data
> From these training animals, I removed 39 Norwegian, 6 Dutch, and 10 German animals, that had a reliability < 15% for their EBV.
 This is just to remove those few animals with very little (or even no) information. A cut-off of 15% is still low, but taking a higher value quickly removes more German animals.

> I then deregressed the EBVs within the remaining training animals, within each country. Attached are the resulting files (same format as before). (from Mario)

## Validation set / EBV
>Note that the files:

- Include all EBVs (training + validation).
- The EBV are on the original country scale, whereas the DRP are on a standardized scale with mean 0 and genetic SD of 1. For computation of the accuracy, the scale of the EBV, however, does not matter.
- Have (per country) the same order of columns as the DRP files. Where the DRP files have columns with ERC, these EBV files have columns with reliabilities (all on a 0-100 scale).

> Maybe for the validation animals you want to use a minimum reliability? For the training, their reliability should be 15% or greater. That may be a sensible cut-off for the validation animals as well? (Mario)
