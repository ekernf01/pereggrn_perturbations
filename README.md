This is a companion package to our [collection of perturbation data](https://github.com/ekernf01/perturbation_data). It finds, loads, and validates the format of perturbation datasets.

### Installation
    
    pip install git+https://github.com/ekernf01/pereggrn_perturbations

### Usage

Using this module involves setting the data path, loading metadata, and validating the datasets.

```
import pereggrn_perturbations
# Set this to point to the "perturbations" folder in the perturbation data collection. 
pereggrn_perturbations.set_data_path("path/to/perturbation_data/perturbations")
# What datasets are available?
pereggrn_perturbations.load_perturbation_metadata()
# Grab a dataset
nakatake_et_al = pereggrn_perturbations.load_perturbation("nakatake") 
# Validate the format
pereggrn_perturbations.check_perturbation_dataset(ad = nakatake_et_al)
```

### Data format

- Files
    - each dataset is required to have a perturb-seq style experiment stored as an h5ad file `<path>/dataset_name/test.h5ad`.
    - Data may optionally include a time-series experiment stored at `<path>/dataset_name/train.h5ad`.
    - Here, `path` is the same path provided to `pereggrn_perturbations.set_data_path`, and `dataset_name` is the same string provided to `pereggrn_perturbations.load_perturbation`.
- Format (general):
    - In `.var`, there must be a column `highly_variable_rank` ranking the genes in terms of variance or overdispersion. The most interesting genes should be ranked lowest. Positive integer. This must match exactly between `train.h5ad` and `test.h5ad` if `train.h5ad` is present. This allows `pereggrn` users to experiment with different amounts of highly variable genes. The gene names themselves must also match exactly between `train.h5ad` and `test.h5ad`.
    - In `.obs`, there must be columns:
        - `is_control` (boolean): Whether the sample is a control or is perturbed.
        - `perturbation`: (string) The gene(s) perturbed. For multi-gene perturbations, this may contain comma-separated gene names.
        - `expression_level_after_perturbation` (numeric): Expression levels after perturbations. For multi-gene perturbations, this may contain comma-separated numbers stored in strings.
        - `perturbation_type` (string): The type of perturbation (e.g., "overexpression", "knockout", "knockdown"). This is only checked for `test.h5ad`.
    - In `.raw`, there must be raw expression data. 
    - In `.uns`, there must be iterables "perturbed_and_measured_genes" and "perturbed_but_not_measured_genes".
- Format (`train.h5ad` only):
    - `timepoint` column: This column should contain numeric values representing the time points at which data was collected.
    - `cell_type` column: This column should contain string values indicating the cell types in the dataset.
- Contents:
    - Each anndata object must be at least one control sample.
    - `test.h5ad` must include non-control samples.
    - The number of perturbations must match the number of reported expression levels.
    - The reported post-perturbation expression must closely match the data in the matrix, except in knockout samples.
    - All measured data must actually be measured.
    - All perturbed data must actually be perturbed.
    - Expression in `.X` must be normalized and natural-log-transformed.
    - Contents of `.raw` should be raw in the sense of having non-normalized integer counts (needed by e.g. Dictys) and having all genes, not a subset (needed by GeneFormer). This is not checked yet.

