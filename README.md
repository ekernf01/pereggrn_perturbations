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

- Files: each dataset is required to have a perturb-seq style experiment stored as an h5ad file `<path>/dataset_name/test.h5ad`. Data may optionally include a time-series experiment stored at `<path>/dataset_name/train.`. Here, `path` is the same path provided to `pereggrn_perturbations.set_data_path`, and `dataset_name` is the same string provided to `pereggrn_perturbations.load_perturbation`.
- Format (general): 
    - 🚧 These docs are under construction 🚧
- Format (perturbations): Any perturbation dataset must contain specific columns in .obs:
        - `perturbation`: The gene(s) perturbed. String. For multi-gene perturbations, this may contain comma-separated gene names.
        - `expression_level_after_perturbation` column: Expression levels after perturbations. Numeric. For multi-gene perturbations, this may contain comma-separated numbers stored in strings.
        - `perturbation_type` column: Types of perturbations (e.g., "overexpression", "knockout", "knockdown").
- Format (time-series): any time-series dataset must contain specific columns in .obs:
    - `timepoint` column: This column should contain numeric values representing the time points at which data was collected.
    - `cell_type` column: This column should contain string values indicating the cell types in the dataset.
