This is a companion package to our [collection of perturbation data](https://github.com/ekernf01/perturbation_data).

```python
import pereggrn_perturbations
# Set this to point to the "perturbations" folder in the perturbation data collection. 
pereggrn_perturbations.set_data_path("path/to/perturbation_data/perturbations")
# What datasets are available?
pereggrn_perturbations.load_perturbation_metadata()
# Grab one
nakatake_et_al = pereggrn_perturbations.load_perturbation("nakatake") 
```

## Perturbation Datasets Checker

This module is designed to validate the format and integrity of perturbation datasets. It performs the following tasks:

- Path Management: Set and retrieve paths for perturbation data using environment variables.
- Data Loading: Load perturbation datasets and their metadata from specified directories.
- Validation: Check the integrity and format of perturbation data using custom criteria, including verifying the presence and format of required files, ensuring consistency between datasets, and confirming the correctness of dataset contents.

The module utilizes the pandas, scanpy, anndata, numpy, and os libraries to perform these tasks. It is intended for use with specific data structures commonly found in biological datasets, particularly those involving gene perturbations.

### Installation

To use this module, ensure that you have the required libraries installed:
    
    pip install pandas scanpy anndata numpy

### Usage
Using this module involves setting the data path, loading metadata, and validating the datasets as described below:

- `set_data_path(path: str)` sets the path for perturbation data and verifies the existence of required files.
#### Example usage
    set_data_path('/path/to/your/data')

- `get_data_path()` retrieves the path to the perturbation data from an environment variable.
#### Example usage
    data_path = get_data_path()

- `load_perturbation_metadata()` loads metadata about perturbations from a CSV file at the designated path returned from `get_data_path()`.
#### Example usage
    metadata = load_perturbation_metadata()

- `load_perturbation(dataset_name: str, is_timeseries: bool = False)` loads a specific perturbation dataset from an AnnData file.
    - The `is_timeseries` parameter is a boolean flag that indicates whether the dataset being validated includes time-series data. This parameter affects the types of checks and validations performed on the dataset.
    - When `is_timeseries` is `False`:
        - The dataset is treated as a static dataset, meaning it contains data from a single time point or independent of time.
        - Validation checks do not include time-specific metadata.
        - The dataset is expected to contain perturbation data, but not necessarily any information about the timing of observations or changes over time.
    - When `is_timeseries` is `True`:
        - The dataset is treated as a time-series dataset, meaning it includes data across multiple time points.
        - Additional validation checks are performed to ensure the presence and correctness of time-specific metadata.
        - The dataset must contain specific columns in the observation metadata:
            - `timepoint` column: This column should contain numeric values representing the time points at which data was collected.
            - `cell_type` column: This column should contain string values indicating the cell types in the dataset.
#### Example usage
    dataset = load_perturbation('nakatake')

- `check_perturbation_dataset(dataset_name: str = None, ad: anndata.AnnData = None, is_timeseries = False, do_full = False, is_perturbation = True)` validates the format and integrity of a loaded perturbation dataset.
    - `dataset_name`
        - Default (None):
            - If `dataset_name` is not provided, the function expects `ad` (AnnData object) to be provided instead.
            - If both `dataset_name` and `ad` are `None`, the function raises a `ValueError`.
        - When Provided:
            - The function loads the dataset with the given `dataset_name` and performs validation checks.
            - The dataset is loaded using the `load_perturbation` function.
    - `ad`
        - Default (None):
            - If `ad` is not provided, the function expects `dataset_name` to be provided instead.
            - If both `dataset_name` and `ad` are `None`, the function raises a `ValueError`.
        - When Provided:
            - The function uses the provided AnnData object `ad` for validation checks without loading a dataset from a file.
    - `is_timeseries`
        - Default (False):
            - The dataset is treated as a static dataset, meaning it contains data from a single time point or independent of time.
            - Time-series-specific checks (such as verifying the presence of `timepoint` and `cell_type columns`) are not performed.
        - When True:
            - The dataset is treated as a time-series dataset, meaning it includes data across multiple time points.
            - Additional checks ensure the presence and correctness of time-specific metadata:
                - `timepoint` column: Numeric values representing the time points at which data was collected.
                - `cell_type` column: String values indicating the cell types in the dataset.
    - `do_full`
        - Default (False):
            - The function performs a basic validation, which includes necessary checks but skips more expensive or comprehensive checks.
            - This mode is faster but may miss some detailed issues.
        - When `True`:
            - The function performs a full validation, including all comprehensive checks.
            - This mode is thorough but may be slower due to the additional checks performed.
    - `is_perturbation`
        Default (True):
        - The function assumes the dataset is a perturbation dataset, which includes additional metadata related to perturbations.
        - The validation checks ensure the presence of perturbation-specific metadata, such as:
            - `perturbation column`: Describes the perturbations applied.
            - `expression_level_after_perturbation` column: Expression levels after perturbations.
            - `perturbation_type` column: Types of perturbations (e.g., overexpression, knockout, knockdown).
        - When `False`:
            - The function treats the dataset as a non-perturbation dataset, which does not include perturbation-specific metadata.
            - Perturbation-specific checks are skipped.
#### Example usage
    is_valid = check_perturbation_dataset(ad=dataset)
    print("Dataset validation result:", is_valid)