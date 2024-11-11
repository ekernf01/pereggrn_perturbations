import os            # Import the os module to interact with the operating system
import pandas as pd  # Import the pandas library for data manipulation and analysis
import scanpy as sc  # Import the scanpy library for analyzing single-cell sequencing data
import anndata       # Import the anndata library to handle annotated data matrices in biology
import numpy as np   # Import the numpy library for numerical operations

def get_data_path():
    """
    Retrieve the path to the perturbation data from an environment variable.

    Returns:
        str: The path stored in the 'PERTURBATION_PATH' environment variable.
    """
    return os.environ['PERTURBATION_PATH'] # Return the path stored in the 'PERTURBATION_PATH' environment variable

def set_data_path(path: str):
    """
    Set the path for perturbation data in an environment variable and check for the existence of required files.

    Args:
        path (str): Path to the directory containing perturbation datasets.

    Raises:
        FileNotFoundError: If the required files are not found at the specified path.
    """
    # Check if the 'perturbations.csv' file exists at the specified path
    if not os.path.isfile(os.path.join(path, "perturbations.csv")): 
        raise FileNotFoundError("perturbations.csv should be a file in the folder whose name is provided to this function.")
        # Check if the 'test.h5ad' file exists within the 'nakatake' sub-directory of the specified path
    if not os.path.isfile(os.path.join(path, "nakatake", "test.h5ad")):
        raise FileNotFoundError("There should be an AnnData file at <your_input>/nakatake/test.h5ad (and others like it for the other datasets).")
    os.environ['PERTURBATION_PATH'] = path # Set the 'PERTURBATION_PATH' environment variable to the validated path
    return

def load_perturbation_metadata():
    """
    Load metadata about perturbations from a CSV file at the designated path.

    Returns:
        pandas.DataFrame: DataFrame containing perturbation metadata.

    Raises:
        KeyError: If the required environment variable is not set.
    """
    # Attempt to load and return the CSV file containing perturbation metadata
    try:
        return pd.read_csv(os.path.join(get_data_path(), "perturbations.csv"))
    # Raise a KeyError if the required environment variable is not set
    except KeyError as e:
        raise(KeyError("Before using the data you must call set_data_path('path/to/collection') to point to the perturbation data collection."))

def load_perturbation(dataset_name: str, is_timeseries: bool = False, is_screen: bool = False):
    """
    Load a perturbation dataset from an AnnData file.

    Args:
        dataset_name (str): Name of the dataset, taken from the metadata row names.
        is_timeseries (bool, optional): If True, loads the training data without perturbation (usually a timecourse). Defaults to False. Ignored if is_screen is True.
        is_screen (bool, optional): If True, loads CRISPR screen or literature review data.

    Returns:
        anndata.AnnData: Perturbation data in a uniform format.

    Raises:
        KeyError: If the dataset cannot be found at the specified path.
    """
    try:
        if is_screen:
            return  pd.read_csv(os.path.join(get_data_path(), dataset_name, "screen.csv"))
        else:
            t = "train" if is_timeseries else "test" # Determine the dataset type based on whether it is part of a time series
            return sc.read_h5ad(os.path.join(get_data_path(), dataset_name, f"{t}.h5ad"))
    except KeyError as e:
        raise(FileNotFoundError(f"Could not find data at this path. Before using the data you must call set_data_path('path/to/collection') to point to the perturbation data collection. Error was {e}."))
    except FileNotFoundError as e:
        raise(FileNotFoundError(f"Could not find data at this path. Before using the data you must call set_data_path('path/to/collection') to point to the perturbation data collection. Error was {e}."))

# Function to validate the format and integrity of a loaded perturbation dataset
def check_perturbation_dataset(dataset_name: str = None, ad: anndata.AnnData = None, is_timeseries = False, do_full = False, is_perturbation = True):
    """
    Validate the format and integrity of a loaded perturbation dataset.

    Args:
        dataset_name (str, optional): Name of the dataset. Provide exactly one of `dataset_name` or `ad`. This specifies the dataset to be checked.
        ad (anndata.AnnData, optional): AnnData object containing perturbation data. Provide exactly one of `dataset_name` or `ad`. This specifies the dataset to be checked.
        is_timeseries (bool, optional): If True, performs checks specific to time-series data. Defaults to False.
        do_full (bool, optional): If True, performs a full validation, including more time-consuming checks. Defaults to False. This mode is faster but may miss some detailed issues.
        is_perturbation (bool, optional): If True, treats the data as a perturbation dataset with additional metadata. Defaults to True.
    Returns:
        bool: True if the input data are correctly formatted.

    Raises:
        ValueError: If both or neither `dataset_name` and `ad` are provided.
        AssertionError: If various expected conditions are not met.
    """
    # Ensure that exactly one of 'ad' or 'dataset_name' is provided
    if ad is None and dataset_name is None:
        raise ValueError("Provide exactly one of ad and dataset_name")
    if not ad is None and not dataset_name is None:
        raise ValueError("Provide exactly one of ad and dataset_name")
    if ad is None and dataset_name is not None: 
        # A tiny bit of recursion helps us check a dataset with separate train and test folds. 
        # The base-case: AnnData input.
        try:
            assert all(load_perturbation(dataset_name, is_timeseries = True).var_names == 
                       load_perturbation(dataset_name, is_timeseries = False).var_names), "Gene names do not match between train and test data."
            try:
                print("Checking FACS-based CRISPR screen data...", flush = True)
                screen = load_perturbation(dataset_name, is_screen = True)
                assert screen.columns[0] == "perturbation", "If present, FACS-based CRISPR screen data must have a 'perturbation' column as the first column."
                assert len(set(screen["perturbation"]).intersection(load_perturbation(dataset_name, is_timeseries = False).var_names)) > 10, "If present, CRISPR screen data must have at least 10 genes in common with the test data feature set."
            except FileNotFoundError:
                print("...no FACS-based CRISPR screen data found, which is fine.")
                screen = None
            check_perturbation_dataset(ad=load_perturbation(dataset_name, is_timeseries = True), is_timeseries = True, is_perturbation = False)
            check_perturbation_dataset(ad=load_perturbation(dataset_name, is_timeseries = False), is_timeseries = True, is_perturbation = True)
        except FileNotFoundError:
            # If only test data is found, check its validity. It's allowed to lack timeseries metadata. 
            check_perturbation_dataset(ad=load_perturbation(dataset_name, is_timeseries = False), is_timeseries = False, is_perturbation = True)
        return
    
    # We will later select a variable number of genes based on this ranking. 
    print("Checking gene metadata...", flush = True)
    # Check gene metadata for ranking and ensure all necessary columns are present.
    assert "highly_variable_rank" in set(ad.var.columns), "Genes must be ranked in .var['highly_variable_rank']"
    assert all(~ad.var["highly_variable_rank"].isnull()), "Gene rankings should not be missing for any genes."
    assert all(ad.var["highly_variable_rank"]>=0), "Gene rankings must be positive integers"

    # Time   
    if is_timeseries: 
        print("Checking celltype and timepoint labels...", flush = True)
        assert "timepoint" in set(ad.obs.columns), "Time-series data must have a numeric 'timepoint' column"
        assert "cell_type" in set(ad.obs.columns), "Time-series data must have a string 'cell_type' column"

    # Names of genes perturbed
    print("Checking perturbation labels...", flush = True)
    # Validate perturbation labels and their corresponding expression levels
    assert "perturbation" in set(ad.obs.columns), "No 'perturbation' column"
    
    # Level of those genes after perturbation
    assert "expression_level_after_perturbation" in set(ad.obs.columns), "No 'expression_level_after_perturbation' column"
    iter = 0
    for i in ad.obs.index: 
        iter = iter + 1
        p = ad.obs.loc[i, "perturbation"] 
        elap = ad.obs.loc[i, "expression_level_after_perturbation"] 
        n_levels = len(str(elap).split(",")) 
        n_perts =  len(str(p   ).split(",")) 
        assert n_levels==n_perts, f"Too many or too few expression_level_after_perturbation entries in sample {i}: {p} has {n_perts} and {elap} has {n_levels}"
        # This check takes a long time, so by default we will apply it only to the first 1000 samples.
        if ~ad.obs.loc[i, "is_control"] and (ad.obs.loc[i, "perturbation_type"] != "knockout") and (do_full or iter < 1000):
            for x,g in zip(str(elap).split(","), str(p   ).split(",")):
                if g in ad.var_names:
                    assert np.isnan(float(x)) or np.abs(float(x) - float(ad[i,g].X[0,0])) < 0.0001, f"For observation {i}, post-perturbation expression is given in .obs as {x} but the value in .X is {ad[i,g].X[0,0]}."

    print("Checking control labels...", flush = True)
    assert "is_control"   in set(ad.obs.columns), "No 'is_control' column" 
    assert bool==ad.obs["is_control"].dtype, "non-boolean 'is_control' column"
    assert       ad.obs["is_control"].any(), "no controls found" 

    # Validate perturbation types
    if is_perturbation:
        # Overexpression / knockout / knockdown
        assert "perturbation_type" in set(ad.obs.columns), "No 'perturbation_type' column"    
        assert all(
            [pt in {"overexpression", "knockout", "knockdown"} 
            for pt in ad.obs["perturbation_type"]]
        ),  "Invalid 'perturbation_type' column"

        assert not ad.obs["is_control"].all(), "only controls found in test data"

        # if it says it's (not) measured, make sure it's (not) measured.
        print("Checking which genes are measured...", flush = True) 
        assert all( [    g in ad.var_names for g in ad.uns["perturbed_and_measured_genes"]] ),     "perturbed_and_measured_genes"    " not all measured"
        assert all( [not g in ad.var_names for g in ad.uns["perturbed_but_not_measured_genes"]] ), "perturbed_and_not_measured_genes sometimes measured"
    
        # If it says it's perturbed, make sure it's perturbed. 
        has_multiple_genes_hit = "perturbations_overlap" in ad.uns.keys() and ad.uns["perturbations_overlap"]
        if has_multiple_genes_hit:
            all_genes_hit = set.union(*[set(p.split(",")) for p in ad.obs["perturbation"]])      
        else:
            all_genes_hit = set(ad.obs["perturbation"]) 
        assert all( [g     in all_genes_hit for g in ad.uns["perturbed_and_measured_genes"]] ),     "perturbed_and_measured_genes"  " not perturbed"
        assert all( [g     in all_genes_hit for g in ad.uns["perturbed_but_not_measured_genes"]] ), "perturbed_and_not_measured_genes not perturbed"
    
    # Expression in `.X` should be normalized and natural-log-transformed. 
    print("Checking for log-transform and raw data...", flush = True)
    if "skip_log_check" in ad.uns.keys() and ad.uns["skip_log_check"]:
        pass
    else:
        assert ad.X.max() < 15, "Expression values too big -- did you log them?" #exp(15) is about 3 million -- too big to be a transcript count.
    # Raw data should be present in `raw`.
    assert ad.raw is not None, "raw data are missing"
    print("... done.")
    return True 


# Requirements from https://github.com/openproblems-bio/task_perturbation_prediction
#
# AnnData object
#  obs: 'dose_uM', 'timepoint_hr', 'raw_cell_id', 'hashtag_id', 'well', 'container_format', 'row', 'col', 'plate_name', 'cell_id', 'cell_type', 'split', 'donor_id', 'sm_name'
#  obsm: 'HTO_clr', 'X_pca', 'X_umap', 'protein_counts'
#  layers: 'counts'
#
# | Slot                      | Type      | Description                                    |
# |:--------------------------|:----------|:-----------------------------------------------|
# | `obs["dose_uM"]`          | `integer` | Dose in micromolar.                            |
# | `obs["timepoint_hr"]`     | `float`   | Time point measured in hours.                  |
# | `obs["raw_cell_id"]`      | `string`  | Original cell identifier.                      |
# | `obs["hashtag_id"]`       | `string`  | Identifier for hashtag oligo.                  |
# | `obs["well"]`             | `string`  | Well location in the plate.                    |
# | `obs["container_format"]` | `string`  | Format of the container (e.g., 96-well plate). |
# | `obs["row"]`              | `string`  | Row in the plate.                              |
# | `obs["col"]`              | `integer` | Column in the plate.                           |
# | `obs["plate_name"]`       | `string`  | Name of the plate.                             |
# | `obs["cell_id"]`          | `string`  | Unique cell identifier.                        |
# | `obs["cell_type"]`        | `string`  | Type of cell (e.g., B cells, T cells CD4+).    |
# | `obs["split"]`            | `string`  | Dataset split type (e.g., control, treated).   |
# | `obs["donor_id"]`         | `string`  | Identifier for the donor.                      |
# | `obs["sm_name"]`          | `string`  | Name of the small molecule used for treatment. |
# | `obsm["HTO_clr"]`         | `matrix`  | Corrected counts for hashing tags.             |
# | `obsm["X_pca"]`           | `matrix`  | Principal component analysis results.          |
# | `obsm["X_umap"]`          | `matrix`  | UMAP dimensionality reduction results.         |
# | `obsm["protein_counts"]`  | `matrix`  | Count data for proteins.                       |
# | `layers["counts"]`        | `matrix`  | Raw count data for each gene across cells.     |

def exportToOpenProblemsFormat(dataset_name: str):
    ad = load_perturbation(dataset_name)
    metadata = load_perturbation_metadata().query("name==@dataset_name")
    ad.layers["counts"] = ad.raw.X.copy()
    del ad.X
    del ad.raw
    ad.obs["dose_uM"] = np.nan
    for i in ad.obs_names:
        dose = ""
        for gene in ad.obs.loc[i, "perturbation"].split(","):
            try:
                dose_this_gene = ad[i, gene].layers["counts"]
            except KeyError:
                dose_this_gene = np.nan
            dose = dose + "," + str(dose_this_gene)
        ad.obs.loc[i, "dose_uM"] = dose
    if metadata["type"].values[0]=="timeseries":
        raise ValueError("Sorry -- timeseries data export is not supported yet because the time labels are too irregular.")
    else:
        ad.obs["timepoint_hr"] = metadata["measured_at"].str.replace("days", "").replace(" ", "").astype(float)*24
        try:
            ad.obs["cell_type"] = metadata["cell_type"]
        except KeyError:
            print("No cell type information found -- you may need to fix it manually until we update the data on Zenodo.")
    ad.obs["raw_cell_id"] = np.nan
    ad.obs["hashtag_id"] = np.nan
    ad.obs["well"] = np.nan
    ad.obs["container_format"] = np.nan
    ad.obs["row"] = np.nan
    ad.obs["col"] = np.nan
    ad.obs["plate_name"] = np.nan
    ad.obs["cell_id"] = ad.obs.index
    ad.obs["split"] = ["control" if is_control else "treated" for is_control in ad.obs["is_control"]]
    try:
        ad.obs["donor_id"] = ad.obs["donor"]
    except KeyError:
        ad.obs["donor_id"] = np.nan
    ad.obs["sm_name"] = ad.obs["perturbation"]
    assert "X_pca" in ad.obsm.keys(), "X_pca missing"
    assert "X_umap" in ad.obsm.keys(), "X_umap missing"
    return ad