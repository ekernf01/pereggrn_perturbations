import os
import pandas as pd
import scanpy as sc
import anndata


def get_data_path():
  return os.environ['PERTURBATION_PATH']
   
def set_data_path(path: str):
  if not os.path.isfile(os.path.join(path, "perturbations.csv")):
    raise FileNotFoundError("perturbations.csv should be a file in the folder whose name is provided to this function.")
  if not os.path.isfile(os.path.join(path, "nakatake", "test.h5ad")):
    raise FileNotFoundError("There should be an AnnData file at <your_input>/nakatake/test.h5ad (and others like it for the other datasets).")
  os.environ['PERTURBATION_PATH'] = path
  return

def load_perturbation_metadata():
    try:
        return pd.read_csv(os.path.join(get_data_path(), "perturbations.csv"))
    except KeyError as e:
        raise(KeyError("Before using the data you must call set_data_path('path/to/collection') to point to the perturbation data collection."))

def load_perturbation(dataset_name: str, is_train: bool = False):
    """Load a perturbation dataset. 

    Args:
        dataset_name (str): Taken from the metadata rownames.
        is_train (bool, optional): If True, this will return separate training data with no 
            perturbation (usually a timecourse). Defaults to False.

    Returns:
        anndata.AnnData: Perturbation data in a uniform format as described by `check_perturbation_dataset` or the README. 
    """
    t = "train" if is_train else "test"
    try:
        return sc.read_h5ad(os.path.join(get_data_path(), dataset_name, f"{t}.h5ad"))
    except KeyError as e:
        raise(KeyError("Dataset not found at this path. Before using the data you must call set_data_path('path/to/collection') to point to the perturbation data collection."))

def check_perturbation_dataset(dataset_name: str = None, ad: anndata.AnnData = None, is_train = False):
    """Enforce expectations on a perturbation dataset.

    Args:
        h5ad_file (str): Path to file containing perturbation data.
        ad (anndata.AnnData): AnnData object containing perturbation data.
        is_train (bool): If True, this is treated as a training dataset, which is expected 
            to contain time-series data and also extra metadata such as "timepoint".

    Raises: 
        ValueError or AssertionError for various problems with the input
    Returns: 
        True if the input data are correctly formatted
    """
    if ad is None and dataset_name is None:
        raise ValueError("Provide exactly one of ad and dataset_name")
    if not ad is None and not dataset_name is None:
        raise ValueError("Provide exactly one of ad and dataset_name")
    if ad is None and dataset_name is not None:
        # This tiny bit of recursion helps us check a dataset with separate train and test folds. 
        # The base-case: dataset_name = None and ad is not None.
        check_perturbation_dataset(ad=load_perturbation(dataset_name, is_train = False), is_train = False)
        try:
            check_perturbation_dataset(ad=load_perturbation(dataset_name, is_train = True), is_train = True)
        except FileNotFoundError:
            pass
        return
    
    # We will later select a variable number of genes based on this ranking. 
    assert "highly_variable_rank" in set(ad.var.columns), "Genes must be ranked in .var['highly_variable_rank']"
    assert all(~ad.var["highly_variable_rank"].isnull()), "Gene rankings should not be missing for any genes."
    assert all(ad.var["highly_variable_rank"]>=0), "Gene rankings must be positive integers"

    # Time
    if is_train:
        assert "timepoint" in set(ad.obs.columns), "Time-series data must have a numeric 'timepoint' column"
        assert "cell_type" in set(ad.obs.columns), "Time-series data must have a string 'cell_type' column"

    # Names of genes perturbed
    assert "perturbation" in set(ad.obs.columns), "No 'perturbation' column"
    
    # Level of those genes after perturbation
    assert "expression_level_after_perturbation" in set(ad.obs.columns), "No 'expression_level_after_perturbation' column"
    for i in ad.obs.index:
        p = ad.obs["perturbation"].astype(str)[i]
        elap = ad.obs["expression_level_after_perturbation"].astype(str)[i]
        n_levels = len(str(elap).split(","))
        n_perts =  len(str(p   ).split(","))
        assert n_levels==n_perts, f"Too many or too few expression_level_after_perturbation entries in sample {i}: {p} has {n_perts} and {elap} has {n_levels}"

    # Overexpression / knockout / knockdown
    assert "perturbation_type" in set(ad.obs.columns), "No 'perturbation_type' column"    
    assert all(
        [pt in {"overexpression", "knockout", "knockdown"} 
        for pt in ad.obs["perturbation_type"]]
    ),  "Invalid 'perturbation_type' column"

    # Boolean column is_control with both T and F
    assert "is_control"   in set(ad.obs.columns), "No 'is_control' column"
    assert bool==ad.obs["is_control"].dtype, "non-boolean 'is_control' column"
    assert       ad.obs["is_control"].any(), "no controls found"
    assert is_train or not   ad.obs["is_control"].all(), "only controls found in test data"

    # if it says it's (not) measured, make sure it's (not) measured.
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
    if "skip_log_check" in ad.uns.keys() and ad.uns["skip_log_check"]:
        pass
    else:
        assert ad.X.max() < 15, "Expression values too big -- did you log them?" #exp(15) is about 3 million -- too big to be a transcript count.
            
    # Raw data should be present in `raw`.
    assert ad.raw is not None, "raw data are missing"
    
    return True
