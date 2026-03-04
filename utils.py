from typing import List
import json
import glob
import numpy as np
import awkward as ak
from skimage.measure import block_reduce
import yaml
import glob
from sklearn.metrics import roc_curve, auc


def get_file_dict(yaml_file_path: str) -> dict:
    """
    Build a dictionary mapping process names to list of file names.
    - look recursively for all .root files in the specified path
    """
    # Get the directory of the YAML file
    with open(yaml_file_path, 'r') as f:
        raw_dict = yaml.safe_load(f)
    
    prefix = raw_dict["prefix"]
    path_dict = raw_dict["paths"]

    return {process: f"{prefix}/{path}" for process, path in path_dict.items()}
    
    # d = {}
    # for process, path in path_dict.items():
    #     # file_list = glob.glob(f"{prefix}/{path}/**/*.root", recursive=True)
        
    #     d[process] = glob.glob(f"{prefix}/{path}/**/*.root", recursive=True)
    
    # return d


def get_file_dict_old(json_file_path: str) -> dict:
    """
    Build a dictionary mapping process names to list of file names.
    - look recusrively for all .root files in the specified path
    """
    # Get the directory of the JSON file
    with open(json_file_path) as f:
        raw_dict = json.load(f)

    prefix = raw_dict["prefix"]
    path_dict = raw_dict["paths"]
    
    d = {}
    for process, path in path_dict.items():
        d[process] =  glob.glob(f"{prefix}/{path}/**/*.root", recursive=True)

    return d


def get_l1_dict(json_file_path: str) -> dict:
    # Get the directory of the JSON file
    with open(json_file_path) as f:
        raw_list = json.load(f)
    return raw_list


def get_fractions_above_threshold(scores):
    thresholds = scores.flatten()
    thresholds.sort()
    fractions = np.linspace(1, 0, len(thresholds))
    return thresholds, fractions


def get_rounded_str(value) -> str:
    """
    Applies relative rounding to a number based on its magnitude:
    - 0 to 10: 2 decimal places
    - 10 to 100: 1 decimal place
    - 100 and above: integer (0 decimal places)
    
    Parameters:
    -----------
    value : float or int
        The number to round
        
    Returns:
    --------
    float or int
        The rounded number
    """
    # Handle edge cases
    if value is None or np.isnan(value):
        return value
    
    # Take absolute value for comparison (to handle negative numbers)
    abs_value = abs(value)
    
    # Apply rounding based on magnitude
    if abs_value < 10:
        # Round to 2 decimal places
        return str(round(value, 2))
    elif abs_value < 100:
        # Round to 1 decimal place
        return str(round(value, 1))
    else:
        # Round to integer
        return str(int(round(value, 0)))


def get_dense_tower_deposits(
    tower_ieta: ak.Array, tower_iphi: ak.Array, tower_iet: ak.Array
) -> np.ndarray:
    """
    Get the dense tower deposits in cicada's eta-range from the given ragged tower arrays
    returns a numpy array of shape (n_events, 72, 56)
    """
    n_events = len(tower_ieta)

    mask = (tower_ieta >= -28) & (tower_ieta <= 28)
    tower_ieta, tower_iphi, tower_iet = tower_ieta[mask], tower_iphi[mask], tower_iet[mask]
    tower_ieta = ak.where(tower_ieta < 0, tower_ieta, tower_ieta - 1)

    tower_ieta = tower_ieta + 28
    tower_iphi = (tower_iphi + 1) % 72

    ids = np.arange(n_events, dtype=np.int32)
    ids = ak.flatten(ak.broadcast_arrays(ids, tower_ieta)[0])

    tower_ieta = ak.flatten(tower_ieta).to_numpy()
    tower_iphi = ak.flatten(tower_iphi).to_numpy()
    tower_iet = ak.flatten(tower_iet).to_numpy()

    et_tower = np.zeros((n_events, 72, 56), dtype=int)
    et_tower[ids, tower_iphi, tower_ieta] = tower_iet

    return et_tower


def get_dense_region_deposits(
    region_ieta: ak.Array, region_iphi: ak.Array, region_et: ak.Array
) -> np.ndarray:
    """
    Get the dense region deposits in cicada's eta-range from the given ragged region arrays
    returns a numpy array of shape (n_events, 18, 14)
    """
    n_events = len(region_ieta)

    ids = np.arange(n_events, dtype=np.int32)
    ids = ak.flatten(ak.broadcast_arrays(ids, region_ieta)[0])

    region_ieta = ak.flatten(region_ieta).to_numpy()
    region_iphi = ak.flatten(region_iphi).to_numpy()
    region_et = ak.flatten(region_et).to_numpy()

    et_region = np.zeros((n_events, 18, 14), dtype=int)
    et_region[ids, region_iphi, region_ieta] = region_et

    return et_region


def get_region_deposits(
    tower_ieta: ak.Array, tower_iphi: ak.Array, tower_iet: ak.Array
) -> np.ndarray:
    """
    Get the dense region deposits in cicada's eta-range from the given ragged tower arrays
    returns a numpy array of shape (n_events, 18, 14)
    """
    et_tower = get_dense_tower_deposits(tower_ieta, tower_iphi, tower_iet)
    et_region = block_reduce(et_tower, (1, 4, 4), np.sum)
    return et_region


def get_region_deposits_from_ntuple_et_array(ntuple_et_array: np.ndarray) -> np.ndarray:
    """
    Get the dense region deposits from the flat ntuple et_region array.
    This requires reshaping and permuting rows and columns because
    the et_region is stored in a different order than expected.
    The output is identical to the output of get_dense_region_deposits.
    @input ntuple_flat_et_array: np.ndarray
        A flat array of shape (n_events, 252) containing the et_region data in weird order
    @returns: np.ndarray
        A 2D array of shape (n_events, 18, 14) representing the region deposits as expected by CICADA.
    """
    row_indices = [9, 10, 11, 0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 6, 7, 8]
    col_indices = [6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 10, 11, 12, 13]
    return ntuple_et_array.reshape(-1, 18, 14)[:, row_indices, :][:, :, col_indices]


def get_anomaly_scores_ae(
        inputs: np.ndarray, outputs: np.ndarray, 
) -> np.ndarray:
    """
    Calculate anomaly score from in- and outputs of an autoencoder.
    """
    mse = np.mean((inputs - outputs) ** 2, axis=(1,2))
    score = np.log(mse * 32)
    return score


def quantize(arr: np.ndarray, precision: tuple = (16, 8)) -> np.ndarray:
    word, int_ = precision
    decimal = word - int_
    step = 1 / 2**decimal
    max_ = 2**int_ - step
    arrq = step * np.round(arr / step)
    arrc = np.clip(arrq, 0, max_)
    return arrc


def get_roc_from_scores(
        bg_scores: np.ndarray, sig_scores: np.ndarray, bkg_weights: np.ndarray = None, sig_weights: np.ndarray = None
) -> tuple[np.ndarray, np.ndarray]:
    y_trues = np.concatenate([np.ones_like(sig_scores), np.zeros_like(bg_scores)])
    y_preds = np.concatenate([sig_scores, bg_scores])
    if bkg_weights is None:
        bkg_weights = np.ones_like(bg_scores)
    if sig_weights is None:
        sig_weights = np.ones_like(sig_scores)
    weights = np.concatenate([sig_weights, bkg_weights])
    fpr, tpr, _ = roc_curve(y_trues, y_preds, sample_weight=weights)
    return fpr, tpr


def get_roc_dict(
        score_dict: dict, bg_label: str, sig_labels: List[str], weight_dict: dict = None
) -> dict:
    """
    Get a dictionary of the form
    {
        proc_0: (fpr, tpr),
        ...,
        proc_n: (fpr, tpr),
    }
    """
    d = {}
    bg_score = score_dict[bg_label]
    for proc in sig_labels:
        sig_score = score_dict[proc]
        if weight_dict is not None:
            bkg_weight = weight_dict.get(bg_label, None)
            sig_weight = weight_dict.get(proc, None)
            d[proc] = get_roc_from_scores(bg_score, sig_score, bkg_weight, sig_weight)
        else:
            d[proc] = get_roc_from_scores(bg_score, sig_score)
    return d
