import bvbabel as bv
from scipy import ndimage
from scipy.stats import t
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os 
import pickle
import glob
import datetime
import os

T_DF = 42

def calc_cluster_sizes(data, threshold):
    positive_cluster_map, _ = ndimage.label(data>threshold, structure=np.ones((3,3,3))) # count adjacent non-zero voxels in all 26 neighbors
    negative_cluster_map, _ = ndimage.label(data<-1*threshold, structure=np.ones((3,3,3))) 
    pos_cluster_sizes = np.bincount(positive_cluster_map.ravel())[1:] # remove clusters of 0s
    neg_cluster_sizes = np.bincount(negative_cluster_map.ravel())[1:] # remove clusters of 0s
    cluster_sizes = np.concatenate([pos_cluster_sizes, neg_cluster_sizes])
    cluster_sizes = np.sort(cluster_sizes)
    return cluster_sizes

def permute_non_zero(org_vmp_data):
    vmp_data_copy = org_vmp_data.copy()
    # shuffle non-zero data
    vmp_nonzero = vmp_data_copy[vmp_data_copy!=0]
    np.random.shuffle(vmp_nonzero)
    # create new VMP
    vmp_data_copy[org_vmp_data!=0] = vmp_nonzero
    return vmp_data_copy

def permutation_test(org_vmp_data, p_threshold=0.025, permutations=5_000, whole_brain=False):
    threshold = t.ppf(1-p_threshold, T_DF)
    all_cluster_sizes = []
    for perm_i in range(permutations):
        if whole_brain:
            # shuffle all voxels
            vmp_data_copy = org_vmp_data.copy()
            np.random.shuffle(vmp_data_copy)
        else:
            # permute only roi voxels (all others are zero)
            vmp_data_copy = permute_non_zero(org_vmp_data)
        # count adjacent non-zero voxels
        cluster_sizes = calc_cluster_sizes(vmp_data_copy, threshold)
        all_cluster_sizes.append(cluster_sizes)
    all_cluster_sizes = np.concatenate(all_cluster_sizes)
    if len(all_cluster_sizes) > 0:
        cluster_thresh = np.percentile(all_cluster_sizes, 95)
    else:
        cluster_thresh = np.inf
    real_cluster_sizes = calc_cluster_sizes(org_vmp_data, threshold)
    significant_clusters = sorted(real_cluster_sizes[real_cluster_sizes>=cluster_thresh])
    print(f'Significant clusters: {significant_clusters}')
    if len(significant_clusters) > 0:
        significant_clusters = significant_clusters[-1]
    else:
        significant_clusters = None
    return significant_clusters, cluster_thresh

def split_file_name(file):
    file_parts = file.split('_') # for example:.VMPs/mmi_sma.vmp
    measure, roi = file_parts[0:2]
    path_separator = os.sep
    measure = measure.split(path_separator)[1]
    roi = roi.split('.')[0] # for example: VMPs/mmi_sma.vmp
    if len(file_parts)>2:
        task = file_parts[2].split('.')[0] # for example: VMPs/meanVel_sma_motorTask.vmp
    else:
        task = 'mainTask'
    return measure, roi, task
    
def save_as_pkl(file, org_vmp_data):
    pkl_file = file.lower()
    pkl_file = pkl_file.replace('.vmp', '.pkl')
    pkl_file = pkl_file.replace('vmps', 'pkls')
    if not os.path.exists(pkl_file):
        print(f'Pickling {pkl_file}')
        with open(pkl_file, 'wb') as f:
            pickle.dump(org_vmp_data, f)
        
def load_pkl(file):
    with open(file, 'rb') as f:
        org_vmp_data =  pickle.load(f)
    return org_vmp_data

def parse_args():
    parser = argparse.ArgumentParser(description="Permutation test for VMP data")
    parser.add_argument('-p', '--pvalue', type=float, default=0.025) 
    parser.add_argument('-n', '--n_permutations', type=int, default=5_000)
    parser.add_argument('-f', '--file', type=str, default=None)
    args = parser.parse_args()
    if args.pvalue:
        p_threshold = args.pvalue
    if args.n_permutations:
        permutations = args.n_permutations
    if args.file:
        files = [args.file]
    else:   
        vmp_files = glob.glob('VMPs/*.vmp')
        pkl_files = glob.glob('pkls/*.pkl')
        if len(vmp_files) == len(pkl_files):
            # if all VMPs were saved as pkl, run on pkl
            files = pkl_files
        else:
            files = vmp_files
    return files, p_threshold, permutations

def main():
    files, p_threshold, permutations = parse_args()
    all_results = pd.DataFrame()
    for file in files:
        print(f'Working on {file}')
        measure, roi, task = split_file_name(file)
        # read data with bvbabel or pickle
        if '.vmp' in file:
            _, org_vmp_data = bv.vmp.read_vmp(file)
            save_as_pkl(file, org_vmp_data)
        elif '.pkl' in file:
            org_vmp_data = load_pkl(file)
        
        if 'wholebrain' in roi.lower():
            is_whole_brain = True
            roi_pvalue = p_threshold / 10 # use p=0.005 instead of p=0.05 (two-sided)
            roi_size = org_vmp_data.size
            threshold = t.ppf(1-roi_pvalue, T_DF)
        else:
            is_whole_brain = False
            roi_pvalue = p_threshold
            roi_size = np.sum(org_vmp_data!=0)
            threshold = t.ppf(1-roi_pvalue, T_DF) 
        sig_cluster_size, cluster_thresh = permutation_test( org_vmp_data, p_threshold=roi_pvalue, 
                                                            permutations=permutations, whole_brain=is_whole_brain)
        real_cluster_sizes = calc_cluster_sizes(org_vmp_data, threshold=threshold)
        max_cluster = np.max(real_cluster_sizes) if len(real_cluster_sizes)!=0 else 0
        results_df = pd.Series({'measure': measure, 
                                'roi': roi, 
                                'task': task,
                                't_thresh': str(threshold)[:6], 
                                'p_thresh': roi_pvalue, 
                                'permutations': permutations,
                                'largest_real_significant_cluster': sig_cluster_size, 
                                'cluster_thresh': cluster_thresh, 
                                'roi_size': roi_size, 
                                'largest_real_cluster':max_cluster})
        all_results = pd.concat([all_results, results_df], axis=1)
    all_results = all_results.T
    all_results = all_results.reset_index(drop=True)
    today = datetime.datetime.now()
    p_str = str(p_threshold).split(".")[-1]
    all_results.to_csv(f'results_p{p_str}_{today.year}_{today.month}_{today.day}_{today.hour}{today.minute}.csv')

if __name__ == "__main__":
    main()
