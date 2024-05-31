import scanpy as sc
import numpy as np
import argparse
import csv
import pandas as pd
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from scipy.sparse import csr_matrix

def compute_metrics(clust_path):
    clust_labels = pd.read_csv(clust_path)
    clusts = clust_labels.iloc[:,2]
    labels = clust_labels["true_labels"]
    print(clusts)
    print(labels)
    ari = metrics.adjusted_rand_score(labels, clusts)
    mi = metrics.adjusted_mutual_info_score(labels, clusts)
    hs = metrics.homogeneity_score(labels, clusts)
    cs = metrics.completeness_score(labels, clusts)
    vm = metrics.v_measure_score(labels, clusts)
    fm = metrics.fowlkes_mallows_score(labels, clusts)
    scores = {"ARI": ari, "mutual_info": mi, "homogeneity" : hs, "completeness" : cs, "vmeasure" : vm, "FM" : fm, "run": clust_labels.columns[2]}
    return scores

def save_metrics(clustering_metrics, clust_path):
    csv_name = clust_path.split("/")[3]
    csv_name = csv_name.split("-labelled.csv")[0]
    combo_name = clust_path.split("/")[4]
    combo_name = combo_name.split("-labelled.csv")[0]
    print(combo_name)
    print(csv_name)
    
    csv_name = "/home/campbell/cfang/bb-rebuttal/results/supervised_metrics/singleR/" + csv_name +"/" + combo_name + ".csv"
    print(csv_name)
    with open(csv_name, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=clustering_metrics.keys())
        clustering_metrics = [clustering_metrics]
        writer.writeheader()
        writer.writerows(clustering_metrics)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--clust", action="store", type=str)
    args = parser.parse_args()
    clust_metrics = compute_metrics(args.clust)
    save_metrics(clust_metrics, args.clust)





