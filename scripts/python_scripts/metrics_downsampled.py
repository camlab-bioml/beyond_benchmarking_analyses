import scanpy as sc
import numpy as np
import argparse
import csv
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from scipy.sparse import csr_matrix


def compute_metrics(sce_path, clusters_path):
    try:
        adata = sc.read_h5ad(sce_path)
    except:
        cluster_values = []
        with open(clusters_path, newline="") as csv_file:
            reader = csv.reader(csv_file, delimiter=" ", quotechar="|")
            for row in reader: # for some reason can't open the csv file?
                row = row[0].split(",")
                cluster_values.append(row[1])
        print(cluster_values)
        pipeline_name = cluster_values[1]
        print(pipeline_name)
        silhouette = None
        ch_index = None
        db_index = None
        print("Encountered NA clusters.")
        clust_metrics = {"pipeline": pipeline_name, "silhouette":
            silhouette, "ch_index": ch_index,
                         "db_index": db_index}
        return clust_metrics
    try:
        x = adata.layers["logcounts"].toarray()
    except:
        x = adata.layers["logcounts"]
    print(adata.obs)
    cluster_values = []
    with open(clusters_path, newline="") as csv_file:
        reader = csv.reader(csv_file, delimiter=" ", quotechar="|")
        for row in reader:
            row = row[0].split(",")
            cluster_values.append(row[1])
    pipeline_name = cluster_values[0]
    cluster_values = cluster_values[1:]
    for i in range(len(cluster_values)):
        cluster_values[i] = int(cluster_values[i].strip("\""))
    cluster_values = np.array(cluster_values, dtype=int)
    try:
        silhouette = metrics.silhouette_score(x, cluster_values, metric='euclidean')
        ch_index = metrics.calinski_harabasz_score(x, cluster_values)
        db_index = metrics.davies_bouldin_score(x, cluster_values)
    except:
        silhouette = None
        ch_index = None
        db_index = None
        print("Encountered singular cluster.")
    clust_metrics = {"pipeline": pipeline_name, "silhouette":
        silhouette, "ch_index": ch_index,
                     "db_index": db_index}
    return clust_metrics


def save_metrics(clustering_metrics, csv_name):

    with open(csv_name, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=clustering_metrics.keys())
        clustering_metrics = [clustering_metrics]
        writer.writeheader()
        writer.writerows(clustering_metrics)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sce", action="store", type=str)
    parser.add_argument("--clust", action="store", type=str)
    args = parser.parse_args()
    sce_path = args.sce
    csv_name = sce_path.split("/")[-1]
    csv_name = csv_name.split(".h5ad")[0]
    sce_name = sce_path.split("/")[-3]
    print(csv_name)
    downsampleType = sce_path.split("/")[-2]
    csv_name = "/home/campbell/cfang/automl_scrna/results/metrics/" + \
               sce_name + "/" + downsampleType + "/metrics-" + \
               sce_name + "-" + csv_name + ".csv"
    print(csv_name)
    clust_metrics = compute_metrics(args.sce, args.clust)
    save_metrics(clust_metrics, csv_name)
