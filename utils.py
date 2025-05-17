import scanpy as sc
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter, defaultdict
from scipy.sparse import csr_matrix
from scipy import sparse
import itertools as it
import anndata as ad
import tqdm

def get_leiden(adata):
    sc.pp.calculate_qc_metrics(adata, qc_vars=[], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata.X.sum(axis = 1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes = 500)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_pcs = 30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution = 0.5)
    # sc.pl.umap(adata, color=['leiden'])
    return adata.obs.leiden

def reclassify(X, N):
    X = np.array(X)
    data = np.array(X[X >= 1]).flatten()
    c = Counter(data)
    # chunk = len(data) // N
    # vals = [[] for _ in range(N)]
    # i = 0
    # curr = 0
    # # for k in sorted(c):
    # #     vals[i].append(k)
    # #     curr += c[k]
    # #     if curr >= chunk:
    # #         i += 1
    # #         curr = 0
    body = 1
    while body in c:
        body += 1
    
    chunk = (body + 1) // N + 1
    vals = [[] for _ in range(N)]
    i = 0
    curr = 0
    for k in range(1, body + 1):
        vals[i].append(k)
        if len(vals[i]) >= chunk:
            i += 1
    vals.append([body + 1, max(c) + 1])
    print(vals)
    for val, chunk_boundaries in enumerate(vals): # tqdm.tqdm(enumerate(vals)):
        if len(chunk_boundaries) == 0:
            continue

        l, r = chunk_boundaries[0], chunk_boundaries[-1]
        v = 1 if val == 0 else (l + r) // 2
        # print("chunking", l, r, val + 1)
        X = np.where((l <= X) & (X <= r), v, X)
    return X

def compare_leidens(l1, l2):
    groups1 = defaultdict(list)
    groups2 = defaultdict(list)
    for i, a in enumerate(l1):
        groups1[a].append(i)
    for i, a in enumerate(l2):
        groups2[a].append(i)
    total_matches = 0
    total = 0
    scores = {}
    for k1, v1 in groups1.items():
        subscores = {}
        for k2, v2 in groups2.items():
            s = sum(1 for i in v1 if i in v2)
            # print(s, k1, k2, len(v1), len(v2))
            subscores[k2] = (s / len(v1), s, len(v1))
            # if best_score is None or s > best_score:
            #     best_score = s 
            #     best_group = k2
        # total_matches += best_score
        # total += len(v1)
        scores[k1] = subscores
        # print(f"Best group for {k1}")
    # print(f"Overall score: {total_matches / total}")
    # print(scores)
    # st = {}
    # for k2, v2 in groups2.items():
    #     subscores = {k1: v1[k2] for k1, v1 in scores.items()}
    #     st[k2] = subscores
    taken = []
    order = sorted(list(groups1.keys()), key=lambda i: -max([j[2] for j in list(scores[i].values())]))
    for k1 in order:
        values = sorted(scores[k1].keys(), key=lambda i: -scores[k1][i][0])
        # print(values, taken)
        while values and values[0] in taken:
            values.pop(0)
        if not values or scores[k1][values[0]][0] < 0.25:
            # print(f"No match for {k1}")
            total += len(groups1[k1])
        else:
            match = values[0]
            taken.append(match)
            total_matches += scores[k1][match][1]
            total += scores[k1][match][2]
            # print(f"Best group for {k1}: {match}, score={scores[k1][match][0]}")
    # print(f"Overall score: {total_matches / total}")
    # for k1, v1 in sorted(list(groups1.items():
        
    # return scores
    return total_matches / total
