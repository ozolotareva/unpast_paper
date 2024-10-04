import pandas as pd

def index_values(l):
    return [(i, x) for i, x in enumerate(l)]

def sort_indexed_values(l):
    # sort descending 
    return sorted(l, key = lambda y: y[1], reverse=True)

def format_sklearn_output(M, n_cluster, labels, transposed=True):
    # select samples based on this: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-78#Sec9
    # Operationally, this was achieved by sorting the genes in descending order by their coefficients in a given column of W (column j) and selecting only the first consecutive genes from the sorted list whose highest entry in W was the coefficient in column j
    group = {}
    if n_cluster is None:
        if transposed:
            n_cluster = len(M[:, 1])
        else:
            n_cluster = len(M[1, :])
    for j in range(n_cluster):
        # get column in W and mark index --> (index, value)
        sample_factors_j = index_values(M[:, j]) if transposed else index_values(M[j, :])
        # sort descending
        sample_factors_j = sort_indexed_values(sample_factors_j)
        samples_in_group = set()
        for i, x in sample_factors_j:
            row = M[i, :] if transposed else M[:, i]
            # include test for only max
            if row.max() == x and (len(row[row==row.max()]) == 1):
                # highest entry in M is in row j
                samples_in_group.add(i)
            else:
                # not consecutive anymore, stop adding samples
                break
            
        group[j] = {
            'samples': {labels[i] for i in samples_in_group},
            'n_samples': len(samples_in_group)
        }
    return pd.DataFrame(group).T