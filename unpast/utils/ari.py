from __future__ import annotations
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import sparse

from collections.abc import Sequence
from typing import Any, Literal


def get_similarity_ari(
    binarized_data,
    min_ARI=0.1,
    block_size=256,
    n_jobs=-1,
    verbose=True,
):
    """
    Calculate pairwise Adjusted Rand Index between binary columns.

    Parameters
    ----------
    binarized_data : pandas.DataFrame, scipy sparse matrix, or numpy.ndarray
        Binary matrix with shape (n_samples, n_features).
        Values must be 0/1 or bool.

    min_ARI : float, default=0.1
        Store only off-diagonal ARI values greater than or equal to
        this threshold. The diagonal is always equal to 1.

    block_size : int, default=256
        Number of columns in one multiplication block.

        The main dense temporary array for each worker has approximately
        shape (block_size, block_size), rather than
        (n_features, n_features).

    n_jobs : int, default=1
        Number of worker threads.

        n_jobs=1:
            sequential execution.

        n_jobs=-1:
            use all available CPU cores.

        Threads are used instead of processes so that the sparse input
        matrix is shared rather than copied for every worker.

    verbose : bool, default=True
        Show progress using tqdm when tqdm is installed.

    Returns
    -------
    ari_matrix : scipy.sparse.csr_matrix
        Symmetric sparse ARI matrix with shape
        (n_features, n_features).

        Missing off-diagonal entries represent ARI values below min_ARI,
        not necessarily exact zeroes.

    feature_names : numpy.ndarray
        DataFrame column names. For array/sparse input, integer indices.
    """
    if not -1.0 <= min_ARI <= 1.0:
        raise ValueError("min_ARI must be between -1 and 1.")

    if block_size < 1:
        raise ValueError("block_size must be a positive integer.")

    if n_jobs == 0 or n_jobs < -1:
        raise ValueError("n_jobs must be 1, 2, ..., or -1.")

    # Convert input to a CSC sparse matrix.
    # int32 is important: uint8 matrix multiplication can overflow.
    if isinstance(binarized_data, pd.DataFrame):
        feature_names = binarized_data.columns.to_numpy(copy=True)

        is_sparse_dataframe = all(
            isinstance(dtype, pd.SparseDtype)
            for dtype in binarized_data.dtypes
        )

        if is_sparse_dataframe:
            X = binarized_data.sparse.to_coo().tocsc()
            X.sum_duplicates()

            if X.data.size and not np.all(
                (X.data == 0) | (X.data == 1)
            ):
                raise ValueError(
                    "binarized_data must contain only 0/1 values."
                )

            X.eliminate_zeros()
            X = X.astype(np.int32, copy=False)

        else:
            array = binarized_data.to_numpy(copy=False)

            if not np.all((array == 0) | (array == 1)):
                raise ValueError(
                    "binarized_data must contain only 0/1 values "
                    "and no missing values."
                )

            X = sparse.csc_matrix(array, dtype=np.int32)

    elif sparse.issparse(binarized_data):
        X = binarized_data.tocsc(copy=True)
        X.sum_duplicates()

        if X.data.size and not np.all(
            (X.data == 0) | (X.data == 1)
        ):
            raise ValueError(
                "binarized_data must contain only 0/1 values."
            )

        X.eliminate_zeros()
        X = X.astype(np.int32, copy=False)
        feature_names = np.arange(X.shape[1])

    else:
        array = np.asarray(binarized_data)

        if array.ndim != 2:
            raise ValueError(
                "binarized_data must be a two-dimensional array."
            )

        if not np.all((array == 0) | (array == 1)):
            raise ValueError(
                "binarized_data must contain only 0/1 values "
                "and no missing values."
            )

        X = sparse.csc_matrix(array, dtype=np.int32)
        feature_names = np.arange(X.shape[1])

    n_samples, n_features = X.shape

    if n_samples < 2:
        raise ValueError("ARI requires at least two samples.")

    if n_features == 0:
        return (
            sparse.csr_matrix((0, 0), dtype=np.float32),
            feature_names,
        )

    # Number of ones in each binary column.
    positive_counts = np.asarray(
        X.sum(axis=0)
    ).ravel().astype(np.int64)

    def choose2(x):
        """Calculate x choose 2 in float64."""
        x = np.asarray(x, dtype=np.float64)
        return x * (x - 1.0) * 0.5

    total_pairs = choose2(n_samples)

    # For each column:
    # C(number of ones, 2) + C(number of zeroes, 2)
    marginal_pairs = (
        choose2(positive_counts)
        + choose2(n_samples - positive_counts)
    )

    blocks = [
        (start, min(start + block_size, n_features))
        for start in range(0, n_features, block_size)
    ]

    def process_block_row(block_index):
        """
        Process one row of upper-triangular matrix blocks.
        """
        i0, i1 = blocks[block_index]
        Xi = X[:, i0:i1]

        row_parts = []
        col_parts = []
        value_parts = []

        # Only calculate upper-triangular block pairs.
        for j0, j1 in blocks[block_index:]:
            Xj = X[:, j0:j1]

            # n11[a, b] is the number of samples in which
            # both binary columns are equal to 1.
            n11_block = (Xi.T @ Xj).toarray()

            diagonal_block = i0 == j0

            for local_i, i in enumerate(range(i0, i1)):
                # Exclude the diagonal and lower triangle.
                local_j_start = (
                    local_i + 1 if diagonal_block else 0
                )

                if local_j_start >= j1 - j0:
                    continue

                js = np.arange(
                    j0 + local_j_start,
                    j1,
                    dtype=np.int64,
                )

                n11 = n11_block[
                    local_i,
                    local_j_start:
                ].astype(np.int64, copy=False)

                a = positive_counts[i]
                b = positive_counts[js]

                # Binary 2 × 2 contingency table:
                #
                #             column j
                #             1       0
                # column i 1  n11     n10
                #          0  n01     n00
                n10 = a - n11
                n01 = b - n11
                n00 = n_samples - a - b + n11

                observed_index = (
                    choose2(n11)
                    + choose2(n10)
                    + choose2(n01)
                    + choose2(n00)
                )

                expected_index = (
                    marginal_pairs[i]
                    * marginal_pairs[js]
                    / total_pairs
                )

                denominator = (
                    0.5
                    * (
                        marginal_pairs[i]
                        + marginal_pairs[js]
                    )
                    - expected_index
                )

                ari = np.zeros(js.size, dtype=np.float64)

                regular = np.abs(denominator) > 1e-15

                ari[regular] = (
                    observed_index[regular]
                    - expected_index[regular]
                ) / denominator[regular]

                # Handle 0/0 cases, for example two constant
                # partitions or two equivalent partitions for n=2.
                if np.any(~regular):
                    equal_vectors = (
                        (n11 == a)
                        & (b == a)
                    )

                    complementary_vectors = (
                        (n11 == 0)
                        & ((a + b) == n_samples)
                    )

                    same_partition = (
                        equal_vectors | complementary_vectors
                    )

                    ari[~regular] = same_partition[
                        ~regular
                    ].astype(np.float64)

                # Remove tiny numerical deviations beyond [-1, 1].
                np.clip(ari, -1.0, 1.0, out=ari)

                keep = ari >= min_ARI

                if np.any(keep):
                    kept_js = js[keep]

                    row_parts.append(
                        np.full(
                            kept_js.size,
                            i,
                            dtype=np.int64,
                        )
                    )
                    col_parts.append(kept_js)
                    value_parts.append(
                        ari[keep].astype(np.float32)
                    )

        if not value_parts:
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.float32),
            )

        return (
            np.concatenate(row_parts),
            np.concatenate(col_parts),
            np.concatenate(value_parts),
        )

    if n_jobs == -1:
        worker_count = os.cpu_count() or 1
    else:
        worker_count = n_jobs

    worker_count = min(worker_count, len(blocks))

    executor = None

    if worker_count == 1:
        iterator = map(
            process_block_row,
            range(len(blocks)),
        )
    else:
        executor = ThreadPoolExecutor(
            max_workers=worker_count
        )
        iterator = executor.map(
            process_block_row,
            range(len(blocks)),
        )

    if verbose:
        try:
            from tqdm.auto import tqdm

            iterator = tqdm(
                iterator,
                total=len(blocks),
                desc="Pairwise ARI",
            )
        except ImportError:
            pass

    upper_rows = []
    upper_cols = []
    upper_values = []

    try:
        for rows, cols, values in iterator:
            if values.size:
                upper_rows.append(rows)
                upper_cols.append(cols)
                upper_values.append(values)
    finally:
        if executor is not None:
            executor.shutdown()

    if upper_values:
        rows = np.concatenate(upper_rows)
        cols = np.concatenate(upper_cols)
        values = np.concatenate(upper_values)

        upper = sparse.coo_matrix(
            (values, (rows, cols)),
            shape=(n_features, n_features),
            dtype=np.float32,
        ).tocsr()

        ari_matrix = upper + upper.T

    else:
        ari_matrix = sparse.csr_matrix(
            (n_features, n_features),
            dtype=np.float32,
        )

    ari_matrix = (
        ari_matrix
        + sparse.eye(
            n_features,
            format="csr",
            dtype=np.float32,
        )
    )

    ari_matrix.eliminate_zeros()
    ari_matrix.sort_indices()
    
    return ari_matrix, feature_names




def cluster_ari_modules(
    ari_matrix: sparse.spmatrix | np.ndarray,
    feature_names: Sequence[Any],
    *,
    algorithm: Literal["leiden", "louvain"] = "leiden",
    partition_type: Literal[
        "RBConfiguration",
        "CPM",
        "Modularity",
    ] = "RBConfiguration",
    resolution: float = 1.0,
    edge_threshold: float | None = None,
    min_module_size: int = 2,
    use_weights: bool = True,
    seed: int | None = 0,
    initial_membership: Sequence[int] | None = None,
    leiden_n_iterations: int = -1,
    leiden_max_comm_size: int = 0,
    partition_kwargs: dict[str, Any] | None = None,
) -> tuple[list[list[Any]], list[Any]]:
    """
    Cluster a sparse ARI similarity graph using Leiden or Louvain.

    Parameters
    ----------
    ari_matrix
        Square sparse ARI similarity matrix.

        The diagonal is ignored. The graph is treated as undirected.
        For asymmetric input, max(A[i, j], A[j, i]) is used.

        Only positive edge weights are supported by this wrapper.

    feature_names
        Feature names in the same order as rows and columns of
        ``ari_matrix``.

    algorithm
        Clustering algorithm:

        - "leiden": package ``leidenalg``
        - "louvain": package ``louvain``

    partition_type
        Optimised quality function:

        - "RBConfiguration":
          modularity-like objective with adjustable resolution.
          Usually a practical default.

        - "CPM":
          Constant Potts Model with adjustable resolution.

        - "Modularity":
          standard modularity. ``resolution`` is ignored.

    resolution
        Resolution for RBConfiguration or CPM.

        Larger values generally produce more and smaller modules.
        Smaller values generally produce fewer and larger modules.

    edge_threshold
        Optional additional threshold applied to the ARI matrix.

        Only edges with weight >= edge_threshold are retained.
        If None, all positive nonzero entries already stored in the
        sparse matrix are retained.

    min_module_size
        Minimum number of features in a returned module.

        Isolated features and features assigned to smaller communities
        are returned in ``not_clustered_features``.

    use_weights
        If True, ARI values are used as edge weights.

        If False, every retained edge has equal importance.

    seed
        Random seed passed to Leiden or Louvain.

    initial_membership
        Optional initial community label for every feature.

    leiden_n_iterations
        Leiden only.

        Number of optimisation iterations. A negative value runs until
        an iteration produces no further improvement.

    leiden_max_comm_size
        Leiden only.

        Maximum allowed module size. Zero means no limit.

    partition_kwargs
        Additional arguments passed to the selected partition class,
        for example:

            {"node_sizes": node_sizes}

        Parameters already represented explicitly in this function
        cannot be passed here.

    Returns
    -------
    modules
        List of modules. Each module is a list of feature names.

        Modules are sorted from largest to smallest. Feature order inside
        each module follows the original feature order.

    not_clustered_features
        Features that:

        - have no retained ARI edges; or
        - belong to communities smaller than ``min_module_size``.

        Original feature order is preserved.
    """
    if min_module_size < 1:
        raise ValueError("min_module_size must be at least 1.")

    if edge_threshold is not None and edge_threshold < 0:
        raise ValueError(
            "edge_threshold must be non-negative or None."
        )

    if resolution <= 0:
        raise ValueError("resolution must be positive.")

    if leiden_max_comm_size < 0:
        raise ValueError(
            "leiden_max_comm_size must be non-negative."
        )

    names = np.asarray(feature_names, dtype=object)

    if names.ndim != 1:
        raise ValueError(
            "feature_names must be one-dimensional."
        )

    # Make a sparse working copy so that the input matrix is not modified.
    if sparse.issparse(ari_matrix):
        A = ari_matrix.tocsr(copy=True)
    else:
        array = np.asarray(ari_matrix)

        if array.ndim != 2:
            raise ValueError(
                "ari_matrix must be two-dimensional."
            )

        A = sparse.csr_matrix(array)

    if A.shape[0] != A.shape[1]:
        raise ValueError(
            "ari_matrix must be square; "
            f"received shape {A.shape}."
        )

    n_features = A.shape[0]

    if len(names) != n_features:
        raise ValueError(
            "feature_names length must equal the dimensions "
            f"of ari_matrix: {len(names)} != {n_features}."
        )

    if n_features == 0:
        return [], []

    # Normalise sparse representation.
    A.sum_duplicates()

    # ARI = 1 on the diagonal must not produce self-loops.
    A.setdiag(0)
    A.eliminate_zeros()

    if A.data.size:
        if not np.all(np.isfinite(A.data)):
            raise ValueError(
                "ari_matrix contains NaN or infinite values."
            )

        # The previous ARI function normally stores only ARI >= 0.1.
        # Negative weights are rejected to keep behaviour consistent
        # across the supported partition types.
        if np.any(A.data < 0):
            raise ValueError(
                "Negative ARI edges are not supported by this "
                "wrapper. Filter the matrix so that only positive "
                "similarities are stored."
            )

        if np.any(A.data > 1.0 + 1e-7):
            raise ValueError(
                "ari_matrix contains values greater than 1."
            )

    # Apply an optional stricter graph threshold.
    if edge_threshold is not None and A.data.size:
        A.data[A.data < edge_threshold] = 0
        A.eliminate_zeros()

    algorithm_key = str(algorithm).strip().lower()

    if algorithm_key not in {"leiden", "louvain"}:
        raise ValueError(
            "algorithm must be 'leiden' or 'louvain'."
        )

    partition_aliases = {
        "rbconfiguration":
            "RBConfigurationVertexPartition",
        "rb":
            "RBConfigurationVertexPartition",
        "cpm":
            "CPMVertexPartition",
        "modularity":
            "ModularityVertexPartition",
    }

    partition_key = (
        str(partition_type)
        .replace("_", "")
        .replace("-", "")
        .replace(" ", "")
        .lower()
    )

    if partition_key not in partition_aliases:
        raise ValueError(
            "partition_type must be 'RBConfiguration', "
            "'CPM', or 'Modularity'."
        )

    if initial_membership is not None:
        initial_membership = list(initial_membership)

        if len(initial_membership) != n_features:
            raise ValueError(
                "initial_membership must contain one label "
                "per feature."
            )

    extra_kwargs = dict(partition_kwargs or {})

    reserved_arguments = {
        "weights",
        "seed",
        "initial_membership",
        "resolution_parameter",
        "n_iterations",
        "max_comm_size",
    }

    conflicting_arguments = (
        reserved_arguments.intersection(extra_kwargs)
    )

    if conflicting_arguments:
        raise ValueError(
            "Pass reserved parameters through the explicit "
            "function arguments, not partition_kwargs: "
            f"{sorted(conflicting_arguments)}"
        )

    try:
        import igraph as ig
    except ImportError as exc:
        raise ImportError(
            "python-igraph is required. Install it with:\n"
            "pip install igraph"
        ) from exc

    # mode="max":
    # for each pair, use max(A[i, j], A[j, i]).
    #
    # Weighted_Adjacency accepts a scipy sparse matrix and does not
    # create a dense n_features × n_features adjacency matrix.
    graph = ig.Graph.Weighted_Adjacency(
        A,
        mode="max",
        attr="weight",
        loops=False,
    )

    # The sparse working copy is no longer required.
    del A

    # No retained relationships: every feature is unclustered.
    if graph.ecount() == 0:
        return [], names.tolist()

    if algorithm_key == "leiden":
        try:
            import leidenalg as backend
        except ImportError as exc:
            raise ImportError(
                "Leiden clustering requires leidenalg. "
                "Install it with:\n"
                "pip install leidenalg"
            ) from exc

    else:
        try:
            import louvain as backend
        except ImportError as exc:
            raise ImportError(
                "Louvain clustering requires the louvain "
                "package. Install it with:\n"
                "pip install louvain"
            ) from exc

    partition_class_name = partition_aliases[partition_key]
    partition_class = getattr(
        backend,
        partition_class_name,
    )

    # ModularityVertexPartition has no resolution parameter.
    if partition_class_name != "ModularityVertexPartition":
        extra_kwargs["resolution_parameter"] = float(
            resolution
        )

    find_partition_kwargs: dict[str, Any] = {
        "initial_membership": initial_membership,
        "weights": "weight" if use_weights else None,
        "seed": seed,
        **extra_kwargs,
    }

    # These arguments exist in leidenalg.find_partition(), but not in
    # louvain.find_partition().
    if algorithm_key == "leiden":
        find_partition_kwargs["n_iterations"] = (
            leiden_n_iterations
        )
        find_partition_kwargs["max_comm_size"] = (
            leiden_max_comm_size
        )

    partition = backend.find_partition(
        graph,
        partition_class,
        **find_partition_kwargs,
    )

    membership = np.asarray(
        partition.membership,
        dtype=np.int64,
    )

    if membership.shape != (n_features,):
        raise RuntimeError(
            "The clustering backend returned an invalid "
            "membership vector."
        )

    degree = np.asarray(
        graph.degree(),
        dtype=np.int64,
    )

    communities: dict[int, list[int]] = {}

    for vertex_index, community_id in enumerate(membership):
        # An isolated feature is never considered part of a module,
        # even if the backend assigned it a community identifier.
        if degree[vertex_index] == 0:
            continue

        communities.setdefault(
            int(community_id),
            [],
        ).append(vertex_index)

    module_indices = [
        indices
        for indices in communities.values()
        if len(indices) >= min_module_size
    ]

    # Largest modules first.
    module_indices.sort(
        key=len,
        reverse=True,
    )

    clustered_mask = np.zeros(
        n_features,
        dtype=bool,
    )

    modules: list[list[Any]] = []

    for indices in module_indices:
        clustered_mask[indices] = True

        modules.append(
            [names[index] for index in indices]
        )

    not_clustered_features = names[
        ~clustered_mask
    ].tolist()

    return modules, not_clustered_features


def fclust_leiden_ari(df_bin_data, # binarized data
			algorithm="leiden",
                        min_ari= 0.25 # ignore weak overlaps with low ARI
                       ):
    #from unpast.utils.method import get_similarity_jaccard
    # calculate feature similarity matrix (Jaccard) for louvain clustering
    ari_matrix, feature_names = get_similarity_ari(df_bin_data,
                                           min_ARI=min_ari,
                                           block_size=256,
                                           n_jobs=-1,
                                           verbose=True
                                          )
    # feature clustering
    modules, not_clustered_features = cluster_ari_modules(ari_matrix,feature_names,algorithm=algorithm, resolution=1)
    print("modules: %s"%len(modules))
    print("not clustered features: %s"%len(modules))
    return modules, not_clustered_features
