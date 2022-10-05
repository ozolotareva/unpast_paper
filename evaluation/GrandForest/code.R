# Install packages
for (pkg in c("grandforest", "geomnet")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        invisible(switch(pkg,
            "grandforest" = devtools::install_github("SimonLarsen/grandforest"),
            "geomnet" = devtools::install_version(
                "geomnet",
                version = "0.3.1", repos = "http://cran.us.r-project.org"
            )
        ))
    }
}

data("HumanBioGRIDInteractionOfficial", package = "simpIntLists")

# convert edge lists into two-column data frame
edges <- lapply(HumanBioGRIDInteractionOfficial, function(x) {
    data.frame(
        source = as.character(x$name),
        target = as.character(x$interactors),
        stringsAsFactors = FALSE
    )
})
edges <- data.table::rbindlist(edges)

read_data <- function(dataset) {
    file_name <- ifelse(dataset == "GDC",
        "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv",
        "METABRIC_1904_17Kgenes.log2_exprs_z_v6.tsv"
    )

    expr_df <- read.table(paste0("/root/projects/data/real_data/", file_name),
        sep = "\t", header = TRUE, row.names = 1
    )

    # NOTE sanity checks
    # NA
    stopifnot(!anyNA(expr_df))
    # only primary cancer sites
    if (dataset == "GDC") {
        stopifnot(sum(sapply(colnames(expr_df), function(w) {
            unlist(strsplit(w, "\\."))[4] != "01"
        })) == 0)
    }

    dir.create(paste0("/root/projects/data/outputs/", dataset),
        showWarnings = FALSE
    )

    # read ref table
    file_name <- ifelse(dataset == "GDC",
        "TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv",
        "METABRIC_1904_17Kgenes.subtypes_and_signatures_v6.tsv"
    )
    labels_df <- read.table(paste0("/root/projects/data/real_data/", file_name),
        sep = "\t", header = TRUE, row.names = 1
    )

    # survival data
    file_name <- ifelse(dataset == "GDC",
        "TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv",
        "METABRIC_1904.annotation_v6.tsv"
    )
    surv_df <- read.table(paste0("/root/projects/data/real_data/", file_name),
        sep = "\t", header = TRUE, row.names = 1
    )

    data_object <- list(expr = expr_df, labels = labels_df, survival = surv_df)
    return(data_object)
}

create_analysis <- function(dataset, data, par_tibble, comb_n, n = 30) {
    dir.create(file.path("/root/projects/data/outputs", dataset, comb_n),
        showWarnings = FALSE
    )

    model <- grandforest::grandforest_unsupervised(
        data = t(data$expr),
        graph_data = edges,
        num.trees = as.numeric(par_tibble["num.trees"]),
        importance = as.character(par_tibble["importance"]),
        subgraph = as.character(par_tibble["subgraph"]),
        num.threads = 1,
        verbose = TRUE,
        write.forest = FALSE,
        seed = 2022
    )

    # gene importance estimates
    values <- sort(grandforest::importance(model), decreasing = TRUE)[1:n]
    top_n <- tibble::rownames_to_column(as.data.frame(values), var = "gene")
    top_n["label"] <- top_n$gene

    png(file.path(
        "/root/projects/data/outputs", dataset, comb_n, "gene_importance.png"
    ))
    plot <- ggplot2::ggplot(top_n, ggplot2::aes(
        reorder(gene, -values), values
    )) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(x = "gene", y = "importance")
    print(plot)
    dev.off()

    subnetwork <- as.data.frame(
        edges[edges$source %in% top_n$gene & edges$target %in% top_n$gene, ]
    )

    net_df <- na.omit(ggplot2::fortify(geomnet::as.edgedf(subnetwork), top_n))
    net_df <- net_df[apply(net_df[, c("from_id", "to_id")], 1, function(i) {
        i[1] != i[2]
    }), ]

    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "gene_module_network.png"
    ))
    p <- ggplot2::ggplot(net_df, ggplot2::aes(
        from_id = from_id, to_id = to_id
    )) +
        geomnet::geom_net(ggplot2::aes(colour = values, label = label),
            layout.alg = "circle", directed = FALSE,
            colour = "lightblue", size = 15,
            labelon = TRUE, labelcolour = "black", vjust = 0.5, fontsize = 3
        ) +
        geomnet::theme_net()
    print(p)
    dev.off()

    kmeans_obj <- create_endophenotypes(dataset, data, top_n, comb_n)

    results_object <- list(top_n = top_n, model = model, cluster = kmeans_obj)
    return(results_object)
}

create_endophenotypes <- function(dataset, data, top_n, comb_n) {
    set.seed(2022)

    # Extract module features and scale/center values.
    d_scaled <- scale(t(data$expr)[, top_n$gene], center = TRUE, scale = TRUE)
    colnames(d_scaled) <- top_n$gene

    # Elbow Method
    tot_withinss <- sapply(1:10, function(k) {
        kmeans(d_scaled, centers = k)$tot.withinss
    })
    elbow_df <- data.frame(k = 1:10, tot_withinss = tot_withinss)
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "elbow_method.png"
    ))
    plot <- ggplot2::ggplot(elbow_df, ggplot2::aes(x = k, y = tot_withinss)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(breaks = 1:10)
    print(plot)
    dev.off()

    # Silhouette Analysis
    sil_width <- sapply(2:10, function(k) {
        cluster::pam(d_scaled, k = k)$silinfo$avg.width
    })
    sil_df <- data.frame(k = 2:10, sil_width = sil_width)
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "silhouette_analysis.png"
    ))
    p <- ggplot2::ggplot(sil_df, ggplot2::aes(x = k, y = sil_width)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_continuous(breaks = 2:10)
    print(p)
    dev.off()

    # Gap Statistic
    gap_stat <- cluster::clusGap(d_scaled,
        FUN = kmeans, nstart = 25,
        K.max = 10, B = 50
    )
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "gap_statistic.png"
    ))
    print(factoextra::fviz_gap_stat(gap_stat))
    dev.off()

    # TODO pick the right k value automatically
    k <- 5

    # Cluster into groups using k-means
    cl <- kmeans(d_scaled, centers = k)
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "k_means.png"
    ))
    print(factoextra::fviz_cluster(cl, geom = "point", data = d_scaled))
    dev.off()

    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "heatmap.png"
    ))
    ht <- ComplexHeatmap::Heatmap(d_scaled,
        split = cl$cluster, name = "expression",
        show_row_names = FALSE
    )
    ComplexHeatmap::draw(ht)
    dev.off()

    # Label overlay
    pca_res <- prcomp(d_scaled)
    pca_df <- data.frame(
        cluster = as.factor(kmeans(d_scaled, centers = k)$cluster)
    )
    pca_df$PC1 <- pca_res$x[, 1]
    pca_df$PC2 <- pca_res$x[, 2]
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "pca.png"
    ))
    plot <- ggplot2::ggplot(ggplot2::aes(
        x = PC1, y = PC2, col = cluster
    ), data = pca_df) +
        ggplot2::geom_point() +
        ggplot2::facet_grid(. ~ data$labels[, "PAM50"])
    print(plot)
    dev.off()

    # survival analysis
    cl_survival <- data.frame(
        data$survival[, c("OS", "OS.time")],
        cluster = cl$cluster
    )
    png(file.path(
        "/root/projects/data/outputs",
        dataset, comb_n, "survival.png"
    ))
    plot <- survminer::ggsurvplot(
        survival::survfit(
            survival::Surv(OS.time, OS) ~ cluster,
            data = cl_survival
        ), cl_survival,
        pval = TRUE
    )
    print(plot$plot)
    dev.off()

    return(cl)
}

# NOTE Main
data_tcga <- read_data("GDC")
data_mbr <- read_data("Mbr")

par_df <- purrr::cross_df(list(
    num.trees = seq(1000, 10000, length.out = 3),
    importance = c("impurity", "impurity_corrected", "permutation"),
    subgraph = c("bfs", "dfs", "random")
))

par_comb <- apply(par_df, 1, function(x) {
    paste(x, collapse = ";")
})

write.table(par_df, "/root/projects/data/outputs/par_df.tsv",
    sep = "\t",
    row.names = FALSE
)

library(doMC)

doMC::registerDoMC(nrow(par_df))

time_tcga <- data.frame(par_comb = par_comb, time = rep(0, nrow(par_df)))
res_tcga <- foreach::foreach(i = seq_len(nrow(par_df))) %dopar% {
    start <- Sys.time()
    create_analysis("GDC", data_tcga, par_df[i, ], i)
    end <- Sys.time()
    time_tcga[i, "time"] <- as.numeric(end - start)
}
write.csv(time_tcga, "/root/projects/data/outputs/time_tcga.tsv", sep = "\t")

save.image("/root/projects/data/outputs/grandforest1.RData")

time_mbr <- data.frame(par_comb = par_comb, time = rep(0, nrow(par_df)))
res_mbr <- foreach::foreach(i = seq_len(nrow(par_df))) %dopar% {
    start <- Sys.time()
    create_analysis("Mbr", data_mbr, par_df[i, ], i)
    end <- Sys.time()
    time_mbr[i, "time"] <- as.numeric(end - start)
}
write.csv(time_mbr, "/root/projects/data/outputs/time_mbr.tsv", sep = "\t")

save.image("/root/projects/data/outputs/grandforest2.RData")

# TODO export data to calc metrics
