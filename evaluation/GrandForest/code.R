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

head(edges)

expr_df <- read.table(
    "./data/real_data/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.log2_exprs_z_v6.tsv",
    sep = "\t", header = TRUE, row.names = 1
)

# NOTE sanity checks

# NA
stopifnot(!anyNA(expr_df))

# only primary cancer sites
stopifnot(sum(sapply(colnames(expr_df), function(w) {
    unlist(strsplit(w, "\\."))[4] != "01"
})) == 0)

# NOTE train model
# sink("./DESMOND2/evaluation/GrandForest/grandforest.log")
model <- grandforest::grandforest_unsupervised(
    data = t(expr_df),
    graph_data = edges,
    # 10000
    num.trees = 1000,
    # "none", "impurity_corrected", "permutation". ’impurity’ == Gini index
    importance = "impurity",
    # "dfs", "random", "random.adjusted"
    subgraph = "bfs",
    # FALSE
    replace = TRUE,
    num.threads = 1,
    verbose = TRUE,
    write.forest = FALSE,
    seed = 2022
)
# sink()

# obtain gene importance estimates of the 30 most important genes
values <- sort(grandforest::importance(model), decreasing = TRUE)[1:30]
top30 <- tibble::rownames_to_column(as.data.frame(values), var = "gene")

# mapping Entrez IDs to gene names
# top30["entrezid"] <- AnnotationDbi::mapIds(
#     org.Hs.eg.db::org.Hs.eg.db, top30$gene, "ENTREZID", "SYMBOL"
# )

top30["label"] <- top30$gene

png("./DESMOND2/evaluation/GrandForest/gene_importance.png")
plot <- ggplot2::ggplot(top30, ggplot2::aes(reorder(gene, -values), values)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "gene", y = "importance")
print(plot)
dev.off()

subnetwork <- as.data.frame(
    edges[edges$source %in% top30$gene & edges$target %in% top30$gene, ]
)

net_df <- na.omit(ggplot2::fortify(geomnet::as.edgedf(subnetwork), top30))
net_df <- net_df[apply(net_df[, c("from_id", "to_id")], 1, function(i) {
    i[1] != i[2]
}), ]

ggplot2::ggplot(net_df, ggplot2::aes(from_id = from_id, to_id = to_id)) +
    geomnet::geom_net(ggplot2::aes(colour = values, label = label),
        layout.alg = "circle", directed = FALSE,
        colour = "lightblue", size = 15,
        labelon = TRUE, labelcolour = "black", vjust = 0.5, fontsize = 3
    ) +
    geomnet::theme_net()

# NOTE stratify the patients into new endophenotypes

# Extract module features and scale/center values.
d_scaled <- scale(t(expr_df)[, top30$gene], center = TRUE, scale = TRUE)
colnames(d_scaled) <- top30$gene

# Cluster into two groups using k-means
cl <- kmeans(d_scaled, centers = 2, nstart = 20)
print(cl$cluster)

ComplexHeatmap::Heatmap(d_scaled, split = cl$cluster, name = "expression")

# NOTE real labels

# read ref table
ref_df <- read.table("./data/real_data/TCGA-BRCA_1079_17Kgenes.Xena_TCGA_PanCan.subtypes_and_signatures_v6.tsv",
    sep = "\t", header = TRUE, , row.names = 1
)
labels_df <- ref_df[, c("PAM50", "claudin_low")]
labels_df["PAM50+"] <- paste(labels_df$PAM50, labels_df$claudin_low, sep = "_")
head(labels_df)


# NOTE survival analysis

# read survival data
surv_df <- read.table("./data/real_data/TCGA-BRCA_1079.Xena_TCGA_PanCan.annotation_v6.tsv",
    sep = "\t", header = TRUE, , row.names = 1
)
head(surv_df)

cl_survival <- data.frame(surv_df[, c("OS", "OS.time")], cluster = cl$cluster)
survminer::ggsurvplot(
    survival::survfit(
        survival::Surv(OS, OS.time) ~ cluster,
        data = cl_survival
    ),
    pval = TRUE
)$plot
