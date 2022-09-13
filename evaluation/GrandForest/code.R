# Install packages
devtools::install_github("SimonLarsen/grandforest")

data("HumanBioGRIDInteractionEntrezId", package = "simpIntLists")

# convert edge lists into two-column data frame
edges <- lapply(HumanBioGRIDInteractionEntrezId, function(x) {
    data.frame(
        source = as.character(x$name),
        target = as.character(x$interactors),
        stringsAsFactors = FALSE
    )
})
edges <- data.table::rbindlist(edges)

head(edges)

expr_df <- read.table("", sep = "\t")

model <- grandforest::grandforest_unsupervised(
    data = expr_df,
    graph_data = edges,
    num.trees = 10000
)

# mapping Entrez IDs to gene names
top25 <- grandforest::importance(model) %>%
    sort(decreasing = TRUE) %>%
    head(25) %>%
    as_data_frame() %>%
    rownames_to_column(var = "gene") %>%
    mutate(label = AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, gene, "SYMBOL", "ENTREZID"))

print(top25)

ggplot2::ggplot(top25, ggplot2::aes(reorder(label, -value), value)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "gene", y = "importance")

library(geomnet)

subnetwork <- filter(edges, source %in% top25$gene & target %in% top25$gene)

net_df <- ggplot2::fortify(as.edgedf(subnetwork), top25)

ggplot2::ggplot(net_df, ggplot2::aes(from_id = from_id, to_id = to_id)) +
    geomnet::geom_net(ggplot2::aes(colour = importance, label = label),
        layout.alg = "circle", directed = FALSE,
        colour = "lightblue", size = 15,
        labelon = TRUE, labelcolour = "black", vjust = 0.5, fontsize = 3
    ) +
    geomnet::theme_net()

# Extract module features and scale/center values.
d_scaled <- scale(D[, top25$gene], center = TRUE, scale = TRUE)
colnames(d_scaled) <- top25$label

# Cluster into two groups using k-means
cl <- kmeans(d_scaled, centers = 2, nstart = 20)
print(cl$cluster)

ComplexHeatmap::Heatmap(d_scaled, split = cl$cluster, name = "expression")

cl_survival <- data.frame(survival, cluster = cl$cluster)
survminer::ggsurvplot(
    survival::survfit(
        survival::Surv(
            os_time, os_event) ~ cluster, data = cl_survival),
            pval = TRUE)$plot
