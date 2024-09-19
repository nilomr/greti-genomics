# Import configuration file # nolint
config <- config::get()
box::use(R / plot[titheme, base_breaks_x, base_breaks_y])
box::use(patchwork)

progressr::handlers(progressr::handler_pbcol(
    adjust = 1.0,
    complete = function(s) cli::bg_red(cli::col_black(s)),
    incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))

# parameters
percent.train <- 0.7
num.threads <- 44
num.eigenvectors <- 15
set.seed(555)

# Input files
immigrant.info <- file.path(config$path$data, "600K_immigrants.fam")
immigrant.vcf <- file.path(config$path$data, "600K_immigrants.vcf")
resident.vcf <- file.path(config$path$data, "600K_residents.vcf")

immigrant.gds <- file.path(config$path$data, "immigrant_recode.gds")
resident.gds <- file.path(config$path$data, "resident_recode.gds")
all.gds <- file.path(config$path$data, "all_recode.gds")


# Load all the data
all.maf <- SNPRelate::snpgdsOpen(all.gds, allow.duplicate = TRUE)

# Load labels
immigrant.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
    SNPRelate::snpgdsOpen(immigrant.gds),
    "sample.id"
))
resident.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
    SNPRelate::snpgdsOpen(resident.gds),
    "sample.id"
))

all.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(all.maf, "sample.id"))
all.label <- ifelse(all.id %in% immigrant.id, "immigrant", "resident")


# Function to split data into training and testing sets
split_data <- function(all.id, all.label, percent.train) {
    # Split data by class
    labels <- unique(all.label)
    class_0 <- all.id[all.label == labels[1]]
    class_1 <- all.id[all.label == labels[2]]

    # Determine the size of the smaller class
    min_class_size <- min(length(class_0), length(class_1))

    # Undersample the majority class
    if (length(class_0) > min_class_size) {
        class_0 <- sample(class_0, min_class_size)
    } else {
        class_1 <- sample(class_1, min_class_size)
    }

    # Combine the balanced classes
    balanced_ids <- c(class_0, class_1)
    balanced_labels <- c(
        rep(labels[1], length(class_0)),
        rep(labels[2], length(class_1))
    )

    # Calculate the number of samples for training
    num_train <- round(length(balanced_ids) * percent.train)

    # Ensure equal representation of both classes in train and test sets
    num_train_per_class <- num_train %/% 2

    # Randomly select samples for training from each class
    train_indices <- c(
        sample(1:min_class_size, num_train_per_class),
        sample((min_class_size + 1):(2 * min_class_size), num_train_per_class)
    )

    train.id <- balanced_ids[train_indices]
    test.id <- balanced_ids[-train_indices]
    train.label <- balanced_labels[train_indices]
    test.label <- balanced_labels[-train_indices]

    list(
        train.id = train.id, test.id = test.id,
        train.label = train.label, test.label = test.label
    )
}
# Function to perform PCA and prepare data frames
prepare_data <- function(
    all.maf, train.id, test.id,
    train.label, test.label, num.eigenvectors, num.threads) {
    pca <- SNPRelate::snpgdsPCA(
        all.maf,
        sample.id = train.id, num.thread = num.threads, verbose = FALSE
    )
    snp_loadings <- SNPRelate::snpgdsPCASNPLoading(
        pca, all.maf,
        num.thread = num.threads, verbose = FALSE
    )
    sample_loadings <- SNPRelate::snpgdsPCASampLoading(
        snp_loadings, all.maf,
        sample.id = test.id, num.thread = num.threads, verbose = FALSE
    )
    train_eigenvects <- data.frame(
        sample.id = pca$sample.id, pca$eigenvect[, 1:num.eigenvectors]
    )
    test_eigenvects <- data.frame(
        sample.id = sample_loadings$sample.id, sample_loadings$eigenvect[, 1:num.eigenvectors]
    )
    train_df <- data.frame(sample.id = train.id, pop = train.label, train_eigenvects[, -1])
    test_df <- data.frame(sample.id = test.id, pop = test.label, test_eigenvects[, -1])
    list(train_df = train_df, test_df = test_df)
}

# Function to train a random forest with optional randomized labels
train_random_forest <- function(
    train_df, test_df, num.eigenvectors, randomize_labels = FALSE) {
    if (randomize_labels) {
        train_df$pop <- sample(train_df$pop)
    }

    return(randomForest::randomForest(
        train_df[, 3:(2 + num.eigenvectors)],
        y = as.factor(train_df$pop),
        data = train_df,
        xtest = test_df[, 3:(2 + num.eigenvectors)],
        ytest = as.factor(test_df$pop),
        ntree = 500,
        replace = FALSE,
        proximity = TRUE
    ))
}

# Function to run a single iteration
run_iteration <- function(
    all.id, all.label, all.maf, percent.train, num.eigenvectors, num.threads, randomize_labels) {
    # Split data
    split <- split_data(all.id, all.label, percent.train)

    # Prepare data
    data <- prepare_data(
        all.maf, split$train.id, split$test.id, split$train.label, split$test.label,
        num.eigenvectors, num.threads
    )

    # Train random forest
    rf <- train_random_forest(data$train_df, data$test_df, num.eigenvectors, randomize_labels)

    # Extract OOB error rates
    oob_data <- rfPermute::plotTrace(rf, plot = FALSE)$data

    # Return both OOB data, the random forest model, and the train/test labels
    return(list(
        oob_data = oob_data, rf_model = rf,
        train_labels = split$train.label, test_labels = split$test.label
    ))
}

# Main function to run multiple iterations
run_analysis <- function(
    all.id, all.label, all.maf, percent.train,
    num.eigenvectors, num.threads,
    n_iterations = 1) {
    set.seed(123) # For reproducibility

    # Initialize progress bar
    progressr::with_progress({
        p <- progressr::progressor(along = 1:(2 * n_iterations))

        # Run iterations with original labels
        original_results <- replicate(n_iterations,
            {
                p()
                run_iteration(
                    all.id, all.label, all.maf, percent.train,
                    num.eigenvectors, num.threads, FALSE
                )
            },
            simplify = FALSE
        )

        # Run iterations with randomized labels
        random_results <- replicate(n_iterations,
            {
                p()
                run_iteration(
                    all.id, all.label, all.maf, percent.train,
                    num.eigenvectors, num.threads, TRUE
                )
            },
            simplify = FALSE
        )
    })

    # Separate OOB data, RF models, and train/test labels
    original_oob <- lapply(original_results, function(x) x$oob_data)
    original_rf_models <- lapply(original_results, function(x) x$rf_model)
    original_train_labels <- lapply(original_results, function(x) x$train_labels)
    original_test_labels <- lapply(original_results, function(x) x$test_labels)

    random_oob <- lapply(random_results, function(x) x$oob_data)
    random_rf_models <- lapply(random_results, function(x) x$rf_model)
    random_train_labels <- lapply(random_results, function(x) x$train_labels)
    random_test_labels <- lapply(random_results, function(x) x$test_labels)

    # Combine OOB results
    oob_df <- dplyr::as_tibble(do.call(rbind, original_oob))
    random_oob_df <- dplyr::as_tibble(do.call(rbind, random_oob))

    list(
        oob_df = oob_df,
        random_oob_df = random_oob_df,
        original_rf_models = original_rf_models,
        random_rf_models = random_rf_models,
        original_train_labels = original_train_labels,
        original_test_labels = original_test_labels,
        random_train_labels = random_train_labels,
        random_test_labels = random_test_labels
    )
}

calculate_mean_error <- function(data) {
    data |>
        dplyr::group_by(trees, class) |>
        dplyr::summarise(mean_error = mean(error), .groups = "drop") |> # nolint: object_usage_linter.
        dplyr::filter(class != "OOB")
}


# ----------------------------------- main ----------------------------------- #

# TODO:
# - RF error rate plot
# - Summary statistics for random and original models
# - MDS plot
# - Confusion matrix
# - Calculate some sort of distance between points and cluster centroids

# Check if results file exists
results_file <- file.path(config$path$data, "rf_results.rds")

if (file.exists(results_file)) {
    # Load the results
    results <- readRDS(results_file)
} else {
    # Run the analysis
    results <- run_analysis(all.id, all.label, all.maf, percent.train,
        num.eigenvectors, num.threads,
        n_iterations = 20
    )

    # Save the results
    saveRDS(results, results_file)
}


# Calculate mean error for all iterations
mean_oob_df <- calculate_mean_error(results$oob_df)
mean_random_oob_df <- calculate_mean_error(results$random_oob_df)

# Create the plot data
plot_data <- list(
    oob_df = results$oob_df |>
        dplyr::filter(class != "OOB"),
    random_oob_df = results$random_oob_df |>
        dplyr::filter(class != "OOB"),
    mean_oob_df = mean_oob_df,
    mean_random_oob_df = mean_random_oob_df |> dplyr::group_by(trees) |>
        dplyr::summarise(mean_error = mean(mean_error), .groups = "drop")
)

# Create the plot
oobplot <- ggplot2::ggplot() +
    ggplot2::geom_line(
        data = plot_data$random_oob_df,
        ggplot2::aes(x = trees, y = error, color = "randomised"),
        alpha = 0.1
    ) +
    ggplot2::geom_line(
        data = plot_data$mean_random_oob_df,
        ggplot2::aes(x = trees, y = mean_error, color = "randomised"),
        linewidth = 1.5
    ) +
    ggplot2::geom_line(
        data = plot_data$oob_df,
        ggplot2::aes(x = trees, y = error, color = class),
        alpha = 0.1
    ) +
    ggplot2::geom_line(
        data = plot_data$mean_oob_df,
        ggplot2::aes(x = trees, y = mean_error, color = class),
        linewidth = 1.5
    ) +
    ggplot2::geom_hline(
        yintercept = 50,
        linetype = "dotted",
        color = "black"
    ) +
    ggplot2::labs(
        title = "OOB error rate (Solid: Original, Dashed: Randomized)",
        x = "Number of trees",
        y = "Percent correct"
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "immigrant" = "#5d8566",
            "resident" = "#c29007",
            "randomised" = "#919191"
        )
    ) +
    base_breaks_y(40:100, expand = ggplot2::expansion(mult = .05)) +
    base_breaks_x(0:500, expand = ggplot2::expansion(mult = .05)) +
    titheme() +
    ggplot2::theme(
        aspect.ratio = 1.8,
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.95),
        legend.justification = c("right", "top"),
        legend.title = ggplot2::element_text(hjust = 1)
    ) +
    ggplot2::guides(
        color = ggplot2::guide_legend(
            title = "Class",
            label.position = "left",
            label.hjust = 1
        )
    )

# Save the plot
ggplot2::ggsave(file.path(config$path$figures, "oob_plot.jpg"),
    plot =
        oobplot,
    width = 8,
    height = 6
)

# get data from tree 100 to tree 500 and plot the distribution of percent correct


oob_df <- dplyr::bind_rows(
    list(
        oob_df = plot_data$oob_df,
        random_oob_df = plot_data$random_oob_df
    ),
    .id = "source"
) |>
    dplyr::filter(trees >= 100 & trees <= 500)

# get the MDS plot data
# run a single random forest model:

set.seed(123)
split <- split_data(all.id, all.label, percent.train)
data <- prepare_data(
    all.maf, split$train.id, split$test.id, split$train.label, split$test.label,
    num.eigenvectors, num.threads
)

# Use all available data as here we are not interested in the model's performance
data <- data$train_df |> dplyr::bind_rows(data$test_df)

set.seed(123)
rf <- randomForest::randomForest(
    data[, 3:(2 + num.eigenvectors)],
    y = as.factor(data$pop),
    ntree = 300,
    replace = FALSE,
    proximity = TRUE,
    oob.prox = TRUE
)

# cmd the proximity matrix from the model and plot
prox <- 1 - rf$proximity
mds <- stats::cmdscale(prox, k = 2)

# plot the MDS
mds_df <- data.frame(mds, class = rf$y, predicted = rf$predicted)


# plot the distribution of percent correct by class

oob_dist_plot <- ggplot2::ggplot() +
    ggdist::stat_slab(
        data = oob_df |> dplyr::filter(source == "oob_df"),
        ggplot2::aes(y = error, fill = class),
        position = "identity", alpha = 0.65,
        density = ggdist::density_bounded(bandwidth = 1, trim = FALSE)
    ) +
    ggdist::stat_slab(
        data = oob_df |> dplyr::filter(source == "random_oob_df"),
        ggplot2::aes(y = error, fill = "randomised"),
        position = "identity", alpha = 0.65,
        density = ggdist::density_bounded(bandwidth = 1, trim = FALSE),
    ) +
    ggdist::stat_pointinterval(
        data = oob_df |> dplyr::filter(source == "oob_df"),
        ggplot2::aes(y = error, color = class),
        position = ggplot2::position_dodge(
            width = .5, preserve = "single"
        )
    ) +
    ggdist::stat_pointinterval(
        data = oob_df |> dplyr::filter(source == "random_oob_df"),
        ggplot2::aes(y = error, color = "randomised"),
        position = ggplot2::position_dodge(
            width = .5, preserve = "single"
        )
    ) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = "black") +
    ggplot2::labs(
        title = "Distribution of percent correct by class",
        y = "Percent correct", x = "Density"
    ) +
    ggplot2::scale_fill_manual(
        values = c(
            "immigrant" = "#5d8566", "resident" = "#c29007", "randomised" = "#716e6e"
        ), labels = c("immigrant", "resident", "random\nperformance")
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "immigrant" = "#5d8566", "resident" = "#c29007", "randomised" = "grey"
        )
    ) +
    base_breaks_y(40:100, expand = ggplot2::expansion(mult = .05)) +
    base_breaks_x(0:1, expand = ggplot2::expansion(add = c(-0.1, 0))) +
    titheme() +
    ggplot2::theme(
        aspect.ratio = 2.3,
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.95),
        legend.justification = c("right", "top"),
        legend.title = ggplot2::element_text(hjust = 1)
    ) +
    ggplot2::guides(
        fill = ggplot2::guide_legend(
            title = "Class",
            label.position = "left",
            label.hjust = 1
        ),
        color = "none"
    )


mds_plot <-
    ggplot2::ggplot(mds_df, ggplot2::aes(x = X1, y = X2, color = predicted)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = class, stroke = NA), # nolint
        shape = 21,
        size = 1,
        alpha = 0.7
    ) +
    ggplot2::scale_fill_manual(values = c(
        "immigrant" = "#5d8566", "resident" = "#c29007"
    )) +
    # add points (real label) and circles (predicted label)
    ggplot2::geom_point(
        ggplot2::aes(color = predicted),
        shape = 21,
        size = 3,
        stroke = 0.5,
        alpha = 0.7,
        show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(
        "immigrant" = "#5d8566", "resident" = "#c29007"
    )) +
    titheme() +
    ggplot2::labs(
        title = "MDS plot of random forest proximity matrix",
        x = "MDS1", y = "MDS2"
    ) +
    ggplot2::guides(
        color = ggplot2::guide_legend(
            title = "Predicted",
            label.position = "left",
            label.hjust = 1,
            order = 1,
            override.aes = list(
                color = c("#5d8566", "#c29007"),
                fill = c(NA, NA),
                shape = 21,
                size = 3.5,
                stroke = 0.5
            )
        ),
        fill = ggplot2::guide_legend(
            title = "Actual",
            label.position = "left",
            label.hjust = 1,
            order = 2,
            override.aes = list(
                color = c("#5d8566", "#c29007"),
                fill = c("#5d8566", "#c29007"),
                shape = 21,
                size = 1.2,
                alpha = 0.9
            )
        )
    ) +
    base_breaks_x(range(-0.2, 0.6), expand = ggplot2::expansion(mult = .05)) +
    base_breaks_y(range(-0.2, 0.6), expand = ggplot2::expansion(mult = .05)) +
    ggplot2::theme(
        aspect.ratio = 1.5,
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.title = ggplot2::element_text(hjust = 1)
    )


rf_combi_plot <- oob_dist_plot +
    ggplot2::theme(plot.margin = ggplot2::margin(r = 20)) +
    mds_plot +
    patchwork::plot_annotation(
        tag_levels = "A"
    ) &
    ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        plot.tag = ggplot2::element_text(face = "plain")
    )

# save the plot
ggplot2::ggsave(
    file.path(config$path$figures, "mds_combi_plot.jpg"),
    plot = rf_combi_plot,
    width = 15,
    height = 11,
    units = "cm",
    bg = "white"
)
