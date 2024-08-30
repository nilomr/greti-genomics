# Import configuration file
config <- config::get()

# parameters
percent.train <- 0.7
num.threads <- 44
num.eigenvectors <- 5
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
# TODO balence so that there are .5 per class in both train and test sets
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
    balanced_labels <- c(rep(labels[1], length(class_0)), rep(labels[2], length(class_1)))

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
prepare_data <- function(all.maf, train.id, test.id, train.label, test.label, num.eigenvectors, num.threads) {
    pca <- SNPRelate::snpgdsPCA(all.maf, sample.id = train.id, num.thread = num.threads)
    snp_loadings <- SNPRelate::snpgdsPCASNPLoading(pca, all.maf, num.thread = num.threads)
    sample_loadings <- SNPRelate::snpgdsPCASampLoading(snp_loadings, all.maf, sample.id = test.id, num.thread = num.threads)
    train_eigenvects <- data.frame(sample.id = pca$sample.id, pca$eigenvect[, 1:num.eigenvectors])
    test_eigenvects <- data.frame(sample.id = sample_loadings$sample.id, sample_loadings$eigenvect[, 1:num.eigenvectors])
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
    return(list(oob_data = oob_data, rf_model = rf, train_labels = split$train.label, test_labels = split$test.label))
}

# Main function to run multiple iterations
run_analysis <- function(
    all.id, all.label, all.maf, percent.train,
    num.eigenvectors, num.threads,
    n_iterations = 1) {
    set.seed(123) # For reproducibility

    # Run iterations with original labels
    original_results <- replicate(n_iterations,
        run_iteration(
            all.id, all.label, all.maf, percent.train,
            num.eigenvectors, num.threads, FALSE
        ),
        simplify = FALSE
    )

    # Run iterations with randomized labels
    random_results <- replicate(n_iterations,
        run_iteration(
            all.id, all.label, all.maf, percent.train,
            num.eigenvectors, num.threads, TRUE
        ),
        simplify = FALSE
    )

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


# Run the analysis
results <- run_analysis(all.id, all.label, all.maf, percent.train, num.eigenvectors, num.threads)


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
    mean_random_oob_df = mean_random_oob_df
)

# Create the plot
oobplot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = plot_data$oob_df, ggplot2::aes(
        x = trees, y = error, color = class
    ), alpha = 0.1) +
    ggplot2::geom_line(data = plot_data$random_oob_df, ggplot2::aes(
        x = trees, y = error, color = class
    ), alpha = 0.1, linetype = "dashed") +
    ggplot2::geom_line(data = plot_data$mean_oob_df, ggplot2::aes(
        x = trees, y = mean_error, color = class
    ), linewidth = 1.5) +
    ggplot2::geom_line(data = plot_data$mean_random_oob_df, ggplot2::aes(
        x = trees, y = mean_error, color = class
    ), linewidth = 1.5, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
        title = "OOB error rate (Solid: Original, Dashed: Randomized)",
        x = "Number of trees", y = "Percent correct"
    ) +
    ggplot2::scale_color_manual(values = c(
        "OOB" = "#4e4e4e",
        "immigrant" = "#5d8566", "resident" = "#c29007"
    )) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = "red") +
    ggplot2::theme(aspect.ratio = 1, text = ggplot2::element_text(size = 12)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1))

# Save the plot
ggplot2::ggsave(file.path(config$path$figures, "oob_plot.jpg"),
    plot =
        oobplot,
    width = 8,
    height = 6
)


# Calculate summary statistics for the random and original models



# randomForest::MDSplot(results$original_rf_models[[1]], results$original_train_labels[[1]])

# rfPermute::confusionMatrix(results$original_rf_models[[1]])
# # rfPermute::plotPredictedProbs(output.forest, bins = 20, plot = TRUE)
# # rfPermute::plotProximity(output.forest)
# # summary(output.forest)
# # rfPermute::plotTrace(output.forest)

# results$original_rf_models[[1]]$importance
# results$original_rf_models[[1]]$confusion
