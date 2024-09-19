# Import configuration file
config <- config::get()

# parameters
percent.train <- 0.7
num.threads <- 44
num.eigenvectors <- 20
set.seed(555)

# Input files
immigrant.info <- file.path(config$path$data, "600K_immigrants.fam")
immigrant.vcf <- file.path(config$path$data, "600K_immigrants.vcf")
resident.vcf <- file.path(config$path$data, "600K_residents.vcf")

immigrant.gds <- file.path(config$path$data, "immigrant_recode.gds")
resident.gds <- file.path(config$path$data, "resident_recode.gds")
all.gds <- file.path(config$path$data, "all_recode.gds")

# read in kinship CSV
kinship <- readr::read_csv(file.path(config$path$data, "kinship.csv")) |>
    dplyr::mutate(kinship_degree = dplyr::case_when(
        kinship_degree == "unrelated" ~ 0,
        kinship_degree == "third_degree" ~ 1,
        kinship_degree == "second_degree" ~ 2,
        kinship_degree == "first_degree" ~ 3
    )) |>
    # add 'Wytham_UK_' to the beginning of each IID column
    dplyr::mutate(dplyr::across(dplyr::starts_with("IID"), ~ paste("Wytham_UK_", ., sep = "")))

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
split_data <- function(all.id, all.label, percent.train, kinship, kinship_threshold = 3) {
    # Randomly select training samples
    num_train <- round(length(all.id) * percent.train)
    train_indices <- sample(1:length(all.id), num_train)
    train.id <- all.id[train_indices]
    train.label <- all.label[train_indices]

    # Assign the rest to the test set and ensure that none in the test set are above the kinship threshold
    test.id <- all.id[-train_indices]
    test.label <- all.label[-train_indices]

    kinship_filtered <- kinship[kinship$kinship_degree <= kinship_threshold, ]
    test.id <- test.id[test.id %in% kinship_filtered$IID1 & test.id %in% kinship_filtered$IID2]
    test.label <- all.label[all.id %in% test.id]

    # Subsample the training so that it has exactly 100 samples from each class and the test set has 20 samples from each class
    train_immigrant_indices <- sample(which(train.label == "immigrant"), 100)
    train_resident_indices <- sample(which(train.label == "resident"), 100)
    test_immigrant_indices <- sample(which(test.label == "immigrant"), 20)
    test_resident_indices <- sample(which(test.label == "resident"), 20)

    train.id <- c(train.id[train_immigrant_indices], train.id[train_resident_indices])
    train.label <- c(train.label[train_immigrant_indices], train.label[train_resident_indices])
    test.id <- c(test.id[test_immigrant_indices], test.id[test_resident_indices])
    test.label <- c(test.label[test_immigrant_indices], test.label[test_resident_indices])

    return(
        list(
            train.id = train.id, test.id = test.id, train.label = train.label,
            test.label = test.label
        )
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


calculate_mean_error <- function(data) {
    data |>
        dplyr::group_by(trees, class) |>
        dplyr::summarise(mean_error = mean(error), .groups = "drop") |> # nolint: object_usage_linter.
        dplyr::filter(class != "OOB")
}


# ----------------------------------- main ----------------------------------- #


# Calculate summary statistics for the random and original models


# Initialize empty lists to store the prederror values
prederror_threshold_0 <- list()
prederror_threshold_3 <- list()


for (i in 1:20) {
    # Split data
    split <- split_data(all.id, all.label, percent.train, kinship, kinship_threshold = 0)

    # Prepare data
    data <- prepare_data(
        all.maf, split$train.id, split$test.id, split$train.label, split$test.label,
        num.eigenvectors, num.threads
    )

    # Train random forest
    rf <- train_random_forest(data$train_df, data$test_df, num.eigenvectors, FALSE)

    # Store the prederror values
    prederror_threshold_0[[i]] <- rfPermute::plotTrace(rf, plot = FALSE)$data |> dplyr::as_tibble()
}

# Run the analysis 5 times with threshold 3
for (i in 1:20) {
    # Split data
    split <- split_data(all.id, all.label, percent.train, kinship, kinship_threshold = 3)

    # Prepare data
    data <- prepare_data(
        all.maf, split$train.id, split$test.id, split$train.label, split$test.label,
        num.eigenvectors, num.threads
    )

    # Train random forest
    rf <- train_random_forest(data$train_df, data$test_df, num.eigenvectors, FALSE)

    # Store the prederror values
    prederror_threshold_3[[i]] <- rfPermute::plotTrace(rf, plot = FALSE)$data |> dplyr::as_tibble()
}

# Prepare the data for plotting
prederror_threshold_0_df <- dplyr::bind_rows(prederror_threshold_0) |>
    dplyr::mutate(threshold = 0)

prederror_threshold_3_df <- dplyr::bind_rows(prederror_threshold_3) |>
    dplyr::mutate(threshold = 3)

prederror_df <- dplyr::bind_rows(prederror_threshold_0_df, prederror_threshold_3_df) |>
    dplyr::mutate(threshold = as.factor(threshold)) |>
    dplyr::filter(class != "OOB")

# Add mean error values
mean_prederror_df <- dplyr::bind_rows(
    calculate_mean_error(prederror_threshold_0_df) |>
        dplyr::mutate(threshold = 0),
    calculate_mean_error(prederror_threshold_3_df) |>
        dplyr::mutate(threshold = 3)
) |>
    dplyr::mutate(threshold = as.factor(threshold))


# Plot the prediction error rate for different kinship thresholds
kin_error_plot <-
    ggplot2::ggplot(prederror_df, ggplot2::aes(
        x = trees, y = error, color = class, linetype = threshold
    )) +
    ggplot2::geom_line(alpha = 0.1) +
    ggplot2::geom_line(data = mean_prederror_df, ggplot2::aes(
        x = trees, y = mean_error, color = class, linetype = threshold
    ), size = 1.5) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dotted", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
        title = "Correct classification for different kinship thresholds",
        x = "Number of trees", y = "Correct classification (%)"
    ) +
    ggplot2::scale_color_manual(
        values = c(
            "immigrant" = "#5d8566",
            "resident" = "#c29007"
        )
    ) +
    ggplot2::scale_linetype_manual(
        values = c(3, 1),
        labels = c(
            "0" = "unrelated",
            "3" = "all birds"
        )
    ) +
    # y axis to percentage
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0)) +
    # also expand the x axis
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        # panel border 1px black
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1)
    )

ggplot2::ggsave(
    file.path(config$path$figures, "kinship_error_plot.png"),
    plot = kin_error_plot,
    width = 8,
    height = 6,
    dpi = 300,
    # white background
    bg = "white"
)
