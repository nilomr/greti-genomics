#' Custom ggplot2 theme
#'
#' This function returns a custom ggplot2 theme with a minimalistic design.
#'
#' @return A ggplot2 theme object.
#' @import ggplot2
#' @export
#' @examples
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point() +
#'     labs(title = "Custom ggplot2 theme", subtitle = "A minimalistic design") +
#'     titheme()
#'
#' @export
titheme <- function() {
    ggplot2::theme(
        text = ggplot2::element_text(
            size = 12, family = "Roboto Condensed",
            colour = "#272727"
        ),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        # axis.ticks = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
            fill = "transparent", colour = NA
        ),
        panel.background = ggplot2::element_rect(
            fill = "transparent", colour = NA
        ),
        aspect.ratio = .8,
        axis.title.y = ggplot2::element_text(
            size = 12,
            margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)
        ),
        axis.title.x = ggplot2::element_text(
            size = 12,
            margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)
        ),
        strip.background = ggplot2::element_rect(
            fill = "transparent", colour = NA
        ),
        strip.text = ggplot2::element_text(
            size = 12, margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 10)
        ),
        plot.subtitle = ggplot2::element_text(
            margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)
        ),
        legend.background = ggplot2::element_rect(
            fill = "transparent", colour = NA
        ),
        legend.box.background = ggplot2::element_rect(
            fill = "transparent", colour = NA
        ),
        legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.spacing.y = ggplot2::unit(0.3, "lines"),
        # title to roboto condensed, size 12
        plot.title = ggplot2::element_text(
            size = 13, face = "bold", colour = "black",
            margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)
        ),
        plot.background = ggplot2::element_rect(
            fill = "#f3f3f3", colour = NA
        )
    )
}




#' Custom color palette
#'
#' This function returns a custom color palette with 2 or 3 colors.
#'
#' @param n An integer specifying the number of colors in the palette.
#' @return A named vector of color hex codes.
#' @examples
#' # get a 2-color palette
#' titpalette(2)
#'
#' # get a 3-color palette
#' titpalette(3)
#'
#' # error: n must be 2 or 3
#' titpalette(4)
#'
#' @export
titpalette <- function(n = 3, order = NULL) {
    # Define the lookup table
    lookup <- list(
        `2` = c("#459395", "#FDA638"),
        `3` = c("#459395", "#d36044", "#fdae38"),
        `4` = c("#459395", "#d35a4a", "#FDA638", "#674A40")
    )

    # Check if n is a valid key in the lookup table
    if (!as.character(n) %in% names(lookup)) {
        stop("n must be one of ", paste(names(lookup), collapse = ", "))
    }

    # Get the corresponding palette
    palette <- lookup[[as.character(n)]]

    # Reorder the palette if the order argument is specified
    if (!is.null(order)) {
        if (length(order) != length(palette)) {
            stop("order must have the same length as the palette")
        }
        palette <- palette[order]
    }

    # Return the palette
    return(palette)
}

#' @export
reds <- c("#d35a4a", "#db7b6e", "#e49c92", "#edbdb6")
#' @export
blues <- c("#539ca1", "#386e72")
#' @export
yellows <- c("#fdae38", "#c58439")
#' @export
persian <- c("#3b6b7e", "#577383")


#' Custom Breaks for ggplot2 Axes
#'
#' These functions create custom breaks for the x and y axes in ggplot2 plots,
#' along with corresponding axis lines.
#'
#' @param x A numeric vector of data points from which to calculate breaks.
#' @param expand A ggplot2 expansion object to control the expansion of the axis.
#'               Default is `ggplot2::expansion(add = c(0, 0))` for `base_breaks_x`
#'               and `ggplot2::expansion(add = c(-0.1, 0))` for `base_breaks_y`.
#'
#' @return A list containing a `geom_segment` object for the axis line and a
#'         `scale_x_continuous` or `scale_y_continuous` object for the axis breaks.
#'
#' @examples
#' # Example usage with ggplot2
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point() +
#'     base_breaks_x(mtcars$wt) +
#'     base_breaks_y(mtcars$mpg)
#' print(p)
#'
#' @import ggplot2
#' @export
base_breaks_x <- function(x, expand = ggplot2::expansion(mult = .05)) {
    b <- pretty(x)
    d <- data.frame(y = -Inf, yend = -Inf, x = min(b), xend = max(b))
    list(
        ggplot2::geom_segment(data = d, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE),
        ggplot2::scale_x_continuous(breaks = b, expand = expand)
    )
}

#' Custom Breaks for ggplot2 Axes
#'
#' These functions create custom breaks for the x and y axes in ggplot2 plots,
#' along with corresponding axis lines.
#'
#' @param x A numeric vector of data points from which to calculate breaks.
#' @param expand A ggplot2 expansion object to control the expansion of the axis.
#'               Default is `ggplot2::expansion(add = c(0, 0))` for `base_breaks_x`
#'               and `ggplot2::expansion(add = c(-0.1, 0))` for `base_breaks_y`.
#'
#' @return A list containing a `geom_segment` object for the axis line and a
#'         `scale_x_continuous` or `scale_y_continuous` object for the axis breaks.
#'
#' @examples
#' # Example usage with ggplot2
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'     geom_point() +
#'     base_breaks_x(mtcars$wt) +
#'     base_breaks_y(mtcars$mpg)
#' print(p)
#'
#' @import ggplot2
#' @export
base_breaks_y <- function(x, expand = ggplot2::expansion(mult = .05)) {
    b <- pretty(x)
    d <- data.frame(x = -Inf, xend = -Inf, y = min(b), yend = max(b))
    list(
        ggplot2::geom_segment(data = d, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE),
        ggplot2::scale_y_continuous(breaks = b, expand = expand)
    )
}
