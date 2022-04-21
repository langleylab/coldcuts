#' Plot the ontology graph
#'
#' Creates a `ggplot` plot of the ontology graph/tree
#'
#' @param segmentation a \code{segmentation} class object
#' @param circular a logical indicating whether the layout should be circular
#'
#' @return a `ggplot` object with the graph
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggraph ggraph geom_edge_elbow geom_node_point geom_node_text
#' @importFrom ggplot2 aes_string theme_bw theme element_blank ggtitle
#'    scale_fill_manual coord_fixed
#'
#' @export

ontology_plot <- function(segmentation, circular = TRUE) {
  o <- ontology(segmentation)
  cols <- o$col
  names(cols) <- o$acronym
  
  onto_graph <- data.frame("from" = o[as.character(o$parent_structure_id), "acronym"], "to" = o$acronym)
  
  if (any(is.na(onto_graph[, 1]))) {
    onto_graph[which(is.na(onto_graph[, 1])), 1] <- onto_graph[which(is.na(onto_graph[, 1])), 2]
  } else if (any(is.na(onto_graph[, 2]))) {
    onto_graph[which(is.na(onto_graph[, 2])), 2] <- onto_graph[which(is.na(onto_graph[, 2])), 1]
  }
  
  g <- graph_from_data_frame(onto_graph)
  g <- as_tbl_graph(g)
  
  p <- ggraph(g, "dendrogram", circular = circular) +
    geom_edge_elbow() +
    geom_node_point(aes_string(fill = "name"), shape = 21, color = "black", size = 8) +
    geom_node_text(aes_string(label = "name"), color = "black", size = 2) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )
  
  return(p)
}


#' Plot slices from a segmentation
#'
#' Creates a `ggplot` plot of up to 3 different slices (one per anatomical plane) of the segmentation
#'
#' @param segmentation a \code{segmentation} class object
#' @param s_slice numeric, the slice number (index) on the sagittal plane
#' @param c_slice numeric, the slice number (index) on the coronal plane
#' @param a_slice numeric, the slice number (index) on the axial plane
#' @param smooth logical, should polygons be smoothed before drawing? Default is \code{TRUE}. Uses kernel smoothing.
#' @param smoothness numeric, the argument passed to the kernel smoothing function.
#' @param show_axis_rulers logical, should rulers for each axis be shown in each view? Default is \code{TRUE}.
#' @param show_outline logical, should the maximum segmentation outline be plotted? Default is \code{TRUE}.
#' @param show_labels logical, should structure acronym labels be shown on the plot? Default is \code{FALSE}.
#' @param label_size numeric, the size of the labels.
#' @param minsize numeric, the minimum number of points to create a polygon. Any structure in the slice with fewer points than this number will be discarded. Default is 10.
#' @param wrap_options numeric, the number of rows for the facet wrapping. Default is 1.
#'
#' @return a `ggplot` object with the graph
#'
#' @importFrom polylabelr poi
#' @importFrom ggplot2 ggplot geom_polygon aes_string theme_bw
#'    scale_fill_gradientn labs theme element_blank ggtitle scale_fill_manual
#'    coord_fixed scale_x_reverse facet_wrap geom_vline geom_hline
#' @importFrom gridExtra grid.arrange
#' @importFrom ggrepel geom_text_repel
#'
#' @export

seg_plot <- function(segmentation, s_slice = NULL, c_slice = NULL,
                     a_slice = NULL, smooth = TRUE, smoothness = 3,
                     show_axis_rulers = TRUE, show_outline = TRUE,
                     show_labels = FALSE, label_size = 2, minsize = 10,
                     wrap_options = 1) {
  if (smooth & !"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")
  if (class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")
  
  planes <- c("sagittal", "coronal", "axial")
  if (is.null(s_slice)) s_slice <- NULL
  if (is.null(c_slice)) c_slice <- NULL
  if (is.null(a_slice)) a_slice <- NULL
  selected_slices <- list(s_slice, c_slice, a_slice)
  planes <- planes[!sapply(selected_slices, is.null)]
  selected_slices <- unlist(selected_slices[!sapply(selected_slices, is.null)])
  names(selected_slices) <- planes
  
  # Select the slices and add axis information (important for faceting the plot later)
  slicelist <- lapply(planes, function(p) {
    sagittal_list <- lapply(segmentation@slices[[p]][[selected_slices[p]]], function(x) {
      polylist <- lapply(x, function(y) {
        if (nrow(y@coords) < minsize) {
          return(NULL)
        } else {
          return(poly_build(y))
        }
      })
      return(polylist[!unlist(lapply(polylist, is.null))])
    })
  })
  
  names(slicelist) <- planes
  
  dflist <- lapply(names(slicelist), function(s) {
    df <- as.data.frame(do.call(rbind, lapply(slicelist[[s]], function(x) do.call(rbind, x))))
    df$axis <- paste0(s, " slice ", unique(df$slice))
    return(df)
  })
  
  names(dflist) <- names(slicelist)
  
  if (show_labels) {
    # Centers for the sagittal plane
    centerlist <- lapply(names(dflist), function(d) {
      centersdf <- as.data.frame(do.call(
        rbind,
        lapply(
          unique(dflist[[d]]$structure),
          function(x) {
            chosen_poly <- dflist[[d]][dflist[[d]]$structure == x, ]
            chosen_poly <- chosen_poly[chosen_poly$id == names(which.max(table(chosen_poly$id))), ]
            return(unlist(poi(chosen_poly[, 1:2], precision = 0.01))[1:2])
          }
        )
      ))
      centersdf$axis <- unique(dflist[[d]]$axis)
      centersdf$acronym <- segmentation@ontology[unique(dflist[[d]]$structure), "acronym"]
      return(centersdf)
    })
    
    centers_axes <- do.call(rbind, centerlist)
  }
  
  # Smoothing
  if (smooth) {
    slicelist <- lapply(dflist, function(x) poly_smooth(x, smoothness = smoothness))
  } else {
    slicelist <- dflist
  }
  
  # Create the data frame containing all slices and their (selected) polygons
  
  all_str_polys_axes <- do.call(rbind, slicelist)
  all_str_polys_axes$structure <- all_str_polys_axes$structure
  
  all_str_polys_axes$col <- ontology(segmentation)[as.character(all_str_polys_axes$structure), "col"]
  
  all_str_polys_axes$acronym <- ontology(segmentation)[as.character(all_str_polys_axes$structure), "acronym"]
  
  # Color palette wrangling - several polygons have the same color, so this is necessary
  cols <- sapply(
    unique(all_str_polys_axes$acronym),
    function(x) unique(all_str_polys_axes[all_str_polys_axes$acronym == x, "col"])
  )
  
  cols <- as.character(cols[levels(factor(all_str_polys_axes$acronym))])
  
  # Plot
  p <- ggplot() +
    geom_polygon(
      data = all_str_polys_axes,
      aes_string(x = "x", y = "y", group = "id", subgroup = "subid", fill = "acronym"),
      color = "black"
    ) +
    scale_fill_manual(values = c(cols, "white")) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    facet_wrap(~axis, nrow = wrap_options)
  
  # Axis rulers
  if (is.null(c_slice) | is.null(s_slice) | is.null(a_slice)) show_axis_rulers <- FALSE
  if (show_axis_rulers) {
    
    # We show how the coronal and axial planes are sliced on the sagittal plot;
    # how the sagittal and axial planes are sliced on the coronal plot;
    # how the coronal and sagittal planes are sliced on the axial plot
    
    rulers <- data.frame(
      axis = c(
        unique(dflist$sagittal$axis),
        unique(dflist$coronal$axis),
        unique(dflist$axial$axis)
      ),
      xi = c(
        as.numeric(names(segmentation@slices$coronal)[c_slice]),
        as.numeric(names(segmentation@slices$sagittal)[s_slice]),
        as.numeric(names(segmentation@slices$sagittal)[s_slice])
      ),
      yi = c(
        as.numeric(names(segmentation@slices$axial)[a_slice]),
        as.numeric(names(segmentation@slices$axial)[a_slice]),
        as.numeric(names(segmentation@slices$coronal)[c_slice])
      )
    )
    
    p <- p +
      geom_vline(
        data = rulers, aes_string(xintercept = "xi"),
        linetype = "dashed",
        color = "gray",
        size = 0.5
      ) +
      geom_hline(
        data = rulers, aes_string(yintercept = "yi"),
        linetype = "dashed",
        color = "gray",
        size = 0.5
      )
  }
  
  # Include outline
  if (show_outline) {
    outlines <- lapply(planes, function(x) {
      if (smooth) {
        outline_poly <- poly_smooth(segmentation@outlines[[x]],
                                    by = "cluster",
                                    smoothness = smoothness
        )
      } else {
        outline_poly <- segmentation@outlines[[x]]
      }
      outline_poly$axis <- unique(dflist[[x]]$axis)
      return(outline_poly)
    })
    all_outlines <- do.call(rbind, outlines)
    p <- p + geom_polygon(
      data = all_outlines,
      aes_string(x = "x", y = "y", group = "cluster"),
      col = "black",
      fill = "NA"
    )
  }
  
  # Adding labels
  
  if (show_labels) {
    p <- p + geom_text_repel(
      data = centers_axes,
      aes_string(x = "x", y = "y", label = "acronym"),
      color = "white",
      segment.color = "black",
      segment.size = 0.3,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = 0,
      size = label_size,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
  }
  return(p + coord_fixed())
}


#' Plot numeric data on a segmentation
#'
#' Creates a `ggplot` plot where numeric data, mapped to structures, fill the segmentation shapes in a colour scale
#'
#' @param segmentation a \code{segmentation} class object.
#' @param feature character, feature to be plotted from the \code{assay} slot.
#' @param assay character, the name of the \code{assay} slot in the \code{segmentation}.
#' @param projection Character, name of the maximum projection slot to be used? Default is \code{NULL}, which uses the same name as selected in \code{assay}.
#' @param plane character, the plane for the maximum projection. Default is "sagittal".
#' @param by character with two elements, [1] the name of the column in the \code{sampledata} of the \code{assay} to be used for subsetting and aggregating and [2] the value to which the column in [1] should be equal to
#' @param aggr_fun function, the function to be used to aggregate samples across rows when using \code{by}
#' @param rng numeric with two elements with lower and upper bound for the color scale. Default is \code{NULL}
#' @param smooth logical, should shapes be smoothed? Default is \code{TRUE}.
#' @param smoothness numeric, the smoothing to be used. Default is 3.
#' @param minsize numeric, minimum number of vertices to draw a polygon. Default is 10.
#' @param color_pal character, the color palette to be used. The default is the `Sunset` palette from \code{colorspace}
#' @param show_side character, the side of the projection to be plotted. one of "first" ,"second", or "both". Default is "both".
#' @param show_labels logical, should segmentation labels be shown? Default is \code{TRUE}
#' @param labelsize numeric, the size of the labels. Default is 2.
#' @param remove_axes logical, should axes be shown? Default is \code{FALSE}
#'
#' @return a `ggplot` plot in which structures are coloured according to a numeric value
#'
#' @importFrom polylabelr poi
#' @importFrom ggplot2 ggplot geom_polygon aes_string theme_bw
#'    scale_fill_gradientn labs theme element_blank ggtitle
#'    coord_fixed scale_x_reverse
#' @importFrom gridExtra grid.arrange
#' @importFrom ggrepel geom_text_repel
#' @importFrom utils installed.packages
#'
#' @export

seg_feature_plot <- function(segmentation, feature, assay, projection = NULL,
                             plane = "sagittal", by = NULL, aggr_fun = NULL,
                             rng = NULL, smooth = TRUE, smoothness = 3,
                             minsize = 10, color_pal = NULL, show_side = "both",
                             show_labels = TRUE, labelsize = 2,
                             remove_axes = TRUE) {
  if (smooth & !"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")
  if (class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")
  
  if (is.null(projection)) {
    if (length(segmentation@projections) == 0) stop("The segmentation must include a projection to plot assay data. Run `seg_projection_add()` first.")
    projection <- assay
  } else {
    projection <- projection
  }
  if (length(feature) > 1) stop("You can only one plot one feature at a time.")
  if (class(feature) != "character") stop("feature must be a character.")
  if (!feature %in% rownames(segmentation@assays[[assay]]@values)) stop(paste0("Feature ", feature, " cannot be found in the row names of assay ", assay, "."))
  
  if (!is.null(by) & length(by) != 2) stop("Argument `by` must contain exactly 2 elements: name of the column in `sampledata` and value of the column")
  if (!is.null(by) & is.null(aggr_fun)) stop("If `by` is not NULL you must specify an aggregation function for `aggr_fun`.")
  
  if (is.null(color_pal)) cpal <- colorspace::sequential_hcl(palette = "Sunset", n = 25) else cpal <- color_pal
  
  selected <- segmentation@projections[[projection]][[plane]][[1]][which(lapply(segmentation@projections[[projection]][[plane]][[1]], function(x) nrow(x@coords)) > minsize)]
  
  proj_1 <- do.call(rbind, lapply(selected, poly_build))
  
  centers <- as.data.frame(do.call(rbind, lapply(unique(proj_1$structure), function(x) {
    df <- proj_1[proj_1$structure == x, ]
    df <- df[df$id == df$id[which.max(table(df$id))], ]
    return(unlist(poi(df[, 1:2], precision = 0.01))[1:2])
  })))
  
  centers$structure <- unique(proj_1$structure)
  centers$acronym <- ontology(segmentation)[as.character(centers$structure), "acronym"]
  centers$col <- ontology(segmentation)[as.character(centers$structure), "col"]
  
  if (smooth) proj_1 <- poly_smooth(proj_1, by = "subid", smoothness = smoothness)
  
  proj_1$dir <- names(projections(segmentation, projection)[[plane]])[1]
  proj_1$acronym <- ontology(segmentation)[as.character(proj_1$structure), "acronym"]
  centers$dir <- unique(proj_1$dir)
  
  selected <- segmentation@projections[[projection]][[plane]][[2]][which(lapply(segmentation@projections[[projection]][[plane]][[2]], function(x) nrow(x@coords)) > minsize)]
  
  proj_2 <- do.call(rbind, lapply(selected, poly_build))
  
  centers2 <- as.data.frame(do.call(rbind, lapply(unique(proj_2$structure), function(x) {
    df <- proj_2[proj_2$structure == x, ]
    df <- df[df$id == df$id[which.max(table(df$id))], ]
    return(unlist(poi(df[, 1:2], precision = 0.01))[1:2])
  })))
  
  centers2$structure <- unique(proj_2$structure)
  centers2$acronym <- ontology(segmentation)[as.character(centers2$structure), "acronym"]
  centers2$col <- ontology(segmentation)[as.character(centers2$structure), "col"]
  
  if (smooth) proj_2 <- poly_smooth(proj_2, by = "subid", smoothness = smoothness)
  proj_2$dir <- names(projections(segmentation, projection)[[plane]])[2]
  proj_2$acronym <- ontology(segmentation)[as.character(proj_2$structure), "acronym"]
  centers2$dir <- unique(proj_2$dir)
  
  if (!is.null(by) & !is.null(aggr_fun)) {
    column_by <- by[1]
    value_by <- by[2]
    
    cdata <- segmentation@assays[[assay]]@sampledata
    samples_keep <- rownames(cdata)[cdata[, column_by] == value_by]
    cdata_keep <- cdata[samples_keep, ]
    
    values_keep <- segmentation@assays[[assay]]@values[feature, samples_keep, drop = FALSE]
    
    if (length(samples_keep) > 1) {
      values_agg_by <- do.call(cbind, lapply(unique(cdata_keep$structure_acronym), function(y) {
        samples_aggregate <- rownames(cdata_keep)[cdata_keep$structure_acronym == y]
        if (length(samples_aggregate) > 1) {
          values_aggregate <- apply(values_keep[, samples_aggregate, drop = FALSE], 1, aggr_fun)
        } else {
          values_aggregate <- values_keep[, samples_aggregate, drop = FALSE]
        }
        return(values_aggregate)
      }))
    } else {
      values_agg_by <- values_keep
    }
    
    colnames(values_agg_by) <- unique(cdata_keep$structure_acronym)
    
    values_plot <- values_agg_by
  } else {
    values_plot <- segmentation@assays[[assay]]@values[feature, , drop = FALSE]
  }
  
  struct_available <- intersect(names(segmentation@assays[[assay]]@mapping), colnames(values_plot))
  
  struct_df <- data.frame(
    "structures" = rep(struct_available, lengths(segmentation@assays[[assay]]@mapping[struct_available])),
    "id" = sapply(unlist(segmentation@assays[[assay]]@mapping[struct_available]), function(x) ontology(segmentation)$id[ontology(segmentation)$acronym == x])
  )
  
  struct_df$gene_expression <- values_plot[feature, struct_df$structures]
  
  struct_df <- struct_df[struct_df$id %in% seg_metadata(segmentation)$structures, , drop = FALSE]
  
  rownames(struct_df) <- struct_df$id
  
  proj_1$gene_exp <- as.numeric(struct_df[proj_1$structure, "gene_expression"])
  proj_2$gene_exp <- as.numeric(struct_df[proj_2$structure, "gene_expression"])
  
  p1 <- ggplot(proj_1, aes_string(x = "x", y = "y", group = "subid", fill = "gene_exp")) +
    geom_polygon(col = "black", size = 0.3)
  
  if (smooth) {
    p1 <- p1 + geom_polygon(data = poly_smooth(segmentation@outlines[[plane]], by = "cluster"), aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  } else {
    p1 <- p1 + geom_polygon(data = segmentation@outlines[[plane]], aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  }
  
  p2 <- ggplot(proj_2, aes_string(x = "x", y = "y", group = "subid", fill = "gene_exp")) +
    geom_polygon(col = "black", size = 0.3)
  
  if (smooth) {
    p2 <- p2 + geom_polygon(data = poly_smooth(segmentation@outlines[[plane]], by = "cluster"), aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  } else {
    p2 <- p2 + geom_polygon(data = segmentation@outlines[[plane]], aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  }
  
  if (is.null(by)) {
    plot_title_1 <- unique(proj_1$dir)
    plot_title_2 <- unique(proj_2$dir)
  } else {
    plot_title_1 <- paste0(unique(proj_1$dir), " - ", by[1], ": ", by[2])
    plot_title_2 <- paste0(unique(proj_2$dir), " - ", by[1], ": ", by[2])
  }
  
  if (is.null(rng)) limits <- range(struct_df$gene_expression) else limits <- rng
  
  p1 <- p1 +
    theme_bw() +
    scale_fill_gradientn(
      na.value = "white",
      colours = cpal,
      limits = limits
    ) +
    labs(fill = feature) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    
    ggtitle(plot_title_1) +
    coord_fixed()
  
  p2 <- p2 +
    theme_bw() +
    scale_fill_gradientn(
      na.value = "white",
      colours = cpal,
      limits = limits
    ) +
    labs(fill = feature) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(plot_title_2) +
    coord_fixed()
  
  if (remove_axes) {
    p1 <- p1 + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
    p2 <- p2 + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  }
  
  if (plane == "sagittal") {
    p1 <- p1 + scale_x_reverse()
  } else if (plane == "coronal") {
    p2 <- p2 + scale_x_reverse()
  } else if (plane == "axial") {
    p1 <- p1 + scale_x_reverse()
  }
  
  if (show_labels) {
    p1 <- p1 + geom_text_repel(
      data = centers,
      aes_string(x = "x", y = "y", label = "acronym"),
      color = "white",
      segment.color = "black",
      segment.size = 0.2,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = 0.6,
      size = labelsize,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
    
    p2 <- p2 + geom_text_repel(
      data = centers2,
      aes_string(x = "x", y = "y", label = "acronym"),
      color = "white",
      segment.color = "black",
      segment.size = 0.2,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = 0.6,
      size = labelsize,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
  }
  
  if (show_side == "first") {
    return(p1)
  } else if (show_side == "second") {
    return(p2)
  } else if (show_side == "both") {
    return(grid.arrange(p1, p2, ncol = 2))
  }
}

#' Plot numeric data on a segmentation from a complex assay
#'
#' Creates a `grid` plot where numeric data, mapped to structures,
#'   fill the segmentation shapes in a colour scale and a heatmap is displayed
#'   at the bottom, showing the whole range of values across an experimental
#'   variable.
#'
#' @param segmentation a \code{segmentation} class object.
#' @param feature character, feature to be plotted from the \code{assay} slot.
#' @param assay character, the name of the \code{assay} slot in the \code{segmentation}.
#' @param projection character, name of the maximum projection slot to be used? Default is \code{NULL}, which uses the same name as selected in \code{assay}.
#' @param plane character, the plane for the maximum projection. Default is "sagittal".
#' @param by character with two elements, [1] the name of the column in the \code{sampledata} of the \code{assay} to be used for subsetting and aggregating and [2] the value to which the column in [1] should be equal to
#' @param aggr_fun function, the function to be used to aggregate samples across rows when using \code{by}
#' @param scale logical, should the values in the heatmap be scaled? Default is FALSE.
#' @param smooth logical, should shapes be smoothed? Default is \code{TRUE}.
#' @param smoothness numeric, the smoothing to be used. Default is 3.
#' @param minsize numeric, minimum number of vertices to draw a polygon. Default is 10.
#' @param color_pal character, the color palette to be used. The default is the `Sunset` palette from \code{colorspace}
#' @param show_labels logical, should segmentation labels be shown? Default is \code{TRUE}
#' @param labelsize numeric, the size of the labels. Default is 2.
#'
#' @return a `grid` plot in which structures are coloured according to a numeric
#'   value in a segmentation projection and in a heatmap
#'
#' @importFrom stringr str_count
#' @importFrom ComplexHeatmap rowAnnotation Heatmap draw anno_mark
#' @importFrom grid grid.grabExpr unit gpar
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 theme theme_void element_blank
#'
#' @export

seg_feature_complex_plot <- function(segmentation, feature, assay,
                                     projection = NULL, plane = "sagittal",
                                     by = NULL, aggr_fun = mean, scale = FALSE,
                                     smooth = TRUE, smoothness = 3,
                                     minsize = 10, color_pal = NULL,
                                     show_labels = TRUE, labelsize = 2) {
  if (smooth & !"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")
  if (class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")
  
  if (is.null(projection)) {
    if (length(segmentation@projections) == 0) stop("The segmentation must include a projection to plot assay data. Run `seg_projection_add()` first.")
    projection <- assay
  } else {
    projection <- projection
  }
  if (length(feature) > 1) stop("You can only one plot one feature at a time.")
  if (class(feature) != "character") stop("feature must be a character.")
  if (!feature %in% rownames(segmentation@assays[[assay]]@values)) stop(paste0("Feature ", feature, " cannot be found in the row names of assay ", assay, "."))
  
  if (length(by) != 2) stop("Argument `by` must contain exactly 2 elements: name of the column in `sampledata` and value of the column")
  if (!is.null(by) & is.null(aggr_fun)) stop("If `by` is not NULL you must specify an aggregation function for `aggr_fun`.")
  
  if (is.null(color_pal)) cpal <- colorspace::sequential_hcl(palette = "Sunset", n = 25) else cpal <- color_pal
  
  if (is.null(plane)) plane <- names(segmentation@projections[[projection]])[1]
  
  cdata <- segmentation@assays[[assay]]@sampledata
  values <- segmentation@assays[[assay]]@values[feature, , drop = FALSE]
  
  exp_list <- lapply(unique(cdata[, by[1]]), function(x) values[, cdata$sample_id[cdata[, by[1]] == x], drop = FALSE])
  
  for (i in seq_len(length(exp_list))) {
    cdata_current <- cdata[cdata[, by[1]] == unique(cdata[, by[1]])[i], , drop = FALSE]
    exp_list[[i]] <- exp_list[[i]][, order(cdata_current[, "structure_acronym"]), drop = FALSE]
    colnames_data <- table(cdata_current[cdata_current[, by[1]] == unique(cdata_current[, by[1]]), "structure_acronym"])
    colnames_data <- paste0(rep(names(colnames_data), colnames_data), "_", unlist(lapply(colnames_data, seq_len)))
    colnames(exp_list[[i]]) <- colnames_data
  }
  
  
  n_by <- apply(table(cdata[, by[1]], cdata$structure_acronym), 2, max)
  colnames_by <- paste0(rep(names(n_by), n_by), "_", unlist(lapply(n_by, seq_len)))
  
  exp_mat <- matrix(NA, ncol = length(colnames_by), nrow = length(unique(cdata[, by[1]])))
  rownames(exp_mat) <- unique(cdata[, by[1]])
  colnames(exp_mat) <- colnames_by
  
  for (i in 1:nrow(exp_mat)) {
    exp_mat[i, colnames(exp_list[[i]])] <- as.numeric(exp_list[[i]][1, ], nrow = 1)
  }
  
  exp_agg_mat <- do.call(cbind, lapply(names(n_by), function(x) {
    rowMeans(exp_mat[, paste0(x, "_", seq_len(n_by[x]))], na.rm = TRUE)
  }))
  
  exp_agg_mat[is.nan(exp_agg_mat)] <- NA
  colnames(exp_agg_mat) <- names(n_by)
  
  if (scale) exp_agg_mat <- t(scale(t(exp_agg_mat)))
  
  label_list <- lapply(segmentation@assays[[assay]]@mapping[colnames(exp_agg_mat)], function(x) {
    paste0(strwrap(paste0(x, collapse = ", "), width = 12), collapse = "\n")
  })
  
  line_pad <- max(unlist(lapply(label_list, function(x) str_count(x, "\n")))) / 2
  
  ha <- rowAnnotation(foo = anno_mark(
    side = "left",
    extend = 1.5,
    labels_gp = gpar(cex = 0.5),
    at = seq_len(ncol(exp_agg_mat)),
    padding = unit(line_pad, "lines"),
    link_width = unit(12, "mm"),
    labels = label_list
  ))
  
  hm <- Heatmap(t(exp_agg_mat),
                row_title = "Structure acronym",
                column_title = by[1],
                heatmap_legend_param = list(title = feature),
                cluster_columns = FALSE,
                cluster_rows = FALSE, left_annotation = ha,
                col = cpal,
                na_col = "gray"
  )
  
  hm <- grid.grabExpr(draw(hm, padding = unit(c(5, 20, 20, 20), "mm")))
  
  if (scale) {
    this_rng <- NULL
  } else {
    this_rng <- range(exp_agg_mat, na.rm = TRUE)
  }
  
  fp <- seg_feature_plot(segmentation,
                         feature = feature,
                         projection = projection,
                         by = by,
                         assay = assay,
                         plane = plane,
                         smooth = smooth,
                         smoothness = smoothness,
                         color_pal = color_pal,
                         show_labels = show_labels,
                         labelsize = labelsize,
                         aggr_fun = aggr_fun,
                         show_side = "first",
                         rng = this_rng
  ) +
    theme_void()
  
  if (!scale) {
    fp <- fp +
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(2, 2, 2, 0), "mm"))
  } else {
    fp <- fp +
      theme(plot.margin = unit(c(2, 2, 2, 0), "mm"))
  }
  
  sp <- seg_feature_plot(segmentation,
                         feature = feature,
                         projection = projection,
                         by = by,
                         assay = assay,
                         plane = plane,
                         smooth = smooth,
                         smoothness = smoothness,
                         color_pal = color_pal,
                         show_labels = show_labels,
                         labelsize = labelsize,
                         aggr_fun = aggr_fun,
                         show_side = "second",
                         rng = this_rng
  ) +
    theme_void()
  
  if (!scale) {
    sp <- sp +
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(2, 2, 2, 0), "mm"))
  } else {
    sp <- sp +
      theme(plot.margin = unit(c(2, 2, 2, 0), "mm"))
  }
  
  
  grobs_toplot <- list(fp, sp, hm)
  
  grid.arrange(
    grobs = grobs_toplot,
    heights = c(1, 1, 0.2),
    layout_matrix = rbind(
      c(1, 2),
      c(3, 3),
      c(4, NA)
    )
  )
}



#' Remove a projection from a segmentation
#'
#' Removes a named projection from a segmentation
#'
#' @param segmentation a \code{segmentation} class object.
#' @param name character, name of the projection to be removed from the \code{segmentation} slot.
#'
#' @return a segmentation object without the named projection
#'
#' @export

seg_projection_remove <- function(segmentation, name) {
  if (class(segmentation) != "segmentation") stop("Must provide a segmentation class object")
  if (!name %in% names(segmentation@projections)) stop(paste0("The projection named ", name, " was not found in this segmentation."))
  segmentation@projections[name] <- NULL
  return(segmentation)
}

#' Plot a projection from a segmentation
#'
#' Plots a named projection from a segmentation using the ontology color scheme
#'
#' @param segmentation a \code{segmentation} class object.
#' @param name character, name of the projection to be plotted from the \code{segmentation} slot.
#' @param plane character, name of the plane (one of `sagittal`, `coronal`, or `axial`)
#' @param minsize numeric, minimum number of vertices to draw a polygon. Default is 10.
#' @param smooth logical, should polygons be smoothed? Default is \code{TRUE}
#' @param smoothness numeric, the kernel bandwidth for kernel smoothing. Default is 3.
#' @param show_labels logical, should structure acronyms be shown as labels? Default is \code{FALSE}
#' @param remove_axes logical, should axes be removed? Default is \code{TRUE}
#'
#' @return a `ggplot` plot of the projection in both directions
#'
#' @importFrom polylabelr poi
#' @importFrom ggplot2 ggplot geom_polygon aes_string theme_bw
#'    scale_fill_gradientn labs theme element_blank ggtitle
#'    coord_fixed scale_x_reverse
#' @importFrom gridExtra grid.arrange
#'
#' @export

seg_projection_plot <- function(segmentation, name, plane, minsize = 10,
                                smooth = TRUE, smoothness = 3,
                                show_labels = FALSE, remove_axes = TRUE) {
  if (smooth & !"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")
  if (class(segmentation) != "segmentation") stop("Must provide a segmentation class object")
  if (!name %in% names(segmentation@projections)) stop(paste0("The projection named ", name, " was not found in this segmentation."))
  if (!plane %in% names(segmentation@projections[[name]])) stop(paste0("The plane named ", plane, " was not found in this segmentation's projection named ", name, "."))
  
  selected <- segmentation@projections[[name]][[plane]][[1]][which(lapply(segmentation@projections[[name]][[plane]][[1]], function(x) nrow(x@coords)) > minsize)]
  
  proj_1 <- do.call(rbind, lapply(selected, poly_build))
  
  centers <- as.data.frame(do.call(rbind, lapply(unique(proj_1$structure), function(x) {
    df <- proj_1[proj_1$structure == x, ]
    df <- df[df$id == df$id[which.max(table(df$id))], ]
    return(unlist(poi(df[, 1:2], precision = 0.01))[1:2])
  })))
  
  centers$structure <- unique(proj_1$structure)
  centers$acronym <- ontology(segmentation)[as.character(centers$structure), "acronym"]
  centers$col <- ontology(segmentation)[as.character(centers$structure), "col"]
  
  if (smooth) proj_1 <- poly_smooth(proj_1, by = "subid", smoothness = smoothness)
  proj_1$dir <- names(projections(segmentation, name)[[plane]])[1]
  proj_1$acronym <- ontology(segmentation)[as.character(proj_1$structure), "acronym"]
  centers$dir <- unique(proj_1$dir)
  
  selected <- segmentation@projections[[name]][[plane]][[2]][which(lapply(segmentation@projections[[name]][[plane]][[2]], function(x) nrow(x@coords)) > minsize)]
  
  proj_2 <- do.call(rbind, lapply(selected, poly_build))
  
  centers2 <- as.data.frame(do.call(rbind, lapply(unique(proj_2$structure), function(x) {
    df <- proj_2[proj_2$structure == x, ]
    df <- df[df$id == df$id[which.max(table(df$id))], ]
    return(unlist(poi(df[, 1:2], precision = 0.01))[1:2])
  })))
  
  centers2$structure <- unique(proj_2$structure)
  centers2$acronym <- ontology(segmentation)[as.character(centers2$structure), "acronym"]
  centers2$col <- ontology(segmentation)[as.character(centers2$structure), "col"]
  
  if (smooth) proj_2 <- poly_smooth(proj_2, by = "subid", smoothness = smoothness)
  proj_2$dir <- names(projections(segmentation, name)[[plane]])[2]
  proj_2$acronym <- ontology(segmentation)[as.character(proj_2$structure), "acronym"]
  centers2$dir <- unique(proj_2$dir)
  
  proj_all <- rbind(proj_1, proj_2)
  
  centers_all <- rbind(centers, centers2)
  
  cols_1 <- sapply(
    unique(centers$acronym),
    function(x) unique(centers[centers$acronym == x, "col"])
  )
  
  cols_1 <- as.character(cols_1[levels(factor(centers$acronym))])
  
  cols_2 <- sapply(
    unique(centers2$acronym),
    function(x) unique(centers2[centers2$acronym == x, "col"])
  )
  
  cols_2 <- as.character(cols_2[levels(factor(centers2$acronym))])
  
  p1 <- ggplot(proj_1, aes_string(x = "x", y = "y", group = "subid", fill = "acronym")) +
    geom_polygon(col = "black", size = 0.3)
  
  if (smooth) {
    p1 <- p1 + geom_polygon(data = poly_smooth(segmentation@outlines[[plane]], by = "cluster"), aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  } else {
    p1 <- p1 + geom_polygon(data = segmentation@outlines[[plane]], aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  }
  
  p1 <- p1 + theme_bw() +
    scale_fill_manual(values = c(cols_1, "white")) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(unique(proj_1$dir)) +
    coord_fixed()
  
  p2 <- ggplot(proj_2, aes_string(x = "x", y = "y", group = "subid", fill = "acronym")) +
    geom_polygon(col = "black", size = 0.3)
  
  if (smooth) {
    p2 <- p2 + geom_polygon(data = poly_smooth(segmentation@outlines[[plane]], by = "cluster"), aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  } else {
    p2 <- p2 + geom_polygon(data = segmentation@outlines[[plane]], aes_string(x = "x", y = "y", group = "cluster"), inherit.aes = FALSE, col = "black", fill = "NA")
  }
  
  p2 <- p2 + theme_bw() +
    scale_fill_manual(values = c(cols_2, "white")) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    ggtitle(unique(proj_2$dir)) +
    coord_fixed()
  
  if (remove_axes) {
    p1 <- p1 + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
    p2 <- p2 + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  }
  
  if (plane == "sagittal") {
    p1 <- p1 + scale_x_reverse()
  } else if (plane == "coronal") {
    p2 <- p2 + scale_x_reverse()
  } else if (plane == "axial") {
    p1 <- p1 + scale_x_reverse()
  }
  
  if (show_labels) {
    p1 <- p1 + geom_text_repel(
      data = centers,
      aes_string(x = "x", y = "y", label = "acronym"),
      color = "white",
      segment.color = "black",
      segment.size = 0.2,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = 0.6,
      size = 2,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
    
    p2 <- p2 + geom_text_repel(
      data = centers2,
      aes_string(x = "x", y = "y", label = "acronym"),
      color = "white",
      segment.color = "black",
      segment.size = 0.2,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = 0.6,
      size = 2,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
  }
  
  return(grid.arrange(p1, p2, ncol = 2))
}
