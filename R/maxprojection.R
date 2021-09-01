#' Draw segmentation
#'
#' Creates a `segmentation` object from an array or a `NIfTI` file
#'
#' @param nifti_file a character string pointing to the NIfTI file address
#' @param ontology_file a character string pointing to the ontology .csv file address
#' @param array a 3D array containing structure annotations (optional if \code{nifti_file} and \code{molten_array} are not specified). Default is \code{NULL}
#' @param molten_array an array in molten form, such as the result of \code{reshape2::melt()} (optional if \code{nifti_file} and \code{array} are not specified). Default is \code{NULL}
#' @param outliers a numeric vector indicating outlier points to be eliminated before drawing polygons. Default is \code{NULL}
#' @param verbose logical, should messages on the process be displayed? Default is \code{TRUE}
#' @param directions_from character, indicates the original orientation of the array. Default is \code{"RAS"} (Left to Right, Posterior to Anterior, Inferior to Superior).
#' @param directions_to character, indicates the final orientation of the array. The array is rotated only if it is different from \code{directions_from}. Default is \code{"RAS"}, which results in the proper assignment of axis labels.
#' @param planes character, either \code{"all"} (default) or any subset of \code{"sagittal"}, \code{"coronal"}, and \code{"axial"}. Determines along which axes the array wil be sliced.
#' @param draw_outline logical, should the outline of the whole array be drawn? Default is `TRUE`
#' @param subset_sagittal numeric vector, indicates a subset of slices along the sagittal plane that will be drawn. Default is \code{NULL}
#' @param subset_coronal numeric vector, indicates a subset of slices along the coronal plane that will be drawn. Default is \code{NULL}
#' @param subset_axial numeric vector, indicates a subset of slices along the axial plane that will be drawn. Default is \code{NULL}
#' @param reference_space character, what reference space does this segmentation use, e.g. \code{"MNI152"}, \code{"CCFv3"}, \code{"JFRC2010"}. Default is \code{NULL}
#' @param citation list of character strings, what is the citation for this segmentation? Default is \code{NULL}.
#' @param parallel logical, should the polygons be drawn using parallel processing? Uses \code{future.apply::future_lapply()}, so the proper setup must be done beforehand. Default is \code{FALSE}.
#'
#' @return a `segmentation` S4 object containing the following slots:
#'
#' @export

drawSegmentation <- function(nifti_file = NULL,
                             ontology_file,
                             array = NULL,
                             molten_array = NULL,
                             outliers = NULL,
                             verbose = TRUE,
                             directions_from = "RAS",
                             directions_to = "RAS",
                             planes = "all",
                             subset_structures = NULL,
                             draw_outline = TRUE,
                             subset_sagittal = NULL,
                             subset_coronal = NULL,
                             subset_axial = NULL,
                             reference_space = NULL,
                             citation = NULL,
                             parallel = TRUE) {
  #Sanity checks
  if (is.null(nifti_file) & is.null(array) & is.null(molten_array)) stop("Must provide at least one 3D array or one NIfTI file.")
  if (planes != "all" & !any(planes %in% c("sagittal", "coronal", "axial"))) stop("You must choose at least one plane among sagittal, coronal or axial.")
  if (!is.null(subset_sagittal) & class(subset_sagittal) != "numeric") stop("You must provide numeric indices to subset planes")
  if (!is.null(subset_coronal) & class(subset_coronal) != "numeric") stop("You must provide numeric indices to subset planes")
  if (!is.null(subset_axial) & class(subset_axial) != "numeric") stop("You must provide numeric indices to subset planes")

  if (verbose) cat("Adding ontology...")
  ontology <- read.csv(ontology_file, header = TRUE, sep = ",")
  colnames(ontology)[grep("id|ID|Id|iD", fixed = TRUE, colnames(ontology))] <- "id"
  ontology$id <- as.character(ontology$id)
  rownames(ontology) <- ontology$id
  if(any(subset_structures %nin% ontology$id)) {
    stop(paste0("Structures ", paste(setdiff(subset_structures, ontology$id), collapse = ", "), " were not found in the ontology. Make sure you are using the correct ontology for this array." ))
  }
  if (verbose) cat("done.\n")

  if (!is.null(subset_sagittal) | !is.null(subset_coronal) | !is.null(subset_axial)) {
    subsets <- list(
      "sagittal" = subset_sagittal,
      "coronal" = subset_coronal,
      "axial" = subset_axial
    )
    subsets <- subsets[!is.null(subsets)]
  }

  if (!is.null(nifti_file) & is.null(array) & is.null(molten_array)) {
    if (verbose) cat("Reading NIfTI file...")
    nifti <- oro.nifti::readNIfTI(nifti_file)
    if (verbose) cat("done.\n")
    n_image <- oro.nifti::img_data(nifti)
    pixdims <- oro.nifti::pixdim(nifti)[2:4]
    units <- oro.nifti::xyzt_units(nifti)
    filename = nifti_file
  } else if (is.null(nifti_file) & !is.null(array) & is.null(molten_array)) {
    n_image <- array
    pixdims <- "Unknown"
    units <- "Unknown"
    filename = paste0("Array named \"", deparse(substitute(molten_array)), "\".")
  } else if (is.null(nifti_file) & is.null(array) & !is.null(molten_array)) {
    pixdims <- "Unknown"
    units <- "Unknown"
    filename = paste0("Molten array named \"", deparse(substitute(molten_array)), "\".")
  }

  if (directions_from != directions_to) {
    dirs = c(directions_from, directions_to)
    if (verbose) cat("Changing directions...")
    n_image <- changeDirections(n_image, directions_from, directions_to)
    if (verbose) cat("done.\n")
  } else {
    dirs = paste0(directions_from, " (unchanged)")
  }

  if (!is.null(outliers)) n_image[t(outliers)] <- 0
  ndims <- dim(n_image)

  if (planes == "all") planes_chosen <- c("sagittal", "coronal", "axial") else planes_chosen <- planes
  slices <- outlines <- list("sagittal" = NULL, "coronal" = NULL, "axial" = NULL)

  if (!is.null(array) & is.null(molten_array) & (is.null(subset_sagittal) & is.null(subset_coronal) & is.null(subset_axial))) {
    if (verbose) cat("Melting array...")
    M <- as.data.table(n_image)
    M <- M[M$value > 0, ]
    if (verbose) cat("done.\n")
  } else if (is.null(array) & !is.null(molten_array)) {
    if (any(M$value == 0)) M <- M[M$value > 0, ]
    M <- molten_array
  } else if (is.null(array) & is.null(molten_array) & !is.null(nifti_file)) {
    if (verbose) cat("Melting array...")
    M <- as.data.table(n_image)
    M <- M[M$value > 0, ]
    if (verbose) cat("done.\n")
  }

  if(!is.null(subset_structures)) {
    if (verbose) cat("Subsetting structures...")
    M <- M[M$value %in% subset_structures,]
    if (verbose) cat("done.\n")
  }

  if(any(unique(M$value) %nin% ontology$id)) {
    stop(paste0("Structures ", paste(setdiff(unique(M$value), ontology$id), collapse = ", "), " were found in the array but not in the ontology. Make sure you are using the correct ontology for this array." ))
  }

  if (draw_outline) {
    for (i in planes_chosen) {
      if (verbose) cat("Drawing", i, "outline...")
      ol <- drawOutline(M, plane = i)
      outlines[[i]] <- ol
      if (verbose) cat("done.\n")
    }
  }

  if (!is.null(subset_sagittal) | !is.null(subset_coronal) | !is.null(subset_axial)) {
    if (verbose) cat("Subsetting slices...")

    M_s <- reshape2::melt(n_image[subset_sagittal, , ])
    M_s$V1 <- subset_sagittal
    M_s <- M_s[, c(4, 1, 2, 3)]
    M_s <- as.data.table(M_s[M_s$value > 0, ])

    M_c <- reshape2::melt(n_image[, subset_coronal, ])
    M_c$V2 <- subset_coronal
    M_c <- M_c[, c(1, 4, 2, 3)]
    M_c <- as.data.table(M_c[M_c$value > 0, ])

    M_a <- reshape2::melt(n_image[, , subset_axial])
    M_a$V3 <- subset_axial
    M_a <- M_a[, c(1, 2, 4, 3)]
    M_a <- as.data.table(M_a[M_a$value > 0, ])

    colnames(M_s) <- colnames(M_c) <- colnames(M_a) <- c("V1", "V2", "V3", "value")
    if (verbose) cat("done.\n")
  }

  if (is.null(subset_sagittal) & is.null(subset_coronal) & is.null(subset_axial)) {
    structure_ids <- unique(M$value)
  } else {
    structure_ids <- union(union(unique(M_s$value), unique(M_c$value)), unique(M_a$value))
  }

  if (is.null(subset_sagittal) & is.null(subset_coronal) & is.null(subset_axial)) {
    for (i in planes_chosen) {
      if (verbose) cat("Slicing along the", i, "plane...")
      sl <- makeSlices(M, plane = i)
      if (verbose) cat("done.\n")

      if (verbose) cat("Drawing polygons...")
      sp <- makePolygonSets(sl, parallel = parallel, verbose = verbose)
      slices[[i]] <- sp
      if (verbose) cat("done.\n")
    }
  } else if (!is.null(subset_sagittal) | !is.null(subset_coronal) | !is.null(subset_axial)) {
    sl_s <- makeSlices(M_s, plane = "sagittal")
    sl_c <- makeSlices(M_c, plane = "coronal")
    sl_a <- makeSlices(M_a, plane = "axial")

    if (verbose) cat("Drawing polygons...")
    sp_s <- makePolygonSets(sl_s, parallel = parallel)
    slices[["sagittal"]] <- sp_s
    sp_c <- makePolygonSets(sl_c, parallel = parallel)
    slices[["coronal"]] <- sp_c
    sp_a <- makePolygonSets(sl_a, parallel = parallel)
    slices[["axial"]] <- sp_a
    if (verbose) cat("done.\n")
  }

  dims_eff = lengths(slices)

  if (verbose) cat("Compiling structure tables...")
  str_table_all <- lapply(planes_chosen, function(n) {
    structure_slices <- as.data.table(cbind(unlist(lapply(unlist(slices[[n]]), function(x) x@structure)),
                                            unlist(unlist(lapply(unlist(slices[[n]]), function(x) x@slice)))))
    structure_slices <- structure_slices[!duplicated(structure_slices[,1:2]),]
    colnames(structure_slices) <- c("structure", "slice")
    return(structure_slices)
  })
  names(str_table_all) <- planes_chosen
  for(i in names(slices)) {
         for(j in names(slices[[i]])) {
              names(slices[[i]][[j]]) <- str_table_all[[i]][str_table_all[[i]]$slice == j, structure]
           }
     }
  if (verbose) cat("done.\n")

  if (verbose) cat("Compiling metadata...")
  metadata <- list(
    "filename" = filename,
    "dirs" = dirs,
    "planes" = planes_chosen,
    "dims_original" = ndims,
    "dims_effective" = dims_eff,
    "pixdim" = pixdims,
    "units" = units,
    "reference_space" = reference_space,
    "structures" = structure_ids,
    "citation" = citation
  )
  if (verbose) cat("done.\n")

  return(new("segmentation",
    slices = slices,
    ontology = ontology,
    outlines = outlines,
    metadata = metadata,
    structure_tables = str_table_all,
    projections = list()
  ))
}


changeDirections <- function(array,
                             from = "PIR",
                             to = "RAS") {
  if (from == "PIR" & to == "RAS") {
    array_rot <- freesurferformats::rotate3D(freesurferformats::rotate3D(array, axis = 1, degrees = 90), axis = 3, degrees = 90)
  }
  if (from == "RAS" & to == "PIR") {
    array_rot <- freesurferformats::rotate3D(freesurferformats::rotate3D(array, axis = 3, degrees = 270), axis = 1, degrees = 270)
  }
  if (from == "LPS" & to == "RAS") {
    array_rot <- freesurferformats::rotate3D(array, axis = 3, degrees = 180)
  }
  if (from == "RAS" & to == "LPS") {
    array_rot <- freesurferformats::rotate3D(array, axis = 3, degrees = 180)
  }
  if (from == "LPS" & to == "PIR") {
    array_rot <- freesurferformats::rotate3D(freesurferformats::rotate3D(array, axis = 1, degrees = 270), axis = 3, degrees = 90)
  }
  if (from == "PIR" & to == "LPS") {
    array_rot <- freesurferformats::rotate3D(freesurferformats::rotate3D(array, axis = 3, degrees = 90), axis = 1, degrees = 270)
  }
  return(array_rot)
}

#' Get visual center
#'
#' Finds the visual center (pole of inaccessibility) for a polygon.
#'
#' @param poly a polygon in array or matrix form, in which the first two columns contain x and y vertex coordinates
#'
#' @return a vector containing the coordinates of the visual center of the input polygon

getCenter <- function(poly) {
  return(unlist(polylabelr::poi(poly[, 1:2]))[1:2])
}

#' Check slices
#'
#' Checks which slices in a polygon list contain a given structure id
#'
#' @param structures a character containing one or more structure ids
#' @param segmentation an object of class \code{segmentation}
#' @param planes a vector indicating the planes in which the structure is looked up
#'
#' @return a vector of slices in each plane containing the structure(s) in \code{str}.
#'
#' @export

sliceCheck <- function(structures,
                       segmentation,
                       planes) {
  structures <- as.character(structures)
  if(any(structures %nin% unique(do.call(rbind, segmentation@structure_tables)$structure))) {
    not_found <- setdiff(structures, unique(do.call(rbind, segmentation@structure_tables)$structure))
    stop("Structure(s) ", paste(not_found, collapse = ", "), " not found in this segmentation.
  Perhaps there is a typo? Check `metaData(segmentation)` and/or `segmentation@structure_tables` to see the available structures.")
  }

  structure_by_plane <- lapply(planes, function(x) {
      unique(segmentation@structure_tables[[x]][structure %in% structures,slice])
    })
  names(structure_by_plane) <- planes
  return(structure_by_plane)
}

#' Smooth polygon sets
#'
#' Uses kernel smoothing to smooth polygon sets
#'
#' @param polygon_set a data frame containing polygons, separated by `subid`
#' @param smoothness a numeric indicating the smoothing parameter, passed to `smooth_ksmooth`
#' @param min_points a numeric indicating the minimum number of points to smooth. If a polygon has less than this number of points, it is not smoothed.
#'
#' @return a polygon set with smoothed coordinates
#'
#' @export

smoothPolygons <- function(polygon_set,
                           by = "subid",
                           smoothness = 3,
                           min_points = 5) {
  smoothing_list <- lapply(unique(polygon_set[, by]), function(x) {
    poly <- as.data.frame(polygon_set[polygon_set[, by] == x, ])

    if (nrow(poly) < min_points) {
      poly_sm <- poly[, 1:2]
    } else {
      poly_sm <- as.data.frame(smoothr::smooth_ksmooth(as.matrix(poly[, 1:2]),
        smoothness = smoothness,
        wrap = TRUE
      ))
    }

    # Restore/add original columns
    colnames(poly_sm) <- c("x", "y")

    for (i in colnames(poly)[3:ncol(poly)]) {
      poly_sm[, i] <- unique(poly[, i])
    }

    return(poly_sm)
  })
  smooth_df <- as.data.frame(do.call(rbind, smoothing_list))
  return(smooth_df)
}

#' Not in
#'
#' The opposite of `%in%`
#'
#' @return a logical value indicating which elements are NOT included in the input

`%nin%` <- Negate(`%in%`)

#' Draw maximum outline
#'
#' Draws a polygon obtained by projecting all polygons on the plane
#'
#' @param M a molten array (of class `data.table`) containing coordinates
#'
#' @return a data frame containing ordered coordinates for the outline polygon
#'
#' @export

drawOutline <- function(M,
  plane) {
  switch(plane,
         "sagittal" = {
           uc <- 1
           xc <- 2
           yc <- 3
         },
         "coronal" = {
           uc <- 2
           xc <- 1
           yc <- 3
         },
         "axial" = {
           uc <- 3
           xc <- 1
           yc <- 2
         },
         stop("Must select a `plane` out of \"sagittal\", \"coronal\", \"axial\"")
  )

  columns <- c(xc, yc)
  df <- M[, ..columns]
  df <- df[!duplicated(df[, 1:2]), ]
  colnames(df) <- c("x", "y")
  outline <- makePolygon(df)
  return(outline)
}

#' Plot the ontology graph
#'
#' Creates a `ggplot` plot of the ontology graph/tree
#'
#' @param o a data frame containing the ontology/tree
#' @param circular a logical indicating whether the layout should be circular
#'
#' @return a `ggplot` object with the graph
#'
#' @export

makeOntologyGraph <- function(o,
                              circular = TRUE) {
  cols <- o$col
  names(cols) <- o$acronym

  onto_graph <- data.frame("from" = o[as.character(o$parent_structure_id), "acronym"], "to" = o$acronym)

  if (any(is.na(onto_graph[, 1]))) {
    onto_graph[which(is.na(onto_graph[, 1])), 1] <- onto_graph[which(is.na(onto_graph[, 1])), 2]
  } else if (any(is.na(onto_graph[, 2]))) {
    onto_graph[which(is.na(onto_graph[, 2])), 2] <- onto_graph[which(is.na(onto_graph[, 2])), 1]
  }

  g <- igraph::graph_from_data_frame(onto_graph)
  g <- tidygraph::as_tbl_graph(g)


  p <- ggraph(g, "dendrogram", circular = circular) +
    geom_edge_elbow() +
    geom_node_point(aes(fill = name), shape = 21, color = "black", size = 8) +
    geom_node_text(aes(label = name), color = "black", size = 2) +
    scale_fill_manual(values = cols) +
    theme(legend.position = "none")

  return(p)
}

plotSegmentation <- function(segmentation,
                             s_slice = NULL,
                             c_slice = NULL,
                             a_slice = NULL,
                             smooth = TRUE,
                             smoothness = 3,
                             show_axis_rulers = TRUE,
                             show_outline = TRUE,
                             show_labels = FALSE,
                             label_size = 2,
                             grep_structures = NULL,
                             minsize = 10,
                             # autoslice = TRUE,
                             invert = FALSE,
                             wrap_options = 1) {

  planes = c("sagittal", "coronal", "axial")
  if(is.null(s_slice)) s_slice <- NULL
  if(is.null(c_slice)) c_slice <- NULL
  if(is.null(a_slice)) a_slice <- NULL
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
          return(buildPolygon(y))
        }
      })
      return(polylist[!unlist(lapply(polylist, is.null))])
    })
  })

  names(slicelist) <- planes

  dflist <- lapply(names(slicelist), function(s) {
    df <-  as.data.frame(do.call(rbind, lapply(slicelist[[s]], function(x) do.call(rbind, x))))
    df$axis <- paste0(s, " slice ", unique(df$slice))
    return(df)
  })

  names(dflist) <- names(slicelist)

  if (show_labels) {
    # Centers for the sagittal plane
    centerlist <- lapply(names(dflist), function(d) {
      centersdf<- as.data.frame(do.call(
        rbind,
        lapply(
          unique(dflist[[d]]$structure),
          function(x) {
            chosen_poly <- dflist[[d]][dflist[[d]]$structure == x, ]
            chosen_poly <- chosen_poly[chosen_poly$id == names(which.max(table(chosen_poly$id))), ]
            return(unlist(polylabelr::poi(chosen_poly[, 1:2]))[1:2])
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
  if (smooth) slicelist <- lapply(dflist, function(x) smoothPolygons(x, smoothness = smoothness))

  # Create the data frame containing all slices and their (selected) polygons

  all_str_polys_axes <- do.call(rbind, slicelist)
  all_str_polys_axes$structure <- all_str_polys_axes$structure

  all_str_polys_axes$col <- segmentation@ontology[as.character(all_str_polys_axes$structure), "col"]

  all_str_polys_axes$acronym <- segmentation@ontology[as.character(all_str_polys_axes$structure), "acronym"]

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
      aes(x = x, y = y, group = id, subgroup = subid, fill = acronym),
      color = "black"
    ) +
    scale_fill_manual(values = c(cols, "white")) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    facet_wrap(~axis, nrow = wrap_options)

  # Axis rulers
  if(is.null(c_slice) | is.null(s_slice) | is.null(a_slice)) show_axis_rulers = FALSE
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
        data = rulers, aes(xintercept = xi),
        linetype = "dashed",
        color = "gray",
        size = 0.5
      ) +
      geom_hline(
        data = rulers, aes(yintercept = yi),
        linetype = "dashed",
        color = "gray",
        size = 0.5
      )
  }

  # Include outline
  if (show_outline) {

    outlines <- lapply(planes, function(x) {
      if (smooth) { outline_poly <- smoothPolygons(segmentation@outlines[[x]],
                                               by = "cluster",
                                               smoothness = smoothness)
        } else { outline_poly = segmentation@outlines[[x]] }
        outline_poly$axis <- unique(dflist[[x]]$axis)
        return(outline_poly)
        })
    all_outlines <- do.call(rbind, outlines)
    p <- p + geom_polygon(
      data = all_outlines,
      aes(x = x, y = y, group = cluster),
      col = "black",
      fill = "NA"
    )
  }

  # Adding labels
  if (!is.null(grep_structures)) padding <- 1 else padding <- 0

  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(
      data = centers_axes,
      aes(x = x, y = y, label = acronym),
      color = "white",
      segment.color = "black",
      segment.size = 0.3,
      bg.color = "black",
      bg.r = 0.15,
      alpha = 1,
      box.padding = padding,
      size = label_size,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
  }
  return(p  + coord_fixed())
}




#' Plot numeric data on a segmentation
#'
#' Creates a `ggplot` plot where numeric data, mapped to structures, fill the segmentation shapes in a colour scale
#'
#' @param segmentation a \code{segmentation} class object.
#' @param feature character, feature to be plotted from the \code{assay} slot.
#' @param assay character, the name of the \code{assay} slot in the \code{segmentation}.
#' @param max_projection logical, should the maximum projection be computed beforehand? Default is \code{TRUE}.
#' @param plane character, the plane for the maximum projection. Default is "sagittal".
#' @param smooth logical, should shapes be smoothed? Default is \code{TRUE}.
#' @param smoothness numeric, the smoothing to be used. Default is 3.
#' @param min_points numeric, the minimum points to smooth polygons. Default is 5.
#'
#' @return a `ggplot` plot in which structures are coloured according to a numeric value
#'
#' @export


plotBrainMap <- function(segmentation,
                         feature,
                         assay,
                         plane = "sagittal",
                         smooth = TRUE,
                         smoothness = 3,
                         min_points = 5){

  mapping = segmentation@assays[[assay]]@mapping

  included_strs <- intersect(as.character(unique(unlist(mapping))), metaData(segmentation)$structures)

    dmp_df <- do.call(rbind, lapply(segmentation@projections[[assay]][[plane]][[1]], function(x) cbind(x@coords[,1:2], "slice" = x@slice, "structure" = x@structure, "id" = x@id, "subid" = x@subid)))
    dmp2_df <- do.call(rbind, lapply(segmentation@projections[[assay]][[plane]][[2]], function(x) cbind(x@coords[,1:2], "slice" = x@slice, "structure" = x@structure, "id" = x@id, "subid" = x@subid)))

  if(smooth) {
    dmp_df <- smoothPolygons(dmp_df, smoothness = smoothness, min_points = min_points)
    dmp2_df <- smoothPolygons(dmp2_df, smoothness = smoothness, min_points = min_points)
  }

  centers <- as.data.frame(do.call(rbind, lapply(unique(dmp_df$structure), function(x) {
    df <- dmp_df[dmp_df$structure == x,]
    df <- df[df$id == df$id[which.max(table(df$id))],]
    return(unlist(polylabelr::poi(df[, 1:2]))[1:2])
  })))

  centers$structure <- unique(dmp_df$structure)
  centers$acronym <- segmentation@ontology[centers$structure, "acronym"]


  centers2 <- as.data.frame(do.call(rbind, lapply(unique(dmp2_df$structure), function(x) {
    df <- dmp2_df[dmp2_df$structure == x,]
    df <- df[df$id == df$id[which.max(table(df$id))],]
    return(unlist(polylabelr::poi(df[, 1:2]))[1:2])
  })))

  centers2$structure <- unique(dmp2_df$structure)
  centers2$acronym <- segmentation@ontology[centers2$structure, "acronym"]

  dmp_df$dir <-  centers$dir <- names(segmentation@projections[[assay]][[plane]])[1]
  dmp2_df$dir <- centers2$dir <-  names(segmentation@projections[[assay]][[plane]])[2]

  dmp_all <- rbind(dmp_df, dmp2_df)
  centers_all <- rbind(centers, centers2)

  struct_df <- data.frame("structures" = rep(names(segmentation@assays[[assay]]@mapping), lengths(segmentation@assays[[assay]]@mapping)),
                          "id" = unlist(segmentation@assays[[assay]]@mapping))

  struct_df$gene_expression <- segmentation@assays[[assay]]@values[feature, struct_df$structures]

  struct_df <- struct_df[struct_df$id %in% metaData(segmentation)$structures,]
  rownames(struct_df) <- struct_df$id

  dmp_all$gene_exp <- as.numeric(struct_df[dmp_all$structure, "gene_expression"])

  p <- ggplot(dmp_all, aes(x = x, y = y, group = id, fill = gene_exp)) +
    geom_polygon(col = "black", size = 0.3) +
    facet_wrap(~dir) +
    coord_fixed() +
    geom_polygon(data = smoothPolygons(segmentation@outlines[[plane]], by = 'cluster'), aes(x = x, y = y, group = cluster), inherit.aes = FALSE, col = "black", fill = "NA") +
    theme_bw() +
    scale_fill_gradientn(na.value = "white",
                         colours =  colorspace::sequential_hcl(palette = "Sunset", n = 25)) +
    labs(fill = paste0(feature, "\nexpression (TPM)")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  p <- p + ggrepel::geom_text_repel(
    data = centers_all,
    aes(x = x, y = y, label = acronym),
    color = "white",
    segment.color = "black",
    segment.size = 0.3,
    bg.color = "black",
    bg.r = 0.15,
    alpha = 1,
    box.padding = 0.2,
    size = 2,
    max.overlaps = Inf,
    inherit.aes = FALSE
  )

  return(p)
}
