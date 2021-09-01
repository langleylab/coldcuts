
#' Build a polygon
#'
#' Creates a polygon in a data frame from a \code{segPointSet} object
#'
#' @param seg_point_set A  \code{segPointSet} object
#'
#' @return A data frame containing polygon coordinates, \code{structure}, \code{id}, and \code{subid}
#'
#' @export

buildPolygon <- function(seg_point_set) {
  if (class(seg_point_set) != "segPointSet") stop("Must provide a `segPointSet` class object")

  return(as.data.frame(
    cbind(seg_point_set@coords[,1:2],
          "slice" = seg_point_set@slice,
          "structure" = seg_point_set@structure,
          "id" = seg_point_set@id,
          "subid" = seg_point_set@subid)))
}

#' Make slices
#'
#' Divides a molten 3D array in slices along a selected plane
#'
#' @param M a data frame molten from a 3D array in long format, with x, y and z columns as \code{"Var1"}, \code{"Var2"}, \code{"Var3"} and the value column containing the voxel annotation
#' @param plane character, one of \code{"sagittal"}, \code{"coronal"}, or \code{"axial"}. The plane along which to slice the array. Default is \code{"sagittal"}
#' @param return character, one of \code{"planes"} or \code{"structures"}. Decides whether the returning object is a list of planes, or a nested list of structures (i.e. data frames divided by unique \code{value} values). Default is \code{"structures"}.
#'
#' @return A list of planes, one per slice (if `return` is set to \code{"planes"}), or a list of lists of structures, one list per slice
#'
#' @export

makeSlices <- function(M,
                       plane = "sagittal",
                       return = "structures",
                       parallel = FALSE) {
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

  axis_plane <- data.table(split(M, by = colnames(M)[uc]))[[1]]

  columns <- c(xc, yc, uc, 4)

  for (i in seq_len(length(axis_plane))) {
    axis_plane[[i]] <- axis_plane[[i]][, ..columns]
    colnames(axis_plane[[i]]) <- c("x", "y", "slice", "structure")
  }

  # Naming slices and ordering
  names(axis_plane) <- unique(M[, ..uc])[[1]]
  axis_plane <- axis_plane[order(as.numeric(names(axis_plane)))]

  if (return == "planes") {
    return(axis_plane)
  }

  axis_strlist <- lapply(axis_plane, function(x) split(x, by = "structure"))
  # Naming structures using their codes
  for (i in seq_len(length(axis_strlist))) {
    names(axis_strlist[[i]]) <- unlist(lapply(
      axis_strlist[[i]],
      function(x) unique(x$structure)
    ))
  }
  names(axis_strlist) <- names(axis_plane)
  if (return == "structures") {
    return(axis_strlist)
  }
}

#' Create a bounding box
#'
#' Generates a rectangular bounding box around a point cloud, filling it with 0 wherever there are no points
#'
#' @param point_set a data frame containing point coordinates. Columns must be named as \code{"x"} and \code{"y"}.
#' @param step_size a numeric indicating how much bigger the bounding box should be compared to the point cloud. Default is 1 unit.
#'
#' @return A data frame containing the original points and the bounding box coordinates
#'
#' @export

boundingBox <- function(point_set,
                        step_size = 1) {

  # Make the bounding box larger than the structure by step_size (1 for non-transformed voxel arrays)
  bounding_box <- c(range(point_set[, 1]), range(point_set[, 2])) + c(-step_size, step_size, -step_size, step_size)

  # Fill the box with 0s
  bounding_box_filled <- expand.grid(
    c(seq(bounding_box[1], bounding_box[2], by = step_size)),
    c(seq(bounding_box[3], bounding_box[4], by = step_size))
  )

  colnames(bounding_box_filled) <- c("x", "y")

  # Join box and structure voxels
  bounding_box_all <- rbind(bounding_box_filled, point_set[, 1:2])

  # Only structure voxels are 1s, the rest are 0s
  bounding_box_all$value <- 1
  bounding_box_all$value[1:nrow(bounding_box_filled)] <- 0
  bounding_box_all <- bounding_box_all[-1 * which(duplicated(bounding_box_all[, 1:2], fromLast = TRUE)), ]

  return(bounding_box_all)
}



#' Fill a polygon with points
#'
#' Fills a polygon with regularly spaced points.
#'
#' @param point_set a polygon point set in array or matrix form, in which the first two columns contain x and y vertex coordinates.
#' @param step_size a numeric indicating the spacing of the points
#'
#' @return A data frame containing coordinates for points within the polygon, including polygon vertices, inheriting all the other columns from the input polygon point set.
#'
#' @export

fillPolygon <- function(point_set,
                        step_size = 1) {
  point_set <- point_set[!duplicated(point_set[, 1:2]), ]

  box <- c(range(point_set[, 1]), range(point_set[, 2]))
  if (diff(box[1:2]) <= 1 | diff(box[3:4]) <= 1) {
    return(point_set)
  }

  box <- expand.grid(
    c(seq(box[1], box[2], by = step_size)),
    c(seq(box[3], box[4], by = step_size))
  )

  colnames(box) <- c("x", "y")
  box <- box[which(splancs::inout(box, point_set, bound = NULL)), ]
  box <- rbind(box, point_set[, 1:2])
  box <- box[!duplicated(box), ]


  if (length(colnames(point_set)) > 2) {
    for (i in colnames(point_set)[3:ncol(point_set)]) {
      if (length(unique(point_set[, i])) > 1) {
        box[, i] <- NA
      } else {
        box[, i] <- unique(point_set[, i])
      }
    }
  }

  return(as.data.frame(box))
}

#' Make polygon
#'
#' Makes a polygon from a point set
#'
#' @param point_set a polygon in array or matrix form, in which the first two columns contain x and y vertex coordinates.
#' @param step_size a numeric indicating the spacing of the points
#'
#' @return A data frame containing coordinates for polygon vertices and the slice it originated from (if present in the original point set).
#'
#' @export


makePolygon <- function(point_set,
                        step_size = 1) {
  if ("x" %nin% colnames(point_set) | "y" %nin% colnames(point_set)) stop("point_set must have columns named \"x\" and \"y\"")

  m <- t(reshape2::acast(boundingBox(point_set, step_size = step_size), formula = x ~ y))

  ib <- isoband::isobands(seq_len(ncol(m)), seq_len(nrow(m)), m, levels_low = 0, levels_high = 1)

  ib_df <- data.frame(
    "x" = ib[[1]]$x,
    "y" = ib[[1]]$y,
    "cluster" = ib[[1]]$id
  )

  if ("slice" %in% colnames(point_set)) ib_df$slice <- unique(point_set$slice)

  ib_df <- ib_df[which(ib_df$x %nin% range(ib_df$x) & ib_df$y %nin% range(ib_df$y)), ]

  ib_df$x <- ib_df$x + min(point_set$x) - 2
  ib_df$y <- ib_df$y + min(point_set$y) - 2

  return(ib_df)
}

#' Subset by structures
#'
#' Subset a \code{segmentation} class object to only contain certain structures
#'
#' @param segmentation an object of class \code{segmentation}
#' @param structures a vector of character structure IDs
#' @param planes_chosen a vector of planes to be further subset. Defaults to \code{NULL}, meaning the 3 anatomical planes will be considered.
#'
#' @return A \code{segmentation} class object containing only the structures defined in \code{structures} and, optionally, the planes defined in \code{planes_chosen}. The metadata are recompiled based on the new composition of the object.
#'
#' @export


subsetByStructures <- function(segmentation,
                               structures,
                               planes_chosen = NULL){

  if(is.null(planes_chosen)) planes_chosen = metaData(segmentation)$planes

  structures <- as.character(structures)

  which_slice_ok <- sliceCheck(structures = structures,
                               segmentation = segmentation,
                               planes = planes_chosen)

  segmentation@slices <- lapply(planes_chosen, function(x) {
    segmentation@slices[[x]] <- segmentation@slices[[x]][as.character(which_slice_ok[[x]])]
    segmentation@slices[[x]] <- lapply(segmentation@slices[[x]], function(y) y[structures])
  })

  names(segmentation@slices) <- planes_chosen

  segmentation@metadata$dims_effective <- lengths(segmentation@slices)
  segmentation@metadata$dims_effective <- lengths(segmentation@slices)
  segmentation@metadata$structures <- structures

  return(segmentation)
}


#' Create a maximum projection
#'
#' Creates a point set or a polygon (ordered set of xy coordinates) resulting as the projection on the same plane of polygons in all slices.
#'
#'@param name a character containing the name of the maximum projection slot
#' @param structures a character containing one or more structure IDs
#' @param segmentation a \code{segmentation} class object
#' @param planes_chosen a character vector indicating which planes should be used to draw the projection. Default is \code{NULL}, using all planes
#' @param make_polygon logical, should the resulting polygon be created? Default is \code{TRUE}
#'
#' @return a \code{segmentation} object containing a maximum projection slot as a \code{segPointSet} object, with 2 views per selected plane.
#'
#' @export

addMaxProjection <- function(name,
                             structures,
                             segmentation,
                             planes_chosen = NULL,
                             make_polygon = TRUE) {

    if(is.null(planes_chosen)) planes_chosen = names(segmentation@slices)

    sub_segmentation = subsetByStructures(segmentation = segmentation, structures = structures, planes_chosen = planes_chosen)

    fill_list_AB <- list()

    for(i in planes_chosen) {

      all_polygons <- lapply(unlist(sub_segmentation@slices[[i]]), function(x) buildPolygon(x))

      filled_polygons <- lapply(all_polygons, fillPolygon)
      fill_df <- do.call(rbind, filled_polygons)
      fill_df <- fill_df[!duplicated(fill_df[, 1:2]), ]
      structure_projection <- fill_df[, c(1,2,4)]
      structure_projection$value <- 1
      colnames(structure_projection)[1:2] <- c("x", "y")
      if (make_polygon) ib_df <- do.call(rbind, lapply(unique(structure_projection$structure), function(x) {
      ib_temp <- makePolygon(structure_projection[structure_projection$structure == x,])
      ib_temp$structure <- x
      ib_temp$id <- paste0(ib_temp$structure, "_", ib_temp$cluster)
      return(ib_temp)
      })) else ib_df <- structure_projection

    ib_df <- ib_df[!duplicated(ib_df[,1:2]),]



    fill_list_AB[[i]] <- lapply(unique(ib_df$id), function(s)
      new("segPointSet",
          coords = ib_df[ib_df$id == s, c("x", "y")],
          slice = "maximum projection",
          structure = as.character(ib_df[ib_df$id == s, "structure"]),
          id =  as.character(s),
          subid =  as.character(s)))
    }

    fill_list_BA <- list()

    for(i in planes_chosen) {

      all_polygons <- lapply(unlist(sub_segmentation@slices[[i]]), function(x) buildPolygon(x))
      all_polygons <- all_polygons[rev(seq_len(length(all_polygons)))]
      filled_polygons <- lapply(all_polygons, fillPolygon)
      fill_df <- do.call(rbind, filled_polygons)
      fill_df <- fill_df[!duplicated(fill_df[, 1:2]), ]
      structure_projection <- fill_df[, c(1,2,4)]
      structure_projection$value <- 1
      colnames(structure_projection)[1:2] <- c("x", "y")
      if (make_polygon) ib_df <- do.call(rbind, lapply(unique(structure_projection$structure), function(x) {
        ib_temp <- makePolygon(structure_projection[structure_projection$structure == x,])
        ib_temp$structure <- x
        ib_temp$id <- paste0(ib_temp$structure, "_", ib_temp$cluster)
        return(ib_temp)
      })) else ib_df <- structure_projection

      ib_df <- ib_df[!duplicated(ib_df[,1:2]),]

      fill_list_BA[[i]] <- lapply(unique(ib_df$id), function(s)
        new("segPointSet",
            coords = ib_df[ib_df$id == s, c("x", "y")],
            slice = "maximum projection",
            structure = as.character(ib_df[ib_df$id == s, "structure"]),
            id =  as.character(s),
            subid =  as.character(s)))
    }

    dirs <- list("RAS" = data.frame(row.names = c("sagittal", "coronal", "axial"),
                                    "dir1" = c("LR", "PA", "IS"),
                                    "dir2" = c("RL", "AP", "SI")),
                 "PIR" = data.frame(row.names  = c("sagittal", "coronal", "axial"),
                                    "dir1" = c("LR", "AP", "SI"),
                                    "dir2" = c("RL", "PA", "IS")),
                 "LPS" = data.frame(row.names  = c("sagittal", "coronal", "axial"),
                                    "dir1" = c("RL", "AP", "IS"),
                                    "dir2" = c("LR", "PA", "SI")))


    if(name %in% names(segmentation@projections)) warning(paste0("The name ", name, " for the projection is already present, so the slot will be overwritten."))
    segmentation@projections[[name]] <- list()

    for(i in planes_chosen){
      segmentation@projections[[name]][[i]] <- list()
      segmentation@projections[[name]][[i]][[dirs[[segmentation@metadata$dirs]][i, 1]]] <- fill_list_AB[[i]]
      segmentation@projections[[name]][[i]][[dirs[[segmentation@metadata$dirs]][i, 2]]] <- fill_list_BA[[i]]
    }
    return(segmentation)
}

#' Make polygon sets
#'
#' Creates a list of polygons for every structure in every slice
#'
#' @param structure_list a list of point sets divided by structure
#'
#' @return A nested list containing ordered coordinates of points for every polygon for every structure in every slice, each point with an identifier for its grouping and hole status.
#'
#' @export

makePolygonSets <- function(structure_list,
                            parallel = FALSE,
                            verbose = FALSE) {

  # Make a polygon for every structure in every slice
  if (parallel) applyfun <- future.apply::future_lapply else applyfun <- lapply

  polygon_list <- lapply(names(structure_list), function(n) {
    applyfun(names(structure_list[[n]]), function(y) {
      if(verbose) {
        cat("\r", "Drawing from slice " , n, ", #", which(names(structure_list) == n), " of ", length(structure_list), ", structure ", y, "...", sep = "")
      }
      ib_df <- makePolygon(structure_list[[n]][[y]])
      ib_df$structure <- unique(structure_list[[n]][[y]]$structure)

      # If there is more than one cluster there could be a hole in the polygon set
      if (length(unique(ib_df$cluster)) > 1) {
        hole_list <- list()

        # For every cluster, we see whether it is contained in any other clusters
        for (i in unique(ib_df$cluster)) {
          hole_list[[as.character(i)]] <- lapply(
            unique(ib_df$cluster)[unique(ib_df$cluster) != i],
            function(x) {
              any(splancs::inout(
                pts = ib_df[ib_df$cluster == x, 1:2],
                poly = ib_df[ib_df$cluster == i, 1:2]
              ))
            }
          )
        }
        # Name the list according to the cluster that is being tested
        for (i in seq_len(length(hole_list))) {
          names(hole_list[[i]]) <- unique(ib_df$cluster)[unique(ib_df$cluster) != names(hole_list)[i]]
        }

        # Convert to a matrix for easier identification of clusters with holes
        hmat <- do.call(rbind, hole_list)
        hole_per_cluster <- apply(hmat, 1, function(x) colnames(hmat)[x == TRUE])
        hole_per_cluster <- hole_per_cluster[which(lengths(hole_per_cluster) > 0)]
        hole_per_cluster <- data.frame(
          "cluster_outer" = rep(names(hole_per_cluster), lengths(hole_per_cluster)),
          "cluster_inner" = unlist(hole_per_cluster)
        )

        # Cluster sizes are calculated to avoid multiple holes (we always take the bigger cluster to be the outer one)
        hole_per_cluster$outer_size <- table(ib_df$cluster)[hole_per_cluster$cluster_outer]
        hole_per_cluster$inner_size <- table(ib_df$cluster)[hole_per_cluster$cluster_inner]
        hole_per_cluster <- hole_per_cluster[which(hole_per_cluster$inner_size < hole_per_cluster$outer_size), ]

        # Handle duplicates (polygons that appear to be holes for more than one polygon)
        if(any(duplicated(hole_per_cluster$cluster_inner))) {
          dups <- unique(hole_per_cluster$cluster_inner[which(duplicated(hole_per_cluster$cluster_inner))])
          discard_list <- list()
          for(i in seq_len(length(dups))) {
            check_size <- hole_per_cluster[hole_per_cluster$cluster_inner == dups[[i]],]
            discard_list[[i]] <- rownames(check_size)[-which.max(check_size$outer_size)]
          }
          hole_per_cluster <- hole_per_cluster[setdiff(rownames(hole_per_cluster), unlist(discard_list)),]
        }
        # Initialize every point as being outer
        ib_df$hole <- "outer"

        # Clusters are now going to change (if cluster A is within cluster B, its points become B as well)
        ib_df$cluster_hole <- ib_df$cluster

        for (i in unique(hole_per_cluster$cluster_inner)) {
          ib_df$hole[ib_df$cluster == i] <- "inner"
          ib_df$cluster_hole[ib_df$cluster == i] <- hole_per_cluster[hole_per_cluster$cluster_inner == i, 1]
        }

        ib_df$id_hole <- paste0(ib_df$structure, "_", ib_df$cluster_hole, "_", ib_df$cluster)
      } else {

        # No need to do this if there is only one polygon in this structure
        ib_df$hole <- "outer"
        ib_df$cluster_hole <- ib_df$cluster
        ib_df$id_hole <- paste0(ib_df$structure, "_", ib_df$cluster_hole, "_", ib_df$cluster)
      }

      # Columns necessary for aesthetics drawing
      # group aesthetic
      ib_df$id <- paste0(ib_df$structure, "_", ib_df$cluster_hole)
      # subgroup aesthetic
      ib_df$subid <- paste0(ib_df$id_hole, "_", ib_df$hole)

      # We return a list splitting on the cluster
      # this way reassigning values to the smoothed polygon is easier

      return(lapply(unique(ib_df$subid), function(x) {
        ib_subid <- ib_df[ib_df$subid == x, ]

        return(new("segPointSet",
                   coords = ib_subid[, c("x", "y")],
                   slice = n,
                   structure = as.character(unique(ib_subid$structure)),
                   id = as.character(unique(ib_subid$id)),
                   subid = as.character(unique(ib_subid$subid))
        ))
      }))
    }) # end of second lapply
  }) # end of first lapply
  names(polygon_list) <- names(structure_list)
  return(polygon_list)
}


#' Create a palette for an ontology/tree
#'
#' Generates a maximally separated color palette for tree nodes taking into account the level at which colors should be initially resolved
#'
#' @param o a data frame containing the ontology/tree
#' @param graph_depth a numeric defining the graph depth at which the initial colouring should be resolved. Can only be used if \code{start_nodes} is NULL
#' @param start_nodes a character vector with the acronyms of the nodes that are being selected for initial colouring
#' @param colorspace the color space to be used as input for \code{qualpal} - see \code{?qualpal} for more information. Defaults to \code{"pretty"}.
#' @param default_col a character indicating the color to be used for nodes that will not be colored
#'
#' @return a vector of colors ordered by the same order as the nodes
#'
#' @export

makeOntologyPalettes <- function(o,
                                 graph_depth = NULL,
                                 start_nodes = NULL,
                                 colorspace = "pretty",
                                 default_col = "white") {

  # Start the colors with the default
  rownames(o) <- o$id
  o$col <- default_col
  o$visited <- FALSE

  if (is.null(start_nodes) & is.null(graph_depth)) stop("Must choose one between start_nodes or graph_depth")
  if (!is.null(start_nodes) & !is.null(graph_depth)) stop("Must choose one between start_nodes or graph_depth")

  # Maximally spaced color palette for the graph depth we will group on
  if (is.null(start_nodes) & !is.null(graph_depth)) {
    st_pal <- qualpalr::qualpal(n = sum(o$graph_depth == graph_depth & o$has_children == TRUE), colorspace = colorspace)
    st_pal_hex <- st_pal$hex
    st_pal_HSL <- st_pal$HSL
    names(st_pal_hex) <- o$id[o$graph_depth == graph_depth & o$has_children == TRUE]
  } else if (!is.null(start_nodes) & is.null(graph_depth)) {
    st_pal <- qualpalr::qualpal(n = sum(o$acronym %in% start_nodes), colorspace = colorspace)
    st_pal_hex <- st_pal$hex
    st_pal_HSL <- st_pal$HSL
    names(st_pal_hex) <- o$id[o$acronym %in% start_nodes]
  }

  # Assign colors to selected depth/nodes
  o[names(st_pal_hex), "col"] <- st_pal_hex
  o[names(st_pal_hex), "visited"] <- TRUE

  # Assign color palettes to structures that are drawn (no children)
  for (i in names(st_pal_hex)) {
    colvec <- o$col[grepl(i, o$parent_path) & o$has_children == FALSE & o$visited == FALSE]

    if (length(colvec) < 2) {
      colvec <- st_pal_hex[i]
    } else {
      H_factor <- 15 / (length(colvec)^2)

      colvec <- hex(colorspace::HLS(
        seq(st_pal_HSL[st_pal_hex[i], 1], st_pal_HSL[st_pal_hex[i], 1] + (H_factor * length(colvec)), length.out = length(colvec)),
        scales::rescale(seq(st_pal_HSL[st_pal_hex[i], 3], st_pal_HSL[st_pal_hex[i], 3] - (0.045 * length(colvec)), length.out = length(colvec)), to = c(st_pal_HSL[st_pal_hex[i], 3], 0.2)),
        scales::rescale(seq(st_pal_HSL[st_pal_hex[i], 2], st_pal_HSL[st_pal_hex[i], 2] + (0.015 * length(colvec)), length.out = length(colvec)), to = c(st_pal_HSL[st_pal_hex[i], 2], 0.8))
      ))

      o$col[grepl(i, o$parent_path) & o$has_children == FALSE & o$visited == FALSE] <- colvec
      o$visited[grepl(i, o$parent_path) & o$has_children == FALSE & o$visited == FALSE] <- TRUE

      # Color the rest of the nodes with children but not root using the root node color
      o$col[grepl(i, o$parent_path) & o$has_children == TRUE & o$visited == FALSE] <- st_pal_hex[i]
      o$visited[grepl(i, o$parent_path) & o$has_children == TRUE & o$visited == FALSE] <- TRUE
    }
  }

  # Warn the user about nodes they did not pick/will not be coloured
  if (any(o$col == default_col)) warning(paste0("The following nodes will not be coloured: ", paste(o$acronym[o$col == default_col], collapse = ", ")))

  return(o$col)
}
