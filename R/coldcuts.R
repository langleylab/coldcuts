#' Draw segmentation
#'
#' Creates a `segmentation` object from an array or a `NIfTI` file
#' 
#' @param nifti_file a character string pointing to the NIfTI file address. Default is \code{NULL}
#' @param nrrd_file a character string pointing to the NRRD file address. Default is \code{NULL}
#' @param array a 3D array containing structure annotations (optional if \code{nifti_file} and \code{molten_array} are not specified). Default is \code{NULL}
#' @param molten_array an array in molten form, such as the result of \code{reshape2::melt()} (optional if \code{nifti_file} and \code{array} are not specified). Default is \code{NULL}
#' @param ontology an object containing the segmentation ontology (e.g. from another segmentation). Default is \code{NULL}
#' @param ontology_file a character string pointing to the ontology .csv or .txt file address. Default is \code{NULL}
#' @param outliers a numeric vector indicating outlier points to be eliminated before drawing polygons. Default is \code{NULL}
#' @param verbose logical, should messages on the process be displayed? Default is \code{TRUE}
#' @param directions_from character, indicates the original orientation of the array. Default is \code{"RAS"} (Left to Right, Posterior to Anterior, Inferior to Superior).
#' @param directions_to character, indicates the final orientation of the array. The array is rotated only if it is different from \code{directions_from}. Default is \code{"RAS"}, which results in the proper assignment of axis labels.
#' @param planes character, either \code{"all"} (default) or any subset of \code{"sagittal"}, \code{"coronal"}, and \code{"axial"}. Determines along which axes the array wil be sliced.
#' @param subset_structures character vector, structures to be subset for the segmentation. Default is \code{NULL}
#' @param draw_outline logical, should the outline of the whole array be drawn? Default is \code{TRUE}
#' @param subset_sagittal numeric vector, indicates a subset of slices along the sagittal plane that will be drawn. Default is \code{NULL}
#' @param subset_coronal numeric vector, indicates a subset of slices along the coronal plane that will be drawn. Default is \code{NULL}
#' @param subset_axial numeric vector, indicates a subset of slices along the axial plane that will be drawn. Default is \code{NULL}
#' @param reference_space character, what reference space does this segmentation use, e.g. \code{"MNI152"}, \code{"CCFv3"}, \code{"JFRC2010"}. Default is \code{NULL}
#' @param citation list of character strings, what is the citation for this segmentation? Default is \code{NULL}.
#' @param parallel logical, should the polygons be drawn using parallel processing? Uses \code{future.apply::future_lapply()}, so the proper setup must be done beforehand. Default is \code{FALSE}.
#'
#' @return a \code{segmentation} S4 object.
#'
#' @details This function takes a volume file such as a NIfTI, NRRD or 3D array
#'    and performs all the slicing and contour drawing operations, loading metadata,
#'    subsetting structures, etc.
#'
#'    The output of this function is a \code{segmentation}
#'    object that can be passed to plotting functions and can be accessed with
#'    several auxiliary methods.
#'
#'    Ideally, the user only needs to run this function
#'    once, and save the resulting object (as a RDS file, for example) for distribution
#'    and future usage. All the operations performed by this function assume the
#'    orientation of the volume is RAS.
#'
#'    Additional arguments such as
#'    \code{directions_from} and \code{directions_to} can be supplied to rotate
#'    the volume if its native orientation is different from RAS. Outlier points
#'    can be removed by specifying the index in the volume array in the \code{outliers}
#'    argument.
#'
#'    The user can choose to limit the final \code{segmentation} object
#'    to a subset of either structures, planes, or specific slices (indexed by plane).
#'
#' @examples
#' array_cubes <- dummy_cubes(10)
#' ontology <- dummy_ontology()
#' 
#' seg <- seg_draw(
#'   array = array_cubes,
#'   ontology = ontology,
#'   verbose = FALSE
#' )
#' 
#' # Keeping only the sagittal plane
#'  seg_sagittal <- seg_draw(
#'   array = array_cubes,
#'   ontology = ontology,
#'   planes = "sagittal",
#'   verbose = FALSE
#' )
#' 
#' # Keeping only a subset of structures
#'  seg_sub <- seg_draw(
#'   array = array_cubes,
#'   ontology = ontology,
#'   subset_structures = c("A", "B"),
#'   verbose = FALSE
#' )
#' 
#'
#' @importFrom data.table data.table as.data.table
#' @importFrom oro.nifti readNIfTI img_data pixdim xyzt_units
#' @importFrom utils head read.csv read.table
#' @importFrom methods new
#'
#' @export

seg_draw <- function(nifti_file = NULL, nrrd_file = NULL, array = NULL,
                     molten_array = NULL, ontology = NULL, ontology_file = NULL,
                     outliers = NULL, verbose = TRUE, directions_from = "RAS",
                     directions_to = "RAS", planes = "all",
                     subset_structures = NULL, draw_outline = TRUE,
                     subset_sagittal = NULL, subset_coronal = NULL,
                     subset_axial = NULL, reference_space = NULL,
                     citation = NULL, parallel = FALSE) {
  # Sanity checks
  if (!is.null(nrrd_file) & !"nat" %in% rownames(installed.packages())) stop("In order to read a NRRD file you must first install the package `nat`.")
  if (is.null(nifti_file) & is.null(nrrd_file) & is.null(array) & is.null(molten_array)) stop("Must provide at least a NIfTI/NRRD file, or an array/molten array.")
  if (is.null(ontology) & is.null(ontology_file)) stop("Must provide at least an ontology or an ontology file.")
  if (planes != "all" & !any(planes %in% c("sagittal", "coronal", "axial"))) stop("You must choose at least one plane among sagittal, coronal or axial.")
  if (!is.null(subset_sagittal) & is(subset_sagittal,"numeric")) stop("You must provide numeric indices to subset planes")
  if (!is.null(subset_coronal) & is(subset_coronal, "numeric")) stop("You must provide numeric indices to subset planes")
  if (!is.null(subset_axial) & is(subset_axial, "numeric")) stop("You must provide numeric indices to subset planes")

  if (verbose) cat("Adding ontology...")
  if (!is.null(ontology_file) & is.null(ontology)) {
    if (grepl("\\.csv$", ontology_file)) {
      read_fun <- read.csv
      read_sep <- ","
    } else if (grepl("\\.txt$", ontology_file)) {
      read_fun <- read.table
      read_sep <- "\t"
    }

    ontology <- read_fun(ontology_file, header = TRUE, sep = read_sep)
  } else if (is.null(ontology_file) & !is.null(ontology)) {
    ontology <- ontology
  }

  missing_fields <- setdiff(c("id", "name", "acronym", "parent_structure_id", "structure_id_path", "col"), colnames(ontology))

  if (length(missing_fields) > 0) stop(paste0("Columns ", paste(missing_fields, collapse = ", ")), " were not found in the ontology. \n Please provide a complete ontology table.")

  colnames(ontology)[grep("id|ID|Id|iD", fixed = TRUE, colnames(ontology))] <- "id"
  ontology$id <- as.character(ontology$id)
  rownames(ontology) <- ontology$id

  if (any(!subset_structures %in% ontology$acronym)) {
    stop(paste0("Subset structures ", paste(setdiff(subset_structures, ontology$acronym), collapse = ", "), " were not found in the ontology. \n Make sure you are using the correct ontology for this volume, or subsetting the right structures."))
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

  if (!is.null(nifti_file) & is.null(array) & is.null(molten_array) & is.null(nrrd_file)) {
    if (verbose) cat("Reading NIfTI file...")
    nifti <- readNIfTI(nifti_file)
    if (verbose) cat("done.\n")
    n_image <- img_data(nifti)
    if (sum(pixdim(nifti)[2:5] > 0) > 3) {
      warning("This file has 4 recorded dimensions, but only the first 3 will be used.", immediate. = TRUE)
      n_image <- array(n_image, dim = dim(n_image)[1:3])
    }
    pixdims <- pixdim(nifti)[2:4]
    units <- xyzt_units(nifti)
    filename <- nifti_file
  } else if (is.null(nifti_file) & is.null(array) & is.null(molten_array) & !is.null(nrrd_file)) {
    if (verbose) cat("Reading NRRD file...")
    n_image <- nat::read.nrrd(file = nrrd_file, ReadData = TRUE)
    filename <- nrrd_file
    pixdims <- "Unknown"
    units <- "Unknown"
    if (verbose) cat("done.\n")
  } else if (is.null(nifti_file) & !is.null(array) & is.null(molten_array) & is.null(nrrd_file)) {
    n_image <- array
    pixdims <- "Unknown"
    units <- "Unknown"
    filename <- paste0("Array named \"", deparse(substitute(molten_array)), "\".")
  } else if (is.null(nifti_file) & is.null(array) & !is.null(molten_array) & is.null(nrrd_file)) {
    pixdims <- "Unknown"
    units <- "Unknown"
    filename <- paste0("Molten array named \"", deparse(substitute(molten_array)), "\".")
  }

  if (directions_from != directions_to) {
    dirs <- directions_to
    if (verbose) cat("Changing directions...")
    n_image <- dir_change(n_image, directions_from, directions_to)
    if (verbose) cat("done.\n")
  } else {
    dirs <- directions_from
  }

  if (is.null(molten_array)) {
    if (!is.null(outliers)) n_image[t(outliers)] <- 0
    ndims <- dim(n_image)
  }

  if (planes == "all") planes_chosen <- c("sagittal", "coronal", "axial") else planes_chosen <- planes
  slices <- outlines <- list("sagittal" = NULL, "coronal" = NULL, "axial" = NULL)

  if (!is.null(array) & is.null(molten_array) & (is.null(subset_sagittal) & is.null(subset_coronal) & is.null(subset_axial))) {
    if (verbose) cat("Melting array...")
    M <- as.data.table(n_image)
    M <- M[M$value > 0, ]
    if (verbose) cat("done.\n")
  } else if (is.null(array) & !is.null(molten_array)) {
    M <- molten_array
    if (any(M$value == 0)) M <- M[M$value > 0, ]
  } else if (is.null(array) & is.null(molten_array) & (!is.null(nifti_file) | !is.null(nrrd_file))) {
    if (verbose) cat("Melting array...")
    M <- as.data.table(n_image)
    M <- M[M$value > 0, ]
    if (verbose) cat("done.\n")
  }

  if (!is.null(subset_structures)) {
    if (verbose) cat("Subsetting structures...")
    subset_structures_id = ontology$id[ontology$acronym %in% subset_structures]
    M <- M[M$value %in% subset_structures_id, ]
    if (verbose) cat("done.\n")
  }

  if (any(!unique(M$value) %in% ontology$id)) {
    missing_strs <- setdiff(unique(M$value), ontology$id)
    max <- min(c(10, length(missing_strs)))
    if (length(missing_strs) > 10) error_add_str <- paste0(" and ", length(missing_strs) - 10, " more ") else error_add_str <- ""
    stop(paste0("Structures ", paste(setdiff(unique(M$value)[1:max], ontology$id), collapse = ", "), error_add_str, "were found in the array but not in the ontology.\n Make sure you are using the correct ontology for this volume."))
  }

  if (draw_outline) {
    for (i in planes_chosen) {
      if (verbose) cat("Drawing", i, "outline...")
      ol <- outline_draw(M, plane = i)
      outlines[[i]] <- ol
      if (verbose) cat("done.\n")
    }
  }

  if (!is.null(subset_sagittal) | !is.null(subset_coronal) | !is.null(subset_axial)) {
    if (verbose) cat("Subsetting slices...")

    M_s <- reshape2::melt(n_image[subset_sagittal, , ])
    M_s$V1 <- subset_sagittal
    M_s <- M_s[, c(4, 1, 2, 3)]
    M_s <- data.table(M_s[M_s$value > 0, ])

    M_c <- reshape2::melt(n_image[, subset_coronal, ])
    M_c$V2 <- subset_coronal
    M_c <- M_c[, c(1, 4, 2, 3)]
    M_c <- data.table(M_c[M_c$value > 0, ])

    M_a <- reshape2::melt(n_image[, , subset_axial])
    M_a$V3 <- subset_axial
    M_a <- M_a[, c(1, 2, 4, 3)]
    M_a <- data.table(M_a[M_a$value > 0, ])

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
      sl <- slice_make(M, plane = i)
      if (verbose) cat("done.\n")

      if (verbose) cat("Drawing polygons...")
      sp <- poly_set_make(sl, parallel = parallel, verbose = verbose)
      slices[[i]] <- sp
      if (verbose) cat("done.\n")
    }
  } else if (!is.null(subset_sagittal) | !is.null(subset_coronal) | !is.null(subset_axial)) {
    sl_s <- slice_make(M_s, plane = "sagittal")
    sl_c <- slice_make(M_c, plane = "coronal")
    sl_a <- slice_make(M_a, plane = "axial")

    if (verbose) cat("Drawing polygons...")
    sp_s <- poly_set_make(sl_s, parallel = parallel)
    slices[["sagittal"]] <- sp_s
    sp_c <- poly_set_make(sl_c, parallel = parallel)
    slices[["coronal"]] <- sp_c
    sp_a <- poly_set_make(sl_a, parallel = parallel)
    slices[["axial"]] <- sp_a
    if (verbose) cat("done.\n")
  }

  dims_eff <- lengths(slices)

  if (verbose) cat("Compiling structure tables...")
  str_table_all <- lapply(planes_chosen, function(n) {
    structure_slices <- data.table(cbind(
      unlist(lapply(unlist(slices[[n]]), function(x) x@structure)),
      unlist(unlist(lapply(unlist(slices[[n]]), function(x) x@slice)))
    ))
    structure_slices <- structure_slices[!duplicated(structure_slices[, 1:2]), ]
    colnames(structure_slices) <- c("structure", "slice")
    return(structure_slices)
  })
  names(str_table_all) <- planes_chosen
  for (i in names(slices)) {
    for (j in names(slices[[i]])) {
      names(slices[[i]][[j]]) <- str_table_all[[i]][str_table_all[[i]]$slice == j, structure]
    }
  }
  if (verbose) cat("done.\n")

  if (verbose) cat("Compiling metadata...")

  if (!is.null(molten_array)) {
    ndims <- c(max(M$Var1), max(M$Var2), max(M$Var3))
  }

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


#' Change array directions
#'
#' Rotates the array so that it goes from one orientation to another. Currently supports PIR, RAS and LPS. Depends on \code{freesurferformats}.
#'
#' @param array a 3D array containing voxel data.
#' @param from character, a 3-word direction code that the array is in. One of \code{PIR}, \code{RAS}, or \code{LPS}.
#' @param to character, a 3-word direction code that the array will be rotated to. One of \code{PIR}, \code{RAS}, or \code{LPS}.
#'
#' @return a 3D array rotated according to original and destination directions
#'
#' @importFrom freesurferformats rotate3D
#'
#' @export

dir_change <- function(array, from = "PIR", to = "RAS") {
  if (!from %in% c("PIR", "RAS", "LPS", "LAS")) stop("Must provide one of the following values to argument from: PIR, RAS, LPS, LAS")
  if (!to %in% c("PIR", "RAS", "LPS", "LAS")) stop("Must provide one of the following values to argument to: PIR, RAS, LPS, LAS")

  if (from == "PIR" & to == "RAS") {
    array_rot <- rotate3D(rotate3D(array, axis = 1, degrees = 90), axis = 3, degrees = 90)
  }
  if (from == "RAS" & to == "PIR") {
    array_rot <- rotate3D(rotate3D(array, axis = 3, degrees = 270), axis = 1, degrees = 270)
  }
  if (from == "LPS" & to == "RAS") {
    array_rot <- rotate3D(array, axis = 3, degrees = 180)
  }
  if (from == "RAS" & to == "LPS") {
    array_rot <- rotate3D(array, axis = 3, degrees = 180)
  }
  if (from == "LPS" & to == "PIR") {
    array_rot <- rotate3D(rotate3D(array, axis = 1, degrees = 270), axis = 3, degrees = 90)
  }
  if (from == "PIR" & to == "LPS") {
    array_rot <- rotate3D(rotate3D(array, axis = 3, degrees = 90), axis = 1, degrees = 270)
  }
  if (from == "LAS" & to == "RAS") {
    array_rot <- array[dim(array)[1]:1, , ]
  }
  if (from == "LAS" & to == "LPS") {
    array_rot <- rotate3D(array[dim(array)[1]:1, , ], axis = 3, degrees = 180)
  }
  return(array_rot)
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

seg_slice_check <- function(structures, segmentation, planes) {
  structures <- as.character(structures)
  if (any(!structures %in% unique(do.call(rbind, segmentation@structure_tables)$structure))) {
    not_found <- setdiff(structures, unique(do.call(rbind, segmentation@structure_tables)$structure))
    stop("Structure(s) ", paste(not_found, collapse = ", "), " not found in this segmentation.
  Perhaps there is a typo? Check `seg_metadata(segmentation)` and/or `segmentation@structure_tables` to see the available structures.")
  }

  structure_by_plane <- lapply(planes, function(x) {
    unique(segmentation@structure_tables[[x]][structure %in% structures, "slice"])
  })
  names(structure_by_plane) <- planes
  return(structure_by_plane)
}

#' Smooth polygon sets
#'
#' Uses kernel smoothing to smooth polygon sets
#'
#' @param polygon_set a data frame containing polygons, separated by `subid`
#' @param by a character containing the name of the column grouping polygons to be smoothed. Default is \code{subid}
#' @param smoothness a numeric indicating the smoothing parameter, passed to `smooth_ksmooth`
#' @param min_points a numeric indicating the minimum number of points to smooth. If a polygon has less than this number of points, it is not smoothed.
#'
#' @return a polygon set with smoothed coordinates
#'
#' @importFrom utils installed.packages
#'
#' @export

poly_smooth <- function(polygon_set, by = "subid", smoothness = 3,
                        min_points = 5) {
  if (!"smoothr" %in% rownames(installed.packages())) stop("In order to use smoothing you must first install the package `smoothr`.")

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


#' Draw maximum outline
#'
#' Draws a polygon obtained by projecting all polygons on the plane
#'
#' @param M a molten array (of class `data.table`) containing coordinates
#' @param plane character, one of "sagittal", "coronal", "axial"
#'
#' @return a data frame containing ordered coordinates for the outline polygon
#'
#' @importFrom data.table data.table
#'
#' @export

outline_draw <- function(M, plane) {
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
  df <- data.table(M[, columns, with = FALSE])
  df <- data.table(df[!duplicated(df[, 1:2]), ])
  colnames(df) <- c("x", "y")
  df$x = as.numeric(unlist(df$x))
  df$y = as.numeric(unlist(df$y))
  outline <- poly_make(df)
  return(outline)
}

#' Get a specific slice from a segmentation
#'
#' Renders a slice from a plane of a segmentation as a data frame
#'
#' @param segmentation a \code{segmentation} class object.
#' @param plane character, name of the plane (one of `sagittal`, `coronal`, or `axial`)
#' @param slice numeric or character, the slice number. If numeric, it will be used as index of the slice. If character, it will be the slice as named in the segmentation.
#' @param fill logical, should the function return just the polygon vertices or all voxels? Default is \code{FALSE}.
#'
#' @return a `data.frame` containing the x y coordinates for all points within the slice (filled polygons)
#'
#' @export

seg_get_slice <- function(segmentation, plane, slice, fill = FALSE) {
  if (!plane %in% c("sagittal", "coronal", "axial")) stop("The plane argument must be one of \"sagittal\", \"coronal\" or \"axial\".")
  if (!plane %in% names(segmentation@slices)) stop(paste0("the ", plane, " plane was not found in this segmentation."))
  if (is(slice, "numeric") & length(segmentation@slices[[plane]]) < slice) stop(paste0("Slice ", slice, " is higher than total amount of slices for this plane"))
  if (is(slice, "character") & !slice %in% segmentation@structure_tables[[plane]]$slice) stop(paste0("Slice ", slice, " not found in this plane"))

  if (fill) {
    df <- do.call(rbind, lapply(segmentation@slices[[plane]][[slice]], function(x) do.call(rbind, lapply(x, function(y) poly_fill(poly_build(y))))))
  } else {
    df <- do.call(rbind, lapply(segmentation@slices[[plane]][[slice]], function(x) do.call(rbind, lapply(x, function(y) poly_build(y)))))
  }
  return(df[!duplicated(df[, 1:3]), ])
}


#' Get a set of indices according to the ontology
#'
#' Retrieves all sets of indices belonging to a particular acronym
#'
#' @param segmentation a \code{segmentation} class object.
#' @param structures character, vector containing one or more structure acronyms
#' @param group_first logical, should priority be given to the group of structures or to the structure with the same acronym? Default is TRUE (prioritizing the group)
#'
#' @return a vector containing the indices of the structures
#'
#' @export

seg_select_str <- function(segmentation, structures, group_first = TRUE) {
  if (length(setdiff(structures, ontology(segmentation)$acronym)) > 0) stop("Some structures were not found in the ontology. Check spelling and/or case.")

  subset_str <- data.frame("str" = structures)
  indices <- list()

  ontology(segmentation)$has_children <- sapply(ontology(segmentation)$id, function(x) x %in% ontology(segmentation)$parent_structure_id)
  ambiguous_id <- as.character(seg_metadata(segmentation)$structures)[ontology(segmentation)[as.character(seg_metadata(segmentation)$structures), "has_children"]]
  ambiguous <- ontology(segmentation)[as.character(ambiguous_id), "acronym"]

  subset_str$id <- sapply(structures, function(x) ontology(segmentation)$id[ontology(segmentation)$acronym == x])
  subset_str$ambiguous <- structures %in% ambiguous
  subset_str$children <- ontology(segmentation)[subset_str$id, "has_children"]

  for (i in subset_str$str[!subset_str$ambiguous & !subset_str$children]) indices[[i]] <- i

  if (any(subset_str$ambiguous & subset_str$children) & group_first) {
    select_amb <- subset_str$str[subset_str$ambiguous]

    for (i in select_amb) {
      str_with_children <- ontology(segmentation)$id[which(ontology(segmentation)$acronym %in% i & ontology(segmentation)$has_children)]
      str_without_children <- ontology(segmentation)$acronym[which(ontology(segmentation)$acronym %in% i & !ontology(segmentation)$has_children)]
      preselection <- ontology(segmentation)[which(unlist(sapply(str_with_children, function(x) grepl(x, ontology(segmentation)$structure_id_path)))), ]
      preselection <- preselection[preselection$id %in% seg_metadata(segmentation)$structures, ]

      indices[[i]] <- union(preselection$acronym, i[i %in% ambiguous])
    }
  } else if (any(subset_str$ambiguous & subset_str$children) & !group_first) {
    select_amb <- subset_str$str[subset_str$ambiguous]

    for (i in select_amb) indices[[i]] <- i
  }

  if (any(!subset_str$ambiguous & subset_str$children)) {
    select_noamb_withkids <- subset_str$str[!subset_str$ambiguous & subset_str$children]

    for (i in select_noamb_withkids) {
      str_with_children <- ontology(segmentation)$id[which(ontology(segmentation)$acronym %in% i & ontology(segmentation)$has_children)]
      str_without_children <- ontology(segmentation)$acronym[which(ontology(segmentation)$acronym %in% i & !ontology(segmentation)$has_children)]
      preselection <- ontology(segmentation)[which(unlist(sapply(str_with_children, function(x) grepl(x, ontology(segmentation)$structure_id_path)))), ]
      preselection <- preselection[preselection$id %in% seg_metadata(segmentation)$structures, ]

      indices[[i]] <- preselection$acronym
    }
  }

  indices <- unlist(indices)

  return(indices)
}


#' Build dummy cubes
#'
#' Builds a 3D array with 8 dummy cubes, for demo/testing purposes only.
#'
#' @param n the number of voxels in the side of each cube 
#'
#' @return a 3D array containing 8 different cubes valued 1 to 8
#'
#' @export

dummy_cubes = function(n) {
  coords_1 = as.data.frame(expand.grid(seq_len(n), seq_len(n), seq_len(n)))
  coords_2 = as.data.frame(expand.grid(seq_len(n)+n, seq_len(n), seq_len(n)))
  coords_3 = as.data.frame(expand.grid(seq_len(n), seq_len(n)+n, seq_len(n)))
  coords_4 = as.data.frame(expand.grid(seq_len(n), seq_len(n), seq_len(n)+n))
  coords_5 = as.data.frame(expand.grid(seq_len(n)+n, seq_len(n)+n, seq_len(n)))
  coords_6 = as.data.frame(expand.grid(seq_len(n), seq_len(n)+n, seq_len(n)+n))
  coords_7 = as.data.frame(expand.grid(seq_len(n)+n, seq_len(n), seq_len(n)+n))
  coords_8 = as.data.frame(expand.grid(seq_len(n)+n, seq_len(n)+n, seq_len(n)+n))
  coords = rbind(coords_1, coords_2, coords_3, coords_4, coords_5, coords_6, coords_7, coords_8)
  colnames(coords) = c("x", "y", "z")
  coords$value = rep(seq_len(8), each = n^3)
  arr = reshape2::acast(coords, formula = x ~ y ~ z, value.var = "value")
  return(arr)
}


#' Build dummy ontology
#'
#' Builds an ontology table to go with the output of \code{dummy_cubes()}. For demo/testing purposes only.
#'
#' @return an ontology table for 8 cubes named A to H.
#'
#' @importFrom grDevices palette
#' 
#' @export


dummy_ontology = function(){
  ids = seq_len(9)
  onto = data.frame(id = ids, 
                    name = LETTERS[ids],
                    acronym = LETTERS[ids],
                    parent_structure_id = c(rep(9, 8), NA),
                    structure_id_path = c(rep("/9/", 8), NA),
                    col = c(palette()[ids[1:8]], "white"),
                    graph_depth = c(rep(2, 8), 1))
  
  return(onto)             
}

