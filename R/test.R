#' Make slices test
#'
#' Divides a molten 3D array in slices along a selected plane
#'
#' @param M a data frame molten from a 3D array in long format, with x, y and z columns as \code{"Var1"}, \code{"Var2"}, \code{"Var3"} and the value column containing the voxel annotation
#' @param plane character, one of \code{"sagittal"}, \code{"coronal"}, or \code{"axial"}. The plane along which to slice the array. Default is \code{"sagittal"}
#' @return A list of lists of structures, one list per slice
#'
#' @export

slice_make_test <- function(M,
                       plane = "sagittal") {
  
  if(plane == "sagittal"){
    uc <- 1
    xc <- 2
    yc <- 3
  } else if (plane == "coronal"){
    uc <- 2
    xc <- 1
    yc <- 3
  } else if (plane == "axial") {
    uc <- 3
    xc <- 1
    yc <- 2
  } else stop("Must select a `plane` out of \"sagittal\", \"coronal\", \"axial\"")
  
  
  axis_plane <- data.table::data.table(split(M, by = colnames(M)[uc]))[[1]]
  
  columns <- c(xc, yc, uc, 4)
  
  for (i in seq_len(length(axis_plane))) {
    axis_plane[[i]] <- data.table::as.data.table(axis_plane[[i]][, columns, with = FALSE])
    colnames(axis_plane[[i]]) <- c("x", "y", "slice", "structure")
  }
  
  names(axis_plane) <- unique(M[, uc, with = FALSE])[[1]]
  axis_plane <- axis_plane[order(as.numeric(names(axis_plane)))]
  
  
  axis_strlist <- lapply(axis_plane, function(x) split(x, by = "structure"))
  
  # Naming structures using their codes
  for (i in seq_len(length(axis_strlist))) {
    names(axis_strlist[[i]]) <- unlist(lapply(
      axis_strlist[[i]],
      function(x) unique(x$structure)
    ))
  }
  names(axis_strlist) <- names(axis_plane)
  
  return(axis_strlist)
  
}
