#' Build a structure mesh
#'
#' Uses marching cubes and mesh reduction to render the segmentation in 3 dimensions in a single mesh.
#'
#' @param segmentation a `segmentation` class object
#' @param subset_str a character vector indicating the structure(s) to be rendered, by acronym.
#' @param pct_reduce numeric, target percentage of vertices for downsampling and remeshing. Closer to 0 gives a simpler mesh. Default is 0.1.
#' @param verbose logical, show progress of the meshing? Default is \code{FALSE}.
#' @return a reduced 3d mesh (as `mesh3d` object) from one or more structures
#'
#' @importFrom utils installed.packages
#' 
#' @export

seg_mesh_build <- function(segmentation,
                          subset_str,
                          pct_reduce = 0.1,
                          verbose = FALSE) {

  if(!"rgl" %in% rownames(installed.packages()) | !"Rvcg" %in% rownames(installed.packages())) stop("In order to use 3D rendering you must first install `rgl` and `Rvcg`.")
  if(class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")
  if(is.null(subset_str)) stop("You must specify a structure acronym")
  if(!subset_str %in% ontology(segmentation)[as.character(segmentation@metadata$structures), "acronym"]) stop(paste0("Structure ", subset_str, " was not found in this segmentation."))
  if(pct_reduce <= 0 | pct_reduce > 1) stop("pct_reduce must be between 0 (excluded) and 1 (included).")

  segmentation <- seg_sub_str(segmentation, structures = subset_str)

  if(verbose) cat("Calculating holes...")

  poly_list <- lapply(segmentation@slices[["sagittal"]], function(x)
    lapply(x, function(y) lapply(y, poly_build)))


  poly_holed <- lapply(poly_list, function(x) {

      keep_borders <- do.call(rbind, lapply(x, function(y) do.call(rbind, lapply(y, function(z) z[grepl("inner", z$subid),]))))

      full_df <- do.call(rbind, lapply(x, function(y) do.call(rbind, lapply(y, poly_fill))))
      reordered_df <- rbind(full_df[grepl("outer", full_df$subid),], full_df[grepl("inner", full_df$subid),])
      reordered_df <- reordered_df[!duplicated(reordered_df[,1:2],fromLast = TRUE),]
      reordered_df <- reordered_df[grepl("outer", reordered_df$subid),]
      reordered_df <- rbind(reordered_df, keep_borders)

      return(reordered_df)
    })

     poly_table <- do.call(rbind, poly_holed)

    poly_table <- poly_table[,1:4]
    poly_table$x <- as.numeric(poly_table$x)
    poly_table$y <- as.numeric(poly_table$y)
    poly_table$slice <- as.numeric(poly_table$slice)
    poly_table$structure <- as.numeric(poly_table$structure)

    if(verbose) cat("done.\n")

    if(verbose) cat("Adding bounding box...")

    cube = array(0, dim = c(diff(range(poly_table$x)) + 10, diff(range(poly_table$y)) + 10, diff(range(poly_table$slice)) + 10))
    mcube <- reshape2::melt(cube)
    colnames(mcube) <- c("x", "y", "slice")

    mcube$x <- mcube$x + min(range(poly_table$x)) - 5
    mcube$y <- mcube$y + min(range(poly_table$y)) - 5
    mcube$slice <- mcube$slice + min(range(poly_table$slice)) - 5

    colnames(mcube) <- c("x", "y", "slice", "structure")
    mcube <- mcube[,colnames(poly_table)]
    poly_table_bound <- rbind(poly_table, mcube)
    poly_table_bound <- poly_table_bound[!duplicated(poly_table_bound[,1:3]),]

    if(verbose) cat("done.\n")

    if(verbose) cat("Rebuilding array...")

    reconstructed_array <- reshape2::acast(poly_table_bound, formula = x ~ y ~ slice, value.var = "structure", fill = 0)

    reconstructed_array[!reconstructed_array %in% unique(poly_table_bound$structure)] <- 0

    if(verbose) cat("done.\n")

    if(verbose) cat("Creating mesh...")

    ijk2ras.90x <- structure(c(1, 0, 0, 0,
                               0, 0, -1, 0,
                               0, 1, 0, 0,
                               0, 0 ,0, 1),
                             .Dim = c(4,4))

    seg_mesh <- Rvcg::vcgIsosurface(reconstructed_array, threshold = 1, IJK2RAS = ijk2ras.90x)

    if(verbose) cat("done.\n")

    if(verbose) cat(paste0("Reducing mesh at ", pct_reduce * 100, "%..."))

    reduced.mesh <- Rvcg::vcgQEdecim(mesh = seg_mesh, percent = pct_reduce, silent = !verbose)

    reduced.mesh$vb[1,] <- reduced.mesh$vb[1,] + min(poly_table[poly_table$structure != 0, "x"])
    reduced.mesh$vb[2,] <- reduced.mesh$vb[2,] + min(poly_table[poly_table$structure != 0, "slice"])
    reduced.mesh$vb[3,] <- reduced.mesh$vb[3,] - min(poly_table[poly_table$structure != 0, "y"])

    if(verbose) cat("done.\n")

    return(reduced.mesh)

}

#' Render meshes from a segmentation
#'
#' Renders a list of mesh3d objects using ontology-defined colors
#'
#' @param segmentation a `segmentation` class object
#' @param subset_str a character vector indicating the structure(s) to be rendered, by acronym. Default is \code{NULL}, meaning the whole segmentation will be rendered as a single mesh.
#' @param iterations numeric, iterations for HC smoothing of the meshes
#' @param group_first logical, in the case of ambiguous acronyms, should the rendering privilege plotting the group or the single structure? Default is \code{TRUE}, meaning the group will be given priority.
#' @param style character, one of "matte", "shiny" or "monochrome" as a rendering style
#' @param mono_color character, the color of meshes when character = "monochrome"
#'
#' @return a plot of all meshes selected from the segmentation `meshes` slot, smoothed and rendered in a \code{rgl} window.
#'
#' @importFrom utils installed.packages
#' 
#' @export

seg_meshlist_render <- function(segmentation,
                            subset_str = NULL,
                            iterations = 4,
                            group_first = TRUE,
                            style = "matte",
                            mono_color = "gray") {

  if(!"rgl" %in% rownames(installed.packages()) | !"Rvcg" %in% rownames(installed.packages())) stop("In order to use 3D rendering you must first install `rgl` and `Rvcg`.")

  if(class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")

  if(class(iterations) != "numeric") stop("Iterations must be an integer")
  if(iterations %% 1 != 0) stop("Iterations must be an integer or coercible to an integer")

  #if(!is.null(subset_str) & any(!subset_str %in% names(segmentation@meshes))) stop("Some structures were not found in the meshes slot.")
  iterations <- as.integer(iterations)

  if(is.null(subset_str)) {

      indices = names(segmentation@meshes)

    } else {

      indices = intersect(seg_select_str(segmentation, structures = subset_str, group_first = group_first), names(segmentation@meshes))
    }

  style <- match.arg(style, choices = c("shiny","matte", "monochrome"))

  for(i in indices) {

    mesh_to_plot <- segmentation@meshes[[i]]

    if(style == "matte") {
      mesh_to_plot$material$specular <- "black"
    }
    if(style == "monochrome") {
      mesh_to_plot$material$specular <- "black"
      mesh_material = mono_color
    } else {
      mesh_material = segmentation@ontology[segmentation@ontology$acronym == i, "col"]
    }
          rgl::shade3d(rgl::rotate3d(
                 Rvcg::vcgSmooth(mesh_to_plot,
                                 "HC",
                                 iteration = iterations),
                 angle = pi, x = 0, y = 1, z = 0),
                 material = mesh_material,
                 meshColor = "vertices",
                 depth_mask = TRUE)
    }
}


#' Add mesh list to a segmentation
#'
#' Builds a list of structure meshes and stores it in the meshes slot.
#'
#' @param segmentation a `segmentation` class object
#' @param structures a character vector of structure acronyms, indicating the structure(s) to be rendered. Default is \code{NULL}, meaning all structures will be rendered.
#' @param pct_reduce numeric, target percentage of vertices for downsampling and remeshing. Closer to 0 gives a simpler mesh. Default is 0.1.
#' @param verbose logical, show progress of the meshing? Default is \code{FALSE}.
#' @return a `segmentation` class object with a list of triangular meshes (as `mesh3d` objects) in the `meshes` slot.
#'
#' @importFrom utils installed.packages
#' 
#' @export

seg_meshlist_add <- function(segmentation,
                             structures = NULL,
                             pct_reduce = 0.1,
                             verbose = FALSE){

  if(!"rgl" %in% rownames(installed.packages()) | !"Rvcg" %in% rownames(installed.packages())) stop("In order to use 3D rendering you must first install `rgl` and `Rvcg`.")
  if(class(segmentation) != "segmentation") stop("You must provide a segmentation class object.")
  if(class(pct_reduce) != "numeric") stop("pct_reduce must be a numeric.")
  if(pct_reduce <= 0 | pct_reduce > 1) stop("pct_reduce must be between 0 (excluded) and 1 (included).")

  if(is.null(structures)) {

    subset_str = ontology(segmentation)[as.character(seg_metadata(segmentation)$structures), "acronym"]

    } else {

      if(any(!structures %in% ontology(segmentation)$acronym)) stop("Some structures were not found in this segmentation.")
      if(any(!structures %in% ontology(segmentation)[as.character(seg_metadata(segmentation)$structures), "acronym"])){
        ontology(segmentation)$has_children <- sapply(ontology(segmentation)$id, function(x) x %in% ontology(segmentation)$parent_structure_id)
        str_with_children <- ontology(segmentation)$id[which(ontology(segmentation)$acronym %in% structures & ontology(segmentation)$has_children)]
        str_without_children <- ontology(segmentation)$acronym[which(ontology(segmentation)$acronym %in% structures & !ontology(segmentation)$has_children)]
        subset_str <- c(ontology(segmentation)[unlist(sapply(str_with_children, function(x) grepl(x, ontology(segmentation)$structure_id_path))) &
                                             !ontology(segmentation)$has_children, "acronym"], str_without_children)
        }
    }

  meshlist <- list()
  for(i in subset_str) {
    if(verbose) cat("Rendering structure ", i, ", ", which(subset_str == i), " of ", length(subset_str), "...\n", sep = "")
    meshlist[[i]] <- seg_mesh_build(segmentation = segmentation, subset_str = i, pct_reduce = pct_reduce, verbose = verbose)
  }

  segmentation@meshes <- meshlist

  if(verbose) cat("All meshes built!\n")

  return(segmentation)
}


