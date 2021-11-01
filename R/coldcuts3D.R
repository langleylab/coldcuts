#' Build segmentation meshes in 3D
#'
#' Uses marching cubes and mesh reduction to render the segmentation in 3 dimensions
#'
#' @param segmentation a `segmentation` class object
#' @param subset_str a character vector indicating the structure(s) to be rendered. Default is \code{NULL}, meaning the whole segmentation will be rendered as a single mesh.
#' @param pct_reduce numeric, target percentage of vertices for downsampling and remeshing. Closer to 0 gives a simpler mesh. Default is 0.1.
#' @param verbose logical, show progress of the meshing? Default is \code{FALSE}.
#' @return a smoothed 3d mesh (as `mesh3d` object) from one or more structures
#'
#' @export

seg_build_meshes <- function(segmentation, 
                          subset_str = NULL,
                          pct_reduce = 0.1,
                          verbose = FALSE) {
  
  if(!"rgl" %in% rownames(installed.packages()) | !"Rvcg" %in% rownames(installed.packages())) stop("In order to use 3D rendering you must first install `rgl` and `Rvcg`.")
  if(is.null(subset_str)) stop("You must specify a structure acronym")
  if(!subset_str %in% ontology(segmentation)[as.character(segmentation@metadata$structures), "acronym"]) stop(paste0("Structure ", subset_str, " was not found in this segmentation."))
  
  subset_str_id <- ontology(segmentation)$id[ontology(segmentation)$acronym == subset_str]
  segmentation <- seg_sub_str(segmentation, structures = subset_str_id)
  
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

#' Renders mesh list
#'
#' Renders a list of mesh3d objects using ontology-defined colors
#'
#' @param segmentation a `segmentation` class object
#' @param subset_str a character vector indicating the structure(s) to be rendered. Default is \code{NULL}, meaning the whole segmentation will be rendered as a single mesh.
#' @param iterations numeric, iterations for HC smoothing of the meshes
#' @param style character, one of "matte" or "shiny" as a rendering style
#' 
#' @return a plot of all meshes selected from the meshlist in an \code{rgl} window.
#'
#' @export

seg_render_meshes <- function(segmentation, 
                            subset_str = NULL,
                            iterations = 4,
                            style = "matte") {
  
  if(!"rgl" %in% rownames(installed.packages()) | !"Rvcg" %in% rownames(installed.packages())) stop("In order to use 3D rendering you must first install `rgl` and `Rvcg`.")
  if(!is.null(subset_str) & any(!subset_str %in% names(segmentation@meshes))) stop("Some structures were not found in the meshes slot.")
  
  if(is.null(subset_str)) {
      indices = names(segmentation@meshes)
    } else {
      indices = subset_str 
    }
  
  style <- match.arg(style, choices = c("shiny","matte"))
  
  for(i in indices) {
    
    mesh_to_plot <- segmentation@meshes[[i]]
    
    if(style == "matte") {
      mesh_to_plot$material <- list()
      mesh_to_plot$material$specular <- "black"
    } 
           rgl::shade3d(rgl::rotate3d(
                 Rvcg::vcgSmooth(mesh_to_plot, 
                                 "HC", 
                                 iteration = iterations), 
                 angle = pi, x = 0, y = 1, z = 0), 
                 material = segmentation@ontology[segmentation@ontology$acronym == i, "col"], 
                 meshColor = "faces")
    } 
}

