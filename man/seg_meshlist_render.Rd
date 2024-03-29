% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coldcuts3D.R
\name{seg_meshlist_render}
\alias{seg_meshlist_render}
\title{Render meshes from a segmentation}
\usage{
seg_meshlist_render(
  segmentation,
  subset_str = NULL,
  iterations = 4,
  group_first = TRUE,
  style = "matte",
  mono_color = "gray"
)
}
\arguments{
\item{segmentation}{a `segmentation` class object}

\item{subset_str}{a character vector indicating the structure(s) to be rendered, by acronym. Default is \code{NULL}, meaning the whole segmentation will be rendered as a single mesh.}

\item{iterations}{numeric, iterations for HC smoothing of the meshes}

\item{group_first}{logical, in the case of ambiguous acronyms, should the rendering privilege plotting the group or the single structure? Default is \code{TRUE}, meaning the group will be given priority.}

\item{style}{character, one of "matte", "shiny" or "monochrome" as a rendering style}

\item{mono_color}{character, the color of meshes when character = "monochrome"}
}
\value{
a plot of all meshes selected from the segmentation `meshes` slot, smoothed and rendered in a \code{rgl} window.
}
\description{
Renders a list of mesh3d objects using ontology-defined colors
}
