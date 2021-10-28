<!-- badges: start -->
  [![](https://img.shields.io/badge/devel%20version-0.0.0.900-orange.svg)](https://github.com/langleylab/coldcuts)
  [![R-CMD-check](https://github.com/langleylab/coldcuts/workflows/R-CMD-check/badge.svg)](https://github.com/langleylab/coldcuts/actions)
<!-- badges: end -->

<img src="https://user-images.githubusercontent.com/21171362/131806935-ed890016-a845-4274-8cb3-cd78c16aeb00.png" align="right" alt="" width="200" />

# coldcuts 

**`coldcuts`** is an R package that allows you to **draw and plot automatically** segmentations from 3D voxel arrays. 

The name is inspired by one of Italy's best products.

🎓 You can find the documentation and a tutorial to get started at the package's page: langleylab.github.io/coldcuts

🗂 You can find additional segmentation files, ontologies and other information at https://langleylab.github.io/coldcuts/articles/segmentations.html

## Motivation

When dealing with neuroimaging data, or any other type of numerical data derived from brain tissues, it is important to situate it in its anatomical and structural context. Various authors provide parcellations or segmentations of the brain, according to their best interpretation of which functional and anatomical boundaries make sense for our understanding of the brain. There are several stand-alone tools that allow to visualize and manipulate segmentations. However, neuroimaging data, together with other functional data such as transcriptomics, is often manipulated in a statistical programming
language such as R which does not have trivial implementations for the visualization of segmentations.

To bridge this gap, some R packages have been recently published:

 - [ggseg](https://github.com/ggseg/ggseg) by Athanasia Mo Mowinckel and Didac Vidal-Piñeiro
 - [cerebroViz](https://github.com/ethanbahl/cerebroViz) by Ethan Bahl, Tanner Koomar, and Jacob J Michaelson
 - [fsbrain](https://github.com/dfsp-spirit/fsbrain) by Tim Schäfer and Christine Ecker

**`ggseg`** and **`cerebroviz`** offer 2D (and 3D in the case of **`ggseg3d`**) visualizations of human brain segmentations, with the possibility of integration with external datasets. These segmentations are manually curated, which means that new datasets must be manually inserted, and they are limited to the human brain in scope. **`ggseg`** in particular has made available several segmentations of human cortical surface atlases. 
**`fsbrain`** focuses on 3D visualization of human MRI data with external data integration and visualization in both native space and transformed spaces. It does not depend on manually curated datastes (beyond segmentations).

While these tools provide a wealth of beautiful visualization interfaces, we felt the need to implement a tool to systematically create 2D (and potentially 3D) objects that are easily shared and manipulated in R, with the addition of labels, external datasets and simple operations such as subsetting and projecting, with minimal need for manual curation and without limiting ourselves to a particular species. 

Thus, **`coldcuts`** is our attempt at bridging the gap between imaging/high throughput brain data and R through data visualization.


## Installing the package

⬇️ You can install this package using `devtools::install_github()`:

```{r}
devtools::install_github("langleylab/coldcuts")
```

## Getting started

🏃🏽‍♀️ You can find a small example to get started [here](https://langleylab.github.io/coldcuts/articles/coldcuts.html)



