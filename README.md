<!-- badges: start -->
  [![](https://img.shields.io/badge/devel%20version-0.0.0.900-orange.svg)](https://github.com/langleylab/coldcuts)
  [![R-CMD-check](https://github.com/langleylab/coldcuts/workflows/R-CMD-check/badge.svg)](https://github.com/langleylab/coldcuts/actions)
<!-- badges: end -->

<img src="https://user-images.githubusercontent.com/21171362/131806935-ed890016-a845-4274-8cb3-cd78c16aeb00.png" align="right" alt="" width="200" />

# coldcuts 

**`coldcuts`** is an R package that allows you to **draw and plot automatically** segmentations from 3D voxel arrays. 

The name is inspired by one of Italy's best products.

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

You can install this package using `devtools::install_github()`:

```{r}
devtools::install_github("langleylab/coldcuts")
```

## Finding additional segmentations (WIP)

Many institutions have published their own segmentations as NIfTI, NRRD or RAW files.
In this table we curate a few segmentations including their source(s).

| Name | Organism | Source | Format | Voxel dimension | Orientation | Ontology | RDS | Citation | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ABA Human Half | _H. sapiens_ | [source](https://community.brain-map.org/t/allen-human-reference-atlas-3d-2020-new/405) | NIfTI | 500 &mu;m | RAS | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Ding et al.](https://pubmed.ncbi.nlm.nih.gov/27418273/) _J. Comp. Neurol._ 2016 |
| ABA Human Full | _H. sapiens_ | [source](https://community.brain-map.org/t/allen-human-reference-atlas-3d-2020-new/405) | NIfTI | 500 &mu;m | RAS | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Ding et al.](https://pubmed.ncbi.nlm.nih.gov/27418273/) _J. Comp. Neurol._ 2016 |
| Hammersmith Atlas | _H. sapiens_ | [source](http://brain-development.org/brain-atlases/adult-brain-atlases/adult-brain-maximum-probability-map-hammers-mith-atlas-n30r83-in-mni-space/) | NIfTI | 1 mm | RAS | custom ontology | [RDS]() | [Hammers et al.](http://www.ncbi.nlm.nih.gov/pubmed/12874777) _Human Brain Mapping_ 2003, [Gousias et al.](http://www.ncbi.nlm.nih.gov/pubmed/18234511) _NeuroImage_ 2007, [Faillenot et al.](http://doi.org/10.1016/j.neuroimage.2017.01.073) _NeuroImage_ 2017, [Wild et al.](http://doi.org/10.1371/journal.pone.0180866) PLoS ONE 2017, | 
| ABA Mouse CCF2017 | _M. musculus_ | [source](http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/) | NRRD | 50 &mu;m | PIR | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Lein et al.](https://www.nature.com/articles/nature05453) _Nature_ 2007 |
| Drosophila JRC2018 | _D. melanogaster_ | [source](https://www.janelia.org/open-science/jrc-2018-brain-templates) | NRRD | 38 &mu;m | RAS | custom ontology | [RDS]() |  [Bogovic et al.](https://journals.plos.ocoldcutsrg/plosone/article?id=10.1371/journal.pone.0236495) _PlOS ONE_ 2020, [Ito et al.](https://www.cell.com/neuron/fulltext/S0896-6273(13)01178-1) _Neuron_ 2014 |
| Chimpanzee Davi130 Juna.Chimp | _P. troglodytes _ | [source](https://www.chimpanzeebrain.org/node/2347) | NIfTI | 500 &mu;m | RAS | custom ontology | [RDS]() | [Vickery et al.](https://pubmed.ncbi.nlm.nih.gov/33226338/) _eLife_ 2020 |
