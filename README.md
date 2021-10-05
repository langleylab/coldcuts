[![](https://img.shields.io/badge/devel%20version-0.0.0.900-orange.svg)](https://github.com/langleylab/coldcuts)

# coldcuts

<img src = "https://user-images.githubusercontent.com/21171362/131806935-ed890016-a845-4274-8cb3-cd78c16aeb00.png" width=200>

**_disclaimer: this package is still in development, will have bugs and will undergo additional changes_**

## Introduction
**`coldcuts`** is an R package that allows you to **draw and plot automatically** segmentations from 3D voxel arrays. 

The name is inspired by one of Italy's best products.

Voxel arrays can come in different forms, either as NIfTI files, NRRD files, RAW files or directly imported in R as such.

These arrays host a value of 0 for "empty space", and another numerical value corresponding to a structural ID in any other position. 

Once the user has added the relevant information, such as orientation, units of measurement, voxel size, etc (all optional), the package "slices" the array along the 3 anatomical axes. For the time being, we will always assume that the orientation of the array is **RAS**, i.e. from Left to **R**ight (x axis), from Posterior to **A**nterior (y axis) and from Inferior to **S**uperior (z axis). If you are importing an array that was created in a different orientation, such as PIR or LPS, you can use an argument in the functions to reorient the array to RAS. RAS being the default orientation means that the anatomical planes will always be in this order: **sagittal**, **coronal** and **axial**. 

Once the array has been sliced, the package computes the polygon(s) for each structure, i.e. for every group of pixels in each slice bearing the same numerical value. Polygons are computed through a marching squares algorithm (as implemented in the `isoband` package) which assigns non-consecutive polygons to different groups. 

Since some of the polygons can have holes in a specific slice - as is the case, for instance, of the forebrain white matter in the human segmentation - the package functions automatically detect which polygons within the same structure are contained by bigger polygons, and treat them as holes in draw time.

Once the drawing procedure is finished, a `segmentation` object is created which can be plotted and accessed in different ways. 

## Installing the package

You can install this package using `devtools::install_github()`:

```{r}
devtools::install_github("langleylab/coldcuts")
```

## Drawing a segmentation

In this package we refer to 2 different types of operations: **drawing** is the procedure which effectively constructs polygons, whereas **plotting** consists in creating plots where polygons are selected and coloured according to user-defined rules.

Therefore the creation of a `segmentation` object is achieved by calling `drawSegmentation()`. 

To use this function we need to have access to 2 files: the voxel array and the ontology table. We then have to instruct the function on whether array rotations are required, and add other types of metadata on the fly.
In this example, **after cloning the git repository**, we use the `annotation.nii.gz` file from the Allen Institute for Brain Sciences, containing a 500&mu;-spaced voxel array with 141 annotated structures. We add the ontology file `allen_human_ontology.csv` as well, both of which can be found in the `data` folder of this repository.
We add the MNI152 reference space to the metadata. Since this NIfTI file has already some relevant information in its header, the function will add it for us.

```{r}
seg <- drawSegmentation(nifti_file = "./data/annotation.nii.gz", 
                        ontology_file = "./data/allen_human_ontology.csv",
                        reference_space = "MNI152",
                        verbose = FALSE)
                        
seg

Segmentation from file:  ./data/annotation.nii.gz 
                    Plane(s):  sagittal, coronal, axial 
                  Directions:  RAS 
       Dimensions (original):  394 x 466 x 378 
      Dimensions (effective):  148 x 362 x 310 
            Pixel dimensions:  0.5 x 0.5 x 0.5 
             Dimension units:  mm 
             Reference space:  MNI152 
              Structures (n):  141 
             Structures (ID):  10368, 12805, 10557, 12369, 10460, 10602 and 135 more. 
       Structures (acronyms):  BL, 4V, FWM, Aq, Pin, 3V and 135 more. 
Type showCitation("name_of_segmentation") to display the citation(s) for this segmentation.
```

## Inspecting a segmentation ontology

We can access the `ontology` slot using the `ontology()` function: 

```{r}
head(ontology(seg))

         id atlas_id                       name acronym st_level ontology_id hemisphere_id weight parent_structure_id depth graph_id graph_order
10153 10153       NA               neural plate      NP       NA          11             3   8660                  NA     0       16           0
10154 10154       NA                neural tube      NT       NA          11             3   8660               10153     1       16           1
10155 10155       NA                      brain      Br       NA          11             3   8660               10154     2       16           2
10156 10156       NA forebrain (prosencephalon)       F       NA          11             3   8660               10155     3       16           3
10157 10157       NA   gray matter of forebrain     FGM       NA          11             3   8660               10156     4       16           4
10158 10158       NA              telencephalon     Tel       NA          11             3   8660               10157     5       16         657
                          structure_id_path color_hex_triplet neuro_name_structure_id neuro_name_structure_id_path failed sphinx_id structure_name_facet
10153                               /10153/            D7D8D8                      NA                           NA      f      6600           3041346888
10154                         /10153/10154/            D0D0D1                      NA                           NA      f      6601           1172862311
10155                   /10153/10154/10155/            E0E0E0                      NA                           NA      f      6602           3016132225
10156             /10153/10154/10155/10156/            EBD6D0                      NA                           NA      f      6603           2015023488
10157       /10153/10154/10155/10156/10157/            EBD6D0                      NA                           NA      f      6604            545890574
10158 /10153/10154/10155/10156/10157/10158/            EBD6D0                      NA                           NA      f      7257           2480649783
      failed_facet                  safe_name     col
10153    734881840               neural plate #D7D8D8
10154    734881840                neural tube #D0D0D1
10155    734881840                      brain #E0E0E0
10156    734881840 forebrain (prosencephalon) #EBD6D0
10157    734881840   gray matter of forebrain #EBD6D0
10158    734881840              telencephalon #EBD6D0
```

This ontology was downloaded from the Allen Brain Atlas API, and it contains several pieces of information that we will not need. However, we do need the `id`, `name`, `acronym`, `parent_structure_id`, `structure_id_path` and `col` fields.

Ontologies may contain additional IDs, which are part of how the authors of the segmentation have classified structures into higher order groupings (e.g. "hipothalamus" as a higher order grouping for several hipothalamic nuclei). This ontology can be visualized as a tree using `plotOntologyGraph()`:

```{r}
plotOntologyGraph(seg)
```

<img width="1010" alt="ontology_graph" src="https://user-images.githubusercontent.com/21171362/131797720-1f92d2ff-4b02-4671-a869-1787f8c160c1.png">

Acronyms can be compared to their names looking at the ontology.

## Plotting a segmentation

The basic functionality with no external datasets allows to plot all structures within specific slices. We pick 3 random slices with the `s_slice`, `c_slice` and `a_slice` arguments for sagittal, coronal and axial respectively. 

```{r}
plotSegmentation(seg, s_slice = 50, c_slice = 60, a_slice = 30)
```

<img width="1144" alt="slices_1" src="https://user-images.githubusercontent.com/21171362/131649114-0c9aeee4-4f54-46f9-b0ed-d440af18bbe3.png">

We can also choose only one axis, and we can show structural labels (as acronyms from the ontology table):

```{r}
plotSegmentation(seg, s_slice = 50, show_labels = TRUE)

```

<img width="717" alt="slices_2" src="https://user-images.githubusercontent.com/21171362/131649611-c4a10103-d1fb-4eed-9eec-6e9320341a11.png">

Additionally, we can subset a `segmentation` object by structure(s). Here is an example in which we isolate the Forebrain White Matter (FWM):

```{r}
seg_sub <- subsetByStructures(segmentation = seg, planes_chosen = "sagittal", structures = "10557")

plotSegmentation(seg_sub, s_slice = 50, show_labels = TRUE)
```

<img width="716" alt="slices_3" src="https://user-images.githubusercontent.com/21171362/131650045-91b1d568-f49f-498e-b6fa-059e1c0ea49a.png">

## Creating a `segmentationAssay`

The `segmentationAssay` S4 class was created to host an external dataset, adding numerical tables and other information to a segmentation. 
For instance, in the case of GTEx brain data, a `segmentationAssay` can host the median TPM count table (in the `values` slot), additional information on the samples (in the `sampledata` slot) and, most importantly, the mapping from samples to structures (`mapping` slot). 

Depending on the segmentation, the assay may require a one-to-many mapping of samples to structures. For instance, the GTEx brain dataset contains "caudate basal ganglia", which has to be mapped to 3 different structural IDs in the Allen Human Brain Atlas segmentation. 

This mapping has to be manually curated by the user, as there is no way to detect automatically this sort of correspondences. 

Here we will show how to create this mapping for the GTEx bulk RNA-seq data, using median gene-level TPM values. 

These values can be found [here](https://www.gtexportal.org/home/datasets), in particular the _GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz_ file (available in the `data` folder of this repository as `gtex_tpm.gct`).

```{r}
# read the table
gtex <- read.table("./data/gtex_tpm.gct", header = TRUE, sep = "\t")

# subsetting and wrangling names
gtex <- gtex[,c(1, 2, grep("Brain", colnames(gtex)))]
colnames(gtex) <- gsub(colnames(gtex), pattern = "Brain...", replacement = "")
colnames(gtex) <- gsub(colnames(gtex), pattern = "[.]", replacement = "_")
colnames(gtex) <- gsub(colnames(gtex), pattern = "__", replacement = "_")
colnames(gtex) <- gsub(colnames(gtex), pattern = "_$", replacement = "")

#removing duplicates (mostly 0 TPM genes)
gtex <- gtex[!duplicated(gtex$Description),]
rownames(gtex) <- gtex$Description
```

Once the table is subsetted and readable, we map samples to structural IDs creating a list:

```{r}
# create mapping
brain_regions <- list()

brain_regions[["Amygdala"]] <- ontology(seg)[grep(pattern  ="amyg",
                                                         x = ontology(test2)$name), "id"]

brain_regions[["Cortex"]] <- ontology(seg)[grep(pattern  = "frontal",
                                                       x = ontology(seg)$name), "id"]

brain_regions[["Frontal_Cortex_BA9"]] <- ontology(seg)[grep(pattern  = "cingulate gyrus, rostral",
                                                                   x = ontology(seg)$name), "id"]

brain_regions[["Caudate_basal_ganglia"]] <- ontology(seg)[grep(pattern  = "caudate",
                                                                      x = ontology(seg)$name), "id"]

brain_regions[["Putamen_basal_ganglia"]] <- ontology(seg)[grep(pattern  = "putamen",
                                                                      x = ontology(seg)$name), "id"]

brain_regions[["Hypothalamus"]] <- ontology(seg)$id[grep(ontology(seg)[which(ontology(seg)$acronym == "HTH"),"id"], ontology(seg)$parent_structure_id)]

brain_regions[["Substantia_nigra"]] <- ontology(seg)$id[grep(pattern  = "SN", x = ontology(seg)$acronym)]

brain_regions[["Hippocampus"]] <- ontology(seg)[grep(pattern  = "hipp",
                                                            x = ontology(seg)$name), "id"]

brain_regions[["Nucleus_accumbens_basal_ganglia"]] <- ontology(seg)[grep(pattern  = "accumbens",
                                                                                x = ontology(seg)$name), "id"]

brain_regions[["Cerebellum"]] <- ontology(seg)[grep(pattern  = "cerebell",
                                                           x = ontology(seg)$name), "id"]

brain_regions[["Hippocampus"]] <- brain_regions[["Hippocampus"]][2:length(brain_regions[["Hippocampus"]])]

brain_regions[["Hippocampus"]] <- brain_regions[["Hippocampus"]][brain_regions[["Hippocampus"]] != "10377"]
```

We now have everything we need to create a `segmentationAssay` object and add it to our segmentation:

```{r}
gtexAssay <- new("segmentationAssay",
                 values = as.matrix(gtex[,3:ncol(gtex)]),
                 mapping = brain_regions)

# placeholder until the getter and setter functions are written
seg@assays <- list("gtex" = gtexAssay)

```

## Creating a maximum projection and plotting external data

When plotting data such as expression values in the segmentation, not all structures will be visible within the same slice, making the choice of a specific slice hard. For this reason we can create a maximum projection of structure slices on every plane, in both directions: slice-level polygons for every structure are joined together into single polygons and displayed in the order that they appear from both points of view, e.g. in the sagittal plane, from left to right (LR) and from right to left (RL). 

Structures can be subset for maximum projections, since using all structures may just result in producing a side view of the segmentation (which may be a desired behaviour in some cases). Since maximum projections are specifically important for plotting datasets, we create a projection that has the same name as the dataset using the `addMaxProjection()` function, only for the sagittal plane:

```{r}

structures_chosen <- unlist(seg@assays$gtex@mapping)
structures_chosen <- intersect(structures_chosen, metaData(seg)$structures)

seg <- addMaxProjection(name = "gtex", 
                        segmentation = seg, 
                        structures = structures_chosen, 
                        planes_chosen = "sagittal")
```
At this point we can plot data from one of our assays (GTEx) in the sagittal projection, using `plotBrainMap()`, specifying the `assay` name and the gene (`feature` argument):

```{r}
plotBrainMap(segmentation = seg,
             assay = "gtex",
             feature = "APP")
```

<img width="1033" alt="gtex_map" src="https://user-images.githubusercontent.com/21171362/131799896-cdede3a8-fb22-4558-906f-45db8eef2b13.png">

## Finding additional segmentations (WIP)

Many institutions have published their own segmentations as NIfTI, NRRD or RAW files.
In this table we curate a few segmentations including their source(s).

| Name | Organism | Source | Format | Voxel dimension | Orientation | Ontology | RDS | Citation | 
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ABA Human Half | _H. sapiens_ | [source](https://community.brain-map.org/t/allen-human-reference-atlas-3d-2020-new/405) | NIfTI | 500 &mu;m | RAS | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Ding et al.](https://pubmed.ncbi.nlm.nih.gov/27418273/) _J. Comp. Neurol._ 2016 |
| ABA Human Full | _H. sapiens_ | [source](https://community.brain-map.org/t/allen-human-reference-atlas-3d-2020-new/405) | NIfTI | 500 &mu;m | RAS | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Ding et al.](https://pubmed.ncbi.nlm.nih.gov/27418273/) _J. Comp. Neurol._ 2016 |
| ABA Mouse CCF2017 | _M. musculus_ | [source](http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_masks/) | NRRD | 50 &mu;m | PIR | [ontology](http://help.brain-map.org/display/api/Downloading+an+Ontology%27s+Structure+Graph) | [RDS]() | [Lein et al.](https://www.nature.com/articles/nature05453) _Nature_ 2007 |
| Drosophila JRC2018 | _D. melanogaster_ | [source](https://www.janelia.org/open-science/jrc-2018-brain-templates) | NRRD | 38 &mu;m | RAS | custom ontology | [RDS]() |  [Bogovic et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236495) _PlOS ONE_ 2020, [Ito et al.](https://www.cell.com/neuron/fulltext/S0896-6273(13)01178-1) _Neuron_ 2014 |


# Additional details
## The `segmentation` class and its components

The main workhorse of this package is its own S4 class, `segmentation`, which hosts many components necessary for plotting and integrating with external datasets.

<img width="844" alt="segmentation_class" src="https://user-images.githubusercontent.com/21171362/131645071-e3942aeb-fca8-426a-be9f-45a771904796.png">


`segmentation` class objects hold the following slots:
- **slices** is a list with 3 nested lists, named **sagittal**, **coronal** and **axial** for the three anatomical planes. Every anatomical plane is a list of slices along the plane spaced by 1 voxel (SL1...SLN in the figure), every slice is a list of structures (ST1...STN in the figure), with every structure being a list of `segPointSet` objects, another S4 class that holds polygon information in a lightweight format (P1...PN). 
- **ontology** is a data frame containing all relevant information for structures: their IDs, names, acronyms, hierarchical relationship if present, colour codes etc. 
- **projections** contains maximum projections for a series of structures in each plane. More than one maximum projection can be hosted in each `segmentation` object to accommodate for more assays. Every projection can be differentiated by the subset of structures that are being projected, which is usually done to match between an assay (external dataset) and a segmentation.
- **assays** is a list of `segmentationAssay` objects, another S4 class that holds numerical `values` - such as gene counts - sample metadata (`sampledata`) and most importantly a one-to-many `mapping` of structures to the relevant samples. 
- **metadata** is a list of various types of information regarding the segmentation, such as the file name, orientation, pixel/voxel dimension, original citation, reference space, etc. Most of these are manually compiled by the user, although an attempt is made at extracting some of them from a file containing this information in its header (e.g. NiFTi).
- **structure_tables** is a list of data.tables, one per anatomical axis, which is used by some internal functions to check which slices contain which structures.


## TODO
- ontology palettes
- other species and access to their segmentations
