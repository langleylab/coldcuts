# coldcuts

## Introduction
This is an R package that allows you to **draw and plot automatically** segmentations from 3D voxel arrays. 

Voxel arrays can come in different forms, either as NIfTI files, NRRD files, RAW files or directly imported in R as such.

These arrays host a value of 0 for "empty space", and another numerical value corresponding to a structural ID in any other position. 

Once the user has added the relevant information, such as orientation, units of measurement, voxel size, etc (all optional), the package "slices" the array along the 3 anatomical axes. For the time being, we will always assume that the orientation of the array is **RAS**, i.e. from Left to **R**ight (x axis), from Posterior to **A**nterior (y axis) and from Inferior to **S**uperior (z axis). If you are importing an array that was created in a different orientation, such as PIR or LPS, you can use an argument in the functions to reorient the array to RAS. RAS being the default orientation means that the anatomical planes will always be in this order: **sagittal**, **coronal** and **axial**. 

Once the array has been sliced, the package computes the polygon(s) for each structure, i.e. for every group of pixels in each slice bearing the same numerical value. Polygons are computed through a marching squares algorithm (as implemented in the `isoband` package) which assigns non-consecutive polygons to different groups. 

Since some of the polygons can have holes in a specific slice - as is the case, for instance, of the forebrain white matter in the human segmentation - the package functions automatically detect which polygons within the same structure are contained by bigger polygons, and treat them as holes in draw time.

Once the drawing procedure is finished, a `segmentation` object is created which can be plotted and accessed in different ways. 


## Drawing a segmentation

In this package we refer to 2 different types of operations: **drawing** is the procedure which effectively constructs polygons, whereas **plotting** consists in creating plots where polygons are selected and coloured according to user-defined rules.

Therefore the creation of a `segmentation` object is achieved by calling `drawSegmentation()`. 

To use this function we need to have access to 2 files: the voxel array and the ontology table. We then have to instruct the function on whether array rotations are required, and add other types of metadata on the fly.
In this example we use the `annotation.nii.gz` file from the Allen Institute for Brain Sciences, containing a 500-micron-spaced voxel array with 141 annotated structures. We add the ontology file `allen_human_ontology.csv` as well, and we add the MNI152 reference space to the metadata. Since this NIfTI file has already some relevant information in its header, the function will add it for us.

```{r}
seg <- drawSegmentation(nifti_file = "annotation.nii.gz", 
                        ontology_file = "allen_human_ontology.csv",
                        reference_space = "MNI152")
                        
seg

Segmentation from file:  annotation.nii.gz 
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

## The `segmentation` class and its components

The main workhorse of this package is its own S4 class, `segmentation`, which hosts many components necessary for plotting and integrating with external datasets.

<img width="844" alt="segmentation_class" src="https://user-images.githubusercontent.com/21171362/131645071-e3942aeb-fca8-426a-be9f-45a771904796.png">


`segmentation` class objects hold the following slots:
- **slices** is a list with 3 nested lists, named **sagittal**, **coronal** and **axial** for the three anatomical planes. Every anatomical plane is a list of slices along the plane spaced by 1 voxel (SL1...SLN in the figure), every slice is a list of structures (ST1...STN in the figure), with every structure being a list of `segPointSet` objects, another S4 class that holds polygon information in a lightweight format (P1...PN). 
- **ontology** is a data frame containing all relevant information for structures: their IDs, names, acronyms, hierarchical relationship if present, colour codes etc. 
- **projections** contains maximum projections for a series of structures in each plane. More than one maximum projection can be hosted in each `segmentation` object to accommodate for more assays. Every projection can be differentiated by the subset of structures that are being projected, which is usually done to match between an assay (external dataset) and a segmentation.
- **assays** is a list of `segmentationAssay` objects, another S4 class that holds numerical `values` - such as gene counts - sample metadata (`sampledata`) and most importantly a one-to-many `mapping` of structures to the relevant samples. For instance, the GTEx brain dataset contains "Caudate" which has to be mapped to 3 different structural IDs in the Allen Human Brain Atlas segmentation. 
- **metadata** is a list of various types of information regarding the segmentation, such as the file name, orientation, pixel/voxel dimension, original citation, reference space, etc. Most of these are manually compiled by the user, although an attempt is made at extracting some of them from a file containing this information in its header (e.g. NiFTi).
- **structure_tables** is a list of data.tables, one per anatomical axis, which is used by some internal functions to check which slices contain which structures.
