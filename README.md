# CytoCompare: an R Package for Computational Comparisons of Cytometry Profiles



Cytometry is an experimental technique used to measure molecules expressed by cells at a single cell resolution.
Recently, several technological improvements have made possible to greatly increase the number of cell markers that can be simultaneously measured.
Many computational methods have been proposed to identify clusters of cells having similar phenotypes.

Nevertheless, no computational method exists to compare the phenotypes of cell clusters identified by different clustering approaches.
These phenotypic comparisons are necessary to choose the appropriate clustering methods and settings.
Because of this lack of computational tools, comparisons of cell cluster phenotypes are often performed manually, a highly biased and time-consuming process.

CytoCompare is an R package that performs comparisons between the phenotypes of cell clusters with the purpose of identifying similar ones.
For each comparison of two cell clusters, CytoCompare provides a similarity measure as well as a p-value asserting the statistical significance of the similarity.
CytoCompare can generate parallel coordinates, parallel heatmaps, multidimensional scaling or circular graph representations to visualize easily cell cluster phenotypes and the comparison results.
CytoCompare can import clustering results from various algorithms including SPADE, viSNE/ACCENSE, and Citrus, the most current widely used algorithms.


# Table of Contents
1. [Package overview](#package_overview)
2. [Package installation](#package_installation)
3. [Importation of cytometry cell clusters](#object_importation)
    3. [Importation of cell clusters from viSNE/ACCENSE results](#visne_importation)
    4. [Importation of cell clusters from SPADE results](#spade_importation)
    5. [Importation of cell clusters from CITRUS results](#citrus_importation)
    6. [Importation of cell clusters from text files](#cluster_importation)
4. [Manipulation of CLUSTER objects](#object_manipulation)
    1. [Extraction of CLUSTER objects](#object_extraction)
    2. [Combination of CLUSTER objects](#object_combination)
5. [Representations of cytometry profiles](#object_representation)
    1. [Visualization of cluster profiles](#cluster_visualization)
    2. [Pairwise visualization of cytometry profiles](#pairwise_visualization)
6. [Comparisons of cytometry profiles](#object_comparison)
7. [Manipulation of comparison results](#res_object)
    1. [Structure of comparison results](#res_structure)
    2. [Summarization of comparison results](#res_summarization)
    3. [Extraction of comparison results](#res_extraction)
    4. [Combination of comparison results](#res_combination)
    5. [Visualization of single comparison results](#res_visualization)
    6. [Visualization of multiple comparison results](#res_d3js)
        1. [Circular graph visualization](#res_d3js_circular)
        2. [MDS visualization](#res_d3js_mds)
8. [Miscellaneous functions](#miscellaneous)
    1. [Biplot representations of cell profiles](#biplot)
    2. [Density heatmap representations of cell cluster profiles](#dheatmap)
9. [Template of the compare() function](#compare_template)
10. [Structures of the CytoCompare objects](#object_structure)
    1. [Overview of CytoCompare objects](#object_structure_uml)
    2. [Structure of CELL object](#cell_structure)
    3. [Structure of CLUSTER object](#cluster_structure)
    4. [Structure of GATE object](#gate_structure)
    5. [Structure of MWEIGHTS object](#mweights_object)
        1. [Structure of MWEIGHTS object](#mweights_structure)
        2. [Summarization of MWEIGHTS object](#mweights_summarization)
        3. [Extraction of MWEIGHTS object](#mweights_extraction)
        4. [Combination of MWEIGHTS object](#mweights_set)
        5. [Visualization of MWEIGHTS object](#mweights_visualization)
    6. [Structure of DENSITY object](#density_object)
        1. [Structure of DENSITY object](#density_structure)
        2. [Summarization of DENSITY object](#density_summarization)
        3. [Visualization of DENSITY object](#density_visualization)
11. [License](#license)
12. [References](#references)

# <a name="package_overview"></a> 1. Package overview
Flow and mass cytometry are experimental techniques used to characterize cells at a single-cell level.
Both techniques use labeled antibodies to measure cell marker expressions.
Flow cytometry uses antibodies conjugated with fluorochromes to quantify stained cells with a laser detection system.
The recently introduced mass cytometry ([CyTOF](https://www.fluidigm.com/products/cytof) [[1](http://www.ncbi.nlm.nih.gov/pubmed/21551058)]) uses antibodies conjugated with metals to quantify stained cells with a mass spectrometer.
While flow cytometry can currently handle up to 18 fluorochromes, mass cytometry can currently handle up to 40 metals.
Mass cytometry offers important perspectives as it can potentially evaluate more than 100 metals.
Another recently introduced technique, called hyper-spectral cytometry [[2](http://www.ncbi.nlm.nih.gov/pubmed/24271566)], combines ultrafast optical spectroscopy with flow cytometry and can also increase up to 40 the number of simultaneously usable fluorochromes.
Such significant technological improvements make flow and mass cytometry suitable for massive data-mining and bioinformatics developments, in order to better explore and analyze complex cell interaction systems.

Cytometry data can be manually analyzed using a gating strategy (hierarchical gating usually followed by Boolean gating) or using automatic gating/clustering algorithms.
In both approaches, the aim is to identify cell populations having similar phenotypes, based on the expressions of selected markers.
Hierarchical gating is done using biplot representations, where each axis represents the intensity of a marker of interest and where dots represent the cells.
Through iterative steps, operators select cell subsets and explore their phenotypes.
Boolean gating is done by quantifying cells present within all possible combinations of gates for a population of interest.
Even if Boolean gating can be done in an automatic way, gates still have to be delineated for each cell marker.
FlowJo [[3](http://www.flowjo.com)] and CytoBank [[4](https://www.cytobank.org)] are among the most common software for manual gating.
Both hierarchical and Boolean gating can be very fastidious and biased. 
The difficulty mainly comes from the high dimensionality of data, that is to say from the high number of markers and cell populations.
This task becomes even more complicated as the number of markers increases or as the studies become more complex. 
On the other hand, automatic gating methods (such as [viSNE](https://www.c2b2.columbia.edu/danapeerlab/html/cyt.html) [[5](http://www.ncbi.nlm.nih.gov/pubmed/23685480)] /[ACCENSE](http://www.cellaccense.com/) [[6](http://www.ncbi.nlm.nih.gov/pubmed/24344260)], [SPADE](http://cytospade.org/) [[7](http://www.ncbi.nlm.nih.gov/pubmed/21964415)] or [Citrus](https://github.com/nolanlab/citrus/wiki) [[8](http://www.ncbi.nlm.nih.gov/pubmed/24979804)]) use various algorithms to computationally identify cell populations having similar phenotypes and provides then less biased results. 

Once identified cell populations, also named cell clusters, need to be deeper characterized or associated with known phenotypes for further investigations. 

While all these tools have been developed, no computational strategy currently exists to compare the phenotypes of identified cell clusters.
These phenotypic comparisons are crucial to explore cell clustering results obtained from different computational approaches and to recognize cell populations with common or marginal phenotypes.
Thereby, such comparisons are often performed manually through multiple comparisons of marker expression densities, a highly biased and time-consuming process.
Indeed, comparisons of cell cluster phenotypes based on the mean or median of marker expressions can be biased as they do not consider all the characteristics of marker expression distributions.
Proper phenotypic comparisons of cell marker expressions must be based on marker expression density distributions.
Additionally, no visualization methods have been proposed to represent such phenotypic comparisons in easily interpretable ways.

Cytometry experts can face at least five different situations that can necessitate comparing the cell clusters based on their phenotypes:

* The first situation occurs when selecting the cell clustering algorithm to perform the analysis, which implies to compare clustering results from different algorithms. 
* The second situation happens when selecting the settings for one given automatic gating algorithm, which implies to compare clustering results obtained using different parameters. 
* The third situation takes place when ascertaining the reproducibility of the clustering results by repeating the analysis using the same settings. 
* The fourth situation arises when users aim to organize the identified cell clusters based on their phenotypes. 
* The fifth situation occurs when comparing the phenotypes of cell populations from different biological studies.

CytoCompare is an R package that performs comparisons between the phenotypes of cell clusters with the purpose of identifying similar ones.
For each comparison of two cell clusters, CytoCompare provides a similarity measure as well as a p-value asserting the statistical significance of the similarity.
CytoCompare can generate parallel coordinates, parallel heatmaps, multidimensional scaling or circular graph representations to visualize easily cell cluster phenotypes and the comparison results.
CytoCompare can import clustering results from various algorithms including viSNE/ACCENSE, SPADE, and Citrus, the most current widely used algorithms.

We designed CytoCompare in a way that it can be easily used by non-bioinformatician experts, but can also be easily customizable by users with more expertise in bioinformatics.
Through the multiple cell cluster manipulation methods, CytoCompare is also a powerful analysis pipeline for high-dimensional cytometry data.


# <a name="package_installation"></a> 2. Package installation
The `ggplot2`, `ggrepel`, `grid`, `igraph`, `MASS`, `RJSONIO`, and `XML` R packages as well as the `flowCore` and `flowUtils` Bioconductor packages are required for running CytoCompare. These packages can be installed using the following commands:

```r
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('grid')
install.packages('igraph')
install.packages('MASS')
install.packages('RJSONIO')
install.packages('XML')

source("http://bioconductor.org/biocLite.R")
biocLite(suppressUpdates=TRUE)
biocLite("flowCore",suppressUpdates=TRUE)
biocLite("flowUtils",suppressUpdates=TRUE)
```

CytoCompare is available on [GitHub](https://github.com/), at https://github.com/tchitchek-lab/CytoCompare. Its installation can be done via the `devtools` package using the following commands:

```r
install.packages('devtools')
library("devtools")
install_github('tchitchek-lab/CytoCompare')
```

Once installed, CytoCompare can be loaded using the following command:

```r
library("CytoCompare")
```



An example dataset, obtained from public mass cytometry data [[9](http://www.ncbi.nlm.nih.gov/pubmed/21964415)], is available in CytoCompare.
This example dataset consists of three healthy human bone marrow samples, unstimulated or stimulated by BCR-inductor or IL-7, and measured using a mass cytometry panel of more than 30 cell markers. 
This panel has been designed to identify a large spectrum of immune cell types. 
A SPADE analysis has been performed to identify cell clusters. 
Then, these cell clusters have been manually labeled based on their marker expressions [[9](http://www.ncbi.nlm.nih.gov/pubmed/21964415)]. 
SPADE cell clusters corresponding to six majors cell types (B, CD33+ monocytes, naive or memory CD4+ T cells, and naive or memory CD8+ T cells) have been extracted to create this example dataset.

The raw data and the R object (`CytoCompareExample.rdata`) corresponding to this example dataset are available on a public ftp server: [ftp://ftp.cytocompare.org/public/](ftp://cytocompare:cytocompare@ftp.cytocompare.org/public/) (username: cytocompare, password: cytocompare).

The R object containing the example dataset can be directly downloaded and loaded in the current R session using the following command:



```r
 # downloads the 'CytoCompareExample.rdata' file on the current folder, if no such file exists, and load it into the current R session
load.examples()
```

Once downloaded, the following objects will be available:

* `bm_example.clusters`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters for all the different cell populations
* `bm_example.clusters.b`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters associated to the B cell populations
* `bm_example.clusters.mono`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters associated to monocyte cell populations
* `bm_example.clusters.tCD4naive`, a `CLUSTER` object containing the phenotypes SPADE of cell clusters associated to the naive CD4+ T cell populations
* `bm_example.clusters.tCD8naive`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters associated to the naive CD8+ T cell populations identified by SPADE
* `bm_example.clusters.tCD4mem`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters associated to memory CD4+ T cell populations
* `bm_example.clusters.tCD8mem`, a `CLUSTER` object containing the phenotypes of SPADE cell clusters associated to memory CD8+ T cell populations

# <a name="object_importation"></a> 3. Importation of cytometry cell clusters 

## <a name="visne_importation"></a> 3.1. Importation of cell clusters from viSNE/ACCENSE results
The `import.VISNE_ACCENSE()` function imports cell clusters from viSNE/ACCENSE results into a `CLUSTER` object. 
The viSNE algorithm is a visualization algorithm for high-dimensional cytometry data [[7](http://www.ncbi.nlm.nih.gov/pubmed/23685480)] that aims to represent cell in a 2-dimensional space using a dimensionality reduction process.
For each cell profile a 2-dimensional coordinate is calculated.
The 2 dimensions are named tSNE1 and tSNE2 and the results are presented via FCS files having two extra markers corresponding to these dimensions. 
The ACCENSE algorithm combine the viSNE results with a k-means clustering algorithm,  for both visualization and automated cell classification into subpopulations. 

An import of a viSNE/ACCENSE results can be done using the following command:


```r
# imports viSNE/ACCENSE results containing several cell clusters 
imported.visne <- import.VISNE_ACCENSE('./CLUSTER.viSNE/accense_output.csv')
```

Some parameters can be specified to apply numeric transformations on the marker expression values, or to exclude some markers in the import or transformation procedures:

* the `exclude` parameter is a character vector containing the marker names to be excluded in the import procedure
* the `trans` parameter is a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation 
* the `trans.para` parameter is a named list containing parameters for the transformation
* the `trans.exclude` parameter is a character vector containing the marker names for which no transformation must be applied on


The available transformation functions are (`trans` parameter):

* `"arcsinh"` for an arc sine hyperbolic transformation of the data (default choice)
* `"log"` for a logarithmic transformation of the data
* `"none"` for no transformation of the data

The available transformation parameters are (`trans.para` parameter):

* for the `"arcsinh"` transformation, the scale (cofactor) can be specified using the `arcsinh.scale` parameter
* for the `"log"` transformation, the base and the shift can be specified using the `log.base` and `log.shift` parameters. The 'log.shift' parameter allows the value "auto" which automatically identify the log shift avoiding to apply log transformations on negative values.


A summary of the `CLUSTER` object, imported by the `import.VISNE_ACCENSE()` function, can be done using the following command:

```r
print(imported.visne)
```

```
## Object class: CLUSTER
## Object name: cell_cluster
## Number of clusters profiles: 6
## Markers: 
## Time
## Cell Length
## 191-DNA
## 193-DNA
## 103-Viability
## 115-CD45
## 110-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 139-CD45RA
## 141-pPLCgamma2
## 142-CD19
## 144-CD11b
## *145-CD4
## *146-CD8
## *148-CD34
## *150-pSTAT5
## *147-CD20
## *152-Ki67
## *154-pSHP2
## 151-pERK1/2
## 153-pMAPKAPK2
## 156-pZAP70/Syk
## 158-CD33
## 160-CD123
## 159-pSTAT3
## 164-pSLP-76
## 165-pNFkB
## 166-IkBalpha
## 167-CD38
## 168-pH3
## 170-CD90
## 169-pP38
## 171-pBtk/Itk
## 172-pS6
## 174-pSrcFK
## 176-pCREB
## 175-pCrkL
## 110_114-CD3
## EventNum
## density
## Number of markers: 42
## Clustering markers: 
## none
## Number of clustering markers: 0
## Density bin width: 0.05
## Cluster profile names and number of associated cells:
##  1: 617 cells
##  2: 752 cells
##  3: 456 cells
##  4: 573 cells
##  5: 373 cells
## and 1 more...
```


## <a name="spade_importation"></a> 3.2. Importation of cell clusters from SPADE results
The `import.SPADE()` function imports cell clusters identified by the SPADE algorithm into a `CLUSTER` object. 
The Spanning Tree Progression of Density Normalized Events ([SPADE](http://cytospade.org/) [5]) algorithm is a visualization and analysis algorithm for high-dimensional cytometry data, available in [R](https://github.com/nolanlab/spade/) or in [Cytobank](https://www.cytobank.org/) [4].
SPADE was designed to analyze mass cytometry data but can also handle flow cytometry data.

SPADE identifies clusters of cells having similar expression for selected markers using an agglomerative hierarchical clustering-based algorithm combined with a density-based down-sampling procedure. 
SPADE represents its clustering results using a tree representation where each node represents a cell cluster and where similar cell clusters are linked using a minimal spanning tree approach.
In SPADE tree representations, the sizes of the nodes are proportional to the amount of cells present within each cluster. 
In order to characterize cell populations and to associate them with known phenotypes, nodes can be gradient-colored based on theirs mean expression intensities for a specific marker.

An import of SPADE clustering results can be done using the following command:


```r
# imports a SPADE folder containing several cell clusters 
imported.spade <- import.SPADE('./CLUSTER.SPADE/')
```

Some parameters can be specified to apply numeric transformations on the marker expression values, to exclude some markers in the import or transformation procedures, or to specify how to extract a SPADE result archive:

* the `exclude` parameter is a character vector containing the marker names to be excluded in the import procedure
* the `trans` parameter is a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation 
* the `trans.para` parameter is a named list containing parameters for the transformation
* the `trans.exclude` parameter is a character vector containing the marker names for which no transformation must be applied on
* the `bin.width` parameter is a numeric value indicating the width of the bins for the marker expression densities computations (0.05 by default)
* the `extract.folder` parameter is a folder path for extracting the SPADE zip archive (temporary folder by default)
* the `extract.folder.del` parameter is a logical value indicating if the extracted SPADE results should be removed after the extraction
* the `zip` parameter is a logical value that specify if the path indicates a zip file

The available transformation functions are (`trans` parameter):

* `"arcsinh"` for an arc sine hyperbolic transformation of the data (default choice)
* `"log"` for a logarithmic transformation of the data
* `"none"` for no transformation of the data

The available transformation parameters are (`trans.para` parameter):

* for the `"arcsinh"` transformation, the scale (cofactor) can be specified using the `arcsinh.scale` parameter
* for the `"log"` transformation, the base and the shift can be specified using the `log.base` and `log.shift` parameters. The 'log.shift' parameter allows the value "auto" which automatically identify the log shift avoiding to apply log transformations on negative values.


A summary of the `CLUSTER` object, imported by the `import.SPADE()` function, can be done using the following command:

```r
#print(imported.spade)
```


## <a name="citrus_importation"></a> 3.3. Importation of cell clusters from CITRUS results
The `import.CITRUS()` function imports cell cluster identified by the Citrus algorithm into a `CLUSTER` object. 
[Citrus](https://github.com/nolanlab/citrus/wiki) [[6](http://www.ncbi.nlm.nih.gov/pubmed/24979804)] is an algorithm developed similarly to SPADE, which can furthermore identify cell clusters associated with different biological condition phenotypes.

An import of Citrus clustering results can be done using the following command:


```r
# imports a Cirtus result file containing several cell clusters 
imported.citrus <- import.CITRUS('./CLUSTER.Citrus/citrusClustering.rData')
```

Some parameters can be specified to apply numeric transformations on the marker expression values, to exclude some markers in the import or transformation procedures, or to specify clusters to import:

* the `exclude` parameter is a vector containing the marker names to be excluded in the import procedure
* the `bin.width` parameter is a numeric value indicating the width of the bins for the marker expression densities computations (0.05 by default)
* the `minimumClusterSizePercent` parameter is a numeric value indicating the minimal ratio of cells per cluster to import
* the `cluster.selection` parameter is a character vector containing the names of the clusters to import

A summary of the `CLUSTER` object, imported by the `import.CITRUS()` function, can be done using the following command:

```r
print(imported.citrus)
```


## <a name="cluster_importation"></a> 3.4. Importation of cell clusters from text files
The `import.CLUSTER()` function imports one or several cell clusters from a tab separated file into a `CLUSTER` object.
In this case, the marker expressions of each cluster are assumed to be normally distributed.
Therefore, the tab separated file must contain for each cell cluster the means and the standard deviations of the expression markers and must be formatted as the following:

* each row must represent a cell cluster
* each column must represent a marker
* each cell in the table must contain the marker expression means and the standard deviations for a given cell cluster separated by a semicolon

The first column must contain the cell cluster names and the first row must contain the marker names.

A typical tab separated file to import must look like the following:

|        | marker_1 | marker_2 | marker_3 | marker_4 | marker_5 | marker_i | marker_n |
|--------|---------|---------|---------|---------|---------|---------|---------|
| cluster_1 | 5.23;0.94 | 0.42;0.88 | 4.26;0.33 | 3.23;2.73 | 1.57;1.69 | 3.14;1.71 | 6.01;2.40 | 
| cluster_2 | 2.56;1.94 | 2.34;2.93 | 5.45;1.88 | 2.34;0.97 | 1.27;0.97 | 1.70;0.48 | 6.32;1.67 | 
| cluster_3 | 2.67;0.36 | 4.10;0.74 | 4.56;2.80 | 3.56;1.29 | 0.89;0.53 | 0.42;0.54 | 4.32;1.39 | 
| cluster_4 | 4.33;2.76 | 4.32;2.76 | 3.89;0.55 | 3.45;0.76 | 3.45;0.84 | 5.07;0.28 | 5.65;2.28 | 
| cluster_5 | 5.24;1.75 | 2.76;0.63 | 5.08;2.17 | 2.56;0.93 | 5.23;1.99 | 4.98;1.30 | 7.98;0.67 | 


An import of a such file can be done using the following command:


```r
# imports a tab separated file containing several cell clusters 
imported.clusters <- import.CLUSTER('./CLUSTER.TXT/clusters.txt')
```

The `exclude` parameter allows to specify some markers to exclude in the import procedure. This parameter takes a character vector containing the marker names to be excluded.

A summary of the `CLUSTER` object, imported by the `import.CLUSTER()` function, can be done using the following command:

```r
print(imported.clusters)
```

```
## Object class: CLUSTER
## Object name: clusters.txt
## Number of clusters profiles: 5
## Markers: 
## marker_1
## marker_2
## marker_3
## marker_4
## marker_5
## marker_x
## marker_n
## Number of markers: 7
## Clustering markers: 
## none
## Number of clustering markers: 0
## No Densities associated to marker expressions
## Cluster profile names and number of associated cells:
##  cluster_1: no associated number of cells
##  cluster_2: no associated number of cells
##  cluster_3: no associated number of cells
##  cluster_4: no associated number of cells
##  cluster_5: no associated number of cells
```

It is to note that `CLUSTER` objects constructed via the `import.CLUSTER()` function do not contain the densities of expression markers.


# <a name="object_manipulation"></a> 4. Manipulation of CLUSTER objects

## <a name="object_extraction"></a> 4.1. Extraction of CLUSTER objects
The extract function `[i,j]` can be used to extract subsets of a `CLUSTER` object. The parameter `i` represents a vector of cell cluster to extract and the parameter `j` represents a vector of markers to extract. Both `i` and `j` can be numeric, logical or character vectors.

For example, subsets of a `CLUSTER` object can be extracted using the following commands:

```r
# extracts some specific cell clusters from a given CLUSTER object
monoclusters.subset1 <- bm_example.clusters.mono[c(1,3,5,7,11)]
# extracts the first 10 cell clusters and extracts some markers from a given CLUSTER object 
monoclusters.subset2 <- bm_example.clusters.mono[1:10,5:10]
```


## <a name="object_combination"></a> 4.2. Combination of CLUSTER objects
The combine function `c()` can be used to combine two or several `CLUSTER` objects.

This function is especially useful when combining cell clusters obtained from different cell clustering results into one single `CLUSTER` object. 

A combination of two or more `CLUSTER` objects can be done using the following commands:

```r
# combines a set of CLUSTER objects into a single CLUSTER object
tcells.combine1 <- c(bm_example.clusters.tCD4mem,
                           bm_example.clusters.tCD4naive,
                           bm_example.clusters.tCD8mem,
                           bm_example.clusters.tCD8naive)

# combines a list of CLUSTER objects into a single CLUSTER object
tcells.list <- list(bm_example.clusters.tCD4mem,
                 bm_example.clusters.tCD4naive,
                 bm_example.clusters.tCD8mem,
                 bm_example.clusters.tCD8naive)
tcells.combine2 <- do.call("c",tcells.list)
```

# <a name="object_representation"></a> 5. Representations of cytometry cell clusters

## <a name="cluster_visualization"></a> 5.1. Visualization of the phenotype of a cell cluster
The `plot()` function can be used to plot the cell cluster phenotype of a `CLUSTER` object, via parallel coordinates.
In such parallel coordinate representation the x-axis represents the different markers and the y-axis represents the marker expressions.
Error bars indicate the marker expression standard deviations.
Representations are generated using the `ggplot2` library [10] and can be modified by users.
 
A `CLUSTER` object, which contains one or several cell clusters, can be plotted using the following command:

```r
# plots the phenotype of the second cell cluster from a CLUSTER object
plot(bm_example.clusters[2])
```

<img src="README.figures/single_cluster_plot-1.png" style="display: block; margin: auto;" />

If the `CLUSTER` object contains several cell clusters, then all the different cell clusters will be plotted separately using the following command:

```r
# plots the phenotypes of first 3 cell clusters in a CLUSTER object
plot(bm_example.clusters[1:3])
```

<img src="README.figures/single_cluster_multiple_plot-1.png" style="display: block; margin: auto;" />

## <a name="pairwise_visualization"></a> 5.2. Pairwise visualization of cell clusters
Two cytometry cell clusters can be visualized in one single representation using parallel coordinates via the `plot()` function.
Such representations are especially useful to better visualize the marker expression similarities between two cell clusters. 



```r
# plots two cell clusters in one single parallel coordinate representation
plot(bm_example.clusters[1],bm_example.clusters[2])
```

<img src="README.figures/cluster_cluster_plot-1.png" style="display: block; margin: auto;" />



# <a name="object_comparison"></a> 6. Comparisons of cytometry cell clusters 

The phenotypes of cell cluster contained in `CLUSTER` objects can be compared using the `compare()` function. 
Comparison results are stored in a `RES` object (please refer to section 7 for more details). 

For each phenotypic comparison of two cell clusters, a similairty measure is provided as well as a p-value asserting the statistical significance of the similiarity. 

Comparisons can be performed based on the whole set of common markers between the two profiles, or based on a subset of markers specified by the user. Moreover, markers can be weighted in the comparison procedure, via a `MWEIGHTS` object.

When comparing two cell clusters, CytoCompare first computes a similarity measure for each pair of cell markers (step 1). 
In our statistical approach, the Kolmogorov-Smirnov distance is used to qualify the similarity between marker distribution densities [22]. 
Here, the distance corresponds to the maximal difference between the two cumulative distribution functions. 

![](./README.figures/CytoCompare-manuscript-Figure1.png)

In our approach, marker similarity distances below a threshold, defined by the user, model a success in a Bernoulli experiment (steps 2).
The p-value, computed by an exact right-tailed binomial test, asserts if the proportion of marker successes is statistically higher than a proportion of expected marker successes defined by the user.
Using this p-value, it is then possible to assert whether the two cell clusters have a phenotype statistically similar or not (step 3).
The marker similarity measures are also aggregated using a weighted mean to obtain one single similarity measure for each comparison of two cell clusters. 

If only one object is provided to the `compare()` function then the comparisons will be performed between all clusters of this object. 
If two objects are provided to the `compare()` function then the comparisons will be performed between all possible pairs of cluster between these two objects.

Importantly, users can define their own function to perform the statistical comparisons of the profiles, using the `method` parameter (please refer to section 9 for more details). 

Different parameters can be defined, via the `method.params` named list, to specify the behavior of such kind of comparisons:

* the `D.th` parameter indicates the similarity threshold to use for comparison between two cell clusters (default values is set to 0.30);
* the `P` parameter indicates the proportion of marker successes to statistically consider the similarity significant (default values is set to 0.75).

For example, comparisons between different cell clusters can be done using the following commands:

```r
# compares a set of cell clusters from the same CLUSTER object, using a specific MWEIGHTS object
res.clusters_b.small <- compare(bm_example.clusters.b[c("205","208","147")],mweights=bm_example.mweights)
print(res.clusters_b.small)
```

```
## Object class: RES
## Number of comparisons: 9
## Markers: 
## 103-Viability
## 110-CD3
## 110_114-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 115-CD45
## 139-CD45RA
## 141-pPLCgamma2
## 142-CD19
## 144-CD11b
## 145-CD4
## 146-CD8
## 147-CD20
## 148-CD34
## 150-pSTAT5
## 151-pERK1/2
## 152-Ki67
## 153-pMAPKAPK2
## 154-pSHP2
## 156-pZAP70/Syk
## 158-CD33
## 159-pSTAT3
## 160-CD123
## 164-pSLP-76
## 165-pNFkB
## 166-IkBalpha
## 167-CD38
## 168-pH3
## 169-pP38
## 170-CD90
## 171-pBtk/Itk
## 172-pS6
## 174-pSrcFK
## 175-pCrkL
## 176-pCREB
## 191-DNA
## 193-DNA
## Cell Length
## Number of markers: 39
## Profiles present in the comparisons:
## CLUSTER:b:205
## CLUSTER:b:208
## CLUSTER:b:147
```


```r
# compares all cell clusters from the same CLUSTER object, using a specific MWEIGHTS object
res.clusters_b <- compare(bm_example.clusters.b,mweights=bm_example.mweights)
```



```r
# compares cell clusters from different CLUSTER objects, using a specific MWEIGHTS object
bm_example.clusters.b@name         <- "b"
bm_example.clusters.mono@name      <- "mono"
bm_example.clusters.tCD4naive@name <- "tCD4naive"
bm_example.clusters.tCD4mem@name   <- "tCD4mem"
bm_example.clusters.tCD8naive@name <- "tCD8naive"
bm_example.clusters.tCD8mem@name   <- "tCD8mem"
clusters  <- list(bm_example.clusters.b,
                  bm_example.clusters.mono,
                  bm_example.clusters.tCD4naive,
                  bm_example.clusters.tCD4mem,
                  bm_example.clusters.tCD8naive,
                  bm_example.clusters.tCD8mem)
res.clusters <- RES()
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    res.clusters <- c(res.clusters,compare(clusters[[i]],clusters[[j]],mweights=bm_example.mweights))
  }
}
```


# <a name="res_object"></a> 7. Manipulation of comparison results

## <a name="res_structure"></a> 7.1. Structure of comparison results
The `RES` object is a S4 object containing one or several comparison results. 
This object mainly stores for each comparison result: the aggregated distance and the marker distances, the associated similarity p-value and the marker similarity successes.

Different slots are available for a given `RES` object:

* the slot `comparisons` is a data.frame containing for each comparison: the profile names,the similarity distance, the distance threshold used and the associated p-value
* the slot `comparisons.nb` is an integer indicating the number of comparisons 
* the slot `markers` is a character vector containing the marker names used in the comparisons
* the slot `marker.weights` is a numeric vector containing the weigths associated to each marker involded in the comparisons
* the slot `marker.distances` is a data.frame containing the marker similarity for each comparison
* the slot `marker.successes` is a data.frame containing the marker successes for each comparison

## <a name="res_summarization"></a> 7.2. Summarization of comparison results
The `print()` or `show()` functions can be used to display a summary of the `RES` objects.

A textual representation of a `RES` object can be displayed using the following commands:

```r
# prints the summary of a given RES object
print(res.clusters_b.small)
```

```
## Object class: RES
## Number of comparisons: 9
## Markers: 
## 103-Viability
## 110-CD3
## 110_114-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 115-CD45
## 139-CD45RA
## 141-pPLCgamma2
## 142-CD19
## 144-CD11b
## 145-CD4
## 146-CD8
## 147-CD20
## 148-CD34
## 150-pSTAT5
## 151-pERK1/2
## 152-Ki67
## 153-pMAPKAPK2
## 154-pSHP2
## 156-pZAP70/Syk
## 158-CD33
## 159-pSTAT3
## 160-CD123
## 164-pSLP-76
## 165-pNFkB
## 166-IkBalpha
## 167-CD38
## 168-pH3
## 169-pP38
## 170-CD90
## 171-pBtk/Itk
## 172-pS6
## 174-pSrcFK
## 175-pCrkL
## 176-pCREB
## 191-DNA
## 193-DNA
## Cell Length
## Number of markers: 39
## Profiles present in the comparisons:
## CLUSTER:b:205
## CLUSTER:b:208
## CLUSTER:b:147
```

## <a name="res_extraction"></a> 7.3. Extraction of comparison results
The extract function `[i]` can be used to extract comparisons of `RES` objects. 
The parameter `i` represents a vector of the comparisons to extract. 
This vector can be a numeric, logical or character.

For example, a subset of comparisons of a `RES` object can be extracted using the following command:

```r
# extracts the two first comparison of a given RES object
res.clusters_b.small_subset <- res.clusters_b.small[1:2]
print(res.clusters_b.small_subset)
```

```
## Object class: RES
## Number of comparisons: 2
## Markers: 
## 103-Viability
## 110-CD3
## 110_114-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 115-CD45
## 139-CD45RA
## 141-pPLCgamma2
## 142-CD19
## 144-CD11b
## 145-CD4
## 146-CD8
## 147-CD20
## 148-CD34
## 150-pSTAT5
## 151-pERK1/2
## 152-Ki67
## 153-pMAPKAPK2
## 154-pSHP2
## 156-pZAP70/Syk
## 158-CD33
## 159-pSTAT3
## 160-CD123
## 164-pSLP-76
## 165-pNFkB
## 166-IkBalpha
## 167-CD38
## 168-pH3
## 169-pP38
## 170-CD90
## 171-pBtk/Itk
## 172-pS6
## 174-pSrcFK
## 175-pCrkL
## 176-pCREB
## 191-DNA
## 193-DNA
## Cell Length
## Number of markers: 39
## Profiles present in the comparisons:
## CLUSTER:b:205
## CLUSTER:b:208
```

## <a name="res_combination"></a> 7.4. Combination of comparison results
The combine function `c()` can be used to combine two or several `RES` objects.
RES objects can be combined to an empty RES object (i.e. `RES()`).

A combination of two or more `RES` objects can be done using the following commands:

```r
# combines a set of RES objects into a single RES object
res.clusters <- c(res.clusters_b.small_subset[1],res.clusters_b.small_subset[2])
print(res.clusters)
```

```
## Object class: RES
## Number of comparisons: 2
## Markers: 
## 103-Viability
## 110-CD3
## 110_114-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 115-CD45
## 139-CD45RA
## 141-pPLCgamma2
## 142-CD19
## 144-CD11b
## 145-CD4
## 146-CD8
## 147-CD20
## 148-CD34
## 150-pSTAT5
## 151-pERK1/2
## 152-Ki67
## 153-pMAPKAPK2
## 154-pSHP2
## 156-pZAP70/Syk
## 158-CD33
## 159-pSTAT3
## 160-CD123
## 164-pSLP-76
## 165-pNFkB
## 166-IkBalpha
## 167-CD38
## 168-pH3
## 169-pP38
## 170-CD90
## 171-pBtk/Itk
## 172-pS6
## 174-pSrcFK
## 175-pCrkL
## 176-pCREB
## 191-DNA
## 193-DNA
## Cell Length
## Number of markers: 39
## Profiles present in the comparisons:
## CLUSTER:b:205
## CLUSTER:b:208
```

It is to note that the comparisons to combine must be unique.

## <a name="res_visualization"></a> 7.5. Visualization of single comparison results
The `plot()` function can be used to plot the `RES` objects, via bar plots.
In such representation each bar corresponds to a marker with a height proportional to the marker similarity measures, and where the bars are colored if they model a success in the Bernoulli experiment. 
Representations are generated using the `ggplot2` library [10] and can be modified by users.

Two `RES` objects, which contains one or several comparison results, can be plotted using the following commands:

```r
# plots the second comparison results based on distance measurements
plot(res.clusters_b.small[2])
```

```
## [1] 0.35
```

<img src="README.figures/res_single_plot1-1.png" style="display: block; margin: auto;" />
*The cell cluster profiles are similar*

If the `RES` object contains several comparison results, then all the different comparison results will be plotted separately using the following command:

```r
# plots three comparison results
plot(res.clusters_b.small[c(2,3,4)])
```

```
## [1] 0.35
## [1] 0.35
## [1] 0.35
```

<img src="README.figures/res_multiple_plot-1.png" style="display: block; margin: auto;" />

## <a name="res_d3js"></a> 7.6. Visualization of comparison results with D3.js
The comparison results (i.e. a `RES` object) can also be visualized using circular graphs or using Multidimensional scaling (MDS) representations [11]. 
Both representations are generated as Scalable Vector Graphics (SVG) elements using the Data-Driven Documents library ([D3.js](http://d3js.org/)) [12]. 
Representation results are interactive and saved as HTML files. 
It is to note that, you can use the `webshot` R package to automatically convert these HTML files to png, jpg or pdf files.

### <a name="res_d3js_circular"></a> 7.6.a. Visualization of multiple comparison results
In the case of a circular graph, each node in the representation is a cell, a cell cluster, or a gate profile. Links between nodes represent similar profiles, based on the aggregated p-values. 
Only comparison with a p-value below a specific threshold are represented.

Through two sliders, users can specify the link tension in the representation and can change the p-value cutoff. 
In the HTML file, users can obtained more specific details about the significant similarities or inclusions by moving the mouse cursor over the links.

Circular graph representations can be generated and saved as a HTML files using the following commands:

```r
# generates a circular graph representation based on comparison results between all cell clusters of the B cell population
res.graph(res.clusters_b,filename="clusters_b.html")
```

![](./README.figures/graph_clusters_b.png)
 


```r
# generates a circular graph representation based on comparison results between all cell clusters
res.graph(res.clusters,filename="graph_clusters.html")
```

![](./README.figures/graph_clusters.png)
 


```r
# generates a circular graph representation based a more complex example
random.cell.b    <- sample(bm_example.cells.b@profiles.nb,10)
random.cell.mono <- sample(bm_example.cells.mono@profiles.nb,10)
profiles <- list(bm_example.clusters.b,
                 bm_example.gates.b, 
                 bm_example.cells.b[random.cell.b],
                 bm_example.cells.mono[random.cell.mono],
                 bm_example.gates.mono,
                 bm_example.clusters.mono, 
                 bm_example.clusters.tCD4naive,
                 bm_example.clusters.tCD4mem,
                 bm_example.clusters.tCD8naive,
                 bm_example.clusters.tCD8mem)
res_graph <- RES()
for(i in 1:(length(profiles))){
  for(j in i:length(profiles)){
    res_graph <- c(res_graph,compare(profiles[[i]],profiles[[j]],mweights=bm_example.mweights))
  }
}
res.graph(res_graph,filename="graph_full.html")
```

![](./README.figures/graph_full.png)
 


### <a name="res_d3js_mds"></a> 7.6.b. MDS visualization
MDS methods aim to represent the similarities and differences among high dimensionality objects into a space of a lower dimensions, generally in two or three dimensions for visualization purposes [11].
In MDS representations, the Kruskal Stress (KS) indicates the amount of information lost during the dimensionality reduction process.

In the case of a MDS representation, each dot in the visualization represents a cell cluster.
Distances between the dots are then proportional to comparison measures calculated between the cluster.

A MDS representation can be generated and saved as a HTML file using the following commands:

```r
# generates a MDS representation based on small portion of B cells and monocytes profiles
res.mds(res.cells,filename="mds_cells.html")
```

![](./README.figures/mds_cells.png)
 



```r
# generates a MDS representation based on a set of cell cluster profiles
res.mds(res.clusters,filename="mds_clusters.html")
```

![](./README.figures/mds_clusters.png)
 



```r
# generates a MDS representation based on a set of cell cluster profiles and use different colors
profiles_names <- sort(unique(c(res.clusters@comparisons$profile1,res.clusters@comparisons$profile2)))
cols <- profiles_names
cols[grep("CLUSTER:b",cols)]         <- "purple"
cols[grep("CLUSTER:mono",cols)]      <- "green"
cols[grep("CLUSTER:tCD4naive",cols)] <- "orange"
cols[grep("CLUSTER:tCD4mem",cols)]   <- "red"
cols[grep("CLUSTER:tCD8naive",cols)] <- "cyan"
cols[grep("CLUSTER:tCD8mem",cols)]   <- "blue"
names(cols) <- profiles_names
res.mds(res.clusters,"mds_clusters_colors.html",cols=cols)
```

![](./README.figures/mds_clusters_colors.png)



# <a name="miscellaneous"></a> 8. Miscellaneous functions

# <a name="dheatmap"></a> 8.2. Density heatmap representations of cell cluster profiles
Density heatmap of a `CLUSTER` object can be generated and plotted using the `dheatmap()` function.
In such representation, each bar corresponds to a marker and the color gradient is proportional to the marker expression density. 

A marker density heatmap representation of a `CLUSTER` object can be displayed using the following command:

```r
# plots a marker density heatmap for a given CLUSTER object
dheatmap(bm_example.clusters.b[1])
```

<img src="README.figures/single_dheatmap_plot-1.png" style="display: block; margin: auto;" />

If the `CLUSTER` object contains several cell cluster profiles, then all the different cell cluster profiles will be plotted separately using the following command:

```r
# plots marker density heatmaps for the first 3 cell cluster profiles of a CLUSTER object
dheatmap(bm_example.clusters[1:3])
```

<img src="README.figures/multiple_dheatmap_plot-1.png" style="display: block; margin: auto;" />


# <a name="compare_template"></a> 9. Template of the compare() function
The `compare()` function can be customized via the `method` parameter, in order to specify a function that can handle the different statistical comparisons. 

Such function template should look like the following:


```r
compare_default <- function(profile1.type,
                            profile2.type,
                            profile1.intensities,
                            profile1.mean,
                            profile1.sd,
                            profile1.density,
                            profile1.nbcells,
                            profile1.range,
                            profile2.intensities,
                            profile2.mean,
                            profile2.sd,
                            profile2.density,
                            profile2.nbcells,
                            profile2.range,
                            ...){ 
  
# code for handling the different profile comparisons
  
  res <- list(measure      = measure,
          pvalue           = pvalue,
          marker.distances = measures,
          marker.successes = successes)
      
  return(res)
}
```

where: 

* `profile1.type` is a character specifying the type of the first profile (`CELL`, `CLUSTER` or `GATE`)
* `profile2.type` is a character specifying the type of the second profile (`CELL`, `CLUSTER` or `GATE`)
* `profile1.intensities` a numeric vector containing the marker intensities for the first cell profile
* `profile1.mean` is a numeric vector containing the marker expression means of the first cluster profile
* `profile1.sd` a numeric vector containing the marker expression standard deviations of the first cluster profile
* `profile1.density` is a named list containing the marker expression densities (`DENSITY` objects) of the first cluster profile
* `profile1.nbcells` a numeric value containing the number of cells associated with the first cluster profile
* `profile1.range` is a numeric array containing the marker intensity ranges of the second gate profile
* `profile2.intensities` is a numeric vector containing the marker intensities for the second cell profile
* `profile2.mean` is a numeric vector containing the marker expression means of the second cluster profile
* `profile2.sd` is a numeric vector containing the marker expression standard deviations of the second cluster profile
* `profile2.density` is a named list containing the marker expression densities (`DENSITY` objects) of the second cluster profile
* `profile2.nbcells` is a numeric value containing the number of cells associated with the second cluster profile
* `profile2.range` is a numeric array containing the marker intensity ranges of the second gate profile
* `...` are other parameters for the custom comparison function

and where the returned named list has the following elements:

* `measure` corresponds to the aggregated distance between the two profiles (or NA in case of inclusion assessments)
* `pvalue` corresponds to the similarity or inclusion p-value between the two profiles
* `marker.distances` corresponds to the marker distances (or vector of NA in case of inclusion assessments)
* `marker.successes` corresponds to the marker similarity or inclusion successes

The parameterization of the `compare()` function with a custom comparison function can be done using the following command:

```r
# compares two cytometry profiles using a customized function and a named parameter list
compare(profile1,profile2,method="compare.sub",method.para=list(para1=1,para2=2))
```

# <a name="object_structure"> 10. Structure of the cytometry objects


##<a name="cluster_structure"></a> 10.3. Structure of CLUSTER object
The `CLUSTER` object is a S4 object containing one or several cell cluster profiles. 
This object mainly stores for each cell cluster profile: the means, the standard deviations and the densities of each marker.

Different slots are available for a given `CLUSTER` object:

* the slot `name` is a character indicating the internal name of the `CLUSTER` object
* the slot `profiles` is a character vector containing the names of the cell cluster profiles
* the slot `profiles.nb` is an integer value indicating the number of cell cluster profiles
* the slot `profiles.sizes` is an integer vector indicating the number of cells associated to each cluster profile
* the slot `markers` is a character vector containing the marker names
* the slot `markers.nb` is an integer value indicating the number of markers
* the slot `markers.clustering` is a logical vector specifying the markers used as clustering markers
* the slot `means` is a numeric matrix containing the means of each marker for each cluster profile
* the slot `sd` is a numeric matrix containing the standard deviations of each marker for each cluster profile
* the slot `densities` is a matrix of `DENSITY` objects containing the density of each marker for each cluster profile
* the slot `overview.function` is a character specifying the name of a function to call when plotting the `CLUSTER` object overview (please refer to the documentation of the `plot()` function)
* the slot `graph` is an object that can be used to store a visual representation of the cell clusters (e.g. a SPADE tree)
* the slot `graph.layout` is a numeric matrix that can be used to store the positions of cell clusters in a 2-dimensional space (e.g. a SPADE tree layout)


## <a name="mweights_object"></a> 10.5. Structure of MWEIGHTS object

### <a name="mweights_structure"></a> 10.5.a. Structure of MWEIGHTS object
The `MWEIGHTS` object is a S4 object containing the marker weights to use in the comparison computations. This object mainly stores for each marker: the markers names and marker weights.  

Different slots are available for a given `MWEIGHTS` object:

* the slot `markers` is a character vector containing the marker names
* the slot `weights` is a numeric vector containing the marker weights


### <a name="mweights_summarization"></a> 10.5.b. Summarization of MWEIGHTS object
The `print()` or `show()` functions can be used to display a summary of the `MWEIGHTS` objects.

A textual representation of a `MWEIGHTS` object can be displayed using the following command:

```r
# prints the summary of a given MWEIGHTS object
print(bm_example.mweights)
```

```
## Object class: MWEIGHTS
## Number of markers: 39
## Markers with associated weights:
##  103-Viability: 0
##  110-CD3: 0
##  110_114-CD3: 1
##  111-CD3: 0
##  112-CD3: 0
##  114-CD3: 0
##  115-CD45: 1
##  139-CD45RA: 1
##  141-pPLCgamma2: 0
##  142-CD19: 1
##  144-CD11b: 1
##  145-CD4: 1
##  146-CD8: 1
##  147-CD20: 1
##  148-CD34: 1
##  150-pSTAT5: 0
##  151-pERK1/2: 0
##  152-Ki67: 0
##  153-pMAPKAPK2: 0
##  154-pSHP2: 0
##  156-pZAP70/Syk: 0
##  158-CD33: 1
##  159-pSTAT3: 0
##  160-CD123: 1
##  164-pSLP-76: 0
##  165-pNFkB: 0
##  166-IkBalpha: 0
##  167-CD38: 1
##  168-pH3: 0
##  169-pP38: 0
##  170-CD90: 1
##  171-pBtk/Itk: 0
##  172-pS6: 0
##  174-pSrcFK: 0
##  175-pCrkL: 0
##  176-pCREB: 0
##  191-DNA: 0
##  193-DNA: 0
##  Cell Length: 0
```


### <a name="mweights_extraction"></a> 10.5.c. Extraction of MWEIGHTS object
The extract function `[i]` can be used to extract marker weights of `MWEIGHTS` objects. The parameter `i` represents a vector of the markers to extract. This vector can be a numeric, logical or character.

For example, a subset of markers of a `MWEIGHTS` object can be extracted using the following command:

```r
# extracts 5 markers of a given `MWEIGHTS` object
bm_example.mweights.sub1 <- bm_example.mweights[5:10]
print(bm_example.mweights.sub1)
```

```
## Object class: MWEIGHTS
## Number of markers: 6
## Markers with associated weights:
##  112-CD3: 0
##  114-CD3: 0
##  115-CD45: 1
##  139-CD45RA: 1
##  141-pPLCgamma2: 0
##  142-CD19: 1
```

### <a name="mweights_set"></a> 10.5.d. Combination of MWEIGHTS object
The set function `[i] <- value` can be used to set marker weights of `MWEIGHTS` objects.

For example, the markers weights of a `MWEIGHTS` object can be set using the following command:

```r
# sets some weights of a given `MWEIGHTS` object
bm_example.mweights.sub1[1:3] <- c(0.3,0,0.7)
print(bm_example.mweights.sub1)
```

```
## Object class: MWEIGHTS
## Number of markers: 6
## Markers with associated weights:
##  112-CD3: 0.3
##  114-CD3: 0
##  115-CD45: 0.7
##  139-CD45RA: 1
##  141-pPLCgamma2: 0
##  142-CD19: 1
```


### <a name="mweights_visualization"></a> 10.5.e. Visualization of MWEIGHTS object
The `plot()` function can be used to plot the `MWEIGHTS` objects, via bar plots where each bar corresponds to a marker and where the bar heights are proportional to the marker weights. Representations are generated using the `ggplot2` library [10] and can be modified by users.

A `MWEIGHTS` object can be plotted using the following command:

```r
# plots a given MWEIGHTS object
plot(bm_example.mweights)
```

<img src="README.figures/mweights_single_plot-1.png" style="display: block; margin: auto;" />

## <a name="density_object"></a> 10.6 Structure of DENSITY object

### <a name="density_structure"></a> 10.6.a. Structure of DENSITY object
The `DENSITY` object is a S4 object containing the marker expression densities. This object mainly stores for each marker: the bin characteristics, the negative and positive marker densities values, and the number of cells used in the density estimation.

The `DENSITY` objects are present in `CLUSTER` objects, and can be accessed via the `densities` slot (please note the `[[1]]`):

```r
# accesses to the DENSITY object of the first marker of the first cell cluster profile
density <- bm_example.clusters.b[1]@densities[1,1][[1]]
```

Different slots are available for a given `DENSITY` object:

* the slot `name` is a character indicating the internal name of the `DENSITY` object
* the slot `bin.interval` is a numeric vector of two values specifying the density boundaries
* the slot `bin.nb` is a numeric vector of two values specifying the numbers of negative and positive bins
* the slot `values.pos` is a numeric vector containing the positive density bins
* the slot `values.neg` is a numeric vector containing the negative density bins
* the slot `point.nb` is a numeric value indicating the number of point used to compute the expression density
* the slot `bin.width` is a numeric value indicating the width of the bins used in the density estimation

### <a name="density_summarization"></a> 10.6.b. Summarization of DENSITY object
The `print()` or `show()` functions can be used to display a summary of the `DENSITY` objects.

A textual representation of a `DENSITY` object can be displayed using the following command:

```r
# prints the summary of a given DENSITY object
density <- bm_example.clusters.b[1]@densities[1,1][[1]]
print(density)
```

```
## Object class: DENSITY
## Bin widths: 0.05
## Number of bins: 133
## Minimal bin: -1.5
## Maximal bin: 5.14999999999999
## Number of negative bins: 30
## Number of positive bins: 103
## Number of values used in the density estimation: 1859
```

### <a name="density_visualization"></a> 10.6.c. Visualization of DENSITY object
The `plot()` function can be used to plot the `DENSITY` objects, via histogram plots where each bar corresponds to a density bin and where a smooth line represents the average estimation of the marker expression density. Representations are generated using the `ggplot2` library [10] and can be modified by users.

A `DENSITY` object can be plotted using the following commands:

```r
# plots the expression density for the first marker of the first cluster profile
density <- bm_example.clusters.b[1]@densities[1,1][[1]]
plot(density)
```

<img src="README.figures/multiple_density_plot-1.png" style="display: block; margin: auto;" />

```r
# plots the expression density for the second marker of the first cluster profile
density <- bm_example.clusters.b@densities[1,2][[1]]
plot(density)
```

<img src="README.figures/multiple_density_plot-2.png" style="display: block; margin: auto;" />

```r
# plots the expression density for the first marker of the second cluster profile
density <- bm_example.clusters.b@densities[2,1][[1]]
plot(density)
```

<img src="README.figures/multiple_density_plot-3.png" style="display: block; margin: auto;" />

# 13. <a name="license"></a> License
CytoCompare is freely distributed under the GLP-3 license.


# 14. <a name="references"></a> References 
[1] - Bendall, S. C., Simonds, E. F., Qiu, P., Amir, E. D., Krutzik, P. O., Finck, R., Nolan, G. P. (2011). Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. Science (New York, N.Y.), 332(6030), 687-96.

[2] - Gregori, G., Rajwa, B., Patsekin, V., Jones, J., Furuki, M., Yamamoto, M., & Paul Robinson, J. (2014). Hyperspectral cytometry. Current Topics in Microbiology and Immunology, 377, 191-210.

[3] - FlowJo | Single Cell Analysis Software - 
http://www.flowjo.com/

[4] - Cytobank - 
https://www.cytobank.org

[5] - Qiu, P., Simonds, E. F., Bendall, S. C., Gibbs, K. D., Bruggner, R. V, Linderman, M. D., Plevritis, S. K. (2011). Extracting a cellular hierarchy from high-dimensional cytometry data with SPADE. Nature Biotechnology, 29(10), 886-91.

[6] - Bruggner, R. V, Bodenmiller, B., Dill, D.L., Tibshirani, R.J., and Nolan, G.P. (2014). Automated identification of stratifying signatures in cellular subpopulations. Proceedings of the National Academy of Sciences of the United States of America 111, E2770-E2777.

[7] - Amir ED, Davis KL, Tadmor MD, Simonds EF, Levine JH, Bendall SC, Shenfeld DK, Krishnaswamy S, Nolan GP, Pe'er D: viSNE enables visualization of high dimensional single-cell data and reveals phenotypic heterogeneity of leukemia. Nat Biotechnol 2013, 31:545-52.

[8] - Flow Cytometry Data File Format Standards. http://isac-net.org/Resources-for-Cytometrists/Data-Standards/Data-File-Standards/Flow-Cytometry-Data-File-Format-Standards.aspx

[9] - Ellis B, Haaland P, Hahne F, Meur NL, Gopalakrishnan N, Spidlen J and Jiang M. flowCore: flowCore: Basic structures for flow cytometry data. R package version 1.34.7.

[10] - Grammar of Graphics library -
http://ggplot2.org/

[11] - Kruskal, J. B., and Wish, M. (1978), Multidimensional Scaling, Sage University Paper series on Quantitative Application in the Social Sciences, 07-011. Beverly Hills and London: Sage Publications.

[12] - A JavaScript visualization library for HTML and SVG -
http://d3js.org
