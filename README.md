# CytoCompare: an R Package for Computational Comparisons of Cytometry Profiles



# Table of Contents
1. [Package overview](#package_overview)
2. [Package installation](#package_installation)
3. [Importation of cytometry profiles](#object_importation)
    1. [Importation of cell profiles from FCS files](#fcs_importation)
    2. [Importation of cell profiles from viSNE FCS files](#visne_importation)
    3. [Importation of cell profiles from text files](#cell_importation)
    4. [Importation of cell cluster profiles from SPADE results](#spade_importation)
    5. [Importation of cell cluster profiles from CITRUS results](#citrus_importation)
    6. [Importation of cell cluster profiles from text files](#cluster_importation)
    7. [Importation of gate profiles from Gating-ML XML files](#gatingml_importation)
    8. [Importation of gate profiles from text files](#gate_importation)
4. [Manipulation of cytometry profiles](#object_manipulation)
    1. [Extraction of cytometry profiles](#object_extraction)
    2. [Combination of cytometry profiles](#object_combination)
    3. [Transformation of cytometry profiles](#object_transformation)
    4. [Exportation of cytometry profiles](#object_export)
5. [Representations of cytometry profiles](#object_representation)
    1. [Visualization of cell profiles](#cell_visualization)
    2. [Visualization of cluster profiles](#cluster_visualization)
    3. [Visualization of gate profile](#gate_visualization)
    4. [Pairwise visualization of cytometry profiles](#pairwise_visualization)
6. [Comparisons of cytometry profiles](#object_comparison)
    1. [Overview of the comparison approach](#compare_function)
    2. [Comparisons between two cell profiles](#cell_cell_compare)
    3. [Comparisons between two cell cluster profiles](#cluster_cluster_compare)
    4. [Comparisons between two gate profiles](#gate_gate_compare)
    5. [Comparisons between a cell profile and a gate profile](#cell_gate_compare)
    6. [Comparisons between a cell profile and a cell cluster profile](#cell_cluster_compare)
    7. [Comparisons between a cell cluster profile and a gate profile](#cluster_gate_compare)
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

Flow and mass cytometry are experimental techniques used to characterize cells at a single-cell level. Both techniques use labeled antibodies to measure cell marker expressions. Flow cytometry uses antibodies conjugated with fluorochromes to quantify stained cells with a laser detection system. The recently introduced mass cytometry ([CyTOF](https://www.fluidigm.com/products/cytof) [[1](http://www.ncbi.nlm.nih.gov/pubmed/21551058)]) uses antibodies conjugated with metals to quantify stained cells with a mass spectrometer. While flow cytometry can currently handle up to 18 fluorochromes, mass cytometry can currently handle up to 40 metals. Mass cytometry offers important perspectives as it can potentially evaluate more than 100 metals. Another recently introduced technique, called hyper-spectral cytometry [[2](http://www.ncbi.nlm.nih.gov/pubmed/24271566)], combines ultrafast optical spectroscopy with flow cytometry and can also increase up to 40 the number of simultaneously usable fluorochromes. Such significant technological improvements make flow and mass cytometry suitable for massive data-mining and bioinformatics developments, in order to better explore and analyze complex cell interaction systems.

Cytometry data can be manually analyzed using a gating strategy (hierarchical gating usually followed by Boolean gating) or using automatic gating/clustering algorithms.
In both approaches, the aim is to identify cell populations having similar profiles, based on the expressions of selected markers.
Hierarchical gating is done using biplot representations, where each axis represents the intensity of a marker of interest and where dots represent the cells. Through iterative steps, operators select cell subsets and explore their marker expression profiles. Boolean gating is done by quantifying cells present within all possible combinations of gates for a population of interest. Even if Boolean gating can be done in an automatic way, gates still have to be delineated for each cell marker. FlowJo [[3](http://www.flowjo.com)] and CytoBank [[4](https://www.cytobank.org)] are among the most common software for manual gating. Both hierarchical and Boolean gating can be very fastidious and biased. 
The difficulty mainly comes from the high dimensionality of data, that is to say from the high number of markers and cell populations. This task becomes even more complicated as the number of markers increases or as the studies become more complex. 
On the other hand, automatic gating methods (such as [SPADE](http://cytospade.org/) [[5](http://www.ncbi.nlm.nih.gov/pubmed/21964415)] or [Citrus](https://github.com/nolanlab/citrus/wiki) [[6](http://www.ncbi.nlm.nih.gov/pubmed/24979804
)]) use various algorithms to computationally identify cell populations having similar profiles and provides then less biased results. Moreover, several dimensionality reduction  tools, such as viSNE [[7](http://www.ncbi.nlm.nih.gov/pubmed/23685480)], have been proposed to better visualize cells in cytometry profiles. Once identified cell populations, also named cell clusters, need to be deeper characterized or associated with known phenotypes for further investigations. 

Characterization of identified cell populations is a challenge becoming increasingly dominant in analysis of high-dimensional cytometry data.

A typical task in this characterization process is to identify similar cell populations within a study or among different studies. Indeed, automatic gating methods can generate a large number of cell clusters, with overlapping marker expression phenotypes, which need to be organized further. Such organization consists then to gather identified cell clusters by families of similar phenotypes. Comparisons of cell cluster profiles among different studies can also reveal the effects of the clustering parameters. For instance, researchers can be interested to compare different automatic gating analyses and to highlight cell clusters found in the different analyses. In a complementary manner, researchers can be interested to compare automatic analyses performed in different experimental conditions.

Another typical task in this characterization process is to identify included cytometry profiles, such as cell or cell cluster profiles included in gate profiles. This task is trivial when identifying cells included within cytometry gates and corresponds to regular hierarchical gating procedure. The introduction of automatic gating approaches transposed this task to the identification of cells included in cell clusters (instead gates). Moreover, the identification of cell cluster profiles included in gate profiles can also be relevant for the characterization of the cytometry profiles. This last comparison can be interested to assess the effectiveness of automatic gating with manual gating.

There is currently a lack of computational methods allowing to compare cytometry profiles, in order to identify *similar* profiles or to identify *included* profiles. These absences are critical as these operations have to be done manually and are intrinsically subject to bias.

To answer this need we designed CytoCompare, an R package allowing to compare cytometry profiles obtained from various sources, based on the characteristics of cell markers specified by the user. CytoCompare is then able to identify similar or included cytometry profiles allowing to better interpret them.

Three types of cytometry profiles can be handled by CytoCompare: 

  1. the `CELL` object can handles cell profiles, modeled by intensities of expression markers;
  2. the `CLUSTER` object can handles cell cluster profiles, modeled by densities of expression markers;
  3. the `GATE` object can handles gate profiles, modeled by intensity ranges of expression markers.

Each of these cytometry profiles (handled by `CELL`, `CLUSTER`, or `GATE` objects) can contain one or several cytometry profiles. 

For instance, a set of cell profiles can be formalized by the following matrix: 
<img src="README.figures/table_cell-1.png" style="display: block; margin: auto;" />
where a, b and c are different cells and a_n, b_n and c_n corresponds to their respective expression intensities for marker_n.

For instance, a set of cell cluster profiles can be formalized by the following matrix: 
<img src="README.figures/table_cluster-1.png" style="display: block; margin: auto;" />
where A, B and C are different cell clusters and A_n, B_n and C_n corresponds to their respective expression densities for marker_n.

For instance, a set of gate profiles can be formalized by the following matrix: 
<img src="README.figures/table_gate-1.png" style="display: block; margin: auto;" />
where alpha, beta and gamma are different gates and [min_{alpha n},max_{alpha n}], [min_{beta n},max_{beta n}] and [min_{gamma n},max_{gamma n}] corresponds to their respective expression boundaries for marker_n.

The following diagram summarizes the comparisons available in CytoCompare, depending on the types of the cytometry profiles:

![](README.figures/ComparisonTypes.png)
 

Depending on the type of cytometry profile to compare, different strategies are used. The following table summarizes those strategies.

Cytometry profile 1 | Cytometry profile 2 | Type of comparisons
:-------------:|:-------------:|------------------------------------------
    `CELL`     |   `CELL`      |   **similarity** of cell marker expressions based on the Euclidean distance (D_i for each marker_i)
    `CLUSTER`  |   `CLUSTER`   |   **similarity** of cluster marker expression densities based on the Kolmogorov-Smirnov distance (D_i for each marker_i)
    `GATE`     |   `GATE`      |   **similarity** of gate marker expression boundaries based on the Kolmogorov-Smirnov distance (D_i for each marker_i)
    `CELL`     |   `GATE`      |   **inclusion** of the cell marker expressions within the gate expression boundaries 
    `CELL`     |   `CLUSTER`   |   **inclusion** of the cell marker expressions within the cluster expression ranges\*
    `CLUSTER`  |   `GATE`      |   **inclusion** of the cluster expression ranges* within the gate expression boundaries

\* Cluster expression ranges are calculated based quantiles of marker expression densities

For each comparison of two cytometry profiles, CytoCompare computes a p-value asserting the significance of the profile similarity or inclusion. Moreover, a similarity measure is provided when comparing cytometry profile of similar types.

The following chart summarizes the approach used in CytoCompare.

![](README.figures/ComparisonMethods.png)
 


In the context of a similarity comparison (step 1.a. and 1.b.), a distance is computed between each marker of the two profiles (D_i). Marker distances below a distance threshold, specified by the user, will correspond to a marker similarity success. Euclidean distance will be used when comparing cell profiles while the Kolmogorov-Smirnov distance will be used when comparing cluster or gate profiles. A weight can be associated to each marker, in order to modulate their importance. An aggregation (step 3.) of marker distances is performed using an exact binomial test where marker successes are considered as successful Bernoulli experiments. Thereby, the proportion of marker successes is compared to a probability of success (P) specified by the user. A aggregated distance (D), corresponding to the weigthed mean of marker distances, is additionaly returned.

In the context of an inclusion comparison (step 2.a. and 2.b.), an inclusion assessment is performed for each marker of the two profiles. As illustrated below, a cell profile marker is considered as included in a gate profile when its expression value is within the range of the marker boundaries. Similarly, a cell profile marker is considered as included in a cell cluster profile when its expression value is within the range of the marker cluster defined based on quantiles of marker expression densities. Finally, a cell cluster profile marker is considered as included in a gate profile when its expression boundaries is within the range of the marker gate. As for similarity comparisons, weights associated to each marker. The aggregation (step 3.) of marker inclusion is also performed using an exact binomial test.

![](README.figures/inclusion_assessments.png)
 

We designed CytoCompare in a way that it can be easily used by non bioinformatician experts, but can also be easily customizable by users with more expertise in bioinformatics. CytoCompare offers many visualization representations to make comparison results and intermediary objects easily understandable (such as parallel coordinates, Multidimensional scaling or circular graphs representations). Through the multiple profiles manipulation methods, CytoCompare is also a powerful analysis pipeline for high-dimensional cytometry data. Finally, users can also create the cytometry profiles based on their own file formats and define their own statistical methods for the comparisons of the different types of profiles. 

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

Alternatively, the `install.requiredpackages()` function of the CytoCompare package can be used to install these packages.

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



An example dataset, obtained from public mass cytometry data [[5](http://www.ncbi.nlm.nih.gov/pubmed/21964415)], is available in CytoCompare. This example dataset consists of three cytometry profiles from healthy human bone marrow samples, unstimulated or stimulated by BCR-inductor or IL-7, and measured using a mass cytometry panel of more than 30 cell markers. This panel has been designed to identify a large spectrum of immune cell types. A SPADE analysis has been performed to identify cell clusters profiles. Then, these cell clusters have been manually labeled based on their marker expressions [[5](http://www.ncbi.nlm.nih.gov/pubmed/21964415)]. SPADE cell clusters corresponding to six majors cell types (B, CD33+ monocytes, naive or memory CD4+ T cells, and naive or memory CD8+ T cells) have been extracted. Furthermore, a set of gates have been constructed based on these six major cell types.

The raw data and the R object (`CytoCompareExample.rdata`) corresponding to this example dataset are available on a public ftp server: [ftp://ftp.cytocompare.org/public/](ftp://cytocompare:cytocompare@ftp.cytocompare.org/public/) (username: cytocompare, password: cytocompare).

The R object containing the example dataset can be directly downloaded and loaded in the current R session using the following command:



```r
 # downloads the 'CytoCompareExample.rdata' file on the current folder, if no such file exists, and load it into the current R session
load.examples()
```

Once downloaded, the following objects will be available:

* `bm_example.cells.b`, a `CELL` object containing the cell profiles of the B cell populations
* `bm_example.cells.mono`, a `CELL` object containing the cell profiles of the monocyte cell populations
* `bm_example.cells.tCD4naive`, a `CELL` object containing the cell profiles of the naive CD4+ T cell populations
* `bm_example.cells.tCD8naive`, a `CELL` object containing the cell profiles of the naive CD8+ T cell populations
* `bm_example.cells.tCD4mem`, a `CELL` object containing the cell profiles of the memory CD4+ T cell populations
* `bm_example.cells.tCD8mem`, a `CELL` object containing the cell profiles of the memory CD8+ T cell populations
* `bm_example.clusters`, a `CLUSTER` object containing the cell cluster profiles for all the different cell populations identified by SPADE
* `bm_example.clusters.b`, a `CLUSTER` object containing the cell cluster profiles of the B cell populations identified by SPADE
* `bm_example.clusters.mono`, a `CLUSTER` object containing the cell cluster profiles of the monocyte cell populations  identified by SPADE
* `bm_example.clusters.tCD4naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD4+ T cell populations identified by SPADE
* `bm_example.clusters.tCD8naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD8+ T cell populations identified by SPADE
* `bm_example.clusters.tCD4mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD4+ T cell populations identified by SPADE
* `bm_example.clusters.tCD8mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD8+ T cell populations identified by SPADE
* `bm_example.gates`, a `GATE` object containing the gate profiles constructed based on the 6 main cell populations identified by SPADE
* `bm_example.gates.b`, a `GATE` object containing the gate profiles constructed based on the B cell populations
* `bm_example.gates.mono`, a `GATE` object containing the gate profiles constructed based on the monocyte cell populations
* `bm_example.gates.tCD4naive`, a `GATE` object containing the gate profiles constructed based on the naive CD4+ T cell populations
* `bm_example.gates.tCD8naive`, a `GATE` object containing the gate profiles constructed based on the naive CD8+ T cell populations
* `bm_example.gates.tCD4mem`, a `GATE` object containing the gate profiles constructed based on the memory CD4+ T cell populations
* `bm_example.gates.tCD8mem`, a `GATE` object containing the gate profiles constructed based on the memory CD8+ T cell populations
* `bm_example.mweights`, a `MWEIGHTS` object containing markers that can be used for the comparison computations
* `bm_example.visne`, a list of three `CELL` objects containing the viSNE cell profiles for each biological sample


# <a name="object_importation"></a> 3. Importation of the cytometry profiles 

## <a name="fcs_importation"></a> 3.1. Importation of cell profiles from FCS files
The `import.FCS()` function imports cell profiles from a Flow Cytometry Standard ([FCS](http://isac-net.org/Resources-for-Cytometrists/Data-Standards/Data-File-Standards/Flow-Cytometry-Data-File-Format-Standards.aspx) [8]) file into a `CELL` object (see section 10.2. for more details). The FCS is a standard file format used to store flow and mass cytometry data. These files mainly store the intensities of each marker for each cell profile.
CytoCompare can handles for now the same FCS file versions as supported in the [flowCore](http://www.bioconductor.org/packages/release/bioc/html/flowCore.html) [9] Bioconductor package (FCS versions 1.0, 2.0, 3.0 and 3.1).

An import of a FCS file can be done using the following command:


```r
# imports a FCS containing several cell profiles 
imported.fcs <- import.FCS('./CELL.FCS/Bcells/Marrow1_01_Basal1_Singlets_B cells.fcs')
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

If a set of files is specified, then the files are merged in the import procedure.

An import of a FCS file using a `"log"` transformation in base 10 with a shift of 1 and using the `exclude` parameter can be made using the following command:


```r
# imports a tab separated file containing several cell profiles with some import specifications
imported.fcs <- import.FCS('./CELL.FCS/Bcells/Marrow1_01_Basal1_Singlets_B cells.fcs', trans="log", trans.para=list(log.base=10,log.shift=1), exclude=c("Cell_length"))
```

A summary of a `CELL` object, imported by the `import.FCS()` function, can be done using the following command:

```r
print(imported.fcs)
```

```
## Object class: CELL
## Object name: Marrow1_01_Basal1_Singlets_B cells.fcs
## Number of cell profiles: 14198
## Number of markers: 43
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
## 145-CD4
## 146-CD8
## 148-CD34
## 150-pSTAT5
## 147-CD20
## 152-Ki67
## 154-pSHP2
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
## cluster
```


## <a name="visne_importation"></a> 3.2. Importation of cell profiles from viSNE FCS files
The `import.VISNE()` function imports cell profiles from viSNE FCS files into a `CELL` object (see section 10.2. for more details). viSNE is a visualization algorithm for high-dimensional cytometry data [[7](http://www.ncbi.nlm.nih.gov/pubmed/23685480)]. ViSNE aims to represent cell profiles in a 2-dimensional space using a dimensionality reduction process. For each cell profile a 2-dimensional coordinate is calculated. The 2 dimensions are named tSNE1 and tSNE2 and the results are presented via FCS files having two extra markers corresponding to these dimensions.  

An import of a viSNE FCS file can be done using the following command:


```r
# imports a viSNE FCS file containing several cell profiles 
imported.visne <- import.VISNE('./CELL.viSNE/Marrow1_01_Basal1_Singlets_viSNE.fcs')
```

Some parameters can be specified to apply numeric transformations on the marker expression values, or to exclude some markers in the import or transformation procedures:

* the `exclude` parameter is a character vector containing the marker names to be excluded in the import procedure
* the `trans` parameter is a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation 
* the `trans.para` parameter is a named list containing parameters for the transformation
* the `trans.exclude` parameter is a character vector containing the marker names for which no transformation must be applied on
* the `tSNE1` parameter is a character indicating the marker name of the first viSNE dimension (tSNE1)
* the `tSNE2` parameter is a character indicating the marker name of the second viSNE dimension (tSNE2)


The available transformation functions are (`trans` parameter):

* `"arcsinh"` for an arc sine hyperbolic transformation of the data (default choice)
* `"log"` for a logarithmic transformation of the data
* `"none"` for no transformation of the data

The available transformation parameters are (`trans.para` parameter):

* for the `"arcsinh"` transformation, the scale (cofactor) can be specified using the `arcsinh.scale` parameter
* for the `"log"` transformation, the base and the shift can be specified using the `log.base` and `log.shift` parameters. The 'log.shift' parameter allows the value "auto" which automatically identify the log shift avoiding to apply log transformations on negative values.


A summary of the `CELL` object, imported by the `import.VISNE()` function, can be done using the following command:

```r
print(imported.visne)
```

```
## Object class: CELL
## Object name: Marrow1_01_Basal1_Singlets_viSNE.fcs
## Number of cell profiles: 10000
## Number of markers: 44
## Markers: 
## Time
## Cell Length
## 191-DNA
## 193-DNA
## 103-Viability
## 115-CD45 (v)
## 110-CD3
## 111-CD3
## 112-CD3
## 114-CD3
## 139-CD45RA (v)
## 141-pPLCgamma2
## 142-CD19 (v)
## 144-CD11b (v)
## 145-CD4 (v)
## 146-CD8 (v)
## 148-CD34
## 150-pSTAT5
## 147-CD20 (v)
## 152-Ki67
## 154-pSHP2
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
## tSNE1
## tSNE2
```


## <a name="cell_importation"></a> 3.3. Importation of cell profiles from text files
The `import.CELL()` function imports one or several cell profiles from a tab separated file into a `CELL` object (see section 10.2. for more details). Such tab separated file must contain for each cell profile the intensities of the marker and must be formatted as the following:

* each row must represent a cell profile
* each column must represent a marker
* each cell in the table must contain the marker expression intensities for a given cell profile

The first column must contain the cell names and the first row must contain the marker names.

A typical tab separated file to import must look like the following:

|         | marker_1 | marker_2 | marker_3 | marker_4 | marker_5 | marker_i | marker_n |
|---------|--------|--------|--------|--------|--------|--------|--------|
| cell_1 |  5.23  |  0.42  |  4.26  |  3.23  |  1.57  |  3.14  |  6.01  |
| cell_2 |  2.56  |  2.34  |  5.45  |  2.34  |  1.27  |  1.70  |  6.32  |
| cell_3 |  2.67  |  4.10  |  4.56  |  3.56  |  0.89  |  0.42  |  4.32  |
| cell_4 |  4.33  |  4.32  |  3.89  |  3.45  |  3.45  |  5.07  |  5.65  |
| cell_5 |  5.24  |  2.76  |  5.08  |  2.56  |  5.23  |  4.98  |  7.98  |


An import of such file can be done using the following command:


```r
# imports a tab separated file containing several cell profiles 
imported.cells <- import.CELL('./CELL.TXT/cells.txt')
```

The `exclude` parameter allows to specify some markers to exclude in the import procedure. This parameter takes a character vector containing the marker names to be excluded.

A summary of the `CELL` object, imported by the `import.CELL()` function, can be done using the following command:

```r
print(imported.cells)
```

```
## Object class: CELL
## Object name: cells.txt
## Number of cell profiles: 5
## Number of markers: 7
## Markers: 
## marker_1
## marker_2
## marker_3
## marker_4
## marker_5
## marker_i
## marker_n
```


## <a name="spade_importation"></a> 3.4. Importation of cell cluster profiles from SPADE results
The `import.SPADE()` function imports cell cluster profiles identified by the SPADE algorithm into a `CLUSTER` object (see section 10.3. for more details). The Spanning Tree Progression of Density Normalized Events ([SPADE](http://cytospade.org/) [5]) algorithm is a visualization and analysis algorithm for high-dimensional cytometry data, available in R or in [Cytobank](https://www.cytobank.org/) [4]. SPADE was designed to analyze mass cytometry data but can also handle flow cytometry data.

SPADE identifies clusters of cells having similar expression profiles for selected markers using an agglomerative hierarchical clustering-based algorithm combined with a density-based down-sampling procedure. Given a set of FCS files (usually one file per sample), SPADE identifies cell clusters based on the whole dataset and provides then for each sample the amount of cells present within each cluster.
The number of desired cell clusters to obtain needs to be specified by the user as well as some down-sampling parameters. 
SPADE represents its clustering results using a tree representation where each node represents a cell cluster and where similar cell clusters are linked using a minimal spanning tree approach.
In SPADE tree representations, the sizes of the nodes are proportional to the amount of cells present within each cluster. In order to characterize cell populations and to associate them with known phenotypes, nodes can be gradient-colored based on theirs mean expression intensities for a specific marker.

An import of SPADE clustering results can be done using the following command:


```r
# imports a SPADE folder containing several cell cluster profiles 
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
print(imported.spade)
```

```
## Object class: CLUSTER
## Object name: CLUSTER.SPADE
## Number of clusters profiles: 253
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
## 145-CD4
## 146-CD8
## 148-CD34
## 150-pSTAT5
## 147-CD20
## 152-Ki67
## 154-pSHP2
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
## 115-CD45
## 139-CD45RA
## 142-CD19
## 144-CD11b
## 145-CD4
## 146-CD8
## 148-CD34
## 147-CD20
## 158-CD33
## 160-CD123
## 167-CD38
## 170-CD90
## 110_114-CD3
## Number of clustering markers: 13
## Density bin width: 0.05
## Cluster profile names and number of associated cells:
##  1: 3140 cells
##  2: 951 cells
##  3: 2056 cells
##  4: 1375 cells
##  5: 30 cells
## and 248 more...
```


## <a name="citrus_importation"></a> 3.5. Importation of cell cluster profiles from CITRUS results
The `import.CITRUS()` function imports cell cluster profiles identified by the Citrus algorithm into a `CLUSTER` object (see section 10.3. for more details). [Citrus](https://github.com/nolanlab/citrus/wiki) [[6](http://www.ncbi.nlm.nih.gov/pubmed/24979804
)] is an algorithm developed similarly to SPADE, which can furthermore identify cell clusters associated with different biological condition phenotypes.

An import of Citrus clustering results can be done using the following command:


```r
# imports a Cirtus result file containing several cell cluster profiles 
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


## <a name="cluster_importation"></a> 3.6. Importation of cell cluster profiles from text files
The `import.CLUSTER()` function imports one or several cell cluster profiles from a tab separated file into a `CLUSTER` object (see section 10.3. for more details). In this case, the marker expressions of each cluster are assumed to be normally distributed. Therefore, the tab separated file must contain for each cell cluster profile the means and the standard deviations of the expression markers and must be formatted as the following:

* each row must represent a cell cluster profile
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
# imports a tab separated file containing several cell cluster profiles 
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

It is to note that `CLUSTER` objects constructed via the `import.CLUSTER()` function do not contain the densities of expression markers (please refer to the documentation of the `compare()` function).


## <a name="gatingml_importation"></a> 3.7. Importation of gate profiles from Gating-ML XML files
The `import.GATINGML()` function imports gate profiles from a GatingML-XML file into a `GATE` object (see section 10.4. for more details). GatingML is a standard file format mainly used in cytometry to store definitions of several types of gate (rectangle, ellipses, polygons). CytoCompare allows for now to import 2-dimensional range gates, modeled in GatingML by the entry type 'RectangleGate'.

An import of a GatingML file can be done using the following command:


```r
# imports a GatingML-XML file containing several gate profiles 
imported.gatingml <- import.GATINGML('./GATE.GatingML/cytobank_gates.xml')
```

A parameter can be specified to select specific gates in the import procedure:

* the `filterId` parameter is a character vector indicating the identifiers of the gates to import

A summary of the `GATE` object, imported by the `import.GATINGML()` function, can be done using the following command:

```r
print(imported.gatingml)
```

```
## Object class: GATE
## Object name: cytobank_gates.xml
## Number of gate profiles: 2
## Markers: 
## (Tb159)Di
## (Sm149)Di
## (Nd142)Di
## (Pr141)Di
## Number of markers: 4
## Gate profile names: 
##  Gate_1002709_Z2F0ZTE.
##  Gate_1002710_Z2F0ZTI.
```


## <a name="gate_importation"></a> 3.8. Importation of gate profiles from text files
The `import.GATE()` function imports one or several range gate profiles from a tab separated file into a `GATE` object (see section 10.4. for more details). Such tab separated file must contain for each gate profile the ranges of the expression markers and must be formatted as the following:

* each row must represent a gate profile
* each column must represent a marker
* each cell in the table must contain the marker expression lower and upper bounds for a given gate profile separated by a semicolon

The first column must contain the gate names and the first row must contain the marker names.

A typical gate file to import must look like the following:

|      | marker_1 | marker_2 | marker_3 | marker_4 | marker_5 | marker_i | marker_n |
|------|---------|---------|---------|---------|---------|---------|---------|
| gate_1 | 1.34:2.72 | 2.77:3.80 | 4.73:6.07 | 1.7:2.44  | 3.00:3.33 | 4.21:5.55 | 5.44:6.15 |
| gate_2 | 0.20:0.81 | 1.34:2.86 | 3.63:4.88 | 4.38:5.10 | 2.96:4.66 | 3.34:4.19 | 1.21:1.68 |
| gate_3 | 4.18:5.61 | 3.48:1.53 | 2.86:4.52 | 4.59:5.94 | 0.78:2.31 | 3.61:5.61 | 1.39:3.28 |
| gate_4 | 1.51:3.45 | 4.24:4.37 | 1.72:2.64 | 1.62:3.54 | 0.62:1.97 | 1.84:3.71 | 0.67:2.31 |
| gate_5 | 0.51:1.27 | 5.46:6.92 | 5.07:6.91 | 0.71:1.23 | 5.24:5.59 | 5.12:5.87 | 5.21:5.56 |

An import of a such file can be done using the following command:


```r
# imports a tab separated file containing several gate profiles 
imported.gates <- import.GATE('./GATE.TXT/gates.txt')
```

The `exclude` parameter allows to specify some markers to exclude in the import procedure. This parameter takes a character vector containing the marker names to be excluded.

A summary of the `GATE` object, imported by the `import.GATE()` function, can be done using the following command:

```r
print(imported.gates)
```

```
## Object class: GATE
## Object name: gates.txt
## Number of gate profiles: 5
## Markers: 
## marker_1
## marker_2
## marker_3
## marker_4
## marker_5
## marker_i
## marker_n
## Number of markers: 7
## Gate profile names: 
##  gate_1
##  gate_2
##  gate_3
##  gate_4
##  gate_5
```


# <a name="object_manipulation"></a> 4. Manipulation of cytometry profiles

## <a name="object_extraction"></a> 4.1. Extraction of cytometry profiles
The extract function `[i,j]` can be used to extract subsets of `CELL`, `CLUSTER`, or `GATE` objects. The parameter `i` represents a vector of profiles to extract and the parameter `j` represents a vector of markers to extract. Both `i` and `j` can be numeric, logical or character vectors.

For example, subsets of a `CELL` object can be extracted using the following commands:

```r
# extracts the first 10 cell profiles of a given CELL object
bcells.subset1 <- bm_example.cells.b[1:10]
# extracts the first 10 cell profiles and extracts 3 markers of a given CELL object
bcells.subset2 <- bm_example.cells.b[1:10,c("110-CD3"," 110_114-CD3","111-CD3")]
```

For example, subsets of a `CLUSTER` object can be extracted using the following commands:

```r
# extracts some specific cell cluster profiles of a given CLUSTER object
monoclusters.subset1 <- bm_example.clusters.mono[c(1,3,5,7,11)]
# extracts the first 10 cell profiles and extracts some markers of a given CLUSTER object 
monoclusters.subset2 <- bm_example.clusters.mono[1:10,5:10]
```

For example, subsets of a `GATE` object can be extracted using the following commands:

```r
# extracts a specific gate profile of a given GATE object
gates.subset1 <- bm_example.gates["tCD4mem"]
# extracts two gate profiles and extracts some markers of a given GATE object 
gates.subset2 <- bm_example.gates[c("tCD4naive","tCD8naive"),c(1,2,3,4,5)]
```

## <a name="object_combination"></a> 4.2. Combination of cytometry profiles
The combine function `c()` can be used to combine two or several `CELL`, `CLUSTER`, or `GATE` objects. All the different objects to combine must be of the same type.

This function is especially useful when combining cell profiles obtained from different FCS files into one single `CELL` object. 

A combination of two or more cytometry objects can be done using the following commands:

```r
# combines a set of CELL objects into a single CELL object
tcells.combine1 <- c(bm_example.cells.tCD4mem,
                           bm_example.cells.tCD4naive,
                           bm_example.cells.tCD8mem,
                           bm_example.cells.tCD8naive)

# combines a list of CELL objects into a single CELL object
tcells.list <- list(bm_example.cells.tCD4mem,
                 bm_example.cells.tCD4naive,
                 bm_example.cells.tCD8mem,
                 bm_example.cells.tCD8naive)
tcells.combine2 <- do.call("c",tcells.list)
```

## <a name="object_transformation"></a> 4.3. Transformation of cytometry profiles
The `as.CELL()` function can be used to coerce a given numeric matrix into a `CELL` object. It is to note that the numeric matrix must have its columns named (which correspond to the marker names). Such coercion can be done using the following commands:

```r
# generates a random numeric matrix
cells <- matrix(runif(100*20,min=0,max=5),nrow=100,ncol=20)
# names the column names of the numeric matrix
colnames(cells) <- paste0("marker_",as.character(1:ncol(cells))) 
# coerces a given numeric matrix into a CELL object
cell.frommatrix <- as.CELL(cells)
```

The `as.CLUSTER()` function can be used to coerce a given `CELL` object into a `CLUSTER` object. 
This function transforms the cell profiles into cell cluster profiles by computing the means, standard deviations, and densities for each marker. Such coercion can be done using the following command:

```r
# coerces a given CELL object into a CLUSTER object
cluster.fromcell <- as.CLUSTER(bm_example.cells.b)
```

The `as.GATE()` function can be used to coerce a given `CLUSTER` or a given `CELL` object into a `GATE` object.
This function transforms the cell or cell cluster profiles into gate profiles by computing the intensity ranges for each marker. By default the intensity ranges are computed based on the 0.01 and 0.99 quantiles, but can be specified by the user. Such coercions can be done using the following commands:

```r
# coerces a given CELL object into a GATE object based on the 0.01 and 0.99 quantiles
gate.fromcell1 <- as.GATE(bm_example.cells.mono)
# coerces a given CELL object into a GATE object based on the 0.05 and 0.95 quantiles
gate.fromcell2 <- as.GATE(bm_example.cells.mono, quantiles=c(0.05,0.95))
# coerces a given CELL object into a GATE object based on the minimal and maximal intensities
gate.fromcell3 <- as.GATE(bm_example.cells.mono, quantiles=c(0,1))
```

## <a name="object_export"></a> 4.4. Exportation of cytometry profiles
The `export()` function can be used to export a `CELL` object into a FCS file or export a `GATE` object into a GatingML-XML file. The exportation restore the marker expression to their value before importation procedure.

An export of a `CELL` object can be done using the following command:

```r
# exports a given CELL object to a FCS file
export(bm_example.cells.tCD4mem,"cell.fcs")
```

An export of a `GATE` object can be done using the following command:

```r
# exports the summary of a given CELL object to an GatingML-XML file
export(bm_example.gates.b,"gate.xml")
```

# <a name="object_representation"></a> 5. Representations of cytometry profiles

## <a name="cell_visualization"></a> 5.1. Visualization of cell profiles
The `plot()` function can be used to plot cell profiles of a `CELL` object, via parallel coordinates where the x-axis represents the different markers and where the y-axis represents the marker expressions. Representations are generated using the `ggplot2` library [10] and can be modified by users.

A `CELL` object, which contains one or several cell profiles, can be plotted using the following command:

```r
# plots the first cell profile of a CELL object
plot(bm_example.cells.tCD4mem[1])
```

<img src="README.figures/single_cell_plot-1.png" style="display: block; margin: auto;" />

If the `CELL` object contains several cell profiles, then all the different cell profiles will be plotted separately using the following command:

```r
# plots the first 3 cell profiles of a CELL object
plot(bm_example.cells.tCD4mem[1:3])
```

<img src="README.figures/single_cell_multiple_plot-1.png" style="display: block; margin: auto;" />

To plot different cell profiles in a same representation please to the section 5.4.

If `CELL` object has been imported via a viSNE FCS file, then biplot overview representations of the cell objects is possible using the following commands:

```r
# plots overviews of CELL objects imported from viSNE FCS files
plot(bm_example.visne[[1]],overview=TRUE)
```

<img src="README.figures/cell_overview_plot-1.png" style="display: block; margin: auto;" />

```r
plot(bm_example.visne[[2]],overview=TRUE)
```

<img src="README.figures/cell_overview_plot-2.png" style="display: block; margin: auto;" />

```r
plot(bm_example.visne[[3]],overview=TRUE)
```

<img src="README.figures/cell_overview_plot-3.png" style="display: block; margin: auto;" />

ViSNE representation aims to visualize cell profiles in a 2-dimensional space using a dimensionality reduction process. For each cell profile a 2-dimensional coordinate is calculated. The 2 dimensions are named tSNE1 and tSNE2.


## <a name="cluster_visualization"></a> 5.2. Visualization of cluster profiles
The `plot()` function can be used to plot cell cluster profiles of a `CLUSTER` object, via parallel coordinates where the x-axis represents the different markers, where the y-axis represents the marker expressions, and where error bars indicate the marker expression standard deviations. Representations are generated using the `ggplot2` library [10] and can be modified by users.
 
A `CLUSTER` object, which contains one or several cell cluster profiles, can be plotted using the following command:

```r
# plots the second cell cluster profile of a CLUSTER object
plot(bm_example.clusters[2])
```

<img src="README.figures/single_cluster_plot-1.png" style="display: block; margin: auto;" />

If the `CLUSTER` object contains several cell cluster profiles, then all the different cell cluster profiles will be plotted separately using the following command:

```r
# plots the first 3 cell cluster profiles of a CLUSTER object
plot(bm_example.clusters[1:3])
```

<img src="README.figures/single_cluster_multiple_plot-1.png" style="display: block; margin: auto;" />

To plot different cell cluster profiles in a same representation please to the section 5.4.

If a `CLUSTER` object has been imported via SPADE result files, then a SPADE tree overview representation of the cluster object is possible using the following command:

```r
# plots an overview of a CLUSTER object imported from SPADE result files
plot(bm_example.clusters,overview=TRUE)
```

<img src="README.figures/cluster_overview_plot-1.png" style="display: block; margin: auto;" />

In a SPADE tree representation, each node represents a cell cluster and similar cell clusters are linked using a minimal spanning tree approach. The sizes of the nodes are proportional to the amount of cells present within each cluster.


## <a name="gate_visualization"></a> 5.3. Visualization of gate profile
The `plot()` function can be used to plot gate profiles of a `GATE` object, via ribbons where the x-axis represents the different markers and where the y-axis represents the marker intensity ranges. Representations are generated using the `ggplot2` library [10] and can be modified by users.

A `GATE` object, which contains one or several gates profiles, can be plotted using the following command:

```r
# plots the first gate profile of a GATE object and restricts it to the 30 first markers
plot(bm_example.gates[1,1:30])
```

<img src="README.figures/single_gate_plot-1.png" style="display: block; margin: auto;" />

If the `GATE` object contains several gate profiles, then all the different gate profiles will be plotted separately using the following command:

```r
# plots the first 3 gate profiles of a GATE object
plot(bm_example.gates[1:3])
```

<img src="README.figures/single_gate_multiple_plot-1.png" style="display: block; margin: auto;" />

To plot different gate profiles in a same representation please to the section 5.4.

## <a name="pairwise_visualization"></a> 5.4. Pairwise visualization of cytometry profiles
Two cytometry profiles can be visualized in one single representation using the `plot()` function. Such representations are especially useful to better visualize the marker similarities or inclusions between two cytometry profiles. All combinations of representations, between two cytometry profiles, are possible and can be performed using the following commands: 


```r
# plots two cell profiles in one single representation
plot(bm_example.cells.b[1],bm_example.cells.mono[2])
```

<img src="README.figures/cell_cell_plot-1.png" style="display: block; margin: auto;" />


```r
# plots two cell cluster profiles in one single representation
plot(bm_example.clusters[1],bm_example.clusters[2])
```

<img src="README.figures/cluster_cluster_plot-1.png" style="display: block; margin: auto;" />


```r
# plots two gate profiles in one single representation
plot(bm_example.gates[1],bm_example.gates[2])
```

<img src="README.figures/gate_gate_plot-1.png" style="display: block; margin: auto;" />


```r
# plots a cell profile and a cell cluster profile in one single representation
plot(bm_example.cells.b[1],bm_example.clusters[2])
```

<img src="README.figures/cell_cluster_plot-1.png" style="display: block; margin: auto;" />


```r
# plots a cell profile and a gate profile in one single representation
plot(bm_example.cells.b[1],bm_example.gates[2])
```

<img src="README.figures/cell_gate_plot-1.png" style="display: block; margin: auto;" />


```r
# plots a cell cluster profile and a gate profile in one single representation
plot(bm_example.clusters[1],bm_example.gates[2])
```

<img src="README.figures/cluster_gate_plot-1.png" style="display: block; margin: auto;" />


# <a name="object_comparison"></a> 6. Comparisons of cytometry profiles 

## <a name="compare_function"></a> 6.1. Overview of the comparison approach
Cytometry profiles contained in `CELL`, `CLUSTER`, or `GATE` objects can be compared using the `compare()` function. Comparison results are stored in a `RES` object (please refer to section 7 for more details). Comparisons can be performed between profiles of same types or between profiles of different types:

* in the default statistical approach, if the comparisons are performed on profiles of same types then profiles will be compared to identify similar profiles;
* in the default statistical approach, if the comparisons are performed on profiles of different types then profiles will be compared to identify included profiles.

For each comparison of two cytometry profiles, a p-value asserting the statistical significance is provided. In the case of a comparison between profiles of same types, a similarity distance is calculated. Comparisons can be performed based on the whole set of common markers between the two profiles, or based on a subset of markers specified by the user. Moreover, markers can be weighted in the comparison procedure, via a `MWEIGHTS` object.

If only one object is provided to the `compare()` function then the comparisons will be performed between all profiles of this object. If two objects are provided to the `compare()` function then the comparisons will be performed between all possible pairs of profiles between these two objects.

Importantly, users can define their own function to perform the statistical comparisons of the profiles, using the `method` parameter (please refer to section 9 for more details). 

Different parameters can be defined, via the `method.params` named list, to specify the behavior of such kind of comparisons:

* the `D.th` parameter indicates the similarity threshold to use for comparison between profiles of same types (default values are precised in sections below);
* the `P` parameter indicates the proportion of marker successes to statistically overtake to consider the similarity or inclusion as significant.

## <a name="cell_cell_compare"></a> 6.2. Comparisons between two cell profiles
In the case of comparisons between two cell profiles, the marker distances are calculated based on the Euclidean distance (by default, the similarity threshold parameter `D.th` is here equals to 1.5).

For example, comparisons between different cell profiles can be done using the following commands:

```r
# compares some cell profiles from same CELL objects, using a specific MWEIGHTS object and using a similarity threshold of 1
res.cells_mono <- compare(bm_example.cells.mono[1:30],mweights=bm_example.mweights,method.params=list(D.th=1))
print(res.cells_mono)
```

```
## Object class: RES
## Number of comparisons: 900
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
## CELL:mono:1
## CELL:mono:2
## CELL:mono:3
## CELL:mono:4
## CELL:mono:5
## and 25 more...
```

```r
res.cells_b <- compare(bm_example.cells.b[1:30],mweights=bm_example.mweights,method.params=list(D.th=1)) 
print(res.cells_b)
```

```
## Object class: RES
## Number of comparisons: 900
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
## CELL:b:1
## CELL:b:2
## CELL:b:3
## CELL:b:4
## CELL:b:5
## and 25 more...
```


```r
# compares cell profiles from two different CELL objects, using a specific MWEIGHTS object and using a similarity threshold of 1
res.cells_mono.vs.b <- compare(bm_example.cells.mono[1:30],bm_example.cells.b[1:30],mweights=bm_example.mweights,method.params=list(D.th=1))
print(res.cells_mono.vs.b)
```

```
## Object class: RES
## Number of comparisons: 900
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
## CELL:mono:1
## CELL:mono:2
## CELL:mono:3
## CELL:mono:4
## CELL:mono:5
## and 55 more...
```


## <a name="cluster_cluster_compare"></a> 6.3. Comparisons between two cell cluster profiles
In the case of comparisons between two cell cluster profiles, the marker distances are calculated based on the Kolmogorov-Smirnov distance (by default, the similarity threshold parameter `D.th` is here equals to 0.3). The Kolmogorov-Smirnov distance corresponds to the maximum absolute distance between two cumulative distribution functions.

Additionally to `D.th` and `P`, different advanced parameters can be defined, via the `method.params` named list. Please refer to the CytoCompare R documentation to obtain more details.

For example, comparisons between different cell cluster profiles can be done using the following commands:

```r
# compares a set of cell cluster profiles from the same CLUSTER object, using a specific MWEIGHTS object
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
# compares all cell cluster profiles from the same CLUSTER object, using a specific MWEIGHTS object
res.clusters_b <- compare(bm_example.clusters.b,mweights=bm_example.mweights)
```



```r
# compares cell cluster profiles from different CLUSTER objects, using a specific MWEIGHTS object
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


## <a name="gate_gate_compare"></a> 6.4. Comparisons between two gate profiles
In the case of comparisons between two gate profiles, gates are modeled by uniform distributions, and the marker distances are calculated based on the Kolmogorov-Smirnov distance (by default, the similarity threshold parameter `D.th` is here equals to 0.3).

For example, comparisons between different gate profiles can be done using the following commands:

```r
# compares two gate profiles from different GATE objects, using a specific MWEIGHTS object
res.gates_b.vs.mono <- compare(bm_example.gates.b,bm_example.gates.mono,mweights=bm_example.mweights)
print(res.gates_b.vs.mono)
```

```
## Object class: RES
## Number of comparisons: 1
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
## GATE:b:b
## GATE:mono:mono
```


```r
# compares gate profiles from different GATE objects, using a specific MWEIGHTS object
gates     <- list(bm_example.gates.tCD4naive,bm_example.gates.mono,bm_example.gates.tCD8naive,bm_example.gates.tCD8mem)
res.gates <- RES()
for(i in 1:length(gates)){
  for(j in 1:length(gates)){
    res.gates <- c(res.gates,compare(gates[[i]],gates[[j]],mweights=bm_example.mweights))
  }
}
```


## <a name="cell_gate_compare"></a> 6.5. Comparisons between a cell profile and a gate profile
In the case of comparisons between cell profiles and gate profiles, an inclusion assessment is performed for each marker of the two profiles. A cell profile marker is considered as included in a gate profile when its expression value is within the range of the marker boundaries.

For example, a comparison between different cell profiles and gate profiles can be done using the following command:

```r
# compares a set of cell with a gate profile, using a specific MWEIGHTS object
res.gates_mono.vs.b <- compare(bm_example.cells.b[1:1000],bm_example.gates.mono,mweights=bm_example.mweights)
print(res.gates_mono.vs.b)
```

```
## Object class: RES
## Number of comparisons: 1000
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
## CELL:b:1
## CELL:b:2
## CELL:b:3
## CELL:b:4
## CELL:b:5
## and 996 more...
```


## <a name="cell_cluster_compare"></a> 6.6. Comparisons between a cell profile and a cell cluster profile
In the case of comparisons between cell cluster profiles and cluster profiles, an inclusion assessment is performed for each marker of the two profiles. A cell profile marker is considered as included in a cell cluster profile when its expression value is within the range of the marker cluster defined based on quantiles of marker expression densities.

Additionally to `P`, the `cluster.quantiles` parameter indicates the quantiles that will define the cluster profile ranges (set to c(0.10,0.90) by default).

For example, a comparison between different cell profiles and cluster profiles can be done using the following command:

```r
# compares a set of cell with cell cluster profiles, using a specific MWEIGHTS object
res.cells_b.vs.cluster_b <- compare(bm_example.cells.b[c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)],bm_example.clusters.b,mweights=bm_example.mweights)
print(res.cells_b.vs.cluster_b)
```

```
## Object class: RES
## Number of comparisons: 1175
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
## CELL:b:2
## CELL:b:3
## CELL:b:5
## CELL:b:7
## CELL:b:11
## and 67 more...
```


## <a name="cluster_gate_compare"></a> 6.7. Comparisons between a cell cluster profile and a gate profile
In the case of comparisons between cell cluster profiles and gate profiles, an inclusion assessment is performed for each marker of the two profiles. A cell cluster profile marker is considered as included in a gate profile when its expression boundaries is within the range of the marker gate.

Additionally to `D.th` and `P`, the `cluster.quantiles` parameter indicates the quantiles that will define the cluster profile ranges (set to c(0.10,0.90) by default). Different advanced parameters can also be defined, via the `method.params` named list. Please refer to the CytoCompare R documentation to obtain more details.

For example, a comparison between different cell cluster profiles and gate profiles can be done using the following command:

```r
# compares all the cell cluster profiles with gates profiles, using a specific MWEIGHTS object
res.cluster_mono.vs.gate_mono <- compare(bm_example.clusters.mono,bm_example.gates.mono,mweights=bm_example.mweights)
print(res.cluster_mono.vs.gate_mono)
```

```
## Object class: RES
## Number of comparisons: 48
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
## CLUSTER:mono:33
## CLUSTER:mono:229
## CLUSTER:mono:241
## CLUSTER:mono:30
## CLUSTER:mono:166
## and 44 more...
```


# <a name="res_object"></a> 7. Manipulation of comparison results

## <a name="res_structure"></a> 7.1. Structure of comparison results
The `RES` object is a S4 object containing one or several comparison results. This object mainly stores for each comparison result: the aggregated distance and the marker distances (in the case of similarity comparison), the associated similarity or inclusion p-value and the marker similarity or inclusion successes.

Different slots are available for a given `RES` object:

* the slot `comparisons` is a data.frame containing for each comparison: the profile names, the type of the comparison (similarity or inclusion), the similarity distance (or NA in case of inclusion), the distance threshold used (or NA in case of inclusion) and the associated p-value
* the slot `comparisons.nb` is an integer indicating the number of comparisons 
* the slot `comparison.type` is an character indicating if the comparisons have been performed by evaluating similarities or inclusions
* the slot `markers` is a character vector containing the marker names used in the comparisons
* the slot `marker.weights` is a numeric vector containing the weigths associated to each marker involded in the comparisons
* the slot `marker.distances` is a data.frame containing the marker similarity for each comparison (or NA in case of inclusion)
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
The extract function `[i]` can be used to extract comparisons of `RES` objects. The parameter `i` represents a vector of the comparisons to extract. This vector can be a numeric, logical or character.

For example, a subset of comparisons of a `RES` object can be extracted using the following command:

```r
# extracts the two first comparison of a given RES object
res.cells_mono.vs.b_subset <- res.cells_mono.vs.b[1:2]
print(res.cells_mono.vs.b_subset)
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
## CELL:mono:1
## CELL:b:1
## CELL:b:2
```

## <a name="res_combination"></a> 7.4. Combination of comparison results
The combine function `c()` can be used to combine two or several `RES` objects. RES objects can be combined to an empty RES object (i.e. `RES()`).

A combination of two or more `RES` objects can be done using the following commands:

```r
# combines a set of RES objects into a single RES object
res.cells <- c(res.cells_mono,res.cells_b,res.cells_mono.vs.b)
print(res.cells)
```

```
## Object class: RES
## Number of comparisons: 2700
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
## CELL:mono:1
## CELL:mono:2
## CELL:mono:3
## CELL:mono:4
## CELL:mono:5
## and 55 more...
```

It is to note that the comparisons to combine must be unique.

## <a name="res_visualization"></a> 7.5. Visualization of single comparison results
The `plot()` function can be used to plot the `RES` objects, via bar plots where each bar corresponds to a marker with a height proportional to the marker distances or inclusion assessments, and where the bars are colored if they model a success. Representations are generated using the `ggplot2` library [10] and can be modified by users.

Two `RES` objects, which contains one or several comparison results, can be plotted using the following commands:

```r
# plots a named comparison results based on distance measurements
plot(res.cells_b["CELL:b:3/vs/CELL:b:7"])
# plots the second comparison results based on distance measurements
plot(res.clusters_b.small[2])
# plots the second comparison results based on distance measurements
plot(res.gates[2])
```

```
## [1] 1
## [1] 0.3
## [1] 0.3
```

<img src="README.figures/res_single_plot1-1.png" style="display: block; margin: auto;" />
*In the first and second plot the cell profile and the cluster profiles are similar in contrary to the third plot where gate profiles seems to be strongly dissimilar*

A `RES` object, which contains the result of inclusion assessments between a cell profile and a gate profile, can be plotted using the following command:

```r
# plots the 8th comparison result based on inclusion assessments
plot(res.gates_mono.vs.b[8])
# plots a named comparison result based on inclusion assessments
plot(res.cells_b.vs.cluster_b["CELL:b:61/vs/CLUSTER:b:60"])
#plots the first comparison result based on inclusion assessments
plot(res.cluster_mono.vs.gate_mono[1])
```
<img src="README.figures/res_single_plot2-1.png" style="display: block; margin: auto;" />
*In the first plot, the cell profile is not included in the gate profile in contrary to the second and third plot where the cytometry profiles seems to be included in the other ones*

If the `RES` object contains several comparison results, then all the different comparison results will be plotted separately using the following command:

```r
# plots three comparison results
plot(res.cells_mono[c(2,3,4)])
```

```
## [1] 1
## [1] 1
## [1] 1
```

<img src="README.figures/res_multiple_plot-1.png" style="display: block; margin: auto;" />

## <a name="res_d3js"></a> 7.6. Visualization of comparison results with D3.js
The comparison results (i.e. a `RES` object) can also be visualized using circular graphs or using Multidimensional scaling (MDS) representations [11]. Both representations are generated as Scalable Vector Graphics (SVG) elements using the Data-Driven Documents library ([D3.js](http://d3js.org/)) [12]. Representation results are interactive and saved as HTML files. It is to note that, you can use the `webshot` R package to automatically convert these HTML files to png, jpg or pdf files.


### <a name="res_d3js_circular"></a> 7.6.a. Visualization of multiple comparison results
In the case of a circular graph, each node in the representation is a cell, a cell cluster, or a gate profile. Links between nodes represent similar profiles, based on the aggregated p-values. Only comparison with a p-value below a specific threshold are represented.

Through two sliders, users can specify the link tension in the representation and can change the p-value cutoff. In the HTML file, users can obtained more specific details about the significant similarities or inclusions by moving the mouse cursor over the links.

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
MDS methods aim to represent the similarities and differences among high dimensionality objects into a space of a lower dimensions, generally in two or three dimensions for visualization purposes [11]. In MDS representations, the Kruskal Stress (KS) indicates the amount of information lost during the dimensionality reduction process.

In the case of a MDS representation, each dot in the visualization represents a cell, a cell cluster or a gate profile. Distances between the dots are then proportional to comparison measures calculated between the profiles. It is to note that only complete `RES` objects of the same type can be represented in a MDS representation.

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

# <a name="biplot"></a> 8.1. Biplot representations of cell profiles
Biplot representation of a `CELL` object can be generated and plotted using the `biplot()` function. In such representation, each dot corresponds to a cell profile and dot are plotted in a 2-dimensional space corresponding to the marker expressions. 

A biplot representation of a `CELL` object can be displayed using the following command:

```r
# plots a biplot representation for a given CELL object using the "115-CD45" and "145-CD4" markers 
biplot(bm_example.cells.b,"110_114-CD3","139-CD45RA")
```

<img src="README.figures/biplot_plot-1.png" style="display: block; margin: auto;" />

# <a name="dheatmap"></a> 8.2. Density heatmap representations of cell cluster profiles
Density heatmap of a `CLUSTER` object can be generated and plotted using the `dheatmap()` function. In such representation, each bar corresponds to a marker and the color gradient is proportional to the marker expression density. 

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
## <a name="object_structure_uml"/> 10.1. Overview of CytoCompare objects

The following UML diagram summarizes the structure of the package:

![](README.figures/UMLDiagram.png)
 

## <a name="cell_structure"></a> 10.2. Structure of CELL object
The `CELL` object is a S4 object containing one or several cell profiles. This object mainly stores for each cell profile: the intensities of each marker.

Different slots are available for a given `CELL` object:

* the slot `name` is a character indicating the internal name of the `CELL` object
* the slot `profiles` is a character vector containing the names of the cell profiles 
* the slot `profiles.nb` is an integer value indicating the number of cell profiles
* the slot `markers` is a character vector containing the marker names
* the slot `markers.nb` is an integer value indicating the number of markers
* the slot `intensities` is a numeric matrix containing the intensities of each marker for each cell profile
* the slot `trans` is a character specifying the name of a transformation function applied on the marker expression intensities.
* the slot `trans.para` is a named list containing parameters for the transformation.
* the slot `trans.exclude` is a character vector containing the marker names for which no transformation has been applied on
* the slot `overview.function` is a character specifying the name of a function to call when plotting the `CELL` object overview (please refer to the documentation of the `plot()` function)
* the slot `layout` is a numeric matrix that can be used to store the positions of cells in a 2-dimensional space (e.g. tSNE1 and tSNE2 dimensions provided by viSNE)



##<a name="cluster_structure"></a> 10.3. Structure of CLUSTER object
The `CLUSTER` object is a S4 object containing one or several cell cluster profiles. This object mainly stores for each cell cluster profile: the means, the standard deviations and the densities of each marker.

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



## <a name="gate_structure"></a> 10.4. Structure of GATE object
The `GATE` object is a S4 object containing one or several gate profiles. This object mainly stores for each gate profile: the intensity ranges of each marker. 

Different slots are available for a given `GATE` object:

* the slot `name` is a character indicating the internal name of the `GATE` object
* the slot `profiles` is a character vector containing the names of the gate profiles
* the slot `profiles.nb` is an integer value indicating the number of cell gate profiles
* the slot `markers` is a character vector containing the marker names
* the slot `markers.nb` is an integer value indicating the number of markers
* the slot `ranges` is a 3-dimensional numeric array containing the intensity ranges of each marker for each gate profile


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
