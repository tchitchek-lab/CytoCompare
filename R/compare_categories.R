

compare_categories <- function(spade1, spade2, spade1.cluster=NULL, spade2.cluster=NULL, spade1.rename=NULL, spade2.rename=NULL, markers, method="spearman", Rlim=0.8, plim=0.05, num=5, th.min_cells = 50){
  
	results_1 <- importResultsFromSPADE(spade1,th.min_cells = th.min_cells)
	results_2 <- importResultsFromSPADE(spade2,th.min_cells = th.min_cells) 

	if (is.null(spade1.cluster)){
		spade1.cluster <- results_1@cluster.names
	}
	if (is.null(spade2.cluster)){
		spade2.cluster <- results_2@cluster.names
	}
	if (!is.character(spade1.cluster)){
		stop("Error: spade1.cluster must be a character")
	}
	if (!is.character(spade2.cluster)){
		stop("Error: spade2.cluster must be a character")
	}


	if (is.null(spade1.rename)){
		spade1.rename <- rep("1:",length(spade1.cluster))
	}
	if (length(spade1.rename) != length(spade1.cluster)){
		stop("Error: size of spade1.rename is not equal to the size of spade1.cluster")
	}
	if (is.null(spade2.rename)){
		spade2.rename <- rep("2:",length(spade2.cluster))
	} 
	if (length(spade2.rename) != length(spade2.cluster)){
		stop("Error: size of spade2.rename is not equal to the size of spade2.cluster")
	}
	if (!is.character(spade1.rename)){
		stop("Error: spade1.rename must be a character")
	}
	if (!is.character(spade2.rename)){
		stop("Error: spade2.rename must be a character")
	}

	mat1 <- computePhenoTable(results_1@cluster.phenotypes, bounds = results_1@bounds, num = num)
	mat2 <- computePhenoTable(results_2@cluster.phenotypes, bounds = results_2@bounds, num = num)

	mat1 <- reshape2::dcast(mat1, cluster ~ marker)
	mat2 <- reshape2::dcast(mat2, cluster ~ marker)

	check1 <- setdiff(markers,colnames(mat1))
	check2 <- setdiff(markers,colnames(mat2))

	if (length(check1)>0){
		stop("Error: all markers are not present within spade1")
	}

	if (length(check2)>0){
		stop("Error: all markers are not present within spade2")
	}
	mat1 <- mat1[spade1.cluster,markers]
	mat2 <- mat2[spade2.cluster,markers]

	namer1 <- c()
	namer2 <- c()
	for (k in 1:length(spade1.rename)){
		namer1[k]<- paste0(spade1.rename[k],"_",spade1.cluster[k])
	}
	for (k in 1:length(spade2.rename)){
		namer2[k]<- paste0(spade2.rename[k],"_",spade2.cluster[k])
	}

	pvalue <- c()
	correlation <- c()
	distance <- 1:length(markers)
	profil1 <- c()
	profil2 <- c()
	if ((spade1==spade2)&(spade1.cluster==spade2.cluster)){
		for (i in 1:(length(spade1.rename)-1)){
			print(paste0(i," out of ",length(spade1.rename)))
			for (j in (i+1):length(spade2.rename)){
				m1 <- unlist(mat1[i,])
				m2 <- unlist(mat2[j,])
				profil1 <- c(profil1,namer1[i])
				profil2 <- c(profil2,namer2[j])
				if (is.na(m1[1]) | is.na(m2[1])){
					correlation <-c(correlation,0)
					pvalue <- c(pvalue,1)
					distance <- rbind(distance,rep(5,length(markers)))
				} else {
					test <- stats::cor.test(m1,m2,method=method)
					manhattan <- abs(m1-m2)
					correlation <-c(correlation,test$estimate)
					pvalue <- c(pvalue,test$p.value)
					distance <- rbind(distance,manhattan)
				}
			}
		}
	} else {
		for (i in 1:length(spade1.rename)){
			print(paste0(i," out of ",length(spade1.rename)))
			for (j in 1:length(spade2.rename)){
				m1 <- unlist(mat1[i,])
				m2 <- unlist(mat2[j,])
				profil1 <- c(profil1,namer1[i])
				profil2 <- c(profil2,namer2[j])
				if (is.na(m1[1]) | is.na(m2[1])){
					correlation <-c(correlation,0)
					pvalue <- c(pvalue,1)
					distance <- rbind(distance,rep(5,length(markers)))
				} else {
					test <- stats::cor.test(m1,m2,method=method)
					manhattan <- abs(m1-m2)
					correlation <- c(correlation,test$estimate)
					pvalue <- c(pvalue,test$p.value)
					distance <- rbind(distance,manhattan)
				}
			}
		}
	}

	data <- data.frame("profile1"=profil1,"profile2"=profil2,"R"=correlation,"pvalue"=pvalue)
	if(!is.null(plim)){
		data$pvalue[data$pvalue>plim] <- 1
	}
	if(!is.null(Rlim)){
		data$pvalue[data$R<Rlim] <- 1
	}
	distance <- distance[2:nrow(distance),]
	data$profile1 <- paste0("CLUSTER:",data$profile1)
	data$profile2 <- paste0("CLUSTER:",data$profile2)
	data$measure <- rowSums(distance)
	for (i in 1:length(data$profile1)){
		rownames(distance)[i] <- paste0(data$profile1[i],"/vs/",data$profile2[i])
	}

	marker.distances = distance
	marker.successes = distance

	colnames(marker.distances) <- markers
	colnames(marker.successes) <- markers

	res <- RES(comparisons = data[,c(1,2,4,5)],
		comparisons.nb     = nrow(data),
		markers            = markers,
		marker.distances   = as.data.frame(marker.distances),
		marker.successes   = as.data.frame(marker.successes),
		marker.weights     = rep(1,length(markers)))

	return(res)
}


# @title Internal - Generates marker expression scores describing phenotypes
# 
# @description 
# This function is used internally to generate a melted numeric dataframe of discrete expression categories for each marker of each cluster.
# 
# @details 
# The function calculates the mean of median expressions between samples (NA values are removed).
# The function assign calculated mean of median expressions to a category (the 'num' parameter provides the number of categories) between bounds of marker expressions.
# The bound parameter must contain the marker name in colnames and two rows. 
# For each marker, first row is the lower bound value et second row is the upper bound value.
# The resulting matrix of this function contains 3 columns: "cluster", "marker" and "value".
# 
# @param cluster.phenotypes a dataframe containing the marker median expressions for each cluster of each sample
# @param bounds a dataframe containing the bounds for each marker
# @param num a numeric value specifying the number of markers expression categories
#  
# @return a numeric matrix of expression scores
# 
#' @import plyr
computePhenoTable <- function(cluster.phenotypes, bounds, num = 5){

    cluster.phenotypes        <- stats::na.omit(cluster.phenotypes)
    cluster.phenotypes.melted <- reshape2::melt(cluster.phenotypes, id.vars = c("sample", "cluster"))
    
    colnames(cluster.phenotypes.melted) <- c("sample", "cluster", "marker", "value")
    cluster.phenotypes.melted$marker    <- as.vector(cluster.phenotypes.melted$marker)
    means                               <- plyr::ddply(cluster.phenotypes.melted,c("cluster", "marker"),function(df){mean(df$value, na.rm = TRUE)})
    
    colnames(means)       <- c("cluster", "marker", "value")

    for (i in seq_len(nrow(means))) {
        
        cluster <- means[i, "cluster"]
        value   <- means[i, "value"]

        min     <- bounds[1, means[i, "marker"]]
        max     <- bounds[2, means[i, "marker"]]

        seq               <- seq(from = min, to = max, length.out = num)
        means[i, "value"] <- which.min(abs(value - seq))
        
    }

    return(means)
}


# @title Importation of clustering results generated by SPADE
#
# @description 
# The 'importResultsFromSPADE()' function imports SPADE cell clustering results from a specified path.
# This function imports the cluster phenotype matrix and count matrix as well as the SPADE tree.
# This function apply an hyperbolic sine transformation to imported FCS data (unless 'use.raw.medians' = TRUE) and compute the marker range quantiles.
# 
# @details
# This function returns a 'SPADEResults' object including 'flowset', 'fcs.files', 'graph' and 'graph.layout' slots. 
# The computation of marker range quantiles can be approximated using 'quantile.approximation' parameter which is more efficient in term of loading time and memory usage.
#  
# @param path a character specify the path of SPADE results folder
# @param exclude.markers a character vector of markers to exclude (case insensitive)
# @param probs a vector of probabilities with 2 values in [0,1] to compute marker range quantiles. First is the lower bound and second is the upper bound.
# @param use.raw.medians a logical specifying if arcsinh transformed or raw medians will be used in the cluster expression matrix (FALSE by default)
# @param quantile.approximation a logical specifying if marker range quantiles are computed using all cells (FALSE), or is the means of the quantile of each samples (TRUE)
# @param th.min_cells a numeric specifying the minimum number of cell in a cluster of a sample to take its phenotype in account
# @param load.phenotype a logical specifying if the phenotype matrix and fcs file will be loaded
#
# @return a S4 object of class 'SPADEResults'
#
# @export 
#
#' @import igraph
importResultsFromSPADE <- function(path,
                                   exclude.markers        = c("cell_length", "FileNum", "density", "time"),
                                   probs                  = c(0.05, 0.95),
                                   use.raw.medians        = FALSE,
                                   quantile.approximation = FALSE,
                                   th.min_cells           = 0,
                                   load.phenotype         = TRUE){
    
    message("[START] - importing SPADE clustering results")
    path <- normalizePath(path, "/", mustWork = TRUE)
    
    message(paste0(basename(path), "\n"))
    
    if (typeof(exclude.markers) != "character") {
        stop("Error in importResultsFromSPADE: The 'exclude.markers' parameter must be a character vector")
    }
    
    if (length(probs) != 2) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must only 2 numeric values")
    } else if (probs[1] > probs[2]) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must contain a first value greather than the second value")
    } else if (probs < 0 || probs > 1) {
        stop("Error in importResultsFromSPADE: The 'probs' parameter must contain values included in the domain: [0;1]")
    }
    
    if (!is.logical(use.raw.medians)) { stop("Error in importResultsFromSPADE: The 'use.raw.medians' parameter must be a logical") }
    if (!is.logical(quantile.approximation)) { stop("Error in importResultsFromSPADE: The 'quantile.approximation' parameter must be a logical") }
    
    if (th.min_cells < 0) { stop("Error in importResultsFromSPADE: The 'th.min_cells' parameter must be a stricly positive integer")  }
    
    if (load.phenotype) {
    
        fcs.files  <- dir(path, full.names = TRUE, pattern = ".fcs.density.fcs.cluster.fcs$")
        list       <- load.flowSet(fcs.files = fcs.files, exclude.markers = exclude.markers, use.raw.medians = use.raw.medians)
        flowset    <- list$flowset
        dictionary <- list$dictionary
        
        message("\tcompute quantiles bounds...")
        
        if (quantile.approximation) {
            quantiles <- computeQuantile.approximation(flowset, probs)
        } else {
            quantiles <- computeQuantile(flowset, probs)
        }
        gc()
    } else {
        fcs.files <- character(0)
        flowset   <- NULL
        quantiles <- data.frame()
    }
    message("\treading SPADE results...")
    
    files <- dir(paste(path, "/tables/bySample/", sep = ""), full.names = TRUE)
    if(length(files)==0)
		stop("Error when importing cell cluster abundances. Please check that subfolder \"./tables/bySample/\" is well existing")
 
    cluster.phenotypes <- data.frame(stringsAsFactors = FALSE)
    cluster.abundances <- data.frame()
    sample.names       <- c()
    for (file in files) {
        
        name                 <- gsub(".fcs.density.fcs.cluster.fcs.anno.Rsave_table.csv$", "", basename(file))
        sample.names         <- c(sample.names, name)
        SPADE.matrix         <- utils::read.table(file, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
        SPADE.matrix[, "ID"] <- as.character(SPADE.matrix[, "ID"])
        
        cluster.abundances.sample <- SPADE.matrix[, "count"]
        
        if (nrow(cluster.abundances)) {
            cluster.abundances <- cbind(cluster.abundances, cluster.abundances.sample)
            samples.headers    <- append(samples.headers, name)
        } else {
            cluster.abundances     <- data.frame(row.names = SPADE.matrix[, "ID"], cluster.abundances.sample)
            samples.headers <- name
        }
        if (load.phenotype) {
            cluster.phenotypes.sample <- SPADE.matrix[, grep("count|percenttotal", colnames(SPADE.matrix), invert = TRUE)]
            
            cluster.phenotypes.sample[cluster.abundances.sample < th.min_cells, 2:ncol(cluster.phenotypes.sample)] <- rep(NA, ncol(cluster.phenotypes.sample) - 1)
            
            cluster.phenotypes.sample <- cbind(name = rep(name, nrow(cluster.phenotypes.sample)), cluster.phenotypes.sample)
            cluster.phenotypes <- rbind(cluster.phenotypes, cluster.phenotypes.sample)
        }
    }
    
    colnames(cluster.abundances) <- samples.headers
    
    if (load.phenotype) {
        cluster.phenotypes.header    <- colnames(cluster.phenotypes)
        
        cluster.phenotypes.header[1] <- "sample"
        cluster.phenotypes.header[2] <- "cluster"
        colnames(cluster.phenotypes) <- cluster.phenotypes.header
        
        cluster.phenotypes        <- filter.medians(cluster.phenotypes, use.raw.medians)
        
        cluster.phenotypes.header <- colnames(cluster.phenotypes)
        clustering.markers.index  <- grep("_clust", cluster.phenotypes.header)
        cluster.phenotypes.header <- gsub("_clust", "", cluster.phenotypes.header)
        
        colnames(cluster.phenotypes) <- rename.markers(cluster.phenotypes.header, dictionary = dictionary)
        clustering.markers           <- colnames(cluster.phenotypes)[clustering.markers.index]
        
        if (!is.null(exclude.markers)) {
            cluster.phenotypes <- exclude.markers.spade(cluster.phenotypes, exclude.markers)
            clustering.markers <- setdiff(clustering.markers, exclude.markers)
        }
        
    	graph        <- igraph::read.graph(paste(path, "/mst.gml", sep = ""), format = "gml")
        graph.layout <- as.matrix(utils::read.table(paste0(path, "/layout.table"), sep = " ", quote = "", stringsAsFactors = FALSE))
        
        markers.names <- colnames(cluster.phenotypes[, -c(1, 2)])
    } else {
        clustering.markers <- character(0)
        markers.names      <- character(0)
        graph.layout       <- NULL
        graph              <- NULL
    }
    
    res <- methods::new("SPADEResults",
                        cluster.phenotypes = cluster.phenotypes,
                        use.raw.medians    = use.raw.medians,
                        cluster.abundances = cluster.abundances,
                        cluster.names      = rownames(cluster.abundances),
                        sample.names       = sample.names,
                        marker.names       = markers.names,
                        clustering.markers = clustering.markers,
                        cluster.number     = nrow(cluster.abundances),
                        flowset            = flowset,
                        fcs.files          = fcs.files,
                        bounds             = quantiles,
                        th.min_cells       = th.min_cells,
                        graph.layout       = graph.layout,
                        graph              = graph)
    
    message("[END] - importing SPADE clustering results")
    
    return(res)
    
}


# @title Internal - Computation of quantile with FCS flowset marker by marker 
#
# @description 
# This function is used internally to compute the marker range quantiles.
# 
# @details 
# This function performs the exact calculation of quantiles with all cells but needs more resources (time and memory usage) than 'computeQuantile.approximation'.
# 
# @param flowset a flowCore flowset
# @param probs a numeric vector of 2 values specifying the quantiles to compute
# 
# @return a numeric matrix containing the quantiles of each marker
computeQuantile <- function(flowset, probs = c(0.05,0.95)){

    bounds  <- data.frame()
    markers <- flowset@colnames
    markers <- setdiff(markers, "cluster")
    
    for (marker in markers) {
        temp <- c()
        for (sample in 1:length(flowset)) {    
            frame <- flowset[[sample]]@exprs
            temp  <- c(temp, frame[, marker])
        }
        if (nrow(bounds) > 0) {
            bounds <- cbind(bounds, stats::quantile(temp, probs = probs))
        }else{
            bounds <- as.data.frame(stats::quantile(temp, probs = probs))
        }
    }
    
    colnames(bounds) <- markers
    rownames(bounds) <- probs
    
    return(bounds)
}


# @title Internal - Computation of quantiles with FCS flowset sample by sample
#
# @description 
# This function is used internally to compute bound using the 2 selected percentiles from each sample. 
# For the first percentile, the mean of this percentile from all samples is compute (same for the second percentiles). 
# This approach avoid to compute percentile based the whole dataset and then seed up computation.
# 
# @details 
# This function performs an approximate calculation of quantiles using less memory than computeQuantile.
# 
# @param flowset a flowCore flowset 
# @param probs a numeric vector of 2 values specifying the quantiles to compute
# 
# @return a numeric matrix containing the quantiles of each marker
# 
#' @importFrom flowCore fsApply
computeQuantile.approximation <- function(flowset,probs = c(0.05,0.95)){
    
    bounds.by.sample <- flowCore::fsApply(flowset[, flowset@colnames != "cluster"], flowCore::each_col, stats::quantile, probs = probs)
    
    lower.bounds <- bounds.by.sample[seq(from = 1, to = nrow(bounds.by.sample), by = 2), ]
    upper.bounds <- bounds.by.sample[seq(from = 2, to = nrow(bounds.by.sample), by = 2), ]
    
    lower.bounds <- apply(lower.bounds, 2, mean)
    upper.bounds <- apply(upper.bounds, 2, mean)
    
    bounds <- data.frame(row.names = probs, rbind(lower.bounds, upper.bounds), check.names = FALSE)
    
    return(bounds)
    
}


# @title Internal - Filtering of medians from a SPADE result matrix
#
# @description 
# This function is used internally to remove raw or transform medians from a SPADE result matrix. CVS medians are always removed.
# 
# @param data a SPADE matrix
# @param use.raw.medians a logical specifying if "transformed" or "raw" medians will be use (FALSE by default)
# 
# @return a numeric matrix without the cell markers to exclude
filter.medians <- function(data,use.raw.medians = FALSE){
    
    if (use.raw.medians) {
        exclude <- "^medians|^cvs"
    }else{
        exclude <- "^raw_medians|^cvs"
    }
    
    data           <- data[,grep(exclude, colnames(data), invert = TRUE, ignore.case = TRUE)]
    colnames(data) <- gsub("^medians|^cvs|^raw_medians", "", colnames(data))
    
    return(data)
    
}


# @title Loading of FCS files object into a 'Results' object
#
# @description 
# This function loads the FCS files to the 'flowset' slot of the 'Results' object.
# 
# @details
# If a 'Results' object is provided, others parameters ('fcs.files', 'exclude.markers' and 'use.raw.medians') will be ignored (they will be retrieved from the 'Results' object). 
# 
# @param Results a Results object (with 'fcs.files' slot not null) (optional)
# @param fcs.files a character vector containing the absolute path of the original FCS files
# @param exclude.markers a character vector of markers to exclude (case insensitive)
# @param use.raw.medians a logical specifying if the arcsinh transformation must be performed or not
# @param pattern a character specifying the pattern of the FCS file 
# 
# @return a S4 'flowSet' object
#
# @export 
# 
#' @importFrom flowCore read.flowSet arcsinhTransform transformList transform
load.flowSet <- function(Results = NULL,
                         fcs.files,
                         exclude.markers,
                         use.raw.medians,
                         pattern = ".fcs.density.fcs.cluster.fcs") {
    
    message("FCS files loading:")
    
    if (!is.null(Results)) {
        if (!is.null(Results@fcs.files)) {
            fcs.files  <- Results@fcs.files
        } else {
            stop("Error in load.flowSet: The 'Results' parameter required a 'Results' object with a not null 'fcs.files' slot")
        }
    }
    
    flowset <- flowCore::read.flowSet(fcs.files, emptyValue = TRUE)
    
    samples.names                  <- gsub(pattern, "", basename(fcs.files))
    flowCore::sampleNames(flowset) <- samples.names
    
    dictionary       <- flowset[[1]]@parameters@data[, c(1, 2)]
    dictionary[, 1]  <- make.names(dictionary[, 1])
    
    dictionary[is.na(dictionary[, 2]),2] <- dictionary[is.na(dictionary[, 2]),1] 
	flowset@colnames <- rename.markers(flowset@colnames, dictionary = dictionary)

    if (!is.null(Results)) {
        exclude.markers  <- setdiff(flowset@colnames, c(Results@marker.names, "cluster"))
    }
    
	if (!is.null(exclude.markers)) {
        flowset <- exclude.markers.spade(flowset, exclude.markers, colnames.FCS = flowset@colnames)
    }
    
	if ((is.null(Results) && !use.raw.medians) || ((!is.null(Results)) && !Results@use.raw.medians)) {
        message("\tarchsin transform...")
        
        transform.arcsinh <- flowCore::arcsinhTransform(a = 0, b = 0.2)
        
        marker.toTransform <- setdiff(flowset@colnames, "cluster")
		transformations    <- flowCore::transformList(marker.toTransform, transform.arcsinh)
		flowset            <- flowCore::transform(flowset, transformations)
    }
    
    return(list(flowset = flowset, dictionary = dictionary))
    
}


# @title Internal - Removing of cell markers to exclude from a matrix
#
# @description 
# This function is used internally to remove one or several cell markers.
# 
# @details 
# If the data parameter is a dataframe the colnames.FCS parameter is ignored but if the data parameter is a flowset, the colnames.FCS parameter is required.
# 
# @param data a numeric matrix or flowset
# @param exclude a character vector containing the cell markers to be excluded (case intensive)
# @param colnames.FCS a character vector containing column names if data is a FCS flowset
# 
# @return a numeric matrix without the cell markers to exclude
exclude.markers.spade <- function(data, exclude, colnames.FCS = NULL){
    
    if (!is.null(colnames.FCS)) {
        column <- colnames.FCS
    } else {
        column <- colnames(data) 
    }
    
    exclude.flags <- toupper(exclude) %in% toupper(column)
    
    if (any(!(exclude.flags))) {
        warning(paste0("Unknown marker to exclude: ", paste(exclude[!exclude.flags], collapse = ", ")))
    }
    
    data    <- data[ , -which(toupper(column) %in% toupper(exclude))]
    
    return(data)
}


#' @title SPADEResults class definition
#' 
#' @description 
#' The SPADEResults object is a S4 object containing cell clustering results. 
#' 
#' This object mainly stores the cluster abundance matrix (i.e. the number of cells associated to each sample for each cluster) and the cluster phenotypes matrix (i.e. the median expressions for each marker of each cluster).
#'  
#' In addition, this object can contain information about clustering results, such a SPADE tree. 
#' 
#' @details 
#' The 'cluster.abundances' dataframe contains the number of cells associated to each sample for each cluster.
#' This dataframe stores the clusters in rows and the samples in columns.
#' 
#' The 'cluster.phenotypes' dataframe stores the median expressions for each marker of each cluster.
#' This dataframe stores in the first column the sample names, in the second column the cluster names, and in the others columns the maker median expressions.
#' 
#' The 'bounds' dataframe contains the marker expressions boundaries (minimum and maximum, or specific percentiles) for each marker.
#' 
#' The 'print()' and 'show()' can be used to display a summary of this object. 
#' 
#' @slot cluster.abundances a dataframe containing the number of cells associated to each sample for each cluster
#' @slot cluster.phenotypes a dataframe containing the median expressions for each marker of each cluster
#' @slot sample.names a character vector containing the sample names
#' @slot cluster.names a character vector containing the cluster names
#' @slot cluster.number a numeric specifying the number of clusters
#' @slot marker.names a character vector containing the marker names
#' @slot clustering.markers a character vector specifying the markers that have been used by the clustering algorithms
#' @slot bounds a numeric data.frame containing the marker expressions boundaries for each marker
#' @slot use.raw.medians a logical specifying if the marker expressions correspond to raw or transformed data
#' @slot flowset a flowSet object containing the imported SPADE FCS files
#' @slot fcs.files a character vector containing the location of the imported FCS files
#' @slot graph a igraph object containing the SPADE tree structure
#' @slot graph.layout a numeric matrix containing the SPADE tree layout
#' @slot assignments a dataframe containing annotations for each sample samples such as a biological condition ("bc"), a timepoint condition ("tp") or an individual ("ind") assignment 
#' @slot th.min_cells a numeric specifying the minimal number of cells that a cluster for a given samples needs to have to be taken into consideration in its phenotypical characterization
#' 
#' @import igraph methods
#' 
#' @name SPADEResults-class
#' @rdname SPADEResults-class
#' @exportClass SPADEResults
SPADEResults <- setClass("SPADEResults",
    slots = c(cluster.abundances = "data.frame",
        cluster.phenotypes = "data.frame",
        sample.names       = "character",
        cluster.names      = "character",
        cluster.number     = "numeric",
        marker.names       = "character",
        clustering.markers = "character",
        bounds             = "data.frame",
        use.raw.medians    = "logical",
        flowset            = "ANY",
        fcs.files          = "character",
        graph              = "ANY",
        graph.layout       = "ANY",
        assignments        = "ANY",
        th.min_cells       = "numeric"),
    validity = function(object){
        if ((length(object@marker.names) != 0) && (length(object@marker.names) + 2) != ncol(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: marker.names length (",length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistent with cluster.phenotypes size (number of columns : ",ncol(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if (nrow(object@cluster.abundances) != object@cluster.number){
            message(paste0("Error in Results object: cluster.number (",object@cluster.number,") is inconsistent with cluster.abundances matrix size (",nrow(object@cluster.abundances),")"))
            return(FALSE)
            }
        if ((length(object@marker.names) != 0) && nrow(object@bounds) != 2) {
            message(paste0("Error in Results object: bounds number of rows (",nrow(object@bounds),") is incorrect (only 2 rows accepted)"))
            return(FALSE)
            }
        if ((length(object@marker.names) != 0) && ncol(object@bounds) != length(object@marker.names)) {
            message(paste0("Error in Results object: bounds number of columns (",ncol(object@bounds),") is inconsistent with marker.names length (",
            length(object@marker.names), ")"))
			message("It is likely that automatic gating results were generated using FCS containing different cell markers.")
			message("Please use FCS files containing identical set of cell markers.")
            return(FALSE)
        }
        if (object@th.min_cells < 0) {
            message(paste0("Error in Results object: th.min_cells must be positif"))
            return(FALSE)
        }
        if (!is.null(object@assignments) && !is.data.frame(object@assignments)) {
            message(paste0("Error in Results object: assignments must be a dataframe"))
            return(FALSE)
        }
        if (!is.null(object@assignments) && (!all(rownames(object@assignments) %in% object@sample.names))) {
            message(paste0("Error in Results object: assignments must contains all samples in rownames.\n These ones are missing:",paste(setdiff(object@sample.names, rownames(object@assignments)), collapse = " ")))
            return(FALSE)
        }
        if (!is.null(object@flowset) && (class(object@flowset)[1] != "flowSet")) {
            message("Error in Results object: flowset must be of class flowSet or null")
            return(FALSE)
        }
        if (!is.null(object@fcs.files) && (class(object@fcs.files)[1] != "character")) {
            message("Error in Results object: fcs.files must be a character vector or null")
            return(FALSE)
        }
        if (!is.null(object@graph) && (class(object@graph)[1] != "igraph")) {
            message("Error in Results object: graph must be of class igraph or null")
            return(FALSE)
        }
        if (!is.null(object@graph.layout) && (class(object@graph.layout)[1] != "matrix")) {
            message("Error in Results object: graph.layout must be a matrix or null")
            return(FALSE)
        }
        if ((length(object@marker.names) != 0) && (length(object@sample.names) * object@cluster.number) != nrow(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: sample.names length (",length(object@sample.names),") and cluster.number (",object@cluster.number,") are inconsistent with cluster.phenotypes size (number of row : ",nrow(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if ((length(object@marker.names) != 0) && (length(object@marker.names) + 2) != ncol(object@cluster.phenotypes)) {
            message(paste0("Error in Results object: marker.names length (",length(object@marker.names)," + 2 for cluster IDs and sample names) are inconsistent with cluster.phenotypes size (number of columns : ",ncol(object@cluster.phenotypes),")"))
            return(FALSE)
        }
        if (nrow(object@cluster.abundances) != object@cluster.number) {
            message(paste0("Error in Results object: cluster.number (",object@cluster.number,") is inconsistent with cluster.abundances matrix size (",nrow(object@cluster.abundances),")"))
            return(FALSE)
        }
        if (ncol(object@cluster.abundances) != length(object@sample.names)) {
            message(paste0("Error in Results object: number of samples (",length(object@sample.names),") is inconsistent with cluster.abundances matrix size (",ncol(object@cluster.abundances),")"))
            return(FALSE)
        }
        if (length(object@clustering.markers) > length(object@marker.names)) {
            message(paste0("Error in Results object: clustering.markers length (",length(object@clustering.markers),") can not be higher than marker.names length (",length(object@marker.names),")"))
            return(FALSE)
        }
        if (!all(object@clustering.markers %in% object@marker.names)) {
            message(paste0("Error in Results object: clustering.markers must contains markers included in marker.names (",setdiff(object@clustering.markers, object@marker.names),")"))
            return(FALSE)
        }
        if (length(object@clustering.markers) == 0) {
            warning("Warning in Results object: clustering.markers length is 0")
        }
        for (fcs.file in object@fcs.files) {
            if (!file.exists(fcs.file)) {
                message(paste0("Error in Results object: FCS file not exist :", fcs.file))
                return(FALSE)
            }
        }
        return(TRUE)
    }
)
