# title Internal - Renaming cell markers
#
# @description This function is used internally to rename the cell markers based on a dictionary.
#
# @details dictionary is a data.frame used to rename the marker names. The first column must correspond to the original marker names, the second column must correspond to the new marker names. 
#
# @param header a character vector containing the original maker names
# @param dictionary a character vector containing a correspondence between the original and the new marker names
#
# @return a character vector containing the renamed marker names
rename.markers <- function(header,dictionary){
    header <- make.names(header)
    dictionary[,1] <- as.vector(dictionary[,1])
    dictionary[,2] <- as.vector(dictionary[,2])
    if(length(unique(dictionary[,1]))!=length(dictionary[,1])){
        stop("Duplicate in dictionary 'original marker names'")
    }
    occurences <- table(dictionary[,2])
    occurences <- occurences[occurences>1]
    redondant.names <- names(occurences)
    for (i in 1:nrow(dictionary)) {
        if (any(dictionary[i, 2] %in% redondant.names)) {
            temp             <- gsub("X\\.(.*)\\.Di(_clust$|$)","(\\1)Di\\2",dictionary[i,1])
            dictionary[i,2] <- paste0(dictionary[i,2],"-",temp)
        }
        header[which(header == dictionary[i,1])[1]] <- dictionary[i,2]
    }
    return(header)
}


# title Internal - Removing of cell markers to exclude from a matrix
#
# @description This function is used internally to remove one or several cell markers from a numeric matrix.
#
# @param data a numeric matrix
# @param exclude a character vector containing the cell markers to be excluded
#
# @return a numeric matrix without the cell markers to exclude
exclude.markers <- function(data,exclude){
    exclude.flags <- exclude %in% colnames(data)
    if(any(!(exclude.flags))){
        warning(paste0("Unknown marker to exclude: ",paste(exclude[!exclude.flags],collapse=", ")))
    }
    data    <- data[,!(colnames(data) %in% exclude)]
    return(data)
}
 

#' @title Importation of cell profiles from a tab separated file
#'
#' @description Imports one or several cell profiles from a tab separated file into a CELL object.
#'
#' @details Tab separated file to import must contain for each cell profile the intensities of the marker and must be formatted as the following:\cr
#' * each row must represent a cell profile;\cr
#' * each column must represent a marker;\cr
#' * each cell in the table must contain the marker expression intensities for a given cell profile;\cr
#' The first column must contain the cell names and the first row must contain the marker names.
#'
#' @param file a character indicating the location of a tab separated file to import
#' @param dictionary a two-column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#'
#' @return a S4 object of class CELL
#'
#' @export
import.CELL <- function(file,
                        dictionary = NULL,
                        exclude    = NULL){
    
    message(paste0("Importing ",file))
    if(!file.exists(file))
        stop(paste0("Error in import.CELL: ",file," does not exist"))
    
    data     <- utils::read.table(file,header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
    profiles <- rownames(data)
    data     <- as.matrix(data)
    
    if(!is.null(dictionary)){
        colnames(data) <- rename.markers(colnames(data),dictionary)
    }
    if(!is.null(exclude)){
        data <- exclude.markers(data,exclude)
    }
    markers        <- colnames(data)
    dimnames(data) <- NULL
    
    cell <- CELL(name = basename(file),
        profiles      = profiles,
        profiles.nb   = nrow(data),
        markers       = markers,
        markers.nb    = length(markers),
        intensities   = data)
        
    return(cell)
}

#' @title Importation of cell profiles from one or several FCS files
#' 
#' @description Imports one or several cell profiles from a FCS file or from a set of FCS files into a CELL object.
#'
#' @details If a set of files is specified, then the files are merged in the import procedure.
#'
#' Several transformations can be applied on the expression marker intensities via the 'trans' parameter:
#' 
#' The transformation functions can be parametrized using the named list 'trans.para'. The scale (cofactor) of the arcsinh transformation function can be parametrized using the 'arcsinh.scale' value. The shift of the log transformation function can be parametrized using the 'log.shift' value and the base of the log transformation function can be parametrized using the 'log.base' value. The 'log.shift' allows the value "auto" which automatically identify the log shift avoiding to apply log transformations on negative values.
#' 
#' The 'rescale' parameter can be used to rescale the marker expression intensities from 0 to 1 (with respect to the distribution proportion). The rescaling can be performed based on minimal and maximal expression values or based on specified quantiles using the 'rescale.quantiles' parameter. This strategy is especialy usefull when comparing cell or cell cluster profiles obtained from different experimental/staining conditions. 
#'    
#' @importFrom flowCore exprs read.FCS
#'
#' @param path a character vector indicating the location to a FCS file or to a set of FCS files
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#' @param trans a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation
#' @param trans.para a named list containing parameters for the transformation. Please refer to the details section for more details
#' @param trans.exclude a character vector containing the marker names for which no transformation must be applied on (including the rescaling transformation)
#' @param rescale a logical specifying if marker expression intensities must be rescale between 0 and 1
#' @param rescale.quantiles a numeric vector of two values specifying the quantiles of the marker expression intensities used to rescale
#'
#' @return a S4 object of class CELL
#'
#' @export
import.FCS <- function(path,
                        exclude           = NULL,
                        trans             = "arcsinh",
						trans.para	      = switch(trans,
												   "arcsinh" = list(arcsinh.scale=5),
												   "log"     = list(log.shift="auto",log.base=10),
												   "none"    = NULL),
                        trans.exclude     = "cluster",
                        rescale           = FALSE,
                        rescale.quantiles = c(0,1)){
    
    if(length(path)==0){
        stop("Error in import.FCS: path has a length of 0")
    }else if(length(path)==1 && is.na(file.info(path)$isdir)){
        stop(paste0("Error in import.FCS: ",path," do not exist"))
    }else if(length(path)>1 && !any(file.exists(path))){
        stop(paste0("Error in import.FCS: ",file[!any(file.exists(path))]," do not exist"))
    }
    
    if(length(path)==1 && file.info(path)$isdir==TRUE){
        path <- list.files(path,full.names = TRUE,pattern=".+\\.fcs$")
    }
    
    markers    <- c()
    cells      <- list()
    dictionary <- extract.dictionary(path[1])
    for(i in 1:length(path)){
        message(paste0("Importing ",path[i]))
        suppressWarnings(data <- flowCore::exprs(flowCore::read.FCS(path[i],trans = FALSE)))
        colnames(data)        <- rename.markers(colnames(data),dictionary)
        if(!is.null(exclude)){
            data <- exclude.markers(data,exclude)
        }
        markers <- colnames(data)[!(colnames(data) %in% trans.exclude)]
        
        if(trans=="arcsinh"){
            for(marker in markers){
                data[,marker] <- arcsinh(data[,marker],trans.para$arcsinh.scale)
            }    
        }else if(!is.element(trans,c("log","none"))){
            stop("Unknown transformation to apply")
        }
        
        markers  <- colnames(data)
        
        if(nrow(data)>1)
            profiles <- as.character(1:nrow(data))
        else    
            profiles <- as.character(NA)
        
        rownames(data) <- profiles

        cells[[i]]     <- as.CELL(data)
        
    }
    cell  <- do.call("c",cells)
    
    if(trans=="log"){
        to.transform <- match(setdiff(markers,trans.exclude),markers)
        if(trans.para$log.shift=="auto"){
            min                  <- min(cell@intensities[,to.transform])
            trans.para$log.shift <- ifelse(min>0,0,abs(floor(min-1)))
        }
        cell@intensities[,to.transform] <- log(cell@intensities[,to.transform] + trans.para$log.shift,trans.para$log.base)
    }
	
	if(rescale){
	    to.transform     <- match(setdiff(markers,trans.exclude),markers)
		cell@intensities[,to.transform] <- apply(cell@intensities[,to.transform], 2, rescale.intensities, rescale.quantiles)
	}

	cell@name          <- paste(basename(path),sep="",collapse=";")
    cell@trans         <- trans
    cell@trans.para    <- trans.para
    cell@trans.exclude <- trans.exclude
    return(cell)
}


#' @title Importation of cell profiles from viSNE FCS files
#'
#' @description Imports viSNE FCS files containing one or several cell profiles into a CELL object.
#'
#' @details ViSNE is a dimensionality reduction algorithm designed for analysis and visualization of high-dimensionality cytometry data (PMID:23685480). In a viSNE map, each dot of the representation corresponds to a cell profile in a two-dimensional space (tSNE1 and tSNE2 dimensions).
#' 
#' The 'rescale' parameter can be used to rescale the marker expression intensities from 0 to 1 (with respect to the distribution proportion). The rescaling can be performed based on minimal and maximal expression values or based on specified quantiles using the 'rescale.quantiles' parameter. This strategy is especialy usefull when comparing cell or cell cluster profiles obtained from different experimental/staining conditions. 
#' 
#' @param path a character vector indicating the location to a viSNE FCS file or to a set of viSNE FCS files
#' @param dictionary a two-column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#' @param trans a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation
#' @param trans.para a named list containing parameters for the transformation. Please refer to the details section for more details
#' @param trans.exclude a character vector containing the marker names for which no transformation must be applied on
#' @param tSNE1 a character indicating the marker name of the first viSNE dimension (tSNE1)
#' @param tSNE2 a character indicating the marker name of the second viSNE dimension (tSNE2)
#' @param rescale a logical specifying if marker expression intensities must be rescale between 0 and 1
#' @param rescale.quantiles a numeric vector of two values specifying the quantiles of the marker expression intensities used to rescale
#'
#' @return a S4 object of class CELL
#'
#' @export
import.VISNE <- function(path,
                        dictionary        = NULL,
                        exclude           = NULL,
                        trans             = "arcsinh",
						trans.para	      = switch(trans,
												   "arcsinh" = list(arcsinh.scale=5),
												   "log"     = list(log.shift="auto",log.base=10),
												   "none"    = NULL),
                        trans.exclude     = NULL,
                        tSNE1             = "tSNE1",
                        tSNE2         	  = "tSNE2",
                        rescale           = FALSE,
                        rescale.quantiles = c(0,1)){
                         
    cell <- import.FCS(path,
        trans             = trans,
        trans.para        = trans.para,
        trans.exclude     = unique(c(trans.exclude,"cluster")),
        exclude           = exclude,
		rescale           = rescale,
		rescale.quantiles = rescale.quantiles)
    
    cell@layout   <- cell[,c(tSNE1,tSNE2)]@intensities
    cell@overview.function  <- "ggplot.viSNEmap"
    
    return(cell)
}


#' @title Importation of cell cluster profiles from a tab separated file
#'
#' @description Imports one or several cell cluster profiles from a tab separated file into a CLUSTER object. In this case, the marker expressions of each cluster are assumed to be normally distributed.
#'
#' @details Tab separated file to import must contain for each cell cluster profile the means and the standard deviations of the expression markers and must be formatted as the following:\cr
#' * each row must represent a cell cluster profile;\cr
#' * each column must represent a marker;\cr
#' * each cell in the table must contain the marker expression means and the standard deviations for a given cell cluster separated by a semicolon;\cr
#' The first column must contain the cell cluster names and the first row must contain the marker names.
#'
#' It is to note that `CLUSTER` objects constructed via the `import.CLUSTER()` function do not contain the densities of expression markers (please refer to the documentation of the `compare()` function).
#'
#' @importFrom igraph graph.empty layout.auto
#'
#' @param file a character specifying the location of a tab separated file to import
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#'
#' @return a S4 object of class CLUSTER
#'
#' @export
import.CLUSTER <- function(file,
                            exclude    = NULL){
    
    message(paste0("Importing ",file))
    
    if(!file.exists(file))
        stop(paste0("Error in import.CLUSTER: ",file," does not exist"))
        
    data <- utils::read.table(file,header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
    
    if(!is.null(exclude)){
        data <- exclude.markers(data,exclude)
    }
    markers        <- colnames(data)
    
    means          <- matrix(0,nrow=nrow(data),ncol=ncol(data))
    sd             <- matrix(0,nrow=nrow(data),ncol=ncol(data))
    for(i in 1:nrow(data)){
        for(j in 1:ncol(data)){
            tab        <- unlist(strsplit(data[i,j],split=";",fixed=TRUE))
            means[i,j] <- as.numeric(tab[1])
            sd[i,j]    <- as.numeric(tab[2])
        }
    }
    
    profiles               <- as.character(rownames(data))
    profiles.sizes         <- as.integer(rep(NA,nrow(data)))
    markers.clustering     <- rep(FALSE,ncol(data))
    densities              <- matrix(list())
    graph                  <- igraph::graph.empty(0,directed=FALSE)
    graph.layout           <- igraph::layout.auto(graph)
    colnames(graph.layout) <- NULL
    
    cluster <- CLUSTER(name = basename(file),
        profiles            = profiles,
        profiles.nb         = nrow(data),
        profiles.sizes      = profiles.sizes,
        markers             = markers,
        markers.nb          = length(markers),
        markers.clustering  = markers.clustering,
        means               = means,
        sd                  = sd,
        densities           = densities,
        graph               = graph,
        graph.layout        = graph.layout)
    return(cluster)
}


#' @title Importation of cell cluster profiles from SPADE results
#'
#' @description Imports one or several cell cluster profiles identified by the SPADE algorithm into a CLUSTER object.
#'
#' @details SPADE is a popular visualization and analysis algorithm that identifies clusters of cells having similar expression profiles for selected markers using an agglomerative hierarchical clustering-based algorithm combined with a density-based down-sampling procedure (PMID:21964415). Given a set of FCS files (usually one file per sample), SPADE identifies cell clusters based on the whole dataset and provides then for each sample the amount of cells present within each cluster.
#' 
#' The 'rescale' parameter can be used to rescale the marker expression intensities from 0 to 1 (with respect to the distribution proportion). The rescaling can be performed based on minimal and maximal expression values or based on specified quantiles using the 'rescale.quantiles' parameter. This strategy is especialy usefull when comparing cell or cell cluster profiles obtained from different experimental/staining conditions. 
#' 
#' @importFrom igraph read.graph
#'
#' @param path a character indicating the location to a zip or a folder containing the SPADE results
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#' @param trans a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation
#' @param trans.para a named list containing parameters for the transformation. Please refer to the details section for more details
#' @param trans.exclude a character vector containing the marker names for which no transformation must be applied on
#' @param bin.width a numeric value indicating the width of the bins for the marker expression densities computations
#' @param extract.folder a folder path for extracting the SPADE zip archive (temporary folder by default)
#' @param extract.folder.del a logical value indicating if the extracted SPADE results should be removed after the extraction
#' @param zip a logical value that specify if the path specifies a zip file
#' @param rescale a logical specifying if marker expression intensities must be rescale between 0 and 1
#' @param rescale.quantiles a numeric vector of two values specifying the quantiles of the marker expression intensities used to rescale
#'
#' @return a S4 object of class CLUSTER
#'
#' @export
import.SPADE <- function(path,
                            exclude            = NULL,
                            trans              = "arcsinh",
                            trans.para         = switch(trans,
                                                    "arcsinh" = list(arcsinh.scale=5),
                                                    "log"     = list(log.shift="auto",log.base=10),
                                                    "none"    = NULL),
                            trans.exclude      = NULL,
                            bin.width          = 0.05,
                            extract.folder     = NULL,
                            extract.folder.del = FALSE,
                            zip                = FALSE,
							rescale            = FALSE,
							rescale.quantiles  = c(0,1)){
    
    if(is.na(file.info(path)$isdir))
        stop(paste0("Error in import.SPADE: ",path," does not exist"))
    
    if(zip){
        if(is.null(extract.folder))
            extract.folder <- paste0(tempdir(),"/extract_spade-",basename(tempfile()))
        suppressWarnings(dir.create(extract.folder))
        message(paste0("Unzipping archive into: ",extract.folder,"..."),appendLF=FALSE)
        utils::unzip(path,exdir = extract.folder)
        message("done")
        path <- extract.folder
    }
    
    fcs.files <- list.files(path,full.names = TRUE,pattern=".+\\.cluster\\.fcs$")
    dictionary <- extract.dictionary(fcs.files[1])
    data <- import.FCS(fcs.files,
        trans             = trans,
        trans.para        = trans.para,
        trans.exclude     = unique(c(trans.exclude,"cluster")),
        exclude           = exclude,
		rescale           = rescale,
		rescale.quantiles = rescale.quantiles)
    
    cluster                    <- as.CLUSTER(object=data,cluster="cluster",bin.width=bin.width)
    
    cluster@name               <- basename(path)
    markers.clustering         <- unlist(strsplit(gsub("\\\"","",readLines(paste(path,"clusters.table",sep="/"),n=1))," "))
    markers.clustering         <- rename.markers(markers.clustering,dictionary)
    cluster@markers.clustering <- cluster@markers %in% markers.clustering
    cluster@graph              <- igraph::read.graph(paste(path,"mst.gml",sep="/"),format="gml")
    cluster@graph.layout       <- as.matrix(utils::read.table(paste0(path,"/layout.table")))
    colnames(cluster@graph.layout) <- c("x","y")
    rownames(cluster@graph.layout) <- 1:nrow(cluster@graph.layout)
    cluster@overview.function      <- "ggplot.SPADEtree"
    
    if(!is.null(extract.folder) && extract.folder.del==TRUE){
        unlink(path,recursive = TRUE)
    }
    return(cluster)
}


#' @title Importation of cell cluster profiles from a Citrus result
#'
#' @description Imports one or several cell cluster profiles identified by the Citrus algorithm into a CLUSTER object.
#'
#' @details Citrus is an algorithm that clusters cells using a hierarchical clustering procedure (similarly to SPADE) and then identifies the cell clusters that are significantly associated with different biological condition phenotypes (PMID:24979804). 
#'
#' @param file a character indicating the location of the citrusClustering.Rdata file
#' @param dictionary a two-column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)
#' @param exclude a vector containing the marker names to be excluded in the import procedure
#' @param bin.width a numeric value indicating the width of the bins for the marker expression densities computations
#' @param minimumClusterSizePercent a numeric value indicating the minimal ratio of cells per cluster to import
#' @param cluster.selection a character vector containing the names of the clusters to import
#'
#' @return a S4 object of class CLUSTER
#'
#' @export
import.CITRUS <- function(file,
                          dictionary                = NULL,
                          exclude                   = NULL,
                          bin.width                 = 0.05,
                          minimumClusterSizePercent = 0.05,
                          cluster.selection         = NULL){
    
    message(paste0("Importing ",file))
    if(is.na(file.exists(file)))
        stop(paste0("Error in import.CITRUS: ",file," does not exist"))
        
    env <- new.env()
    load(file,env)
    citrus.combinedFCSSet <- env[["citrus.combinedFCSSet"]]
    citrus.foldClustering <- env[["citrus.foldClustering"]]
    citrus.clustering     <- citrus.foldClustering$allClustering
    
    # from citrus.selectClusters.minimumClusterSize
    citrus.allCluster   <- citrus.foldClustering$allClustering
    cluster.sizes       <- sapply(citrus.clustering$clusterMembership,length)
    size.min            <- (length(cluster.sizes)+1) * minimumClusterSizePercent
    largeEnoughClusters <- which(cluster.sizes >= size.min)
    # from citrus.selectClusters.minimumClusterSize
    
    if(!is.null(cluster.selection)){
        largeEnoughClusters <- largeEnoughClusters[largeEnoughClusters%in%cluster.selection]
        if(length(largeEnoughClusters) == 0)
            stop("Selected clusters do not have enough cell")
    }
    data <- c()
    for(clusterId in largeEnoughClusters){
        
        # from citrus.exportCluster - begin
        clusterData = citrus.combinedFCSSet$data[citrus.clustering$clusterMembership[[clusterId]],]
        if (!is.null(citrus.combinedFCSSet$scaleColumns)) {
            clusterData[, citrus.combinedFCSSet$scaleColumns] = t((t(clusterData[,
            citrus.combinedFCSSet$scaleColumns]) * citrus.combinedFCSSet$scaleColumns.SD) +
            citrus.combinedFCSSet$scaleColumns.mean)
        }
        if (!is.null(citrus.combinedFCSSet$transformColumns)) {
            clusterData[, citrus.combinedFCSSet$transformColumns] = sinh(clusterData[,
            citrus.combinedFCSSet$transformColumns]) * citrus.combinedFCSSet$transformCofactor
        }
        # from citrus.exportCluster - end
        
        data <- rbind(data,clusterData)
    }
    data <- citrus.combinedFCSSet$data
    
    
    if(!is.null(dictionary)){
        colnames(data) <- rename.markers(colnames(data),dictionary)
    }
    if(!is.null(exclude)){
        data <- exclude.markers(data,exclude)
    }
    markers        <- colnames(data)
    
    cluster <- as.CLUSTER(data,cluster = "fileEventNumber",bin.width = bin.width)
    cluster@name <- basename(file)
    
    
    return(cluster)
}


#' @title Importation of range gate profiles from a tab separated file
#'
#' @description Imports one or several range gate profiles from a tab separated file into a GATE object.
#'
#' @details Tab separated file to import must contain for each gate profile the ranges of the expression markers and must be formatted as the following:\cr
#' * each row must represent a gate profile;\cr
#' * each column must represent a marker;\cr
#' * each cell in the table must contain the marker expression lower and upper bounds for a given gate profile separated by a semicolon.\cr
#' The first column must contain the gate names and the first row must contain the marker names.
#'
#' @param file a character indicating the location of a tab separated file to import
#' @param exclude a character vector containing the marker names to be excluded in the import procedure
#'
#' @return a S4 object of class GATE
#'
#' @export
import.GATE <- function(file,
                        exclude    = NULL){
    
    message(paste0("Importing ",file))
    if(!file.exists(file))
        stop(paste0("Error in import.GATE: ",file," does not exist"))
    
    data <- utils::read.table(file,header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
    
    if(!is.null(exclude)){
        data <- exclude.markers(data,exclude)
    }
    markers        <- colnames(data)
    
    ranges_inf     <- matrix(0,nrow=nrow(data),ncol=ncol(data))
    ranges_sup     <- matrix(0,nrow=nrow(data),ncol=ncol(data))
    for(i in 1:nrow(data)){
        for(j in 1:ncol(data)){
            tab <- unlist(strsplit(data[i,j],split=":",fixed=TRUE))
            ranges_inf[i,j] <- as.numeric(tab[1])
            ranges_sup[i,j] <- as.numeric(tab[2])
        }
    }
    
    ranges           <- cbind(ranges_inf,ranges_sup)
    dim(ranges)      <- c(nrow(data),ncol(data),2)
    dimnames(ranges) <- list(NULL,NULL,c("inf","sup"))
    
    gate <- GATE(name  = basename(file),
        markers        = markers,
        markers.nb     = ncol(data),
        profiles       = as.character(rownames(data)),
        profiles.nb    = nrow(data),
        ranges         = ranges)
    return(gate)
}


#' @title Importation of  range gate profiles from a Gating-ML file
#'
#' @description Imports one or several range gate profiles from a Gating-ML file into a GATE object.
#'
#' @details Gating-ML is a standard file format for gate definitions developed to facilitate the interchange between different analysis software. Gating-ML rectangular gates, that are defined by ranges of expression markers, can be imported in CytoCompare using this function.\cr
#'
#' @import XML
#'
#' @param file a character indicating the location of a Gating-ML xml file to import
#' @param filterId a character vector indicating the identifiers of the gates to import
#'
#' @return a S4 object of class GATE
#'
#' @export
import.GATINGML <- function(file,filterId=NULL){
    
    message(paste0("Importing ",file))
    if(!file.exists(file))
        stop(paste0("Error in import.GATE: ",file," does not exist"))
    
    # the flowUtils::read.gatingML function seems to be not correctly exported in the flowUtils package and cannot be used
    # alternative is to parse the Gating-ML file via the XML library
    xml     <- XML::xmlParse(file) 
    xmltop  <- XML::xmlRoot(xml)

    gate    <- list()
    count   <- 1
    for(i in 1:length(XML::xmlChildren(xmltop))){
        element <- xmltop[[i]]
        
        if(XML::xmlName(element)=="RectangleGate"){
            
            id <- XML::xmlAttrs(element)[["id"]]
            
            if(is.null(filterId) || (!is.null(filterId) && id %in% filterId)){
                marker1          <- element[[2]]
                marker1_min      <- as.numeric(XML::xmlAttrs(marker1)[["min"]])
                marker1_max      <- as.numeric(XML::xmlAttrs(marker1)[["max"]])
                marker1_name     <- XML::xmlAttrs(marker1[[1]])[["name"]]
                
                marker2          <- element[[3]]
                marker2_min      <- as.numeric(XML::xmlAttrs(marker2)[["min"]])
                marker2_max      <- as.numeric(XML::xmlAttrs(marker2)[["max"]])
                marker2_name     <- XML::xmlAttrs(marker2[[1]])[["name"]]
                
                ranges_inf       <- c(marker1_min,marker2_min)
                ranges_sup       <- c(marker1_max,marker2_max)
                markers          <- c(marker1_name,marker2_name)
                    
                ranges           <- c(ranges_inf,ranges_sup)
                dim(ranges)      <- c(1,length(markers),2)
                dimnames(ranges) <- list(NULL,NULL,c("inf","sup"))
                
                gate[[count]] <- GATE(profiles = id,
                    profiles.nb                = as.integer(1),
                    markers                    = markers,
                    markers.nb                 = length(markers),
                    ranges                     = ranges)
                count <- count+1
            }
        }
    }
    
    gate      <- do.call("c",gate)
    gate@name <- basename(file)
    return(gate)
}

# title Internal - Extract a dictionary from a FCS file
#
# @description This function is used internally to extract the correspondence between the original marker names (first column) and the true marker names (second column)
#
# @param fcs.file a character indicating the location of the fcs file containing the correspondences
#
# @return a two-column data.frame providing the correspondence between the original marker names (first column) and the true marker names (second column)
extract.dictionary <- function(fcs.file){
    flowframe       <- flowCore::read.FCS(fcs.file[1],trans = FALSE)
    dictionary      <- flowframe@parameters@data[, c(1, 2)]
    dictionary[, 1] <- make.names(dictionary[, 1])
    return(dictionary)
}