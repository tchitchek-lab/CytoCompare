# @title CELL class definition
#
# @description CELL is a S4 object containing one or several cell profiles.
#
# @details This object mainly stores for each cell profile, the intensities of each marker.
#
# The slot 'trans.para' is a named list contains different parameters depending of the transformation applied on the marker expression intensities. The scale (cofactor) of the arcsinh transformation function is parametrized using the 'arcsinh.scale' value. The shift of the log transformation function is parametrized using the 'log.shift' value and the base of the log transformation function is parametrized using the 'log.base' value. If no transformation function have been applied, the 'trans.para' slot is set to NULL.
#
# @slot name a character indicating the internal name of the CELL object
# @slot profiles a character vector containing the names of the cell profiles 
# @slot profiles.nb an integer value indicating the number of cell profiles
# @slot markers a character vector containing the marker names
# @slot markers.nb an integer value indicating the number of markers
# @slot intensities a numeric matrix containing the intensities of each marker for each cell profile
# @slot trans a character specifying the name of a transformation function applied on the marker expression intensities. Possible values are "arcsinh" for arc sin hyperbolic transformation, "log" for logarithmic transformation, or "none" for no transformation
# @slot trans.para a named list containing parameters of the transformation. Please refer to the details section for more details
# @slot trans.exclude a character vector containing the marker names for which no transformation has been applied on
# @slot overview.function a character specifying the name of a function to call when plotting the CELL object overview (please refer to the documentation of the 'plot()' function)
# @slot layout a numeric matrix that can be used to store the positions of cells in a 2-dimensional space (e.g. tSNE1 and tSNE2 dimensions provided by viSNE)
#
# @import methods
#
# @name CELL-class
# @rdname CELL-class
# @exportClass CELL
CELL <- setClass("CELL",
    slots=c(name      = "character",
		profiles          = "character",
		profiles.nb       = "integer",
		markers           = "character",
		markers.nb        = "integer",
		intensities       = "matrix",
		trans             = "character",
		trans.para        = "ANY",
		trans.exclude     = "ANY",
		overview.function = "character",
		layout            = "matrix"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@profiles.nb!=length(object@profiles))
            stop("Error in profiles.nb slot: profiles.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        if(!is.element(object@trans,c("arcsinh","log","none")))
            stop("Error in trans slot: trans do not contain allowed value (allowed value are \"arcsinh\", \"log\", \"none\")")
        if(!is.null(object@trans.para) && class(object@trans.para)[1]!="list")
            stop("Error in trans.para slot: trans.para must be a of type list or NULL")
        if(!is.null(object@trans.exclude) && class(object@trans.exclude)[1]!="character")
            stop("Error in trans.exclude slot: trans.exclude must be a of type character or NULL")
        return(TRUE)
    }
)
setMethod("initialize",c("CELL"),
    function(.Object,
            name              = "",
            profiles          = "",
            profiles.nb       = 0,
            markers           = "",
            markers.nb        = 0,
            intensities       = as.matrix(0),
            overview.function = "",
            layout            = as.matrix(0)){
        if(profiles.nb==0)
            stop("Error can not create a CELL object with no profile")
        .Object@name              = name     
        .Object@profiles          = profiles
        .Object@profiles.nb       = profiles.nb
        .Object@markers           = markers
        .Object@markers.nb        = markers.nb
        .Object@trans             = "none"
		.Object@trans.para        = NULL
		.Object@trans.exclude     = NULL
		.Object@intensities       = intensities
        .Object@overview.function = overview.function
        .Object@layout            = layout
        validObject(.Object)
        return(.Object)
    }
)


#' @title CLUSTER class definition
#'
#' @description CLUSTER is a S4 object containing one or several cell cluster profiles.
#'
#' @details This object mainly stores for each cell cluster profile, the means, the standard deviations and the densities of each marker.
#'
#' @slot name a character indicating the internal name of the CLUSTER object
#' @slot profiles a character vector containing the names of the cell cluster profiles
#' @slot profiles.nb an integer value indicating the number of cell cluster profiles
#' @slot profiles.sizes an integer vector indicating the number of cells associated to each cluster profile
#' @slot markers a character vector containing the marker names
#' @slot markers.nb an integer value indicating the number of markers
#' @slot markers.clustering a logical vector specifying the makers used as clustering markers
#' @slot means a numeric matrix containing the means of each maker for each cluster profile
#' @slot sd a numeric matrix containing the standard deviations of each maker for each cluster profile
#' @slot densities a matrix of DENSITY objects containing the densities of each marker for each cluster profile
#' @slot overview.function a character specifying the name of a function to call when plotting the CLUSTER object overview (please refer to the documentation of the 'plot()' function)
#' @slot graph an object that can be used to store a visual representation of the cell clusters (e.g. a SPADE tree)
#' @slot graph.layout a numeric matrix that can be used to store the positions of cell clusters in a 2-dimensional space (e.g. a SPADE tree layout)
#' 
#' @import methods
#'
#' @name CLUSTER-class
#' @rdname CLUSTER-class
#' @exportClass CLUSTER
#' @export CLUSTER
CLUSTER <- setClass("CLUSTER",
    slots=c(name           = "character",
        profiles           = "character",
        profiles.nb        = "integer",
        profiles.sizes     = "integer",
        markers            = "character",
        markers.nb         = "integer",
        markers.clustering = "logical",
        means              = "matrix",
        sd                 = "matrix",
        overview.function  = "character",
        graph              = "ANY",
        graph.layout       = "matrix",
        densities          = "matrix"),
    validity=function(object){
        if(length(object@profiles)!=length(unique(object@profiles)))
            stop("Error in profiles slot: profile names are not unique")
        if(length(object@markers)!=length(unique(object@markers)))
            stop("Error in markers slot: marker names are not unique")
        if(object@profiles.nb!=length(object@profiles))
            stop("Error in profiles.nb slot: profiles.nb do not correspond to the number of profile names")
        if(object@markers.nb!=length(object@markers))
            stop("Error in markers.nb slot: markers.nb do not correspond to the number of marker names")
        if(object@profiles.nb==0)
            stop("Error can not create a CLUSTER object with no profile")
        return(TRUE)
    }
)
setMethod("initialize",c("CLUSTER"),
    function(.Object,
            name               = "",
            profiles           = "",
            profiles.nb        = 0,
            profiles.sizes     = 0,
            markers            = "",
            markers.nb         = 0,
            markers.clustering = FALSE,
            means              = as.matrix(0),
            sd                 = as.matrix(0),
            overview.function  = "",
            graph              = NULL,
            graph.layout       = as.matrix(0),
            densities          = as.matrix(0)){
        if(profiles.nb==0)
            stop("Error can not create a CLUSTER object with no profile")
        .Object@name               <- name     
        .Object@profiles           <- profiles
        .Object@profiles.nb        <- profiles.nb
        .Object@profiles.sizes     <- profiles.sizes
        .Object@markers            <- markers
        .Object@markers.nb         <- markers.nb
        .Object@markers.clustering <- markers.clustering
        .Object@means              <- means
        .Object@sd                 <- sd
        .Object@overview.function  <- overview.function
        .Object@graph              <- graph
        .Object@graph.layout       <- graph.layout
        .Object@densities          <- densities
        validObject(.Object)
        return(.Object)
    }
)


#' @title DENSITY class definition
#'
#' @description DENSITY is a S4 object used to stores a marker expression density.
#'
#' @details This object mainly stores for each marker: the bin characteristics, the negative and positive marker densities values, and the number of cells used in the density estimation. Densities are stored using two numeric vectors: 'values.neg' for the negative densities and 'values.pos' for the positive densities. This strategy allows to compute and store densities without defining an absolute minimal value or an absolute maximal value.
#'
#' @slot name a character indicating the internal name of the CELL object
#' @slot bin.interval a numeric vector of two values specifying the density boundaries
#' @slot bin.nb a numeric vector of two values specifying the numbers of negative and positive bins
#' @slot values.pos a numeric vector containing the positive density bins
#' @slot values.neg a numeric vector containing the negative density bins
#' @slot point.nb a numeric value indicating the number of point used to compute the expression density
#' @slot bin.width a numeric value indicating the width of the bins used in the density estimation
#'
#' @import methods
#'
#' @name DENSITY-class
#' @rdname DENSITY-class
#' @exportClass DENSITY
#' @export DENSITY
DENSITY <- setClass("DENSITY",
    slots=c(name         = "character",
            bin.interval = "numeric",
            bin.nb       = "numeric",
            values.pos   = "numeric",
            values.neg   = "numeric",
            point.nb     = "numeric",
            bin.width    = "numeric")
)
setMethod("initialize",c("DENSITY"),
    function(.Object,bin.width=0.05,bin.interval,point.nb,values.neg,values.pos,values=NULL,name=""){
        if(is.null(values)){
            .Object@bin.interval <- bin.interval
            .Object@bin.nb       <- c(length(values.neg),length(values.pos))
            .Object@point.nb     <- point.nb
            .Object@values.pos   <- values.pos
            .Object@values.neg   <- values.neg
        }else{
            .Object@point.nb     <- length(values)
            density              <- compute.density(values,bin.width)
            .Object@bin.interval <- density$bin.interval
            .Object@bin.nb       <- density$bin.nb
            .Object@values.pos   <- density$values.pos
            .Object@values.neg   <- density$values.neg 
        }
        .Object@bin.width <- bin.width
        .Object@name      <- name
        validObject(.Object)
        return(.Object)
    }
)


#' @title MWEIGHTS class definition
#'
#' @description MWEIGHTS is a S4 object containing the marker weights to use in the comparison computations. 
#'
#' @details This object mainly stores for each marker: the markers names and marker weights. 
#'
#' @slot markers a character vector containing the marker names
#' @slot weights a numeric vector containing the marker weights
#'
#' @import methods
#'
#' @name MWEIGHTS-class
#' @rdname MWEIGHTS-class
#' @exportClass MWEIGHTS
#' @export MWEIGHTS
MWEIGHTS <- setClass("MWEIGHTS",
    slots = c(markers = "character",
              weights = "numeric")
)


#' @title RES class definition
#'
#' @description RES is a S4 object containing one or several comparison results.
#'
#' @details This object mainly stores for each comparison result: the associated distances (with the associated distance threshold) or inclusion p-value, and the marker successes. In the case of similarity comparisons, the aggregated distance between the two profiles and the marker distances are also stored in this object.
#'
#' @slot comparisons a data.frame containing for each comparison: the profile names, the type of the comparison (similarity or inclusion), the aggregated distance (or NA in case of inclusion), the distance threshold used (\code{$D_th$}) and the associated p-value
#' @slot comparisons.nb is an integer indicating the number of comparisons  
#' @slot markers a character vector containing the marker names used in the comparisons
#' @slot marker.weights a numeric vector containing the weigths associated to each marker involded in the comparisons
#' @slot marker.distances a data.frame containing the marker distances for each comparison (or NA in case of inclusion assessments)
#' @slot marker.successes a data.frame containing the marker successes for each comparison
#' 
#' @import methods
#'
#' @name RES-class
#' @rdname RES-class
#' @exportClass RES
#' @export RES
RES <- setClass("RES",
    slots=c(comparisons   = "data.frame",
        comparisons.nb    = "numeric",
        markers           = "character",
		marker.weights    = "numeric",
        marker.distances  = "data.frame",
        marker.successes  = "data.frame"),
    validity=function(object){ 
        if(object@comparisons.nb!=nrow(object@comparisons))
            stop("Error in comparisons.nb slot: comparisons.nb do not correspond to the number of comparisons")  
        if(nrow(object@marker.distances)!=nrow(object@comparisons))
            stop("Error in marker.distances slot: marker.distances do not correspond to the number of comparisons")
        if(nrow(object@marker.successes)!=nrow(object@comparisons))
            stop("Error in marker.successes slot: marker.successes do not correspond to the number of comparisons")
        if(length(object@markers)!=length(object@marker.weights))
            stop("Error in marker.weights slot: marker.weights do not correspond to the number of markers")
        return(TRUE)
    }
)
setMethod("initialize",c("RES"),
          function(.Object,
                   comparisons      = data.frame(),
                   comparisons.nb   = 0,
                   markers          = character(0),
                   marker.distances = data.frame(),
                   marker.successes = data.frame(),
                   marker.weights   = integer(0)){
            .Object@comparisons      <- comparisons     
            .Object@comparisons.nb   <- comparisons.nb
            .Object@markers          <- markers
            .Object@marker.distances <- marker.distances          
            .Object@marker.successes <- marker.successes
            .Object@marker.weights   <- marker.weights
            validObject(.Object)
            return(.Object)
          }
)


# title Internal - Common markers between two cytometry objects
#
# @description Identifies the common markers between two cytometry objects (CELL, CLUSTER or GATE objects) and store the results in a MWEIGHTS object.
#
# @details The weight of each marker is set to '1' but can be changed afterwards using the 'set()' function.
#
# @param object1 a CELL, CLUSTER or GATE object
# @param object2 another CELL, CLUSTER or GATE object
#
# @return a S4 object of class MWEIGHTS.
inner.intersect <- function(object1,object2){
    cmarkers <- intersect(object1@markers,object2@markers)
    weights  <- rep(1,length(cmarkers))
    mweights <- MWEIGHTS(markers=cmarkers,weights=weights)
    return(mweights)
}


#' @title Identification of common markers between two cytometry objects
#'
#' @description Identifies the common markers between two cytometry objects (CELL, CLUSTER or GATE objects) and store the results in a MWEIGHTS object.
#'
#' @details The weight of each marker is set to 1 but can be changed afterwards using the function set.
#'
#' @param x a CELL, CLUSTER or GATE object
#' @param y a CELL, CLUSTER or GATE object
#'
#' @return a S4 object of class MWEIGHTS
#' @name intersect
#' @rdname intersect-methods
NULL


#' @rdname intersect-methods
#' @export
setMethod("intersect",c("CLUSTER","CLUSTER"),
    function(x,y){
        inner.intersect(x,y)
    }
)


#' @title Extraction of subsets of data from CytoCompare objects
#'
#' @description Extracts subsets of CELL, CLUSTER, GATE, MWEIGHTS or RES object.
#'
#' @details For cytometry objects (CELL, CLUSTER, or GATE objects), the parameter 'i' represents a vector of profiles to extract and the parameter 'j' represents a vector of markers to extract.
#'
#' For MWEIGHTS objects, the parameter 'i' represents a vector of markers to extract.
#'
#' For RES objects, the parameter 'i' represents a vector of comparisons to extract.
#'
#' @param x a CELL, CLUSTER, GATE, MWEIGHTS or RES object
#' @param i a numeric, logical or character vector
#' @param j a numeric, logical or character vector
#'
#' @return a S4 object of class CELL, CLUSTER, GATE, MWEIGHTS or RES
#'
#' @name extract
#' @rdname extract-methods
NULL


#' @rdname extract-methods
#' @export
setMethod("[",c("CLUSTER","ANY","ANY"),
    function(x,i,j){
    
        if(!missing(i) && length(i)>x@profiles.nb)
            stop("Too many cluster profiles to extract")
        if(!missing(j) && length(j)>length(x@markers))
            stop("Too many markers to extract")
        if(!missing(i) && !(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(!missing(j) && !(is.logical(j) || is.numeric(j) || is.character(j)))
            stop(paste0("Wrong types in j: ",typeof(j)))
            
        if(missing(i)){
            k <- 1:x@profiles.nb
        }else{
            if(is.character(i)){
                k <- which(x@profiles %in% i)
            }else{
                k <- i
            }
        }
        
        if(missing(j)){
            l <- 1:length(x@markers)
        }else{
            if(is.character(j)){
                l <- match(j,x@markers)
            }else{
                l <- j
            }
        }
        
        if(missing(i)){
            overview.function <- x@overview.function    
            graph             <- x@graph
            graph.layout      <- x@graph.layout
        }else{
            overview.function <- ""
            graph             <- igraph::graph.empty(0,directed=FALSE)
            graph.layout      <- igraph::layout.auto(graph)
        }
            
        cluster <- CLUSTER(name = x@name,
            profiles            = x@profiles[k],
            profiles.nb         = length(x@profiles[k]),
            profiles.sizes      = as.integer(x@profiles.sizes[k]),
            markers             = x@markers[l],
            markers.nb          = length(x@markers[l]),
            markers.clustering  = x@markers.clustering[l],
            means               = x@means[k,l,drop=FALSE],
            sd                  = x@sd[k,l,drop=FALSE],
            densities           = x@densities[k,l,drop=FALSE],
            graph               = graph,
            graph.layout        = graph.layout,
            overview.function   = overview.function)
        return(cluster)
    }
)


#' @rdname extract-methods
#' @export
setMethod("[",c("MWEIGHTS","ANY"),
    function(x,i){
    
        if(!(is.character(i) || is.numeric(i) || is.logical(i)))
            stop("Wrong types in parameter i")
    
        if(is.character(i)){
            k <- match(i,x@markers)
        }else{
            k <- i
        }
        
        return(MWEIGHTS(markers=x@markers[k],weights=x@weights[k]))    
    }
)

#' @rdname extract-methods
#' @export
setMethod("[",c("RES","ANY"),
    function(x,i){
        
        if(!(is.logical(i) || is.numeric(i) || is.character(i)))
            stop(paste0("Wrong types in i: ",typeof(i)))
        if(is.logical(i) && length(i)!=nrow(x@comparisons))
            stop("logical vector i must have the same length as the number of comparisons in the RES object")
            
        if(is.character(i)){
            idx <- rownames(x@comparisons) %in% i
        }else{
            idx <- i
        }
        
        comparisons <- x@comparisons[idx,]
        
        res <- RES(comparisons = comparisons,
            comparisons.nb     = nrow(comparisons),
            markers            = x@markers,
            marker.distances   = x@marker.distances[idx,],
            marker.successes   = x@marker.successes[idx,],
            marker.weights     = x@marker.weights)
        return(res)
    }
)


#' @title Change marker weights in a MWEIGHTS object
#'
#' @description Sets the weights of the markers in a MWEIGHTS object.
#'
#' @details Marker weights can be set based on their indexes or names in the MWEIGHTS object.
#'
#' @param x a MWEIGHTS object
#' @param i a numeric, logical or character vector
#' @param value a numeric vector containing the new marker weight values
#' @return a S4 object of class MWEIGHTS
#'
#' @name set
#' @rdname set-methods
NULL

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","numeric","missing","numeric"),
    function(x,i,value){
    
        if(min(i)<1 || max(i)>length(x@weights))
            stop("Indices out of range")
        if(length(i)>length(x@weights))
            stop("Too many markers to modify")
            
        x@weights[i] <- value
        return(x)
    }
)

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","character","missing","numeric"),
    function(x,i,value){
    
        if(sum(!i %in% x@markers)>0)
            stop("Unknown marker in i")
        if(length(i)>length(x@weights))
            stop("Too many markers to modify")
            
        idx            <- match(i,x@markers)
        x@weights[idx] <- value
        return(x)
    }
)

#' @rdname set-methods
#' @export
setMethod("[<-",c("MWEIGHTS","logical","missing","numeric"),
    function(x,i,value){
    
        if(length(i)!=length(x@weights))
            stop("logical vector i must have the same length as the number of weights in the MWEIGHTS object")
        
        x@weights[i] <- value
        return(x)
    }
)


#' Combination of CytoCompare objects
#'
#' @description Combines two or several CELL, CLUSTER, GATE or RES objects.
#'
#' @details All the different objects to combine must be of the same type.
#'
#' This function is especially useful when combining comparison results from different RES objects into a single RES object (the RES objects must share the same markers and marker weights). RES objects can be combined to an empty RES object (i.e. RES()).
#' 
#' This function is also especially useful when combining cell profiles obtained from different FCS files into one single CELL object. 
#'
#' @param x a first CELL, CLUSTER, GATE or RES object
#' @param ... further objects of the same class as x to be combined
#' @param recursive a logical value indicating if the function recursively descends through lists combining all their elements into a vector. Not implemented and should be set to FALSE
#'
#' @return a S4 object of class CELL, CLUSTER, GATE or RES
#'
#' @name c
#' @rdname c-methods
NULL

#' @rdname c-methods
#' @export
setMethod("c",c("CELL"),
    function(x,...){
    
        other.CELL  <- list(x,...)
        name        <- c()
        markers     <- x@markers
        profiles.nb <- 0
        
        i <- 1
        for(cell in other.CELL){
            if(class(cell) != "CELL")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(cell)))
            name        <- c(name,cell@name)
            markers     <- union(markers,cell@markers)
            profiles.nb <- profiles.nb+cell@profiles.nb
            i <- i+1
        }
        name        <- paste0(name,collapse=";")
        profiles    <- as.character(1:profiles.nb)
        
        intensities <- matrix(NA,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        
        i <- 1
        for(cell in other.CELL){
            nb                                   <- cell@profiles.nb
            intensities[i:(i+nb-1),cell@markers] <- cell@intensities
            i <- i+nb
        }
        dimnames(intensities) <- NULL
        
        cell <- CELL(name = name,
            profiles      = profiles,
            profiles.nb   = length(profiles),
            markers       = markers,
            markers.nb    = length(markers),
            intensities   = intensities)
        return(cell)
    }
)

#' @rdname c-methods
#' @export
setMethod("c",c("CLUSTER"),       
    function(x,...){
    
        other.CLUSTER  <- list(x,...)
        name           <- c()
        markers        <- c()
        profiles       <- c()
        profiles.sizes <- c()
        profiles.nb    <- 0
        
        i <- 1 
        for(cluster in other.CLUSTER){
            if(class(cluster) != "CLUSTER")
                stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(cluster)))
            name           <- c(name,cluster@name)
            markers        <- union(markers,cluster@markers)
            profiles.sizes <- c(profiles.sizes,cluster@profiles.sizes)
            profiles.nb    <- profiles.nb+cluster@profiles.nb
            profiles       <- c(profiles,cluster@profiles)
            i <- i+1
        }
        
        if(length(profiles)!=length(unique(profiles))){
            stop("Error: profile names are not unique")
        }
        
        name <- paste0(name,collapse=";")
        
        means          <- matrix(0,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        sd             <- matrix(0,ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        densities      <- matrix(list(),ncol=length(markers),nrow=profiles.nb,dimnames=list(profiles,markers))
        
        i <- 1
        for(cluster in other.CLUSTER){
            nb                                    <- cluster@profiles.nb
            means[i:(i+nb-1),cluster@markers]     <- cluster@means
            sd[i:(i+nb-1),cluster@markers]        <- cluster@sd
            densities[i:(i+nb-1),cluster@markers] <- cluster@densities
            i <- i+nb
        }
        dimnames(means)     <- NULL
        dimnames(sd)        <- NULL
        dimnames(densities) <- NULL
        
        markers.clustering  <- rep(FALSE,length(markers))
        graph               <- igraph::graph.empty(0,directed=FALSE)
        graph.layout        <- igraph::layout.auto(graph)
        
        cluster <- CLUSTER(name = name,
            profiles            = profiles,
            profiles.nb         = length(profiles),
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
)


#' @rdname c-methods
#' @export
setMethod("c",c("RES"),   
function(x,...){

    other.RES        <- list(x,...)
    comparisons      <- NULL
    distances        <- NULL
    successes        <- NULL
    markers          <- NULL
    marker.weights   <- NULL
    i <- 1
    for(res in other.RES){
        if(class(res) != "RES")
            stop(paste0("Cannot combine objects of different classes (element at position ",i," is of type ",class(res)))
        if(any(rownames(res@comparisons) %in% rownames(comparisons)))
            stop(paste0("Duplicated comparison results found"))

        if(res@comparisons.nb>0){
            if(i==1){
                markers        <- res@markers
                marker.weights <- res@marker.weights
            }
            if(!identical(res@markers,markers)){
                stop("Error: to combine RES object, all markers slot must be identical")
            }
            if(!identical(res@marker.weights,marker.weights)){
                stop("Error: to combine RES object, all marker.weights slot must be identical")
            }
            
            comparisons <- rbind(comparisons,res@comparisons)
            distances   <- rbind(distances,res@marker.distances)
            successes   <- rbind(successes,res@marker.successes)
            
            i <- i+1
        }
        
    }

    rownames(comparisons) <- paste(comparisons[,"profile1"],comparisons[,"profile2"],sep="/vs/")
    comparisons[,"profile1"]  <- as.character(comparisons[,"profile1"])
    comparisons[,"profile2"]  <- as.character(comparisons[,"profile2"])
    
    res <- RES(comparisons = comparisons,
        comparisons.nb     = nrow(comparisons),
        markers            = markers,
        marker.distances   = distances,
        marker.successes   = successes,
        marker.weights     = marker.weights)
    return(res)
    }
)


#' @title Coercion to a CELL object
#'
#' @description Coerces a numeric matrix into a CELL object.
#'
#' This function transforms a numeric matrix into one or several cell profiles.
#'
#' @details The matrix must have its column names corresponding to the cell markers.
#'
#' @param object a numeric matrix
#' @param name a character specifying the internal name of the CELL object to create
#'
#' @return a S4 object of class CELL
#'
#' @name as.CELL
#' @rdname as.CELL-methods
#'
#' @export
setGeneric("as.CELL", function(object,name="cell") { standardGeneric("as.CELL") })

#' @rdname as.CELL-methods
#' @export
setMethod("as.CELL",c("matrix"),
    function(object,name){  
        data           <- object
        dimnames(data) <- NULL
        cell <- CELL(name = name,
            profiles      = as.character(1:nrow(object)),
            profiles.nb   = nrow(object),
            markers       = colnames(object),
            markers.nb    = ncol(object),
            intensities   = data)
        return(cell)
    }
)

#' @title Coercion to a CLUSTER object
#'
#' @description Coerces a CELL object or a numeric matrix into a CLUSTER object.
#'
#' This function transforms the cell profiles from a CELL object or from a numeric matrix into one or several cell cluster profiles by computing the means, the standard deviations, and the densities of each marker.
#'
#' @details The 'cluster' parameter is especially useful when importing FCS files containing the SPADE clustering results (where an additional channel is used to indicate the associations between cells and cell clusters) or FCS files from any other automatic gating algorithm.
#'
#' In the context of a numeric matrix coercion, the matrix must have its column names corresponding to the cell markers.
#'
#' @importFrom igraph graph.empty layout.auto
#'
#' @param object a CELL object or a numeric matrix
#' @param name a character specifying the internal name of the CLUSTER object to create
#' @param cluster a character indicating a channel name that can be used to gather the cell profiles into several cluster profiles. If a channel named is specified then the created CLUSTER object will contain as many profiles as different values present in this channel. If this parameter is NULL then a CLUSTER object with only one profile will be created
#' @param bin.width a numeric value indicating the width of the bins in the density estimation computations (default=0.05)
#'
#' @return a S4 object of class CLUSTER
#'
#' @name as.CLUSTER
#' @rdname as.CLUSTER-methods
#'
#' @export
setGeneric("as.CLUSTER",function(object,name=object@name,cluster=NULL,bin.width=0.05){standardGeneric("as.CLUSTER")})

#' @rdname as.CLUSTER-methods
#' @export
setMethod("as.CLUSTER",c("CELL"),
    function(object,name=object@name,cluster=NULL,bin.width=0.05){
        colnames(object@intensities) <- object@markers
        cluster                      <- as.CLUSTER(object=object@intensities,name=object@name,cluster=cluster,bin.width=bin.width)
        return(cluster)
    }
)

#' @rdname as.CLUSTER-methods
#' @export
setMethod("as.CLUSTER",c("matrix"),
    function(object,name="cell_cluster",cluster=NULL,bin.width=0.05){
        markers <- colnames(object)
        if(is.null(cluster)){
            profiles       <- "1"
            profiles.nb    <- as.integer(1)
            profiles.sizes <- integer(profiles.nb)
            message("Computation of marker means, standard deviations and densities...")
            means          <- matrix(apply(object,2,mean),nrow=profiles.nb,ncol=length(markers))
            sd             <- matrix(apply(object,2,function(x){sqrt(stats::var(x))}),nrow=profiles.nb,ncol=length(markers))
            densities      <- matrix(list(),nrow=1,ncol=length(markers))
            for(i in 1:length(markers)){
                densities[[i]] <- DENSITY(bin.width=bin.width,values=object[,i],name=markers[i])
                message(paste0(round(i/length(markers)*100),"%"))
            }
            message("done")
            profiles.sizes <- nrow(object)
        }else{
            val.cluster    <- sort(unique(object[,cluster]))
            profiles       <- as.character(val.cluster)
            profiles.nb    <- length(val.cluster)
            profiles.sizes <- integer(profiles.nb)
            markers        <- markers[markers!=cluster]
            message("Computation of marker means, standard deviations and densities...")
            means          <- matrix(0,nrow=profiles.nb,ncol=length(markers))
            sd             <- matrix(0,nrow=profiles.nb,ncol=length(markers))
            densities      <- matrix(list(),nrow=profiles.nb,ncol=length(markers))
            for(i in 1:profiles.nb){
				current_profile   <- object[object[,cluster]==val.cluster[i],markers,drop=FALSE]
                profiles.sizes[i] <- nrow(current_profile)
                means[i,]         <- apply(current_profile,2,mean)
                sd[i,]            <- apply(current_profile,2,sd)
                densities[i,]     <- apply(current_profile,2,function(x){ return(DENSITY(bin.width=bin.width,values=x))})
                for(j in 1:length(markers)){
					densities[i,j][[1]]@name <- markers[j]
                }
                message(paste0(round(i/profiles.nb*100),"%"))
            }
            message("done")
        }
        
        markers.clustering     <- rep(FALSE,length(markers))
        graph                  <- igraph::graph.empty(0,directed=FALSE)
        graph.layout           <- igraph::layout.auto(graph)
        colnames(graph.layout) <- NULL
        cluster <- CLUSTER(name = name,
            profiles            = profiles,
            profiles.nb         = profiles.nb,
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
)

