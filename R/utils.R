# title Internal - Inverse hyperbolic sine (arcsinh) transformation
#
# @description Computes the inverse hyperbolic sine (arcsinh) of a real number: log(sqrt((x/coeff)^2+x/coeff)), where 'coeff' is a rescaling coefficient (cofactor). The default value of 'coeff' is set to 5.
#
# @details This transformation is an alternative to the log transformation for cytometry data processing.
#
# @param x a numeric vector of values to transform
# @param coeff a numeric value indicating the rescaling coefficient (cofactor)
#
# @return a numeric vector of transformed values
arcsinh <- function(x,coeff=5){
    res <- log((x/coeff)+sqrt((x/coeff)^2+1))
    return(res)
}


# title Internal - hyperbolic sine (sinh) transformation
#
# @description Computes the hyperbolic sine (sinh) of a real number: (exp(x)-exp(-x))/2*coeff, where 'coeff' is a rescaling coefficient (cofactor). The default value of 'coeff' is set to 5.
#
# @details This transformation is used to convert values from an inverse hyperbolic sine (arcsinh) transformation
#
# @param x a numeric vector of values to transform
# @param coeff a numeric value indicating the rescaling coefficient (cofactor)
#
# @return a numeric vector of transformed values
sinh_coeff <- function(x,coeff=5){
    res <- (exp(x)-exp(-x))/2*coeff
    return(res)
}


# title Internal - Density computation
#
# @description Estimates the probability distribution of a continuous variable. This function is used internally to estimate the expression densities of the cell markers. The 'bin.width' parameter indicates the width of the bins, and is set to 0.05 by default.
#
# @details This function computes one histogram for the negative values and one histogram for the positive values. Please refer to the DENSITY class definition for more details about the provided results.
#
# @param x a numeric vector of values
# @param bin.width a numeric value specifying the width of the bins
#
# @return a named list containing the density computation results
compute.density <- function(x,bin.width=0.05){

	val.neg <- c()
    val.pos <- c()
    x.neg   <- x[x<0]
    x.pos   <- x[x>=0]
    lim.neg <- suppressWarnings(min(x.neg))
    lim.pos <- suppressWarnings(max(x.pos))
    n       <- length(x)
    if(abs(lim.neg) == Inf){
        val.neg <- 0
        inter.neg <- -bin.width
    }else{
        ind <- -bin.width
        while(ind > lim.neg-bin.width){
            val.neg <- c(val.neg,sum(x.neg>=ind)/n)
            x.neg   <- x.neg[x.neg<ind]
            ind     <- ind-bin.width
        }
        inter.neg <- ind+bin.width
    }
	
	if(abs(lim.pos) == Inf){
        val.pos <- 0
        inter.pos <- bin.width
    }else{
        ind <- bin.width
        while(ind < lim.pos+bin.width){
			val.pos <- c(val.pos,sum(x.pos<=ind)/n)
            x.pos   <- x.pos[x.pos>ind]
            ind     <- ind+bin.width
        }
        inter.pos <- ind-bin.width
    }
	
	if(is.null(val.pos))
		val.pos <- 0
    
    res <- list(point.nb = n,
    bin.interval         = c(inter.neg,inter.pos),
    bin.nb               = c(length(val.neg),length(val.pos)),
    values.pos           = val.pos,
    values.neg           = val.neg)
    
    return(res)
}


# title Internal - Display several ggplot objects into a single representation
#
# @description Plots vertically several ggplot objects into a single representation.
#
# @details This function is used when several cell, cell cluster or gate profiles are displayed at the same time.
#
# @param list a list of ggplot objects to be displayed vertically
# 
# @return none
#
#' @importFrom grid grid.newpage pushViewport viewport grid.layout 
multiplot <- function(list=NULL) {
    layout <- matrix(seq(1,length(list)),ncol=1,nrow=length(list))
    if(length(list) == 1){
      display.plot(list[[1]])
    }else{
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout), 1)))
        for(i in 1:length(list)){
            idx <- as.data.frame(which(layout==i, arr.ind=TRUE))
            vp  <- grid::viewport(layout.pos.row=idx$row, layout.pos.col=idx$col)
            display.plot(list[[i]], newpage=FALSE, vp=vp) 
        }
    }
}
#' @importFrom ggplot2 ggplot_build ggplot_gtable
display.plot <- function(x, newpage = is.null(vp), vp = NULL, ...){
    if (newpage) 
        grid::grid.newpage()
    grDevices::recordGraphics(requireNamespace("ggplot2", quietly = TRUE), list(), getNamespace("ggplot2"))
    data   <- ggplot2::ggplot_build(x)
    gtable <- ggplot2::ggplot_gtable(data)
    if (is.null(vp)){
        grid::grid.draw(gtable)
    }
    else{
        if(is.character(vp)) 
            grid::seekViewport(vp)
        else 
            grid::pushViewport(vp)
        grid::grid.draw(gtable)
        grid::upViewport()
    }
    invisible(data)
}


#' @title Cumulative distribution function of a DENSITY object
#'
#' @description Provides the cumulative distribution function (CDF) of a DENSITY object at specific values.
#'
#' @details This function is used internally when comparing the marker expression densities of cluster profiles with the default comparison approach. This function can also be used by users willing to define their own statistical functions when comparing cell, cell cluster or gate profiles.
#'
#' @param density a DENSITY object
#' @param values a numeric value or a numeric vector specifying the densities to compute
#'
#' @return a numeric value of the cumulative distribution values
#'
#' @export
cdf.density <- function(density,values){
    bins.neg    <- seq(from=density@bin.interval[1], to=0, by=density@bin.width)
    bins.pos    <- seq(from=density@bin.width, to=density@bin.interval[2], by=density@bin.width)
    bins        <- c(bins.neg,bins.pos)
    bins.length <- length(bins)
    bins.values <- c(rev(density@values.neg),density@values.pos)
    pdensities  <- NULL
    for(value in values){
        if(value >= bins[bins.length]){
            pdensities <- c(pdensities,1)
        }else if(value <= bins[2]){
            pdensities <- c(pdensities,0)
        }else{
            idx        <- which(bins >= value)[1]-2
            pdensities <- c(pdensities,sum(bins.values[1:idx]))
        }
    }
    return(pdensities)
}


#' @title Cumulative distribution function of a uniform distribution
#'
#' @description Provides the cumulative distribution function (CDF) of a uniform distribution at a specific value.
#'
#' @details This function is used internally when comparing the marker expression densities of gate profiles with the default comparison approach. This function can also be used by users willing to define their own statistical functions for comparing cell, cell cluster or gate profiles.
#'
#' @param bounds a numeric vector indicating the support (lower and upper bounds) of the uniform distribution
#' @param value a numeric value specifying the densities to compute
#'
#' @return a numeric value of the cumulative distribution value
#'
#' @export
cdf.uniform <- function(bounds,value){ 
    return(stats::punif(value,min=bounds[1],max=bounds[2]))
}

#' @title Quantiles of a DENSITY object
#'
#' @description Provides the quantiles of a DENSITY object, at specific values.
#'
#' @details This function is used internally when comparing the expression densities of cell markers with the proposed comparison approach. This function can also be used by users willing to define their own statistical functions for comparing cell, cell cluster or gate profiles.
#'
#' @param density a DENSITY object
#' @param values a numeric value or a numeric vector specifying the quantiles to compute
#'
#' @return a numeric vector of two values or a numeric matrix containing the bin ranges of the quantiles values
#'
#' @export
quantiles.density <- function(density,values){
    bins.neg    <- seq(from  = density@bin.interval[1],
                        to   = 0, 
                        by   = density@bin.width)
    bins.pos    <- seq(from  = density@bin.width, 
                        to   = density@bin.interval[2],
                        by   = density@bin.width)
    bins        <- c(bins.neg,bins.pos)
    bins.length <- length(bins)
    bins.values <- c(rev(density@values.neg), density@values.pos)
    bins.cumsum <- cumsum(bins.values)
    quantiles   <- NULL
    for(value in values){
        if(value <= bins.cumsum[1]){
            quantiles <- rbind(quantiles, rep(bins[1],2))
        }else if(value >= bins.cumsum[bins.length-1]){
            quantiles <- rbind(quantiles, rep(bins[bins.length],2))
        }else{
            idx       <- which(bins.cumsum >= value)[1]
            quantiles <- rbind(quantiles, bins[c(idx-1,idx)])
        }
    }
    return(quantiles)
}

#' @title Quantiles of a uniform distribution
#'
#' @description Provides the quantiles of a uniform distribution
#'
#' @details This function is used internally when comparing the marker expression ranges of gate profiles with the default comparison approach. This function can also be used by users willing to define their own statistical functions for comparing cell, cell cluster or gate profiles.
#'
#' @param bounds a numeric vector indicating the support (lower and upper bounds) of the uniform distribution
#' @param value a numeric value specifying the quantile to compute
#'
#' @return a numeric value of the quantiles value
#'
#' @export
quantiles.uniform <- function(bounds,value){
    return(stats::qunif(value,min=bounds[1],max=bounds[2]))
}


#' @title Creation of a MWEIGHTS object
#'
#' @description Creates a MWEIGHTS object based on a set of marker names, where all marker weights are set to 1.
#'
#' @details This function is a short-cut to the following code:\cr
#' markers        <- c("marker1","marker2","marker3","markeri...","markern")\cr
#' weights        <- rep(1,length(markers))\cr
#' mweights       <- MWEIGHTS(markers=markers,weights=weights)\cr\cr
#' The weight of each marker is set to 1 but can be changed afterwards using the function set.
#'
#' @param markers a character vector specifying the names of the markers
#'
#' @return a MWEIGHTS object containing the marker weights
#'
#' @export
create.MWEIGHTS <- function(markers){
    weights  <- rep(1,length(markers))
    wmarkers <- MWEIGHTS(markers=markers,weights=weights)
    return(wmarkers)
}


#' @title Retrieving of an example dataset of CytoCompare objects
#'
#' @description Downloads and loads an example dataset of CytoCompare objects constructed based on cytometry profiles obtained from healthy human bone marrow unstimulated or stimulated (PMID:21964415).
#'
#' This example dataset consists on three cytometry profiles of healthy human bone marrow, unstimulated or stimulated by BCR-inductor or IL-7, measured using a mass cytometry panel of more than 30 cell markers. This panel has been designed to identify a large spectrum of immune cell types like monocytes, B, or CD4+ and CD8+ T cells. A SPADE analysis has been performed to identify cell clusters, that have been then manually labelled based on theirs profiles. SPADE cell clusters corresponding to 6 majors cell types have been extracted and a set of rectangle gates have been constructed based these cell types.
#'
#' Once downloaded, the following objects will be available:\cr
#' * `bm_example.cells.b`, a `CELL` object containing the cell profiles of the B cell populations;\cr
#' * `bm_example.cells.mono`, a `CELL` object containing the cell profiles of the monocyte cell populations;\cr
#' * `bm_example.cells.tCD4naive`, a `CELL` object containing the cell profiles of the naive CD4+ T cell populations;\cr
#' * `bm_example.cells.tCD8naive`, a `CELL` object containing the cell profiles of the naive CD8+ T cell populations;\cr
#' * `bm_example.cells.tCD4mem`, a `CELL` object containing the cell profiles of the memory CD4+ T cell populations;\cr
#' * `bm_example.cells.tCD8mem`, a `CELL` object containing the cell profiles of the memory CD8+ T cell populations;\cr
#' * `bm_example.clusters`, a `CLUSTER` object containing the cell cluster profiles for all the different cell populations, identified by SPADE;\cr
#' * `bm_example.clusters.b`, a `CLUSTER` object containing the cell cluster profiles of the B cell populations, identified by SPADE;\cr
#' * `bm_example.clusters.mono`, a `CLUSTER` object containing the cell cluster profiles of the monocyte cell cluster profiles, identified by SPADE;\cr
#' * `bm_example.clusters.tCD4naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD4+ T cell populations, identified by SPADE;\cr
#' * `bm_example.clusters.tCD8naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD8+ T cell populations, identified by SPADE;\cr
#' * `bm_example.clusters.tCD4mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD4+ T cell populations, identified by SPADE;\cr
#' * `bm_example.clusters.tCD8mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD8+ T cell populations, identified by SPADE;\cr
#' * `bm_example.gates`, a `GATE` object containing the gate profiles constructed based on the six main cell populations identified by SPADE;\cr
#' * `bm_example.gates.b`, a `GATE` object containing the gate profiles constructed based on the B cell populations;\cr
#' * `bm_example.gates.mono`, a `GATE` object containing the gate profiles constructed based on the monocyte cell populations;\cr
#' * `bm_example.gates.tCD4naive`, a `GATE` object containing the gate profiles constructed based on the naive CD4T cell populations;\cr
#' * `bm_example.gates.tCD8naive`, a `GATE` object containing the gate profiles constructed based on the naive CD8T cell populations;\cr
#' * `bm_example.gates.tCD4mem`, a `GATE` object containing the gate profiles constructed based on the memory CD4T cell populations;\cr
#' * `bm_example.gates.tCD8mem`, a `GATE` object containing the gate profiles constructed based on the memory CD8T cell populations;\cr
#' * `bm_example.mweights`, a `MWEIGHTS` object containing cell markers that can be used in for comparison computations;\cr
#' * `bm_example.visne`, a list of three `CELL` objects containing the viSNE cell profiles for each biological sample.
#'
#' @details This function downloads a CytoCompareExample.rdata file containing the different CytoCompare objects (from a public ftp server "ftp://ftp.cytocompare.org/public/rdata/").
#'
#' @param del.file a logical specifying if the download CytoCompareExample.rdata file must be erased after the loading
#'
#' @return none
#'
#' @export
load.examples <- function(del.file=FALSE){
    rdata        <- "./CytoCompareExample.rdata"
	expected.md5 <- "496d2451cb781dd1ecf7b7dfd245cec6"

    if(!file.exists(rdata) || tools::md5sum(rdata)!=expected.md5){
        if(del.file){
            rdata <- tempfile(fileext = ".rdata")
            rdata <- gsub("\\\\","/",rdata)
        }
        utils::download.file("ftp://cytocompare:cytocompare@ftp.cytocompare.org/public/rdata/CytoCompareExample.rdata",rdata,mode="wb")
    }
    obj.loaded <- load(rdata,.GlobalEnv)
    if(del.file){
        file.remove(rdata)
    }
    message("The following objects have been loaded:")
    for(obj in obj.loaded){
        message("        ",obj)
    }
}


# @title Internal - Creation of a FlowFrame object
#
# @description Creates a FlowFrame object based on a numeric matrix of cell profiles. 
# 
# @details This function is used internally when writing cell profiles contained in a CELL object into a FCS file.
#
# @param intensities a numeric matrix corresponding to the cell profile intensities.
# @param markers a character vector corresponding to marker names.
#
# @return none
createFlowFrame <- function(intensities,markers) {
    
    colnames(intensities) <- markers
    p                     <- c()    
    description           <- list() 
    
    description[["$DATATYPE"]] <- "F"
    
    for (i in 1:ncol(intensities)) {
        name  <- markers[i]
        min   <- min(intensities[,i])
        max   <- max(intensities[,i])
        range <- max-min+1
        
        l           <- matrix(c(name,name,range,min,max),nrow=1)
        colnames(l) <- c("name","desc","range","minRange","maxRange")
        rownames(l) <- paste0("$P",i) 
        p           <- rbind(p,l)
        
        description[[paste("$P",i,"N",sep="")]] <- name;
        description[[paste("$P",i,"S",sep="")]] <- name;
        description[[paste("$P",i,"R",sep="")]] <- toString(range);
        description[[paste("$P",i,"B",sep="")]] <- "32";
        description[[paste("$P",i,"E",sep="")]] <- "0,0";
    }
    
    dataframe <- as(data.frame(p), "AnnotatedDataFrame")
    flowframe <- flowCore::flowFrame(intensities, dataframe, description=description)
    
    return(flowframe)
}


# @title Internal - Assessing if a DENSITY objet is unimodal
#
# @description This function is used internally to report if the distribution of the DENISTY object is unimodal using an Hartigan's dip test. 
# 
# @details H0 hypothesis of the Hartigan's dip test corresponds to unimodality and H1 to non-unimodality (multimodality)
#
# @param DENSITY a DENSITY object to be assessed
# @param th.pvalue a numeric specifying the pvalue threshold of the test
#
# @return a logical specifying if the DENSITY objet is unimodal or not
is.unimodal <- function(DENSITY,th.pvalue){
    expression.vector <- DENSITY.to_vector(DENSITY)
    p.values          <- diptest::dip.test(expression.vector)$p.value
    return(ifelse(p.values>th.pvalue,TRUE,FALSE))
}


# @title Internal - Assessing if a DENSITY objet has low spread 
#
# @description This function is used internally to report if the distribution of the DENISTY object has a spread below a threshold. The spread is mesured using the IQR (Inter Quartile Range). 
#
# @param DENSITY a DENSITY object to be assessed
# @param IQR.th a numeric specifying the maximal IQR threshold
#
# @return a logical specifying if the DENSITY objet has low spread or not
has.low_spread <- function(DENSITY,IQR.th=2){
    quartiles <- quantiles.density(DENSITY,c(0.25,0.75))
    IQR       <- max(quartiles) - min(quartiles)
    return(ifelse(IQR<IQR.th,TRUE,FALSE))
}

# @title Internal - Generating random values following a specified distribution
#
# @description This function is used internally to generate a numerical vector of n values which follow a specific distribution provided by a DENSITY object.
#
# @param DENSITY a DENISTY object containing a distribution
# @param n a numeric specifying the desired number of values 
#
# @return a vector of n values which follow the provided density
DENSITY.to_vector <- function(DENSITY,
                              n=DENSITY@point.nb){
    vector <- c()
    half   <- DENSITY@bin.width/2
    values <- c(rev(DENSITY@values.neg),DENSITY@values.pos)
    i      <- 1
	nmin   <- length(DENSITY@values.pos)+length(DENSITY@values.neg)
	size   <- n

	if(n<nmin)	n <- nmin
	
    for(bin.mean in seq(DENSITY@bin.interval[1]+half,DENSITY@bin.interval[2]-half,by=DENSITY@bin.width)){
        nb     <- round(n * values[i])
        vector <- c(vector,stats::rnorm(nb,mean=bin.mean,sd=half))
        i <- i+1
    }
	
	if(n<nmin)	vector <- sample(x=vector,size=size)
	
    return(vector)
}

# @title Internal - Rescaling of a vector of marker expression intensities
#
# @description This function is used internally to rescale a vector of marker expression intensities to a new interval, by keeping the distribution proportion.
#
# @details By default, the vector is rescaled to have a range going from 0 to 1. Then, the minimal density value will be set to 0 and the maximal one will be set to 1. Instead of minimal and maximal values, it is possible to provide quantiles of the distribution. 
#
# @param x a vector of marker expression intensities to be rescaled
# @param interval.after a numeric vector of two values specifying the range of the new interval
# @param quantiles a numeric vector of two values specifying the quantiles of the vector used to rescale (please refer to the detail section)
#
# @return a new rescaled vector of marker expression intensities
rescale.intensities <- function(x,
								interval.after = c(0,1),
								quantiles      = c(0,1)){

	stopifnot(is.numeric(x))
	stopifnot(is.numeric(interval.after))
    stopifnot(is.numeric(quantiles))
    stopifnot(length(interval.after)==2)
    stopifnot(length(quantiles)==2)
    stopifnot(interval.after[1]<interval.after[2])
    stopifnot(quantiles[1]<quantiles[2])
					
	interval.to_scale <- stats::quantile(x,quantiles)
	
	range.after    <- interval.after[2]  - interval.after[1]
    range.before   <- interval.to_scale[2] - interval.to_scale[1]
	rescale.factor <- range.after/range.before
	shift          <- interval.after[1] - (interval.to_scale[1]*rescale.factor)
	x              <- (x*rescale.factor) + shift
	
	return(x)			
}
