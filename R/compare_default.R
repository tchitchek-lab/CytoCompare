# title Internal - Comparison of two CELL, CLUSTER or gate profiles
#
# @description This function is a wrapper for the default comparison methods
#
# @param profile1.type a character specifying the type of the first profile (CELL, CLUSTER or GATE)
# @param profile2.type a character specifying the type of the second profile (CELL, CLUSTER or GATE)
# @param profile1.intensities a numeric vector containing the marker expressions for the first cell profile
# @param profile1.mean a numeric vector containing the marker expression means of the first cluster profile
# @param profile1.sd a numeric vector containing the marker expression standard deviations of the first cluster profile
# @param profile1.density a named list containing the marker expression densities (DENSITY objects) of the first cluster profile
# @param profile1.nbcells a numeric value containing the number of cells associated with the first cluster profile
# @param profile1.range a numeric matrix containing the marker expression ranges of the second gate profile
# @param profile2.intensities a numeric vector containing the marker expressions for the second cell profile
# @param profile2.mean a numeric vector containing the marker expression means of the second cluster profile
# @param profile2.sd a numeric vector containing the marker expression standard deviations of the second cluster profile
# @param profile2.density a named list containing the marker expression densities (DENSITY objects) of the second cluster profile
# @param profile2.nbcells a numeric value containing the number of cells associated with the second cluster profile
# @param profile2.range a numeric matrix containing the marker expression ranges of the second gate profile
# @param ... other parameters for the 'compare_default' methods
#
# @return a named list containing the marker measures, marker successes, the aggregated measure and the aggregated p-value
compare_default <- function(profile1.type,profile2.type,
                            profile1.intensities,profile1.mean,profile1.sd,profile1.density,profile1.nbcells,profile1.range,
                            profile2.intensities,profile2.mean,profile2.sd,profile2.density,profile2.nbcells,profile2.range,
                            ...){  
    if(profile1.type=="CLUSTER" && profile2.type=="CLUSTER"){
        res <- compare_default_CLUSTER_CLUSTER(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,profile2.mean,profile2.sd,profile2.density,profile2.nbcells,...)
    }else{
        stop("")
    }
    return(res)
}

# title Internal - Comparison of two cluster profiles
#
# @description Performs a similarity comparison between two cluster profiles.
#
# @details This function computes the marker distances as the Kolmogorov-Smirnov statistics between the marker expression densities of the two cell cluster profiles. The aggregation of marker distances is performed using an exact binomial test.
#
# @param profile1.mean a numeric vector containing the marker expression means of the first cluster profile
# @param profile1.sd a numeric vector containing the marker expression standard deviations of the first cluster profile
# @param profile1.density a named list containing the marker expression densities (DENSITY objects) of the first cluster profile
# @param profile1.nbcells a numeric value containing the number of cells in the first cluster profile
# @param profile2.mean a numeric vector containing the marker expression means of the second cluster profile
# @param profile2.sd a numeric vector containing the marker expression standard deviations of the second cluster profile
# @param profile2.density a named list containing the marker expression densities (DENSITY objects) of the second cluster profile
# @param profile2.nbcells a numeric value containing the number of cells in the second cluster profile
# @param weights a numeric vector containing the marker weights
# @param D.th a numeric value specifying a threshold below which a marker distance will be considered as a similarity success
# @param P a numeric value specifying the proportion of marker successes to statistically overtake to consider the similarity as significant
# @param nbcells.th a numeric value specifying a threshold below which the marker expression density will be approximated by a normal distribution
# @param dip.pvalue a numeric specifying the pvalue threshold of the test of unimodality of marker expression densities (H0 meaning unimodal and H1 meaning multimodal)
# @param IQR.th a numeric specifying the maximal IQR threshold of marker expression densities
#
# @return a named list containing the marker distances, the used D.th threshold, the marker similarity successes, the aggregated distance and the aggregated similarity p-value
compare_default_CLUSTER_CLUSTER <- function(profile1.mean,profile1.sd,profile1.density,profile1.nbcells,
                                            profile2.mean,profile2.sd,profile2.density,profile2.nbcells,
                                            weights,
                                            D.th=0.25,P=1,
                                            nbcells.th=50,
                                            dip.pvalue=NULL,
                                            IQR.th    =NULL){
    SD_TO_IQR <- 1.349

    if(missing(profile1.density) || profile1.nbcells <= nbcells.th){
        profile1.fun.cdf    <- function(x,i){ stats::pnorm(x,profile1.mean[i],profile1.sd[i]) }
        profile1.fun.q      <- function(x,i){ stats::qnorm(x,profile1.mean[i],profile1.sd[i]) }
        if(!is.null(IQR.th)){
            profile1.unimodal   <- rep(TRUE,length(profile1.mean))
        }
        if(!is.null(dip.pvalue)){
            profile1.low_spread <- ifelse((profile1.sd*SD_TO_IQR)<IQR.th,TRUE, FALSE)
        }
    }else{
        profile1.fun.cdf    <- function(x,i){ cdf.density(profile1.density[[i]],x) }
        profile1.fun.q      <- function(x,i){ quantiles.density(profile1.density[[i]],x) }
        if(!is.null(IQR.th)){
            profile1.unimodal   <- sapply(profile1.density,FUN=is.unimodal,dip.pvalue)
        }
        if(!is.null(dip.pvalue)){
            profile1.low_spread <- sapply(profile1.density,FUN=has.low_spread,IQR.th)
        }
    }
    if(missing(profile2.density) || profile2.nbcells <= nbcells.th){
        profile2.fun.cdf    <- function(x,i){ stats::pnorm(x,profile2.mean[i],profile2.sd[i]) }
        profile2.fun.q      <- function(x,i){ stats::qnorm(x,profile2.mean[i],profile2.sd[i]) }
        if(!is.null(IQR.th)){
            profile2.unimodal   <- rep(TRUE,length(profile2.mean))
        }
        if(!is.null(dip.pvalue)){
            profile2.low_spread <- ifelse((profile2.sd*SD_TO_IQR)<IQR.th,TRUE, FALSE)
        }
    }else{
        profile2.fun.cdf    <- function(x,i){ cdf.density(profile2.density[[i]],x) }
        profile2.fun.q      <- function(x,i){ quantiles.density(profile2.density[[i]],x) }
        if(!is.null(IQR.th)){
            profile2.unimodal   <- sapply(profile2.density,FUN=is.unimodal,dip.pvalue)
        }
        if(!is.null(dip.pvalue)){
            profile2.low_spread <- sapply(profile2.density,FUN=has.low_spread,IQR.th)
        }
    }

    if(profile1.nbcells>=1 && profile2.nbcells>=1){
        measures   <- KS.stat(profile1.fun.cdf = profile1.fun.cdf,
                            profile1.fun.q     = profile1.fun.q,
                            profile2.fun.cdf   = profile2.fun.cdf,
                            profile2.fun.q     = profile2.fun.q,
                            markers.nb         = length(weights))
    }else{
        measures   <- rep(Inf,length(weights))
    }
    if(!is.null(IQR.th)){
        measures[!(profile1.unimodal & profile2.low_spread)] <- 1
    }
    if(!is.null(dip.pvalue)){
        measures[!(profile1.low_spread & profile2.unimodal)] <- 1
    }
    
    successes  <- (measures <= D.th)
    
    measure   <- stats::weighted.mean(measures,weights,na.rm = TRUE)
    pvalue    <- aggreg.binom(successes,weights,P)

    res <- list(measure  = measure,
                pvalue           = pvalue,
                marker.distances = measures,
                D.th             = D.th,
                marker.successes = successes)
    return(res)
}


# title Internal - Weighted binomial test for computations of profile similarity or inclusion p-values
#
# @description Performs a weighted binomial test to compute the significance of profile similarities or profile inclusions. 
#
# @details In the default statistical approach, each marker similarity or inclusion measure below a specific threshold value models a success in a Bernoulli experiment. The computed p-value represents then the significance of the proportion of similar or included cell markers when comparing two cytometry profiles. Each success or fail is repeated the number of times indicated by the 'weights' parameter. 
# 
# @param successes a boolean vector containing the marker successes (TRUE) or fail (FALSE)
# @param weights an integer vector containing the marker weights 
# @param P a numeric value specifying the expected success probability (set at 0.75 by default)
#
# @return a numeric providing the p-value of the weighted binomial test
aggreg.binom <- function(successes,weights,P=1){
    successes <- successes[!is.na(successes)]
    weights   <- weights[!is.na(successes)]
    pvalue    <- stats::binom.test(sum(successes*weights),sum(weights), P,alternative="less")$p.value
	return(pvalue)
}


# title Internal - Computation the Kolmogorov-Smirnov statistics of cell markers
#
# @description Computes the Kolmogorov-Smirnov statistics for each marker of two cytometry profiles.
#
# @details The Kolmogorov-Smirnov statistics corresponds to the maximal absolute difference between the cumulative distributions of two empirical distribution functions. A large statistics indicates that the two probability distributions are different, while a small statistics indicates that they are similar.
#
# @param profile1.fun.cdf the cumulative density function of the first cell cluster object  
# @param profile1.fun.q the quantile function of the first cell cluster object (used to obtain the marker expression ranges)
# @param profile2.fun.cdf the cumulative density function of the second cell cluster object 
# @param profile2.fun.q the quantile function of the second cell cluster object (used to obtain the marker expression ranges)
# @param markers.nb a numeric indicating the number of markers in both objects
#
# @return a number vector containing the Kolmogorov-Smirnov statistics for each marker
KS.stat <- function(profile1.fun.cdf,profile1.fun.q,profile2.fun.cdf,profile2.fun.q,markers.nb){
    measures <- numeric(markers.nb)
    bounds  <- c(0.0001,0.9999)
    for(i in 1:markers.nb){
        lower_bound <- min(profile1.fun.q(0,i),profile2.fun.q(0,i))
        if(is.infinite(lower_bound))
            lower_bound <- min(profile1.fun.q(bounds[1],i),profile2.fun.q(bounds[1],i))
        upper_bound <- max(profile1.fun.q(1,i),profile2.fun.q(1,i))
        if(is.infinite(upper_bound))
            upper_bound <- max(profile1.fun.q(bounds[2],i),profile2.fun.q(bounds[2],i))
        fun <- function(x){
            return(abs(profile1.fun.cdf(x,i)-profile2.fun.cdf(x,i)))
        }
        opt        <- stats::optimize(f=fun,maximum=TRUE,interval=c(lower_bound,upper_bound))
        measures[i] <- opt$objective
    }
    return(measures)
}

