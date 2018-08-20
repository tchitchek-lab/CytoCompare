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
