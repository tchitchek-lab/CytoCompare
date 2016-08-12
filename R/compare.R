#' @title Compare two cytometry profiles.
#'
#' @description Cytometry profiles contained in CELL, CLUSTER, or GATE objects can be compared using the `compare()` function. Comparison results are stored in a RES object. 
#' 
#' Comparisons can be performed between profiles of same types or between profiles of different types. In the default statistical approach:\cr
#' * if the comparisons are performed on profiles of same type then profiles will be compared to identify similar profiles\cr
#' * if the comparisons are performed on profiles of different types then profiles will be compared to identify included profiles
#'
#' In the context of a similarity comparison, a distance is computed between each marker of the two profiles (\code{$D_i$}). Marker distances below a distance threshold, specified by the user, will correspond to a marker similarity success. Euclidean distance will be used when comparing cell profiles while the Kolmogorov-Smirnov distance will be used when comparing cluster or gate profiles. A weight can be associated to each marker, in order to modulate their importance. An aggregation of marker distances is performed using an exact binomial test where marker successes are considered as successful Bernoulli experiments. Thereby, the proportion of marker successes is compared to a probability of success (\code{$P$}) specified by the user. A aggregated distance (\code{$D$}), corresponding to the weigthed mean of marker distances, is additionaly returned.
#'
#' In the context of an inclusion comparison, an inclusion assessment is performed for each marker of the two profiles. A cell profile marker is considered as included in a gate profile when its expression value is within the range of the marker boundaries. Similarly, a cell profile marker is considered as included in a cell cluster profile when its expression value is within the range of the marker cluster defined based on quantiles of marker expression densities. Finally, a cell cluster profile marker is considered as included in a gate profile when its expression boundaries is within the range of the marker gate. As for similarity comparisons, weights associated to each marker. The aggregation of marker inclusion is also performed using an exact binomial test.
#' 
#' Comparisons can be performed based on the whole set of common markers between the two profiles, or based on a subset of markers specified by the user. Moreover, markers can be weighted in the comparison procedure, via a MWEIGHTS object.
#'
#' If only one object is provided to the `compare()` function then the comparisons will be performed between all profiles of this object. If two objects are provided to the `compare()` function then the comparisons will be performed between all possible pairs of profiles between these two objects.
#'
#' Importantly, users can define their own function to perform the statistical comparisons of the profiles, using the `method` parameter. Please refer to the user tutorial for more details about this feature.
#'
#' @details Different parameters can be defined, via the method.params named list, to specify the behaviour of the comparisons:\cr
#' * the D.th parameter indicates the distance threshold\cr
#' * the P parameter indicates the expected proportion of marker successes\cr
#' * the nbcells.th parameter indicates the number of cells per cluster below which the marker expression density of a cell cluster profile will be approximated by a normal distribution\cr
#' * the cluster.quantiles parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profiles 
#' 
#' In the case of comparisons between two cell profiles, the marker distances are calculated based on the Euclidean distance. The parameter 'D.th' is set to 1.50 by default and the parameter 'P' is set to 0.75 default. 
#'
#' In the case of comparisons between two cell cluster profiles, the marker distances are calculated based on the Kolmogorov-Smirnov distance. The parameter 'D.th' is set to 0.30 by default and the parameter 'P' is set to 0.75 default. The nbcells.th parameter indicates the number of cells per cluster below which the density will be approximated by a normal distribution (set to 50 by default)
#'
#' In the case of comparisons between two gate profiles, gates are modeled by uniform distributions, and the marker distances are calculated based on the Kolmogorov-Smirnov distance. The parameter 'D.th' is set to 0.30 by default and the parameter 'P' is set to 0.75 default.
#'
#' In the case of comparisons between a cell profile and a gate profile, a cell profile marker is considered as included in the gate profile when its expression value is within the range of the marker boundaries. The parameter 'P' is set to 0.75 default.
#'
#' In the case of comparisons between a cell profile and cell cluster profile, a cell profile marker is considered as included in a cell cluster profile when its expression value is within the range of the marker cluster defined based on quantiles of marker expression densities. The parameter 'P' is set to 0.75 default. The 'cluster.quantiles' parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profile (set to 0.10 and 0.90 by default).
#'
#' In the case of comparisons between a cell cluster profile and gate profile, a cell cluster profile marker is considered as included in a gate profile when its expression boundaries is within the range of the marker gate. The parameter 'P' is set to 0.75 default. The 'cluster.quantiles' parameter indicates the quantiles that will define the marker expression ranges for the cell cluster profile (set to 0.10 and 0.90 by default).
#'
#' Importantly, in the case of comparisons involving CLUSTER profiles, Hartigan's dip tests and InterQuartile Ranges (IQR) can be computed in order to estimate if the marker expression densities are unimodales with low spreads. The Hartigan's dip test p-value threshold and IQR threshold can be both parametrized using the 'dip.pvalue' and 'IQR.th' parameters. If a marker density do not respect these constraints, the distance is set to 1.
#'
#' @param object1 a CELL, CLUSTER or GATE object
#' @param object2 a CELL, CLUSTER or GATE object
#' @param mweights a MWEIGHTS object specifying the markers to use in the comparison procedure with theirs associated weights
#' @param method a function or a character specifying the name of a function to use when performing the statistical comparisons between the cytometry profiles
#' @param method.params a named character list used to parametrize the comparison function (please see the details section)
#' @param ... other parameters
#'
#' @return a S4 object of class RES
#'
#' @name compare
#' @rdname compare-methods
#'
#' @export
setGeneric("compare", function(object1,object2,...){ standardGeneric("compare")})

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CELL",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.intensities = object1@intensities[j,object1@markers %in% markers],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params,method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object1@name,object1@profiles[j],sep=":"),type="similarity", check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count  <- count+1
            }
        }
        message("done")
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        object.density_tmp    <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type       = "CLUSTER",
                        profile1.mean       = object1@means[i,object1@markers %in% markers],
                        profile1.sd         = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells    = object1@profiles.sizes[i],
                        profile2.mean       = object1@means[j,object1@markers %in% markers],
                        profile2.sd         = object1@sd[j,object1@markers %in% markers],
                        profile2.nbcells    = object1@profiles.sizes[i],
                        weights             = mweights@weights)
                        
                if(object.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers],profile2.density=object1@densities[j,object1@markers %in% markers]))
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object1@name,object1@profiles[j],sep=":"),type="similarity",check.names = FALSE)

                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","missing"),
    function(object1,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object1)
        comparisons           <- NULL
        markers               <- mweights@markers
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object1@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object1@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object1@name,":",object1@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "GATE",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.range       = object1@ranges[j,object1@markers %in% markers,],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object1@name,object1@profiles[j],sep=":"),type="similarity",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")

        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CELL",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)

                res.sub               <- do.call(method,c(params,method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),type="similarity",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "CLUSTER",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),type="similarity",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "GATE",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),type="similarity",check.names = FALSE)
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")

        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing CELL:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "CLUSTER",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        comparisons           <- NULL
        markers               <- mweights@markers
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
       
       message(paste0("Comparing CLUSTER:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "CELL",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)
                
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                    
                res.sub              <- do.call(method,c(params, method.params))
                comparison           <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CELL","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing CELL:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CELL:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CELL",
                        profile2.type        = "GATE",
                        profile1.intensities = object1@intensities[i,object1@markers %in% markers],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CELL",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","CELL"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        
        message(paste0("Comparing GATE:",object1@name," vs. ","CELL:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","CELL:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "CELL",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.intensities = object2@intensities[j,object2@markers %in% markers],
                        weights              = mweights@weights)
                        
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("CELL",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("CLUSTER","GATE"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        profile1.density_tmp  <- prod(dim(object1@densities) == c(object1@profiles.nb,length(object1@markers)))
        
        message(paste0("Comparing CLUSTER:",object1@name," vs. ","GATE:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% CLUSTER:",object1@name,":",object1@profiles[i]," vs. ","GATE:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "CLUSTER",
                        profile2.type        = "GATE",
                        profile1.mean        = object1@means[i,object1@markers %in% markers],
                        profile1.sd          = object1@sd[i,object1@markers %in% markers],
                        profile1.nbcells     = object1@profiles.sizes[i],
                        profile2.range       = object2@ranges[j,object2@markers %in% markers,],
                        weights              = mweights@weights)
                        
                if(profile1.density_tmp)
                    params <- c(params,list(profile1.density=object1@densities[i,object1@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("CLUSTER",object1@name,object1@profiles[i],sep=":"),profile2=paste("GATE",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)

                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)

#' @rdname compare-methods
#' @export
setMethod("compare",c("GATE","CLUSTER"),
    function(object1,object2,
                mweights      = NULL,
                method        = "compare_default",
                method.params = NULL){
        if(is.null(mweights))
            mweights <- intersect(object1,object2)
        markers               <- mweights@markers
        comparisons           <- NULL
        comb.measure          <- NULL
        comb.pvalue           <- NULL
        comb.D_th             <- NULL
        comb.marker.distances <- NULL
        comb.marker.successes <- NULL
        profile2.density_tmp  <- prod(dim(object2@densities) == c(object2@profiles.nb,length(object2@markers)))
        
        message(paste0("Comparing GATE:",object1@name," vs. ","CLUSTER:",object2@name,"..."))
        count <- 1
        for(i in 1:object1@profiles.nb){
            for(j in 1:object2@profiles.nb){
                percent <- round(sqrt(count)/object1@profiles.nb*sqrt(count)/object2@profiles.nb*100)
                message(paste0(percent,"% GATE:",object1@name,":",object1@profiles[i]," vs. ","CLUSTER:",object2@name,":",object2@profiles[j]))
                
                params <- list(profile1.type = "GATE",
                        profile2.type        = "CLUSTER",
                        profile1.range       = object1@ranges[i,object1@markers %in% markers,],
                        profile2.mean        = object2@means[j,object2@markers %in% markers],
                        profile2.sd          = object2@sd[j,object2@markers %in% markers],
                        profile2.nbcells     = object2@profiles.sizes[j],
                        weights              = mweights@weights)
                        
                if(profile2.density_tmp)
                    params <- c(params,list(profile2.density=object2@densities[j,object2@markers %in% markers]))
                    
                res.sub               <- do.call(method,c(params, method.params))
                comparison            <- data.frame(profile1=paste("GATE",object1@name,object1@profiles[i],sep=":"),profile2=paste("CLUSTER",object2@name,object2@profiles[j],sep=":"),type="inclusion",check.names = FALSE)
                                        
                comparisons           <- rbind(comparisons,comparison)
                comb.measure          <- c(comb.measure,res.sub$measure)
                comb.pvalue           <- c(comb.pvalue,res.sub$pvalue)
                comb.D_th             <- c(comb.D_th,res.sub$D.th)
                comb.marker.distances <- rbind(comb.marker.distances,res.sub$marker.distances)
                comb.marker.successes <- rbind(comb.marker.successes,res.sub$marker.successes)
                
                count <- count+1
            }
        }
        message("done")
        
        res <- RES.format(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,mweights@weights)
        return(res)
    }
)


# title Internal - Creation of a RES object based the comparison results 
#
# @description This function is used internally to create a RES object based on calculated marker measures, marker p-values, aggregated measure and aggregated p-value. 
#
# @param comparisons a data.frame with two-columns containing the names of the profiles to compare
# @param markers a character vector containing the marker names
# @param comb.measure a numeric value containing the comparison measure
# @param comb.pvalue a numeric value containing the comparison p-value
# @param comb.D_th a numeric value containing the distance threshold
# @param comb.marker.distances a numeric vector containing the marker comparison measures
# @param comb.marker.successes a numeric vector containing the marker comparison success
#
# @return a RES object containing
RES.format <- function(comparisons,markers,comb.measure,comb.pvalue,comb.D_th,comb.marker.distances,comb.marker.successes,weights){
    
    comparisons[,"profile1"] <- as.character(comparisons[,"profile1"])
    comparisons[,"profile2"] <- as.character(comparisons[,"profile2"])
    comparisons$measure <- comb.measure
    comparisons$pvalue  <- comb.pvalue
    comparisons$D_th    <- comb.D_th
    
    rownames(comparisons)           <- paste(comparisons[,"profile1"],comparisons[,2],sep="/vs/")
    comb.marker.distances           <- data.frame(comb.marker.distances)
    comb.marker.successes           <- data.frame(comb.marker.successes)
    colnames(comb.marker.distances) <- markers
    rownames(comb.marker.distances) <- rownames(comparisons) 
    colnames(comb.marker.successes) <- markers
    rownames(comb.marker.successes) <- rownames(comparisons) 

    res <- RES(comparisons = comparisons,
        comparisons.nb     = nrow(comparisons),
        markers            = markers,
        marker.distances   = comb.marker.distances,
        marker.successes   = comb.marker.successes,
        marker.weights     = weights)
    return(res)
}
