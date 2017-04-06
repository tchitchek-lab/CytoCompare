#' @title Plot for all S4 CytoCompare objects
#'
#' @description Makes visual representations for CELL, CLUSTER, GATE, MWEIGHTS or DENSITY objects.
#'
#' @details Profiles contained in CELL, CLUSTER, GATE objects can be represented alone (object2=NULL) or in combination with another CELL, CLUSTER, GATE objects (via object2).
#'
#' Cell profiles are represented via parallel coordinates where the x-axis represents the different markers and where the y-axis represents the marker expressions.
#' 
#' Cell cluster profiles are represented via parallel coordinates where the x-axis represents the different markers, where the y-axis represents the marker expressions, and where error bars indicate the marker expression standard deviations.
#' 
#' Gate profiles are represented via ribbons where the x-axis represents the different markers and where the y-axis represents the marker intensity ranges
#'
#' In the case of a single CELL object, the parameter 'overview' indicates if an overview of the CELL object must be represented (e.g. viSNE map). In the case of a single CLUSTER object, the parameter 'overview' indicates if an overview of the CLUSTER object must be represented (e.g. a SPADE tree). In both cases, the plot function will call the function indicated in the 'overview.function' slot of the CELL or CLUSTER objects. 
#' 
#' Marker weights contained in MWEIGHTS objects can be represented via bar plots where each bar corresponds to a marker and where the bar heights are proportional to the marker weights. 
#' 
#' Density profiles contained in DENSITY objects can be represented via histogram plots where each bar corresponds to a density bin and where a smooth line represents the average estimation of the marker expression density.
#'
#' If several profiles are present in the CELL, CLUSTER, GATE objects, all profiles or combination of profiles will be plotted. If several comparison results are present in the RES objects, all comparison results will be plotted.
#'
#' Comparison results contained in RES object can also be ploted using this function. In such representation, each bar corresponds to a marker with a height proportional to the marker distance or inclusion assessment, and where the bars are colored if they model a success.
#'
#' @param object1 a CELL, CLUSTER, GATE, MWEIGHTS or DENSITY object to plot
#' @param object2 another object to plot over the first object (only for CELL, CLUSTER or GATE objects)
#' @param overview a logical indicating is an overview of the object must be plotted (only available for single CELL and CLUSTER)
#' @param return.gg a logical indicating if the function should return a list of ggplot objects
#' @param ... other parameters
#'
#' @return if return.gg is TRUE, the function returns a list of ggplot objects
#' 
#' @name plot
#' @rdname plot-methods
#' @export 
setGeneric("plot", function(object1,object2,...){ standardGeneric("plot") })


#' @rdname plot-methods
#' @export
setMethod("plot",c("CLUSTER","missing"),
    function(object1,object2,overview=FALSE,return.gg=FALSE){
        plots <- list()
        if(overview){
            if(object1@overview.function!="")
                plots[[1]] <- do.call(object1@overview.function,list(object1))
            else
                stop("no function defined for plotting the CLUSTER profile")
        }else{
            for(i in 1:object1@profiles.nb){
                
                frame <- data.frame(profile = rep(object1@profiles[i],object1@markers.nb),
                                    markers = object1@markers,
                                    means   = object1@means[i,],
                                    sd      = object1@sd[i,], stringsAsFactors = FALSE)
                frame$range_min <- frame$means-frame$sd
                frame$range_max <- frame$means+frame$sd
                
                frame$markers[object1@markers.clustering] <- paste0(frame$markers[object1@markers.clustering],"_clust")
                
                frame$markers <- factor(frame$markers,levels = frame$markers)
                
                ymin <- min(frame$range_min)
                ymax <- max(frame$range_max)
                
                plots[[i]] <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="means",ymin="range_min",ymax="range_max",group="profile")) +
                    ggplot2::ggtitle(paste0("CLUSTER:",object1@name,":",object1@profiles[i]," (",object1@profiles.sizes[i]," cells)"))+
                    ggplot2::geom_errorbar(stat="identity",color="deepskyblue3")+
                    ggplot2::geom_line(stat="identity",color="deepskyblue3")+
                    ggplot2::geom_point(stat="identity",colour="deepskyblue3") +
                    ggplot2::xlab("markers") +
                    ggplot2::ylab("marker expression") +
                    ggplot2::scale_y_continuous(limits=c(floor(ymin),ceiling(ymax)),breaks=seq(floor(ymin),ceiling(ymax),by=1)) +
                    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
                        legend.position="bottom",
                        panel.background=ggplot2::element_blank(),
                        panel.grid.major.y=ggplot2::element_line(colour="gray75"),
                        panel.grid.major.x=ggplot2::element_line(colour="gray75"),
                        panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
            }
        }
        
        if(return.gg)
            return(plots)
        else
            multiplot(list = plots)
    }
)


#' @rdname plot-methods
#' @export
setMethod("plot",c("CLUSTER","CLUSTER"),
    function(object1,object2,return.gg=FALSE){
            plots <- list()
            plots.count <- 1
            for(i in 1:object1@profiles.nb){
                for(j in 1:object2@profiles.nb){
                    
                    markers <- intersect(object1@markers,object2@markers)
                    
                    frame.profile1 <- data.frame(profiles = rep(paste0("profile1 - CLUSTER:",object1@name,":",object1@profiles[i]," (",object1@profiles.sizes[i]," cells)"),length(markers)),
                                            markers  = markers,
                                            means    = object1@means[i,object1@markers %in% markers],
                                            sd       = object1@sd[i,object1@markers %in% markers], stringsAsFactors = FALSE)
                    frame.profile2 <- data.frame(profiles = rep(paste0("profile2 - CLUSTER:",object2@name,":",object2@profiles[j]," (",object2@profiles.sizes[j]," cells)"),length(markers)),
                                            markers  = markers,
                                            means    = object2@means[j,object2@markers %in% markers],
                                            sd       = object2@sd[j,object2@markers %in% markers], stringsAsFactors = FALSE)
                    
                    frame.profile1$markers <- factor(frame.profile1$markers,levels = frame.profile1$markers)
                    frame.profile2$markers <- factor(frame.profile2$markers,levels = frame.profile2$markers)
                    
                    frame           <- rbind(frame.profile1,frame.profile2)
                    frame$range_min <- frame$mean-frame$sd
                    frame$range_max <- frame$mean+frame$sd
                    
                    ymin <- min(frame$range_min)
                    ymax <- max(frame$range_max)
                
                    plots[[plots.count]] <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="means",ymin="range_min",ymax="range_max",color="profiles",group="profiles")) +
                        ggplot2::ggtitle(paste0("CLUSTER:",object1@name,":",object1@profiles[i]," (",object1@profiles.sizes[i]," cells)","\n","vs\nCLUSTER:",object2@name,":",object2@profiles[j]," (",object2@profiles.sizes[j]," cells)")) +
                        ggplot2::geom_point(stat="identity",size=2) +
                        ggplot2::geom_line(stat="identity",size=0.7) +
                        ggplot2::geom_errorbar() +
                        ggplot2::scale_color_manual(values=c("deepskyblue3","firebrick3")) +
                        ggplot2::xlab("markers") +
                        ggplot2::ylab("marker expression") +
                        ggplot2::scale_y_continuous(limits=c(floor(ymin),ceiling(ymax)),breaks=seq(floor(ymin),ceiling(ymax),by=1)) +
                        ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE)) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
                            plot.title=ggplot2::element_text(hjust=0.5),
                            legend.position="bottom",
                            panel.background=ggplot2::element_blank(),
                            panel.grid.major.y=ggplot2::element_line(colour="gray75"),
                            panel.grid.major.x=ggplot2::element_line(colour="gray75"),
                            panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
                    plots.count <- plots.count+1
                }
            }
                    
            if(return.gg)
                return(plots)
            else
                multiplot(list = plots)
    }   
)


#' @rdname plot-methods
#' @export
setMethod("plot",c("MWEIGHTS","missing"),
    function(object1,object2,return.gg=FALSE){
        
        frame <- data.frame(markers = object1@markers,
                            weight  = object1@weights,stringsAsFactors = FALSE)
        
        frame$markers <- factor(frame$markers,levels = frame$markers)
        
        plot <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="weight")) +
            ggplot2::ggtitle("MWEIGHTS") +
            ggplot2::geom_bar(stat="identity",colour="black",fill="darkorchid3") +
            ggplot2::xlab("markers") +
            ggplot2::ylab("weights") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
                legend.position="bottom",
                panel.background=ggplot2::element_blank(),
                panel.grid.major.y=ggplot2::element_line(colour="gray75"),
                panel.grid.major.x=ggplot2::element_line(colour="gray75"),
                panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
                
        if(return.gg)
            return(plot)
        else
            multiplot(list = list(plot))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("DENSITY","missing"),
    function(object1,object2,return.gg=FALSE){
        
        bin <- seq(object1@bin.interval[1],(object1@bin.interval[2]),by=object1@bin.width)

        frame <- data.frame(bin.min = bin[-sum(object1@bin.nb)],
                            bin.max = bin[-1],
                            value   = c(rev(object1@values.neg),object1@values.pos))
        
        frame$x      <- as.numeric(frame$bin.min)
        frame$smooth <- as.vector(stats::smooth(frame$value))
         
        plot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste0("DENSITY - ",object1@name))+
            ggplot2::geom_rect(data=frame,ggplot2::aes_string(xmin="bin.min",xmax="bin.max",ymax="value",ymin=0),colour="black",fill="deepskyblue3") +
            ggplot2::geom_line(data=frame,ggplot2::aes_string(x="x",y="smooth"),size=1,colour="darkblue") +
            ggplot2::xlab("expression") +
            ggplot2::ylab("density") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
                legend.position="bottom",
                panel.background=ggplot2::element_blank(),
                panel.grid.major.y=ggplot2::element_line(colour="gray75"),
                panel.grid.major.x=ggplot2::element_line(colour="gray75"),
                panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
        
        if(return.gg)
            return(plot)
        else
            multiplot(list = list(plot))
    }
)

#' @rdname plot-methods
#' @export
setMethod("plot",c("RES","missing"),
    function(object1,return.gg=FALSE,...){
        plots       <- list()
        plots.count <- 1
        for(i in 1:object1@comparisons.nb){
            
            
            frame <- data.frame(markers          = paste0(object1[i]@markers," (",object1[i]@marker.weights,")"),
                                marker.distances = as.vector(unlist(object1[i]@marker.distances)),
                                success          = as.vector(unlist(object1[i]@marker.successes)),stringsAsFactors = FALSE)
            
            colors            <- as.numeric(as.vector(unlist(object1[i]@marker.successes)))
            colors[colors==1] <- "green"
            colors[colors==0] <- "red"

            frame$markers <- factor(frame$markers,levels=frame$markers)
            
			if(object1@comparisons[i,"type"] == "similarity"){
			    
			    D_th <- object1@comparisons[i,"D_th"]
			    ymin <- min(frame$marker.distances)
			    ymax <- max(max(frame$marker.distances),D_th,1)
			    print(D_th)
				plots[[plots.count]] <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="marker.distances",fill="success")) +
										ggplot2::ggtitle(paste(object1[i]@comparisons[1,"profile1"]," vs. ",object1[i]@comparisons[1,"profile2"],"\n aggregated distance = ",format(object1[i]@comparisons[1,]$measure,nsmall=4,digit=4),"\n similarity p-value = ",format(object1[i]@comparisons[1,]$pvalue,nsmall=4,digit=4))) +
                         			    ggplot2::geom_hline(yintercept = D_th,size=1,linetype="dashed",color="grey20") +
										ggplot2::geom_bar(stat="identity",colour="black") +
										ggplot2::scale_fill_manual(values=c("TRUE"="green","FALSE"="red")) +
										ggplot2::scale_y_continuous(limits=c(floor(ymin),ceiling(ymax)),breaks=c(seq(floor(ymin),ceiling(ymax),by=1),D_th)) +
										ggplot2::xlab("markers") +
										ggplot2::ylab("distances (D)") +
										ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1, colour=colors),
													   legend.position="bottom",
													   panel.background=ggplot2::element_blank(),
													   panel.grid.major.y=ggplot2::element_line(colour="gray75"),
													   panel.grid.major.x=ggplot2::element_line(colour="gray75"),
													   panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
			}else if (object1@comparisons[i,"type"] == "inclusion"){
				plots[[plots.count]] <- ggplot2::ggplot(frame,ggplot2::aes_string(x="markers",y="success",fill="success")) +
										ggplot2::ggtitle(paste(object1[i]@comparisons[1,]$profile1," vs. ",object1[i]@comparisons[1,"profile2"],"\n inclusion p-value = ",format(object1[i]@comparisons[1,"pvalue"],nsmall=4,digit=4))) +
										ggplot2::geom_tile(color = "black") +
										ggplot2::scale_fill_manual(values=c("TRUE"="green","FALSE"="red")) +
				                        ggplot2::scale_y_discrete(expand=c(0,0.01)) +
										ggplot2::xlab("markers") +
										ggplot2::ylab("inclusions") +
										ggplot2::theme(axis.text.x       =ggplot2::element_text(angle=290, hjust=0, vjust=1, colour=colors),
													   legend.position   ="bottom",
													   panel.background  =ggplot2::element_rect(color="grey80", fill = NA),
													   axis.text.y       =ggplot2::element_blank(),
													   axis.ticks.y      =ggplot2::element_blank(),
													   panel.grid.major.y=ggplot2::element_blank(),
													   panel.grid.major.x=ggplot2::element_blank(),
													   panel.grid.minor.x=ggplot2::element_blank())
				}

            plots.count <- plots.count+1
        }
        
        if(return.gg)
            return(plots)
        else
            multiplot(list = plots)

    }
)


# title Internal - Ggplot representation of a SPADE tree
#
# @description Generates a ggplot representation of a SPADE tree stored in a CLUSTER profile which has been imported from a SPADE result.
#
# @details In such SPADE tree representation, each node of the graph corresponds to a SPADE cluster and nodes are linked based on theirs SPADE similarities using a minimum spanning tree approach. Node sizes are proportional to the numbers of cells associated to each SPADE cluster. This function is called when plotting a CLUSTER object imported from a SPADE result.
#
# @param profile a CLUSTER profile containing the SPADE tree structure and layout.
#
# @return a ggplot object of the SPADE tree representation
#
#' @importFrom igraph get.edgelist
#' @import ggrepel
ggplot.SPADEtree <- function(profile){
    
    graph                  <- profile@graph
    layout                 <- profile@graph.layout
    profiles.labels        <- profile@profiles
    profiles.cluster_sizes <- profile@profiles.sizes
    
    edge            <- igraph::get.edgelist(graph,names=FALSE)
    pos.edge        <- data.frame(xstart=layout[edge[,1],1],xend=layout[edge[,2],1],ystart=layout[edge[,1],2],yend=layout[edge[,2],2])
    pos.vertex      <- data.frame(id=profiles.labels,size=profiles.cluster_sizes,x=layout[,1],y=layout[,2])
    
    plot <- ggplot2::ggplot()+
        ggplot2::ggtitle("SPADE CLUSTER object overview")+
        ggplot2::geom_segment(data=pos.edge,ggplot2::aes_string(x="xstart",xend="xend",y="ystart",yend="yend"),colour="black")+
        ggplot2::geom_point(data=pos.vertex,ggplot2::aes_string(x="x",y="y",size="size"),colour="deepskyblue3")+
        ggrepel::geom_text_repel(data=pos.vertex,ggplot2::aes_string(x="x",y="y",label="id"),size=2)+
        ggplot2::theme(panel.grid.major=ggplot2::element_blank(),
            panel.grid.minor=ggplot2::element_blank(),
            panel.background=ggplot2::element_blank(),
            legend.text=ggplot2::element_blank(),
            axis.text=ggplot2::element_blank(),
            axis.title=ggplot2::element_blank(),
            axis.ticks=ggplot2::element_blank(),
            legend.position="none")
            
    return(plot)
}


#' @title Density heatmap of the marker expression
#'
#' @description Generates a marker density heatmap for the cell cluster profiles stored in a `CLUSTER` object.
#'
#' @details In such representation, each bar corresponds to a marker and the color gradient is proportional to the marker expression density. 
#'
#' @param cluster a CLUSTER object containing one or several cell cluster profiles
#' @param density.max a numeric specifying the maximal density gradient value in the heatmap   
#' @param return.gg a logical indicating if the function should return a list of ggplot objects
#'
#' @return if return.gg is TRUE, the function returns a list of ggplot objects
#'
#' @export
dheatmap <- function(cluster, density.max = NULL, return.gg=FALSE){
    stopifnot(density.max>0)
    plots <- list()
    markers <- cluster@markers
    markers[cluster@markers.clustering] <- paste0(markers[cluster@markers.clustering],"_clust")
    for(i in 1:cluster@profiles.nb){
        den.value <- NULL
        for(j in 1:length(markers)){
            temp.den  <- cluster@densities[[i,j]]
            bin.neg   <- seq(from = temp.den@bin.interval[1],to=0,by = temp.den@bin.width)
            bin.pos   <- seq(from = temp.den@bin.width,to=temp.den@bin.interval[2],by = temp.den@bin.width)
            x.lab     <- c(bin.neg,bin.pos)
            den.value <- rbind(den.value,
                data.frame(markers=rep(markers[j],sum(temp.den@bin.nb)),
                    bin.id.min=x.lab[-sum(temp.den@bin.nb)+1],
                    bin.id.max=x.lab[-1],
                    density=c(rev(temp.den@values.neg),temp.den@values.pos)))
        }
		max.value <- max(den.value$density,na.rm=TRUE)*0.65
		if(is.null(density.max)){
		    density.max <- max.value
		}
        density.breaks <- round(c(0,max.value,density.max),2)
        plots[[i]] <- ggplot2::ggplot(den.value)+
            ggplot2::ggtitle(paste0("CLUSTER:",cluster@name,":",cluster@profiles[i]))+
            ggplot2::geom_linerange(ggplot2::aes_string(x="markers",ymax="bin.id.max",ymin="bin.id.min",color="density"),size=5)+
            ggplot2::scale_color_gradientn(colours=c("black","darkgoldenrod3"),values=c(0,density.max/max.value),limits=c(0,max.value),na.value="darkgoldenrod3",breaks=density.breaks) +
            ggplot2::xlab("markers")+
            ggplot2::ylab("marker expression") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle=290, hjust=0, vjust=1),
                legend.position="bottom",
                panel.background=ggplot2::element_blank(),
                panel.grid.major.y=ggplot2::element_line(colour="gray75"),
                panel.grid.major.x=ggplot2::element_line(colour="gray75"),
                panel.grid.minor.x=ggplot2::element_line(colour="gray90"))
    }
    
    if(return.gg)
        return(plots)
    else
        multiplot(list = plots)
            
}