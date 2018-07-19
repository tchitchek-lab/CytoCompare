#' @title Exportation of CELL objects
#'
#' @description Exports a CELL object into a FCS file.
#'
#' @param object a CELL object
#' @param filename a character indicating the location of the FCS to save
#'
#' @return none
#' 
#' @importFrom flowCore write.FCS
#' @import flowUtils
#'
#' @name export
#' @rdname export-methods
setGeneric("export",function(object,filename) { standardGeneric("export") })

#' @rdname export-methods
#' @export
setMethod("export",c("CELL"),
    function(object,filename="cells.FCS"){
        data <- object@intensities
        if(object@trans!="none"){
            to.restore <- match(setdiff(object@markers,object@trans.exclude),object@markers)
            if(object@trans=="arcsinh"){
                data[,to.restore] <- sinh_coeff(data[,to.restore],coeff=object@trans.para$arcsinh.scale)
            }else if(object@trans=="log"){
                data[,to.restore] <- (object@trans.para$log.base ^ data[,to.restore]) - object@trans.para$log.shift
            }
        }
        flowFrame                 <- createFlowFrame(data,object@markers)
        suppressWarnings(filename <- flowCore::write.FCS(flowFrame,filename))
    }
)
