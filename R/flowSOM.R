test <- function() {
    col.names <- c("CD279 BB515", "CD278 BB660", "CD127 BB700", "CD38 BB790", "TCF1 AF647", "Ki67 AF700", "KLRG1 APC-Fire", "CD223 BV480", "CD137 BV605",
                   "CD244 BV650", "CD366 BV711", "CD39 BV750", "CD28 BV786", "CD45RA BUV395", "CD8 BUV496", "CD185 BUV563", "CD25 BUV615", "CD226 BUV661", "CD27 BUV737", "CD4 BUV805",
                   "TIGIT PE", "Eomes PE-EF610", "CD152 PE-Cy5", "FoxP3 PE-Cy55", "tBET PE-Cy7", "CD197 BV421")

    fcs <- read.FCS("PICI0002_A00_K01304KE01_SPB_A01_002_CD3+.fcs")
    tab <- convert_fcs(fcs, asinh.cofactor = 150, negative.values = "truncate")
    m <- tab[, col.names]


    fsom <- list(
        data = as.matrix(m),
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        toTransform = NULL,
        transformFunction = NULL,
        scale = FALSE,
        prettyColnames = colnames(m),
        scaled.center = NULL,
        scaled.scale = NULL
    )
    class(fsom) <- "FlowSOM"

    ret <- FlowSOM::BuildSOM(fsom, colsToUse = 1:ncol(m))

}

#' @export
run_flowsom <- function(tab, col.names, scale = FALSE) {
    m <- as.matrix(tab[, col.names])

    fsom <- list(
        data = m,
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        toTransform = NULL,
        transformFunction = NULL,
        scale = FALSE,
        prettyColnames = setNames(colnames(m), colnames(m)),
        scaled.center = NULL,
        scaled.scale = NULL
    )

    if(scale) {
        fsom$data <- base::scale(x = fsom$data, center = TRUE, scale = TRUE)
        fsom$scaled.center <- attr(fsom$data, "scaled:center")
        attr(fsom$data, "scaled:center") <- NULL
        fsom$scaled.scale <- attr(fsom$data, "scaled:scale")
        attr(fsom$data, "scaled:scale") <- NULL
    }
    class(fsom) <- "FlowSOM"
    ret <- FlowSOM::BuildSOM(fsom, colsToUse = 1:ncol(m))
    return(ret)
}


