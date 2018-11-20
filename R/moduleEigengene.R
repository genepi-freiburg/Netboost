moduleEigengenes <- function (expr, colors, nPC = 1, align = "along average", 
    scale = TRUE, verbose = 0, indent = 0) 
{
    spaces = indentSpaces(indent)
    if (verbose == 1) 
        printFlush(paste(spaces, "moduleEigengenes: Calculating", 
            nlevels(as.factor(colors)), "module eigengenes in given set."))
    if (is.null(expr)) {
        stop("moduleEigengenes: Error: expr is NULL. ")
    }
    if (is.null(colors)) {
        stop("moduleEigengenes: Error: colors is NULL. ")
    }
    if (is.null(dim(expr)) || length(dim(expr)) != 2) 
        stop("moduleEigengenes: Error: expr must be two-dimensional.")
    if (dim(expr)[2] != length(colors)) 
        stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per variable).")
    if (is.factor(colors)) {
        nl = nlevels(colors)
        nlDrop = nlevels(colors[, drop = TRUE])
        if (nl > nlDrop) 
            stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                "Use colors[, drop=TRUE] to get rid of them."))
    }
    alignRecognizedValues = c("", "along average")
    if (!is.element(align, alignRecognizedValues)) {
        printFlush(paste("ModulePrincipalComponents: Error:", 
            "parameter align has an unrecognised value:", align, 
            "; Recognized values are ", alignRecognizedValues))
        stop()
    }
    maxVarExplained = 10
    if (nPC > maxVarExplained) 
        warning(paste("Given nPC is too large. Will use value", 
            maxVarExplained))
    nVarExplained = min(nPC, maxVarExplained)
    modlevels = levels(factor(colors))

    PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
        ncol = length(modlevels)))
    averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
    varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
    validAEs = rep(FALSE, length(modlevels))
    validColors = colors
    names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
        sep = "")
    names(averExpr) = paste("AE", modlevels, sep = "")
    for (i in c(1:length(modlevels))) {
        if (verbose > 1) 
            printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                modlevels[i]))
        modulename = modlevels[i]
        restrict1 = as.character(colors) == as.character(modulename)
        if (verbose > 2) 
            printFlush(paste(spaces, " ...", sum(restrict1), 
                "variables"))
        datModule = as.matrix(t(expr[, restrict1]))
        n = dim(datModule)[1]
        p = dim(datModule)[2]
        pc = try({
            if (verbose > 5) 
                printFlush(paste(spaces, " ...scaling"))
            if (scale) 
                datModule = t(scale(t(datModule)))
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating SVD"))
            svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, 
                p, nPC))
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating PVE"))
            veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], 
                t(datModule), use = "p")
            varExpl[c(1:min(n, p, nVarExplained)), i] = rowMeans(veMat^2, 
                na.rm = TRUE)
            svd1$v[, 1]
        }, silent = TRUE)
        PrinComps[, i] = pc
        ae = try({
             scaledExpr = scale(t(datModule))
             averExpr[, i] = rowMeans(scaledExpr, na.rm = TRUE)
             if (align == "along average") {
                if (verbose > 4) 
                  printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
                corAve = cor(averExpr[, i], PrinComps[, i], 
                  use = "p")
                if (!is.finite(corAve)) 
                  corAve = 0
                if (corAve < 0) 
                  PrinComps[, i] = -PrinComps[, i]
              }
              0
            if (class(ae) == "try-error") {
                if (verbose > 0) {
                  printFlush(paste(spaces, " ..Average expression calculation of module", 
                    modulename, "failed with the following error:"))
                  printFlush(paste(spaces, "     ", ae, spaces, 
                    " ..the returned average expression vector will be invalid."))
                }
                warning(paste("Average expression calculation of module", 
                  modulename, "failed with the following error \n     ", 
                  ae, "The returned average expression vector will be invalid.\n"))
            }
            validAEs[i] = !(class(ae) == "try-error")
        }
    }
    allAEOK = (sum(!validAEs) == 0)
    list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, 
        nPC = nPC, validAEs = validAEs, allAEOK = allAEOK)
}
