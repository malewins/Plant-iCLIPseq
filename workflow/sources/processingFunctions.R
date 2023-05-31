pureClipGlobalFilter <- function(object, # bindingSiteFinder
                                 scoreCol = "score",
                                 cutoff = 0.01# defines the cutoff for which PureCLIP sites to keep; smallest are 1% (0.01); a cutoff of 0.05 means to remove the lowest 5% PureCLIP sites base on the score
) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    # check for correct cutoff value
    possibleCutoffs = seq(from = 0, to = 1, by = 0.01)
    if (! cutoff %in% possibleCutoffs) {
        msg = paste0("Specified cutoff (", cutoff,
                     "), does not refer to a percentage value. Try cutoffs of the form 0.01, 0.05, ... instead.")
        stop(msg)
    }
    # Check for correct score column
    metaNames = colnames(mcols(getRanges(object)))
    if (! scoreCol %in% metaNames) {
        msg = paste0("Specified scoreCol (", scoreCol, "), is present in the ranges. Add respective column. ")
        stop(msg)
    }
    
    # Create matching vectors for score column from input source
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    score = mcols(rngInitial)[match(scoreCol, colnames(mcols(rngInitial)))][[1]]
    
    # ---
    # Store function parameters in list
    optstr = list(cutoff = cutoff, scoreCol = scoreCol)
    object@params$pureClipGlobalFilter = optstr
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # compute quantile based on PureCLIP score
    qCut = quantile(score, probs = seq(0,1, by = 0.01))
    
    # apply user set cutoff
    names(qCut) = as.numeric(sub("%", "", names(qCut))) / 100
    idx = match(cutoff, names(qCut))
    rngFiltered = rngInitial[score >= qCut[idx]]
    
    # ---
    # Store for plotting
    dfPlot = score
    cutoffPoint = qCut[idx]
    object@plotData$pureClipGlobalFilter$data = dfPlot
    object@plotData$pureClipGlobalFilter$cutoffPoint = cutoffPoint
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "pureClipGlobalFilter()", class = "crosslink sites",
        nIn = length(rngInitial), nOut = length(rngFiltered),
        per = paste0(round(length(rngFiltered)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", optstr$cutoff, ", scoreMatchCol=", optstr$scoreCol)
    )
    object@results = rbind(object@results, resultLine)
    
    # return BSF object with only ranges above cutoff
    objectNew = setRanges(object, rngFiltered)
    return(objectNew)
}

pureClipGeneWiseFilter <- function(object, # bindingSiteFinder
                                   annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then genes must be provided
                                   genes = NULL, # ranges of the genes of class GenomicRanges; can be NULL -> but then annoDB must be provided
                                   cutoff = 0.05, # defines the cutoff for which PureCLIP sites to remove; smallest steps are in 1% (0.01); cutoff of 0.7 means that the top 30% of PureCLIP sites per gene are retained. 
                                   overlappingLoci = c("keepSingle", "removeAll", "keepAll") # Options to deal with PureCLIP sites on overlapping loci genes; use keepSingle as default; keepAll will inflate the number of binding sites by the number of overlaps, removeAll will kill all duplicated instances
) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    # Check if none is specified
    if (is.null(annoDB) & is.null(genes)) {
        msg = paste0("None of the required annotation sources annoDB or genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(annoDB) & !is.null(genes)) {
        msg = paste0("Both of the required annotation sources annoDB or genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if annoDB should be used
    if (!is.null(annoDB) & is.null(genes)) {
        stopifnot(is(annoDB, "OrganismDb"))
        if (!is.null(genes)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            genes = genes(annoDB)
        }
    }
    # Checks if genes should be used
    if (is.null(annoDB) & !is.null(genes)) {
        stopifnot(is(genes, "GenomicRanges"))
        if (!is.null(annoDB)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            genes = genes
        }
    }
    
    # handle options (remove, keep)
    overlappingLoci = match.arg(overlappingLoci, choices = c("keepSingle", "removeAll", "keepAll"))
    
    # ---
    # Store function parameters in list
    optstr = list(cutoff = cutoff, overlappingLoci = overlappingLoci,
                  annoDB = ifelse(!is.null(annoDB), "annoDB", ""), genes = ifelse(!is.null(genes),"genes", ""))
    object@params$geneWiseFilter = optstr
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    
    # Gene specific cutoffs
    # ---------------------
    # set cutpoints
    potentialCutoffs = seq(0, 1, by = 0.01)
    # group scores by gene range
    ols = findOverlaps(genes, rngInitial) %>%
        as.data.frame() %>%
        mutate(score = rngInitial$score[subjectHits]) %>%
        group_by(queryHits) 
    # calculate gene range specific quantile cutoff
    quants = ols %>%
        reframe(quantiles = quantile(score, probs = potentialCutoffs)) %>%
        group_by(queryHits) %>%
        mutate(cutoffs = potentialCutoffs) %>%
        filter(as.character(cutoffs) == as.character(cutoff))
    
    # apply cutoff to every PureCLIP peak per gene 
    peakIdToKeep = ols
    peakIdToKeep$quantiles = quants$quantiles[match(peakIdToKeep$queryHits, quants$queryHits)]
    peakIdToKeep = peakIdToKeep %>%
        filter(score >= quantiles)
    # add gene information
    filteredPerRegion = rngInitial[peakIdToKeep$subjectHits]
    names(filteredPerRegion) = gns$gene_name[peakIdToKeep$queryHits]
    
    # ---
    # Store for plotting
    dfPlot = countOverlaps(filteredPerRegion) %>% table() %>% as.data.frame() %>% rename("#N overlaps" = ".")
    object@plotData$geneWiseFilter$data = dfPlot
    
    # Deal with multiple loci
    # -----------------------
    # check for overlapping loci peaks
    potentialOverlaps = findOverlaps(filteredPerRegion)
    overlappingLociHits = potentialOverlaps[queryHits(potentialOverlaps) != subjectHits(potentialOverlaps)]
    # reduce to a single pair
    singlePair = overlappingLociHits %>%
        as.data.frame() %>%
        filter(queryHits > subjectHits)
    overlappingLociHits = overlappingLociHits[queryHits(overlappingLociHits) %in% singlePair$queryHits]
    # get numbers
    overlappingLociGenePairs = paste0(unique(names(filteredPerRegion[queryHits(overlappingLociHits)])),
                                      ",",
                                      unique(names(filteredPerRegion[subjectHits(overlappingLociHits)])), " | ")
    fractionOfRangesRemoved = paste0(round(length(overlappingLociHits) / length(rngInitial) * 100, digits = 2), "%")
    # Handle cases
    if (length(overlappingLociHits) > 0) {
        msg0 = paste0(fractionOfRangesRemoved, " (", length(overlappingLociHits),"/", length(rngInitial), ")",
                      " peaks overlap with multiple genes in the given gene annotation. \n")
        if (overlappingLoci == "keepSingle") {
            msg1 = "A single instancen of each peak is kept. This is recommended. \n "
            filteredPerRegion = filteredPerRegion[!duplicated(filteredPerRegion)]
            message(c(msg0, msg1))
        }
        if (overlappingLoci == "removeAll") {
            msg1 = "All peaks on any duplicated range are removed. This is not recommended, since many peaks are lost. \n"
            filteredPerRegion = filteredPerRegion[countOverlaps(filteredPerRegion) == 1]
            warning(c(msg0, msg1))
        }
        if (overlappingLoci == "keepAll") {
            msg1 = "Duplicated peaks due to overlapping loci are not removed. This is not recommended, since peak numbers are inflated. \n"
            warning(msg0, msg1)
        }
    }
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "pureClipGeneWiseFilter()", class = "crosslink sites",
        nIn = length(rngInitial), nOut = length(filteredPerRegion),
        per = paste0(round(length(filteredPerRegion)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", optstr$cutoff, ", overlapp option=", optstr$overlappingLoci)
    )
    object@results = rbind(object@results, resultLine)
    
    objectOut = setRanges(object, filteredPerRegion)
    return(objectOut)
}

#' @importFrom stats quantile
.selectQuantilesMultipleConditions <-
    function(covDf, userCond, userNreps, userCutoff) {
        # bind locally used variables
        per <- applyTo <- NULL
        
        # construct general cutoff matrix in 1% steps
        q = as.data.frame(apply(covDf, 2, function(x) {
            quantile(x, probs = seq(0, 1, by = 0.01))
        }))
        q$cut = seq(0, 1, by = 0.01)
        q$per = rownames(q)
        
        # select that part of q that was chosen by user
        qSel = q[q$cut %in% userCutoff, ]
        applyDf = data.frame(levels(userCond), userCutoff)
        # applyDf = data.frame(levels(userCond), userCutoff, userCond)
        
        idx = match(qSel$cut, applyDf$userCutoff)
        qSel$applyTo = applyDf$levels.userCond.[idx]
        qSel = qSel %>% pivot_longer(-c(cut, per, applyTo))
        qSel$sel = vapply(strsplit(qSel$name, "_"), `[`, 2,
                          FUN.VALUE = character(1))
        if (length(unique(userCutoff)) > 1) {
            qSel = qSel[qSel$applyTo == qSel$sel, ]
        }
        # add n.reps support to df
        nRepsDf = data.frame(defaultNreps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$defaultNreps = nRepsDf$defaultNreps[idx]
        return(qSel)
    }

#' @importFrom stats quantile
.selectQuantilesSingleCondtion <-
    function(covDf, userCond, userNreps, userCutoff) {
        # bind locally used variables
        per <- NULL
        
        # construct general cutoff matrix in 1% steps
        q = as.data.frame(apply(covDf, 2, function(x) {
            quantile(x, probs = seq(0, 1, by = 0.01))
        }))
        q$cut = seq(0, 1, by = 0.01)
        q$per = rownames(q)
        
        # select that part of q that was chosen by user
        qSel = q[q$cut %in% userCutoff, ]
        qSel = qSel %>% pivot_longer(-c(cut, per))
        qSel$sel = vapply(strsplit(qSel$name, "_"), `[`, 2,
                          FUN.VALUE = character(1))
        
        # add n.reps support to df
        nRepsDf = data.frame(defaultNreps = userNreps, applyTo = levels(userCond))
        idx = match(qSel$sel, nRepsDf$applyTo)
        qSel$defaultNreps = nRepsDf$defaultNreps[idx]
        return(qSel)
    }

reproducibilityFilter <- function(object,
                                  cutoff = NULL,
                                  n.reps = NULL,
                                  min.crosslinks = 1,
                                  returnType = c("BSFDataSet", "data.frame")) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    
    if (length(cutoff) != length(n.reps)) {
        stop("Number of values for 'cutoff' does not match the number of values
             for 'n.reps'. ")
    }
    
    metaData = getMeta(object)
    
    numberOfConditions = length(levels(metaData$condition))
    if (numberOfConditions != length(cutoff) && !is.null(cutoff)) {
        msg = paste0("Reproducibility filter cutoff does not match the number of conditions. You specified: ",
                     length(cutoff), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("The specified cutoff (",  cutoff, ") ",
                      "is applied to all conditions (",
                      paste(as.character(levels(metaData$condition)), collapse = ",") ,") \n")
        warning(paste0(msg, msg2))
        defaultCutoff = rep(cutoff, numberOfConditions)
    }
    if (is.null(cutoff)) {
        msg = paste0("Reproducibility cutoff not defined. Defaults to 0.05 for each condition. \n")
        message(msg)
        defaultCutoff = rep(0.05, numberOfConditions)
    } 
    if (numberOfConditions == length(cutoff) && !is.null(cutoff)) {
        defaultCutoff = cutoff
    }
    # Manage parameter n.reps
    if (numberOfConditions != length(n.reps) && !is.null(n.reps)) {
        msg = paste0("Parameter n.reps does not match the number of conditions. You specified: ",
                     length(n.reps), ", but there are: ", numberOfConditions, "\n")
        msg2 = paste0("n.reps defaults to N-1 for each condition. \n")
        warning(paste0(msg, msg2))
        
        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    }
    if (is.null(n.reps)) {
        msg = paste0("Parameter n.reps not defined. Defaults to N-1 for each condition. \n")
        message(paste0(msg))
        n.conditions = table(metaData$condition) %>% as.data.frame()
        defaultNreps = n.conditions$Freq -1
    } 
    if (numberOfConditions == length(n.reps) && !is.null(n.reps)) {
        defaultNreps = n.reps
    }
    
    # ---
    # Store function parameters in list
    optstr = list(cutoff = defaultCutoff, n.reps = defaultNreps, min.crosslinks = min.crosslinks)
    object@params$reproducibilityFilter = optstr
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    cond = metaData$condition
    # get number of crosslinks per binding site and replicate
    df = as.data.frame(mcols(coverageOverRanges(
        object, returnOptions = "merge_positions_keep_replicates",
        silent = TRUE)))
    
    # Manage single cutoff for single condition
    if (length(defaultCutoff) == 1) {
        if(length(levels(cond)) > 1) {
            stop("Only one cutoff is given for multiple conditions.")
        }
        # calculate sample specific thresholds
        qSel = .selectQuantilesSingleCondtion(
            covDf = df,
            userCond = cond,
            userNreps = defaultNreps,
            userCutoff = defaultCutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < min.crosslinks,
                            min.crosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))
        
        # ---
        # Store data for diagnostic plot in list
        dfPlot = df %>% 
            pivot_longer(everything()) %>%
            group_by(name, value) %>% dplyr::count() %>%
            separate(name, into = c(NA, "condition"), sep = "_", remove = FALSE)
        object@plotData$reproducibilityFilterPlot$data = dfPlot
        object@plotData$reproducibilityFilterPlot$cutoffs = qSel
        
        # calculate replicate support based on quantile cutoff
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0) 
        }) %>%
            t() %>% as.data.frame() 
        
        # --- 
        # Store results for plotting
        object@plotData$reproducibilitySamplesPlot$data = s
        
        # Filter binding sites ranges by replicate support
        support = rowSums(s) >= defaultNreps
        newRanges = getRanges(object)
        newRanges = newRanges[support]
        newObject = setRanges(object, newRanges)
    }
    
    # Manage multiple cutoffs for multiple conditions
    if (length(defaultCutoff) > 1) {
        if (length(levels(cond)) == 1) {
            stop("multiple cutoffs are given but only one condition exists")
        }
        # calculate sample specific thresholds
        qSel = .selectQuantilesMultipleConditions(
            covDf = df,
            userCond = cond,
            userNreps = defaultNreps,
            userCutoff = defaultCutoff
        )
        # apply minimal crosslink threshold
        qSel$value = ifelse(qSel$value < min.crosslinks,
                            min.crosslinks,
                            qSel$value)
        matchIdx = match(qSel$name, colnames(df))
        
        # ---
        # Store data for diagnostic plot in list
        dfPlot = df %>% 
            pivot_longer(everything()) %>%
            group_by(name, value) %>% dplyr::count() %>%
            separate(name, into = c(NA, "condition"), sep = "_", remove = FALSE)
        object@plotData$reproducibilityFilterPlot$data = dfPlot
        object@plotData$reproducibilityFilterPlot$cutoffs = qSel
        
        # calculate reproducibility per condition
        s = apply(df, 1, function(x) {
            ifelse(x > qSel$value[matchIdx], 1, 0)
        }) %>%
            t %>% as.data.frame() 
        
        # --- 
        # Store results for plotting
        object@plotData$reproducibilitySamplesPlot$data = s
        
        # Filter binding sites ranges by replicate support
        sSplit = vapply(levels(cond), function(x) {
            s %>% dplyr::select(contains(x)) %>% rowSums()
        }, FUN.VALUE = numeric(nrow(s))) %>% as.data.frame()
        idx = match(colnames(sSplit), qSel$sel)
        support = apply(sSplit, 1, function(x) {
            x >= qSel$defaultNreps[idx] 
        }) %>%
            t %>% as.data.frame()
        
        supportAll = apply(support, 1, any)
        
        newRanges = getRanges(object)
        mcols(newRanges) = support
        newRanges = newRanges[supportAll]
        newObject = setRanges(object, newRanges)
    }
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "reproducibilityFilter()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(newRanges),
        per = paste0(round(length(newRanges)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Cutoff=", paste(optstr$cutoff, collapse = ", "),
                         ", n.reps=", paste(optstr$n.reps, collapse = ", "),
                         ", min.crosslinks=", paste(optstr$min.crosslinks, collapse = ", "))
    )
    newObject@results = rbind(newObject@results, resultLine)
    
    # Manage return options
    returnType = match.arg(returnType, choices = c("BSFDataSet", "data.frame"))
    if (returnType == "BSFDataSet") {
        retObj = newObject
    }
    if (returnType == "data.frame") {
        retObj = s
    }
    return(retObj)
}


.resolveGeneOverlapsWithRule <- function(rng, ols, rule, selectID, selectName, selectType) {
    # takes a Hits object from the binding site and gene annotation matching
    # and a rule by which binding sites on multiple genes are matched 
    # and the initial ranges 
    # returns the inital binding site ranges with matched gene info
    olsInfo = data.frame(geneIndex = queryHits(ols), bsIndex = subjectHits(ols), 
                         geneID = selectID[queryHits(ols)], 
                         geneName = selectName[queryHits(ols)],
                         geneType = selectType[queryHits(ols)]) %>%
        group_by(bsIndex) %>%
        arrange(bsIndex) %>%
        mutate(choice = rule$idx[match(geneType, rule$gene_type)]) %>%
        arrange(choice, .by_group = TRUE) %>%
        slice_head(n = 1)
    
    rngAssigned = rng[olsInfo$bsIndex]
    mcols(rngAssigned) = cbind(mcols(rngAssigned), geneID = olsInfo$geneID[olsInfo$bsIndex],
                               geneName = olsInfo$geneName[olsInfo$bsIndex], geneType = olsInfo$geneType[olsInfo$bsIndex])
    # Set binding sites names and IDs
    mcols(rngAssigned)$bsName = paste0("BS", seq_along(rngAssigned), "_", rngAssigned$geneName)
    mcols(rngAssigned)$bsID = paste0("BS", seq_along(rngAssigned), "_", rngAssigned$geneID)
    return(rngAssigned)
}

assignToGenes <- function(object, # bindingSiteFinder
                          annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then genes must be provided
                          genes = NULL, # ranges of the genes of class GenomicRanges; can be NULL -> but then annoDB must be provided
                          geneMatchID = "gene_id", geneMatchName = "gene_name", geneMatchType = "gene_type", # columns that must be present in the genes object which are used for the matching
                          overlappingLoci = c("frequency", "hierarchy", "remove", "keep"), # overlapping loci solution options
                          rule = NULL # rule to apply when option hierarchy is selected
) {
    # local vairables
    datasource = NULL
    
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    # Check if none is specified
    if (is.null(annoDB) & is.null(genes)) {
        msg = paste0("None of the required annotation sources annoDB or genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(annoDB) & !is.null(genes)) {
        msg = paste0("Both of the required annotation sources annoDB or genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if annoDB should be used
    if (!is.null(annoDB) & is.null(genes)) {
        stopifnot(is(annoDB, "OrganismDb"))
        if (!is.null(genes)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "annoDB"
            # extract relevant annotation
            genes = genes(annoDB, columns=c("ENSEMBL", "GENETYPE", "SYMBOL", "GENEID"))
            # Create matching vectors for columns from input annotation
            # --------------------------------------------------------------------------
            selectID = as.character(genes$GENEID)
            selectName = as.character(genes$SYMBOL)
            selectType = as.character(genes$GENETYPE)
        }
    }
    # Checks if genes should be used
    if (is.null(annoDB) & !is.null(genes)) {
        stopifnot(is(genes, "GenomicRanges"))
        if (!is.null(annoDB)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "genes"
            # extract relevant annotation
            genes = genes
            # check correct annotation columns
            inNames = c(geneMatchID, geneMatchName, geneMatchType)
            annoColNames = colnames(mcols(genes))
            if (!all(inNames %in% annoColNames)) {
                msg = paste0("One or multiple of the provided matching columns (",
                             geneMatchID, ", ", geneMatchName, ", ", geneMatchType,
                             ") is not present in the provided annotation. \n")
                stop(msg)
            }
            # Create matching vectors for columns from input annotation
            # --------------------------------------------------------------------------
            selectID = mcols(genes)[match(geneMatchID, colnames(mcols(genes)))][[1]]
            selectName = mcols(genes)[match(geneMatchName, colnames(mcols(genes)))][[1]]
            selectType = mcols(genes)[match(geneMatchType, colnames(mcols(genes)))][[1]]
            
        }
    }
    
    # handle options (hierarchy, frequency, remove, keep)
    overlappingLoci = match.arg(overlappingLoci, choices = c("frequency", "hierarchy", "remove", "keep"))
    # Check multiple loci options
    if (overlappingLoci == "hierarchy" & is.null(rule)) {
        msg1 = paste0("Binding sites on genes with overlapping annotaitons is set to be handled by 'hierarchy', but no rule is provided. \n")
        msg2 = paste0("Change the 'overlappingLoci' option or provide a valid rule. \n")
        stop(c(msg1, msg2))
    }
    
    # ---
    # Store function parameters in list
    optstr = list(source = datasource, geneMatchID = geneMatchID,
                  geneMatchName = geneMatchName, geneMatchType = geneMatchType,
                  overlappingLoci = overlappingLoci, rule = rule)
    object@params$assignToGenes = optstr
    
    # Match Binding sites and genes
    # --------------------------------------------------------------------------
    # split ranges in duplicated event and non duplicated events
    rngInitial = getRanges(object)
    ols = findOverlaps(genes, rngInitial)  
    
    # Deal with multiple loci
    # --------------------------------------------------------------------------
    # Count out how many cases there are 
    countOlsSize = as.data.frame(table(duplicated(subjectHits(ols))))
    totalBS = countOlsSize$Freq[1]
    duplicatedBS = countOlsSize$Freq[2]
    duplicatedFraction = paste0(round(duplicatedBS / totalBS * 100, digits = 2), "%")
    
    # --- 
    # Store results for plotting
    # dfPlot = table(table(subjectHits(ols))) %>% as.data.frame() %>% rename("#N overlaps" = "Var1") %>% arrange(desc(Freq))
    dfPlot = data.frame(geneIndex = queryHits(ols), bsIndex = subjectHits(ols), 
                        geneID = selectID[queryHits(ols)], 
                        geneName = selectName[queryHits(ols)],
                        geneType = selectType[queryHits(ols)]) %>%
        mutate(value = 1) %>% 
        select(-geneIndex, -geneID, -geneName) %>% 
        pivot_wider(names_from = geneType, values_from = value, values_fn = length, values_fill = 0) %>%
        select(-bsIndex) %>%
        mutate_all(~ ifelse(. > 1, 1, .))
    object@plotData$assignToGenes$dataOverlaps = dfPlot
    
    # Case where there are no overlaps
    if (duplicatedBS == 0) {
        msg = paste0("No binding sites (", format(totalBS, big.mark = ',', decimal.mark = "."), ") on overlapping gene loci. \n")
        msg2 = paste0("Parameter 'overlappingLoci' set to '", overlappingLoci, "' not in use. \n")
        message(c(msg, msg2))
    }
    # Cases where there are multiple overlaps
    if (length(overlappingLoci) > 0) {
        msg0 = paste0(duplicatedFraction, " (", duplicatedBS,"/", totalBS, ")",
                      " of binding sites overlap with multiple genes in the given gene annotation. \n")
        if (overlappingLoci == "hierarchy") {
            # apply hierarchy solution 
            msg1 = paste0("Apply 'hierarchy' solution")
            message(c(msg0, msg1))
            ruleMod = data.frame(gene_type = rule, idx = seq_along(rule))
            rngResolved = .resolveGeneOverlapsWithRule(rng = rngInitial, ols = ols, rule = ruleMod,
                                                       selectID = selectID, selectName = selectName, selectType = selectType)
            if (length(rngInitial) != length(rngResolved)) {
                # warning(paste0("Hierarchy - Number of binding sites before and after gene assignemnt not equal!: Before=", length(rngInitial), ", After=", length(rngResolved) ))
                message(paste0(format(length(rngInitial)-length(rngResolved), big.mark = ",", decimal.mark = "."), " binding sites could not be assigned to a gene range.\n"))
            }
            
        }
        if (overlappingLoci == "frequency") {
            # apply frequency solution
            msg1 = paste0("Apply 'frequency' solution")
            message(c(msg0, msg1))
            ruleFreq = selectType %>% table() %>% as.data.frame() %>% rename(gene_type = 1) %>% arrange(desc(Freq)) %>% mutate(idx = row_number())  %>% rename(n = 2)
            rngResolved = .resolveGeneOverlapsWithRule(rng = rngInitial, ols = ols, rule = ruleFreq,
                                                       selectID = selectID, selectName = selectName, selectType = selectType)
            if (length(rngInitial) != length(rngResolved)) {
                # warning(paste0("Frequency - Number of binding sites before and after gene assignemnt not equal!: Before=",length(rngInitial), ", After=", length(rngResolved) ))
                message(paste0(format(length(rngInitial)-length(rngResolved), big.mark = ",", decimal.mark = "."), " binding sites could not be assigned to a gene range.\n"))
            }
            
        }
        if (overlappingLoci == "remove") {
            # apply remove option
            msg1 = "Binding sites from overlapping loci are removed. This is not recommended. Please see options 'heirarchy' and 'frequency'. \n"
            warning(c(msg0, msg1))
            rngResolved = rngInitial[subjectHits(ols)]
            mcols(rngResolved) = cbind(mcols(rngResolved), geneID = selectID[queryHits(ols)],
                                       geneType = selectType[queryHits(ols)], geneName = selectName[queryHits(ols)] )
            rngResolved = rngResolved[countOverlaps(rngResolved) == 1]
        }
        if (overlappingLoci == "keep") {
            msg1 = "Binding sites from overlapping loci are kept. This is not recommended. Please see options 'heirarchy' and 'frequency'. \n"
            warning(c(msg0, msg1))
            rngResolved = rngInitial[subjectHits(ols)]
            mcols(rngResolved) = cbind(mcols(rngResolved), geneID = selectID[queryHits(ols)],
                                       geneType = selectType[queryHits(ols)], geneName = selectName[queryHits(ols)] )
        }
    }
    
    # --- 
    # Store results for plotting
    dfPlot = rngResolved$geneType %>% table() %>% as.data.frame() %>% rename("GeneType" = ".") %>% arrange(desc(Freq))
    object@plotData$assignToGenes$dataSpectrum = dfPlot
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "assignToGenes()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rngResolved),
        per = paste0(round(length(rngResolved)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Overlaps=", optstr$overlappingLoci, ", Source=", optstr$source, 
                         ifelse(!rlang::is_empty(optstr$rule), paste0(", rule=", paste(optstr$rule, collapse = ">")), ""))
    )
    object@results = rbind(object@results, resultLine)
    
    objectOut = setRanges(object, rngResolved)
    return(objectOut)
}

assignToTranscriptRegions <- function(object, # bindingSiteFinder
                                      annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then transcriptRegionList must be provided
                                      transcriptRegionList = NULL, # GRangesList object with the transcript features to be used; can be NULL -> but then annoDB must be provided
                                      overlappingLoci = c("frequency", "hierarchy", "flag", "remove"), # overlapping loci solution options
                                      rule = NULL # rule to apply when option hierarchy is selected
) {
    # local variables
    datasource = NULL
    
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    # Check if none is specified
    if (is.null(annoDB) & is.null(transcriptRegionList)) {
        msg = paste0("None of the required annotation sources annoDB or transcriptRegionList was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(annoDB) & !is.null(transcriptRegionList)) {
        msg = paste0("Both of the required annotation sources annoDB or transcriptRegionList are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # AnnoDB should be used
    if (!is.null(annoDB) & is.null(transcriptRegionList)) {
        stopifnot(is(annoDB, "OrganismDb"))
        if (!is.null(transcriptRegionList)) {
            msg = paste0("Parameter annoDB and transcriptRegionList are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "annoDB"
            # extract relevant annotation
            cdseq = cds(annoDB) 
            intrns = unlist(intronsByTranscript(annoDB)) 
            utrs3 = unlist(threeUTRsByTranscript(annoDB)) 
            utrs5 = unlist(fiveUTRsByTranscript(annoDB)) 
            transcriptRegionList = GRangesList(CDS = cdseq, Intron = intrns, UTR3 = utrs3, UTR5 = utrs5)
            names(transcriptRegionList) = toupper(names(transcriptRegionList))
        }
    }
    # TranscriptRegionList should be used
    if (is.null(annoDB) & !is.null(transcriptRegionList)) {
        stopifnot(is(transcriptRegionList, "GenomicRangesList"))
        if (!is.null(annoDB)) {
            msg = paste0("Parameter annoDB and transcriptRegionList are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            datasource = "transcriptRegionList"
            # extract relevant annotation
            transcriptRegionList = transcriptRegionList
            names(transcriptRegionList) = toupper(names(transcriptRegionList))
            # 
        }
    }
    # handle options (hierarchy, frequency, remove, keep)
    overlappingLoci = match.arg(overlappingLoci, choices = c("frequency", "hierarchy", "flag", "remove"))
    # Check multiple loci options
    if (overlappingLoci == "hierarchy" & is.null(rule)) {
        msg1 = paste0("Multiple region problem is set to be handled by 'hierarchy', but no rule is provided. \n")
        msg2 = paste0("Change the 'overlappingLoci' option or provide a valid rule. \n")
        stop(c(msg1, msg2))
    }
    rule = toupper(rule)
    if (! all(rule %in% names(transcriptRegionList))){
        msg = paste0("Regions defined in rule (", paste(rule, collapse = " > "),
                     ") does not match the input region names (", paste(names(transcriptRegionList), collapse = ", "), ")")
        stop(msg)
    }
    if (overlappingLoci == "hierarchy" & is.null(rule)) {
        msg1 = paste0("Overlap are set to be solved by 'hierarchy', but no rule is provided. \n")
        msg2 = paste0("Change the 'overlappingLoci' option or provide a valid rule. \n")
        stop(c(msg1, msg2))
    } 
    
    # ---
    # Store function parameters in list
    optstr = list(source = datasource, overlappingLoci = overlappingLoci, rule = rule)
    object@params$assignToTranscriptRegions = optstr
    
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    # Match transcript annotation with binding sites
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    rng = getRanges(object)
    countDf = lapply(transcriptRegionList, function(x) {
        countOverlaps(rng, x)
    }) %>% as.data.frame() %>% rename_with(toupper)
    # Add intergenic cases
    countDf = countDf %>% mutate(INTERGENIC = ifelse(rowSums(.) == 0, 1, 0))
    
    # --- 
    # Store results for plotting
    plotDf = countDf
    plotDf[plotDf > 1] = 1
    object@plotData$assignToTranscriptRegions$dataOverlaps = plotDf
    
    # Deal with multiple loci
    # --------------------------------------------------------------------------
    # Count out how many cases there are 
    totalBs = nrow(countDf)
    numberOfLociHit = rowSums(countDf != 0)
    ambigousBs = nrow(countDf[numberOfLociHit != 1,])
    ambigousBsFraction = paste0(round(ambigousBs / totalBs * 100, digits = 2), "%")
    
    # Case where there are no overlaps
    if (ambigousBs == 0) {
        msg0 = paste0("No binding sites (", format(totalBs, big.mark = ',', decimal.mark = "."), ") on overlapping transcript loci. \n")
        msg1 = paste0("Parameter 'overlappingLoci' set to '", overlappingLoci, "' not in use. \n")
        message(c(msg0, msg1))
    }
    # Cases where we have to resolve overlapping loci
    if (ambigousBs > 0) {
        msg0 = paste0(ambigousBsFraction, " (", ambigousBs,"/", totalBs, ")",
                      " of binding sites overlap with multiple different transcript regions in the given annotation. \n")
        if (overlappingLoci == "hierarchy") {
            # apply hierarchy solution
            msg1 = paste0("Apply 'hierarchy' solution. \n ")
            message(c(msg0, msg1))
            countDfSorted = countDf[, rule] 
            countDfSorted[countDfSorted > 0] = 1
            transcriptRegion = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")]
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlappingLoci == "frequency") {
            msg1 = paste0("Apply 'frequency' solution. \n ")
            message(c(msg0, msg1))
            countDfSorted = countDf
            # get rule by frequency
            rule = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")] %>% table() %>% as.data.frame() %>% rename(transcript_type = 1) %>%
                arrange(desc(Freq)) %>% pull(transcript_type) %>% as.character()
            # apply frequency solution
            countDfSorted = countDfSorted[, rule]
            transcriptRegion = names(countDfSorted)[max.col(countDfSorted, ties.method = "first")]
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlappingLoci == "flag") {
            msg1 = paste0("Binding sites with multiple different overlapping transcript regions are being flaged with the tag 'Ambiguous'. \n")
            message(c(msg0, msg1))
            countDfSorted = countDf
            countDfSorted[countDfSorted > 0] = 1
            # flag multiple loci binding sites
            col_names = colnames(countDfSorted)
            transcriptRegion = ifelse(rowSums(countDfSorted) == 1, col_names[max.col(countDfSorted, ties.method = "first")], "Ambiguous")
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            # remove mulitple loci binding sitges
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
        }
        if (overlappingLoci == "remove") {
            msg1 = paste0("Binding sites with multiple different overlapping transcript regions are being removed. This is not recommended! \n")
            warning(c(msg0, msg1))
            countDfSorted = countDf
            countDfSorted[countDfSorted > 0] = 1
            # flag multiple loci binding sites
            col_names = colnames(countDfSorted)
            transcriptRegion = ifelse(rowSums(countDfSorted) == 1, col_names[max.col(countDfSorted, ties.method = "first")], "Ambiguous")
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            # remove mulitple loci binding sitges
            mcols(rng) = cbind.data.frame(mcols(rng), transcriptRegion = transcriptRegion)
            rng = rng[rng$transcriptRegion != "Ambiguous"]
        }
    }
    # --- 
    # Store results for plotting
    dfPlot = rng$transcriptRegion %>% table() %>% as.data.frame() %>% rename("TranscriptRegion" = ".") %>% arrange(desc(Freq))
    object@plotData$assignToTranscriptRegions$dataSpectrum = dfPlot
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "assignToTranscriptRegions()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("Overlaps=", optstr$overlappingLoci, ", Source=", datasource, 
                         ifelse(!rlang::is_empty(optstr$rule), paste0(", rule=", paste(optstr$rule, collapse = ">")), ""))
    )
    object@results = rbind(object@results, resultLine)
    
    outObject = setRanges(object, rng)
    return(outObject)
}

annotateWithScore <- function(object, # bindingSiteFinder
                              scoreRanges, # GRanges object with scores in the meta columns that should be matched 
                              MatchColScore = "score", # column name for matching
                              bsMatchOption = c("max", "sum", "mean")
) {
    
    # bind locally used variables
    qHits <- NULL
    
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    stopifnot(is(scoreRanges, "GRanges"))
    
    # handle options
    bsMatchOption = match.arg(bsMatchOption, choices = c("max", "sum", "mean"))
    
    # match score column
    scoreColNames = colnames(mcols(scoreRanges))
    if (!(MatchColScore %in% scoreColNames)) {
        msg = paste0("Matching columns (", MatchColScore,
                     ") is not present in the provided scoreRanges. \n")
        stop(msg)
    }
    
    # ---
    # Store function parameters in list
    optstr = list(MatchColScore = MatchColScore, bsMatchOption = bsMatchOption)
    object@params$annotateWithScore = optstr
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    rngInitial = getRanges(object)
    rng = getRanges(object)
    ol = findOverlaps(rng, scoreRanges)
    
    if (length(ol) == 0) {
        stop("Ranges of 'object' and 'scoreRanges' do not overlap.")
    }
    
    matchDF = data.frame(qHits = queryHits(ol),
                         sHits = subjectHits(ol),
                         score = scoreRanges$score[subjectHits(ol)])
    
    if (bsMatchOption == "max") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = max(score), .groups = "drop") %>%
            as.data.frame()
    }
    if (bsMatchOption == "sum") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = sum(score), .groups = "drop") %>%
            as.data.frame()
    }
    if (bsMatchOption == "mean") {
        score = dplyr::group_by(matchDF, qHits) %>%
            dplyr::summarize(score = mean(score), .groups = "drop") %>%
            as.data.frame()
    }
    
    # --- 
    # Store results for plotting
    object@plotData$annotateWithScore$data = score$score
    
    # ---
    # Store for results
    resultLine = data.frame(
        funName = "annotateWithScore()", class = "binding sites",
        nIn = length(rngInitial), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rngInitial), digits = 2)*100,"%"),
        options = paste0("MatchOption=", optstr$bsMatchOption, ", MatchColScore=", MatchColScore)
    )
    object@results = rbind(object@results, resultLine)
    
    mcols(rng)$score = score$score
    newObject = setRanges(object, rng)
    return(newObject)
}



# ------------------------------------------------------------------------------
# for estimate binding site width option


.mergeCrosslinkSites <- function(rng, # a GRanges object; -> holds PureCLIP sites, is single nt size
                                 sgn, # the crosslink signal merged per replicates
                                 bsSize, # the binding site size to compute
                                 minWidth, # the minimal width sites should be retained after intital merge
                                 computeOption = c("full", "simple")
) {
    
    # handle compute options
    computeOption = match.arg(computeOption, choices = c("full", "simple"))
    
    # summarize signal over all replicates for mergeing
    sgnMergePlus = sgn$signalPlus
    sgnMergeMinus = sgn$signalMinus
    
    rngS1 = rng
    
    ### Merge peaks for given bs size
    rngS2 = reduce(rngS1, min.gapwidth = bsSize - 1)
    
    ### Keep only regions that are larger or equal to minWidth 
    # -> if minWiidth == 3, then the smallest range to consider is 3
    rngS3 = rngS2[width(rngS2) >= minWidth]
    names(rngS3) = seq_along(rngS3)
    
    ### Center detection and extension
    rngCenterPlus <- GRanges()
    rngCenterMinus <- GRanges()
    rngToProcessPlus <- subset(rngS3, strand == "+")
    rngToProcessMinus <- subset(rngS3, strand == "-")
    
    Counter = 0
    while (TRUE) {
        # quit if no more regions to check
        if (length(rngToProcessMinus) == 0 &
            length(rngToProcessPlus) == 0) {
            break
        } else {
            if (length(rngToProcessPlus) != 0) {
                # get max xlink position of each peak
                peaksMaxPosPlus = as.matrix(sgnMergePlus[rngToProcessPlus])
                peaksMaxPosPlus[is.na(peaksMaxPosPlus)] = -Inf
                peaksMaxPosPlus = max.col(peaksMaxPosPlus,
                                          ties.method = "first")
                
                # make new peaks centered arround max position
                currentPeaksPlus = rngToProcessPlus
                start(currentPeaksPlus) =
                    start(currentPeaksPlus) + peaksMaxPosPlus -1
                end(currentPeaksPlus) = start(currentPeaksPlus)
                currentPeaksPlus = currentPeaksPlus + ((bsSize - 1) / 2)
                # store peaks
                rngCenterPlus = c(rngCenterPlus, currentPeaksPlus)
                # remove peak regions from rest of possible regions
                currentPeaksPlus = as(currentPeaksPlus + ((bsSize - 1) / 2),
                                      "GRangesList")
                
                # update peak regions that are left for processing
                rngToProcessPlus = unlist(psetdiff(rngToProcessPlus,
                                                   currentPeaksPlus))
            }
            if (length(rngToProcessMinus) != 0) {
                peaksMaxPosMinus = as.matrix(sgnMergeMinus[rngToProcessMinus])
                peaksMaxPosMinus[is.na(peaksMaxPosMinus)] = -Inf
                peaksMaxPosMinus = max.col(peaksMaxPosMinus,
                                           ties.method = "last")
                
                currentPeaksMinus = rngToProcessMinus
                start(currentPeaksMinus) =
                    start(currentPeaksMinus) + peaksMaxPosMinus -1
                end(currentPeaksMinus) = start(currentPeaksMinus)
                currentPeaksMinus = currentPeaksMinus + ((bsSize - 1) / 2)
                
                rngCenterMinus = c(rngCenterMinus, currentPeaksMinus)
                
                currentPeaksMinus =
                    as(currentPeaksMinus + ((bsSize - 1) /2),
                       "GRangesList")
                
                rngToProcessMinus = unlist(psetdiff(rngToProcessMinus,
                                                    currentPeaksMinus))
            }
            Counter = Counter + 1
        }
        # compute option simple exits the loop after the first iteration, which is 
        # when the first round of merge and extend is done
        if (computeOption == "simple") {
            if (length(rngCenterPlus) > 0 | length(rngCenterMinus) > 0){
                break
            }
        }
    }
    rngS4 = c(rngCenterPlus, rngCenterMinus)
    rngS4 = .sortRanges(rngS4)
    return(rngS4)
}
.sortRanges <- function(rng) {
    rngSort = GenomeInfoDb::sortSeqlevels(rng)
    rngSort = sort(rngSort)
    return(rngSort)
}
.subsetByChr <- function(object, chr) {
    # subset ranges
    rng = getRanges(object)
    rngSub = rng[seqnames(rng) == chr,]
    
    # subset the signal
    sgn = getSignal(object)
    sgnSub = lapply(sgn, function(selStrand) {
        lapply(selStrand, function(chrList) {
            chrList[names(chrList) == chr]
        })
    })
    
    objectNew = setRanges(object, rngSub)
    objectNew = setSignal(objectNew, sgnSub)
    return(objectNew)
}
.reduceSignalToFrame <- function(object, frame) {
    rng = getRanges(object) 
    # extend ranges by desired frame on which signal should be kept
    rng = rng + frame
    newObj = setRanges(object, rng)
    # drop signal outside of protective frame
    newObj = newObj[seq_along(rng), drop=TRUE]
    # put original sized ranges back in place
    rng = rng-frame
    newObj = setRanges(newObj, rng)
    
    return(newObj)
}
.collapseSamples <- function(signal) {
    pSum <- Reduce(`+`, signal$signalPlus)
    names(pSum) <- names(signal$signalPlus[[1]])
    mSum <- Reduce(`+`, signal$signalMinus)
    names(mSum) <- names(signal$signalMinus[[1]])
    
    mergedSignal <- list(signalPlus = pSum, signalMinus = mSum)
    return(mergedSignal)
}
.approximateBindingSites_medium <- function(rng, sgn, bsSize, minWidth) {
    # approximate binding sites by first iteration of merging
    rng = .mergeCrosslinkSites(
        rng = rng,
        sgn = sgn,
        bsSize = bsSize,
        minWidth = 3,
        computeOption = "simple"
    )
    rng$bsSize = bsSize
    return(rng)
}
.approximateBindingSites_coarse <- function(rng, bsSize, minWidth) {
    ### Merge peaks for given bs size
    # rng = reduce(rng, min.gapwidth = bsSize - 1)
    rng = reduce(rng)
    ### Remove peaks smaller than min width
    rng = rng[width(rng) >= minWidth]
    ### Keep peaks that have the exact size
    # rng = rng[width(rng) == bsSize]
    names(rng) = seq_along(rng)
    ### Approximate binding sites by center of merged region
    rng = resize(rng, fix = "center", width = bsSize)
    rng$bsSize = bsSize
    return(rng)
}

estimateBsWidth <- function(object, # BindingSiteFinder object
                            bsFilterSteps = NULL, # the percentage filter steps to test -> for the pureClipGeneWiseFilter()
                            est.geneFilter_resolution = c("medium", "coarse", "fine", "finest"), # the desired resolution of the bsFilter, default = medium
                            bsWidthSteps = NULL, # the steps to test -> for makeBindingSites()
                            est.maxSites = Inf, 
                            max.bs = 13, # the largest binding site width which should considered in the testing; defaults to 13nt
                            est.bs_resolution = c("medium", "fine", "coarse"),
                            est.subChr = "chr1", # define on which chromosome the training should be done 
                            annoDB = NULL, # annotation data base of class OrganismDbi; can be NULL -> but then genes must be provided
                            genes = NULL # ranges of the genes of class GenomicRanges; can be NULL -> but then annoDB must be provided,
) {
    # INPUT CHECKS
    # --------------------------------------------------------------------------
    stopifnot(is(object, "BSFDataSet"))
    
    # Check if none is specified
    if (is.null(annoDB) & is.null(genes)) {
        msg = paste0("None of the required annotation sources annoDB or genes was specified. ")
        stop(msg)
    }
    # Check if both are specified
    if (!is.null(annoDB) & !is.null(genes)) {
        msg = paste0("Both of the required annotation sources annoDB or genes are specified. Please provide only one of the two. ")
        stop(msg)
    }
    # Checks if annoDB should be used
    if (!is.null(annoDB) & is.null(genes)) {
        stopifnot(is(annoDB, "OrganismDb"))
        if (!is.null(genes)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            genes = genes(annoDB)
        }
    }
    # Checks if genes should be used
    if (is.null(annoDB) & !is.null(genes)) {
        stopifnot(is(genes, "GenomicRanges"))
        if (!is.null(annoDB)) {
            msg = paste0("Parameter annoDB and genes are specified at the same time. Use only one of them.")
            stop(msg)
        } else {
            # extract relevant annotation
            genes = genes
        }
    }
    
    # handle gene filter resolution
    est.geneFilter_resolution = match.arg(est.geneFilter_resolution, choices = c("medium", "coarse", "fine", "finest"))
    if (!is.null(bsFilterSteps)) {
        bsFilterSteps = bsFilterSteps
    } else {
        if (est.geneFilter_resolution == "coarse") {
            bsFilterSteps = seq(from = 0, to = 1, by = 0.2)
            bsFilterSteps = bsFilterSteps[1:length(bsFilterSteps)-1]
        }
        if (est.geneFilter_resolution == "medium") {
            bsFilterSteps = seq(from = 0, to = 1, by = 0.1)
            bsFilterSteps = bsFilterSteps[1:length(bsFilterSteps)-1]
        }
        if (est.geneFilter_resolution == "fine") {
            bsFilterSteps = seq(from = 0, to = 1, by = 0.05)
            bsFilterSteps = bsFilterSteps[1:length(bsFilterSteps)-1]
        }
        if (est.geneFilter_resolution == "finest") {
            bsFilterSteps = seq(from = 0, to = 1, by = 0.01)
            bsFilterSteps = bsFilterSteps[1:length(bsFilterSteps)-1]
        }
    }
    
    # handle binding site compute resolution
    est.bs_resolution = match.arg(est.bs_resolution, choices = c("medium", "fine", "coarse"))
    
    # handle max.bs 
    if (!is.null(bsWidthSteps)) {
        bsWidthSteps = bsWidthSteps
    } else {
        bsWidthSteps = seq(from = 3, to = max.bs, by = 2)
    }
    
    # ---
    # Store function parameters in list
    optstr = list(est.bs_resolution = est.bs_resolution,
                  est.geneFilter_resolution = est.geneFilter_resolution, 
                  max.bs = max.bs, est.maxSites = est.maxSites,
                  est.subChr = est.subChr)
    object@params$estimateBsWidth = optstr
    
    # PREPARE TEST RANGES + SIGNAL
    # --------------------------------------------------------------------------

    suppressWarnings({
        checkRng = getRanges(object)
        # check if subset is part of the seqnames from ranges
        if (!est.subChr %in% levels(seqnames(checkRng))){
            msg = paste0("Chromosome to estimate on (", est.subChr, "), is not included in the ranges: ", paste(levels(seqnames(checkRng)), collapse = ", "))
            stop(msg)
        }
        # limit estimation to a specific chromosome
        if (!is.null(est.subChr)) {
            redObj = .subsetByChr(object, chr = est.subChr)
        } else {
            redObj = object
        }
        # limit estimation to a maximum number of sites
        estRng = getRanges(redObj)
        if (length(estRng) > est.maxSites) {
            estRng = head(estRng, est.maxSites)
            redObj = setRanges(redObj, estRng)
        }
        # limit the clip signal to the ranges used for estimation (plus some extra offset)
        maxFrame = ceiling((max(bsWidthSteps) *3))
        redObj = .reduceSignalToFrame(redObj, frame = maxFrame)
        # collapse signal from replicates 
        sgnMerge = .collapseSamples(getSignal(redObj))
    })
    
    # MAIN COMPUTE
    # --------------------------------------------------------------------------
    
    # a counter that shows to how many percent the computation is done
    counterTotalIterations = length(bsFilterSteps) * length(bsWidthSteps)
    counterCurrentIterations = 0
    
    # print start message 
    msg = paste0("Estimation at: ", round((counterCurrentIterations/counterTotalIterations)*100 ), "%"  )
    print(msg)
    suppressWarnings({
        suppressMessages({
            # calculate binding sites for each filter step and width
            scoreAllDf = lapply(bsFilterSteps, function(bsFilterStep){
                # apply current gene-wise filter
                currFilterObj = pureClipGeneWiseFilter(object = redObj, genes = genes, cutoff = bsFilterStep)
                currRng = getRanges(currFilterObj)
                
                # calculate binding sites for current bsWidth
                rngPerWidth = lapply(bsWidthSteps, function(bsWidthStep){
                    if (est.bs_resolution == "fine") {
                        # calculate with full binding sites
                        currBsObj = makeBindingSites(object = currFilterObj, bsSize = bsWidthStep)
                        currBs = getRanges(currBsObj)
                    }
                    if (est.bs_resolution == "medium") {
                        # approximate binding sites by a single merge and extend round
                        currBs = .approximateBindingSites_medium(rng = currRng, sgn = sgnMerge, bsSize = bsWidthStep, minWidth = 3)
                    }
                    if (est.bs_resolution == "coarse") {
                        # approximate binding sites by reduced pureclip sites center
                        currBs = .approximateBindingSites_coarse(rng = currRng, bsSize = bsWidthStep, minWidth = 3)
                    }
                    return(currBs)
                })
                rngPerWidth = unlist(GRangesList(rngPerWidth))
                
                # handle all plus ranges
                cRangePlus = subset(rngPerWidth, strand == "+")
                if (length(cRangePlus) > 0) {
                    bsSumPlus = sum(sgnMerge$signalPlus[cRangePlus])
                    extendedRangePlus = cRangePlus + cRangePlus$bsSize
                    exSumPlus = sum(sgnMerge$signalPlus[extendedRangePlus])
                    extendedRangePlus$score = (bsSumPlus) / ((exSumPlus - bsSumPlus + 0.01) / 2)
                } else {
                    extendedRangePlus = cRangePlus
                }
                
                # handle all minus ranges
                cRangeMinus = subset(rngPerWidth, strand == "-")
                if (length(cRangeMinus) > 0) {
                    bsSumMinus = sum(sgnMerge$signalMinus[cRangeMinus])
                    extendedRangeMinus = cRangeMinus + cRangeMinus$bsSize
                    exSumMinus = sum(sgnMerge$signalMinus[extendedRangeMinus])
                    extendedRangeMinus$score = (bsSumMinus) / ((exSumMinus - bsSumMinus + 0.01) / 2)
                } else {
                    extendedRangeMinus = cRangeMinus
                }
                
                
                # combine both and make results dataframe
                df = rbind(as.data.frame(mcols(extendedRangeMinus)), as.data.frame(mcols(extendedRangePlus))) %>%
                    group_by(bsSize) %>%
                    summarize(score = median(score), .groups = "keep") %>%
                    mutate(geneWiseFilter = bsFilterStep)
                
                # update chunk counter
                counterCurrentIterations <<- counterCurrentIterations + length(bsWidthSteps)
                msg = paste0("Estimation at: ", round((counterCurrentIterations/counterTotalIterations)*100 ), "%"  )
                print(msg)
                return(df)
            })    
        })
    })
    
    df = do.call("rbind", scoreAllDf)
    
    # --- 
    # Store results for plotting
    object@plotData$estimateBsWidth$data = df
    
    # Estimate binding site width and gene filter based on mean over all iterations
    dfMean = df %>%
        group_by(bsSize) %>%
        summarise(ms = mean(score), sd = sd(score)) %>%
        mutate(geneWiseFilter = "mean") 
    
    # bsSize is estimated by the highest mean score with the lowest standard deviation
    est.bsSize = dfMean %>% 
        slice_max(ms, n = 1, with_ties = TRUE) %>%
        slice_min(sd, n = 1, with_ties = FALSE) %>%
        pull(bsSize)
    
    # GeneFilter is estimated from the selected bsSize
    est.geneFilter = df %>%
        filter(bsSize == est.bsSize) %>%
        filter(score >= max(dfMean$ms[dfMean$bsSize == est.bsSize])) %>%
        head(1) %>%
        pull(geneWiseFilter)
    
    object@params$bsSize = est.bsSize
    object@params$geneFilter = est.geneFilter
    
    # ---
    # Store for results
    rng = getRanges(object)
    resultLine = data.frame(
        funName = "estimateBsWidth()", class = "estimate",
        nIn = length(rng), nOut = length(rng),
        per = paste0(round(length(rng)/ length(rng), digits = 2)*100,"%"),
        options = paste0("est.bs_resolution=", est.bs_resolution, ", est.geneFilter_resolution=", est.geneFilter_resolution,
                         ", max.bs=", max.bs, ", est.maxSites=", est.maxSites, ", est.subChr=", est.subChr)
    )
    object@results = rbind(object@results, resultLine)
    
    return(object)
}
