###############################
#                             #
#Ensembl specific functions   #
###############################

checkWrapperArgs = function(id, type, mart){
    if(missing(type))
        stop("Specify the type of identifier you are using, see ?getGene for details. Valid values for the type argument can be found with the listFilters function.") 
    if(!type %in% listFilters(mart)[,1])
        stop(paste("Invalid identifier type:",type," see ?getGene for details. Use the listFilters function to get the valid value for the type argument.", sep=""))   
    if(missing(id))
        stop("No identifiers specified.  Use the id argument to specify a vector of identifiers for which you want to retrieve the annotation.")
}

getGene <- function( id, type, mart){
    martCheck(mart,"ensembl") 
    checkWrapperArgs(id, type, mart)
    symbolAttrib = switch(strsplit(martDataset(mart), "_", fixed = TRUE, useBytes = TRUE)[[1]][1],
                          hsapiens = "hgnc_symbol",
                          mmusculus = "mgi_symbol",
                          "external_gene_id")
    typeAttrib = switch(type,affy_hg_u133a_2 = "affy_hg_u133a_v2",type)
    attrib = c(typeAttrib, symbolAttrib,
               "description", "chromosome_name", "band", "strand",
               "start_position", "end_position", "ensembl_gene_id")
    table = getBM(attributes = attrib, filters = type, values = id, mart=mart)
    return(table)
}

.checkSequenceArgs <- function(seqType, type, id, chromosome, upstream, downstream) {
    
    validSeqTypes <- c("cdna","peptide","3utr","5utr", "gene_exon", "transcript_exon",
                    "transcript_exon_intron", "gene_exon_intron","coding","coding_transcript_flank",
                    "coding_gene_flank","transcript_flank","gene_flank")
    if(missing(seqType) || !seqType %in% validSeqTypes) {
        stop("Please specify the type of sequence that needs to be retrieved when using biomaRt in web service mode\n.",
             "Valid options are: ", paste(validSeqTypes, collapse = ', '),
             call. = FALSE)
    }
    
    if(missing(type))
        stop("Please specify the type argument.\n",
            "If you use chromosomal coordinates to retrieve sequences ",
            "then the type argument will specify the type of gene indentifiers that you will retrieve with the sequences.\n", 
            "If you use a vector of identifiers to retrieve the sequences ",
            "the type argument specifies the type of identifiers you are using.")
    
    ## must use one and only one of 'id' and 'chromosome'
    if(missing(id) && missing(chromosome))
        stop("You must provide either the 'id' or 'chromosome' argument.",
             call. = FALSE)
    if(!missing(chromosome) && !missing(id))
        stop("You must provide only one of the 'id' and 'chromosome' arguments.",
             call. = FALSE)
    
    if(grepl(pattern = "flank", x = seqType) && 
       (missing(upstream) && missing(downstream))) {
      stop("You must provide either the 'upstream' or 'downstream' ",
           "argument when requesting flanking sequences.", call. = FALSE)
    }
    
}

.getSequenceFromCoords <- function(chromosome, start, end, type, seqType, 
                                   upstream, downstream, mart, 
                                   useCache = TRUE, verbose=FALSE) {
    
    if(missing(start) || missing(end))
        stop("You must specify both a start and end position.")
    
    start = as.integer(start)
    end = as.integer(end)
    
    if(!missing(upstream) && !missing(downstream)) {
        stop("getSequence() only allows specifying either the 'upstream' or 'downstream' argument but not both.")
    }
    
    filters <- list("chromosome_name" = chromosome,
                    "start" = start,
                    "end" = end)
    
    if(!missing(upstream)) {
        filters[[ "upstream_flank" ]] <- upstream
    } else if(!missing(downstream)) {
        filters[[ "downstream_flank" ]] <- downstream
    } 
    
    sequence <- getBM(attributes = c(seqType, type), 
                      filters = filters, 
                      mart = mart, 
                      checkFilters = FALSE, 
                      verbose = verbose,
                      useCache = useCache)
    return(sequence)
}

.getSequenceFromId <- function(id, type, seqType, 
                               upstream, downstream, mart, 
                               useCache = TRUE, verbose=FALSE) {
    
    if(missing(type)) {
        stop("Type argument is missing. ",
             "This will be used to retrieve an identifier along with the sequence so one knows which gene it is from. ", 
             "Use the listFilters() function to select a valid type argument.")
    }
    if(!type %in% listFilters(mart, what="name")) {
        stop("Invalid 'type' argument.  Use the listFilters() function to select a valid type argument.")
    }
    
    if(missing(upstream) && missing(downstream)){
        filters <- list(id)
        names(filters) <- type
    } else if(!missing(upstream) && missing(downstream)) {
        filters <- list(id, upstream)
        names(filters) <- c(type, "upstream_flank")
    } else if (!missing(downstream) && missing(upstream)) {
        filters <- list(id, downstream)
        names(filters) <- c(type, "downstream_flank")
    } else {
        stop("Currently getSequence only allows the user to specify either an upstream of a downstream argument but not both.")
    }
    
    if(!type %in% listAttributes(mart, page = "sequences", what = "name")) {
        mapping_id <- getBM(attributes = c(type, "ensembl_gene_id"), 
                               filters = type, 
                               values = id,
                               mart = mart, 
                               useCache = useCache)
        
        filters[[1]] <- mapping_id$ensembl_gene_id
        names(filters)[1] <- "ensembl_gene_id"
        
        mapping_seq <- getBM(attributes = c(seqType, "ensembl_gene_id"), 
                          filters = filters,
                          mart = mart, 
                          checkFilters = FALSE,
                          verbose = verbose,
                          useCache = useCache)
        
        ## merge data.frames and keep rows for any id that doesn't have a match
        sequence <- merge(mapping_seq, mapping_id, 
                          by = "ensembl_gene_id", all.y = TRUE)
        sequence <- sequence[, !(names(sequence) %in% "ensembl_gene_id")]
        
    } else {
        sequence <- getBM(attributes = c(seqType, type), 
                          filters = filters,
                          mart = mart, 
                          checkFilters = FALSE,
                          verbose = verbose,
                          useCache = useCache)
    }
    return(sequence)    
}

getSequence <- function(chromosome, start, end, id, type, seqType, 
                        upstream, downstream, mart, 
                        useCache = TRUE, verbose = FALSE) {
    
    martCheck(mart,c("ensembl","ENSEMBL_MART_ENSEMBL"))
    
    .checkSequenceArgs(seqType, type, chromosome, id, upstream, downstream)

    if(!missing(chromosome)){
        sequence <- .getSequenceFromCoords(chromosome, start, end, type, seqType, 
                                           upstream, downstream, mart, 
                                           useCache = useCache, verbose=verbose)
    }
    
    if(!missing(id)){
        sequence <- .getSequenceFromId(id, type, seqType, 
                                       upstream, downstream, mart, 
                                       useCache = useCache, verbose=verbose)
    }
    return(sequence)
}
