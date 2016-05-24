

###############################
#                             #
#Ensembl specific functions   #
###############################

## create the host name if we asked for a specific version
.useEnsemblVersion <- function(version = NULL, GRCh = NULL) {
    
    if(!is.null(version)){
        host = paste0("e",version,".ensembl.org")
    }
    
    if(!is.null(GRCh)){
        if(GRCh == 37){ 
            host = paste0("grch",GRCh,".ensembl.org")	
        }
        else{
            warning("Only 37 can be specified for GRCh version.  Using default host.")
            host <- "www.ensembl.org"
        }
    }
    
    return(host)
}

## convert name of mirror into full host name
.useEnsemblMirror <- function(mirror) {
    
    host <- switch(mirror, 
                   uswest = "uswest.ensembl.org",
                   useast = "useast.ensembl.org",
                   asia = "asia.ensembl.org",
                   "www.ensembl.org")
    ## warn if something invalid was specified and we're still using the default
    if(host == "www.ensembl.org") {
        warning("Invalid mirror.  Select a mirror from [uswest,useast,asia]. 
                The default when no mirror is specified points to main ensembl hosted in the UK")
    }
    
    return(host)
}

listEnsembl <- function(mart = NULL, host="www.ensembl.org", version = NULL, GRCh = NULL, mirror = NULL,verbose = FALSE){
    
    if( !is.null(version) || !is.null(GRCh) ) {
        host <- .useEnsemblVersion(version, GRCh)
        if(!is.null(mirror)) {
            warning("version or GRCh arguments can not be used together with the mirror argument.  
                    Will ignore the mirror argument and connect to default ensembl host") 
        }
    } else if(!is.null(mirror)) {
            host <- .useEnsemblMirror(mirror = mirror)
    }
    
    marts = listMarts(mart = mart, host = host, verbose = verbose)
    ## convert to easier to read names
    idx <- match(c("ENSEMBL_MART_ENSEMBL", "ENSEMBL_MART_SNP", "ENSEMBL_MART_FUNCGEN", "ENSEMBL_MART_VEGA"), marts$biomart)
    marts$biomart[ idx[!is.na(idx)] ] <- c("ensembl", "snp", "regulation", "vega")[ !is.na(idx) ]
    return(marts)
}

useEnsembl <- function(biomart, dataset,host = "www.ensembl.org", version = NULL, GRCh = NULL, mirror = NULL ,verbose = FALSE){
    
    if( !is.null(version) || !is.null(GRCh) ) {
        host <- .useEnsemblVersion(version, GRCh)
        if(!is.null(mirror)) {
            warning("version or GRCh arguments can not be used together with the mirror argument.  
                    Will ignore the mirror argument and connect to default ensembl host") 
        }
    } else if(!is.null(mirror)) {
        host <- .useEnsemblMirror(mirror = mirror)
    }
    
    biomart <- switch(biomart, 
                      ensembl = "ENSEMBL_MART_ENSEMBL", 
                      snp = "ENSEMBL_MART_SNP", 
                      regulation = "ENSEMBL_MART_FUNCGEN", 
                      vega = "ENSEMBL_MART_VEGA")
    
    ens = useMart(biomart = biomart, dataset = dataset, host = host, verbose = verbose)	   
    return(ens)
}

getGene <- function( id, type, mart){
    martCheck(mart,"ensembl") 
    checkWrapperArgs(id, type, mart)
    symbolAttrib = switch(strsplit(martDataset(mart), "_", fixed = TRUE, useBytes = TRUE)[[1]][1],
                          hsapiens = "hgnc_symbol",
                          mmusculus = "mgi_symbol",
                          "external_gene_id")
    typeAttrib = switch(type,
                        affy_hg_u133a_2 = "affy_hg_u133a_v2",
                        type)
    attrib = c(typeAttrib,symbolAttrib,"description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id")
    table = getBM(attributes = attrib,filters = type, values = id, mart=mart)
    return(table)
}

getSequence <- function(chromosome, start, end, id, type, seqType, upstream, downstream, mart, verbose=FALSE){
    
    martCheck(mart,c("ensembl","ENSEMBL_MART_ENSEMBL"))
    
    if(missing(seqType) || !seqType %in% c("cdna","peptide","3utr","5utr", "gene_exon", "transcript_exon","transcript_exon_intron","gene_exon_intron","coding","coding_transcript_flank","coding_gene_flank","transcript_flank","gene_flank")){
        stop("Please specify the type of sequence that needs to be retrieved when using biomaRt in web service mode.  Choose either gene_exon, transcript_exon,transcript_exon_intron, gene_exon_intron, cdna, coding,coding_transcript_flank,coding_gene_flank,transcript_flank,gene_flank,peptide, 3utr or 5utr")
    }
    if(missing(type))
        stop("Please specify the type argument.  If you use chromosomal coordinates to retrieve sequences, then the type argument will specify the type of gene indentifiers that you will retrieve with the sequences.  If you use a vector of identifiers to retrieve the sequences, the type argument specifies the type of identifiers you are using.")
    if(missing(id) && missing(chromosome) && !missing(type))
        stop("No vector of identifiers given. Please use the id argument to give a vector of identifiers for which you want to retrieve the sequences.")
    if(!missing(chromosome) && !missing(id))
        stop("The getSequence function retrieves sequences given a vector of identifiers specified with the id argument of a type specified by the type argument.  Or alternatively getSequence retrieves sequences given a chromosome, a start and a stop position on the chromosome.  As you specified both a vector of identifiers and chromsomal coordinates. Your query won't be processed.")
    
    if(!missing(chromosome)){
        if(!missing(start) && missing(end))
            stop("You specified a chromosomal start position but no end position.  Please also specify a chromosomal end position.")
        if(!missing(end) && missing(start))
            stop("You specified a chromosomal end position but no start position.  Please also specify a chromosomal start position.")
        if(!missing(start)){ 
            start = as.integer(start)
            end = as.integer(end)
        }
        if(missing(upstream) && missing(downstream)){
            sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end"), values = list(chromosome, start, end), mart = mart, checkFilters = FALSE, verbose=verbose)
        }
        else{
            if(!missing(upstream) && missing(downstream)){
                sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end","upstream_flank"), values = list(chromosome, start, end, upstream), mart = mart, checkFilters = FALSE, verbose=verbose)
            }
            if(!missing(downstream) && missing(upstream)){
                sequence = getBM(c(seqType,type), filters = c("chromosome_name","start","end","downstream_flank"), values = list(chromosome, start, end, downstream), mart = mart, checkFilters = FALSE, verbose = verbose)
            }
            if(!missing(downstream) && !missing(upstream)){
                stop("Currently getSequence only allows the user to specify either an upstream of a downstream argument but not both.")
            }
        }
    }
    
    if(!missing(id)){
        if(missing(type)) stop("Type argument is missing.  This will be used to retrieve an identifier along with the sequence so one knows which gene it is from.  Use the listFilters function to select a valid type argument.")
        if(!type %in% listFilters(mart, what="name")) 
            stop("Invalid type argument.  Use the listFilters function to select a valid type argument.")
        
        valuesString = paste(id,"",collapse=",",sep="")
        if(missing(upstream) && missing(downstream)){
            sequence = getBM(c(seqType,type), filters = type, values = id, mart = mart, verbose=verbose)
        }
        else{
            if(!missing(upstream) && missing(downstream)){
                sequence = getBM(c(seqType,type), filters = c(type, "upstream_flank"), values = list(id, upstream), mart = mart, checkFilters = FALSE, verbose=verbose)
            }
            if(!missing(downstream) && missing(upstream)){
                sequence = getBM(c(seqType,type), filters = c(type, "downstream_flank"), values = list(id, downstream), mart = mart, checkFilters = FALSE, verbose=verbose)
            }
            if(!missing(downstream) && !missing(upstream)){
                stop("Currently getSequence only allows the user to specify either an upstream of a downstream argument but not both.")
            }
        }
    }
    return(sequence)
}
