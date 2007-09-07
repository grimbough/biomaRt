packageName <- "biomaRt"

setClass("Mart",
         representation(mysql = "logical",
                        connections = "list",
                        mysqldriver = "list",
                        mainTables = "list",
                        biomart = "character",
                        host = "character",
                        port = "character",
                        vschema = "character",
                        dataset = "character",
                        filters = "environment",
                        attributes = "environment",
                        attributePointer = "environment"   
                        ),
         prototype(mysql = FALSE,
                   connections = new("list"),
                   dataset = "",
                   vschema="default"
                   )
         );


setMethod("show","Mart",
  function(object){	
    res = paste("Object of class 'Mart':\n Using the ",object@biomart," BioMart database\n Using the ",object@dataset," dataset\n", sep="")
    cat(res)
})


###############
#messageToUser#
###############

messageToUser <- function(message){
 if(interactive()){
  cat(message)
 }
}

##############################################################
#martCheck                                                   # 
#                                                            #
#This function checks if there is a valid Mart object,       # 
#if a dataset is selected and                                #
#if the correct BioMart database has been selected (optional)# 
##############################################################

martCheck = function(mart, biomart = NULL){
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("You must provide a valid Mart object. To create a Mart object use the function: useMart.  Check ?useMart for more information.")
  }
  if(!is.null(biomart)){
    martcheck = strsplit(mart@biomart,"_")[[1]][1]
    if(martcheck[1] != biomart)stop(paste("This function only works when used with the ",biomart," BioMart.",sep="")) 
   
  }
  if(mart@dataset==""){
    stop("No dataset selected, please select a dataset first.  You can see the available datasets by using the listDatasets function see ?listDatasets for more information.  Then you should create the Mart object by using the useMart function.  See ?useMart for more information");
  }
}

##############
#list marts  #
##############


listMarts <- function( mart, host, user, password, port, includeHosts = FALSE, mysql = FALSE, archive = FALSE){
  
  if(mysql){
    
#MySQL-----------------------------------------------------------    
    require(RMySQL)
    
    if(missing(mart)){
      mart = c("ensembl","vega","snp","msd","uniprot","sequence","wormbase")
    }
    if(missing(host)){
      host = c("martdb.ensembl.org", "martdb.ebi.ac.uk")
      user = c("anonymous","anonymous")
      password = c("","")
      port = c(3316, 3306) 
    }
    
    database = NULL
    driv = dbDriver("MySQL", force.reload = FALSE);
    
    for(i in seq(along=host)){
      connection = dbConnect(driv, user = user[i], host = host[i], password = password[i], port = port[i]);
  
      res = dbGetQuery(connection,"show databases like '%mart%'"); 
                                        #Search latest releases of marts
      if(dim(res)[1] >= 1){
        if(!archive){
         for(j in seq(along=mart)){ 
          matches = grep(mart[j],res[,1]);
          if(length(matches) > 1){
            version = 1;
            latest = 1;
            for(j in seq(along=matches)){
              v = suppressWarnings(as.numeric(strsplit(res[matches[j],1],"_")[[1]][3]));
              if(!is.na(v)){
                if(v > version){
                  latest = j;
                  version = v;
                }
              }
            }
            if(!includeHosts){
              database = c(database,res[matches[latest],1])
            }
            else{
              database = rbind(database,cbind(res[matches[latest],1],host[i]))
            }
          }
          else{
            if(sum(matches)> 0){
              if(!includeHosts){
                database = c(database,res[matches,1]);
              }
              else{
               database = rbind(database,cbind(res[matches,1],host[i]));
              }
            }
          }
          
          dbDisconnect(connection);
         }
        }
        else{
         database = rbind(database,res)
       } 
      }
    }
    return( database );
  }
  
#Webservice-----------------------------------------------------
  else{
    if(missing(host)){
      host = "http://www.biomart.org/biomart/martservice"
    }
    registry = getURL(paste(host,"?type=registry&requestid=biomaRt", sep=""))
    registry = xmlTreeParse(registry)
    registry = registry$doc$children[[1]]
    
    marts = list(biomart = NULL, version = NULL, host = NULL, path = NULL)
    index = 1
    
    for(i in seq(len=xmlSize(registry))){
      if(xmlName(registry[[i]])=="virtualSchema"){
           vschema = xmlGetAttr(registry[[i]],"name")
        for(j in seq(len=xmlSize(registry[[i]]))){
          if(xmlGetAttr(registry[[i]][[j]],"visible") == 1){
            marts$biomart[index] = xmlGetAttr(registry[[i]][[j]],"name")
            marts$version[index] = xmlGetAttr(registry[[i]][[j]],"displayName")
            marts$host[index] = xmlGetAttr(registry[[i]][[j]],"host")
            marts$path[index] = xmlGetAttr(registry[[i]][[j]],"path")
            marts$port[index] = xmlGetAttr(registry[[i]][[j]],"port")
            marts$vschema[index] = vschema
            index=index+1
          }
        }
      }
      if(xmlName(registry[[i]])=="MartURLLocation"){  
        if(xmlGetAttr(registry[[i]],"visible") == 1){
          marts$biomart[index] = xmlGetAttr(registry[[i]],"name")
          marts$version[index] = xmlGetAttr(registry[[i]],"displayName")
          marts$host[index] = xmlGetAttr(registry[[i]],"host")
          marts$path[index] = xmlGetAttr(registry[[i]],"path")
          marts$port[index] = xmlGetAttr(registry[[i]],"port")
          if(!is.null(xmlGetAttr(registry[[i]],"serverVirtualSchema"))){
           marts$vschema[index] =  xmlGetAttr(registry[[i]],"serverVirtualSchema")
          }
          else{
           marts$vschema[index] = vschema
          }
          index=index+1
        }
      }
    }
    if(includeHosts){
     return(marts)
    }
    else{
     ret = data.frame(marts$biomart, marts$version)
     colnames(ret)=c("name","version")
     return(ret)
    } 
  }
}

######################
#Disconnect from mart#
######################

martDisconnect <- function( mart ){
  
  if(missing(mart)){
    stop("No Mart object to disconnect");
  }
  openConnections = names(mart@connections);
  
  if(match("biomart",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$biomart );
  }
}

mapSpeciesToHomologTable <- function(fromESpecies = NULL, toESpecies = NULL) {
  
  table = NULL;
  table = paste(fromESpecies,"_gene_ensembl__homologs_",toESpecies,"__dm",sep="");
  return(table);
}

###########################################################################################
#                          
#Ensembl specific functions
##########################################################################################

######################
#get gene information#
######################

getGene <- function( id, type, mart){
  
  martCheck(mart,"ensembl") 

  if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details.  Note that the array argument is now depricated and should be specified with the type argument.") 
  if(!type %in% ls(mart@filters))stop("Invalid type of identifier see ?getGene for details of use listFilters function to get the valid value for type")   
  symbolAttrib = switch(strsplit(mart@dataset, "_")[[1]][1],hsapiens = "hgnc_symbol",mmusculus = "markersymbol","external_gene_id")
  typeAttrib = switch(type,affy_hg_u133a_2 = "affy_hg_u133a_v2",type)
  attrib = c(typeAttrib,symbolAttrib,"description","chromosome_name","band","strand","start_position","end_position","ensembl_gene_id")
  table = getBM(attributes = attrib,filters = type, values = id, mart=mart)

  return(table)
}

########################
#getFeature            #
########################

getFeature <- function( symbol, OMIMID, GOID, chromosome, start, end, type,  mart){

  martCheck(mart,"ensembl") 

  if( missing( type ))stop("You must provide the identifier type using the type argument.  Values for the type argument are given by the listAttributes function, see ?listAttributes for more information.")
  
  if(!type %in% ls(mart@filters))stop("Invalid type of identifier see ?getGene for details of use listFilters function to get the valid value for type")
  
  if(mart@mysql && !missing(chromosome) && !missing(start)){    
    
    speciesTable <- unique(mart@mainTables$tables[mart@mainTables$keys == "gene_id_key"]);
    
    IDTable <- get(type,mart@filters)$table
    if(IDTable == "main") IDTable = speciesTable
    dbcolID <- get(type,mart@filters)$field;
    query <- paste("select distinct a.",dbcolID,",b.chr_name,b.gene_chrom_start,b.gene_chrom_end from ", IDTable ," as a inner join ",speciesTable," as b on a.gene_stable_id = b.gene_stable_id where b.chr_name ='",chromosome,"' and b.gene_chrom_start >=",start," and b.gene_chrom_end <= ",end," and a.",dbcolID," !='NULL'",sep="");  
    res <- dbGetQuery( conn = mart@connections$biomart,statement = query);
    
    if(dim(res)[1] != 0){
      names(res) = c("id","chromosome", "start", "end");
      table <- as.data.frame(res)     
    }
    else{
      writeLines("No match found")
    }
    return( table )
  }
  
  else{
    
    filter=NULL
    startpos = "start"
    endpos = "end"
    chrname = "chromosome_name"
    attribute = type
    
    if(!missing(symbol)){
      filter = switch(strsplit(mart@dataset, "_")[[1]][1],hsapiens = "hgnc_symbol",mmusculus = "markersymbol","external_gene_id")
      attributes = c(filter,attribute)
      values = symbol
    }
  
    if(!missing(OMIMID)){
      filter = "mim_gene_ac"                      
      attributes = c(filter,attribute)
      values = OMIMID
    }
    
    if(!missing(GOID) && mart@mysql){
      filter = c("go")
      attributes = c("go",attribute)
      values = GOID
    }
    if(!missing(GOID) && !mart@mysql){
      filter = c("go", paste("with_",type,sep=""))
      attributes = c("go",attribute)
      values = list(GOID,TRUE)
    }
    
    if(!missing(chromosome)){
      if(missing(start) && missing(end)){
        if(attribute == "ensembl_gene_id" || attribute == "ensembl_transcript_id"){
          filter = "chromosome_name"
          values = chromosome
        }
        else{
          if(!mart@mysql){
            filter = c("chromosome_name", paste("with_", attribute, sep=""))
          }
          else{
            filter = "chromosome_name"
          }
          values = list(chromosome,TRUE)
        }
        attributes = c(chrname,attribute)
      }
      else{
        if(attribute == "ensembl_gene_id" || attribute == "ensembl_transcript_id"){
          filter = c(chrname,startpos,endpos)
          values = list(chromosome, start, end)
        }
        else{
          filter = c(chrname,startpos,endpos,paste("with_", attribute, sep=""))
          values = list(chromosome, start, end,TRUE)
        }
        attributes = c(chrname,"start_position","end_position",attribute)  
      }
    }
    table = getBM(attributes = attributes, filters = filter, values = values, mart=mart)
    return(table)
  }
}

###################
#get GO annotation#
###################


getGO <- function( id, type, mart){
   
  martCheck(mart,"ensembl")
  typeAttrib = switch(type,affy_hg_u133a_2 = "affy_hg_u133a_v2",type)
  attrib = c(typeAttrib,"go","go_description","evidence_code","ensembl_gene_id")   
  table = getBM(attributes = attrib,filters = type, values = id, mart=mart)
  return(table)
}

##############
#getSequence #
##############


getSequence <- function(chromosome, start, end, id, type, seqType, upstream, downstream, mart, verbose=FALSE){

  martCheck(mart,"ensembl")

  if(mart@mysql){
  
    martdb = "";
    if(!missing(upstream)){
     stop("Use of the upstream argument only works when using biomaRt in webservice mode")
    }
    if(!missing(downstream)){
     stop("Use of the downstream argument only works when using biomaRt in webservice mode")
    }
    if(!missing(seqType)){
     stop("seqType only has to be specified when using biomaRt in webservice mode.  In MySQL mode biomaRt will retrieve genomic sequences only")
    }
    if(!missing(type)){
     stop("Use of the type argument to specify gene identifiers only works when using biomaRt in webservice mode")
    }
    
    if(mart@biomart != "sequence"){
      version = "0"
      marts = listMarts(mysql = TRUE)
      for(i in seq(along=marts)){
        if("sequence" == strsplit(marts[i],"_")[[1]][1]){
          version = strsplit(marts[i],"_")[[1]][3]
        }
      }
      martdb = paste("sequence_mart_",version,sep="")
      mart@biomart="sequence"
      
    }
    mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = "anonymous", host = "martdb.ensembl.org" , port = 3316, dbname = martdb, password = "")
    
    species = strsplit(mart@dataset,"_")[[1]][1]
    sequence <- NULL;  
    speciesTable <- paste( species,"_genomic_sequence__dna_chunks__main",sep="" ); 
    
    if(!missing( chromosome ) && !missing( start ) && !missing( end )){
      for(i in seq(along = chromosome )){
        
        if(end[i] - start[i] > 100000){
          stop("maximum sequence length is 100000 nucleotides, change start and end arguments to make the sequence size smaller")
        }
        
        chunkStart <- (floor((start[i] - 1)/100000)*100000) + 1;
        chunkEnd <- (floor((end[i] - 1)/100000)*100000) + 1;
        
        if(chunkStart == chunkEnd ){  #we only need to get one sequence chunck of 100000 nucleotides
          
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
          chunkseq <- dbGetQuery(conn = mart@connections$biomart, statement = query);
          newstart <- start[i] - (floor((start[i]-1)/100000) * 100000)
          newend <- end[i] - (floor((end[i]-1)/100000) * 100000)
          
          sequence <- c(sequence,substr(as.character(chunkseq), newstart, newend));
          
        }
        
        else{   #query sequence is on 2 sequence chuncks
          
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
          chunkseq1 <- dbGetQuery(conn = mart@connections$biomart, statement = query);
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkEnd,"'",sep="");
          chunkseq2 <- dbGetQuery(conn = mart@connections$biomart, statement = query);
          chunkseq <- paste(as.character(chunkseq1),as.character(chunkseq2), sep=""); 
          
          newstart <- start[i] - (floor((start[i]-1)/100000) * 100000);
          newend <- end[i] - (floor((start[i]-1)/100000) * 100000);
          sequence <- c(sequence,substr(chunkseq, start = newstart, stop = newend));
          
        }
      }
    }
    
    table <- data.frame(chromosome = chromosome, start = start, end = end, sequence = sequence)
    
    return(table)
  }
  else{
    if(missing(seqType) || !seqType %in% c("cdna","peptide","3utr","5utr", "gene_exon", "transcript_exon","transcript_exon_intron","gene_exon_intron","coding","coding_transcript_flank","coding_gene_flank","transcript_flank","gene_flank")){
    stop("Please specify the type of sequence that needs to be retrieved when using biomaRt in web service mode.  Choose either gene_exon, transcript_exon,transcript_exon_intron, gene_exon_intron, cdna, coding,coding_transcript_flank,coding_gene_flank,transcript_flank,gene_flank,peptide, 3utr or 5utr")
    }   
    if(missing(type)) stop("Type argument is missing.  This will be used to retrieve an identifier along with the sequence so one knows which gene it is from.  Use the listFilters function to select a valid type argument.")
    if(!type %in% ls(mart@filters)) stop("Invalid type argument.  Use the listFilters function to select a valid type argument.")

    if(!missing(chromosome)){
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
}

###############################
#getAffyArrays: show affy info#
###############################

getAffyArrays <- function(mart){
  martCheck(mart,"ensembl")
  att=listFilters(mart)
  affy=att[grep("affy",att[,1]),]
  affy=affy[-grep("with",affy[,1]),]
  return(affy)
}

#######################
#getSNP
#######################

getSNP <- function(chromosome, start, end, mart){
  martCheck(mart,"snp")
  if(missing(chromosome) || missing(start) || missing(end) ){
    stop("You have to give chromosome, start and end positions as arguments, see ?getSNP for more information")
  }
  
  if(mart@mysql){ 
    ensemblTable <- unique(mart@mainTables$tables);
    query <- paste("select external_id,  tscid, snp_chrom_start,chrom_strand ,allele, ensemblcoding_bool, ensemblintronic_bool, ensembl5utr_bool, ensembl3utr_bool, ensemblsyn_bool from ",ensemblTable," where chr_name = '",chromosome,"' and snp_chrom_start >= ",start," and snp_chrom_start <= ",end,sep="")
    res <- dbGetQuery( mart@connections$biomart, query );
    if(dim(res)[1] == 0){
      stop("No SNP's found in selected region, check if the input you entered is correct")
    }
    return(res)
    table <- as.data.frame(res)
    colnames(res) = c("refsnp_id", "tscid","chromosome_start","strand" ,"allele","coding" , "intronic", "5utr", "3utr", "syn")
    return( table );
  }
  else{
    tscid="tscid"
    snpstart="chrom_start"
    snpend="chrom_end"
    if(exists("tscid",envir=mart@attributes)){
     attributes = c(tscid,"refsnp_id","allele","chrom_start","chrom_strand")
    }
    else{
      attributes = c("refsnp_id","allele","chrom_start","chrom_strand")
    } 
    table = getBM(attributes = attributes,filters = c("chr_name",snpstart,snpend), values = list(chromosome, start, end), mart=mart)    
    return(table)
  }
}

##################
#getHomolog      #
##################

getHomolog <- function(id, from.type, to.type, from.mart, to.mart) {

  martCheck(to.mart,"ensembl")
  martCheck(from.mart,"ensembl")
  
  if ( missing( from.type )) stop("You must provide the identifier type using the from.type argument")
  if (!from.type %in% ls(from.mart@filters)) stop("Invalid from.type, use the listFilters function on the from.mart to get valid from.type values")

  if ( missing( to.type )) stop("You must provide the identifier type using the to.type  argument")
  if (!to.type %in% ls(to.mart@attributes)) stop("Invalid to.type, use the listAttributes function on the to.mart to get valid to.type values")
  
  if(from.mart@mysql){
    if(!to.mart@mysql)stop("Both mart object should both be using mysql or not.  Your from.mart uses mysql but your to.mart not.")
    
    if( !missing( id )){
      id <- as.character(id);
    }
    else{
      stop("No ids to search for homologs");
    }
    
    fromIDTable <- get(from.type, from.mart@filters)$table
    fromSpeciesTable <- unique(from.mart@mainTables$tables[from.mart@mainTables$keys == "gene_id_key"]);    
    if(fromIDTable == "main") fromIDTable = fromSpeciesTable

    fromCol <- get(from.type, from.mart@filters)$field
    toIDTable <- get(to.type, to.mart@filters)$table
    toSpeciesTable <- unique(to.mart@mainTables$tables[to.mart@mainTables$keys == "gene_id_key"]);    
    if(toIDTable == "main") toIDTable = toSpeciesTable

    toCol <- get(to.type, to.mart@filters)$field;
    
    if(to.type == "ensembl" && from.type != "ensembl"){
      to.species = strsplit(to.mart@dataset,"_")[[1]][1]
      from.species = strsplit(from.mart@dataset,"_")[[1]][1] 
      homolTable <- mapSpeciesToHomologTable(to.species,from.species);
    }
    else{
      to.species = strsplit(to.mart@dataset,"_")[[1]][1]
      from.species = strsplit(from.mart@dataset,"_")[[1]][1]
      homolTable <- mapSpeciesToHomologTable(from.species,to.species);
    }
    
    if (length(id) >= 1) {
      ids <- paste("'",id,"'",sep="",collapse=",")
      res <- NULL
      
      if(from.type == "ensembl" && to.type != "ensembl"){
        query <-  paste("select distinct c.",fromCol,",b.",toCol," from ",
                        homolTable," as c inner join ",
                        toIDTable," as b on b.gene_id_key=c.homol_id ",
                        "where c.",fromCol," in (",ids,") and b.",toCol," != 'NULL'",sep="");
        
      }
      else{
        if(to.type == "ensembl" && from.type != "ensembl"){    
          query <-  paste("select distinct b.",fromCol,",c.",toCol," from ",
                          homolTable," as c inner join ",
                          fromIDTable," as b on b.gene_id_key=c.homol_id ",
                          "where b.",fromCol," in (",ids,") and c.",toCol," != 'NULL'",sep="");
        }
        else{
          if(to.type == "ensembl" && from.type == "ensembl"){
            query <-  paste("select distinct ",fromCol,", homol_stable_id from ",
                            homolTable," where ",fromCol," in (",ids,") and homol_stable_id != 'NULL'",sep="");
          }
          else{
            query <-  paste("select distinct a.",fromCol,",b.",toCol," from  ",
                            fromIDTable," as a inner join ",
                            homolTable," as c on a.gene_id_key=c.gene_id_key inner join ",
                            toIDTable," as b on b.gene_id_key=c.homol_id ",
                            "where a.",fromCol," in (",ids,") and b.",toCol," != 'NULL'",sep="");
          }
        }
      }
      
      res <- dbGetQuery(conn = from.mart@connections$biomart, statement = query);
      
      if (dim(res)[1] == 0) {
        return(NULL)
      }
      
      else {
        
        foundID <- NULL
        MappedID <- NULL
        isNA <- is.na(res[ ,2])
        res <- res[!isNA ,]
        
        for (j in seq(along=id)) {
          
          m <- match(res[, 1], id[j], nomatch = 0)
          
          if (sum(m) == 0) {
            
            foundID <- c(foundID, as.character(id[j]))
            MappedID <- c(MappedID, NA)
          }
          
          else {
            
            foundID <- c(foundID, res[m == 1, 1])
            MappedID <- c(MappedID, res[m == 1, 2])
            
          }
        }
        table <- data.frame(id = foundID, MappedID = MappedID)
        ind = is.na(MappedID)
        table = table[!ind,]
      }
    }
    return(table)
  }
  else{
    if(to.mart@mysql)stop("The Mart objects should either use both mysql or not. Here your from.mart does not use mysql but your to.mart does.")
    result = getLDS(attributes = from.type, filters = from.type, values = id, mart = from.mart, attributesL = to.type, martL = to.mart) 
    return(result)
  }
} 

####################
#export FASTA      #
####################

exportFASTA <- function( sequences, file ){
  if( missing( sequences ) || class( sequences ) != "data.frame"){
    stop("No data.frame given to write FASTA.  The data.frame should be the output of the getSequence function.");
  }
  if( missing(file)){
    stop("Please provide filename to write to");
  }
  if(length(sequences[1,]) == 2){
   for(i in seq(along = sequences[,2])){
     cat(paste(">",sequences[i,2],"\n",sep=""),file = file, append=TRUE);
     cat(as.character(sequences[i,1]),file = file, append = TRUE);
     cat("\n\n", file = file, append = TRUE);
   }
  }
  else{
   for(i in seq(along = sequences[,2])){
     cat(paste(">chromosome_",sequences[i,1],"_start_",sequences[i,2],"_end_",sequences[i,3],"\n",sep=""),file = file, append=TRUE);
     cat(as.character(sequences[i,4]),file = file, append = TRUE);
     cat("\n\n", file = file, append = TRUE);
   }
 }  
}

#################################
# #                           # #
# # Generic BioMart functions # #
# #                           # #
#################################

#######################
#useMart              # 
#######################

useMart <- function(biomart, dataset, host, user, password, port, local = FALSE, mysql = FALSE, archive = FALSE){

  if(mysql){
    
    require(RMySQL)
    driver <- dbDriver("MySQL", force.reload = FALSE);
    
    mart <- new("Mart", biomart = biomart, mysqldriver = list(driver=driver), mysql = TRUE)

    if(! (is.character(biomart) && (length(biomart)==1)))
      stop("'biomart' should be a single character string.")

    if(local){
      if(!missing(host) && !missing(user) && !missing(password) && !missing(port)){
        database <- listMarts(mart = biomart, host = host, user = user, password = password, port = port, mysql=TRUE);
        mart@biomart = strsplit(biomart,"_")[[1]][1]
        mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = user, host = host, dbname = database, password = password)
      }
      else{
        stop(sprintf("Please provide host, user, password and port for using local database '%s'.", biomart))
      }
    }
    else {
      
      version = "0"
      marts=NULL
      
      if(!archive){ 
       marts = listMarts(mysql = TRUE)
       for(i in seq(along=marts)){
        if(biomart == strsplit(marts[i],"_")[[1]][1]){
          version = strsplit(marts[i],"_")[[1]][3]
        }
       }
       martdb=""
       if(version > 0){
         martdb = paste(biomart,"_mart_",version,sep="")
       }
       else{
         martdb = biomart
       }
      }
      else{
        marts = listMarts(mysql = TRUE, archive = archive)
        marts=marts[,1]  
        martdb = biomart    
      }    

      if(!martdb %in% marts) stop("Requested BioMart database is not available please use the function listMarts(mysql=TRUE) to see the valid biomart names you can query using mysql access")
      mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = "anonymous", host = "martdb.ensembl.org" , dbname = martdb, password = "", port = 3316)
      messageToUser(paste("connected to: ",biomart,"\n"))
    }
    if(!missing(dataset)){
      mart = useDataset(mart = mart,dataset = dataset)
    }
    return( mart )
  }
  else{
    if(missing(host)){
      host = "http://www.biomart.org/biomart/martservice"
    }
    marts=listMarts(host = host, includeHosts = TRUE)
    mindex=match(biomart,marts$biomart)
    if(is.na(mindex)){
      stop("Incorrect BioMart name, use the listMarts function to see which BioMart databases are available")
    }

    if(marts$path[mindex]=="")marts$path[mindex]="/biomart/martservice" #temporary to catch bugs in registry
    biomart = sub(" ","%20",biomart)
    mart <- new("Mart", biomart = biomart,vschema = marts$vschema[mindex], host = paste("http://",marts$host[mindex],":",marts$port[mindex],marts$path[mindex],sep=""), mysql= FALSE)
    if(!missing(dataset)){
      mart = useDataset(mart = mart, dataset=dataset)
    }
    return(mart)
  }
}

####################
#listDatasets      #
####################

listDatasets <- function(mart) {
  if(missing(mart) || !is(mart, 'Mart'))
    stop("No Mart object given or object not of class 'Mart'")

  if(mart@mysql){
    
    res = dbGetQuery(mart@connections$biomart,
      "select dataset, version from meta_conf__dataset__main where visible = 1")
    
  } else{
    
    txt = scan(paste(mart@host,"?type=datasets&requestid=biomaRt&mart=",mart@biomart,sep=""),
      sep="\t", blank.lines.skip=TRUE, what="character", quiet=TRUE)

    ## select visible ("1") table sets
    i = intersect(which(txt=="TableSet"), which(txt=="1")-3L)
                                         
    res = data.frame(dataset     = I(txt[i+1L]),
                     description = I(txt[i+2L]),
                     version     = I(txt[i+4L]))

  }

  return(res)
}

############################################
#useDataset: select a BioMart dataset      #             
############################################

useDataset <- function(dataset, mart){

  if(missing(mart) || class(mart)!="Mart") stop("No valid Mart object given, specify a Mart object with the attribute mart")
  if(missing(dataset)) stop("No valid dataset to use given")
  validDatasets=listDatasets(mart)
  if(is.na(match(dataset, validDatasets$dataset)))stop(paste("The given dataset: ",dataset,", is not valid.  Correct dataset names can be obtained with the listDatasets function"))

  filtersEnv = new.env(parent = emptyenv(), hash = TRUE)
  attributesEnv = new.env(parent = emptyenv(), hash = TRUE)
  attributePointerEnv = new.env(parent = emptyenv(), hash = TRUE)
   
  if(mart@mysql){
    res = dbGetQuery(mart@connections$biomart,paste("select xml from meta_conf__dataset__main inner join meta_conf__xml__dm on meta_conf__dataset__main.dataset_id_key = meta_conf__xml__dm.dataset_id_key where dataset = '",dataset,"'",sep=""))
    messageToUser(paste("Reading database configuration of:",dataset,"\n"))
    if(dim(res)[1] == 0) stop("This dataset is not accessible from biomaRt as not xml description of dataset is available")

    xml = xmlTreeParse(res[1,])
    xml = xml$doc$children[[1]]

    messageToUser("Checking attributes and filters ...")   
    parseAttributes(xml, env = attributesEnv, attributePointer = attributePointerEnv)
    parseFilters(xml, filtersEnv)
    messageToUser(" ok\n")
 
    mainTables = getMainTables(xml)

    mart@attributes = attributesEnv
    mart@attributePointer = attributePointerEnv 
    mart@filters = filtersEnv
    mart@dataset = dataset
    mart@mainTables = mainTables
    
    return(mart)
  }
  else{

    config = getURL(paste(mart@host,"?type=configuration&requestid=biomaRt&dataset=",dataset,"&virtualschema=",mart@vschema, sep=""))
    config = xmlTreeParse(config)
    config = config$doc$children[[1]]
 
    messageToUser("Checking attributes and filters ...")
    parseAttributes(config, env = attributesEnv, attributePointer = attributePointerEnv)
    parseFilters(config, filtersEnv)
    messageToUser(" ok\n")
    
    mart@dataset = dataset
    mart@attributes = attributesEnv
    mart@attributePointer = attributePointerEnv 
    mart@filters = filtersEnv
    
    return( mart )
  }
}

###################
#getName          #
###################

getName = function(x, pos) if(is.null(x[[pos]])) NA else x[[pos]] 

#####################
#listAttributes     #
#####################

listAttributes = function( mart , group, category, showGroups = FALSE){

  martCheck(mart)
  summaryA = attributeSummary(mart)
  if(!missing(category) && !category %in% summaryA[,1]) stop(paste("The chosen category: ",category," is not valid, please use the right category name using the attributeSummary function",sep=""))
  if(!missing(group) && !group %in% summaryA[,2]) stop(paste("The chosen group: ",group," is not valid, please use the right group name using the attributeSummary function",sep=""))

  attribList = mget(ls(mart@attributes), env=mart@attributes)
  frame = data.frame(name = names(attribList),description = sapply(attribList, getName, 1), group = sapply(attribList, getName, 5), category = sapply(attribList, getName, 6),  row.names=NULL, stringsAsFactors=FALSE) 
  if(!missing(category)){
   frame = frame[frame[,4] == category,]
  }
  if(!missing(group)){
   frame = frame[frame[,3] == group,]
  }
  ord = order(frame[,4])
  frameOut = frame[ord,]
  if(!showGroups) frameOut = frameOut[,-c(3,4)]
  rownames(frameOut) = seq(len=length(frameOut[,1]))
  return(frameOut)

}

######################
#attributeSummary    #
######################

attributeSummary = function( mart ){

  martCheck(mart)
  attribList = mget(ls(mart@attributes), env=mart@attributes)
  frame = data.frame(category = sapply(attribList, getName, 6), group = sapply(attribList, getName, 5),  row.names=NULL, stringsAsFactors=FALSE) 
  frame = unique(frame)
  ord = order(frame[,1])
  frameOut = frame[ord,]
  rownames(frameOut) = seq(len=length(frameOut[,1]))
  return(frameOut)

}

###################
#listFilters      #
###################

listFilters = function( mart , group, category, showGroups = FALSE, showType = FALSE){

  martCheck(mart)
  summaryF = filterSummary(mart)
  if(!missing(category) && !category %in% summaryF[,1]) stop(paste("The chosen category: ",category," is not valid, please use the right category name using the filterSummary function",sep=""))
  if(!missing(group) && !group %in% summaryF[,2]) stop(paste("The chosen group: ",group," is not valid, please use the right group name using the filterSummary function",sep=""))

  filterList = mget(ls(mart@filters), env=mart@filters)
  frame = data.frame(name = names(filterList),description = sapply(filterList, getName, 1), group = sapply(filterList, getName, 5), category = sapply(filterList, getName, 6),type =  sapply(filterList, getName, 7),  row.names=NULL, stringsAsFactors=FALSE) 

  if(!missing(category)){
   frame = frame[frame[,4] == category,]
  }
  if(!missing(group)){
   frame = frame[frame[,3] == group,]
  }
  ord = order(frame[,4])
  frameOut = frame[ord,]
  if(!showType) frameOut = frameOut[,-5]
  if(!showGroups) frameOut = frameOut[,-c(3,4)]
  rownames(frameOut) = seq(len=length(frameOut[,1]))
  return(frameOut)

}

#################
#filterSummary  #
#################

filterSummary = function( mart ){
  martCheck(mart)
  filterList = mget(ls(mart@filters), env=mart@filters)
  frame = data.frame(category = sapply(filterList, getName, 6), group = sapply(filterList, getName, 5),  row.names=NULL, stringsAsFactors=FALSE) 
  frame = unique(frame)
  ord = order(frame[,1])
  frameOut = frame[ord,]
  rownames(frameOut) = seq(len=length(frameOut[,1]))
  return(frameOut)

}

###################
#filterOptions    #
###################

filterOptions = function(filter, mart){
  martCheck(mart)
  if(!filter %in% ls(mart@filters))stop("Filter not valid")
  options = get(filter, env=mart@filters)$options
  return(options)
}

###################
#filterType       #
###################

filterType = function(filter, mart){
  martCheck(mart)
  if(!filter %in% ls(mart@filters))stop("Filter not valid")
  type = get(filter, env=mart@filters)$type
  return(type)
}

##########################################
#getBM: generic BioMart query function   # 
##########################################

getBM = function(attributes, filters = "", values = "", mart, curl = NULL, output = "data.frame", list.names = NULL, na.value = NA, checkFilters = TRUE, verbose=FALSE, uniqueRows=TRUE){

  martCheck(mart)
  if(missing( attributes ))
    stop("Argument 'attributes' must be specified.")
  if(filters != "" && missing( values ))
    stop("Argument 'values' must be specified.")
  if(output != "data.frame" && output != "list")
    stop("Only data.frame and list are valid output formats for this function")
  if(class(uniqueRows) != "logical")
    stop("Argument 'uniqueRows' must be a logical value, so either TRUE or FALSE")

  invalid = !(attributes %in% ls(mart@attributes))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))
  if(filters[1] != "" && checkFilters){
   invalid = !(filters %in% ls(mart@filters))
   if(any(invalid))
    stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
               "\nPlease use the function 'listFilters' to get valid filter names"))
  }

  for(k in seq(along=attributes)){
   if(attributes[k] %in% ls(mart@attributePointer)) attributes[k] = get(attributes[k], env = mart@attributePointer)
  }  
  
  ## use the mySQL interface
  if(mart@mysql){
    if(output == "data.frame"){
      query = queryGenerator(attributes=attributes, filter=filters, values=values, mart=mart)
      if(verbose) print(query)
      res = dbGetQuery(mart@connections$biomart,query)
      if(dim(res)[1] == 0){
        res = data.frame()
      }
      else{
       colnames(res) = attributes
      }
      return(as.data.frame(res))
    }
    else{
      out = vector("list", length(attributes))
      if(is.null(list.names))
        names(out) = attributes
      else
        names(out) = list.names
      
      for(j in seq(along = attributes)){
        tmp = getBM(c(filters,attributes[j]), filters, values, mart)
        tmp2 = vector("list", length(values))
        names(tmp2) = values
        for(i in seq(along=tmp2)){
          tst = tmp[tmp[,1] %in% values[i],2]
          tst = tst[!is.na(tst)]
          if(length(tst) == 0) tst <- na.value
          tmp2[[i]] = tst
        }
        out[[j]] = tmp2
      }
     return(out)
    }
  }
  else{
    ## use the http/XML interface
    if(output == "data.frame"){
      xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '",as.numeric(uniqueRows),"' count = '0' datasetConfigVersion = '0.6' requestid= \"biomaRt\"> <Dataset name = '",mart@dataset,"'>",sep="")
      attributeXML =  paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
      if(length(filters) > 1){
        if(class(values)!= "list")
          stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on Affymetrix identifier and chromosome, only results that pass both filters will be returned");
        filterXML = NULL
        for(i in seq(along = filters)){
          
          if(filters[i] %in% ls(mart@filters)){
            if(get(filters[i],env=mart@filters)$type == 'boolean' || get(filters[i],env=mart@filters)$type == 'boolean_num'){
              if(!is.logical(values[[i]])) stop(paste("biomaRt error: ",filters[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="")  
              if(!values[[i]]){
                values[[i]] = 1
              }
              else{
               values[[i]] = 0 
              }
              filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' excluded = \"",values[[i]],"\" />", collapse="",sep=""),sep="")
            }
            else{
              valuesString = paste(values[[i]],"",collapse=",",sep="")
              filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
            }
          }
          else{ #used for attributes with values as these are treated as filters in BioMart
            valuesString = paste(values[[i]],"",collapse=",",sep="")
            filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
          } 
        }
      }
      else{
       if(filters != ""){
        if(filters %in% ls(mart@filters)){
          if(get(filters,env=mart@filters)$type == 'boolean' || get(filters,env=mart@filters)$type == 'boolean_num'){
           if(!is.logical(values)) stop(paste("biomaRt error: ",filters," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
           if(!values){
            values = 1
           }
           else{
            values = 0 
           }
           filterXML = paste("<Filter name = '",filters,"' excluded = \"",values,"\" />", collapse="",sep="")
          }
          else{
           valuesString = paste(values,"",collapse=",",sep="")
           filterXML = paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep="")
          }
        }
        else{ #used for attributes with values as these are treated as filters in BioMart
          valuesString = paste(values,"",collapse=",",sep="")
          filterXML = paste(filterXML,paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
        }
       }
       else{
         filterXML=""
       }
      }
      xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset></Query>",sep="")
       
      if(verbose){
       cat(paste(xmlQuery,"\n", sep=""))
      }      

      if(is.null(curl)){
        postRes = postForm(paste(mart@host,"?",sep=""),"query" = xmlQuery)
      }
      else{
        postRes = postForm(paste(mart@host,"?",sep=""),"query" = xmlQuery, curl = curl)
      }
      
      if(postRes != ""){
        if(postRes != "\n" && postRes != "\n\n"){
        ## convert the serialized table into a dataframe
        con = textConnection(postRes)
        result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
        close(con)
        if( substr(as.character(result[1]),1,5) == "ERROR"){  #Will forward the webservice error
          stop(paste("\n",result,"The following query was attempted, use this to report this error\n",xmlQuery ,sep="\n"))
        }

        if(checkFilters){
        stopifnot(ncol(result)==length(attributes))
         if(class(result) == "data.frame"){
           colnames(result) = attributes
         }
        }
       }
       else{
        result = NA
       } 
      } else {
       
        geturl = getURL(mart@host)
        if(geturl == "" || is.null(geturl)){
          stop(paste("The getBM query to BioMart webservice returned no result.  The webservice could be temporarily down, please check if the following URL is active: ",mart@host,".  If this URL is not active then try your query again at a later time when this URL is active.", sep=""))
        }
        else{ 
         result = NULL
        }
      }
      return(result)
    }
    else{
      out <- vector("list", length(attributes))
      if(is.null(list.names))
        names(out) <- attributes
      else
        names(out) <- list.names
      #curl <- getCurlHandle()
      for(j in seq(along = attributes)){
        tmp2 <- vector("list", length(values))
        names(tmp2) <- values
        for(k in seq(along = tmp2)){
       #   tst <- getBM(attributes = attributes[j], filters=filters, values = values[k], mart = mart, curl = curl)
          tst <- getBM(attributes = attributes[j], filters=filters, values = values[k], mart = mart, verbose = verbose)
       
          if(class(tst) == "data.frame"){
            tmp <- unlist(unique(tst[!is.na(tst)]), use.names = FALSE)
            if(length(tmp) > 0)
              tmp2[[k]] <- tmp
            else
              tmp2[[k]] <- na.value
          }else{
            tmp2[[k]] <- na.value
          }
          out[[j]] <- tmp2
        }
      }
      return(out)
    }
  }
}

###################################
#getLDS: Multiple dataset linking #
###################################


getLDS <- function(attributes, filters = "", values = "", mart, attributesL, filtersL = "", valuesL = "", martL, verbose = FALSE, uniqueRows = TRUE) {
  
  martCheck(mart)
  martCheck(martL)

  if(mart@mysql || martL@mysql)stop("This function only works with biomaRt in webservice mode")
  
  invalid = !(attributes %in% ls(mart@attributes))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))

  invalid = !(attributesL %in% ls(martL@attributes))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributesL[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))
 
  if(filters[1] != ""){
   invalid = !(filters %in% ls(mart@filters))
   if(any(invalid))
    stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
               "\nPlease use the function 'listFilters' to get valid filter names"))
  }
  if(filtersL[1] != ""){
   invalid = !(filtersL %in% ls(martL@filters))
   if(any(invalid))
    stop(paste("Invalid filters(s):", paste(filtersL[invalid], collapse=", "),
               "\nPlease use the function 'listFilters' to get valid filter names"))
  }
   
    xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '",as.numeric(uniqueRows),"' count = '0' datasetConfigVersion = '0.6' requestid= \"biomaRt\"> <Dataset name = '",mart@dataset,"'>",sep="")
    attributeXML = paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
    if(length(filters) > 1){
        if(class(values)!= "list")
          stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on affyid and chromosome, only results that pass both filters will be returned");
        filterXML = NULL
        for(i in seq(along=filters)){
          if(get(filters[i],env=mart@filters)$type == 'boolean' || get(filters[i],env=mart@filters)$type == 'boolean_num'){
            if(!is.logical(values[[i]])) stop(paste("biomaRt error: ",filters[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
            if(!values[[i]]){
              values[[i]] = 1
            }
            else{
              values[[i]] = 0 
            }
            filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' excluded = \"",values[[i]],"\" />", collapse="",sep=""),sep="")
          }
          else{
            valuesString = paste(values[[i]],"",collapse=",",sep="")
            filterXML = paste(filterXML,paste("<Filter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
          }
        }
      }
      else{
       if(filters != ""){       
        if(get(filters,env=mart@filters)$type == 'boolean' || get(filters,env=mart@filters)$type == 'boolean_num'){
          if(!is.logical(values)) stop(paste("biomaRt error: ",filters," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
          if(!values){
            values = 1
          }
          else{
            values = 0 
          }
 
         filterXML = paste("<Filter name = '",filters,"' excluded = \"",values,"\" />", collapse="",sep="")
        }
        else{
         valuesString = paste(values,"",collapse=",",sep="")
         filterXML = paste("<Filter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep="")
        }
       }
       else{
         filterXML=""
       }
      }

    xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset>",sep="")
    xmlQuery = paste(xmlQuery, "<Dataset name = '",martL@dataset,"' >", sep="")
    linkedAttributeXML =  paste("<Attribute name = '", attributesL, "'/>", collapse="", sep="")
    
    if(length(filtersL) > 1){
        if(class(valuesL)!= "list")
          stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on affyid and chromosome, only results that pass both filters will be returned");
        linkedFilterXML = NULL
        for(i in seq(along=filtersL)){
          if(get(filtersL[i],env=martL@filters)$type == 'boolean' || get(filtersL[i],env=martL@filters)$type == 'boolean_num'){
            if(!is.logical(valuesL[[i]])) stop(paste("biomaRt error: ",filtersL[i]," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
            if(!valuesL[[i]]){
              valuesL[[i]] = 1
            }
            else{
              valuesL[[i]] = 0 
            } 
            linkedFilterXML = paste(linkedFilterXML,paste("<Filter name = '",filtersL[i],"' excluded = \"",valuesL[[i]],"\" />", collapse="",sep=""),sep="")
          }
          else{
            valuesString = paste(valuesL[[i]],"",collapse=",",sep="")
            linkedFilterXML = paste(linkedFilterXML,paste("<Filter name = '",filtersL[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
          }
        }
      }
      else{
       if(filtersL != ""){       
        if(get(filtersL,env=martL@filters)$type == 'boolean' || get(filtersL,env=martL@filters)$type == 'boolean_num'){
          if(!is.logical(valuesL)) stop(paste("biomaRt error: ",filtersL," is a boolean filter and needs a corresponding logical value of TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or alternatively that all data that not fulfill the requirement should be retrieved."), sep="") 
          if(!valuesL){
            valuesL = 1
          }
          else{
            valuesL = 0 
          }

         linkedFilterXML = paste("<Filter name = '",filtersL,"' excluded = \"",valuesL,"\" />", collapse="",sep="")
        }
        else{
         valuesString = paste(valuesL,"",collapse=",",sep="")
         linkedFilterXML = paste("<Filter name = '",filtersL,"' value = '",valuesString,"' />", collapse="",sep="")
        }
       }
       else{
         linkedFilterXML=""
       }
      }

    xmlQuery = paste(xmlQuery, linkedAttributeXML, linkedFilterXML,"</Dataset></Query>",sep="")

    if(verbose){
      cat(paste(xmlQuery,"\n", sep=""))
    }
    postRes = postForm(paste(mart@host,"?",sep=""),"query"=xmlQuery)
    
    if(postRes != ""){
      con = textConnection(postRes)
      result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
      close(con)
      if(all(is.na(result[,ncol(result)])))
        result = result[,-ncol(result),drop=FALSE]
     } else {
      warning("getLDS returns NULL.")
      result=NULL
    }
    return(result)
} 

##########################################
#
#Configuration Parsers
##########################################


parseAttributes = function(xml, env, attributePointer, group = "", page = ""){
   
   select = TRUE
   
    if(xmlName(xml) == "AttributePage"){
      if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE
      if(!is.null(xmlGetAttr(xml,"hideDisplay")) && xmlGetAttr(xml,"hideDisplay") == "true") select = FALSE
      if(select){
        page = xmlGetAttr(xml,"displayName")
        xmlSApply(xml,"parseAttributes",env, attributePointer = attributePointer, page = page)
      }
    }

    if(xmlName(xml) == "AttributeDescription"){
      if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE
      if(select){
        description = xmlGetAttr(xml,"displayName")
        if(is.null(description)) description = NA
        if(!is.null(xmlGetAttr(xml,"pointerAttribute"))) assign(xmlGetAttr(xml,"internalName"),xmlGetAttr(xml,"pointerAttribute"),env = attributePointer)
        assign(xmlGetAttr(xml,"internalName"),list(description = description, table = xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field"), group = group, page = page),env = env)        
      } 
    }

    if(xmlName(xml)=="AttributeCollection"){
      if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE
      if(select){
         xmlSApply(xml,"parseAttributes",env, attributePointer = attributePointer, group = group, page = page)
      }
    }    
    
    if(xmlName(xml)=="AttributeGroup"){
      if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE
      if(select){
        group = xmlGetAttr(xml,"displayName") 
        xmlSApply(xml,"parseAttributes",env, attributePointer = attributePointer, group = group, page = page)
      }
    } 

   if(xmlName(xml)=="DatasetConfig"){
        xmlSApply(xml,"parseAttributes",env, attributePointer = attributePointer)
   }

   return()
}

parseFilters = function(xml, env, group = "", page = ""){
  
  select = TRUE

  if(xmlName(xml) == "FilterPage"){
    if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE 
    if(!is.null(xmlGetAttr(xml,"hideDisplay")) && xmlGetAttr(xml,"hideDisplay") == "true") select = FALSE 
    if(select){
      page = xmlGetAttr(xml,"displayName")
      xmlSApply(xml,"parseFilters",env, page = page)
    }
  }

  if(xmlName(xml)=="FilterGroup"){
    if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE 
    if(select){
      group = xmlGetAttr(xml,"displayName") 
      xmlSApply(xml,"parseFilters",env, group = group, page = page)
    }
  }

  if(xmlName(xml)=="FilterCollection"){
    if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE 
    if(select){
      xmlSApply(xml,"parseFilters",env, group = group, page = page)
    }
  }

 if(xmlName(xml) == "FilterDescription"){
   hidden = xmlGetAttr(xml,"hidden")
   if(!is.null(hidden) && hidden == "true") select = FALSE
   displayType = xmlGetAttr(xml,"displayType")
   
   if(select && !is.null(displayType) && displayType == "container"){
     xmlSApply(xml,"parseFilters",env, group = group, page = page)
   }

   if((select && !is.null(displayType) && displayType != "container") || (select && is.null(displayType))){
     description = xmlGetAttr(xml,"displayName")
     if(is.null(description))description = NA
     type = xmlGetAttr(xml,"type")
     options = type 
     if(!is.null(displayType) && xmlGetAttr(xml,"displayType") == "list" && type != "boolean") options = parseOptions(xml)
     assign(xmlGetAttr(xml,"internalName"),list(descr = description, table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field"), group = group, page = page, type = type, options = options), env = env)
   }
 } 

 if(xmlName(xml) == "Option"){
   if(!is.null(xmlGetAttr(xml,"hidden")) && xmlGetAttr(xml,"hidden") == "true") select = FALSE
   if(select){
     description = xmlGetAttr(xml,"displayName");
     if(is.null(description))description = NA
     type = xmlGetAttr(xml,"type")
     options = type   
     if( xmlGetAttr(xml,"displayType") == "list" && type != "boolean") options = parseOptions(xml)
     assign(xmlGetAttr(xml,"internalName"),list(descr = description, table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field"), group = group, page = page, type = type, options = options),env = env)
   }
 }
  
  if(xmlName(xml)=="DatasetConfig"){
    xmlSApply(xml,"parseFilters",env)
  }
  return()  
}

parseOptions = function(xml){
 options = NULL
 for(h in seq(len=xmlSize(xml))){
    if(!is.null(xmlGetAttr(xml[[h]],"isSelectable")) && xmlGetAttr(xml[[h]],"isSelectable") == "true"){
      options = c(options,xmlGetAttr(xml[[h]],"value"))
    }
 }
 return(options)
}


getMainTables = function( xml ){
  
  messageToUser("Checking main tables ...")
  names = names(xml)
  
  tableM = NULL
  keyM = NULL
  i=1
  j=1
  for(h in seq(len = xmlSize(xml)-1)){
    if(names[h] == "MainTable"){
      tableM[i] = xmlValue(xml[[h]])
      i = i+1
    }
    if(names[h] == "Key"){
      keyM[j] = xmlValue(xml[[h]])
      j = j+1
    }
  }
  messageToUser(" ok\n")
 
  return(list(tables=tableM,keys=keyM))
}

########################
#MySQL query generator #
########################

queryGenerator <- function(attributes, filter, values, mart){

  ## attributes
  Afield <- Atab <- Akey <- NULL
  for(i in seq(along = attributes)){
    if(!exists(attributes[i],mart@attributes)){
      stop(paste("Attribute: ",attributes[i]," not found, please use the function 'listAttributes' to get valid attribute names",sep=""))
    }
    else{
      Afield = c(Afield,get(attributes[i],mart@attributes)$field)
      if(is.null(get(attributes[i],mart@attributes)$table)){
          stop(paste("Attribute: ",attributes[i]," is not retrievable via BioMart in MySQL mode, please perform query in the default web service mode", sep=""))
      }
      Atab   = c(Atab,get(attributes[i],mart@attributes)$table)
      Akey   = c(Akey,get(attributes[i],mart@attributes)$key)
    }
  }
    
  ## special treatment of main table
  isMain = (Atab=="main")
  m2 = match(Akey, mart@mainTables$keys)
  if(any(is.na(m2) & isMain))
    stop("Internal error: Key does not match MainTable key")
  Atable = ifelse(isMain, mart@mainTables$tables[m2], Atab)
  
  ## filter
  ## FIXME: I assume this is correct - filter can have only length 1?
  ## YES when accessing BioMart databases via MySQL the number of filters is currently limited to 1.  Access via the webservice however allows more filters to be used
  
  if(length(filter)!=1)
    stop(sprintf("'length(filter)' must be 1 when using biomaRt in MySQL mode."))
  
  if(!exists(filter,mart@filters))
    stop(paste("filter(s)",filter," not found, please use the function 'listFilters' to get valid filter names",sep=""))
  
  Ffield = get(filter,mart@filters)$field
  Ftab   = get(filter,mart@filters)$table
  Fkey   = get(filter,mart@filters)$key
  
  ## special treatment of main table
  isMain = (Ftab=="main")
  m2 = match(Fkey, mart@mainTables$keys)
  if(any(is.na(m2) & isMain))
    stop("Internal error: Key does not match MainTable key")
  Ftable = ifelse(isMain, mart@mainTables$tables[m2], Ftab)

query = paste("SELECT DISTINCT",
    paste(Atable, Afield, sep=".", collapse=", "),
    "FROM")
  
  matchGeneKey=match("gene_id_key",c(Akey,Fkey))
  matchTransKey=match("transcript_id_key",c(Akey,Fkey))

  if(!is.na(matchGeneKey) && !is.na(matchTransKey)){ # so keys are mixed and we have to take order into account
     Akey[Akey =="transcript_id_key"] = "gene_id_key"
     Fkey[Fkey =="transcript_id_key"] = "gene_id_key"
  }
  
  ### end Ensembl specific part #######
  
  tables <- unique(c(Atable, Ftable));
  query = paste(query, tables[1])

  if(length(tables) > 1){    
    query = paste(query, " INNER JOIN (",paste("",tables[-1],sep="",collapse=","), ") ON (",paste(paste(tables[-1],Akey[1],sep="."), paste(tables[1], Akey[1], sep="."),sep=" = ", collapse=" AND "))
     query = paste(query, ")", sep="")
  }
  
  query = paste(query, " WHERE ", paste(Ftable, Ffield[1], sep =".")," IN (",paste("'", values, "'", collapse=", ", sep=""), ")", sep="")

  return(query)
}
