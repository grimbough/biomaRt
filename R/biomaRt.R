packageName <- "biomaRt"

setClass("Mart",
         representation(mysql = "logical",
                        connections = "list",
                        mysqldriver = "list",
                        mainTables = "list",
                        biomart = "character",
                        host = "character",
                        vschema = "character",
                        dataset = "character",
                        filters = "environment",
                        attributes = "environment"
                        ),
         prototype(mysql = FALSE,
                   connections = new("list"),
                   dataset = "",
                   vschema="default"
                   )
         );

setClass("martTable",
         representation(id = "character",	
                        table = "list"	
                        ));	

setMethod("show","martTable",
  function(object){	
    nrShow = 5	
    res = paste("Object of class 'martTable' with", length(object@id), "IDs.")
    if(length(object@id) > nrShow)
      res = paste(res, "The first", nrShow, "rows are:")
    cat(strwrap(res, exdent=5), sep="\n")
    n  = min(nrShow, length(object@id))
    df = do.call("data.frame", args=lapply(object@table, "[", 1:n));
    df = cbind(object@id[1:n],df)    
    show(df)
  })


#######################
#list marts
#######################


listMarts <- function( mart, host, user, password, includeHosts = FALSE, mysql = FALSE){
  
  if(mysql){
    
#MySQL-----------------------------------------------------------    
    require(RMySQL)
    
    if(missing(mart)){
      mart <- c("ensembl","vega","snp","msd","uniprot","sequence","wormbase")
    }
    if(missing(host)){
      host <- c("ensembldb.ensembl.org", "martdb.ebi.ac.uk")
      user <- c("anonymous","anonymous")
      password <- c("","")
    }
    
    database <- NULL
    driv <- dbDriver("MySQL", force.reload = FALSE);
    
    for(i in 1:length(host)){
      connection <- dbConnect(driv, user = user[i], host = host[i], password = password[i]);
  
      res <- dbGetQuery(connection,"show databases like '%mart%'"); 
                                        #Search latest releases of marts
      if(dim(res)[1] >= 1){
        for(j in 1:length(mart)){ 
          matches <- grep(mart[j],res[,1]);
          if(length(matches) > 1){
            version <- 1;
            latest <- 1;
            for(j in 1:length(matches)){
              v <- suppressWarnings(as.numeric(strsplit(res[matches[j],1],"_")[[1]][3]));
              if(!is.na(v)){
                
                if(v > version){
                  latest <- j;
                  version <- v;
                }
              }
            }
            if(!includeHosts){
              database <- c(database,res[matches[latest],1])
            }
            else{
              database <- rbind(database,cbind(res[matches[latest],1],host[i]))
            }
          }
          else{
            if(sum(matches)> 0){
              if(!includeHosts){
                database <- c(database,res[matches,1]);
                v <- as.numeric(strsplit(res[matches,1],"_")[[1]][3]);
              }
              else{
                database <- rbind(database,cbind(res[matches,1],host[i]));
              }
            }
          }
          
          dbDisconnect(connection);
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
    registry = getURL(paste(host,"?type=registry", sep=""))
    registry = xmlTreeParse(registry)
    registry = registry$doc$children[[1]]
    
    marts = list(biomart = NULL, version = NULL, host = NULL, path = NULL)
    index = 1
    
    for(i in 1:xmlSize(registry)){
      if(xmlName(registry[[i]])=="virtualSchema"){
           vschema = xmlGetAttr(registry[[i]],"name")
        for(j in 1:xmlSize(registry[[i]])){
          if(xmlGetAttr(registry[[i]][[j]],"visible") == 1){
            marts$biomart[index] = xmlGetAttr(registry[[i]][[j]],"name")
            marts$version[index] = xmlGetAttr(registry[[i]][[j]],"displayName")
            marts$host[index] = xmlGetAttr(registry[[i]][[j]],"host")
            marts$path[index] = xmlGetAttr(registry[[i]][[j]],"path")
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
    stop("no Mart object to disconnect");
  }
  openConnections <- names(mart@connections);
  
  if(match("biomart",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$biomart );
  }
}

####### local functions #############

mapSpeciesToHomologTable <- function(fromESpecies = NULL, toESpecies = NULL) {
  
  table <- NULL;
  table <- paste(fromESpecies,"_gene_ensembl__homologs_",toESpecies,"__dm",sep="");
  return(table);
}

mapFilter <- function(type){
  
  if(!(type %in% c("affy","entrezgene","agilentprobe","refseq","embl","hugo","ensembl","ensemblTrans","flybase", "unigene"))){
    stop("invalid type choose either affy, refseq, embl, hugo, ensembl, ensemblTrans, unigene, agilentprobe or entrezgene");
  }
  mapf = switch(type,
    ensembl =  "ensembl_gene_id",
    ensemblTrans = "ensembl_transcript_id",
    entrezgene = "entrezgene",
    hugo = "hgnc_symbol",
    refseq = "refseq_dna",
    embl = "embl",
    unigene = "unigene",
    agilentprobe = "agilentprobe"
    );
  return(mapf) 
}



######### public functions ##############

###########################################################################################
#                          
#Ensembl specific functions
##########################################################################################

######################
#get gene information#
######################

getGene <- function( id, type, array, mart){
  
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function useMart")
  }
  
  if(mart@biomart != "ensembl"){
    stop("you can only use ensembl for this query.  Please do: useMart('ensembl').  To access VEGA you have to use the more advanced biomaRt function:  useMart. listDatasets, useDataset, listFilters, listAttributes and getBM.");
  }
  
  if(mart@dataset==""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
  
  if(!missing(array)){
    type <- "affy";
  }
  
  if( missing( type )){
    stop("you must provide the identifier type using the type argument")
  }
  
#MySQL-----------------------------------------------------------
  if(mart@mysql){
    
    if("ensembl" != mart@biomart){
      stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
    }
    
    speciesTable <- unique(mart@mainTables$tables[mart@mainTables$keys == "gene_id_key"]);
    
    IDTable = NULL
    dbcolQID = NULL
    
    if(!missing(array)){
      IDTable = get(array,mart@filters)$table
      dbcolQID = get(array,mart@filters)$field
    }
    else{
      IDTable <- get(mapFilter(type),mart@filters)$table
      dbcolQID = get(mapFilter(type),mart@filters)$field
    }
    
    dbcolID <- get(mapFilter("ensembl"),mart@filters)$field;
    
    if (length( id ) >= 1){
      
      ids <- paste("'",id,"'",sep="",collapse=",");
      
      if(type == "ensembl"){
        query <- paste("select distinct ",dbcolID,", display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end, chrom_strand  from ", speciesTable," where ",dbcolQID," in (",ids,")",sep="");
      }
      else{
        query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".gene_stable_id,",speciesTable,".display_id,",speciesTable,".description, ",speciesTable,".band,",speciesTable,".chr_name,",speciesTable,".gene_chrom_start,",speciesTable,".gene_chrom_end,",speciesTable,".chrom_strand  from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,")",sep="");        
      }
      
      res <- dbGetQuery( conn = mart@connections$biomart, statement = query);
      
      if(dim(res)[1] == 0){
        writeLines("Query retrieved no results")
        return(NULL)
      }
      else{
        mt = match(res[,1], id);
        if(any(is.na(mt)))
          stop("Internal error!");
        if(type == "ensembl"){ 
          res = res[order(mt),c(1,2,3,4,5,6,7,8,1)];
        }
        else{
          res = res[order(mt),c(1,3,4,5,6,7,8,9,2)];
        }
        
        names(res) = c("id","symbol", "description", "band", "chromosome", "start", "end","strand" ,"martID");
        table <- as.data.frame(res)
        return(table);     
      }
    }
  }

#Webservice-----------------------------------------------------
  else{
     
      startpos = "start_position"
      endpos = "end_position"
      chrname="chromosome_name"
      geneid="ensembl_gene_id"
      transid="ensembl_transcript_id"
      strand = "strand"
      symbol = switch(mart@dataset, 
                      hsapiens_gene_ensembl = "hgnc_symbol",
                      mmusculus_gene_ensembl = "mgi_symbol",
                      rnorvegicus_gene_ensembl = "mgi_symbol",
                      scerevisiae_gene_ensembl = "sgd",
                      celegans_gene_ensembl = "external_gene_id",
                      cintestinalis_gene_ensembl = "external_gene_id",
                      ptroglodytes_gene_ensembl = "external_gene_id",  
                      frubripes_gene_ensembl = "external_gene_id",
                      agambiae_gene_ensembl = "external_gene_id",
                      ggallus_gene_ensembl = "external_gene_id",
                      xtropicalis_gene_ensembl = "external_gene_id",
                      drerio_gene_ensembl = "external_gene_id",
                      tnigroviridis_gene_ensembl = "external_gene_id",
                      mmulatta_gene_ensembl = "external_gene_id",
                      mdomesticus_gene_ensembl = "external_gene_id",
                      amellifera_gene_ensembl = "external_gene_id",
                      dmelanogaster_gene_ensembl = "external_gene_id",
                      btaurus_gene_ensembl = "external_gene_id",
                      cfamiliaris_gene_ensembl = "external_gene_id",
                      )
   
    if(!missing(array)){
      if(array=="affy_hg_u133a_2"){
        attrib = "affy_hg_u133a_v2"
      }
      else{
        attrib = array
      }
      table = getBM(attributes=c(attrib,symbol,"description",chrname,"band",strand,startpos,endpos,geneid,transid),filters = array, values = id, mart=mart)
    }
    else{
      if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details")
      filter = mapFilter(type)    
      table = getBM(attributes=c(filter,symbol,"description",chrname,"band",strand,startpos, endpos, geneid,transid),filters = filter, values = id, mart=mart)
    }
    if(!is.null(table)){
      colnames(table)=c("ID","symbol", "description", "chromosome","band","strand","chromosome_start","chromosome_end","ensembl_gene_id","ensembl_transcript_id")
    }
    return(table)
  }
}

########################
#get feature           #
########################


getFeature <- function( symbol, OMIM, OMIMID, GO, GOID, array, chromosome, start, end, type,  mart){
  
  if(!missing(array)){
    type <- "affy"
  }
  
  else{
    if( missing( type )){
      stop("you must provide the identifier type using the type argument");
    }
  }
  
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function useMart");
  }
  
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  
  if(mart@dataset == ""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
  
  if(mart@mysql){
    if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl")  && !(type == "hugo") && !(type == "ensembl")  && !(type == "flybase") ){
      stop("invalid type choose either affy, refseq, embl, hugo or entrezgene");
    }
    
    if( type == "affy" ){
      if( !missing( array ) ){
        affyBool <- paste(array,"_bool",sep="");    
      }
    }
    speciesTable <- unique(mart@mainTables$tables[mart@mainTables$keys == "gene_id_key"]);
    
    IDTable = NULL
    dbcolID=NULL
    if(!missing(array)){
      IDTable = get(array,mart@filters)$table
      dbcolID <- get(mapFilter(array),mart@filters)$field;
    }
    else{
      IDTable <- get(mapFilter(type),mart@filters)$table
      dbcolID <- get(mapFilter(type),mart@filters)$field;
    }
    
    if(!missing(symbol)){
      if(type != "ensembl"){
        query <- paste("select distinct ",speciesTable,".display_id, ",IDTable,".",dbcolID,", description from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",speciesTable,".display_id like '%",symbol,"%' and ",IDTable,".",dbcolID," !='NULL'",sep="");
      }
      else{
        query <- paste("select distinct display_id,",dbcolID,", description from ", speciesTable ," where display_id like '%",symbol,"%' and ",dbcolID," !='NULL'",sep="");
      }
    }
    if(!missing(OMIM)){
      OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
      query <- paste("select distinct ",IDTable,".",dbcolID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",OMIMTable,".disease like '%",OMIM,"%'  and ",IDTable,".",dbcolID," != 'NULL'",sep="");  
    }
    if(!missing(OMIMID)){
      OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
      query <- paste("select distinct ",IDTable,".",dbcolID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",OMIMTable,".omim_id in ('",OMIMID,"') and ",IDTable,".",dbcolID," != 'NULL'",sep="");
    }
    if(!missing(GO)){
      GOTable <- get("go",mart@filters)$table;
      query <- paste("select distinct ",IDTable,".",dbcolID,",",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_id_key = ",GOTable,".gene_id_key where ",GOTable,".description like '%",GO,"%' and ",IDTable,".",dbcolID," != 'NULL'",sep="");
    }
    if(!missing(GOID)){
      GOTable <- get("go",mart@filters)$table;    
      query <- paste("select distinct ",IDTable,".",dbcolID,",",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_id_key = ",GOTable,".gene_id_key where ",GOTable,".dbprimary_id in ('",GOID,"') and ",IDTable,".",dbcolID," != 'NULL'",sep="");
    }
    
    if(!missing(chromosome)){
      if(missing(start) && missing(end)){
        query <- paste("select distinct a.",dbcolID,", b.chr_name, b.gene_chrom_start, b.gene_chrom_end from ", IDTable ," as a inner join ",speciesTable," as b on a.gene_stable_id = b.gene_stable_id where b.chr_name = '",chromosome,"' and a.",dbcolID," !='NULL'",sep="");
        
      }
      else{
        query <- paste("select distinct a.",dbcolID,",b.chr_name,b.gene_chrom_start,b.gene_chrom_end from ", IDTable ," as a inner join ",speciesTable," as b on a.gene_stable_id = b.gene_stable_id where b.chr_name ='",chromosome,"' and b.gene_chrom_start >=",start," and b.gene_chrom_end <= ",end," and a.",dbcolID," !='NULL'",sep="");  
        
      }
    }
    
    res <- dbGetQuery( conn = mart@connections$biomart,statement = query);
    
    if(dim(res)[1] != 0){
      if(!missing(symbol)){
        table <- new("martTable", id = as.vector(as.character(res[,2])), table = list(symbol = as.vector(res[,1]), description = as.vector(res[,3])))
      }
      if(!missing(OMIM) || !missing(OMIMID)){
        table <- new("martTable", id = as.vector(as.character(res[,1])), table = list(OMIMID = as.vector(res[,2]), description = as.vector(res[,3])))
      }
      if(!missing(GO) || !missing(GOID)){
        table <- new("martTable", id = as.vector(as.character(res[,1])), table = list(GOID = as.vector(res[,2]), description = as.vector(res[,3]), evidence = as.vector(res[,4])))
        
      } 
      if(!missing(chromosome)){
        table <- new("martTable", id = as.vector(as.character(res[,1])), table = list(chromosome = as.vector(res[,2]), start = as.vector(res[,3]), end = as.vector(res[,4])))
      }
    }
    else{
      writeLines("No match found")
    }
    return( table )
  }
  
#-----webservice------------------------------------
  
  else{
    
    filter=NULL
    attribute=NULL
    values = NULL
    
    startpos = "start"
    endpos = "end"
    chrname="chromosome_name"
    geneid="ensembl_gene_id"
    transid="ensembl_transcript_id"
    strand = "strand"
    attribute = switch(type, hugo="hgnc_symbol",agilentcgh = "agilent_cgh",ensemblTrans="ensembl_transcript_id",agilentprobe="agilent_probe", entrezgene = "entrezgene", locuslink = "entrezgene", embl = "embl", refseq ="refseq_dna", unigene="unigene", affy = array, ensembl="ensembl_gene_id")
    
    if(!missing(symbol)){
     
      filter = switch(mart@dataset, 
                      hsapiens_gene_ensembl = "hgnc_symbol",
                      mmusculus_gene_ensembl = "mgi_symbol",
                      rnorvegicus_gene_ensembl = "mgi_symbol",
                      scerevisiae_gene_ensembl = "sgd",
                      celegans_gene_ensembl = "external_gene_id",
                      cintestinalis_gene_ensembl = "external_gene_id",
                      ptroglodytes_gene_ensembl = "external_gene_id",  
                      frubripes_gene_ensembl = "external_gene_id",
                      agambiae_gene_ensembl = "external_gene_id",
                      ggallus_gene_ensembl = "external_gene_id",
                      xtropicalis_gene_ensembl = "external_gene_id",
                      drerio_gene_ensembl = "external_gene_id",
                      tnigroviridis_gene_ensembl = "external_gene_id",
                      mmulatta_gene_ensembl = "external_gene_id",
                      mdomesticus_gene_ensembl = "external_gene_id",
                      amellifera_gene_ensembl = "external_gene_id",
                      dmelanogaster_gene_ensembl = "external_gene_id",
                      btaurus_gene_ensembl = "external_gene_id",
                      cfamiliaris_gene_ensembl = "external_gene_id",
                      )
 
      attribute = c(filter,attribute)
      values = symbol
    }
    
    if(!missing(OMIM)){
      stop("getFeature currently can only query for omim descriptions using mysql access.  create a mart object using useMart('ensembl',mysql=TRUE)")
    }
    if(!missing(OMIMID)){
      stop("getFeature currently can only query for omim descriptions using mysql access.  create a mart object using useMart('ensembl',mysql=TRUE)")
    }
    if(!missing(GO)){
      stop("getFeature currently can only query for go descriptions using mysql access.  create a mart object using useMart('ensembl',mysql=TRUE)")
    }
    if(!missing(GOID)){
      filter="go"
      attribute = c("go",attribute)
      values = GOID
    }
    
    if(!missing(chromosome)){
      if(missing(start) && missing(end)){
        filter = "chromosome_name"
        attribute = c(transid,chrname,attribute)
        values = chromosome
      }
      else{
        filter=c(chrname,startpos,endpos)
        attribute = c(transid,chrname,"start_position","end_position",attribute)
        values = list(chromosome, start, end)
      }
    }
    table = getBM(attributes = attribute, filters = filter, values = values, mart=mart)  
    return(table)
  }
}

###################
#get GO annotation#
###################


getGO <- function( id, type, array, mart){
  
  table <- NULL;
  go <- NULL;
  
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function useMart")
  }
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(mart@dataset==""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
 
  if(mart@mysql){
    
    if(!missing(array)){
      type <- "affy";
    }
    else{
      if(missing( type )){
        stop("you must provide the identifier type using the type argument")
      }
    }
    GOTable <- get("go",mart@filters)$table;
    IDTable = NULL
    dbcolQID=NULL
    if(!missing(array)){
      IDTable = get(array,mart@filters)$table
      dbcolQID <- get(array,mart@filters)$field;
    }
    else{
      IDTable <- get(mapFilter(type),mart@filters)$table
      dbcolQID <- get(mapFilter(type),mart@filters)$field;
    }
    
    dbcolID <- get(mapFilter("ensembl"),mart@filters)$field;        #get database id col   
    
    if( length( id ) >= 1){
      ids <- paste("'",id,"'",sep="",collapse=",")
      res <- NULL
      if(type == "ensembl"){
        query <- paste("select distinct ",IDTable,".",dbcolQID,",",IDTable,".",dbcolQID,", ",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
      }
      else{
        query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolID,",",GOTable,".dbprimary_id, description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
      }
      res <- dbGetQuery(conn = mart@connections$biomart, statement = query);
      
      if(dim(res)[1] == 0){
        table <- new("martTable", id = id, table = list(GOID = NA, description = NA, evidence = NA, martID = NA))
        
      }
      else{
        mt = match(res[,1], id)
        if(any(is.na(mt)))
          stop("Internal error!")
        if(type == "ensembl"){ 
          res = res[order(mt),c(1,3,4,5,1)];
        }
        else{
          res = res[order(mt),c(1,3,4,5,2)];
        }
        names(res) = c("id","GOID", "description", "evidence", "martID")
        table <- as.data.frame(res)
      }
    }
    return( table );
  }
#--webservice----------------------
  else{
      goid="go"
      geneid="ensembl_gene_id"
      transid="ensembl_transcript_id"
   
    if(!missing(array)){
      if(array=="affy_hg_u133a_2"){
        attrib = "affy_hg_u133a_v2"
      }
      else{
        attrib = array
      }
       
      table = getBM(attributes=c(attrib,goid, "go_description", "evidence_code",geneid,transid),filters = array, values = id, mart=mart)
    }
    else{
      if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details")
      filter = mapFilter(type)
      table = getBM(attributes=c(filter,goid, "go_description", geneid,transid),filters = filter, values = id, mart=mart)
    }
    if(!is.null(table)){
      if(!missing(array)){
        colnames(table)=c("ID","go_id", "go_description", "evidence_code","ensembl_gene_id","ensembl_transcript_id")
      }
      else{
          colnames(table)=c("ID","go_id", "go_description","ensembl_gene_id","ensembl_transcript_id")
      }
    }
    return(table)
  }
}

#####################
#get OMIM annotation#
#####################


getOMIM <- function( id, type, array, mart){
  
  if( missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function useMart")
  }
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(mart@dataset==""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
  
  if(mart@mysql){
    OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
    omim <- NULL;
    table <- NULL;
    species <- "hsapiens";
    if(!missing(array)){
      type <- "affy";
    }
    else{
      if( missing( type )){
        stop("you must provide the identifier type using the type argument")
      }
    }
    
    IDTable = NULL
    dbcolQID = NULL
    if(!missing(array)){
      IDTable = get(array,mart@filters)$table
      dbcolQID = get(array,mart@filters)$field
    }
    else{
      IDTable <- get(mapFilter(type),mart@filters)$table
      dbcolQID = get(mapFilter(type),mart@filters)$field
    }
    
    dbcolID <- get(mapFilter("ensembl"),mart@filters)$field;
    
    if( length( id ) >= 1){
      ids <- paste("'",id,"'",sep="",collapse=",")
      res <- NULL
      if(type == "ensembl"){
        query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolQID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".",dbcolQID," in (",ids,")  and ",OMIMTable,".omim_id != 'NULL'",sep="");
        
      }
      else{
        query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".",dbcolQID," in (",ids,") and ",OMIMTable,".omim_id != 'NULL'",sep="");
      }
      
      res <- dbGetQuery(conn = mart@connections$biomart, statement = query);
      
      if(dim(res)[1] == 0){
        table <- new("martTable", id = id, table = list(OMIMID = NA, disease = NA, martID = NA))
      } 
      else{
        mt = match(res[,1], id)
        if(any(is.na(mt)))
          stop("Internal error!")
        if(type == "ensembl"){ 
          res = res[order(mt),c(1,3,4,1)];
        }
        else{
          res = res[order(mt),c(1,3,4,2)];
        }
        names(res) = c("id","OMIMID", "disease", "martID")
        table <- as.data.frame(res)
      }  
    }
    return( table );
  }
  
#------webservice-------------------------------
  else{
      omimid="omim"
      geneid="ensembl_gene_id"
      transid="ensembl_transcript_id"
   
    if(!missing(array)){
      if(array=="affy_hg_u133a_2"){
        attrib = "affy_hg_u133a_v2"
      }
      else{
        attrib = array
      }
      table = getBM(attributes=c(attrib,omimid, "disease_description",geneid,transid),filters = array, values = id, mart=mart)
    }
    else{
      if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details")
      filter = mapFilter(type)
      table = getBM(attributes=c(filter,omimid, "disease_description",geneid,transid),filters = filter, values = id, mart=mart)
    }
   
    if(!is.null(table)){
      colnames(table)=c("ID","omim_id", "description", "ensembl_gene_id","ensembl_transcript_id")
    }
    return(table)
  }
}

#####################
#getSequence
#####################


getSequence <- function(chromosome, start, end, id, type, seqType, mart){

  if(missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function useMart")
  }
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(mart@dataset==""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }  
  if(mart@mysql){
  
    martdb = "";
    if(mart@biomart != "sequence"){
      version = "0"
      marts = listMarts(mysql = TRUE)
      for(i in 1:length(marts)){
        if("sequence" == strsplit(marts[i],"_")[[1]][1]){
          version = strsplit(marts[i],"_")[[1]][3]
        }
      }
      martdb = paste("sequence_mart_",version,sep="")
      mart@biomart="sequence"
      
    }
    mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = martdb, password = "")
    
    species = strsplit(mart@dataset,"_")[[1]][1]
    sequence <- NULL;  
    speciesTable <- paste( species,"_genomic_sequence__dna_chunks__main",sep="" ); 
    
    if(missing( martTable ) && !missing( chromosome ) && !missing( start ) && !missing( end )){
      for(i in 1:length( chromosome )){
        
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
    
    else{
      if(!missing( martTable )){
                                        #check class!!
        chromosome = martTable@table$chromosome;
        start = martTable@table$start;
        end = martTable@table$end;
        
        for(i in 1:length( chromosome )){
          
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
    }
    
    table <- new("martTable", id = paste(chromosome, start, end, sep = "_") ,table = list(chromosome = chromosome, start = start, end = end, sequence = sequence))
    
    return(table)
  }
  else{
    geneid="gene_stable_id"
    
    species = strsplit(mart@dataset,"_")[[1]][1]
    if(!missing(chromosome)){
      if(seqType %in% c("cdna","peptide","3utr","5utr")){
        query = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0' ><Dataset name = '",species,"_gene_ensembl'><ValueFilter name = 'chromosome_name' value = '",chromosome,"'/><ValueFilter name = 'start' value = '",start,"'/><ValueFilter name = 'end' value = '",end,"'/></Dataset><Links source = '",species,"_gene_ensembl' target = '",species,"_gene_ensembl_structure' defaultLink = '",species,"_internal_transcript_id' /><Dataset name = '",species,"_gene_ensembl_structure'><Attribute name = '",geneid,"'/><Attribute name = 'str_chrom_name'/><Attribute name = 'biotype'/></Dataset><Links source = '",species,"_gene_ensembl_structure' target = '",species,"_genomic_sequence' defaultLink = '",seqType,"' /><Dataset name = '",species,"_genomic_sequence'><Attribute name = '",seqType,"'/></Dataset></Query>",sep="")
      }
      else{
        stop("The type of sequence specified with seqType is not available. Please select from: cdna, peptide, 3utr, 5utr")
      }
    }
    
    if(!missing(id)){
      valuesString = paste(id,"",collapse=",",sep="")
      if(seqType %in% c("cdna","peptide","3utr","5utr")){
        query = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0' ><Dataset name = '",species,"_gene_ensembl'><ValueFilter name = '",mapFilter(type),"' value = '",valuesString,"'/></Dataset><Links source = '",mart@dataset,"' target = '",species,"_gene_ensembl_structure' defaultLink = '",species,"_internal_transcript_id' /><Dataset name = '",species,"_gene_ensembl_structure'><Attribute name = '",geneid,"'/><Attribute name = 'str_chrom_name'/><Attribute name = 'biotype'/></Dataset><Links source = '",species,"_gene_ensembl_structure' target = '",species,"_genomic_sequence' defaultLink = '",seqType,"' /><Dataset name = '",species,"_genomic_sequence'><Attribute name = '",seqType,"'/></Dataset></Query>", sep="")
      }
      else{
        stop("The type of sequence specified with seqType is not available. Please select from: cdna, peptide, 3utr, 5utr")
      }   
    }
    postRes = postForm(paste(mart@host,"?",sep=""),"query"=query)
    
    if(postRes != ""){
      ## convert the serialized table into a dataframe
      con = textConnection(postRes)
      result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
      close(con)
      ## check and postprocess
      if(all(is.na(result[,ncol(result)])))
        result = result[,-ncol(result),drop=FALSE]
      #stopifnot(ncol(result)==length(attributes))
      #if(class(result) == "data.frame"){
      #  colnames(result) = attributes
      #}
      
    } else {
      warning("getBM returns NULL.")
      result=NULL
    }
    return(result)
  }
}

################
#show affy info#
################

getAffyArrays <- function(mart){
  if(missing(mart)){
    stop("no Mart object found, create this first using the function useMart and then give as argument in the function")
  }
  if(is.null(mart@dataset)){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
  att=listFilters(mart)
  affy=att[grep("affy",att[,1]),]
  affy=affy[-grep("with",affy[,1]),]
  return(affy)
}

#######################
#getSNP
#######################

getSNP <- function(chromosome, start, end, mart){
  if( missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function useMart")
  }
  if("snp" != mart@biomart){
    stop("This function only works when using to snp. To use this function use: mart =  useMart('snp')")
  }
  if(missing(chromosome) || missing(start) || missing(end) ){
    stop("you have to give chromosome, start and end positions as arguments, see ?getSNP for more information")
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
    attributes = c(tscid,"refsnp_id","allele","chrom_start","chrom_strand")
    table = getBM(attributes = attributes,filters = c("chr_name",snpstart,snpend), values = list(chromosome, start, end), mart=mart)    
    return(table)
  }
}

##################
#getHomolog      #
##################

getHomolog <- function(id, from.type, to.type, from.array, to.array, from.mart, to.mart) {
  
  if( missing( to.mart )|| class( to.mart ) != 'Mart' || missing( from.mart )|| class( from.mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function useMart")
  }
  
  if("ensembl" != to.mart@biomart || "ensembl" != from.mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  
  if(from.mart@mysql){
    if(!to.mart@mysql)stop("Both mart object should both be using mysql or not.  Your from.mart uses mysql but your to.mart not.")
    
    if( !missing( id )){
      id <- as.character(id);
    }
    else{
      stop("No id's to search for homologs given");
    }
    
    if( !missing( from.array )){
      from.type <- "affy";
    }
    else{
      if ( missing( from.type )) {
        stop("You must provide the identifier type using the from.type argument")
      }
    }
    
    if( !missing( to.array )){
      to.type <- "affy";
    }
    else{
      if ( missing( to.type )) {
        stop("You must provide the identifier type using the to.type  argument")
      } 
    }
    
    fromCol = NULL
    toCol = NULL
    fromIDTable = NULL
    toIDTable = NULL
    
    if(missing(from.array)){
      fromIDTable <- get(mapFilter(from.type), from.mart@filters)$table
      fromCol <- get(mapFilter(from.type), from.mart@filters)$field;
    }
    else{
      fromIDTable <- get(from.array, from.mart@filters)$table
      fromCol <- get(from.array, from.mart@filters)$field;
    }
    if(missing(to.array)){  
      toIDTable <- get(mapFilter(to.type), to.mart@filters)$table
      toCol <- get(mapFilter(to.type), to.mart@filters)$field;
    }
    else{
      toIDTable <- get(to.array, to.mart@filters)$table
      toCol <- get(to.array, to.mart@filters)$field;
    }
    
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
        
        for (j in 1:length(id)) {
          
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
    if(to.mart@mysql)stop("The mart objects should either use both mysql or not. Here your from.mart does not use mysql but you to.mart does.")
    
     # to.attributes=c("ensembl_gene_id","ensembl_transcript_id")
    to.attributes=NULL
    from.attributes = NULL
    filter = NULL
    
    if(!missing(to.array)){
      to.attributes = c(to.attributes,to.array)
    }
    else{
        to.attributes = c(to.attributes,mapFilter(to.type))
    }
    if(!missing(from.array)){
      from.attributes = c(from.attributes,from.array)
      filter = from.array
    }
    else{
      from.attributes="ensembl_gene_id"
      #from.attributes = c(from.attributes,mapFilter(from.type))
      filter = mapFilter(from.type)
    }
    
    xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0'> <Dataset name = '",from.mart@dataset,"'>",sep="")
    attributeXML = ""#paste("<Attribute name = '", from.attributes, "'/>", collapse="", sep="")
    valuesString = paste(id,"",collapse=",",sep="")
    filterXML = paste("<Filter name = '",filter,"' value = '",valuesString,"' />", sep="")
    xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset>",sep="")
   
    species = strsplit(to.mart@dataset,"_")[[1]][1]
    xmlQuery = paste(xmlQuery, "<Links source = '",from.mart@dataset,"' target = '",to.mart@dataset,"' defaultLink = '",species,"_internal_gene_id'/><Dataset name = '",to.mart@dataset,"' >", sep="")
    to.attributeXML =  paste("<Attribute name = '", to.attributes, "'/>", collapse="", sep="") 
    xmlQuery = paste(xmlQuery, to.attributeXML,"</Dataset></Query>",sep="") 
    postRes = postForm(paste(to.mart@host,"?",sep=""),"query"=xmlQuery)
    
    if(postRes != ""){
      ## convert the serialized table into a dataframe
      con = textConnection(postRes)
      result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
      close(con)
      ## check and postprocess
      if(all(is.na(result[,ncol(result)])))
        result = result[,-ncol(result),drop=FALSE]
      #stopifnot(ncol(result)==length(attributes))
      #if(class(result) == "data.frame"){
      #  colnames(result) = attributes
      #}
      
    } else {
      warning("getBM returns NULL.")
      result=NULL
    }
    return(result)
  }
} 


###########################
#Xref functions           #
###########################
#
#Should become getBM

getPossibleXrefs <-  function( mart ) {
  if ( missing( mart ) || class( mart ) != 'Mart') {
    stop('You must supply a valid mart object.  You can use useMart to produce one')
  }
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(!mart@mysql)stop("This function only works via MySQL access to the BioMart.  Create a new Mart object using the function useMart('ensembl',mysql=TRUE)")
  
  res <- dbGetQuery( mart@connections$biomart, "show tables like '%_gene_ensembl__xref%'")
  xref <- strsplit(res[,1],'_gene_ensembl__xref_')
  xref <- lapply(xref,function(s) {
    s[2] <- gsub('__dm','',s[2])
    return(s)
  })
  xrefdf <- do.call('rbind',xref)
  colnames(xrefdf) <- c('species','xref')
  return(xrefdf)
} 

getXref <- function( id, from.species, to.species, from.xref, to.xref, mart) {
  
  if ( missing( mart ) || class( mart )!='Mart'){
    stop('You must specify a valid Mart object which can be created using useMart.');
  }
  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(!mart@mysql)stop("This function only works via MySQL access to the BioMart.  Create a new Mart object using the function useMart('ensembl',mysql=TRUE)")
  if ( missing( id )) {
    stop('You need to give ID(s) to map from');
  }
  if ( missing( from.species )) {
    writeLines('fromspecies is a necessary argument\nValid species for this mart are:');
    #print( getSpecies(mart = mart, db = "ensembl" ));
    stop();
  }
  if ( missing( from.xref ) || missing( to.xref ) ) {
    writeLines('Both fromxref and toxref are required.\nPossible crossreferences are:');
    print(getPossibleXrefs(mart = mart));
    stop();
  }
  
  xp <- paste(id, collapse="','")
  
  if (missing( to.species )) {
    query <- paste( 
                   'SELECT distinct a.dbprimary_id as fromid,b.dbprimary_id as toid,',
                   'a.gene_stable_id as ensemblgene ',
                   'from ',from.species,'_gene_ensembl__xref_',from.xref,'__dm as a ',
                   'left join ',from.species,'_gene_ensembl__xref_',to.xref,'__dm as b ',
                   "on a.gene_id_key=b.gene_id_key where a.dbprimary_id in ('",xp,
                   "') and b.dbprimary_id != 'NULL'",sep="");
  } else {
    if ( from.species == to.species) {
     query <- paste( 
                   'SELECT distinct a.dbprimary_id as fromid,b.dbprimary_id as toid,',
                   'a.gene_stable_id as ensemblgene ',
                   'from ',from.species,'_gene_ensembl__xref_',from.xref,'__dm as a ',
                   'left join ',from.species,'_gene_ensembl__xref_',to.xref,'__dm as b ',
                   "on a.gene_id_key=b.gene_id_key where a.dbprimary_id in ('",xp,
                   "') and b.dbprimary_id != 'NULL'",sep="");
   
    }
    else{
      query <- paste(
                     'SELECT distinct a.dbprimary_id as fromid,b.dbprimary_id as toid,c.* ',
                     'from ',from.species,'_gene_ensembl__xref_',from.xref,'__dm as a ',
                     'left join ',from.species,'_gene_ensembl__homologs_',to.species,'__dm as c ',
                     'on a.gene_id_key=c.gene_id_key ',
                     'left join ',to.species,'_gene_ensembl__xref_',to.xref,'__dm as b ',
                     "on c.homol_id=b.gene_id_key where a.dbprimary_id in ('",xp,
                     "') and b.dbprimary_id != 'NULL'",sep="");
    }
  }
  
  res <- dbGetQuery(con = mart@connections$biomart, query);
  
  if (dim(res)[1]>0){

    table <- new("martTable", id = res[,1], table = list(from.id = res[,1], to.id = res[,2], martID = res[,3]));

  }
  
  else{
    table <- new("martTable", id = id, table = list(from.id = NA, to.id = NA, martID = NA));
  }
  return( table )
}


####################
#export FASTA      #
####################

exportFASTA <- function( martTable, file ){
  if( missing( martTable ) || class( martTable ) != "martTable"){
    stop("No martTable given to write");
  }
  if( missing(file)){
    stop("Please provide filename to write to");
  }
  for(i in 1:length(martTable@id)){
    cat(paste(">",martTable@id[i],"\n",sep=""),file = file, append=TRUE);
    cat(martTable@table$sequence[i],file = file, append = TRUE);
    cat("\n\n", file = file, append = TRUE);
  }  
}

#####################
#INTERPRO functions #
#####################


getINTERPRO <- function( id, type, array, mart){
     
  if( missing( mart )|| class(mart)!='Mart'){
    stop("you must provide a valid Mart object, create with function useMart")
  }

  if("ensembl" != mart@biomart){
    stop("This function only works when using to ensembl. To use this function use: mart =  useMart('ensembl')")
  }
  if(mart@dataset==""){
    stop("Please select a dataset first.  You an see the available datasets by:  listDatasets('ensembl').  Then you should create the mart object with e.g.: mart = useMart(biomart='ensembl',dataset='hsapiens_gene_ensembl')");
  }
  
  if(mart@mysql){ 
    if( !missing( array )){
      type <- "affy";
    }
    else{
      if( missing( type )){
        stop("you must provide the identifier type using the type argument");
      }
    }
    
    uniprotTable = get("interpro_ids",mart@filters)$table
    
    IDTable = NULL
    dbcolID = NULL
    if(!missing(array)){
      IDTable = get(array,mart@filters)$table
      dbcolID = get(array,mart@filters)$field
    }
    else{
      IDTable <- get(mapFilter(type),mart@filters)$table
      dbcolID = get(mapFilter(type),mart@filters)$field
    }
    
    if (length( id ) >= 1){
      
      ids <- paste("'",id,"'",sep="",collapse=",")
      if(type == "ensembl"){
        query <- paste("select distinct ",IDTable,".gene_stable_id,",uniprotTable,".interpro_list, short_description, ",uniprotTable,".description from ", IDTable ," inner join ",uniprotTable," on ",IDTable,".gene_stable_id = ",uniprotTable,".gene_stable_id where ",IDTable,".gene_stable_id in (",ids,") and ",uniprotTable,".interpro_list != 'NULL'",sep="");
      }
      
      else{
        if(type == "ensemblTrans"){
          query <- paste("select distinct ",uniprotTable,".transcript_stable_id,",uniprotTable,".interpro_list, short_description, ",uniprotTable,".description from ", uniprotTable," where ",uniprotTable,".transcript_stable_id in (",ids,") and ",uniprotTable,".interpro_list != 'NULL'",sep="");
          
        }
        else{
          query <- paste("select distinct ",IDTable,".",dbcolID,",",uniprotTable,".interpro_list, short_description,",uniprotTable,".description from ", IDTable ," inner join ",uniprotTable," on ",IDTable,".gene_stable_id = ",uniprotTable,".gene_stable_id where ",IDTable,".",dbcolID," in (",ids,") and ",uniprotTable,".interpro_list != 'NULL'",sep="");
        }
      }
      
      res <- dbGetQuery( conn = mart@connections$biomart,statement = query);
      
      if(dim(res)[1] == 0){
       return(NULL)
      }
      else{
        mt = match(res[,1], id)
        if(any(is.na(mt)))
          stop("Internal error!")
        if(type == "ensembl"){ 
          res = res[order(mt),];
        }
        else{
          res = res[order(mt),];
        }
        names(res) = c("id","INTEPROID", "short description", "description")
        table <- as.data.frame(res)
      }
    }
    return(table)
  }

#--webservice--------------------------------
  
  else{
      interproid="interpro"
      geneid="ensembl_gene_id"
      transid="ensembl_transcript_id"
    
    if(!missing(array)){
      if(array=="affy_hg_u133a_2"){
        attrib = "affy_hg_u133a_v2"
      }
      else{
        attrib = array
      }

      table = getBM(attributes=c(attrib,interproid, "interpro_description",geneid,transid),filters = array, values = id, mart=mart)
    }
    else{
      if(missing(type))stop("Specify the type of identifier you are using, see ?getGene for details")
      filter = mapFilter(type)
          table = getBM(attributes=c(filter,interproid, "interpro_description",geneid,transid),filters = filter, values = id, mart=mart)
    }
    if(!is.null(table)){
      colnames(table)=c("ID","interpro_id", "description", "ensembl_gene_id","ensembl_transcript_id")
    }
    return(table)
  }
}

######################################################################################
#
#BioMart functions
#based on martshell
#
######################################################################################

useMart <- function(biomart, dataset, host, user, password, local = FALSE, mysql = FALSE){

  if(mysql){
    
    require(RMySQL)
    driver <- dbDriver("MySQL", force.reload = FALSE);
    
    mart <- new("Mart", biomart = biomart, mysqldriver = list(driver=driver), mysql = TRUE)

    if(! (is.character(biomart) && (length(biomart)==1)))
      stop("'biomart' should be a single character string.")

    if(local){
      if(!missing(host) && !missing(user) && !missing(password)){
        database <- listMarts(mart = biomart, host = host, user = user, password = password);
        mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = user, host = host, dbname = database, password = password)
        writeLines(paste("connected to: ",database[1,1]))
      }
      else{
        stop(sprintf("Please provide host, user and password for using local database '%s'.", biomart))
      }
    }
    else {
      
      version = "0"
      marts = listMarts(mysql = TRUE)
      for(i in 1:length(marts)){
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

      if(!martdb %in% marts) stop("Requested BioMart database is not available please use the function listMarts(mysql=TRUE) to see the valid biomart names you can query using mysql access")
      mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver$driver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = martdb, password = "")
      writeLines(paste("connected to: ",biomart))
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
    
    mart <- new("Mart", biomart = biomart,vschema = marts$vschema[mindex], host = paste("http://",marts$host[mindex],marts$path[mindex],sep=""), mysql= FALSE)
    if(!missing(dataset)){
      mart = useDataset(mart = mart, dataset=dataset)
    }
    return(mart)
  }
}


listDatasets <- function( mart ){
  if(missing( mart ) || class( mart )!='Mart') stop("No Mart object given or object not of class 'Mart'")
  if(mart@mysql){
    
    res <- dbGetQuery(mart@connections$biomart,"select dataset, version from meta_conf__dataset__main where visible = 1")
    return(res)
    
  }
  else{
    
    datasetsTemp = scan(paste(mart@host,"?type=datasets&mart=",mart@biomart,sep=""),
      sep="\t", blank.lines.skip=TRUE, what="character", quiet=TRUE)
    dataset = NULL
    version = NULL
    index = 0
    for(i in 1:(length(datasetsTemp)-3)){
      if(datasetsTemp[i] == "TableSet"){
        if(datasetsTemp[i+3] == "1"){
          index=index+1;
          dataset[index]=datasetsTemp[i+1]
          version[index]=datasetsTemp[i+4]
        }
      }
    }
    datasets=data.frame(dataset=dataset, version=version)
    return(datasets)
  }
}

##########################################
#
#Configuration Parsers
##########################################


parseAttributes <- function(xml, env){
  
   if(xmlName(xml) == "AttributePage"){
     if(!is.null(xmlGetAttr(xml,"hidden"))){
       if(xmlGetAttr(xml,"hidden") != "true"){
         xmlSApply(xml,"parseAttributes",env)
       }
     }
     else{
       xmlSApply(xml,"parseAttributes",env)
     }
   }
   if(xmlName(xml) == "AttributeDescription"){
     if(!is.null(xmlGetAttr(xml,"hidden"))){
       if(xmlGetAttr(xml,"hidden") != "true"){
         if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
           assign(xmlGetAttr(xml,"internalName"),list(description = xmlGetAttr(xml,"displayName"), table = xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
         }
       }
     }
     else{
       if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
         assign(xmlGetAttr(xml,"internalName"),list(description = xmlGetAttr(xml,"displayName"), table = xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
       }
     }
   }
   if(xmlName(xml)=="AttributeCollection"){
     if(!is.null(xmlGetAttr(xml,"hidden"))){
       if(xmlGetAttr(xml,"hidden")!="true"){
          xmlSApply(xml,"parseAttributes",env)
       }
     }
     else{
       xmlSApply(xml,"parseAttributes",env)
     }
   }
   
   if(xmlName(xml)=="AttributeGroup"){
     if(!is.null(xmlGetAttr(xml,"hidden"))){
       if(xmlGetAttr(xml,"hidden")!="true"){
         xmlSApply(xml,"parseAttributes",env)
       }
     }
     else{
       xmlSApply(xml,"parseAttributes",env)
     }
   }
   if(xmlName(xml)=="DatasetConfig"){
         xmlSApply(xml,"parseAttributes",env)
   }
   return()
}

parseFilters <- function(xml, env){
  
  if(xmlName(xml) == "FilterPage"){
    if(!is.null(xmlGetAttr(xml,"hidden"))){
      if(xmlGetAttr(xml,"hidden")!="true"){
        if(!is.null(xmlGetAttr(xml,"displayName"))){
          if(xmlGetAttr(xml,"displayName")=="FILTERS"){
            xmlSApply(xml,"parseFilters",env)
          }
        }
      }
    }
    else{
      if(!is.null(xmlGetAttr(xml,"displayName"))){
        if(xmlGetAttr(xml,"displayName")=="FILTERS"){
          xmlSApply(xml,"parseFilters",env)
        }
      }
    }
  }
  if(xmlName(xml)=="FilterGroup"){
    if(!is.null(xmlGetAttr(xml,"hidden"))){
      if(xmlGetAttr(xml,"hidden")!="true"){
        xmlSApply(xml,"parseFilters",env)
      }
    }
    else{
      xmlSApply(xml,"parseFilters",env)
    }
  }
  if(xmlName(xml)=="FilterCollection"){
    if(!is.null(xmlGetAttr(xml,"hidden"))){
      if(xmlGetAttr(xml,"hidden")!="true"){
        xmlSApply(xml,"parseFilters",env)
      }
    }
    else{
      xmlSApply(xml,"parseFilters",env)
    }
  }
  if(xmlName(xml) == "FilterDescription"){
    if(!is.null(xmlGetAttr(xml,"hidden"))){
      if(xmlGetAttr(xml,"hidden") != "true"){
        if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
          assign(xmlGetAttr(xml,"internalName"),list(descr= xmlGetAttr(xml,"displayName"), table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
        }
        else{
          xmlSApply(xml,"parseFilters",env)
        }
      }
    }
    else{
      if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
        assign(xmlGetAttr(xml,"internalName"),list(descr= xmlGetAttr(xml,"displayName"), table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
      }
      else{
        xmlSApply(xml,"parseFilters",env)
      }
    }
  }
  
  if(xmlName(xml) == "Option"){
    if(!is.null(xmlGetAttr(xml,"hidden"))){
      if(xmlGetAttr(xml,"hidden") != "true"){
        if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
          assign(xmlGetAttr(xml,"internalName"),list(descr= xmlGetAttr(xml,"displayName"), table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
        }
      }
    }
    else{
      if(!is.null(xmlGetAttr(xml,"tableConstraint"))){
        assign(xmlGetAttr(xml,"internalName"),list(descr= xmlGetAttr(xml,"displayName"), table= xmlGetAttr(xml,"tableConstraint"), key = xmlGetAttr(xml,"key"), field= xmlGetAttr(xml,"field")),env = env)
      }
    }
  }
  
  if(xmlName(xml)=="DatasetConfig"){
    xmlSApply(xml,"parseFilters",env)
  }
  return()  
}


useDataset <- function(dataset, mart){

  if(missing(mart) || class(mart)!="Mart") stop("No valid Mart object given, specify a Mart object with the attribute mart")
  if(missing(dataset)) stop("No valid dataset to use given")
  validDatasets=listDatasets(mart)
  if(is.na(match(dataset, validDatasets$dataset)))stop(paste("The given dataset: ",dataset,", is not valid.  Correct dataset names can be obtained with the listDatasets function"))
  
  if(mart@mysql){
    res <- dbGetQuery(mart@connections$biomart,paste("select xml from meta_conf__dataset__main inner join meta_conf__xml__dm on meta_conf__dataset__main.dataset_id_key = meta_conf__xml__dm.dataset_id_key where dataset = '",dataset,"'",sep=""))
    writeLines(paste("Reading database configuration of:",dataset))
    if(dim(res)[1] == 0) stop("This dataset is not accessible from biomaRt as not xml description of dataset is available")

    xml <- xmlTreeParse(res[1,])
    xml <- xml$doc$children[[1]]
   
    filtersEnv = new.env()
    attributesEnv = new.env()

    parseAttributes(xml, attributesEnv)
    parseFilters(xml, filtersEnv)
    mainTables <- getMainTables(xml)

    writeLines("Checking attributes and filters ...", sep=" ")
    mart@attributes <- attributesEnv
    mart@filters <- filtersEnv
    writeLines("ok")
    mart@dataset = dataset
    mart@mainTables <- mainTables
    
    return(mart)
  }
  else{

    config = getURL(paste(mart@host,"?type=configuration&dataset=",dataset,"&virtualschema=",mart@vschema, sep=""))
    config = xmlTreeParse(config)
    config = config$doc$children[[1]]

    filtersEnv = new.env()
    attributesEnv = new.env()
    
    writeLines("Checking attributes and filters ...", sep=" ")
    parseAttributes(config, attributesEnv)
    parseFilters(config, filtersEnv)
    writeLines("ok")
    
    mart@dataset = dataset
    mart@attributes = attributesEnv
    mart@filters = filtersEnv
    
    return( mart )
  }
}

getMainTables <- function( xml ){
  
  writeLines("Checking main tables ...", sep=" ")
  names <- names(xml)
  
  tableM <- NULL
  keyM <- NULL
  i<-1
  j<-1
  for(h in 1:(xmlSize(xml)-1)){
    if(names[h] == "MainTable"){
      tableM[i] = xmlValue(xml[[h]])
      i <- i+1
    }
    if(names[h] == "Key"){
      keyM[j] = xmlValue(xml[[h]])
      j <- j+1
    }
  }
  writeLines("ok")
 
  return(list(tables=tableM,keys=keyM))
}

listAttributes <- function( mart ){
  if(missing( mart ) || class( mart )!='Mart') stop("No Mart object given or object not of class 'Mart'")
  attribList=mget(ls(mart@attributes), env=mart@attributes)
  mat <- do.call(cbind, attribList)
  frame=data.frame(cbind(colnames(mat),mat[1,]), row.names=NULL) 
  colnames(frame)=c("name","description")
  return(frame)
}

listFilters <- function( mart ){
  if(missing( mart ) || class( mart )!='Mart') stop("No Mart object given or object not of class 'Mart'")
  filterList=mget(ls(mart@filters), env=mart@filters)
  mat <- do.call(cbind, filterList)
  frame=data.frame(cbind(colnames(mat),mat[1,]), row.names=NULL) 
  colnames(frame)=c("name","description")
  return(frame)
}


##------------------------------------------------------------
## getBM
##------------------------------------------------------------

getBM <- function(attributes, filters, values, mart, curl = NULL, output = "data.frame", list.names = NULL, na.value = NA){

  if(missing( mart ) || class( mart )!='Mart')
    stop("Argument 'mart' must be specified and be of class 'Mart'.")
  if(missing( attributes ))
    stop("Argument 'attributes' must be specified.")
  if(missing( filters )) 
    stop("Argument 'filters' must be specified.")
  if(missing( values ))
    stop("Argument 'values' must be specified.")
  if(output != "data.frame" && output != "list")
    stop("Only data.frame and list are valid output formats for this function")

  invalid = !(attributes %in% ls(mart@attributes))
  if(any(invalid))
    stop(paste("Invalid attribute(s):", paste(attributes[invalid], collapse=", "),
               "\nPlease use the function 'listAttributes' to get valid attribute names"))
  invalid = !(filters %in% ls(mart@filters))
  if(any(invalid))
    stop(paste("Invalid filters(s):", paste(filters[invalid], collapse=", "),
               "\nPlease use the function 'listFilters' to get valid filter names"))

  ## use the mySQL interface
  if(mart@mysql){
    if(output == "data.frame"){
      if(length(filters) > 1) stop("biomaRt currently allows only one filter per query, reduce the number of filters to one")
      query <- queryGenerator(attributes=attributes, filter=filters, values=values, mart=mart)
      res <- dbGetQuery(mart@connections$biomart,query)
      if(dim(res)[1] == 0){
        res <- data.frame()
      }
      else{
        mt = match(res[,1], values);
        if(any(is.na(mt)))
          stop("Internal error!");
        res = res[order(mt),];
        names(res) = c(filters, attributes);     
      }  
      return(as.data.frame(res))
    }
    else{
      out <- vector("list", length(attributes))
      if(is.null(list.names))
        names(out) <- attributes
      else
        names(out) <- list.names
      
      for(j in seq(along = attributes)){
        tmp <- getBM(attributes[j], filters, values, mart)
        tmp2 <- vector("list", length(values))
        names(tmp2) <- values
        for(i in seq(along=tmp2)){
          tst <- tmp[tmp[,1] %in% values[i],2]
          tst <- tst[!is.na(tst)]
          if(length(tst) == 0) tst <- na.value
          tmp2[[i]] <- tst
        }
        out[[j]] <- tmp2
      }
     return(out)
    }
  }
  else{
    ## use the http/XML interface
    if(output=="data.frame"){
      xmlQuery = paste("<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0'> <Dataset name = '",mart@dataset,"'>",sep="")
      attributeXML =  paste("<Attribute name = '", attributes, "'/>", collapse="", sep="")
      if(length(filters) > 1){
        if(class(values)!= "list")
          stop("If using multiple filters, the 'value' has to be a list.\nFor example, a valid list for 'value' could be: list(affyid=c('1939_at','1000_at'), chromosome= '16')\nHere we select on affyid and chromosome, only results that pass both filters will be returned");
        filterXML = NULL
        for(i in 1:length(filters)){
         
          filtmp = strsplit(get(filters[i],env=mart@filters)$field, "_")
          if(filtmp[[1]][length(filtmp[[1]])] == 'bool'){
            filterXML = paste(filterXML,paste("<BooleanFilter name = '",filters[i],"' />", collapse="",sep=""),sep="")
          }
          else{
            valuesString = paste(values[[i]],"",collapse=",",sep="")
            filterXML = paste(filterXML,paste("<ValueFilter name = '",filters[i],"' value = '",valuesString,"' />", collapse="",sep=""),sep="")
          }
        }
      }
      else{
       
        filtmp = strsplit(get(filters,env=mart@filters)$field, "_")
        if(filtmp[[1]][length(filtmp[[1]])] == 'bool'){
         filterXML = paste("<BooleanFilter name = '",filters,"' />", collapse="",sep="")
        }
        else{
         valuesString = paste(values,"",collapse=",",sep="")
         filterXML = paste("<ValueFilter name = '",filters,"' value = '",valuesString,"' />", collapse="",sep="")
        }

      }
      xmlQuery = paste(xmlQuery, attributeXML, filterXML,"</Dataset></Query>",sep="")
      
      if(is.null(curl)){
        postRes = postForm(paste(mart@host,"?",sep=""),"query" = xmlQuery)
      }
      else{
        postRes = postForm(paste(mart@host,"?",sep=""),"query" = xmlQuery, curl = curl)
      }
      
      if(postRes != ""){
        ## convert the serialized table into a dataframe
        con = textConnection(postRes)
        result = read.table(con, sep="\t", header=FALSE, quote = "", comment.char = "", as.is=TRUE)
        close(con)
        if( substr(as.character(result[1]),1,5) == "ERROR"){  #Will forward the webservice error
          stop(paste("\n",result,"The following query was attempted, use this to report this error\n",xmlQuery ,sep="\n"))
        }
        ## check and postprocess
        if(all(is.na(result[,ncol(result)])))
          result = result[,-ncol(result),drop=FALSE]
        stopifnot(ncol(result)==length(attributes))
        if(class(result) == "data.frame"){
          colnames(result) = attributes
        }
      } else {
        #stop("The getBM query to BioMart webservice returned no result.  The webservice could be temporarily down, please try query again.")
        ##
        ## FIXME (wh 1.6.2006) - do we really just want to quietly fail, without warning or error? I don't like this,
        ## this will cause hard to trace problems in scripts and processing pipelines.
        result=NULL
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
          tst <- getBM(attributes = attributes[j], filters=filters, values = values[k], mart = mart)
       
          if(!is.null(tst)){
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

#########################################################
#cleanBM: Function to clear getBM output for duplicates #
#########################################################

cleanBM=function(bmresult, query_id_col=1, result_id_col=2 ){
 sel = bmresult[,result_id_col] != ""
 bmresult = bmresult[sel,]
 dupli=duplicated(bmresult[,c(query_id_col,result_id_col)])
 bmresult=bmresult[!dupli,]
 return(bmresult)
}


queryGenerator <- function(attributes, filter, values, mart){

  ## attributes
  Afield <- Atab <- Akey <- NULL
  for(i in 1:length(attributes)){
    if(!exists(attributes[i],mart@attributes)){
      stop(paste("attribute: ",attributes[i]," not found, please use the function 'listAttributes' to get valid attribute names",sep=""))
    }
    else{
      Afield = c(Afield,get(attributes[i],mart@attributes)$field)
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
    stop(sprintf("'length(filter)' must be 1 when using MySQL queries."))
  
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

  if(length(unique(c(Atable,Ftable)))>2)
    stop(paste("This query is currently not possible: your attributes extend over multiple tables (",
               paste(c(Atable,Ftable), collapse=", "), "), please only use attributes from one table per query or do your query via the webservice.", sep=""))
  
  query = paste("SELECT DISTINCT",
    paste(Ftable, Ffield, sep=".", collapse=", "), ",",
    paste(Atable, Afield, sep=".", collapse=", "),
    "FROM")

  ########################################
  # Ensembl Key Order: gene overrides transcript #
  ########################################
  
  matchGeneKey=match("gene_id_key",c(Akey,Fkey))
  matchTransKey=match("transcript_id_key",c(Akey,Fkey))

  if(!is.na(matchGeneKey) && !is.na(matchTransKey)){ # so keys are mixed and we have to take order into account
     Akey[Akey =="transcript_id_key"] = "gene_id_key"
     Fkey[Fkey =="transcript_id_key"] = "gene_id_key"
  }
  
  ### end Ensembl specific part #######
  
  tables <- unique(c(Atable, Ftable));
 
  ## FIXME: why only for tables[1:2] and Akey[1] and not for all?

  query = paste(query, tables[1])

  if(length(tables) > 1){

     query = paste(query, " INNER JOIN ",tables[2], " ON ", tables[1], ".", Akey[1], " = ",tables[2], ".", Fkey[1], sep="")  #old query
    
   # query = paste(query, " INNER JOIN (",paste("",tables[-1],sep="",collapse=","), ") ON (",paste(paste(tables[-1],Akey[1],sep="."), paste(tables[1], Akey[1], sep="."),sep=" = ", collapse=" AND "))
     #query = paste(query, ")", sep="")
  }
  
  query = paste(query, " WHERE ", paste(Ftable, Ffield[1], sep =".")," IN (",
    paste("'", values, "'", collapse=", ", sep=""), ")", sep="")

  #The last part of the query is to avoid returning duplicate results that are the same but where one of the returned values is NULL due to multiple entries in the database
  #By commenting this you can see it's effect
  
  #for(j in 1:length(Afield)){
  # query = paste(query," AND ",paste(Atable[j], Afield[j], sep="."),"!='NULL'")
  #}
  
  return(query)
}
  
##---------TEST FUNCTION FOR WEBSERVICE QUERIES -----

testService <- function(mart){

#query = "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0' ><Dataset name = 'hsapiens_genomic_sequence'><ValueFilter name = 'chr' value = '1'/><ValueFilter name = 'start' value = '10000'/><ValueFilter name = 'end' value = '10020'/><Attribute name = 'raw_sequence'/></Dataset></Query>"

#query = "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0'> <Dataset name = 'hsapiens_snp'><Attribute name = 'tscid'/><Attribute name = 'refsnp_id'/><Attribute name = 'allele'/><Attribute name = 'chrom_start'/><Attribute name = 'chrom_strand'/><ValueFilter name = 'chr_name' value = 'Y' /><ValueFilter name = 'snp_chrom_start' value = '2000' /><ValueFilter name = 'snp_chrom_end' value = '40000' /></Dataset></Query>"

#query = "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0'><Dataset name = 'mmusculus_gene_ensembl'><Attribute name = 'gene_stable_id'/></Dataset><Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'entrezgene'/><ValueFilter name = 'entrezgene' value = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20' /></Dataset><Links source = 'hsapiens_gene_ensembl' target = 'mmusculus_gene_ensembl' defaultLink = 'mmusculus_internal_gene_id'/></Query>"

#query = "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0'><Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'gene_stable_id'/></Dataset><Dataset name = 'mmusculus_gene_ensembl'><ValueFilter name = 'affy_mg_u74av2' value = '95919_at'/></Dataset><Links source = 'mmusculus_gene_ensembl' target = 'hsapiens_gene_ensembl' defaultLink = 'hsapiens_internal_gene_id'/></Query>"

#query = "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' count = '0' ><Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'gene_stable_id' /><Attribute name = 'affy_hg_u95c' /></Dataset><Dataset name = 'mmusculus_gene_ensembl'><ValueFilter name = 'affy_mg_u74av2' value = '95919_at'/></Dataset><Links source = 'mmusculus_gene_ensembl' target = 'hsapiens_gene_ensembl' defaultLink = 'hsapiens_internal_gene_id' /></Query>"
  
result = postForm(paste(mart@host,"?",sep=""),"query"=query)
result = strsplit(result,"\n")[[1]]
result = sapply(result,"strsplit","\t")
return(result)
}
