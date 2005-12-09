.packageName <- "biomaRt"

setClass("Mart",
         representation(
                        connections = "list",
                        mysqldriver = "MySQLDriver",
                        arrayToSpecies = "data.frame",
                        datasets = "list"
                        ),
         prototype(
                   connections = new("list"),
                   arrayToSpecies = data.frame(cbind(x=1, y=1:2))
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


listMarts <- function( mart, host, user, password, includeHosts = FALSE){

  if(missing(mart)){
    mart <- c("ensembl","vega","snp","msd","uniprot")
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
            v <- as.numeric(strsplit(res[matches[j],1],"_")[[1]][3]);
            if(v > version){
              latest <- j;
              version <- v;
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

##################
#Connect to marts#
##################

martConnect <- function(biomarts = "ensembl", host, user, password, mart, local = FALSE){
  
    if(missing(mart)){
      driver <- dbDriver("MySQL", force.reload = FALSE);
      mart <- new("Mart", mysqldriver = driver)
    }
    presentConnections <- length(mart@connections)
    for(i in 1: length(biomarts)){
      if(local){
        if(biomarts[i] == "ensembl" || biomarts[i] == "sequence" || biomarts[i] == "uniprot" || biomarts[i] == "snp" || biomarts[i] == "vega"){
          if(!missing(host) && !missing(user) && !missing(password)){
            database <- listMarts(mart = biomarts[i], host = host[i], user = user[i], password = password[i]);
            mart@connections[[i + presentConnections]] <- dbConnect(drv = mart@mysqldriver,user = user[i], host = host[i] , dbname = database, password = password[i])
            writeLines(paste("connected to: ",database))
            names(mart@connections)[i + presentConnections] <- biomarts[i]
          }
          else{
            stop("you should provide host, user and password when using local databases")
          }
        }
        else{
          stop(paste("mart ",biomarts[i]," is not covered by biomaRt, please choose one of the following BioMarts to connect to: ensembl, vega, sequence, uniprot and snp"))
        }
      }
      else{
        if(biomarts[i] == "ensembl" || biomarts[i] == "sequence" || biomarts[i] == "uniprot" || biomarts[i] == "snp" || biomarts[i] == "vega"){
          database <- switch(biomarts[i],
                             ensembl = listMarts(mart = biomarts[i], host = "ensembldb.ensembl.org", user = "anonymous", password = ""),
                             sequence = listMarts(mart = biomarts[i], host = "ensembldb.ensembl.org", user = "anonymous", password = ""),
                             uniprot = listMarts(mart = biomarts[i], host = "martdb.ebi.ac.uk", user = "anonymous", password = ""),
                             snp = listMarts(mart = biomarts[i], host = "ensembldb.ensembl.org", user = "anonymous", password = ""),
                             vega = listMarts(mart = biomarts[i], host = "ensembldb.ensembl.org", user = "anonymous", password = ""));
          
          #should have a try and catch here ...
          
          mart@connections[[i + presentConnections]] <- switch(biomarts[i],
                                                               ensembl = dbConnect(drv = mart@mysqldriver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = database, password = ""),
                                                               sequence = dbConnect(drv = mart@mysqldriver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = database, password = ""),
                                                               snp = dbConnect(drv = mart@mysqldriver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = database, password = ""),
                                                               vega = dbConnect(drv = mart@mysqldriver,user = "anonymous", host = "ensembldb.ensembl.org" , dbname = database, password = ""),
                                                               uniprot = dbConnect(drv = mart@mysqldriver,user = "anonymous", host = "martdb.ebi.ac.uk" , dbname = database, password = ""));
          writeLines(paste("connected to: ",database))
          names(mart@connections)[i + presentConnections] <- biomarts[i]
        }
        else{
          stop(paste("mart ",biomarts[i]," is not covered by biomaRt, please choose one of the following BioMarts to connect to: ensembl, vega, sequence, uniprot and snp"))
        }
      }
    }
    
    if(match("ensembl",biomarts,nomatch=0) > 0){
      xref <- getPossibleXrefs( mart );
      affySelect <- grep("affy", xref[,2]);
      affyTables <- xref[affySelect,];
      affyTables[,2] <- gsub("affy_","", affyTables[,2]);
      affyPackName <- gsub("_","",affyTables[,2]);
      mart@arrayToSpecies <- data.frame(affyID = affyPackName,EnsemblArrayID = affyTables[,2],species = affyTables[,1]);
    }
    return( mart )
}


######################
#Disconnect from mart#
######################

martDisconnect <- function( mart ){

  if(missing(mart)){
    stop("no Mart object to disconnect");
  }
  openConnections <- names(mart@connections);
  
  if(match("ensembl",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$ensembl );
  }
  if(match("vega",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$vega );
  }
  if(match("sequence",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$sequence );
  }
  if(match("snp",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$snp );
  }
  if(match("uniprot",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$uniprot );
  }
  if(match("biomart",openConnections,nomatch = 0) != 0){
    dbDisconnect( mart@connections$biomart );
  }
}

####### local functions #############


mapArrayToEnsemblTable <- function(array = NULL, species = NULL, mart = NULL, dbtable = NULL){
  
  table <- NULL;
  array <- gsub("_","",array)
  array <- mart@arrayToSpecies$EnsemblArrayID[ match( array, mart@arrayToSpecies$affyID)]; 
  if(dbtable == "xrefdm"){
    table <- paste( species,"_gene_ensembl__xref_affy_",array,"__dm",sep="");
  }
  if(dbtable == "affybool"){
    table <- paste("affy_",array,"_bool",sep="");
  }
  return( table );
}

mapSpeciesToGeneTable <- function( species = NULL, db = NULL){

  table <- NULL;
  if(db == "ensembl"){
    table <- paste( species,"_gene_ensembl__gene__main",sep = "");
  }
  
  if(db == "vega"){
    table <- paste( species,"_gene_vega__gene__main",sep = "");
  }
  return( table );
}

mapSpeciesToGOTable <- function( species = NULL){

  table <- paste(species,"_gene_ensembl__xref_go__dm",sep = "");
  return(table); 
}

mapSpeciesToUNIPROTTable <- function( species = NULL){

#  table <- paste(species,"_gene_ensembl__xref_uniprot_accession__dm",sep = "");
  table <- paste(species,"_gene_ensembl__prot_interpro__dm",sep="");
  return(table); 
}


mapArrayToSpecies <- function( array = NULL, mart = NULL ){

 species = "";
 array <- gsub("_","",array) 
 species <- mart@arrayToSpecies$species[ match( array, mart@arrayToSpecies$affyID)]; 
 return( species )
}

mapSpeciesToEntrezGene <- function( species = NULL, db = NULL ){

  table <- NULL;
  if(db == "ensembl"){
    table <- paste( species, "_gene_ensembl__xref_entrezgene__dm", sep = "");
   }
  if(db == "vega"){
     table <- paste( species, "_gene_vega__xref_entrezgene__dm", sep = "");
  }
  
  return( table );
}

mapSpeciesToHUGO <- function( species = NULL, db = NULL ){

  table <- NULL;
  if(db == "ensembl"){
    table <- paste( species, "_gene_ensembl__xref_hugo__dm", sep = "");
   }
  if(db == "vega"){
    stop("No HUGO annotation currently present in VEGA")
  }
  return( table );
}


mapSpeciesToRefSeq <- function( species = NULL, db = NULL ){

  table <- NULL;
  
  if(db == "ensembl"){
    table <- paste( species,"_gene_ensembl__xref_refseq_dna__dm" , sep = "");
  }
  if(db == "vega"){
    table <- paste( species,"_gene_vega__xref_refseq__dna__dm" , sep = "");
  }
  return( table );
}

mapSpeciesToEMBL <- function( species = NULL ){

  table <- paste( species,"_gene_ensembl__xref_embl__dm" , sep = ""); 
  return( table );
}


mapSpeciesToHomologTable <- function(fromESpecies = NULL, toESpecies = NULL) {
  
  table <- NULL;
  table <- paste(fromESpecies,"_gene_ensembl__homologs_",toESpecies,"__dm",sep="");
  return(table);
}

mapSpeciesToFlybase <- function ( species = NULL){
  table <- NULL;
  if(species == "dmelanogaster"){
    table <- "dmelanogaster_gene_ensembl__xref_flybase_gene_id__dm";
  }
  else{
    stop("Flybase id's can only be used with dmelanogaster as species")
  }
  return(table);
}

getTableColumn <- function(type = NULL){
  tableCol <- switch(type,
                     entrezgene = "dbprimary_id",
                     refseq     = "dbprimary_id",
                     embl       = "dbprimary_id",
                     hugo       = "display_id_list",
                     affy       = "dbprimary_id",
                     flybase    = "dbprimary_id",
                     ensembl    = "gene_stable_id",
                     ensemblTrans = "transcript_stable_id"
                     );
  return(tableCol);
}
######### public functions ##############

######################
#get gene information#
######################


getGene <- function( id, type, array, species, db = "ensembl", mart, output = "data.frame"){

  IDTable <- NULL
  
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if( db != "ensembl" && db != "vega"){
    stop("you can only use ensembl or vega");
  }
     
  if(!missing(array)){
    type <- "affy";
     if( db != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }
  }
  if( missing( type )){
    stop("you must provide the identifier type using the type argument")
  }

  if( type == "affy" ){
    if( missing( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  else{
    if(missing( species ) ){
      stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
    }
  }
  
  speciesTable <- mapSpeciesToGeneTable( species, db = db);
  
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo")  && !(type == "ensembl") && !(type == "ensemblTrans") && !(type == "flybase")){
    stop("invalid type choose either affy, refseq, embl, hugo, ensembl, ensemblTrans or entrezgene");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = species, mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = db ),
                    ensemblTrans =   paste(species,"_gene_ensembl__transcript__main",sep=""),
                    entrezgene = mapSpeciesToEntrezGene( species, db = db),
                    hugo = mapSpeciesToHUGO( species, db = db),
                    refseq = mapSpeciesToRefSeq( species, db = db),
                    embl = mapSpeciesToEMBL( species ),
                    flybase = mapSpeciesToFlybase( species )
                    );
                    

  dbcolID <- getTableColumn("ensembl");
  dbcolQID <- getTableColumn(type);
  
  if (db == "ensembl" || db == "vega"){
    
    if (length( id ) >= 1){
      
      ids <- paste("'",id,"'",sep="",collapse=",");
      
      if(type == "ensembl"){
         query <- paste("select distinct ",dbcolID,", display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end, chrom_strand  from ", speciesTable," where ",dbcolQID," in (",ids,")",sep="");
      }
      else{
          query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".gene_stable_id,",speciesTable,".display_id,",speciesTable,".description, ",speciesTable,".band,",speciesTable,".chr_name,",speciesTable,".gene_chrom_start,",speciesTable,".gene_chrom_end,",speciesTable,".chrom_strand  from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,")",sep="");        
      }

      if(db == "vega"){
        if(match("vega",names(mart@connections), nomatch=0) == 0){
          stop("You need a connection to vega for this query, use martConnect and include 'vega' in your biomarts vector");
        }
        res <- dbGetQuery( conn = mart@connections$vega, statement = query);
      }
      if(db == "ensembl"){
        res <- dbGetQuery( conn = mart@connections$ensembl, statement = query);
      }
      if(dim(res)[1] == 0){
        table <- new("martTable", id = id, table = list(symbol = NA, description = NA, band = NA, chromosome = NA, start = NA, end = NA, martID = NA));
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
        if(output == "martTable"){
         table <- new("martTable", id = as.vector(res[,1]), table = as.list(res[,-1]));
        }
        else{
          table <- as.data.frame(res)
        }
      }
    }
  }  
  return(table);  
}

########################
#get feature           #
########################


getFeature <- function( symbol, OMIM, OMIMID, GO, GOID, array, species, chromosome, start, end, type,  mart){
  
  table <- NULL
  
  if(!missing(array)){
    type <- "affy"
  }
  
  else{
    if( missing( type )){
      stop("you must provide the identifier type using the type argument");
    }
  }
  
  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect");
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector");
  }
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl")  && !(type == "hugo") && !(type == "ensembl")  && !(type == "flybase") ){
    stop("invalid type choose either affy, refseq or entrezgene");
  }
  
  if( type == "affy" ){
    if( !missing( array ) ){
      species <- mapArrayToSpecies( array = array, mart = mart );
      affyBool <- mapArrayToEnsemblTable( array = array, species = species, mart = mart, dbtable = "affybool"  );    
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  else{
    if(missing(species)){
      stop("you must provide the species, valid species names can be found by using the function getSpecies()")
    }
  }
  speciesTable <-  mapSpeciesToGeneTable( species, db = "ensembl" );
  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = species, mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ),
                    flybase = mapSpeciesToFlybase( species )
                    );
    
  dbcolID <- getTableColumn(type);        #get database id col
  
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
    GOTable <- mapSpeciesToGOTable( species );
    query <- paste("select distinct ",IDTable,".",dbcolID,",",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_id_key = ",GOTable,".gene_id_key where ",GOTable,".description like '%",GO,"%' and ",IDTable,".",dbcolID," != 'NULL'",sep="");
  }
  if(!missing(GOID)){
    GOTable <- mapSpeciesToGOTable( species );
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
  
  res <- dbGetQuery( conn = mart@connections$ensembl,statement = query);

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

###################
#get GO annotation#
###################


getGO <- function( id, type, array, species, mart, output="data.frame"){

  table <- NULL;
  go <- NULL;

  if( missing( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
 
  if(!missing(array)){
    type <- "affy";
  }
  else{
    if(missing( type )){
      stop("you must provide the identifier type using the type argument")
    }
  }
  
  if( type == "affy" ){
    if( is.null( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  else{
    if(missing( species ) ){
      stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
    }
  }

  GOTable <- mapSpeciesToGOTable( species );

  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo") && !(type == "ensembl") && !(type == "flybase")){
    stop("invalid type choose either affy,  entrezgene, hugo, ensembl, refseq or embl");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = species, mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ),
                    flybase = mapSpeciesToFlybase( species )
                    );
    
  dbcolID <- getTableColumn("ensembl");        #get database id col   
  dbcolQID <- getTableColumn(type);        #get database id col
 
  if( length( id ) >= 1){
    ids <- paste("'",id,"'",sep="",collapse=",")
    res <- NULL
    if(type == "ensembl"){
      query <- paste("select distinct ",IDTable,".",dbcolQID,",",IDTable,".",dbcolQID,", ",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
    }
    else{
      query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolID,",",GOTable,".dbprimary_id, description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".",dbcolQID," in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
  }
    res <- dbGetQuery(conn = mart@connections$ensembl, statement = query);
    
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
      if(output == "martTable"){
        table <- new("martTable", id = as.vector(res[,1]), table = as.list(res[,-1]))
      } 
      else{
        table <- as.data.frame(res)
      }
    }
  }
  return( table );
}

#####################
#get OMIM annotation#
#####################


getOMIM <- function( id, type, array, mart, output="data.frame"){
  
  OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
  omim <- NULL;
  table <- NULL;
  species <- "hsapiens";

  if( missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if(!missing(array)){
    type <- "affy";
  }
  else{
    if( missing( type )){
      stop("you must provide the identifier type using the type argument")
    }
  }

  if( type == "affy" ){
    if( missing( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  else{
    if( missing( species ) ){
      stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
    }
  }
  
  if(!(type == "affy") && !(type == "refseq") && !(type == "embl") && !(type == "entrezgene") && !(type == "hugo") && !(type == "ensembl")){
    stop("invalid type choose either affy, refseq. embl, hugo, ensembl or entrezgene");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = "hsapiens", mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ));

  dbcolID <- getTableColumn("ensembl");        #get database id col   
  dbcolQID <- getTableColumn(type);        #get database id col
 
  if( length( id ) >= 1){
    ids <- paste("'",id,"'",sep="",collapse=",")
    res <- NULL
    if(type == "ensembl"){
      query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolQID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".",dbcolQID," in (",ids,")  and ",OMIMTable,".omim_id != 'NULL'",sep="");
      
    }
    else{
      query <- paste("select distinct ",IDTable,".",dbcolQID,", ",IDTable,".",dbcolID,",",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".",dbcolQID," in (",ids,") and ",OMIMTable,".omim_id != 'NULL'",sep="");
  }
    res <- dbGetQuery(conn = mart@connections$ensembl, statement = query);
    
    
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
      if(output == "martTable"){ 
        table <- new("martTable", id = as.vector(res[,1]), table = as.list(res[,-1]))
      }
      else{
        table <- as.data.frame(res)
      }
    }  
  }
  return( table );  
}

#####################
#getSequence
#####################


getSequence <- function(species, chromosome, start, end, martTable, mart){
  if(missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("sequence",names(mart@connections),nomatch=0) == 0){
    stop("You are missing a database connection to sequence BioMart for this query.  Add connection to a sequence BioMart to your Mart object via the function martConnect. Use the following command to do this:  martConnect(biomarts='sequence',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password=''). ")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }

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
        chunkseq <- dbGetQuery(conn = mart@connections$sequence, statement = query);
        newstart <- start[i] - (floor((start[i]-1)/100000) * 100000)
        newend <- end[i] - (floor((end[i]-1)/100000) * 100000)
        
        sequence <- c(sequence,substr(as.character(chunkseq), newstart, newend));
        
      }
      
      else{   #query sequence is on 2 sequence chuncks
        
        query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
        chunkseq1 <- dbGetQuery(conn = mart@connections$sequence, statement = query);
        query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkEnd,"'",sep="");
        chunkseq2 <- dbGetQuery(conn = mart@connections$sequence, statement = query);
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
          chunkseq <- dbGetQuery(conn = mart@connections$sequence, statement = query);
          newstart <- start[i] - (floor((start[i]-1)/100000) * 100000)
          newend <- end[i] - (floor((end[i]-1)/100000) * 100000)
          
          sequence <- c(sequence,substr(as.character(chunkseq), newstart, newend));
          
        }
        
        else{   #query sequence is on 2 sequence chuncks
          
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
          chunkseq1 <- dbGetQuery(conn = mart@connections$sequence, statement = query);
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkEnd,"'",sep="");
          chunkseq2 <- dbGetQuery(conn = mart@connections$sequence, statement = query);
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

################
#show affy info#
################


getAffyArrays <- function(mart){
  if(missing(mart)){
    stop("no Mart object found, create this first using the function martConnect and then give as argument in the function")
  }
  print( mart@arrayToSpecies );
}

#######################
#getSNP
#######################


getSNP <- function(species, chromosome, start, end, mart){
  if( missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if(match("snp",names(mart@connections),nomatch=0) == 0){
    stop("You are missing a database connection to snp for this query.  Add connection to snp a BioMart to your Mart object via the function martConnect.   Use the following command to do this:  martConnect(biomarts='snp',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password='')")
  }
  if(missing(chromosome) || missing(start) || missing(end) || missing(species)){
    stop("you have to give chromosome, start and end positions and species as arguments, see ?getSNP for more information")
  }
  
  table <- NULL
  ensemblTable <- paste(species,"_snp__snp__main", sep ="");
  query <- paste("select snp_id_key,snp_chrom_start, allele, tscid, ensemblcoding_bool, ensemblintronic_bool, ensembl5utr_bool, ensembl3utr_bool, ensemblsyn_bool from ",ensemblTable," where chr_name = '",chromosome,"' and snp_chrom_start >= ",start," and snp_chrom_start <= ",end,sep="")
  res <- dbGetQuery( mart@connections$snp, query );
  if(dim(res)[1] == 0){
    stop("No SNP's found in selected region, check if the input you entered is correct")
  }
  else{
    table <- new("martTable", id = res$tscid ,table = list(snpStart = res[,2], allele = res$allele, coding = res[,5], intronic =  res[,6], syn =  res[,9],utr5 = res[,7], utr3= res[,8]))
  }
  return( table ); 
}

##################
#getHomolog      #
##################

getHomolog <- function(id, from.type, to.type, from.array, to.array, from.species, to.species, mart, output = "martTable") {
    
  db = "ensembl";
  if( missing( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }

  if( !missing( id )){
    id <- as.character(id);
  }
  else{
    stop("No id's to search for homologs given");
  }
  
  if( !missing( from.array )){
    from.type <- "affy";
    from.species <- mapArrayToSpecies( array = from.array, mart = mart );
  }
  else{
    if( missing( from.species ) ) {
      stop("You must provide the species to map FROM using the  from.species argument")
    }
    if ( missing( from.type )) {
      stop("You must provide the identifier type using the from.type  argument")
    }
  }

  if( !missing( to.array )){
    to.type <- "affy";
    to.species <- mapArrayToSpecies( array = to.array, mart = mart );
  }
  else{
    if ( missing( to.species )) {
      stop("You must provide the species to map TO using the to.species  argument");
    }
    if ( missing( to.type )) {
      stop("You must provide the identifier type using the to.type  argument")
    } 
  }
  
  
  fromIDTable <-
    switch(from.type,
           entrezgene = mapSpeciesToEntrezGene(from.species,db=db),
           refseq     = mapSpeciesToRefSeq(from.species,db=db),
           embl       = mapSpeciesToEMBL(from.species),
           hugo       = mapSpeciesToHUGO(from.species,db=db),
           affy       =  mapArrayToEnsemblTable( from.array, species = from.species, mart = mart, dbtable = "xrefdm"),
           ensembl    = mapSpeciesToGeneTable(from.species,db=db),
           flybase    = mapSpeciesToFlybase( from.species )
           );
  
  toIDTable <-
    switch(to.type,
           entrezgene = mapSpeciesToEntrezGene(to.species,db=db),
           refseq     = mapSpeciesToRefSeq(to.species,db=db),
           embl       = mapSpeciesToEMBL(to.species),
           hugo       = mapSpeciesToHUGO(to.species,db=db),
           affy       =  mapArrayToEnsemblTable( to.array, species = to.species, mart = mart, dbtable = "xrefdm"),
           ensembl    = mapSpeciesToGeneTable(to.species,db=db),
           flybase    = mapSpeciesToFlybase( to.species )
           );

  fromCol <- getTableColumn(from.type);
  toCol <- getTableColumn(to.type);

  if(to.type == "ensembl" && from.type != "ensembl"){
    homolTable <- mapSpeciesToHomologTable(to.species,from.species);
  }
  else{
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
    
    res <- dbGetQuery(conn = mart@connections$ensembl, statement = query);
    
    if (dim(res)[1] == 0) {
      table <- new("martTable", id = id, table = list( MappedID = rep(NA,length(id))))
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
      table <- new("martTable", id = foundID, table = list(MappedID = MappedID))
    }
  }
  return(table)
} 

###########################################################################################
#                          
#Ensembl specific functions
##########################################################################################


###########################
#Xref functions           #
###########################


getPossibleXrefs <-  function( mart ) {
  if ( missing( mart ) || class( mart ) != 'Mart') {
    stop('You must supply a valid mart object.  You can use martConnect to produce one')
  }
  if( match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
 
  res <- dbGetQuery( mart@connections$ensembl, "show tables like '%_gene_ensembl__xref%'")
  xref <- strsplit(res[,1],'_gene_ensembl__xref_')
  xref <- lapply(xref,function(s) {
    s[2] <- gsub('__dm','',s[2])
    return(s)
  })
  xrefdf <- do.call('rbind',xref)
  colnames(xrefdf) <- c('species','xref')
  return(xrefdf)
} 


getSpecies <- function(mart, db = c("ensembl")) {
  if ( missing( mart ) || class( mart )!='Mart') {
    stop('You must supply a valid mart object.  You can use martConnect to produce one')
  }
  if( match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
 
  #db <- match.arg(db)
  query = "show tables like '%gene__main'";

  if(match(db, names(mart@connections),nomatch=0) == 0){
    stop(paste("You are missing a database connection to ", db," for this query.  Add connection to ",db," a BioMart to your Mart object via the function martConnect.   Use the following command to do this:  martConnect(biomarts='",db,"',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password='')"),sep="")
  }
 
  res <- switch(db,
                ensembl  = dbGetQuery(mart@connections$ensembl,query),
                vega     = dbGetQuery(mart@connections$vega,query),
                snp  = dbGetQuery(mart@connections$snp,query),
                sequence = dbGetQuery(mart@connections$sequence,query)
                );
  
  resvec <- as.vector(as.character(res[,1]));
  speciessub <- grep('gene__main',resvec);
  specieslist <- strsplit(resvec[speciessub],'_');
  species <- sapply(specieslist,function(s) {s[1]},simplify=T);
  return( unique( species ) );
} 

getXref <- function( id, from.species, to.species, from.xref, to.xref, db = "ensembl", mart) {
  
  if ( missing( mart ) || class( mart )!='Mart'){
    stop('You must specify a valid Mart object which can be created using martConnect().');
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  } 
  if ( missing( id )) {
    stop('You need to give ID(s) to map from');
  }
  if ( missing( from.species )) {
    writeLines('fromspecies is a necessary argument\nValid species for this mart are:');
    print( getSpecies(mart = mart, db = "ensembl" ));
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
    if (db == 'vega'){
      stop('VEGA supports only hsapiens')
    }
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
  
  res <- dbGetQuery(con = mart@connections$ensembl, query);
  
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


getINTERPRO <- function( id, type, array, species, mart, output="data.frame"){
  
  genes <- NULL
  speciesTable <- NULL
  db <- "ensembl";
  results <- FALSE;
    
  if( missing( mart )|| class(mart)!='Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
 # if(match("uniprot",names(mart@connections), nomatch=0) == 0){
 #   stop("You need a connection to uniprot for this query, use martConnect and include 'uniprot' in your biomarts vector")
 # }
  if( !missing( array )){
    type <- "affy";
  }
  else{
    if( missing( type )){
      stop("you must provide the identifier type using the type argument");
    }
  }

  if( type == "affy" ){
    if( missing( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  else{
    if( missing( species ) ){
      stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
    }
  }

  uniprotTable <- mapSpeciesToUNIPROTTable( species );
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo")  && !(type == "ensembl")  && !(type == "ensemblTrans")  && !(type == "flybase")){
    stop("invalid type choose either affy, refseq, embl, hugo, ensembl or entrezgene");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = species, mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    ensemblTrans =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ),
                    flybase = mapSpeciesToFlybase( species )
                    );

  dbcolID <- getTableColumn(type);        #get database id col
 
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

    res <- dbGetQuery( conn = mart@connections$ensembl,statement = query);
    
    if(dim(res)[1] == 0){
      table <- new("martTable", id = id, table = list(INTERPROID = NA, shortdescription = NA, description = NA))
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
      if(output=="martTable"){
        table <- new("martTable", id = as.vector(res[,1]), table = as.list(res[,-1]))
      }
      else{
        table <- as.data.frame(res)
      }
    }
  }
  return(table)
}

######################################################################################
#
#BioMart functions
#based on martshell
#
######################################################################################

useMart <- function(biomart, host, user, password, local = FALSE){

  
  driver <- dbDriver("MySQL", force.reload = FALSE);
  mart <- new("Mart", mysqldriver = driver)
  

  if(local){
    if(!missing(host) && !missing(user) && !missing(password)){
      database <- listMarts(mart = biomart, host = host, user = user, password = password);
      mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver,user = user, host = host, dbname = database, password = password)
      writeLines(paste("connected to: ",database[1,1]))
    }
    else{
      stop("you should provide host, user and password when using local databases")
    }
  } 
  else{
    database <- listMarts(includeHosts=TRUE)
    m<-match(biomart,database[,1],nomatch=0)
    if(m == 0){
      stop(paste("BioMart ",biomart, " does not exist or can not be found", sep =""))
    }
    else{
      mart@connections[["biomart"]] <- dbConnect(drv = mart@mysqldriver,user = "anonymous", host = database[m,2] , dbname = database[m,1], password = "")
      writeLines(paste("connected to: ",biomart))
    }
  }  
  return( mart )
}




listDatasets <- function( mart ){
  res <- dbGetQuery(mart@connections$biomart,"select dataset, version from meta_configuration where visible = 1")
  return(res)
}

getAttributes<-function( xml ){
  names <- names(xml)
  writeLines("Checking attributes ...", sep =" ")
  attrib <- table <- fieldA <- keyA <- NULL
  
  for(h in 1:(xmlSize(xml)-1)){
    if(names[h] == "AttributePage"){  #fix xmlSize ...first test if size > 0!!
      if(!is.null(xmlGetAttr(xml[[h]],"hidden"))){
        if(xmlGetAttr(xml[[h]],"hidden") != "true"){                
          for(i in 1:xmlSize(xml[[h]])){
            for(j in 1:xmlSize(xml[[h]][[i]])){
              for(k in 1:xmlSize(xml[[h]][[i]][[j]])){
                if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"hidden"))){
                  if(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"hidden") != "true"){                
                    if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))){
                      table<-c(table,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))
                      attrib<-c(attrib,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"internalName"))
                      fieldA <- c(fieldA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"field"))
                      keyA <- c(keyA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"key"))
                    }
                  }
                }
                else{
                  if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))){
                    table<-c(table,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))
                    attrib<-c(attrib,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"internalName"))
                    fieldA <- c(fieldA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"field"))
                    keyA <- c(keyA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"key"))
                    
                  }
                }
              }
            }
          }
        }
      }
      else{
        for(i in 1:xmlSize(xml[[h]])){
          for(j in 1:xmlSize(xml[[h]][[i]])){
            for(k in 1:xmlSize(xml[[h]][[i]][[j]])){
              if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"hidden"))){
                if(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"hidden") != "true"){                
                  if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))){
                    table<-c(table,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))
                    attrib<-c(attrib,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"internalName"))
                    fieldA <- c(fieldA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"field"))
                    keyA <- c(keyA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"key"))
                   
                  }
                }
              }
              else{
                if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))){
                  table<-c(table,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))
                  attrib<-c(attrib,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"internalName"))
                  fieldA <- c(fieldA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"field"))
                  keyA <- c(keyA,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"key"))
                   
                }
              }
            }
          }
        }
      }
    }
  }
  attributes <- list(attributes = attrib, field = fieldA, table = table, key = keyA)
  writeLines("ok")
  return(attributes)
}


getFilters <- function(xml){
  writeLines("Checking filters ...", sep=" ")

  names <- names(xml)
  
  filter <- NULL
  tableF <- NULL
  fieldF <- NULL
  keyF <- NULL
  ## FIXME:
  ## 1. perhaps could this be simplified (and made faster, b/c currently it's quite slow) by using
  ##    "xmlSApply". One could perhaps avoid the awkward  v = c(v, x)  constructs  to append vectors.
  ## 2. Should make sure that we don't fall into the "for(i in 1:length(x))" trap - consider
  ##    what happens when length(x)=0. Usually "for(i in seq(along=x))" is the way that matches
  ##    better what is intended.

  for(h in 1:(xmlSize(xml)-1)){
    if(names[h] == "FilterPage"){
      if(!is.null(xmlGetAttr(xml[[h]],"displayName"))){
        if(xmlAttrs(xml[[h]])[["displayName"]]=="FILTERS"){
          for(i in 1:xmlSize(xml[[h]])){
            for(j in 1:xmlSize(xml[[h]][[i]])){
              if(xmlSize(xml[[h]][[i]][[j]]) > 0){
                for(k in 1:xmlSize(xml[[h]][[i]][[j]])){
                  if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))){
                    tableF<-c(tableF,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"tableConstraint"))
                    filter<-c(filter,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"internalName"))
                    fieldF <- c(fieldF,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"field"))
                    keyF <- c(keyF,xmlGetAttr(xml[[h]][[i]][[j]][[k]],"key"))
                    
                  }
                  else{
                    if(xmlSize(xml[[h]][[i]][[j]][[k]]) > 0){
                      for(s in 1:xmlSize(xml[[h]][[i]][[j]][[k]])){
                        if(!is.null(xmlGetAttr(xml[[h]][[i]][[j]][[k]][[s]],"tableConstraint"))){
                          tableF<-c(tableF,xmlGetAttr(xml[[h]][[i]][[j]][[k]][[s]],"tableConstraint"))
                          filter<-c(filter,xmlGetAttr(xml[[h]][[i]][[j]][[k]][[s]],"internalName"))
                          fieldF<-c(fieldF,xmlGetAttr(xml[[h]][[i]][[j]][[k]][[s]],"field"))
                          keyF <- c(keyF,xmlGetAttr(xml[[h]][[i]][[j]][[k]][[s]],"key"))
                        }
                      }
                    }
                  }  
                }
              }
            }
          }
        }
      }
    }
  }

  filters <- list(filter = filter,field = fieldF, table = tableF, key = keyF)
  writeLines("ok")
  return(filters)
}

useDataset <- function(dataset, mart){
  res <- dbGetQuery(mart@connections$biomart,paste("select xml from meta_configuration where dataset = '",dataset,"'",sep=""))
  writeLines(paste("Reading database configuration of:",dataset))
  if(dim(res)[1] == 0){
    writeLines("This dataset is not accessible from biomaRt as not xml description of dataset is available")
  }
  else{
    xml <- xmlTreeParse(res[1,])
    xml <- xml$doc$children[[1]]
  }
  #Getting the attributes & filters
 
  attributes <- getAttributes(xml)
  filters <- getFilters(xml)
  mainTables <- getMainTables(xml)
  datasets <- list(name = dataset, xml= xml, mainTables = mainTables, attributes = attributes, filters = filters);
  mart@datasets <- datasets
  
  return(mart)
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
  return(mart@datasets$attributes$attributes)
}

listFilters <- function( mart ){
  return(mart@datasets$filters$filter)

}
    
getBM <- function(attributes, filter, values, mart){
  query <- queryGenerator(attributes=attributes, filter=filter, values=values, mart=mart)
  res <- dbGetQuery(mart@connections$biomart,query)
  if(dim(res)[1] == 0){
    res <- data.frame()
  }
  else{
    mt = match(res[,1], values);
    if(any(is.na(mt)))
      stop("Internal error!");
    res = res[order(mt),];
    names(res) = c(filter, attributes);     
  }  
  return(as.data.frame(res))
}


queryGenerator <- function(attributes, filter, values, mart){

  ## attributes
  if(length(attributes)==0)
    stop("Please specify at least one attribute")

  m = match(attributes, mart@datasets$attributes$attributes)
  if(any(is.na(m)))
    stop(sprintf("attribute(s) %s not found, please use the function 'listAttributes' to get valid attribute names",
                 paste(attributes[is.na(m)], collapse=", ")))
  Afield = mart@datasets$attributes$field[m]
  Atab   = mart@datasets$attributes$table[m]
  Akey   = mart@datasets$attributes$key[m]

  ## special treatment of main table
  isMain = (Atab=="main")
  m2 = match(Akey, mart@datasets$mainTables$keys)
  if(any(is.na(m2) & isMain))
    stop("Internal error: Key does not match MainTable key")
  Atable = ifelse(isMain, mart@datasets$mainTables$tables[m2], Atab)

  if(length(unique(Atable))>1)
    stop(paste("This query is currently not possible: your attributes extend over multiple tables (",
               paste(Atable, collapse=", "), "), please only use attributes from one table per query.", sep=""))
  
  ## filter
  ## FIXME: I assume this is correct - filter can have only length 1?
  if(length(filter)!=1)
    stop(sprintf("'length(filter)' must be 1."))
  
  m = match(filter, mart@datasets$filter$filter, nomatch = 0)
  if(any(is.na(m)))
    stop(sprintf("filter(s) %s not found, please use the function 'listFilters' to get valid filter names",
                 paste(filter[is.na(m)], collapse=", ")))
  Ffield = mart@datasets$filter$field[m]
  Ftab   = mart@datasets$filter$table[m]
  Fkey   = mart@datasets$filter$key[m]

  ## special treatment of main table
  isMain = (Ftab=="main")
  m2 = match(Fkey, mart@datasets$mainTables$keys)
  if(any(is.na(m2) & isMain))
    stop("Internal error: Key does not match MainTable key")
  Ftable = ifelse(isMain, mart@datasets$mainTables$tables[m2], Ftab)

  query = paste("SELECT DISTINCT",
    paste(Ftable, Ffield, sep=".", collapse=", "), ",",
    paste(Atable, Afield, sep=".", collapse=", "),
    "FROM")

  ########################################
  # Key Order: gene overrides transcript #
  ########################################
  ## FIXME: why only for Akey[1] and not for all?
  if(Fkey[1] == "gene_id_key" && Akey[1] == "transcript_id_key")  
    Akey[1] = Fkey[1]

  if(Akey[1] == "gene_id_key" && Fkey[1] == "transcript_id_key")
    Fkey[1] = Akey[1]
  
  tables <- unique(c(Atable, Ftable));
 
  ## FIXME: why only for tables[1:2] and Akey[1] and not for all?
  query = paste(query, tables[1])
  if(length(tables) > 1)
    query = paste(query, " INNER JOIN ",tables[2], " ON ", tables[1], ".", Akey[1], " = ",
      tables[2], ".", Fkey[1], sep="")  
  
  query = paste(query, " WHERE ", 
    paste(Ftab[1], Ffield[1], sep ="."),
    " IN (",
    paste("'", values, "'", collapse=", ", sep=""), ")", sep="")
    
  return(query)
}
