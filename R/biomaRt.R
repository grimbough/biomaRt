.packageName <- "biomaRt"


setClass("Mart",
         representation(
                        connections = "list",
                        mysqldriver = "MySQLDriver",
                        arrayToSpecies = "data.frame"
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
    df = do.call("data.frame", args=lapply(object@table, "[", 1:n))
    show(df)
  })


#setMethod("as.data.frame","martTable",
 #         function( object, row.names = NULL, optional = NULL){
  #          names <- rownames( summary( object@table ))
  #          frame <- new("data.frame", id = as.vector(object@id))
  #          for(i in 1: length(names)){
  #            frame<-cbind(as.vector(object@table[[names[i]]]))
  #          }
  #          return( frame )
  #                   
  #        }
  #        );



#######################
#list marts
#######################

listMarts <- function(mart = NULL, host = "ensembldb.ensembl.org", user = "anonymous", password = ""){
  
  driv <- dbDriver("MySQL", force.reload = FALSE);
  connection <- dbConnect(driv, user = user, host = host, password = password);
  res <- dbGetQuery(connection,"show databases like '%mart%'");
  database <- NULL
  
  #Search latest releases of marts
 # martConf <- NULL  #one mart vs multiple marts      
  if(dim(res)[1] >= 1){
 #   for(i in 1: length(marts)){
      matches <- grep(mart,res[,1]);
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
        database <- res[matches[latest],1];
      #  martConf <- c( martConf, v );
      }
      else{
        databases <- res[matches,1];
        v <- as.numeric(strsplit(res[matches,1],"_")[[1]][3]);
      #  martConf <- c( martConf, v );
      }
  #  }
    
  #  if(!(sum(martConf)/4 == martConf[[1]])){ # mean there is only one new mart
  #    
  #    databases$snp <- databases$ensembl
  #    databases$vega <- databases$ensembl
  #    databases$sequence <- databases$ensembl
      
  #  }
    
  }
  else{
   # if(dim(res)[1] == 1){
   #   if(grep("ensembl",res[1,1]) == 1){
   #     databases$ensembl <- res[1,1]
   #     databases$snp <- databases$ensembl
   #     databases$vega <- databases$ensembl
   #     databases$sequence <- databases$ensembl
   #   }
   #   else{
   #     stop("No Ensembl database found, if you have a local installation, make sure your database follows trhe following naming convention: ensembl_mart_30 where 30 is the release number");
   #   }
      
   # }
   # else{
      stop(paste("No biomaRt database for",mart," found, please check if host is set correctly"));
   # }
  }
  dbDisconnect(connection);
  
  return( database );
  
}

##################
#Connect to marts#
##################

martConnect <- function(biomarts = "ensembl", host = NULL, user = NULL, password = NULL, mart = NULL, local = FALSE){
  
    if(is.null(mart)){
      driver <- dbDriver("MySQL", force.reload = FALSE);
      mart <- new("Mart", mysqldriver = driver)
    }
    presentConnections <- length(mart@connections)
    for(i in 1: length(biomarts)){
      if(local){
        if(biomarts[i] == "ensembl" || biomarts[i] == "sequence" || biomarts[i] == "uniprot" || biomarts[i] == "snp" || biomarts[i] == "vega"){
          database <- listMarts(mart = biomarts[i], host = host[i], user = user[i], password = password[i]);
          mart@connections[[i + presentConnections]] <- dbConnect(drv = mart@mysqldriver,user = user[i], host = host[i] , dbname = database, password = password[i])
          writeLines(paste("connected to: ",database))
          names(mart@connections)[i + presentConnections] <- biomarts[i]
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
      colnames(affyTables) <- c('species','affy array');
      mart@arrayToSpecies <- as.data.frame(cbind(affyTables[,2],affyTables[,1]))
    }
    return( mart )
}


######################
#Disconnect from mart#
######################

martDisconnect <- function( mart = NULL){

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
}

####### local functions #############

mapArrayToEnsemblTable <- function(array = NULL, species = NULL, mart = NULL, dbtable = NULL){
  
  table <- NULL;
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

  table <- paste(species,"_gene_ensembl__xref_uniprot_accession__dm",sep = "");
  return(table); 
}


mapArrayToSpecies <- function( array = NULL, mart = NULL ){

 species = "";
 species <- mart@arrayToSpecies[ match( array, mart@arrayToSpecies[,1]),2 ]; 
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

######### public functions ##############

######################
#get gene information#
######################

getGene <- function( id = NULL, type = NULL, array = NULL, species = NULL, db = "ensembl", mart = NULL){

  genes <- NULL
  IDTable <- NULL
  
  if( is.null( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if( db != "ensembl" && db != "vega"){
    stop("you can only use ensembl or vega");
  }
     
  if(!is.null(array)){
    type <- "affy";
     if( db != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }

  if( type == "affy" ){
    if( is.null( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  
  if(is.null( species ) ){
    stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
  }
  else{
    speciesTable <- mapSpeciesToGeneTable( species, db = db);
  }
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo")  && !(type == "ensembl")){
    stop("invalid type choose either affy, refseq, embl, hugo, ensembl or entrezgene");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = "hsapiens", mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ));
  
  if (db == "ensembl" || db == "vega"){
    conn <- NULL
    
    if (length( id ) >= 1){
      
      ids <- paste("'",id,"'",sep="",collapse=",");
      
      if(type == "ensembl"){
         query <- paste("select distinct gene_stable_id, display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end  from ", speciesTable," where gene_stable_id in (",ids,")",sep="");
      }
      else{
        query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",speciesTable,".display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end  from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,")",sep="");
      }

      if(db == "vega"){
        if(match("vega",names(mart@connections), nomatch=0) == 0){
          stop("You need a connection to vega for this query, use martConnect and include 'vega' in your biomarts vector")
        }
        res <- dbGetQuery( conn = mart@connections$vega, statement = query);
      }
      if(db == "ensembl"){
        res <- dbGetQuery( conn = mart@connections$ensembl, statement = query);
      }
      if(dim(res)[1] == 0){
        
        table <- new("martTable", id = id, table = list(symbol = NA, description = NA, band = NA, chromosome = NA, start = NA, end = NA, martID = NA))
      }
      else{
        
        foundID <- NULL;
        symbol <- NULL;
        description <- NULL;
        band <- NULL;
        chromosome <- NULL;
        start <- NULL;
        end <- NULL;
        martID <- NULL;
        for( j in 1:length(id)){
          m <- match( res[,1],id[j], nomatch=0)
          if(sum(m) == 0){
            foundID <- c(foundID, as.character(id[j]));
            symbol <- c(symbol, NA);
            description <- c(description, NA);
            band <- c(band, NA);
            chromosome <- c(chromosome, NA);
            start <- c(start, NA);
            end <- c(end, NA);
            martID <- c(martID, NA);              
          }
          else{
            
            if(type == "ensembl"){
              foundID <- c(foundID, res[m == 1,1]);
              symbol <- c(symbol, res[m == 1,2]);
              description <- c(description, res[m == 1,3]);
              band <- c(band, res[m == 1,4]);
              chromosome <- c(chromosome, res[m == 1,5]);
              start <- c(start, res[m == 1,6]);
              end <- c(end, res[m == 1,7]);
              martID <- c(martID, res[m == 1,1]);
            }
            else{ 
              foundID <- c(foundID, res[m == 1,1]);
              symbol <- c(symbol, res[m == 1,3]);
              description <- c(description, res[m == 1,4]);
              band <- c(band, res[m == 1,5]);
              chromosome <- c(chromosome, res[m == 1,6]);
              start <- c(start, res[m == 1,7]);
              end <- c(end, res[m == 1,8]);
              martID <- c(martID, res[m == 1,2]);
            }
          }     
        }
        table <- new("martTable", id = foundID, table = list(symbol = symbol, description = description, band = band, chromosome = chromosome, start = start, end = end, martID = martID))
        
      }
    }
  }
  
  return(table);  
}

########################
#get feature           #
########################

getFeature <- function( symbol = NULL, array = NULL, type = NULL, mart = NULL){
  
  table <- NULL
  
  if(!is.null(array)){
    type <- "affy"
  }
  if( is.null( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl")){
    stop("invalid type choose either affy, refseq or entrezgene");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){

      species <- mapArrayToSpecies( array = array, mart = mart );
      speciesTable <- mapSpeciesToGeneTable( species, db = "ensembl" );
      IDTable <- mapArrayToEnsemblTable( array = array, species = species, mart = mart, dbtable = "xrefdm"  );
      affyBool <- mapArrayToEnsemblTable( array = array, species = species, mart = mart, dbtable = "affybool"  );
      
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }

  query <- paste("select distinct gene_stable_id, display_id, description from ",speciesTable," where display_id like '%",symbol,"%' and ",affyBool," = 1",sep="");   
  martID <- dbGetQuery( conn = mart@connections$ensembl,statement = query);
  features <- NULL;
  mID <- NULL;
  sym <- NULL;
  description <- NULL;
  if(dim(martID)[1] != 0){
    
    for(j in 1: dim(martID)[1]){
      feature <- NULL;  
      query <- paste("select distinct display_id_list from ",IDTable," where gene_stable_id = '",martID[j,1],"'",sep="");
      feature <- dbGetQuery( conn = mart@connections$ensembl,statement = query);
      if(dim(feature)[1] != 0){
        feature <- feature[,1];
        select <- is.na(feature);
        features <- c(features,feature[!select]);
        mID <- c(mID,rep(martID[j,1],sum(!select)))
        sym <- c(sym,rep(martID[j,2],sum(!select)))
        description  <- c(description,rep(martID[j,3],sum(!select))) 
      }
      
    }
    
    table <- new("martTable", id = features, table = list(symbol = sym, description = description, martID= mID)) 
  }
  else{
    writeLines("No match found")
  }
  
  return( table )

}

###################
#get GO annotation#
###################

getGO <- function( id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL){

  table <- NULL;
  go <- NULL;

  if( is.null( mart ) || class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
 
  if(!is.null(array)){
    type <- "affy";
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }

  if( type == "affy" ){
    if( is.null( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  
  if(is.null( species ) ){
    stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
  }
  else{
    GOTable <- mapSpeciesToGOTable( species );
  }

  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo") && !(type == "ensembl")){
    stop("invalid type choose either affy,  entrezgene, hugo, ensembl, refseq or embl");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = "hsapiens", mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ));
  
  if( length( id ) >= 1){
    ids <- paste("'",id,"'",sep="",collapse=",")
    res <- NULL
    if(type == "ensembl"){
      query <- paste("select distinct ",IDTable,".gene_stable_id,",IDTable,".gene_stable_id, ",GOTable,".dbprimary_id,",GOTable,".description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".gene_stable_id in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
    }
    else{
      query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",GOTable,".dbprimary_id, description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,") and ",GOTable,".dbprimary_id != 'NULL'",sep="");
  }
    res <- dbGetQuery(conn = mart@connections$ensembl, statement = query);
    
    if(dim(res)[1] == 0){
      table <- new("martTable", id = id, table = list(GOID = NA, description = NA, evidence = NA, martID = NA))
    }
    else{
      foundID <- NULL;
      GOID <- NULL;
      description <- NULL;
      evidence <- NULL;
      martID <- NULL
      for( j in 1:length(id)){
        m <- match( res[,1],id[j], nomatch=0)
        if(sum(m) == 0){
          foundID <- c(foundID, as.character(id[j]));
          GOID  <- c(GOID, NA) 
          description <- c(description, NA);
          evidence <- c(evidence, NA);
          martID <- c(martID, NA);              
        }
        else{
          foundID <- c(foundID, res[m == 1,1]);
          GOID  <- c(GOID, res[m == 1,3])
          description <- c(description, res[m == 1,4]);
          evidence <- c(evidence, res[m == 1,5]);
          martID <- c(martID, res[m == 1,2]);
        }     
      }
      table <- new("martTable", id = foundID, table = list(GOID = GOID, description = description, evidence = evidence, martID = martID));
    }
  }
  
  return( table );
}


#####################
#get OMIM annotation#
#####################

getOMIM <- function( id = NULL, type = NULL, array = NULL, mart = NULL){
  
  OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
  omim <- NULL;
  table <- NULL;
  species <- "hsapiens";

  if( is.null( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if(!is.null(array)){
    type <- "affy";
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }

  if( type == "affy" ){
    if( is.null( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  else{
    if( is.null( species ) ){
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

  if( length( id ) >= 1){
    ids <- paste("'",id,"'",sep="",collapse=",")
    res <- NULL
    if(type == "ensembl"){
      query <- paste("select distinct ",IDTable,".gene_stable_id, ",IDTable,".gene_stable_id,",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".gene_stable_id in (",ids,")",sep="");
      
    }
    else{
      query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".display_id_list in (",ids,")",sep="");
  }
    res <- dbGetQuery(conn = mart@connections$ensembl, statement = query);
    
    
    if(dim(res)[1] == 0){
      table <- new("martTable", id = id, table = list(OMIMID = NA, disease = NA, martID = NA))
    }
    else{
      foundID <- NULL;
      OMIMID <- NULL;
      disease <- NULL;
      martID <- NULL;
      
      for( j in 1:length(id)){
        m <- match( res[,1],id[j], nomatch=0)
        if(sum(m) == 0){
          foundID <- c(foundID, as.character(id[j]));
          OMIMID  <- c(OMIMID, NA) 
          disease <- c(disease, NA);
          martID <- c(martID, NA);              
        }
        else{
          
          foundID <- c(foundID, res[m == 1,1]);
          OMIMID  <- c(OMIMID, res[m == 1,3])
          disease <- c(disease, res[m == 1,4]);
          martID <- c(martID, res[m == 1,2]);
        }     
      }
      table <- new("martTable", id = foundID, table = list(OMIMID = OMIMID, disease = disease, martID = martID))
    }
  }
  return( table );  
}

#####################
#getSequence
#####################


getSequence <- function(species = NULL, chromosome = NULL, start = NULL, end = NULL, martTable = NULL, mart = NULL){
  if( is.null( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("sequence",names(mart@connections),nomatch=0) == 0){
    stop("You are missing a database connection to sequence BioMart for this query.  Add connection to a sequence BioMart to your Mart object via the function martConnect. Use the following command to do this:  martConnect(biomarts='sequence',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password=''). ")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }

  sequence <- NULL;  
  speciesTable <- paste( species,"__dna_chunks__main",sep="" );
  
  if(is.null( martTable ) && !is.null( chromosome ) && !is.null( start ) && !is.null( end )){
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
    if(!is.null( martTable )){
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

getAffyArrays <- function(mart = NULL){
  print( mart@arrayToSpecies );
}


#######################
#getSNP
#######################


getSNP <- function(species = NULL, chromosome = NULL, start = NULL, end = NULL, mart = NULL){
  if( is.null( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if(match("snp",names(mart@connections),nomatch=0) == 0){
    stop("You are missing a database connection to snp for this query.  Add connection to snp a BioMart to your Mart object via the function martConnect.   Use the following command to do this:  martConnect(biomarts='snp',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password='')")
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


getHomolog <- function(id = NULL, from.type = NULL, to.type = NULL, from.species = NULL, to.species = NULL, mart = NULL, output = "martTable") {
  
  #all homology needs to go through ensembl_mart, not vega, etc.
  db = "ensembl";
  if( is.null( mart )|| class( mart ) != 'Mart'){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  
  id <- as.character(id)
  
  if (is.null(from.type)) {
    stop("You must provide the identifier type using the from.type  argument")
  }
  if (is.null(to.type)) {
    stop("You must provide the identifier type using the to.type  argument")
  }
  if (is.null(from.species)) {
    stop("You must provide the species to map FROM using the  from.species argument")
  }
  if (is.null(to.species)) {
    stop("You must provide the species to map TO using the to.species  argument")
  }
  
  fromIDTable <-
    switch(from.type,
           entrezgene = mapSpeciesToEntrezGene(from.species,db=db),
           refseq     = mapSpeciesToRefSeq(from.species,db=db),
           embl       = mapSpeciesToEMBL(from.species),
           hugo       = mapSpeciesToHUGO(from.species,db=db),
           ensembl    = mapSpeciesToGeneTable(to.species,db=db)
           );
  
  toIDTable <-
    switch(to.type,
           entrezgene = mapSpeciesToEntrezGene(to.species,db=db),
           refseq     = mapSpeciesToRefSeq(to.species,db=db),
           embl       = mapSpeciesToEMBL(to.species),
           hugo       = mapSpeciesToHUGO(from.species,db=db),
           ensembl    = mapSpeciesToGeneTable(to.species,db=db)
           );

  homolTable <- mapSpeciesToHomologTable(from.species,to.species);
  
  if (length(id) >= 1) {
    ids <- paste("'",id,"'",sep="",collapse=",")
    res <- NULL
    query <- paste("select distinct a.display_id_list,b.display_id_list from  ",
                   fromIDTable," as a inner join ",
                   homolTable," as c on a.gene_id_key=c.gene_id_key inner  join ",
                   toIDTable," as b on b.gene_id_key=c.homol_id ",
                   "where a.display_id_list in (",ids,")",sep="");
    
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


getPossibleXrefs <-  function( mart = NULL) {
  if ( is.null( mart ) || class( mart ) != 'Mart') {
    stop('You must supply a valid mart object.  You can use martConnect to produce one')
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
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


getSpecies <- function(mart = NULL, db = c("ensembl")) {
  if (is.null(mart) || class(mart)!='Mart') {
    stop('You must supply a valid mart object.  You can use martConnect to produce one')
  }
   if(match("ensembl",names(mart@connections), nomatch=0) == 0){
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

getXref <- function( id  = NULL, from.species = NULL, to.species = NULL, from.xref = NULL, to.xref = NULL, db = c('ensembl'), mart = NULL) {
  
  if (is.null(mart) || class(mart)!='Mart'){
    stop('You must specify a valid Mart object which can be created using martConnect().');
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  } 
  if (is.null(id)) {
    stop('You need to give ID(s) to map from');
  }
  if (is.null(from.species)) {
    writeLines('fromspecies is a necessary argument\nValid species for this mart are:');
    print(getSpecies(mart = mart, db = "ensembl"));
    stop();
  }
  
  if (is.null(from.xref) || is.null(to.xref)) {
    writeLines('Both fromxref and toxref are required.\nPossible crossreferences are:');
    print(getPossibleXrefs(mart = mart));
    stop();
  }
  
  xp <- paste(id, collapse="','")
  
  if (is.null( to.species )) {
    query <- paste( 
                   'SELECT distinct a.display_id_list as fromid,b.display_id_list as toid,',
                   'a.gene_stable_id as ensemblgene ',
                   'from ',from.species,'_gene_ensembl__xref_',from.xref,'__dm as a ',
                   'left join ',from.species,'_gene_ensembl__xref_',to.xref,'__dm as b ',
                   "on a.gene_id_key=b.gene_id_key where a.display_id_list in ('",xp,
                   "')",sep="");
  } else {
    if (db == 'vega'){
      stop('VEGA supports only hsapiens')
    }
    query <- paste(
                   'SELECT distinct a.display_id_list as fromid,b.display_id_list as toid,c.* ',
                   'from ',from.species,'_gene_ensembl__xref_',from.xref,'__dm as a ',
                   'left join ',from.species,'_gene_ensembl__homologs_',to.species,'__dm as c ',
                   'on a.gene_id_key=c.gene_id_key ',
                   'left join ',to.species,'_gene_ensembl__xref_',to.xref,'__dm as b ',
                   "on c.homol_id=b.gene_id_key where a.display_id_list in ('",xp,
                   "')",sep="");
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

#exportFASTA <- function( martTable = NULL ){


#}

#####################
#INTERPRO functions #
#####################



getINTERPRO <- function(id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL){
  
  genes <- NULL
  speciesTable <- NULL
  if(match("uniprot",names(mart@connections),nomatch=0) == 0){
    stop("You are missing a database connection to uniprot for this query.  Add connection to uniprot a BioMart to your Mart object via the function martConnect.  Use the following command to do this:  martConnect(biomarts='uniprot',mart=mart,host='martdb.ebi.ac.uk',user='anonymous',password='')")
  }
  db <- "ensembl";
  results <- FALSE;
    
  if( is.null( mart )|| class(mart)!='Mart'){
    stop("you must provide a valid Mart object, create with function martConnect")
  }
  if(match("ensembl",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to ensembl for this query, use martConnect and include 'ensembl' in your biomarts vector")
  }
  if(match("uniprot",names(mart@connections), nomatch=0) == 0){
    stop("You need a connection to uniprot for this query, use martConnect and include 'uniprot' in your biomarts vector")
  }
  if(!is.null(array)){
    type <- "affy";
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }

  if( type == "affy" ){
    if( is.null( array ) ){
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
    else{
      species <- mapArrayToSpecies( array = array, mart = mart );
    }
  }
  
  if(is.null( species ) ){
    stop( "you must provide the species via the species argument when using this function for Ensembl identifiers" );
  }
  else{
    uniprotTable <- mapSpeciesToUNIPROTTable( species );
  }
  
  if(!(type == "affy") && !(type == "entrezgene") && !(type == "refseq") && !(type == "embl") && !(type == "hugo")  && !(type == "ensembl")){
    stop("invalid type choose either affy, refseq, embl, hugo, ensembl or entrezgene");
  }

  IDTable <- switch(type,
                    affy = mapArrayToEnsemblTable( array, species = "hsapiens", mart = mart, dbtable = "xrefdm"),
                    ensembl =  mapSpeciesToGeneTable( species, db = "ensembl" ),
                    entrezgene = mapSpeciesToEntrezGene( species, db = "ensembl" ),
                    hugo = mapSpeciesToHUGO( species, db = "ensembl"),
                    refseq = mapSpeciesToRefSeq( species, db = "ensembl"),
                    embl = mapSpeciesToEMBL( species ));

  conn <- NULL
  
  if (length( id ) >= 1){
    
    ids <- paste("'",id,"'",sep="",collapse=",")
    if(type == "ensembl"){
      query <- paste("select distinct ",IDTable,".gene_stable_id,",uniprotTable,".display_id_list from ", IDTable ," inner join ",uniprotTable," on ",IDTable,".gene_stable_id = ",uniprotTable,".gene_stable_id where ",IDTable,".gene_stable_id in (",ids,") and ",uniprotTable,".display_id_list != 'NULL'",sep="");
    }
    else{
      query <- paste("select distinct ",IDTable,".display_id_list,",uniprotTable,".display_id_list from ", IDTable ," inner join ",uniprotTable," on ",IDTable,".gene_stable_id = ",uniprotTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,") and ",uniprotTable,".display_id_list != 'NULL'",sep="");
      
    }

    resEns <- dbGetQuery( conn = mart@connections$ensembl,statement = query);

    if(dim(resEns)[1] > 0){
      results <- TRUE
      if (length( resEns[[2]] ) >= 1){
      
        uniprotacc <- paste("'",resEns[[2]],"'",sep="",collapse=",")
        db <- "uniprot" 
        query <- paste("SELECT sptr_ac, ipro_id, name, short_name, type from uniprot__interpro__dm join uniprot__protein__main on uniprot__interpro__dm.protein_id_key = uniprot__protein__main.protein_id_key where sptr_ac in (",uniprotacc,")",sep="");
        resUni <- dbGetQuery( conn = mart@connections$uniprot,statement = query);
      }
    }
    if(dim(resEns)[1] == 0  || dim(resUni)[1] == 0){
      
      table <- new("martTable", id = id, table = list( interproid= NA, name = NA, shortname = NA, type = NA, martID = NA))
    }
    else{
      
      foundID <- NULL;
      interproid <- NULL;
      name <- NULL;
      shortname <- NULL;
      type <- NULL
      
      for( j in 1:length(id)){
        m <- match( resEns[,1],id[j], nomatch=0)
        if(sum(m) == 0){
          foundID <- c(foundID, as.character(id[j]));
          interproid <- c(interproid, NA);
          name <- c(name, NA);
          shortname <- c(shortname, NA);
          type <- c(type, NA);
        }
        else{
          uniprotAcc <- resEns[m == 1,2]
          if(length(uniprotAcc == 1)){ #one to one relation between Ensemble gene and protein in Uniprot
            unimatch <- match( resUni[,1],uniprotAcc, nomatch=0)
            foundID <- c(foundID, rep(resEns[m == 1,1],sum(unimatch)));
            interproid <- c(interproid, resUni[unimatch == 1,2]);
            name <- c(name, resUni[unimatch == 1,3]);
            shortname <- c(shortname, resUni[unimatch == 1,4]);
            type <- c(type, resUni[unimatch == 1,5]);
          }
          else{
            for(h in 1:length(uniprotAcc)){
              unimatch <- match( resUni[,1],uniprotAcc[h], nomatch=0)
              foundID <- c(foundID, rep(resEns[m == 1,1],sum(unimatch)));
              interproid <- c(interproid, resUni[unimatch == 1,2]);
              name <- c(name, resUni[unimatch == 1,3]);
              shortname <- c(shortname, resUni[unimatch == 1,4]);
              type <- c(type, resUni[unimatch == 1,5]);
            }
          }
        }     
      }
      table <- new("martTable", id = foundID, table = list(interproid= interproid, name = name, shortname = shortname, type = type))
      
    }
  }
  return(table)
}
