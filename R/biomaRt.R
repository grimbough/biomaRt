.packageName <- "biomaRt"

####### Class creation #######

setClass("OMIM",
         representation(
                        OMIMID = "vector",
                        disease = "vector"
                        ),
         prototype(
                   OMIMID = NULL,
                   disease = NULL
                   )
         );
                 

setClass("GO",
         representation(
                        GOID = "vector",
                        description = "vector",
                        evidence = "vector"
                        )
         ,
         prototype(
                   GOID = NULL,
                   description = NULL,
                   evidence = NULL
                   )
         );

setClass("Gene",
         representation(
                        id = "character",
                        martID= "character",
                        symbol = "character",
                        description = "character",
                        chromosome = "character",
                        band = "character",
                        start = "numeric",
                        end = "numeric",
                        GO = "GO",
                        OMIM = "OMIM"
                        )
         );
         

setClass("MultiGene","Gene");


setClass("Mart",
         representation(
                        ensembl = "MySQLConnection",
                        vega = "MySQLConnection",
                        sequence = "MySQLConnection",
                        snp = "MySQLConnection",
                        arrayToSpecies = "data.frame",
                        arrayToEnsembl = "data.frame"
                        )
        # ,
        # prototype(
        #           mart = NULL,
                 #  connection = NULL,
                 #  driver = NULL
         #          )
        );


setClass("Feature",
         representation(
                        id = "character",
                        martID= "character",
                        symbol = "character"
                        ));


setClass("martTable",
         representation(id = "character",
                        table = "list"
                        ));

setMethod("show","martTable",
          function(object){
            cat("object of class martTable\n\n")
            cat("slot id\n\n")
            len <- length(object@id)
            if(len < 10){
                print(object@id)
            }
            else{
              print(object@id[1:10])
            }
            cat("\nslot table\n\n")
            names <- rownames(summary(object@table))
            if(len < 10){
              for(i in 1:length(names)){
                cat(names[i],"\n\n")
                print(object@table[[i]])
                cat("\n\n")
              }
            }
            else{
             for(i in 1:length(names)){
                cat(names[i],"\n\n")
                print(object@table[[i]][1:10])
                cat("\n\n")
              }
            }            
          }
          );


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

listMarts <- function(){
  
  driv <- dbDriver("MySQL", force.reload = FALSE);
  connection <- dbConnect(driv, user="anonymous", host="martdb.ebi.ac.uk");
  res <- dbGetQuery(connection,"show databases");
  marts <- c("ensembl_mart_","snp_mart_","sequence_mart_","vega_mart_");
  databases <- list( ensembl = "", snp = "", sequence = "", vega = "" );

  
  #Search latest releases of marts
      
      
  for(i in 1: length(marts)){
    matches <- grep(marts[i],res[,1]);
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
      databases[[i]] <- res[matches[latest],1];
    }
    else{
      databases[[i]] <- res[matches,1];
    }        
  }
  
  dbDisconnect(connection);
  
  return( databases );
  
}

##################
#Connect to marts#
##################

martConnect <- function(){
  
  databases <- listMarts();
  driver <- dbDriver("MySQL", force.reload = FALSE);
  
  tablefile <- system.file("tables/ArrayEnsemblTable.txt",package="biomaRt") ;
  arrayToEnsembl <- read.table(file = tablefile, sep="\t");
  tablefile <- system.file("tables/ArraySpeciesTable.txt",package="biomaRt") ;
  arrayToSpecies <- read.table(file = tablefile, sep="\t");
  
  mart <- new("Mart",
              ensembl = dbConnect(drv = driver,user="anonymous", host="ensembldb.ensembl.org", dbname = databases$ensembl),
              vega = dbConnect(drv = driver,user="anonymous", host="ensembldb.ensembl.org", dbname = databases$vega),
              sequence = dbConnect(drv = driver,user="anonymous", host="ensembldb.ensembl.org", dbname = databases$sequence),
              snp = dbConnect(drv = driver,user="anonymous", host="ensembldb.ensembl.org", dbname = databases$snp),
              arrayToSpecies = arrayToSpecies,
              arrayToEnsembl = arrayToEnsembl
              );
  writeLines(paste("-  Connected to ",databases$ensembl,", ",databases$vega,", ",databases$snp," and ",databases$sequence," -", sep = ""));
 
  return( mart )
}


######################
#Disconnect from mart#
######################

martDisconnect <- function( mart = NULL){

  dbDisconnect( mart@ensembl );
  dbDisconnect( mart@vega );
  dbDisconnect( mart@sequence );
  dbDisconnect( mart@snp );
}



### initialisation ###


####### welcome message #############

####### local functions #############

mapArrayToEnsemblTable <- function(array = NULL, species = NULL, mart = NULL, dbtable = NULL){
  
  table <- "";
  if(dbtable == "xrefdm"){
    table <- paste( species,"_gene_ensembl__xref_",mart@arrayToEnsembl[ match( array, mart@arrayToEnsembl[,1]),2 ],"__dm",sep="");
  }
  if(dbtable == "affybool"){
    table <- paste(mart@arrayToEnsembl[ match( array, mart@arrayToEnsembl[,1]),2 ],"_bool",sep="");
  }

  return( table );
}

mapSpeciesToESpecies <- function( species = NULL,db = NULL , mart = NULL ){

  mSpecies <- NULL;
  
  if( db == "ensembl"){
    query <- paste("select mart_species from meta_release_info where species = '", species, "'",sep = "");
  
    res <- dbGetQuery(conn = mart@ensembl,statement = query );
    mSpecies <- res[1,1]
  }
  if(db == "vega"){
    mSpecies <- "hsapiens";
  }
  if(db == "sequence" || db == "snp"){
  
    if(species == "homo_sapiens"){
      mSpecies <- "hsapiens"
    }
    if(species == "mus_musculus"){
      mSpecies <- "mmusculus"
    }
    if(species == "rattus_norvegicus"){
      mSpecies <- "rnorvegicus"
    }
    if(species == "apis_mellifera"){
      mSpecies <- "amellifera"
    }
    if(species == "caenorhabditis_elegans"){
      mSpecies <- "celegans"
    }
    if(species == "gallus_gallus"){
      mSpecies <- "ggallus"
    }
    
    
  }
  
  return( mSpecies ); 
}

mapESpeciesToGeneTable <- function( species = NULL, db = NULL){

  table <- NULL;
  
  if(db == "ensembl"){
    table <- paste( species,"_gene_ensembl__gene__main",sep = "");
  }
  
  if(db == "vega"){
    table <- paste( species,"_gene_vega__gene__main",sep = "");
  }
  
  return( table );
}

mapESpeciesToGOTable <- function( species = NULL){

  table <- paste(species,"_gene_ensembl__xref_go__dm",sep = "");
  
  return(table); 
}

mapArrayToSpecies <- function( array = NULL, mart = NULL ){

 species = "";
 species <- mart@arrayToSpecies[ match( array, mart@arrayToSpecies[,1]),2 ]; 
 return( species )
}

mapSpeciesToLocusLink <- function( species = NULL, db = NULL ){

  table <- NULL;
  if(db == "ensembl"){
    table <- paste( species, "_gene_ensembl__xref_locuslink__dm", sep = "");
   }
  if(db == "vega"){
     table <- paste( species, "_gene_vega__xref_locuslink__dm", sep = "");
  }
  
  return( table );
}

mapSpeciesToRefSeq <- function( species = NULL, db = NULL ){

  table <- NULL;
  
  if(db == "ensembl"){
    table <- paste( species,"_gene_ensembl__xref_refseq__dm" , sep = "");
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


######### public functions ##############

######################
#get gene information#
######################

getGene <- function( id = NULL, type = NULL, array = NULL, species = NULL, db = "ensembl", mart = NULL, output = "martTable" ){

  genes <- NULL
  IDTable <- NULL
  
  if(!is.null(array)){
   type <- "affy"
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  if( is.null( mart )){
    stop("you must provide a mart connection object, create with function martConnect")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("invalid type choose either affy, refseq, embl or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      if( db != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }
      species <- mapArrayToSpecies( array = array, mart = mart );
      Especies <- mapSpeciesToESpecies( species = species, db = db, mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, db = db);
      IDTable <- mapArrayToEnsemblTable( array = array, species = Especies, mart = mart, dbtable = "xrefdm" );
      
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){
      
      Especies <- mapSpeciesToESpecies( species = species, db = db, mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, db = db);
      speciesTable <- mapESpeciesToGeneTable( Especies, db = db );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
  
  
  if( type == "refseq"){
    if( !is.null( species ) ){
      
      Especies <- mapSpeciesToESpecies( species = species, db = db, mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, db = db );
      speciesTable <- mapESpeciesToGeneTable( Especies, db = db );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
  
  
  if( type == "embl"){
    if( !is.null( species ) ){
      
      Especies <- mapSpeciesToESpecies( species = species, db = db, mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies, db = db );
      speciesTable <- mapESpeciesToGeneTable( Especies, db = db );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  
  if (db == "ensembl" || db == "vega"){
    conn <- NULL
    
    if (length( id ) >= 1){
      
      res <- NULL;
      ids <- paste("'",id[1],"'",sep="")
      if(length(id) >= 2){
        for( i in 2:length(id)){
          ids <- paste(ids,",'",id[i],"'",sep="")
        }
      }
      query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",speciesTable,".display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end  from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,")",sep="");

      if( db == "ensembl"){
        res <- dbGetQuery( conn = mart@ensembl,statement = query);
      }
      
      if( db == "vega"){
        res <- dbGetQuery( conn = mart@vega,statement = query);
      }
      
      if( output == "martTable"){
        
        
        if(dim(res)[1] == 0){
          table <- new("martTable", id = id, table = list(symbol = NA, description = NA, band = NA, chromosome = NA, start = NA, end = NA, martID = NA))
        }
        else{
          foundID <- NULL
          symbol <- NULL
          description <- NULL
          band <- NULL
          chromosome <- NULL
          start <- NULL
          end <- NULL
          martID <- NULL
          
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
      if( output == "geneObject" ){
        if(dim(res)[1] == 0){
          genes[[1]] <- new("Gene", id = as.character(id));
          
        }
        else{
          for( j in 1:length(id)){
            m <- match(res[,1], id[j],nomatch = 0)
            if(sum(m) == 0){
              genes[[j]] <- new("Gene", id = id[j]);
              
            }
            else{
              if(sum(m) == 1){
                genes[[j]] <- new("Gene",id = res[m == 1, 1], martID = res[m == 1,2], symbol = res[ m == 1 ,3], description = res[ m == 1,4], band = res[ m == 1,5], chromosome = res[ m == 1,6], start = as.numeric(res[ m == 1,7]), end = as.numeric(res[m == 1,8]))
              }
              else{
                
                genes[[j]] <- new("MultiGene",id = res[m == 1,1], martID = res[m==1,2],symbol = res[m == 1,3], description = res[m == 1,4], band = res[m == 1,5], chromosome = res[m == 1,6], start = as.numeric(res[m == 1,7]), end = as.numeric(res[ m == 1,8]))
                
              }
            }     
          }
        } 
      }
    }
  } 
  return(table);  
}

########################
#get feature           #
########################

getFeature <- function( symbol = NULL, array = NULL, type = NULL, mart = NULL){

  if(!is.null(array)){
    type <- "affy"
  }
 
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("invalid type choose either affy, refseq or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){

      species <- mapArrayToSpecies( array = array, mart = mart );
      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, db = "ensembl" );
      
      IDTable <- mapArrayToEnsemblTable( array = array, species = Especies, mart = mart, dbtable = "xrefdm"  );
      affyBool <- mapArrayToEnsemblTable( array = array, species = Especies, mart = mart, dbtable = "affybool"  );
      
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
 
  
  query <- paste("select distinct gene_stable_id, display_id, description from ",speciesTable," where display_id like '%",symbol,"%' and ",affyBool," = 1",sep="");   
  martID <- dbGetQuery( conn = mart@ensembl,statement = query);
  features <- NULL;
  mID <- NULL;
  sym <- NULL;
  description <- NULL;
      
  if(dim(martID)[1] != 0){
    
    for(j in 1: dim(martID)[1]){
      feature <- NULL;  
      query <- paste("select distinct display_id_list from ",IDTable," where gene_stable_id = '",martID[j,1],"'",sep="");
      feature <- dbGetQuery( conn = mart@ensembl,statement = query);
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
  
  return( table )

}



###################
#get GO annotation#
###################

getGO <- function( id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL){
  
  go <- NULL

  if(!is.null(array)){
    type <- "affy"
  }
   
  if( is.null( type )){
    stop("ensembl error: you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("invalid type choose either affy,  locuslink refseq or embl");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      
      species <- mapArrayToSpecies( array = array, mart = mart );
      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      GOTable <- mapESpeciesToGOTable( Especies );
      IDTable <-  mapArrayToEnsemblTable( array = array, species = Especies, mart = mart, dbtable = "xrefdm");
        
    
 
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, db = "ensembl");
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
  
  
  if( type == "refseq"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, db = "ensembl" );
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
  
  if( type == "embl"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies );
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  
  if( length( id ) >= 1){
    ids <- paste("'",id[1],"'",sep="")
    if(length(id) >= 2){
      for( i in 2:length(id)){
        ids <- paste(ids,",'",id[i],"'",sep="")
      }
    }
    
    res <- NULL
    query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",GOTable,".dbprimary_id, description, evidence_code from ", IDTable ," inner join ",GOTable," on ",IDTable,".gene_stable_id = ",GOTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,")",sep="");
    
    res <- dbGetQuery(conn = mart@ensembl, statement = query);
    
    if(dim(res)[1] == 0){
      table <- new("martTable", id = id, table = list(GOID = NA, description = NA, evidence = NA, martID = NA))
    }
    else{
      foundID <- NULL
      GOID <- NULL
      description <- NULL
      evidence <- NULL
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
    }
    table <- new("martTable", id = foundID, table = list(GOID = GOID, description = description, evidence = evidence, martID = martID))
  }
  
  return( table );
}


#####################
#get OMIM annotation#
#####################

getOMIM <- function( id = NULL, type = NULL, array = NULL, mart = NULL){
  
  OMIMTable <- "hsapiens_gene_ensembl__disease__dm";
  omim <- NULL
  species <- "homo_sapiens";

  if(!is.null(array)){
    type <- "affy"
  }

  #if(mart@mart != "ensembl"){
  #  stop( "you should connect to ensembl mart for this query" );
  #}
  
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("invalid type choose either affy, refseq. embl or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      IDTable <-  mapArrayToEnsemblTable( array, species = "hsapiens", mart = mart, dbtable = "xrefdm");
 
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, db = "ensembl" );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
 
  if( type == "refseq"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, db = "ensembl");
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
 
  if( type == "embl"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, db = "ensembl", mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  if( length( id ) >= 1){
    ids <- paste("'",id[1],"'",sep="")
    if(length(id) >= 2){
      for( i in 2:length(id)){
        ids <- paste(ids,",'",id[i],"'",sep="")
      }
    }
    
    res <- NULL
    query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",OMIMTable,".omim_id, disease from ", IDTable ," inner join ",OMIMTable," on ",IDTable,".gene_id_key = ",OMIMTable,".gene_id_key where ",IDTable,".display_id_list in (",ids,")",sep="");
    res <- dbGetQuery(conn = mart@ensembl, statement = query);
    
    
    if(dim(res)[1] == 0){
      table <- new("martTable", id = id, table = list(OMIMID = NA, disease = NA, martID = NA))
    }
    else{
      foundID <- NULL
      OMIMID <- NULL
      disease <- NULL
      martID <- NULL
      
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
    }
    table <- new("martTable", id = foundID, table = list(OMIMID = OMIMID, disease = disease, martID = martID))
   
  }
  return( table );  
}

#####################
#getSequence
#####################


getSequence <- function(species = NULL, chromosome = NULL, start = NULL, end = NULL, martTable = NULL, mart = NULL){
  
  sequence <- NULL;
  
  Especies <- mapSpeciesToESpecies( species = species, db= "sequence", mart = mart);
  speciesTable <- paste( Especies,"__dna_chunks__main",sep="" );
  
  if(is.null( martTable ) && !is.null( chromosome ) && !is.null( start ) && !is.null( end )){
    for(i in 1:length( chromosome )){
      
      if(end[i] - start[i] > 100000){
        stop("maximum sequence length is 100000 nucleotides, change start and end arguments to make the sequence size smaller")
      }
      
      chunkStart <- (floor((start[i] - 1)/100000)*100000) + 1;
      chunkEnd <- (floor((end[i] - 1)/100000)*100000) + 1;
      
      if(chunkStart == chunkEnd ){  #we only need to get one sequence chunck of 100000 nucleotides
        
        query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
        chunkseq <- dbGetQuery(conn = mart@sequence, statement = query);
        newstart <- start[i] - (floor((start[i]-1)/100000) * 100000)
        newend <- end[i] - (floor((end[i]-1)/100000) * 100000)
        
        sequence <- c(sequence,substr(as.character(chunkseq), newstart, newend));
        
      }
      
      else{   #query sequence is on 2 sequence chuncks
        
        query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
        chunkseq1 <- dbGetQuery(conn = mart@sequence, statement = query);
        query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkEnd,"'",sep="");
        chunkseq2 <- dbGetQuery(conn = mart@sequence, statement = query);
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
          chunkseq <- dbGetQuery(conn = mart@sequence, statement = query);
          newstart <- start[i] - (floor((start[i]-1)/100000) * 100000)
          newend <- end[i] - (floor((end[i]-1)/100000) * 100000)
          
          sequence <- c(sequence,substr(as.character(chunkseq), newstart, newend));
          
        }
        
        else{   #query sequence is on 2 sequence chuncks
          
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkStart,"'",sep="");
          chunkseq1 <- dbGetQuery(conn = mart@sequence, statement = query);
          query <- paste("select sequence from ", speciesTable ," where chr_name = '", chromosome[i],"' and chr_start = '",chunkEnd,"'",sep="");
          chunkseq2 <- dbGetQuery(conn = mart@sequence, statement = query);
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

###################
#show species info#
###################

getSpecies <- function( db = "ensembl", mart = NULL ){
  
  query <- paste("select species from meta_release_info");
  res <- dbGetQuery( mart@ensembl, query );
  m <- match( "multi_species", res$species );
  sel<- rep(TRUE, length(res$species));
  sel[m] <- FALSE;
  res <- res[sel,]
  return( res ); 
}

#######################
#getSNP
#######################


getSNP <- function(species = NULL, chromosome = NULL, start = NULL, end = NULL, mart = NULL){
  
    
    Especies <- mapSpeciesToESpecies( species = species, db = "snp", mart = mart );
    ensemblTable <- paste(Especies,"_snp__snp__main", sep ="");
    query <- paste("select snp_id_key,snp_chrom_start, allele, tscid, ensemblcoding_bool, ensemblintronic_bool, ensembl5utr_bool, ensembl3utr_bool, ensemblsyn_bool from ",ensemblTable," where chr_name = '",chromosome,"' and snp_chrom_start >= ",start," and snp_chrom_start <= ",end,sep="")
    res <- dbGetQuery( mart@snp, query );
    table <- new("martTable", id = res$tscid ,table = list(snpStart = res[,2], allele = res$allele, coding = res[,5], intronic =  res[,6], syn =  res[,9],utr5 = res[,7], utr3= res[,8]))
    return( table );
 
}

####################
#export FASTA      #
####################

#exportFASTA <- function( martTable = NULL ){


#}


