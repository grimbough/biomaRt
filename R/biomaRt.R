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
                        )#,
        # prototype(
        #           symbol = character(0),
        #           description = character(0),
        #           chromosome = character(0),
        #           band = character(0),
        #           start = numeric(0),
        #           end = numeric(0),
                #   GO = NULL,
                #   OMIM = NULL
         #          )
         );


setClass("MultiGene","Gene");


setClass("Mart",
         representation(
                        mart = "character",
                        connection = "MySQLConnection",
                        driver = "MySQLDriver"
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
####### connect to mart ######
#could make host changable in case people have a local BioMart install

connected <- FALSE;


martConnect <- function( mart = NULL, driver = NULL, host = NULL, dbname = NULL, user = NULL, passwd = NULL ){
  martObj <- NULL;
  driver <- NULL;
  
  if( !is.null( mart ) ){
    if( mart == "ensembl" ){
    
      driver <- dbDriver("MySQL", force.reload = FALSE);
      con <- dbConnect(driver,user="anonymous", host="ensembldb.ensembl.org", dbname="ensembl_mart_28_1");
      connected <- TRUE
      martObj <- new("Mart", mart = "ensembl", connection = con)
      writeLines("-  Connected to Ensembl  -");
      
    }
    
    if( mart == "vega" ){
      
      driver <- dbDriver("MySQL", force.reload = FALSE);
      con <- dbConnect(driver,user="anonymous", host="ensembldb.ensembl.org", dbname="vega_mart_28_1");
      connected <- TRUE
      martObj <- new("Mart", mart = "vega", connection = con)
      writeLines("-  Connected to vega  -");
      
    }

    if( mart == "snp" ){
      
      driver <- dbDriver("MySQL");
      con <- dbConnect(driver,user="anonymous", host="ensembldb.ensembl.org", dbname="snp_mart_28_1");
      connected <- TRUE
      martObj <- new("Mart", mart = "ensembl", connection = con)
      writeLines("-  Connected to snp  -");
      
    }    
  }
  return( martObj )
}

martDisconnect <- function( mart = NULL){

  dbDisconnect( mart@connection )

}



### initialisation ###

tablefile <- system.file("tables/ArrayEnsemblTable.txt",package="biomaRt") ;
arrayToEnsembl <- read.table(file = tablefile, sep="\t");
tablefile <- system.file("tables/ArraySpeciesTable.txt",package="biomaRt") ;
arrayToSpecies <- read.table(file = tablefile, sep="\t");

####### welcome message #############

####### local functions #############

mapArrayToEnsemblTable <- function(array = NULL, species = NULL){
  
  table <- "";
  table <- paste( species,"_gene_ensembl__xref_",arrayToEnsembl[ match( array, arrayToEnsembl[,1]),2 ],"__dm",sep="");

  return( table );
}

mapSpeciesToESpecies <- function( species = NULL, mart = NULL ){
  
  mSpecies <- NULL;
  if( mart@mart == "ensembl"){
    query <- paste("select mart_species from meta_release_info where species = '", species, "'",sep = "");
  
    res <- dbGetQuery(conn = mart@connection,statement = query );
    mSpecies <- res[1,1]
  }
  if(mart@mart == "vega"){
    mSpecies <- "hsapiens";
  }
  return( mSpecies ); 
}

mapESpeciesToGeneTable <- function( species = NULL, mart = NULL){

  table <- NULL;
  
  if(mart@mart == "ensembl"){
    table <- paste( species,"_gene_ensembl__gene__main",sep = "");
  }
  
  if(mart@mart == "vega"){
    table <- paste( species,"_gene_vega__gene__main",sep = "");
  }
  
  return( table );
}

mapESpeciesToGOTable <- function( species = NULL){

  table <- paste(species,"_gene_ensembl__xref_go__dm",sep = "");
  
  return(table); 
}

mapArrayToSpecies <- function( array = NULL ){

 species = "";
 species <- arrayToSpecies[ match( array, arrayToSpecies[,1]),2 ]; 
 return( species )
}

mapSpeciesToLocusLink <- function( species = NULL, mart = NULL ){

  table <- NULL;
  if(mart@mart == "ensembl"){
    table <- paste( species, "_gene_ensembl__xref_locuslink__dm", sep = "");
   }
  if(mart@mart == "vega"){
     table <- paste( species, "_gene_vega__xref_locuslink__dm", sep = "");
  }
  
  return( table );
}

mapSpeciesToRefSeq <- function( species = NULL, mart = NULL ){

  table <- NULL;
  
  if(mart@mart == "ensembl"){
    table <- paste( species,"_gene_ensembl__xref_refseq__dm" , sep = "");
  }
  if(mart@mart == "vega"){
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

getGene <- function( id = NULL, type = NULL, array = NULL, species = NULL, mart= NULL, output = "martTable" ){

  genes <- NULL
  IDTable <- NULL
  
  if(!is.null(array)){
   type <- "affy"
  }
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  if( is.null( mart )){
    stop("you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("ensembl error: invalid type choose either affy, refseq or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      if( mart@mart != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }
      species <- mapArrayToSpecies( array );
      Especies <- mapSpeciesToESpecies( species = species, mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, mart = mart );
      IDTable <- mapArrayToEnsemblTable( array = array, species = Especies );
      
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){
      
      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, mart = mart);
      speciesTable <- mapESpeciesToGeneTable( Especies, mart = mart );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
  
  
  if( type == "refseq"){
    if( !is.null( species ) ){
      
      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, mart = mart );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
  
  
  if( type == "embl"){
    if( !is.null( species ) ){
      if( mart@mart != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }
      
      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies, mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, mart = mart );
      
    }
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  
  if (mart@mart == "ensembl" || mart@mart == "vega"){
    if (length( id ) >= 1){
      
      res <- NULL;
      ids <- paste("'",id[1],"'",sep="")
      if(length(id) >= 2){
        for( i in 2:length(id)){
          ids <- paste(ids,",'",id[i],"'",sep="")
        }
      }
      query <- paste("select distinct ",IDTable,".display_id_list, ",IDTable,".gene_stable_id,",speciesTable,".display_id, description, band,chr_name, gene_chrom_start, gene_chrom_end  from ", IDTable ," inner join ",speciesTable," on ",IDTable,".gene_stable_id = ",speciesTable,".gene_stable_id where ",IDTable,".display_id_list in (",ids,")",sep="");
      res <- dbGetQuery( conn = mart@connection,statement = query);
                                        # if(dim(res)[1]!= 0){
      if( output == "martTable"){
        table <- new("martTable", id = res[,1], table = list(symbol = res[,3], description = res[,4], band = res[, 5], chromosome = res[,6], start = res[,7], end = res[,8], martID = res[,2]))
      }
      if( output == "geneObject" ){
        if(dim(ann)[1] == 0){
          genes[[1]] <- new("Gene", id = as.character(id));
          
        }
        else{
          for( j in 1:length(id)){
            m <- match(ann[,1], id[j],nomatch = 0)
            if(sum(m) == 0){
              genes[[j]] <- new("Gene", id = id[j]);
              
            }
            else{
              if(sum(m) == 1){
                genes[[j]] <- new("Gene",id = ann[m == 1, 1], martID = ann[m == 1,2], symbol = ann[ m == 1 ,3], description = ann[ m == 1,4], band = ann[ m == 1,5], chromosome = ann[ m == 1,6], start = as.numeric(ann[ m == 1,7]), end = as.numeric(ann[m == 1,8]))
              }
              else{
                
                genes[[j]] <- new("MultiGene",id = ann[m == 1,1], martID = ann[m==1,2],symbol = ann[m == 1,3], description = ann[m == 1,4], band = ann[m == 1,5], chromosome = ann[m == 1,6], start = as.numeric(ann[m == 1,7]), end = as.numeric(ann[ m == 1,8]))
                
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

getFeature <- function( symbol = NULL, array = NULL, type = NULL, mart = NULL ){

  if(!is.null(array)){
    type <- "affy"
  }
 
  if( is.null( type )){
    stop("you must provide the identifier type using the type argument")
  }
  if( is.null( mart )){
    stop("you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("ensembl error: invalid type choose either affy, refseq or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      if( mart@mart != "ensembl"){
        stop("you can only use ensembl when working with affy id's");
      }

      species <- mapArrayToSpecies( array );
      Especies <- mapSpeciesToESpecies( species = species, mart = mart );
      speciesTable <- mapESpeciesToGeneTable( Especies, mart = mart );
      
      IDTable <- mapArrayToEnsemblTable( array = array, species = Especies  );
      
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
 if (mart@mart == "ensembl"){
      
      query <- paste("select distinct gene_stable_id, display_id, description from ",speciesTable," where display_id like '%",symbol,"%' and affy_hg_u95av2_bool = 1",sep="");   
      martID <- dbGetQuery( conn = mart@connection,statement = query);
      features <- NULL;
      mID <- NULL;
      sym <- NULL;
      description <- NULL;
      
      if(dim(martID)[1] != 0){
       
        for(j in 1: dim(martID)[1]){
          feature <- NULL;  
          query <- paste("select distinct display_id_list from ",IDTable," where gene_stable_id = '",martID[j,1],"'",sep="");
          feature <- dbGetQuery( conn = mart@connection,statement = query);
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
}



###################
#get GO annotation#
###################

getGO <- function( id = NULL, type = NULL, array = NULL, species = NULL, mart = NULL ){
  
  go <- NULL

  if(!is.null(array)){
    type <- "affy"
  }
  
  if(mart@mart != "ensembl"){
    stop( "you should connect to ensembl mart for this query" );
  }
  
  if( is.null( type )){
    stop("ensembl error: you must provide the identifier type using the type argument")
  }
  
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("ensembl error: invalid type choose either affy or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      
      species <- mapArrayToSpecies( array );
      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      GOTable <- mapESpeciesToGOTable( Especies );
      IDTable <-  mapArrayToEnsemblTable( array = array, species = Especies);
        
    
 
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, mart = mart);
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
  
  
  if( type == "refseq"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, mart = mart );
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
  
  if( type == "embl"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies );
      GOTable <- mapESpeciesToGOTable( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  
  if( length( id ) == 1 ){
    
    ann <- NULL
    query <- paste("select gene_stable_id from ", IDTable ," where display_id_list ='", id,"'",sep="");
    geneID <- dbGetQuery(conn = mart@connection, statement = query);
    query <- paste("select dbprimary_id, description, evidence_code from ", GOTable ," where gene_stable_id ='", geneID[1,1],"'",sep="");
    res <- dbSendQuery(conn = mart@connection, statement = query);
    ann <- fetch(res);
    dbClearResult( res );
    go <- new("GO", GOID = as.vector(ann[,1]), description = as.vector(ann[,2]), evidence = as.vector(ann[,3]))
  }
    
  if( length( id ) > 1){
    for( i in 1:length( id )){
   
      ann <- NULL
      query <- paste("select gene_stable_id from ", IDTable ," where display_id_list ='", id,"'",sep="");
      geneID <- dbGetQuery(conn = mart@connection, statement = query);
      query <- paste("select dbprimary_id, description, evidence_code from ", GOTable ," where gene_stable_id ='", geneID[1,1],"'",sep="");
      res <- dbSendQuery(conn = mart@connection, statement = query);
      ann <- fetch( res );
      dbClearResult( res );
      
      go[[i]] <- new("GO", GOID = as.vector(ann[,1]), description = as.vector(ann[,2]), evidence = as.vector(ann[,3]));
    }
  }
  
  return( go );
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

  if(mart@mart != "ensembl"){
    stop( "you should connect to ensembl mart for this query" );
  }
  
  if( is.null( type )){
    stop("ensembl error: you must provide the identifier type using the type argument")
  }
  if(!(type == "affy") && !(type == "locuslink") && !(type == "refseq") && !(type == "embl")){
    stop("ensembl error: invalid type choose either affy or locuslink");
  }
  
  if( type == "affy" ){
    if( !is.null( array ) ){
      IDTable <-  mapArrayToEnsemblTable( array, species = "hsapiens");
 
    }
    else{
      stop( "you must provide the affymetrix array identifier via the array argument when using this function for affy identifiers" );
    }
  }
  
  
  if( type == "locuslink"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToLocusLink( Especies, mart = mart );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for locuslink identifiers" );
    }
  }
 
  if( type == "refseq"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToRefSeq( Especies, mart = mart);
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for RefSeq identifiers" );
    }
  }
 
  if( type == "embl"){
    if( !is.null( species ) ){

      Especies <- mapSpeciesToESpecies( species = species, mart = mart);
      IDTable <- mapSpeciesToEMBL( Especies );
   
    }
    
    else{
      stop( "you must provide the species via the species argument when using this function for EMBL identifiers" );
    }
  }
  if ( length( id ) == 1 ){

    ann <- NULL
    query <- paste("select gene_id_key from ", IDTable ," where display_id_list = '", id,"'",sep="");
    geneID <- dbGetQuery(conn = mart@connection, statement = query);
    query <- paste("select disease, omim_id from ", OMIMTable ," where gene_id_key ='", geneID[1,1],"'",sep="");
    res <- dbSendQuery(conn = mart@connection, statement = query);
    ann <- fetch( res );
    dbClearResult( res );
    omim <- new("OMIM",OMIMID = as.vector(ann[,2]), disease = as.vector(ann[,1]));
  
  }
  if( length( id ) > 1){
    
    for( i in 1:length(id)){
      
      ann <- NULL
      query <- paste("select gene_id_key from ", IDTable ," where display_id_list = '", id,"'",sep="");
      geneID <- dbGetQuery(conn = mart@connection, statement = query);
      query <- paste("select disease, omim_id from ", OMIMTable ," where gene_id_key ='", geneID[1,1],"'",sep="");
      res <- dbSendQuery(conn = mart@connection, statement = query);
      ann <- fetch( res );
      dbClearResult( res );
      omim[[i]] <- new("OMIM",OMIMID = as.vector( ann[,2] ), disease = as.vector( ann[,1] ));
    }
  }
  return( omim );  
}

################
#show affy info#
################

getAffyArrays <- function(){
  print( arrayToSpecies );
}


getSpecies <- function( mart = NULL ){
  
  query <- paste("select species from meta_release_info");
  res <- dbGetQuery( mart@connection, query );
  m <- match( "multi_species", res$species );
  sel<- rep(TRUE, length(res$species));
  sel[m] <- FALSE;
  res <- res[sel,]
  return( res ); 
}








