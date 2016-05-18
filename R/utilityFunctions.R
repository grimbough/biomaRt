.convertToDataFrame <- function(postRes) {
    ## first break on newline 
    lines <- strsplit(postRes, "\n")[[1]]
    ## any tabs within two quotes should be flagged
    lines <- gsub(pattern = "(\"[[:print:]]+)\\t([[:print:]]+\")", x = lines, replacement = "\\1--t--\\2")
    ## then split each line by tabs
    cols <- data.frame(do.call("rbind", strsplit(lines, "\t")), stringsAsFactors = FALSE)
    for(i in 1:ncol(cols)) { 
        cols[,i] <- .convertNumbers(cols[,i]) 
        if(class(cols[,i]) == "character") {
            cols[,i] <- gsub(pattern = "--t--", x = cols[,i], replacement = "\t")
        }
    }
    return(cols)
}

## R & read.table have an annoying habit of converting any column with just Ts
## in it to a logical - this can be frustrating for returning alleles
## we default to making everything a character and then
## if a column can be cast to integer or numeric we do so,
## otherwise it remains a character vector
.convertNumbers <- function(x) {
    if(class(type.convert(x)) %in% c('integer', 'numeric')) 
        return(type.convert(x))
    else
        return(x)
}

