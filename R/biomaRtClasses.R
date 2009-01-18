setClass("Mart",
         representation(biomart = "character",
                        host = "character",
                        vschema = "character",
                        dataset = "character",
                        filters = "data.frame",
                        attributes = "data.frame",
                        archive = "logical"
                        ),
         prototype(dataset = "",
                   vschema="default",
                   archive = FALSE
                   )
         );

