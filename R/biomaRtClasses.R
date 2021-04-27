Mart <- setClass("Mart",
         representation(biomart = "character",
                        host = "character",
                        vschema = "character",
                        version = "character",
                        dataset = "character",
                        filters = "data.frame",
                        attributes = "data.frame",
                        httr_config = "list"
                        ),
         prototype(dataset = "",
                   vschema="default",
                   version = "",
                   httr_config = list(httr::config())
                   )
         )
