Mart <- setClass("Mart",
         representation(biomart = "character",
                        host = "character",
                        vschema = "character",
                        version = "character",
                        dataset = "character",
                        filters = "data.frame",
                        attributes = "data.frame",
                        http_config = "list"
                        ),
         prototype(dataset = "",
                   vschema="default",
                   version = "",
                   http_config = list()
                   )
         )
