setClass("Mart",
         representation(mysql = "logical",
                        connections = "list",
                        mysqldriver = "list",
                        mainTables = "list",
                        biomart = "character",
                        host = "character",
                        port = "numeric",
                        vschema = "character",
                        dataset = "character",
                        filters = "environment",
                        attributes = "environment",
                        attributePointer = "environment",
                        archive = "logical"
                        ),
         prototype(mysql = FALSE,
                   connections = new("list"),
                   dataset = "",
                   vschema="default",
                   archive = FALSE
                   )
         );

