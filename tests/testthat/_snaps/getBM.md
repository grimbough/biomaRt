# getBM doesn't convert T/F alleles into TRUE/FALSE

    Code
      getBM(attributes = c("refsnp_id", "minor_allele"), filters = c("chr_name",
        "start", "end"), values = list(8, 35127386, 35127386), mart = snp_ensembl)
    Output
        refsnp_id minor_allele
      1 rs1528723            T

