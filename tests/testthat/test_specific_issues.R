# ## the intention here is to test specific issues reported on the bioc support site
# 
# context('Testing specific reported issues')
# mart <- useMart(host="feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_SNP")
# ensembl <- useDataset("hsapiens_snp", mart=mart)
# res <- getBM(attributes=c("chr_name","minor_allele"),
#              filters=c("chr_name", "chrom_start", "chrom_end"), values=list("11", 108202735,108202740), mart=ensembl)
