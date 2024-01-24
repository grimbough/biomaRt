structure(list(method = "POST", url = "https://asia.ensembl.org/info/website/archives/index.html?redirect=no", 
    status_code = 500L, headers = structure(list(Date = "Fri, 12 Jan 2024 15:37:37 GMT", 
        Server = "Apache", `X-Frame-Options` = "SAMEORIGIN", 
        `Content-Security-Policy` = "frame-ancestors 'self'", 
        `Transfer-Encoding` = "chunked", `Content-Type` = "text/plain; charset=utf-8"), class = "httr2_headers"), 
    body = charToRaw("ENSG00000000003\n"), request = structure(list(
        url = "https://asia.ensembl.org/info/website/archives/index.html?redirect=no", 
        method = NULL, headers = list(), body = list(data = list(
            query = structure("%3C%3Fxml%20version%3D%221.0%22%20encoding%3D%22UTF-8%22%3F%3E%0A%3C%21DOCTYPE%20Query%3E%0A%3CQuery%20%20virtualSchemaName%20%3D%20%22default%22%20formatter%20%3D%20%22TSV%22%20header%20%3D%20%220%22%20uniqueRows%20%3D%20%220%22%20count%20%3D%20%22%22%20datasetConfigVersion%20%3D%20%220.6%22%20%3E%0A%09%3CDataset%20name%20%3D%20%22hsapiens_gene_ensembl%22%20interface%20%3D%20%22default%22%20%3E%0A%09%09%3CFilter%20name%20%3D%20%22ensembl_gene_id%22%20value%20%3D%20%22ENSG00000000003%22%2F%3E%0A%09%09%3CAttribute%20name%20%3D%20%22ensembl_gene_id%22%20%2F%3E%0A%09%3C%2FDataset%3E%0A%3C%2FQuery%3E", class = "AsIs")), 
            type = "form", content_type = "application/x-www-form-urlencoded", 
            params = list()), fields = list(), options = list(
            timeout_ms = 10000), policies = list()), class = "httr2_request"), 
    cache = new.env(parent = emptyenv())), class = "httr2_response")
