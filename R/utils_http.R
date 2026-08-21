.guess_port <- function(url) {
  if (startsWith(url, "http://")) {
    return("80")
  }
  return("443")
}
