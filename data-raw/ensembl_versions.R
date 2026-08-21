# It became necessary to create a data frame of Ensembl versions and their
# release dates because the website stopped displaying all versions.

ensembl_versions <- data.frame(
  version = c(
    "116",
    "115",
    "114",
    "113",
    "112",
    "111",
    "110",
    "109",
    "108",
    "107",
    "106",
    "105",
    "104",
    "103",
    "102",
    "101",
    "100",
    "99",
    "98",
    "97",
    "96",
    "95",
    "94",
    "93",
    "92",
    "91",
    "90",
    "89",
    "88",
    "87",
    "86",
    "80",
    "77",
    "75",
    "54"
  ),
  date = c(
    "June 2026",
    "Sep 2025",
    "May 2025",
    "Oct 2024",
    "May 2024",
    "Jan 2024",
    "Jul 2023",
    "Feb 2023",
    "Oct 2022",
    "Jul 2022",
    "Apr 2022",
    "Dec 2021",
    "May 2021",
    "Feb 2021",
    "Nov 2020",
    "Aug 2020",
    "Apr 2020",
    "Jan 2020",
    "Sep 2019",
    "Jul 2019",
    "Apr 2019",
    "Jan 2019",
    "Oct 2018",
    "Jul 2018",
    "Apr 2018",
    "Dec 2017",
    "Aug 2017",
    "May 2017",
    "Mar 2017",
    "Dec 2016",
    "Oct 2016",
    "May 2015",
    "Oct 2014",
    "Feb 2014",
    "May 2009"
  )
)

ensembl_versions$current_release <- ""
ensembl_versions$name <- paste("Ensembl", ensembl_versions$version)
ensembl_versions$url <- sprintf(
  "https://%s.archive.ensembl.org",
  sub(" ", "", tolower(ensembl_versions$date))
)

ensembl_versions <- ensembl_versions[, c(
  "name",
  "date",
  "url",
  "version",
  "current_release"
)]

usethis::use_data(ensembl_versions, overwrite = TRUE)
usethis::use_data(ensembl_versions, overwrite = TRUE, internal = TRUE)
