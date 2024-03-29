## Use this code to initiate your devstuff file
# usethis::use_data_raw("devstuff")

# Put this into your .Rprofile
usethis::use_description(fields = list(
  # Package = "prereviewr", # default, no need to specify
  # Version = "0.0.0.9000", # default, no need to specify
  Title = "Tools for Quantile-Parametrized Distribiutions",
  Description = "Define and sample from quantile and quantile-parameterized distributions.",
  `Authors@R` = 'c(person("Dmytro", "Perepolkin",
                   email = "dperepolkin@gmail.com", role = c("aut", "cre"),
                   comment = c(ORCID = "0000-0001-8558-6183")))',
  # License = "MIT + file LICENSE", # will add via separate command
  # Encoding = "UTF-8", # default, no need to specify
  # LazyData = "true" # default, no need to specify
  URL = "https://github.com/dmi3kno/qpd",
  BugReports = "https://github.com/dmi3kno/qpd/issues",
  Language =  "en" # keep last
)
)
usethis::use_mit_license("Dmytro Perepolkin") #after use_description
usethis::use_git()
usethis::use_github()
usethis::use_r("script")
usethis::use_readme_rmd()
usethis::use_lifecycle_badge("Experimental") # after use_readme_rmd
usethis::use_news_md()
#usethis::use_testthat()
#usethis::use_test() #after use_testthat
#usethis::use_travis()
#usethis::use_appveyor()
#usethis::use_pkgdown()
#usethis::use_pkgdown_travis()
#usethis::use_coverage()
