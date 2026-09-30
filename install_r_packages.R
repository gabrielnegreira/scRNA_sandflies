#this script will install all packages required to run the analysis
#install and initiate renv
install.packages("renv")
renv::init(bare = TRUE, bioconductor = "3.23")

#create a function to read and format the files with the list of packages to install
read_packages <- function(x){
  x <- readLines(x, skipNul = TRUE)
  x <- x[!grepl("#", x)]
  x <- x[x != ""]
  return(x)
}

#get the list of packages to install
cran_packages <- read_packages("r_packages/cran_packages.yml")
bioconductor_packages <- read_packages("r_packages/bioconductor_packages.yml")
bioconductor_packages <- paste0("bioc::", bioconductor_packages) #needed for renv to install them

#install cran packages
renv::install(cran_packages)

#install bioconductor packages
renv::install(bioconductor_packages)