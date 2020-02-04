.cran_packages <- c("data.table", 
                  "futile.logger",
                  "ggplot2",
                  "optparse",
                  "plyr",
                  "readr",
                  "reshape2",
                  "scales",
                  "viridis",
                  "yaml")
                  
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst])
}

sapply(.cran_packages, require, character.only = TRUE)
