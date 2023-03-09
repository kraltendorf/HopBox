# hopbox - correlations between traits
# kayla altendorf

# load packages
library(dplyr)
library(data.table)
library(ggplot2)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(stringr)
library("tidyverse")
library(stringr)
library(corrplot)
library("Hmisc")
library(xtable)

# set working directory
setwd("/users/kayla.altendorf/OneDrive - USDA/Documents/2023/HopBox/")

# read in emmeans from linear models
files <-  list.files("./Output", pattern = c("emmean"), full.names = T)

# set output
output <- list()

for (i in 1:length(files)) {
  dat <- read.csv(files[i]) %>% dplyr::select(1:4) 
  colnames(dat)[4] <- "value"
  dat$trait <- gsub(".csv", "", str_split(files[i],pattern = "_")[[1]][2])
  output[[i]] <- dat
}

output <- do.call("rbind", output) %>% 
  pivot_wider(names_from = trait, values_from = value) %>% 
  rename(weight = avg, 
         cone_density = cone)

colnames(output)
output <- output %>% dplyr::select(length, width, area, perimeter, openness, weight, dgci, cone_density)


# set function for correlation table 
# from: http://sthda.com/english/wiki/elegant-correlation-table-using-xtable-r-package#:~:text=To%20get%20the%20lower%20or%20the%20upper%20part,lower.tri%28x%2C%20diag%20%3D%20FALSE%29%20upper.tri%28x%2C%20diag%20%3D%20FALSE%29
cor.test(x = output$dgci, y = output$area)

corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  # compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  # define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "***", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  # trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  # build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  # remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  # remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 

corstars(output)
write.csv(corstars(output), "./Output/correlation_table.csv")
