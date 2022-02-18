########################
# Agro - 932 Homework 1#
# Author: Mason Lien   #
########################

#load Geno.txt file 
#file is stored in data folder in mutation_01
geno <- read.table(file.choose()) #choose the geno.txt file created from the hw1_arab.sh where the mutation rate is 0.01 
names(geno) <- c("chr", "pos", "ref", "alt", "l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", "l11", "l12", "l13", "l14", "l15", "l16", "l17", "l18", "l19", "l20")
head(geno)

#Calculate the Fst value for each site and visualize the results

for(i in 5:24){
  # replace slash and everything after it as nothing
  geno$newcol <- gsub("/.*", "", geno[,i] )
  # extract the line name
  nm <- names(geno)[i]
  # assign name for this allele
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a1")
  geno$newcol <- gsub(".*/", "", geno[,i] )
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a2")
}

#calculate Fst value for each site and visualize the results
#compute p1, p2, p with a 60/40 

geno$p <- apply(geno[, 25:64], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p <- geno$p/40
geno$p1 <- apply(geno[, 25:47], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p1 <- geno$p1/24
geno$p2 <- apply(geno[, 48:64], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p2 <- geno$p2/16


#theta function
pi <- function(n,p){
  return(n/(n-1)*(1-p^2-(1-p)^2))
}

#compute theta for population
for (i in seq_len(nrow(geno))){
  geno$theta_p1 <- pi(n=20, p=geno$p1)
  geno$theta_p2 <- pi(n=20, p=geno$p2)
}


pdf(file = "graphs/theta_wmutation01.PDF",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(geno$pos, geno$theta, xlab="Physical position", ylab="theta value", main="")
dev.off()

#calculate Fst
geno$fst <- with(geno, ((p1-p)^2 + (p2-p)^2)/(2*p*(1-p)))

#write Fst results
write.table(geno, file = "cache/fst_results_wmutation01.csv", sep = ",", row.names = F, quote = F)

#visualize Fst results with mutation rate of 0.01
fst <- read.csv("cache/fst_results_wmutation01.csv")

pdf(file = "graphs/Fst_wmutation01.PDF",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(fst$pos, fst$fst, xlab="Physical position", ylab="Fst value", main="")
dev.off()


###run same for mutation rate of 10%
geno <- read.table(file.choose()) #choose the geno.txt file created from the hw1_arab.sh where the mutation rate is 0.10
names(geno) <- c("chr", "pos", "ref", "alt", "l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8", "l9", "l10", "l11", "l12", "l13", "l14", "l15", "l16", "l17", "l18", "l19", "l20")
head(geno)

#Calculate the Fst value for each site and visualize the results

for(i in 5:24){
  # replace slash and everything after it as nothing
  geno$newcol <- gsub("/.*", "", geno[,i] )
  # extract the line name
  nm <- names(geno)[i]
  # assign name for this allele
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a1")
  geno$newcol <- gsub(".*/", "", geno[,i] )
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a2")
}


#calculate Fst value for each site and visualize the results
#compute p1, p2, p with a 60/40 

geno$p <- apply(geno[, 25:64], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p <- geno$p/40
geno$p1 <- apply(geno[, 25:47], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p1 <- geno$p1/24
geno$p2 <- apply(geno[, 48:64], 1, function(x) {sum(as.numeric(as.character(x)))})
geno$p2 <- geno$p2/16

#compute theta

pi <- function(n,p){
  return(n/(n-1)*(1-p^2-(1-p)^2))
}

#compute theta for population
for (i in seq_len(nrow(geno))){
  geno$theta_p1 <- pi(n=20, p=geno$p1)
  geno$theta_p2 <- pi(n=20, p=geno$p2)
}

#calculate Fst
geno$fst <- with(geno, ((p1-p)^2 + (p2-p)^2)/(2*p*(1-p)))

#write Fst results
write.table(geno, file = "cache/fst_results_wmutation10.csv", sep = ",", row.names = F, quote = F)

#visualize Fst results
fst <- read.csv("cache/fst_results_wmutation10.csv")

pdf(file = "graphs/Fst_wmutation10.PDF",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
plot(fst$pos, fst$fst, xlab="Physical position", ylab="Fst value", main="")
dev.off()




