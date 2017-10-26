#!/usr/bin/env Rscript
args =  commandArgs(trailingOnly=TRUE)
library(data.table)
library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library("data.table"))

parSapply(cl = cl,args,function(i){
  if (substr(i,nchar(i)-2,nchar(i))==".gz"){
    system(paste("gunzip -k -f",i))
    df <-fread(paste0("cut -f2,6,10,11 ",substr(i,1,nchar(i)-3)))[aminoAcid !=""]
    l=7
  } else {
    df <-fread(paste0("cut -f2,6,10,11 ",i))[aminoAcid !=""]
    l=4
  }
  df[vMaxResolved=="unresolved"]$vMaxResolved <-ifelse(test = df[vMaxResolved=="unresolved"]$vFamilyTies=="", 
							df[vMaxResolved=="unresolved"]$vGeneNameTies,
							df[vMaxResolved=="unresolved"]$vFamilyTies)

  df <- df[!grepl(pattern = "*", df$aminoAcid,fixed = T),]
  df[grep(pattern ="or",x=df$vMaxResolved)]$vMaxResolved <-sapply(strsplit(x =  grep(pattern ="or",x=df$vMaxResolved,value=TRUE), split = ","),"[[",1)
  df[grep(pattern ="or",x=df$vMaxResolved)]$vMaxResolved <- sapply(strsplit(x =  grep(pattern ="or",x=df$vMaxResolved,value=TRUE), split = "-or"),"[[",1)
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = "V0",replacement = "V")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = "-0",replacement = "-")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = '\\*..',replacement = "")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TCRB',replacement = "TRB")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV15-1',replacement = "TRBV15")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV28-1',replacement = "TRBV28")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV1-1',replacement = "TRBV1")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV2-1',replacement = "TRBV2")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV13-1',replacement = "TRBV13")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV14-1',replacement = "TRBV14")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV16-1',replacement = "TRBV16")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV18-1',replacement = "TRBV18")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV19-1',replacement = "TRBV19")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV20',replacement = "TRBV20-1")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV21',replacement = "TRBV21-1")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV23',replacement = "TRBV23-1")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV27-1',replacement = "TRBV27")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV25',replacement = "TRBV25-1")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV30-1',replacement = "TRBV30")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV9-1',replacement = "TRBV9")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBVA-or09_02',replacement = "TRBVA")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = 'TRBV17-1',replacement = "TRBV16")
  df$vMaxResolved <- gsub(x=df$vMaxResolved,pattern = '-1-1',replacement = "-1")

  samp <- strsplit(i,"/")[[1]][length(strsplit(i,"/")[[1]])]
  samp <- substr(samp,1,nchar(samp)-l)
  df$Sample = samp
  df <- df[,c("Sample","aminoAcid","vMaxResolved")]
  path <- paste0(paste0(strsplit(i,"/")[[1]][1:(length(strsplit(i,"/")[[1]])-1)], collapse = "/"),"/immunoseqTidy.txt")
  fwrite(x=df, file=path,quote = F,sep="\t",row.names = F,col.names = F,append = T)
  
  if (substr(i,nchar(i)-2,nchar(i))==".gz") system(paste("rm",substr(i,1,nchar(i)-3)))
  return(NULL)
}
)
stopCluster(cl)
