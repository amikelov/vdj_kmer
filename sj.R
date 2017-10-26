library(pbapply)
library(data.table)
library(ggplot2)
library(parallel)
library(plyr)

mitcr_clean<-function(input_file,out_file){
  SL <- fread(input_file,skip =1,head =T)
  SL <- SL[!grepl("~",SL$`CDR3 amino acid sequence`),]
  SL <- SL[!grepl("*",SL$`CDR3 amino acid sequence`,fixed = T),]
  SL$Percentage <- SL$`Read count`/sum(SL$`Read count`)
  SL<-data.frame(rep(gsub(x = gsub(x=input_file,pattern = "^.*\\/",""),pattern= "_t2.ec2.txt",replacement = ""),nrow(SL)),SL$`CDR3 amino acid sequence`,SL$`V segments`,SL$Percentage)
  fwrite(SL, out_file,sep = "\t",col.names = F)
}


#####################km_list

pair_ints<-function(samples,metafile){
  km_list<-lapply(samples, function(x){
    df<-fread(x,skip = 1,col.names = c("kmer","v","value"))
    df$kmer_v <- paste(df$kmer,df$v,sep = " ")
    df$p<- df$value/sum(df$value)
  
    return(df)
  })
  

  
  names(km_list)<-sapply(strsplit(samples,"/"),tail,1)
  km_list_len<-sapply(samples, function(x){scan(x,n = 1)})
  km_list_kmlen<-sapply(km_list,nrow)
  names(km_list_kmlen)<-names(km_list_len)
  nms <- combn( names(km_list) , 2 , FUN = paste , collapse = " " , simplify = FALSE )



  if (!missing(metafile) ){
    meta<- fread(metafile,header = T,sep2 = ";")
    n_feat<- ncol(meta)
    sample_pairs<-combn(samples,2,FUN =c,simplify=F)
    meta_pairs <- ldply(sample_pairs, function(x){
      df<- cbind(meta[meta$file==x[1],] , meta[meta$file==x[2],])
      names(df) <- c(paste0(names(meta),".x") ,paste0(names(meta),".y"))
      return(df)
    })
    
    for (i in 2:n_feat) {
      meta_pairs[[names(meta)[i]]]<- as.factor(mapply(FUN = function(x,y){paste( sort(c(x,y)),collapse="_")},
                                            as.character(meta_pairs[,i]),as.character(meta_pairs[,i+n_feat])))
    }
  }



  nms <- combn( names(km_list) , 2 , FUN = c , simplify = FALSE )
  ll <- combn( km_list , 2 , simplify = FALSE )


  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  clusterExport(cl,envir = environment(),c("km_list","nms"))
  
  meta_pairs$sj_dist <- parSapply(cl=cl, nms , function(x) {
    df <- merge(km_list[[x[1]]] , km_list[[x[2]]], by = "kmer_v" ) 
    m <- 0.5 * (df$p.x + df$p.y)
    JS <- 0.5 * (sum(df$p.x * log2(df$p.x / m)) + sum(df$p.y * log2(df$p.y / m)))
    return(sqrt(JS))
  })
  stopCluster(cl)
  
  rownames(meta_pairs)<-nms
  meta_pairs$name1<-sapply(strsplit(meta_pairs$file.x,"/"),tail,1)
  meta_pairs$name2<-sapply(strsplit(meta_pairs$file.y,"/"),tail,1)
  meta_pairs$len1<-sapply(meta_pairs$file.x,function(x){km_list_len[x]})
  meta_pairs$len2<-sapply(meta_pairs$file.y,function(x){km_list_len[x]})
  meta_pairs$kmlen1<-sapply(meta_pairs$file.x,function(x){km_list_kmlen[x]})
  meta_pairs$kmlen2<-sapply(meta_pairs$file.y,function(x){km_list_kmlen[x]})
  
  meta_pairs$len_product<-meta_pairs$len1*meta_pairs$len2
  
  meta_pairs$kmlen_product<-as.double(meta_pairs$kmlen1)*as.double(meta_pairs$kmlen2)
  return(res)

}

samples<- system("realpath ~/km_list/kmers/counts/*",intern = T)
km_list_new<-data.frame(sj_dist=pair_ints(samples, sj=T))
km_list_new$name1 <- sapply(strsplit(rownames(km_list_new)," "),"[[",1)
km_list_new$name2 <- sapply(strsplit(rownames(km_list_new)," "),"[[",2)
htable <- read.delim("~/km_list/km_list_TCR1.csv",header=F)
htable$name1 <- sapply(strsplit(as.character(htable$V1)," "),"[[",1)
htable$name2 <- sapply(strsplit(as.character(htable$V1)," "),"[[",2)
names(htable)[2:3] <- c("km_list", "cond")
htable <- htable[,2:5]

names(htable)[3:4] <- c("name2","name1")
km_list_search <- merge(km_list_new,htable)
names(htable)[3:4] <- c("name1","name2")
km_list_search <- rbind(km_list_search, merge(km_list_new,htable))
km_list_search <- unique(km_list_search)

ggplot(data=km_list_search,aes(km_list,sj_dist)) +geom_boxplot()
km_list_search$cond <- as.factor(km_list_search$cond)

ggplot(data=km_list_search,aes(cond,sj_dist)) +geom_boxplot()

ggplot(data=km_list_search,aes(cond,sj_dist)) +geom_boxplot()+geom_jitter()
ggplot(data=km_list_search,aes(km_list,sj_dist)) +geom_boxplot()+geom_jitter()
km_list_withoutLeb<- km_list_search[(km_list_search$name1!= "DL_8_km.txt" | km_list_search$name1!= "DL_8_km.txt"),]

ggplot(data=km_list_search,aes(km_list,sj_dist)) +geom_boxplot()+geom_point()+geom_text(aes(label=name1),hjust=1.3, vjust=0,size=2)+geom_text(aes(label=name2),hjust=-0.5, vjust=0,size=2)
ggplot(data=km_list_search,aes(km_list,sj_dist,fill=cond)) +geom_boxplot()+geom_jitter()+facet_grid(. ~ cond )

################total

# continious
samples<- system("realpath ~/asps/AS/kmers/count/*",intern = T)
ids<-sapply(strsplit(samples,"/"),tail,1)
meta<- data.frame( file=samples,id=ids,cond=c("as","as","h","as","as","h","as","as","as","as","h","h","h","as","as","as","as","as","as","as","as","as","as","as","as","as","as","as","as"),
                   subpop=c("cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd8","cd4","cd4","cd4","cd4","cd8","cd4","cd4","cd4","cd4","cd8","cd4","cd4","cd8","cd4","cd4","cd8"),
                   comp=c(rep("pb",16),"sf","sf","pb","pb","pb","sf","sf","pb","sf","sf","pb","sf","sf"))



df<-data.frame(sj_dist=pair_ints(samples))
df$name1 <- sapply(strsplit(rownames(df)," "),"[[",1)
df$name2 <- sapply(strsplit(rownames(df)," "),"[[",2)

df <- merge(df,meta,by.x='name1',by.y='id')
names(df)[5:7] <- c("cond1","subpop1","comp1")
df <- merge(df,meta,by.x='name2',by.y='id')
names(df)[9:11] <- c("cond2","subpop2","comp2")

df$cond1<- as.character(df$cond1)
df$cond2<- as.character(df$cond2)

df$subpop1<-as.character(df$subpop1)
df$subpop2<-as.character(df$subpop2)

df$comp1<-as.character(df$comp1)
df$comp2<-as.character(df$comp2)


df$cond <- paste(df$cond1,df$cond2,sep = "_")
df[df$cond=="as_h","cond"]<- "h_as"
df$cond<-as.factor(df$cond)
df$subpop <- paste(df$subpop1,df$subpop2,sep = "_")
df[df$subpop=="cd8_cd4","subpop"]<- "cd4_cd8"
df$subpop<-as.factor(df$subpop)
df$comp <- paste(df$comp1,df$comp2,sep = "_")
df[df$comp=="pb_sf","comp"]<- "sf_pb"
df$comp<-as.factor(df$comp)

donors <- data.frame(samp =unique(c(df$name2,df$name1)),
                     donor= c("ASH110","ASH11","Kal","Kal","Kal","Dv","Dv","Dv","Evst","ASH111","Luk","Mikh","Mikh","Mikh","Shep","Shep","Shep","DL","Dv","Evst","IrZv","Kal","Luk","Mikh","Shep","TA","TV","YB","ASH110"))

df <- merge(df,donors, by.x="name2",by.y="samp")
names(df)[7]="donor2"
df <- merge(df,donors, by.x="name1",by.y="samp")
names(df)[8]="donor1"




ggplot(df[df$comp1!="sf" & df$comp2!="sf", ],aes(subpop,sj_dist)) +geom_boxplot()
ggplot(df,aes(comp,sj_dist)) +geom_boxplot()
ggplot(df[df$subpop=="cd8_cd8" & df$comp =="pb_pb", ],aes(cond,sj_dist)) +geom_boxplot() +geom_text(aes(label=name1),hjust=-0.6,size=3)+geom_text(aes(label=name2),hjust=0.6,size=3)
