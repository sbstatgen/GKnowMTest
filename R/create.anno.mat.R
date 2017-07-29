
create.anno.mat<-function(snp, base.path='annotation-', fl.suffix=".txt", ncols=454, start.anno=4, snp.col=3, header=TRUE,index=NULL)
{
 locs.ob<- SNPlocs.Hsapiens.dbSNP142.GRCh37
 snp.gpos=snpsById(locs.ob,snp,ifnotfound="drop")
 df=data.frame(snp.gpos)[,c(1,4)]
 colnames(df)<-c("chrno","rsid")
 df$chrno<-sub("ch","",df$chrno)
 
 nc=(ncols-start.anno)+1

 rn<-character(0)
 fmat<-matrix(NA,ncol=nc,nrow=1) #fmat is annotation matrix

 smat<-matrix(0,ncol=2,nrow=22)
 colnames(smat)<-c("Total SNPs","Matched SNPs")

 ab.snp<-vector("list",22)     # chromosome wise SNPs that are not present in annotation file.
 unm.snp=setdiff(snp,df$rsid)  # un-mapped SNPs. These SNPs are not in dbsnp142. 

 if(is.null(index)) ## Building Index
 {
   anno.snp.list<-vector("list",22)
   message("Building Index File for Input SNPs....")

   for (i in 1:22)
   {
     cat(paste0("chr",i,'\n'))
     tt=fread(paste0(base.path,"chr",i,fl.suffix),select=snp.col)
     tt=data.frame(tt)
     anno.snp.list[[i]]<-tt[,1]
   }
   
  message("Building Index completed...")
  cat("\n")

 }else
 {
   anno.snp.list<-index
 }

 if(header)
 {
   hf<-1
 }else
 {
   hf<-0
 }

 message("Creating Annotation matrix...")
 for(i in 1:22)
 {  
   totsnp<-na.omit(df[df$chrno==i,])$rsid
   cs<-intersect(totsnp,anno.snp.list[[i]]) 
   indx<-sort(match(cs,anno.snp.list[[i]])) 
   len.indx<-length(indx)
   
   smat[i,1]<-length(totsnp)
   smat[i,2]<-len.indx

   ab.snp[[i]]<-setdiff(totsnp,cs)
     
   if(len.indx!=0)
   {
      mat<-matrix(0,nrow=len.indx,ncol=nc)
      path<-paste0(base.path,"chr",i,fl.suffix)
      rn<-c(rn,anno.snp.list[[i]][indx])

      out<-.C("getrow",as.character(path),as.integer(indx),as.integer(c(dim(mat),start.anno,hf)),res=as.vector(as.integer(mat)))

      newmat <- matrix(out$res, dim(mat), byrow=TRUE)
      fmat<-rbind(fmat,newmat)
   }

   cat(paste0("\n Chromosome ",i," (Matched SNPs / Total SNPs): ",len.indx," / ",length(totsnp)))
   cat("\n*****")

 }

 fmat=fmat[-1,]
 rownames(fmat)<-rn

 cat("\n **** Summary ****")
 cat(paste0("\n Total SNPs: ",length(snp)))
 cat(paste0("\n Mapped SNPs: ",(length(snp)-length(unm.snp))))
 cat(paste0("\n No.of Un-Mapped SNPs:",length(unm.snp)))
 cat(paste0("\n Matched SNPs: ",length(rn)))
 cat("\n")

 return(list(anno.mat=fmat,snp.count=smat,ab.snp=ab.snp,unm.snp=unm.snp,index=anno.snp.list))

}

