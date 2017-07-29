
# Merging more than mapping objects into a composite object #
#############################################################

anno.merge<-function(anno.list)
{
  len=length(anno.list)

  if(len>1)
  { 
    all.eq.mat<-vector("list",len)
    all.eq.mat[[1]]<-slot(anno.list[[1]],"eq.mat")[[1]]
    vv=dim(all.eq.mat[[1]])
    vv.row<-vv[1]
    vv.col<-vv[2]

    all.snp.eq<-slot(anno.list[[1]],"snp.eq") 

    rn<-rownames(slot(anno.list[[1]],"snp.df"))
      
    for(j in 2:len)
    {
                                     #all.eq<-paste(all.eq,anno.list[[j]][[1]])
     rn2<-rownames(slot(anno.list[[j]],"snp.df"))

     if(identical(rn,rn2))
     {
       
        all.snp.eq<-paste(all.snp.eq,slot(anno.list[[j]],"snp.eq"))
        
        all.eq.mat[[j]]<-slot(anno.list[[j]],"eq.mat")[[1]]
        vv=dim(all.eq.mat[[j]])
        vv.row<-c(vv.row,vv[1])
        vv.col<-c(vv.col,vv[2])

     }else
     {
       message("Need Identical SNP sets of Mapping Object")
       stop()
     }

    }

    all.snp.eq<-as.factor(all.snp.eq)
    fn.eq<-as.numeric(all.snp.eq)
    tt=as.numeric(unlist(strsplit(levels(all.snp.eq)," ")))
  
    eq.map<-matrix(0,nrow=(length(tt)/len),ncol=len)
    for(i in 1:len)
    {
      eq.map[,i]<-tt[seq(i,length(tt),by=len)]
    }

   
    df<-slot(anno.list[[1]],"snp.df")
    dm.eqmap<-dim(eq.map)
    ob.dm<-c(nrow(df),nrow(eq.map),sum(vv.col),len)
    names(ob.dm)<-c("no_of_snps","no_eq_class","no_of_anno","used_ob")
    dm.ls<-list(NULL,NULL,NULL)
    dm.ls[[1]]<-ob.dm
    dm.ls[[2]]<-vv.row
    dm.ls[[3]]<-vv.col
    names(dm.ls)<-c("ob_dim","eq_class_of_arg","no_of_anno_of_arg")
   
    rr=annotatedSNPset(snp.df=df,snp.eq=fn.eq,eq.mat=all.eq.mat,eq.map=eq.map,dim=dm.ls)
    return(rr)
    #return(list(snp.eq=fn.eq,eq.mat=eq.mat))
  
  }else
  {
    message("Need Atleast two objects for merging")
    return(anno.list)
  } 
 
}
