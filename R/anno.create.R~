
####################################################
# ---- creation of IRangesList from snp index ---- #
####################################################

make.irl<-function(tsindx)
{
  kk=c(1,diff(tsindx))
  kk.start=c(1,which(kk!=1))
  kk.end<-c((kk.start[-1]-1),length(kk))
  return(IRangesList(start=tsindx[kk.start],end=tsindx[kk.end]))

}


##############################################
# ------ creation of nwgene.gr ------------- #
##############################################

## Function to add/subtract upstream span and downstream span to/from genomic postion based on strand 
## if strand is  '+' then upstream increased by d.span*1000 and downstream decreased by u.span*1000
## if strand is  '-' then upstream increased by u.span*1000 and downstream decreased by d.span*1000
span<-function(ip.grl,u.span,d.span)  
{
  sgl<- start(ip.grl)
  st=strand(ip.grl)
  stdf=as.data.frame(st)
  stt<-tapply(stdf$value,stdf$group_name, function(x) as.character(x[1]))
  str.span1=stt
  str.span1[str.span1=="+"]<-u.span*1000
  str.span1[str.span1=="-"]<-d.span*1000
  str.span1=as.numeric(str.span1)
  str.sgl=IntegerList((0/sgl)+str.span1)

  end.span1<-stt
  end.span1[end.span1=="+"]<-d.span*1000
  end.span1[end.span1=="-"]<-u.span*1000
  end.span1=as.numeric(end.span1)
  end.sgl=IntegerList((0/sgl)+end.span1)

  suppressWarnings(start(ip.grl)<-start(ip.grl)-str.sgl)
  suppressWarnings(end(ip.grl)<-end(ip.grl)+end.sgl)
  
  return(ip.grl)
}

create.nwg<-function(ipsnp,pathlist,u.span,d.span)
{
        
  ##### making snp.gr #########
  #############################
  cat(paste0("Obtaining base pair position..",'\r'))    
  orig.pos<- 1:length(ipsnp)

  ipsnp2<-ipsnp
  if(class(pathlist)=="matrix")
  {
    ipsnp2<-intersect(ipsnp,rownames(pathlist))
  }  

  locs.ob<- SNPlocs.Hsapiens.dbSNP142.GRCh37
  snp.gpos=snpsById(locs.ob, ipsnp2, ifnotfound="drop")
  df=data.frame(mcols(snp.gpos))
  os<-setdiff(ipsnp,df[,1])
  mis.pos<-match(os,ipsnp)

  newloc<-start(snp.gpos) #pos(snp.gpos)
  dd=seqnames(snp.gpos)
  newrle<-Rle(sub("ch","chr",rep(runValue(dd),runLength(dd))))

  if(length(mis.pos))
  {
    orig.pos<-orig.pos[-mis.pos]
  }
  cat(paste0("Obtaining base pair position.. 25%",'\r'))  
      
  if(class(pathlist)!="matrix")
  {

    ###### making gn.grl and snp.gr #########
    ##########################################
    y <- org.Hs.egSYMBOL
    mapped_genes <- mappedkeys(y)
    yy <- as.list(y[mapped_genes])

    txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
    txdb <- keepSeqlevels(txdb, seqlevels(txdb)[1:25])
    txbygene = transcriptsBy(txdb, "gene")
    comgene<-intersect(names(txbygene), names(yy))
  
    gn.grl<-txbygene[comgene]
    names(gn.grl)<-yy[comgene]

    gn.grl<-span(gn.grl,u.span,d.span)
    #snp.gr=as(snp.gpos,"GRanges")

    snp.gr<-GRanges(seqnames=newrle,IRanges(start=newloc, end=newloc))#strand=strand(snp.gpos)
 
    cat(paste0("Obtaining base pair position.. 50%",'\r'))

    ##### creation of nw.gene.gr with snp position ########
    ####################################################### 
 
    gn=names(gn.grl)

    res<-findOverlaps(gn.grl,snp.gr)
    resdf<-as.data.frame(res)
    acg<-unique(queryHits(res))
   
    ## Following steps are applied to identify genes with non-consecutive overlapped snps and making a proper granges
    ## example resdf[c(25380:25390),]
     
    all.pos<-rep(FALSE,nrow(resdf))
    dif.subhit<-c(1,diff(resdf$subjectHits))
    dif.qhit<-c(1,diff(resdf$queryHits))
    ch.pos<-which(dif.qhit!=0)
    dif.subhit[ch.pos]<-1
    ch.pos<-c(ch.pos,which(dif.subhit!=1))
    all.pos[ch.pos]<-TRUE
  
    cat(paste0("Obtaining base pair position.. 75%",'\r'))

    res.tp2=tapply(resdf$subjectHits,cumsum(all.pos), function(x) c(min(x), max(x)))#simplify=TRUE)
    utp2<-matrix(unlist(res.tp2), ncol=2,byrow=TRUE)
    all.chrseq<-rep("chr1",nrow(utp2))#length(all.end)
     
    nwgene.gr<-GRanges(seqnames=all.chrseq,IRanges(start=utp2[,1],end=utp2[,2]))

    
    all.gn<-tapply(resdf$queryHits,cumsum(all.pos), function(x) x[1])
    gn.rle<-Rle(all.gn)
    tte<-cumsum(runLength(gn.rle))
    tts<-c(1,(tte[1:(length(tte)-1)]+1))
    names(tte)<-gn[acg]
    names(tts)<-gn[acg]
    # both tts and tte are used to index nwgene.gr(eg:-nwgene.gr[tts["LOC100133331"]:tte["LOC100133331"],]) It gives snps of LOC100133331

    cat(paste0("Obtaining base pair position..completed",'\n'))
   
    return(list(orig.pos,nwgene.gr,tts,tte,newloc,newrle))
     
  }else
  {    
    cat(paste0('\r',"Obtaining base pair position..completed",'\n'))
     
    return(list(orig.pos,newloc,newrle,ipsnp[orig.pos]))
  } 
     
}

#####################################################
#------- creation of pathlist & path2snp map -------#
#####################################################

create.plist<-function(retv,pathlist,ipd.len)
{
  cat(paste0("Mapping Pathways To SNPs..",'\r'))
  
  orig.pos<-retv[[1]]
  used.sindx<-numeric(0)
  pl.id<-numeric(0)
  path2snp.indx<-IRangesList()
  plist<-GRangesList()

  if(class(pathlist)!="matrix")
  {  
    nwgene.gr=retv[[2]]
    tts=retv[[3]]
    tte=retv[[4]] 
    plen<-length(pathlist)
    gname1<-names(tte)   

  }else
  {
    plen<-ncol(pathlist)
    com.snp<-retv[[4]]
    indx<-match(com.snp,rownames(pathlist))   
  }

  j=1

  for( i in 1:plen) 
  {
      
      if(class(pathlist)!="matrix")
      {     
          tgn1<-intersect(gname1,pathlist[[i]]) 
    
          if(length(tgn1)!=0)
          {
             gindx<-unlist(as.vector(IRanges(start=tts[tgn1], end=tte[tgn1])))
           
             tg.gr1<-nwgene.gr[gindx]
             plist<-append(plist,GRangesList(tg.gr1))

             sindx<-orig.pos[sort(unique(unlist(as.vector(ranges(plist[[j]])))))]
             path2snp.indx<-append(path2snp.indx,make.irl(sindx))
             used.sindx<-c(used.sindx,sindx)
             pl.id<-c(pl.id,i)
             j=j+1
          }

      }else ## for annotation data
      {
         vv=intersect(which(pathlist[,i]!=0),indx)
         if(length(vv)!=0)
         {
           vv2= match(vv,indx)
           sindx<-orig.pos[vv2]
           ff=make.irl(sindx)
           gg=GRanges(seqnames=rep("chr1",length(ff[[1]])),ff[[1]])
           plist<-append(plist,GRangesList(gg))
           path2snp.indx<-append(path2snp.indx,ff)
           used.sindx<-c(used.sindx,sindx)
           pl.id<-c(pl.id,i)
         }

      }    
     
      if(i<plen)
      {
        cat(paste0('Mapping Pathways To SNPs.. ',floor((i/plen)*100),'%','\r'))
      }
     
  } ## end of for
  
         
  nopath.sindx<-(1:ipd.len)[-unique(used.sindx)]

  if(length(nopath.sindx)!=0)
  { 
    path2snp.indx<-append(path2snp.indx,make.irl(nopath.sindx))
  }else
  {
    path2snp.indx<-append(path2snp.indx,IRangesList(IRanges(NULL)))
  } 
   
  cat(paste0('\r',"Mapping Pathways To SNPs..completed",'\n'))
  
  return(list(plist,path2snp.indx,pl.id))

}

###################################################################
#------- Creation of Equivalence class of SNPs ----------#
###################################################################

create.pcomb<-function(retv,orig.pos,cls.path,pth.len)
{
  cat(paste0('Creating Equivalence Class..','\r'))
  plist=retv[[1]]
  path2snp.indx=retv[[2]]
  pl.id<-retv[[3]]
  plen<-length(retv[[2]])
  nopath.sindx=unlist(as.vector(unlist(retv[[2]][plen])))

  dj=disjoin(unlist(plist))
  ov.res<-findOverlaps(dj,plist)
   
  fresdf<-as.data.frame(ov.res)
  res2=tapply(fresdf$subjectHits,fresdf$queryHits,function(x) paste(x, collapse=","))
  rm(fresdf)

  restp2<-tapply(names(res2),res2, function(x) paste(x,collapse=","))
  pcomb=names(restp2)
  pcomb.len=length(pcomb)
  
  pcomb.mat<-matrix(0,nrow=pcomb.len+1,ncol=pth.len) 
  snp.eq<-numeric(0)

  for( j in 1:pcomb.len)
  {
       pindx<-pl.id[unique(as.integer(unlist(strsplit(pcomb[j],split=","))))]
       dj.indx<-as.integer(unlist(strsplit(restp2[pcomb[j]],split=",")))
       dj.rng<-ranges(dj[dj.indx])

       if(cls.path!="matrix")
       {
         sindx2=sort(orig.pos[unlist(as.vector(dj.rng))])
       }else
       {
         sindx2=sort(unlist(as.vector(dj.rng)))
       }

       snp.eq[sindx2]<-j
       pcomb.mat[j,pindx]<-1

       cat(paste0('Creating Equivalence Class.. ',floor((j/pcomb.len)*100),'%','\r'))
        
  }
 
  if(length(nopath.sindx))
  {
    snp.eq[nopath.sindx]<-nrow(pcomb.mat)
  }else
  { 
    pcomb.mat<-pcomb.mat[-nrow(pcomb.mat),]
  } 
  
  cat(paste0('\r',"Creating Equivalence Class..completed",'\n'))
  
  return(list(snp.eq,pcomb.mat))
}


###########################################
##-- Creattion of  final input pathway --##
########################################### 

build.path<-function(anno.mat,path.def)
{
  if(is.null(anno.mat) & is.null(path.def))
  {
    message("Provide either anno.mat or path.def!!!!")
    pthl<-FALSE
  
  }else if(!is.null(anno.mat) & !is.null(path.def))
  {
    message("Create object seperately,then use anno.merge to get composite object.")
    pthl<-FALSE

  }else if(!is.null(anno.mat))
  {
    ## removing columns having non-binary values
    nc=ncol(anno.mat)
    tt2=matrix((anno.mat %in% c(0,1)),ncol=nc)
    col.remove<-which(colSums(!tt2)>0)
    pthl<-anno.mat
    if(length(col.remove))
    {
      pthl<-anno.mat[,-col.remove]
    }
 
  }else
  {
    if(class(path.def)=="matrix")
    {
      pthl<-vector("list",ncol(path.def))
      gn.rn<-rownames(path.def)

      for(i in 1: ncol(path.def))
      {
        pthl[[i]]<-gn.rn[which(path.def[,i]==1)]
      }
    
    }else
    {
      pthl<-path.def
    }
  }
  
  return(pthl)
}


#################################
##-- Mapping object creation --##
#################################
 
anno.create <-function(snp.ids, anno.mat=NULL, path.def=NULL,u.span=250,d.span=100)
{   
   pth.len<-numeric(0)
   
   pathlist<-build.path(anno.mat,path.def)
   
   if(!is.logical(pathlist))
   {
     ## Obtaining base pair position..
     ret1<-create.nwg(snp.ids,pathlist,u.span,d.span)

     if(class(pathlist)=="matrix")
     {
       newloc<-ret1[[2]]
       newrle<-ret1[[3]]
       pth.len<-ncol(pathlist)

     }else
     {
       newloc<-ret1[[5]]
       newrle<-ret1[[6]]
       pth.len<-length(pathlist)
     }    
  
     ## Mapping Pathways To SNPs..
     ret2<-create.plist(ret1,pathlist,length(snp.ids))
  
     ## Creating Equivalence Class..
     ret3<-create.pcomb(ret2,ret1[[1]],class(pathlist),pth.len)
    
     chrno<-rep(NA,length(snp.ids)); pos<-rep(0,length(snp.ids))
     chrno[ret1[[1]]]<-as.numeric(sub("chr","",rep(runValue(newrle),runLength(newrle))))
     pos[ret1[[1]]]<-newloc
     pos[-ret1[[1]]]<-NA
  
     op.df<-data.frame(chrno=chrno,chrpos=pos,stringsAsFactors=FALSE)
     rownames(op.df)<-snp.ids
   
     vv<-c(nrow(op.df),nrow(ret3[[2]]),dim(ret3[[2]])[2],1)
     names(vv)<-c("no_snps","no_eq_class","no_anno","used_ob")
     dm.ls<-list(NULL,NULL,NULL)
     dm.ls[[1]]<-vv
     names(dm.ls)<-c("ob_dim","eq_class_arg","no_anno_arg")
  
     cat(paste0("Created Mapping Object Sucessfully",'\n'))
     
     rr=annotatedSNPset(snp.df=op.df,snp.eq=ret3[[1]],eq.mat=ret3[2],eq.map=NULL,dim=dm.ls)
     return(rr)

   }

}


  

    









