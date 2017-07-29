
setClassUnion("matrixORNULL", c("matrix", "NULL"))

## class- defination
setClass("annotatedSNPset",
representation= representation(snp.df="data.frame",snp.eq="numeric",eq.mat="list",eq.map="matrixORNULL",dim="list")
#prototype=list(snp.wpval=data.frame(),path.df=matrix())

)

### constructor ######
annotatedSNPset<-function(snp.df,snp.eq,eq.mat,eq.map,dim) new("annotatedSNPset",snp.df=snp.df,snp.eq=snp.eq,eq.mat=eq.mat,eq.map=eq.map,dim=dim)


#### show method ####
annotatedSNPset.print <- function(arg){

 df=slot(arg,"snp.df")
 cat(paste("$snp.df",".....",nrow(df) ,"X", ncol(df),"data.frame"),"\n")

 print(head(df,3))
 cat("..................\n")
 cat("..................\n")
 print(tail(df,3))
 
 cat("\n======================================================\n")
 
 snp.eq= slot(arg,"snp.eq")
 cat(paste("Equivalence class of SNPs of length",length(snp.eq)),"\n")
 print(head(snp.eq,20))
 
 cat("\n======================================================\n")

 eq.mat= slot(arg,"eq.mat")
 if(length(eq.mat)==1)
 {
   cat(paste("Equivalence class By Path matrix of dimension",nrow(eq.mat[[1]]),"X",ncol(eq.mat[[1]])),"\n")
   if(ncol(eq.mat[[1]])>5)
   { 
     print(eq.mat[[1]][1:4,1:4])

   }else
   {
     print(eq.mat[[1]][1:4,1:ncol(eq.mat[[1]])])
   }

 }else
 {
   cat(paste("List of ",length(eq.mat),"Equivalence class By Path matrix"),"\n")
   cat(paste("First matrix is .. ",nrow(eq.mat[[1]]),"X",ncol(eq.mat[[1]])),"\n") 
   #print(eq.mat[[1]][1:2,1:4])
   if(ncol(eq.mat[[1]])>5)
   { 
     print(eq.mat[[1]][c(1:4),1:4])

   }else
   {
     print(eq.mat[[1]][1:4,1:ncol(eq.mat[[1]])])
   }

   cat("..................\n")
   cat("..................\n")
   ln<-length(eq.mat)
   cat(paste("Last matrix is .. ",nrow(eq.mat[[ln]]),"X",ncol(eq.mat[[ln]])),"\n")
   if(ncol(eq.mat[[ln]])>5)
   { 
     print(eq.mat[[ln]][1:4,1:4])

   }else
   {
     print(eq.mat[[ln]][1:4,1:ncol(eq.mat[[ln]])])
   }
 
   #print(eq.mat[[ln]][1:2,1:4],)
 }
 cat("\n======================================================\n")
 eq.map= slot(arg,"eq.map")
 if(is.null(eq.map))
 {
   cat(paste("Equivalence class map matrix"),"\n")
   print(head(eq.map,4))

 }else
 {
   cat(paste("Equivalence class map matrix of dimension",nrow(eq.map),"X",ncol(eq.map)),"\n")
   print(head(eq.map,4))
 }

 cat("\n======================================================\n")
 dim2= slot(arg,"dim")
 cat(paste("Dimensions of various elments of mapping Object"),"\n")
 print(dim2)

}  

setMethod("show","annotatedSNPset", function(object){annotatedSNPset.print(object)})

### acessor method ######

fn.mat <- function(map.obj)
{
 n.maps <- (map.obj@dim)$ob_dim["used_ob"]
 n.eqc <- (map.obj@dim)$ob_dim["no_eq_class"]

 if(n.maps > 1) 
 {
  
  x.mat<-matrix(1,nrow=n.eqc)
  for(i in 1:n.maps)
  {
   x.mat<- cbind(x.mat, map.obj@eq.mat[[i]][map.obj@eq.map[,i],])
  }
 
  x.mat<-x.mat[,-1]

 }else 
 {
  x.mat <- map.obj@eq.mat[[1]]
 }

 x.mat
}

fn.eq<-function(map.obj)
{
  sn<-rownames(map.obj@snp.df)
  vv=map.obj@snp.eq
  names(vv)<-sn
  return(vv)
}

setGeneric("getEQ",function(x){standardGeneric("getEQ")})
setGeneric("getDF",function(x){standardGeneric("getDF")})
setGeneric("getDIM",function(x){standardGeneric("getDIM")})
setGeneric("getMAT",function(x){standardGeneric("getMAT")})
setMethod("getEQ","annotatedSNPset",function(x){fn.eq(x)})  #signature(x="annotatedSNPset")
setMethod("getDF","annotatedSNPset",function(x){x@snp.df})
setMethod("getDIM","annotatedSNPset",function(x){x@dim[[1]]})
setMethod("getMAT","annotatedSNPset",function(x){fn.mat(x)})

