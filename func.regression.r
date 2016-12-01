func.glmnet.1.cv = function (exFile, meFile)
  
{
  is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
  } 

  for(i in c("glmnet","methods")) {
  if (!is.installed(i)){
    install.packages(i)
    }
  }
  library(glmnet)
  library(methods)
  # get the expression and methylation files
  ex = read.delim(exFile)
  me = read.delim(meFile)
  
  #ex = ex[complete.cases(ex), ]  #omit NA rows
  ex=ex[,colSums(is.na(ex))==0]  #omit NA columns
  
  
  #length(which(rowSums(is.na(test))>0))
  #length(which(colSums(is.na(test))>0))
  
  #me = me[complete.cases(me), ]  #omit NA rows
  me=me[,colSums(is.na(me))==0]  #omit NA columns
  
  
  
  # get the ids
  ex.gene = ex[,1]
  ex.id = names(ex)[-1]
  me.site = me[,1]
  me.id = names(me)[-1]
  n.site = length(me.site)
  # match ids
  m.mch = match(me.id, ex.id)
  idx.mch = which(!is.na(m.mch))
  me.mch =m.mch[idx.mch]
  me.val1 = me[,-1]
  me.val2 =me.val1[,idx.mch]
  ex.val1= as.vector(ex[,-1])
  ex.val2 = ex.val1[me.mch]
  
  genevar=var(as.vector(ex.val2,"numeric"))
  
  
  
  ####standardization
  
 #colMeans(ex.val2)  # faster version of apply(scaled.dat, 2, mean)
 #rowMeans(ex.val2)   ###12.03
 #apply(ex.val2, 1, sd)   #0.32
 #b=scale(me.val2, center = T, scale = TRUE)
 tmp=as.data.frame(t(ex.val2))
 tmp=scale(tmp, center = T, scale = T)
 ex.val2= as.data.frame(t(tmp))
 row.names(ex.val2)="1"
 
 #rowMeans(ex.val2)    #-4.083497e-16
 #apply(ex.val2, 1, sd)   #1
  ############################################################################-------------------------------------------single
  ######single regression, extract beta, p.value, and R2.max, R2.max.var
  beta.single=rep(NA,n.site)
  p.single=rep(NA,n.site)
  R2.single=NULL
  R2.single.max.var=NULL
  R2.single.max=NULL
  
  
  for (i in (1:n.site)){
    lmfit = lm(t(ex.val2) ~ t(me.val2[i,]))
    info = summary(lmfit)
    beta.single[i+1] = info$coefficients[,1][2]
    p.single[i+1]=info$coefficients[,4][2]
    R2.single = c(R2.single,info$"r.squared")
  }
  
  
  R2.single.max=max(R2.single)
  R2.single.max.var=var(R2.single)
  
  
  #cross validation for single regression
  n.sample = length(ex.val2)  
  m = ceiling(n.sample/5)  
  R2.single.cv=NULL
  R2.single.cv.mn=NULL
  R2.single.cv.var=NULL
  R2.single.cv.max=NULL ###next time, delete "mn"
  R2.single.cv.max.var=NULL
  
  for (j in (1:n.site)){
    
    me.val3=me.val2[j,]
    
    for (i in (1:10)){
      idx.out = sample((1:n.sample), m, replace=FALSE)
      idx.in = setdiff((1:n.sample), idx.out)
      newdata = t(as.matrix( me.val3[,idx.out]) ) 
      fit = lm(t(ex.val2[idx.in]) ~ t(me.val3[,idx.in]))
      yhat = newdata %*%(as.vector(coefficients(fit)[-1]))+coefficients(fit)[1]
      resi = ex.val2[idx.out]-yhat
      R2.single.cv =c(R2.single.cv, 
                      (var(as.vector(ex.val2[idx.out], "numeric"), na.rm=TRUE)-
                         var(as.vector(resi, "numeric"), na.rm=TRUE))/
                        var(as.vector(ex.val2[idx.out], "numeric")))
    }
    R2.single.cv.mn = c(R2.single.cv.mn,mean(R2.single.cv))
    R2.single.cv.var = c(R2.single.cv.var,var(R2.single.cv))
    
  }
  
  R2.single.cv.max=max(R2.single.cv.mn)
  R2.single.cv.max.var=max(R2.single.cv.var)
  
  
  ####################################################################################---------------------------------------------multiple
  # multiple regression
  lmfit = lm(t(ex.val2) ~ t(me.val2))
  info = summary(lmfit)
  test = anova(lmfit)
  
  # output of observed data
  beta.multiple = info$coefficients[,1]
  p.multiple = info$coefficients[,4]
  
  p.multiple.overall = test$"Pr(>F)"[1]
  R2.multiple = info$"r.squared"
  R2.multiple.adjust = info$"adj.r.squared"
  sites = c("Intercept", as.character(me.site))
  
  ######## corss validation for multiple
  R2.cv=NULL
  R2.multiple.cv=NULL
  R2.multiple.cv.var=NULL
  
  n.sample = length(ex.val2)  
  m = ceiling(n.sample/5)   ##nFold=5, nTimes=10
  
  for (i in (1:10)){
    idx.out = sample((1:n.sample), m, replace=FALSE)
    idx.in = setdiff((1:n.sample), idx.out)
    newdata = t(as.matrix( me.val2[,idx.out]) ) 
    fit = lm(t(ex.val2[idx.in]) ~ t(me.val2[,idx.in]))
    yhat = newdata %*%(as.vector(coefficients(fit)[-1]))+coefficients(fit)[1]
    resi = ex.val2[idx.out]-yhat
    R2.cv =c(R2.cv, 
             (var(as.vector(ex.val2[idx.out], "numeric"), na.rm=TRUE)-
                var(as.vector(resi, "numeric"), na.rm=TRUE))/
               var(as.vector(ex.val2[idx.out], "numeric")))
  }
  
  # get the mean and variance of R2.cv
  R2.multiple.cv = mean(R2.cv)
  R2.multiple.cv.var = var(R2.cv) 
  
  
  ###########################################################################-------------------------------------------------------glmnet
 R2.glmcv=NULL
  R2.glmnet.cv=NULL
  R2.glmnet.cv.var=NULL
 
 R2.glmnet=NULL

  beta.glmnet=NULL
  
  n.sample = length(ex.val2)  
  m = ceiling(n.sample/5)   ##nFold=5, nTimes=10
 # if (nrow(me.val2)>=2){
    for (i in (1:10)){
      idx.out = sample((1:n.sample), m, replace=FALSE)
      idx.in = setdiff((1:n.sample), idx.out)
      newdata = t(as.matrix( me.val2[,idx.out]) ) 
      
      x=t(me.val2[,idx.in])
      #x[is.na(x)]=0
      y=t(ex.val2[,idx.in])
      #y[is.na(y)]=0
      
      cvfit = cv.glmnet(x, y)  #cvfit here, a list with all the ingredients of the cross-validation fit
      cvfit$lambda.min  # lambda.min is the value of λ that gives minimum mean cross-validated error
      coef(cvfit, s = "lambda.min")
      a=predict(cvfit, newx = t(me.val2[,idx.out]), s = "lambda.min")
      resi=t(ex.val2[idx.out])-a
      R2.glmcv =c(R2.glmcv, 
                  (var(as.vector(ex.val2[idx.out], "numeric"), na.rm=TRUE)-
                     var(as.vector(resi, "numeric"), na.rm=TRUE))/
                    var(as.vector(ex.val2[idx.out], "numeric")))
      
    }
    
    R2.glmnet.cv = mean(R2.glmcv)
    R2.glmnet.cv.var = var(R2.glmcv) 
    
    ####beta for glmnet
    
    
    x=t(me.val2)
    y=t(ex.val2)
    # a=scale(x, center = TRUE, scale = TRUE)
    #  b=scale(y, center = TRUE, scale = TRUE)
    
    cvfit = cv.glmnet(x, y,lambda=NULL,nfolds=10,alpha=1,standardize=T,standardize.response=T) 
    #cvfit here, a list with all the ingredients of the cross-validation fit
    cv.lambda=cvfit$lambda.min
    beta.glmnet=as.data.frame(as.matrix(coef(cvfit,s=cv.lambda)))
 a=predict(cvfit, newx = t(me.val2), s = "lambda.min")
 resi=t(ex.val2)-a
 R2.glmnet = (var(as.vector(ex.val2, "numeric"), na.rm=TRUE)-
                var(as.vector(resi, "numeric"), na.rm=TRUE))/
   var(as.vector(ex.val2, "numeric"))

    
  #} else{
    #R2.glmcv=NA
    #R2.glmnet=NA
    #R2.glmnet.var=NA
   # beta.glmnet=NA
  #}s
 ##############################adding dist
 #folder=strsplit(exFile,"/")[[1]][1]
 folder="run"
 folderGene=strsplit(exFile,"/")[[1]][2]
 gene_anno=read.delim(paste(folder,folderGene,"gene_anno.txt",sep="/"),sep="\t")
 cpg_anno=read.delim(paste(folder,folderGene,"cpg_anno.txt",sep="/"),sep="\t")
 
 #####cpg probes are not consistent between annotation file and the file used to perform linear regression.
 if (length(gene_anno$txStart)!=0){
 

  cpg_profile=data.frame(IlmnID=me[,1],dist=rep(0.001,length(me[,1])),txStart=rep(gene_anno$txStart,length(me[,1])),txEnd=rep(gene_anno$txEnd,length(me[,1])))
  cpg_profile_anno=merge.data.frame(cpg_profile,cpg_anno,by="IlmnID",all.x=T,all.y=F)
  
  #> dim(cpg_profile_anno)
  #[1] 16  7
  for (i in 1:length(me[,1])){
    if (gene_anno$strand=="+"){
      cpg_profile_anno$dist[i]=cpg_profile_anno$MAPINFO[i]-cpg_profile_anno$txStart[i]
    } else{
      cpg_profile_anno$dist[i]=cpg_profile_anno$txEnd[i]-cpg_profile_anno$MAPINFO[i]
    }
    
  }
  
  cpg_profile_anno=cpg_profile_anno[,c(1,2)]
  
  } else {
    cpg_profile_anno=data.frame(IlmnID=me[,1],dist=rep(0.001,length(me[,1])))
  } 
  
  ###########################################################################################------------------------------------final output
  ## organize for output
  out = data.frame(sites,n.site,as.character(ex.gene),beta.single,beta.multiple,beta.glmnet,
                   R2.single.max,R2.single.max.var,R2.single.cv.max,R2.single.cv.max.var,
                   R2.multiple,R2.multiple.adjust,R2.multiple.cv,R2.multiple.cv.var,
                   R2.glmnet,R2.glmnet.cv,R2.glmnet.cv.var,
                   p.single,p.multiple,p.multiple.overall,genevar
                   
  )
  file.name = paste("run/",eval(ex.gene), ".txt",sep="")
  colname = c("CpG","n.site","gene","beta.single","beta.multiple","beta.glmnet",
              "R2.single.max","R2.single.var","R2.single.cv.max","R2.single.cv.max.var",
              "R2.multiple","R2.multiple.adjust","R2.multiple.cv","R2.multiple.cv.var",
              "R2.glmnet","R2.glmnet.cv","R2.glmnet.cv.var",
              "p.single","p.multiple","p.multiple.overall","genevar","dist"
              
  )
  
  out=merge.data.frame(out,cpg_profile_anno,by.x="sites",by.y="IlmnID",all.x=T,all.y=T)
  write.table(out, file=file.name, sep="\t", row.names=F, col.names=colname)
  
  obj=NULL
  obj$sites=sites
  obj$n.site=n.site
  obj$ex.gene=ex.gene
  obj$beta.single=beta.single
  obj$beta.multiple=beta.multiple
  obj$beta.glmnet=beta.glmnet
  obj$R2.single.max=R2.single.max
  obj$R2.single.max.var=R2.single.max.var
  obj$R2.single.cv.max=R2.single.cv.max
  obj$R2.single.cv.max.var=R2.single.cv.max.var
  obj$R2.multiple=R2.multiple
  obj$R2.multiple.adjust=R2.multiple.adjust
  obj$R2.multiple.cv=R2.multiple.cv
  obj$R2.multiple.cv.var=R2.multiple.cv.var
  obj$R2.glmnet=R2.glmnet
  obj$R2.glmnet.cv.var=R2.glmnet.cv.var
  obj$p.single=p.single
  obj$p.multiple=p.multiple
  obj$p.multiple.overall=p.multiple.overall
  obj$genevar=genevar
#  obj$dist=dist 
  obj$R2.glmnet.cv=R2.glmnet.cv

  return ( obj)
  
}

