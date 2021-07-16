########################################################### 
            #### FUNCTIONS - ANNOTATION  ####
########################################################### 

# Annotation
# Add anotation to the read.geo file
add.annotation.bioconductor = function(x, anno.file) {
  set.l = x 
  geo.set = set.l$geo.set
  gse = set.l$gse
  colnames(geo.set) = sub(" ","_",as.vector(gse[[1]]$source_name_ch1))
  geo.set = cbind('PROBEID'=rownames(geo.set), geo.set)
  anno.file = anno.file$bioc_package
  
  if(!is.na(anno.file)) {
    n.install = paste(anno.file,'db',sep='.')
    if(n.install == 'IlluminaHumanMethylation450k.db') {n.install = 'FDb.InfiniumMethylation.hg19'}
    print(n.install)
    if( !n.install %in% rownames(installed.packages())) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(n.install,ask=F, suppressUpdates= T)
      library(n.install, character.only=T)
      if(n.install == 'FDb.InfiniumMethylation.hg19') {  
        z <- get450k()
        set.an = data.frame('PROBEID'=rownames(getNearestGene(z)),"SYMBOL"=getNearestGene(z)$nearestGeneSymbol)
        m.set.an = merge(geo.set, set.an)
        return(m.set.an)
      } else {
        set.an = select(get(n.install), keys=rownames(geo.set), columns=c("PROBEID","SYMBOL","ENTREZID"),  keytype="PROBEID")
        m.set.an = merge(geo.set, set.an) 
        return(m.set.an)
      }
    } else if (n.install == 'IlluminaHumanMethylation27k.db') {
      return(NA)
    } else {
      library(n.install, character.only=T)
      if(n.install == 'FDb.InfiniumMethylation.hg19') {  
        z <- get450k()
        set.an = data.frame('PROBEID'=rownames(getNearestGene(z)),"SYMBOL"=getNearestGene(z)$nearestGeneSymbol)
        m.set.an = merge(geo.set, set.an)
        return(m.set.an)
      } else {
        set.an = select(get(n.install), keys=rownames(geo.set), columns=c("PROBEID","SYMBOL","ENTREZID"),  keytype="PROBEID")
        m.set.an = merge(geo.set, set.an) 
        return(m.set.an)
      }
    }	
  } else { return(NA) }
  
}

# merge exp/methylation per sample.  Uses Mean of the data
merge.gene.expMeth = function(gen.name, dt.val, isEntrez=T){
  if(isEntrez){
    gen.symbol = as.vector(human.annotation[match(gen.name, human.annotation$ENTREZID),2])
  } else {
    gen.symbol = gen.name
  }
  # Expression Values
  exp.v  = grep('EXP', names(dt.val))
  if(length(exp.v) > 0){
    exp.table = NULL #; exp.table = as.data.frame(exp.table)
    exp.dt = lapply(dt.val[exp.v], function(x) x[which(x$SYMBOL == gen.symbol),!colnames(x) %in% c('PROBEID','ENTREZID','SYMBOL')])
    count.v = 0 
    for ( i in exp.dt) {
      if(nrow(i) >= 1) {
        mean.val = unlist(lapply(colnames(i), function(x) mean(as.numeric(as.vector(i[,x])))))
        names(mean.val) = colnames(i)
        if(count.v == 3){ mean.val = log2(mean.val)}
        exp.table = append(exp.table, mean.val )  
        count.v =  count.v + 1
      }
    }} else { exp.table = NULL}
  
  # Methylation Values
  meth.v  = grep('METH', names(dt.val))
  if(length(meth.v) > 0) {
    meth.table = table(NA,NA,NA); meth.table = as.data.frame(meth.table)
    meth.table = NULL #; meth.table = as.data.frame(meth.table)
    meth.dt = lapply(dt.val[meth.v], function(x) x[which(x$SYMBOL == gen.symbol),!colnames(x) %in% c('PROBEID','ENTREZID','SYMBOL')])
    for ( i in meth.dt) {
      if(nrow(i) >= 1) {
        mean.val = unlist(lapply(colnames(i), function(x) mean(as.numeric(as.vector(i[,x])))))
        names(mean.val) = colnames(i)
        meth.table = append(meth.table, mean.val )  
      }
    }
  } else {meth.table = NULL}
  # Return Data 
  ret.vector = list('exp'=exp.table,'meth'=meth.table)
  return(ret.vector)
}

# Plots scatter plot of the merge gene
scatter.expMeth <- function(exp, meth, gen.name) {
    a = exp; print(a)
    b = meth; print(b )
    scatter.smooth(x=b,y=a, xlab = 'Meth', ylab = 'Exp', main=gen.name)
  }

# GEO annotation
getGEO.anno = function(geo.id) {
  # Load required function
  require(GEOquery)
  require(affy)
  library(GEOmetadb)
  library(DBI)
  require(RH2)
  
  # Get by id and Save in a data.frame annotation
  # Change getGPL = F if the firewall make it fails
  gse <- getGEO(geo.id, GSEMatrix = TRUE, getGPL = F)
  # print(length(gse))
  
  # Convert to table adding information
  exprs.data = gse[[1]]@assayData$exprs
  pdata.geo = pData(gse[[1]])
  colnames(exprs.data) = pdata.geo[colnames(exprs.data), ]$title
  
  # Function to check if it's normalized
  # print(as.vector(pData(gse[[1]])$data_processing[1]))
  
  # Save and Return the data
  res.data = list('geo.set' = exprs.data, 'gse'=gse)
  
  # Annotation
  # Query Database
  # sqlfile = paste(getwd(),'GEOmetadb.sqlite',sep='/')
  sqlfile = "/Users/carlos/0_Projects/1_Msc_Biostatistics/7_TFM/data/GEOmetadb.sqlite"
  con = dbConnect(RSQLite::SQLite(), sqlfile)
  plat.id = unique(res.data$gse[[1]]$platform_id)
  all.bioconductor = dbGetQuery(con,paste("select gpl,title,bioc_package from gpl where gpl=",paste(plat.id,"'",sep=''),sep="'"))
  #names(all.bioconductor) = plat.id
  dbDisconnect(con)
  
  an.data.list <- add.annotation.bioconductor(res.data, all.bioconductor)
  has.annotation <-  unlist(lapply(an.data.list, function(x)  length(x) > 1))
  filt.data.list <- an.data.list[has.annotation]
  
  return(filt.data.list)
}

# Comparison data
comp.datasets <- function(set.to.compare, controlPos, samplePos) {
  
  # Filter Data.sets 
  comp.set <- set.to.compare[,!colnames(set.to.compare) %in% c('PROBEID','SYMBOL','ENTREZID')]
  comp.set <- data.frame(lapply(comp.set, as.numeric), stringsAsFactors=FALSE)
  size.data <- ncol(comp.set)
  #control <- grep('mock',colnames(comp.set))
  #sample <- setdiff(1:size.data, control)
  type.d <- rep(NA, size.data)
  type.d[controlPos] <- 'C'; type.d[samplePos] <- 'S'
  d.set <- data.frame('data'=colnames(comp.set),'type'=type.d)[c(controlPos,samplePos),]
  
  # Using Limma to check results
  library(limma)
  design  <-  model.matrix(  ~ 0 + type  ,data=d.set)
  fit <- lmFit(as.matrix(comp.set), design)
  fit2 <- eBayes(fit)
  top <- topTable(fit2, coef=1, adjust="fdr", number=nrow(fit$coefficients), sort.by="none")
  res.top <- data.frame(top, set.to.compare)
  return(res.top)
}

# Function to create the table with all statistics per data
# Only works for the correlation Exp/Methylation
get.shape = function(pos.en, dt.list) {
  a  <-  merge.gene.expMeth(human.annotation$ENTREZID[pos.en], dt.list)
  print(pos.en)
  if(is.null(a$exp) || is.null(a$meth) || length(a$exp) <= 2 ||  length(a$meth) <= 2  ) {
    print(NULL)
  }else{
    return(shape.param(exp=a$exp, meth=a$meth))
  }
}

# Filt values
filt.fun = function(min,max,param,set) {
  return(which(set[,param] >= min & set[,param] <= max))
}

# Functions to be tested
shape.param <-  function(exp, meth) {
  # Libraries
  require(scagnostics)
  require(energy)
  require(entropy)
  require(infotheo)
  require(splines)
  
  # Clean data
  a = as.numeric(exp[!is.na(exp)])
  b = as.numeric(meth[!is.na(meth)])
  
  # Same size , if not random subsample
  if(length(a) < length(b)) {
    b = sample(b,length(a))
  } else if (length(a) > length(b)) {
    a = sample(a,length(b))
  } else { a = a; b = b  }
    
  # Statistics
    corr = cor(a,b,method = "spearman")
    #spline.v = spline(a,b)
    minf = mutinformation(discretize(a),discretize(b))
    dcorr = dcor(a,b)
    sc.vector = scagnostics(a,b)
    res.df = data.frame(corr,minf, dcorr, as.list(sc.vector))
    return(res.df)
}

# Diagnostics using diffrent algoritm . It will return a data.frame with row per gene conteingin all the infomration of thte 
# Auxiliar function wich reproduces the value of a perfect fit
lshape.diagnostics = function(len.check=100, st.x=0.1, st.y=1) {
  l.per.shape  = data.frame('a'=c(st.y:len.check),'b'=rep(st.x,len.check) )
  l.per.shape[c(round(len.check*.5):len.check),]$b <- seq(0,1,length.out = round(len.check*.5)+1)
  l.per.shape[c(round(len.check*.5):len.check),]$a <- rep(st.x,round(len.check*.5)+1)
  val.diag = shape.param(l.per.shape$a, l.per.shape$b)
  with(l.per.shape, plot(y=a, x=b, type='p'))
  with(l.per.shape, lines(spline(x=b,y=a), col='red'))
  return(val.diag)
}

# K-NN clustering takes a data.frame with pos.factor as the factor to classify
knn.algorithm <- function(dataset, sel.seed, k.num, pos.factor) {
  
  # Creamos un set de datos aleatorio
  set.seed(sel.seed)
  
  # Filter table
  dset <- as.data.frame(scale(dataset[,-pos.factor]))
  
  # Creamos los sets de datos aleatorios
  test.set <- sample(1:nrow(dset),round(nrow(dset) * 0.2) )
  train.set <- c(1:nrow(dset))[!c(1:nrow(dset)) %in% test.set]
  
  # Creamos un data.set de training y otro de test. 
  # Para ello dividiremos las muestras para diferentes aplicaciones
  dataset_train <- na.exclude(dset[train.set, ])
  dataset_test <- na.exclude(dset[test.set, ])
  dataset_train_labels <- dataset[train.set, pos.factor]
  dataset_test_labels <- dataset[test.set, pos.factor]
  
  # Aplicamos la función de knn()
  require(class)
  dataset_test_pred <- knn(train = dataset_train, 
                           test = dataset_test,cl = dataset_train_labels, k = k.num)
  
  # Evaluación
  require(gmodels)
  
  # Create the cross tabulation of predicted vs. actual
  dataset_eval <-CrossTable(x = dataset_test_labels, y = dataset_test_pred, prop.chisq = FALSE)
  
  # Devolvemos una lista con elalgoritmo y la evaluación
  return(list('Pred'=dataset_test_pred,'Eval'=dataset_eval))  
}

# Evaluation function of the K-NN clustering
eval.fun <- function(model.set, print.frec=F, pos.name) {
  # Si print.frec = T , mostramos la tabla de frecuencias
  if(print.frec) {
    (round(prop.table(table(model.set$Pred))*100, digits = 2))
  }
  else {
    # Cargamos el set de datos
    eval.set =  as.data.frame(model.set$Eval$t)
    
    # Estimar los parametros
    true.val  = eval.set[eval.set$x == eval.set$y,]
    neg.val  = eval.set[eval.set$x != eval.set$y,]
    TP = true.val[true.val$x==pos.name, 3]; FP=neg.val[neg.val$y!=pos.name, 3]
    TN = true.val[true.val$x!=pos.name, 3]; FN=neg.val[neg.val$y==pos.name, 3]
    
    # Sensibilidad y Especificidad  
    sen.dat = round(100 * (TP / (TP + FN)), digits = 1)
    spc.dat = round(100 * (TN / (TN + FP)), digits = 1)
    
    # Valor predictivo positivo y Valor predictivo negativo
    pos.pred = round(100 * (TP / (FP + TP)), digits = 1)
    neg.pred = round(100 * (TN / (FN + TN)), digits = 1)
    
    # Accuracy y Error
    acc = (TP + TN) / (TP +TN +FP+FN)
    err.class = 1 - acc
    
    return(data.frame(TP,FP, TN, FN, acc, err.class,sen.dat, spc.dat, pos.pred, neg.pred) )
  }
}

# Parameters Evaluation
param.eval <-  function(data.s = shape.tabl, convex=0.1 , skinny=c(0.1,0.5), stringy=0.6 ) {
  # Avoid missed data
  data.s = na.exclude(data.s)
  thres.h = which(data.s$Convex < convex & data.s$Stringy > stringy & data.s$Skinny > skinny[1] & data.s$Skinny < skinny[2] ) 
  
  # We add the factors to 
  data.s$L.shape <- rep('N', nrow(data.s))
  data.s$L.shape[thres.h] <- 'L'
  
#   if(method.class == 'knn') {
#     # Knn Algorithm Evaluation
#     library(kernlab)
#     knn.set = knn.algorithm(dataset = data.s, sel.seed = 123, k.num = 3, pos.factor = ncol(data.s))
# 
#     # Eval the model
#     ev.knn = eval.fun(knn.set,print.frec = F, pos.name = 'L' )
#   }
#   
#   if(method.class == 'svm') {
#     # svm algorithm Evaluation
#     print('test')
#   }
  
  # Return vector
  #ret.dat = list('eval'=ev.knn, 'set.vector'=data.s)
  #print(ev.knn)
  return(data.s)
  
}

# GO ontology enrichment. 
go.enrichment <-  function(dt.res, gen.ass, th.pval) {
  library(goseq)
  geneNames <-  unique(as.vector(human.annotation[as.vector(human.annotation$SYMBOL) %in% rownames(dt.res), 2]))
  int.gene <- rep(0,length(geneNames))
  names(int.gene) <- geneNames
  int.gene[dt.res[rownames(dt.res) %in% geneNames, ]$L.shape == 'L'] <- 1
  pwf=nullp(int.gene, 'hg19','geneSymbol')
  GO.wall=goseq(pwf,gen.ass,"geneSymbol")
  GO.wall = GO.wall[GO.wall$ontology == 'BP' & GO.wall$over_represented_pvalue <= th.pval,]
  return(GO.wall)
}

# Do it with caret
svn.classification = function(param.shape.set, cv=3) {
  # Library
  require(caret)
  # Do crossvalitaion
  train_control<- trainControl(method="cv", number=cv, returnResamp = "all")
  # train the model 
  model<- train( param.shape.set[,-ncol(param.shape.set)], 
                 as.factor(param.shape.set[,ncol(param.shape.set)]), 
                 trControl = train_control, method="svmLinear",trace=FALSE)
  # make predictions
  predictions<- predict(model,param.shape.set[,-ncol(param.shape.set)])
  
  # Append predictions
  mydat<- cbind(param.shape.set, predictions)
  
  # summarize results
  (confusionMatrix<- confusionMatrix( mydat$predictions, param.shape.set$L.shape, positive='L' ))
  
  # List Results 
  res.list = list('conf'=confusionMatrix,'data'=mydat)
  
  # Return
  return(res.list)
  
  
}

# Interpretation
interp.classification <- function(pred.value, list.val.gene) {
  # This is going to handle all the information and print in the screen
  val.pred = pred.value[which(pred.value$L.shape == 'L' & pred.value$predictions == 'L'), ] 
  cor.gene = prop.table(table( list.val.gene$SYMBOL %in% rownames(val.pred)  )) * 100
  #top.go = go.enrichment(val.pred, 'hg19', 0.05) 
  return(list('val'=val.pred, 'cor'=cor.gene))
}

# Aumentar la cv no mejora el SVM algoritm. Otra medida para averiguar este correlación será si son o no negativos en correlación
pred.based.Cor <-  function(param.eval.set) {
  
  # Support vector machine with the predicción
  a = svn.classification(param.eval.set, 2)
  
  # Medida para observar si mejoramos la predicción de valores 
  cor.values = sapply(seq(-1,1,0.1), function(x) length(which(param.eval.set$cor.sp < x )))
  cor.Lshape = sapply(seq(-1,1,0.1), function(x) length(which(param.eval.set$cor.sp < x & param.eval.set$L.shape == 'L' )))
  predCor.Lshape = sapply(seq(-1,1,0.1), function(x) length(which(a$data$cor.sp < x & a$data$predictions == 'L' )))
  propCor.Lshape = sapply(seq(-1,1,0.1), function(x) (length(which(param.eval.set$cor.sp < x &param.eval.set$L.shape == 'L' ) ) /sum(shape.tabl$cor.sp < x )) * 100 )
  names(cor.values) = paste('cor',seq(-1,1,0.1),sep=''); names(cor.Lshape) = paste('cor',seq(-1,1,0.1),sep=''); names(propCor.Lshape) = paste('cor',seq(-1,1,0.1),sep='')
  cor.dataset = data.frame(cor.values, cor.Lshape, propCor.Lshape)
  
  # Show percentage of findings
  print(propCor.Lshape)
  
  # Saturation plot, we convert to log2 scale so we can measuse the point of saturaion
  rng.stringy = paste( "Stringy",
                       paste(format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Stringy)[1],digits = 2),
                             format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Stringy)[2],digits = 2),sep='_'
                       ),sep=':')
  rng.Covex = paste("Convex",
                    paste(format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Convex)[1], digits = 2),
                          format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Convex)[2],digits = 2),sep='_'
                    ),sep=':')
  
  rng.Skinny = paste("Skinny",
                     paste(format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Skinny)[1],digits = 2),
                           format(range(a$data[which(a$data$L.shape == 'L', a$data$predictions == 'L' ),]$Skinny)[2],digits = 2), sep='_'
                     ),sep=':')
  
  name.plot = paste(rng.Covex, rng.Skinny, rng.Skinny, sep = "; ")
  plot(seq(-1,1,0.1), log2(cor.values), type='b', col='darkblue', sub=name.plot); abline(v=seq(-1,1,0.1)[which(log2(cor.values) > 0)][1], col='darkblue')
  lines(seq(-1,1,0.1), log2(cor.Lshape), type='b', col='darkred'); abline(v=seq(-1,1,0.1)[which(log2(cor.Lshape) > 0)][1], col='darkred')
  lines(seq(-1,1,0.1), log2(predCor.Lshape), type='b', col='darkgreen'); abline(v=seq(-1,1,0.1)[which(log2(predCor.Lshape) > 0)][1], col='darkgreen')
  
  # Return
  return(a)
}
