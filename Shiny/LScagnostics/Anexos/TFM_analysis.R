########################################################### 
            #### FUNCTIONS - ANNOTATION  ####
########################################################### 

# Read GEO ID and load information 
read.geo = function(geo.id) {
  # Load required function
  require(GEOquery)
  require(affy)

  # Get by id and Save in a data.frame annotation
  # Change getGPL = F if the firewall make it fails
  gse <- getGEO(geo.id, GSEMatrix = TRUE, getGPL = F)
  print(length(gse))

  # Convert to table adding information
  exprs.data = gse[[1]]@assayData$exprs
  pdata.geo = pData(gse[[1]])
  colnames(exprs.data) = pdata.geo[colnames(exprs.data), ]$title

  # Function to check if it's normalized
  print(as.vector(pData(gse[[1]])$data_processing[1]))

  # Save and Return the data
  res.data = list('geo.set' = exprs.data, 'gse'=gse)
  return(res.data)
}

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
scatter.expMeth <- function(merge.g) {
  if(is.null(merge.g$exp) || is.null(merge.g$meth)) {
    print('Exp or Meth is not in the data')
    } else {
  a = merge.g$exp[order(merge.g$exp)]
  b = merge.g$meth[order(merge.g$exp)]
  scatter.smooth(x=b,y=sample(a, length(b)), xlab = 'Meth', ylab = 'Exp', main=)
  #qplot(x=b,y=sample(a, length(b)), geom='smooth', span =0.5)
  #plot(sample(a, length(b)) ~ b, xlab='meth',ylab='exp') 
  #lines(y=sample(a, length(b)), x=b, col='red')
    }
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
    cor.sp = cor(a,b,method = "spearman")
    spline.v = spline(a,b)
    minf = mutinformation(discretize(a),discretize(b))
    dcor.v = dcor(a,b)
    sc.vector = scagnostics(a,b)
    sp.spline = scagnostics(spline.v$x, spline.v$y)
    res.df = data.frame(cor.sp,minf, dcor.v, as.list(sc.vector), as.list(sp.spline))
    names(res.df)[c(12:20)] <-  sub('1','sp', names(res.df)[c(12:20)])
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
param.eval <-  function(data.s = shape.tabl, convex=0.1 , skinny=c(0.2,0.6), stringy=0.6 ) {
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
  GO.wall=goseq(pwf,gen.ass,"geneSymbol",use_genes_without_cat = TRUE)
  GO.wall = GO.wall[GO.wall$ontology == 'BP' & GO.wall$over_represented_pvalue <= th.pval,]
  return(GO.wall)
}

########################################################### 
      #### DOWNLOAD GEO DATA AND ANNOTATE  ####
###########################################################
# 0. Load All the information 
setwd('/Users/carlos/0_Projects/1_Msc_Biostatistics/7_TFM/data/')
options(stringsAsFactors = F)

# 1. Load GEO from table geo id
data.sets =  read.csv('geo.id.colon-cancer.txt', header=T, sep="\t")
list.data.sets = lapply(1:nrow(data.sets), function(x) read.geo(data.sets$GEOid[x]))
names(list.data.sets) = as.vector(data.sets$GEOid)
#save(list.data.sets, file="geo.id.Rdata")

# 2. Load Annotation from sql file
library(GEOmetadb)
library(DBI)
require(RH2)

# Retrieve information: 
# (a) Load all the information 
#sqlfile = getSQLiteFile() # It download all the information about the bioconductor annotation
load('sqlfile.Rdata')
load('GEOmetadb.sqlite') # File with the path to the sqlLite database 
con = dbConnect(RSQLite::SQLite(), sqlfile)
plat.id = unlist(lapply(1:length(list.data.sets), function(x) levels(list.data.sets[[x]]$gse[[1]]$platform_id)))
all.bioconductor = lapply(plat.id, function(x) dbGetQuery(con,paste("select gpl,title,bioc_package from gpl where gpl=",paste(x,"'",sep=''),sep="'")  )); names(all.bioconductor) = plat.id
dbDisconnect(con)

# (b) annotate information of the list
an.data.list <- lapply(1:length(list.data.sets), function(x) 	add.annotation.bioconductor(list.data.sets[[x]], all.bioconductor[[x]]))
names(an.data.list) <-paste(sub(' ','',data.sets$Type), data.sets$GEOid,sep='.')
#save(an.data.list, file='list.anno.Rdata')

# (c) Filter list only those with annotation values
has.annotation <-  unlist(lapply(an.data.list, function(x)  length(x) > 1))
filt.data.list <- an.data.list[has.annotation]
# save(filt.data.list, file='Filt.Normalized.Rdata')
# filt.data.list2 <- filt.data.list[-6] # This is not well normalized, so I drop it

########################################################### 
#### LOAD PREVIOUS DATASETS  ####
########################################################### 

form.dat = list(
  'meth'=read.csv2('~/OneDrive/TFM_05.10.16/former_results/previous_data/DatosOriginales-Recodificados/DatosMetilacion-Recoded.csv'),
  'exp'=read.csv2('~/OneDrive/TFM_05.10.16/former_results/previous_data/DatosOriginales-Recodificados/DatosMicroarrays-Recoded.csv'))

geo.form.dat = list(
  'meth' = read.csv2('~/OneDrive/TFM_05.10.16/former_results/previous_data/Datos-GEO/GEOMethData.csv'),
  'exp' = read.csv2('~/OneDrive/TFM_05.10.16/former_results/previous_data/Datos-GEO/GEOExpData.csv')
)
former.dat = list('form'=form.dat,'geo.form'=geo.form.dat)
save(former.dat, file='previous.data.Rdata')

########################################################### 
          ### CORRELATION OF THE DATA ###
########################################################### 

# 1. Merge per Gene
# Load Human Information, we will work in human
library(org.Hs.eg.db)
human.annotation  <- data.frame('ENTREZID'=names(as.list(org.Hs.egSYMBOL)), 'SYMBOL'=as.vector(unlist(as.list(org.Hs.egSYMBOL))))
# save(human.annotation, file='human.annotation.Rdata')

# Filter Table only values of Cancer
filt.data.list3 = filt.data.list
filt.data.list3[[1]] = filt.data.list3[[1]][,grep('mock', colnames(filt.data.list3[[1]]), invert=T)]
filt.data.list3[[2]] = filt.data.list3[[2]]
filt.data.list3[[3]] = filt.data.list3[[3]][,grep('control', colnames(filt.data.list3[[3]]), invert=T)]
filt.data.list3[[4]] = filt.data.list3[[4]][,grep('FHC', colnames(filt.data.list3[[4]]), invert=T)]
filt.data.list3[[5]] = filt.data.list3[[5]][,grep('colon', colnames(filt.data.list3[[5]])) ]
filt.data.list3[[5]] = filt.data.list3[[5]][,grep('mock', colnames(filt.data.list3[[5]]), invert = T) ]
filt.data.list3[[6]] = filt.data.list3[[6]][,grep('Control', colnames(filt.data.list3[[6]]), invert = T) ]
filt.data.list3[[7]] = filt.data.list3[[7]][, grep('Control', colnames(filt.data.list3[[7]]), invert=T)]
filt.data.list3[[8]] = filt.data.list3[[8]][, grep('Control', colnames(filt.data.list3[[8]]), invert=T)]
filt.data.list3[[9]] = former.dat$geo.form$meth ; names(filt.data.list3)[9] <- 'METH.geoFormer'
filt.data.list3[[9]]$SYMBOL = rownames(filt.data.list3[[9]])
filt.data.list3[[10]] = former.dat$geo.form$exp ; names(filt.data.list3)[10] <- 'EXP.geoFormer'
filt.data.list3[[10]]$SYMBOL = rownames(filt.data.list3[[10]])

save(filt.data.list3,file='dataset.full.onlyTumor.RData')

# Table with expression/Meth values
values.data.list = lapply(filt.h.anot, function(x) merge.gene.expMeth(x, filt.data.list3) )

# Opt: Plot Correlation expression + Methylation 
scatter.expMeth(merge.gene.expMeth(human.annotation$ENTREZID[10], filt.data.list3))

# 2. Add hsape each gene in the list
# Reduce Presence Gene
list.gene <- unique(unlist(lapply(1:length(filt.data.list3), function(x) unique(unlist(as.vector(filt.data.list3[[x]]$SYMBOL))))))
filt.h.anot <- as.vector(rownames(human.annotation[human.annotation$SYMBOL %in% list.gene, ]))

# Simple
shape.val <-  lapply(filt.h.anot, function(x) get.shape(x, filt.data.list3))

# Parallel computing
library(parallel)
library(snow)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, "get.shape" )
clusterExport(cl, "shape.param" )
clusterExport(cl,"merge.gene.expMeth")
clusterExport(cl,"filt.h.anot")
clusterExport(cl,"human.annotation")
clusterExport(cl,"filt.data.list3")
shape.val <-  parLapply(cl = cl, filt.h.anot, function(x) get.shape(x, filt.data.list3))
stopCluster(cl)

# Give names to the files
names(shape.val) <-  as.vector(human.annotation$SYMBOL[filt.h.anot])
# save(shape.val,file='shape.val.Rdata')

# Test Shape
testShape = lapply(c(10,20,50,100,500,1000), function(x) lshape.diagnostics(x))
table.testShape = table(NA,NA,NA,NA)

for(i in testShape) {
  table.testShape = rbind(table.testShape, i)
}
table.testShape = as.data.frame(table.testShape)
colnames(table.testShape) = colnames(testShape[[1]])
rownames(table.testShape) = paste("N",c(10,20,50,100,500,1000),sep="=")
write.csv2(table.testShape[,c(4:12)], 'tableTest.shape.csv')

# Clean shape table
shape.table <-  table(NA,NA,NA)
for (i in shape.val) {  shape.table = rbind(shape.table, i ) }
is.null.vector = as.vector(unlist(lapply(1:length(shape.val), function(x) is.null(shape.val[[x]])  )))
shape.tabl <-  as.data.frame(shape.table)
rownames(shape.tabl) <- make.names(names(shape.val)[!is.null.vector], unique=TRUE)  
# save(shape.tabl, file='shape.tabl_results.Rdata')
shape.tabl = shape.tabl[, c(1, grep('.[sp1]', colnames(shape.tabl),invert = T))]

# Get some features
shape.range.our.data = rbind(apply(shape.tabl[,c(4:12)],2, min), apply(shape.tabl[,c(4:12)],2, mean), apply(shape.tabl[,c(4:12)],2, max))
rownames(shape.range.our.data) = c('Min','Mean','Max')
plot(c(1:ncol(shape.range.our.data)), shape.range.our.data[1,], 
     type='b',ylim=c(0,max(shape.range.our.data)),xlab='',ylab='', main='Shape Range: Min-Mean-Max',col='darkred')
lines(c(1:ncol(shape.range.our.data)), shape.range.our.data[2,], type='b', col='darkblue')
lines(c(1:ncol(shape.range.our.data)), shape.range.our.data[3,], type='b', col='darkgreen')

# (Opt) Annotate those which represent L-shape
# aux.fun = function(x) { x - lshape.diagnostics(8) }
# diag.table <- abs(as.data.frame(do.call("rbind", apply(shape.tabl, 1, aux.fun)), stringsAsFactors = FALSE))
#save(diag.table, file='diag.table_results.Rdata')

# (Opt) Contrast Genes Related in Colon Cancer
param.shape.set <- param.eval()
cc.gene <-  read.csv('~/OneDrive/TFM_05.10.16/Results/gene_target.csv', stringsAsFactors = F)
prop.table(table(cc.gene$SYMBOL %in% rownames(param.shape.set.clas))) * 100


# (Opt) Machine learning algoritm to classify the shape of the genes, based in all factors
# Al final he decido utilizar las 3 diferentes features de scagnostics. De eso, la mejor manera para clasificar los datos, resulta ser utilizar SVM y dado eso los valores de los parámetros que obtenemos son. 

param.shape.set.clas <-  param.shape.set[param.shape.set$L.shape == "L", ]
scatter.expMeth(merge.gene.expMeth(human.annotation$ENTREZID[6444], filt.data.list3))

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
    top.go = go.enrichment(val.pred, 'hg19', 1) 
    return(list('val'=val.pred, 'cor'=cor.gene, 'go'=top.go))
}


# Valores defaut con mas y menos cross validations 
pred.def.cv2 = svn.classification(param.eval(stringy = 0.63, convex = 0.1, skinny = c(0.28,0.73)), 2)
pred.def.cv3 = svn.classification(param.eval(stringy = 0.63, convex = 0.1, skinny = c(0.28,0.73)), 3)
pred.def.cv4 = svn.classification(param.eval(stringy = 0.63, convex = 0.1, skinny = c(0.28,0.73)), 4)
pred.def.cv10 =svn.classification(param.eval(stringy = 0.63, convex = 0.1, skinny = c(0.28,0.73)), 10)

# No mejora nada cambiando el cross validation: 199 valores. 
# Valores de accuracy y otros parámetros
options(digits = 3)
data.frame(
  'cv2'=c(pred.def.cv2$conf$overall,pred.def.cv2$conf$byClass),
  'cv3'=c(pred.def.cv3$conf$overall,pred.def.cv3$conf$byClass),
  'cv4'=c(pred.def.cv4$conf$overall,pred.def.cv4$conf$byClass),
  'cv10'=c(pred.def.cv10$conf$overall, pred.def.cv10$conf$byClass)
)

# Valores variados de skinny , cross validation 3 
pred.sky0_30.cv3 = svn.classification(param.eval(skinny = c(0,0.3)), 3)
pred.sky60_100.cv3 = svn.classification(param.eval(skinny = c(0.6,1)), 3)
pred.sky0_100.cv3 = svn.classification(param.eval(skinny = c(0,1)), 3)

#Variación  
data.frame(
  'def' = c(pred.def.cv3$conf$overall,pred.def.cv3$conf$byClass),
  'r0.30'=c(pred.sky0_30.cv3$conf$overall,pred.sky0_30.cv3$conf$byClass),
  'r60.100'=c(pred.sky60_100.cv3$conf$overall,pred.sky60_100.cv3$conf$byClass),
  'r0.100'=c(pred.sky0_100.cv3$conf$overall,pred.sky0_100.cv3$conf$byClass)
)

# El accuracy es similiar en todos. Sin embargo la sensibliidad, la predición el valor de positivos predichos y negativos predichos es mayor en el rango de r60.100. Curiosamente, si tenemos en cuenta más factores como correlación, el mutual information, dist cor o splines. Observamos que la sensiblidad y sensitividad aumentan en r60-100 en skinny

# Ocurre también que cuando disminuimos por debajo del valor de 30, la mayoria de muestras desaparecen.
# Conclusión es que si somos muy relajados en skinny será dificil distinguir entre buenos o no buenos resutalado. 
# Pero dado que parece que esta noe s determinante el rango es 0.2 - 1 default

#	Valores variados convex, cv 3 
pred.conv50.cv3 = svn.classification(param.eval(convex = 0.5 ), 3)
pred.conv70.cv3 = svn.classification(param.eval(convex = 0.7), 3)
pred.conv100.cv3 = svn.classification(param.eval(convex = 1), 3)

data.frame(
  'def' = c(pred.def.cv3$conf$overall,pred.def.cv3$conf$byClass),
  'r50'=c(pred.conv50.cv3$conf$overall,pred.conv50.cv3$conf$byClass),
  'r70'=c(pred.conv70.cv3$conf$overall,pred.conv70.cv3$conf$byClass),
  'r100'=c(pred.conv100.cv3$conf$overall,pred.conv100.cv3$conf$byClass)
)

# Si aumentamos mucho el valor de convex todos los elementos desaparecen. Por lo tanto nuestro valor de datos por 0.1 es el correcto


#	Valores variados stringy, cv 3
pred.stry10.cv3 = svn.classification(param.eval( stringy = 0.1), 3)
pred.stry50.cv3 = svn.classification(param.eval( stringy = 0.5), 3)
pred.stry90.cv3 = svn.classification(param.eval( stringy = 0.9), 3)

data.frame(
  'def' = c(pred.def.cv3$conf$overall,pred.def.cv3$conf$byClass),
  'r10'=c(pred.stry10.cv3$conf$overall,pred.stry10.cv3$conf$byClass),
  'r50'=c(pred.stry50.cv3$conf$overall,pred.stry50.cv3$conf$byClass),
  'r100'=c(pred.stry90.cv3$conf$overall,pred.stry90.cv3$conf$byClass)
)

# Una vez determinados los valores más aproximados. 

interp.classification(svn.classification(param.eval(stringy = 0.1, convex = 0.9, skinny = c(0.6,0.8)), 2)$data, cc.gene)$cor

# Las conclusiones finales son las siguienstes: Stringy, convex y skinny son aquellos parámetros que definien el shape. De estos, no he conseguido probar que ser relaciónen entre si. 

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
# We can use this step to remove data, so for instance if not above 0 out # Pero perse ya hay poco, por 
test.stringy <-  lapply(seq(0,1,0.1), function(x) pred.based.Cor(param.eval(stringy = x)) ) 

# (Opt) We can check the top 10 shape data by using a GO ontology 
random10 <- param.shape.set.clas[sample(rownames(param.shape.set.clas), size = 150), ]
go.enrichment(random10, 'hg19', 0.05)

# (Opt) We check the data using our set control 
comp.datasets(set.to.compare = filt.data.list[[1]], controlPos = control, samplePos = sample)

