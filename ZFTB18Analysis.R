DatosMetilacion <- read.csv2("dades/DatosMetilacion.csv")
Met<-DatosMetilacion
rownames(Met)<-Met$X
Met<- as.matrix(Met[,-1])
M1<-Met["ZBTB18",]

arrays <- read.csv2("dades/DatosMicroarrays.csv")
rownames(arrays)<-arrays$X
arrays<- as.matrix(arrays[,-1])
E1<-arrays["ZBTB18",]

plot(M1,E1)
cor(M1,E1)

rnaseq <- read.csv2("C:/Users/asanc/Dropbox (VHIR)/TreballsMeus/2017-01-Selecting_GRM/dades/DatosRNAseq.csv")
rownames(rnaseq)<-rnaseq$X
rnaseq<- as.matrix(rnaseq[,-1])
E2<-rnaseq["ZBTB18",]

plot(M1,E2)
cor(M1,E2)

###################################################

library(readxl)
final_list <- as.data.frame(read_excel("Genes regulated by methylation - final list.xlsx"))
rownames(final_list)<-final_list$gene
final_list<- as.matrix(final_list[,-1])
dim(final_list)
expres<- final_list[,1:30]
mets<- final_list[,31:60]
met2<- mets["ZNF238",]
expres2<-expres["ZNF238",]
plot(met2,expres2)
