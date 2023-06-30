
# Rozwiazanie do zadania z transkryptomiki

setwd("~/EDUGEN_2023/Zadanie_transkryptomika/EDUGEN_2023_zdalne/featurefiles")

############################## WCZYTANIE WYNIKÓW ###############################
counts<-
  read.table("feature_counts_MUT_WT", 
             sep="\t",
             fill=TRUE,
             header=TRUE,
             skip=1,
             stringsAsFactors=FALSE,
             quote = "",
             comment.char = "",
             check.names=FALSE)

########################### MODYFIKACJA NAZW KOLUMN ############################
columns.with.results <- 
  grep(".bam",colnames(counts))

a<-gsub(".bam","",colnames(counts))
colnames(counts)<-a

######################## TABELA Z WYNIKAMI I NAZWAMI GENÓW #####################

mycounts <- counts[, columns.with.results]

rownames(mycounts) <- counts[,1]

library(biomaRt)
mart = useMart("ensembl",
               dataset="hsapiens_gene_ensembl",
               host='jul2018.archive.ensembl.org')

gen<- getBM(values = rownames(mycounts), 
            filters = "ensembl_gene_id", 
            mart = mart, 
            attributes = c("ensembl_gene_id",
                           "hgnc_symbol",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "percentage_gene_gc_content",
                           "gene_biotype",
                           "description",
                           "entrezgene"))

##################################### ANALIZA ##################################

library("DESeq2")

condition <- factor(c(rep("WT",3),
                      rep("MUT",3)))

coldata <- data.frame(row.names=colnames(mycounts), 
                      condition)

dds1 <- DESeqDataSetFromMatrix(countData=mycounts, 
                               colData=coldata, 
                               design=~condition)   

dds1 <- DESeq(dds1)
fr1<-results(dds1,contrast = c("condition","MUT","WT"))

################################### ZADANIE 1 ##################################

###Zadanie1 PODSTAWOWE: Proszę sprawdzić ile mamy obserwacji istotnych po padj
##(P-value z poprawką) poniżej 0.05
###Jaki jest rozkład - histogram wartości pvalue oraz padj i log2FoldChange
    
library(ggplot2)
library(patchwork) 
length(which(fr1$padj<0.05))
g1 <- ggplot(as.data.frame(fr1), aes(x=padj))+geom_histogram()
g2 <- ggplot(as.data.frame(fr1), aes(x=pvalue))+geom_histogram()
g3 <- ggplot(as.data.frame(fr1), aes(x=log2FoldChange))+geom_histogram()
    
g1+g2+g3
    
  
#### Zadanie 1 Dodatkowe 
## Jakie są nazwy NCBI (Hugo names) genów, które są istotnie zmienione,
## co najmniej 5x podwyższone w MUT w skali log2FoldChange

w<-which(fr1$log2FoldChange>5 & fr1$padj<0.05)

geny.nazwy.ncbi <- gen[match(rownames(fr1)[w],gen[,1]),2]
geny.nazwy.ncbi <- geny.nazwy.ncbi[!is.na(geny.nazwy.ncbi)]
geny.nazwy.ncbi <- geny.nazwy.ncbi[grep('\\S', 
                                        geny.nazwy.ncbi, 
                                        invert = FALSE, 
                                        perl = TRUE)] #usunięcie pustych wartości oraz NA

geny.nazwy.ncbi

############################# DALSZA ANALIZA ##################################

setwd("~/EDUGEN_2023/Zadanie_transkryptomika/EDUGEN_2023_zdalne/featurefiles/Second_folder")
      
m<-match(fr1@rownames,gen[,1])

res<-cbind(fr1@rownames,gen[m,],fr1$log2FoldChange,fr1$pvalue,fr1$padj)


colnames(res)<-c("ens", colnames(gen),
                 "log2FoldChange_MUT_WT",
                 "pvalue_MUT_WT",
                 "padj_MUT_WT")

write.table(res,"stats_DESeq2.txt",
            sep="\t",row.names = TRUE, col.names =NA)

df <- as.data.frame(counts(dds1, normalized=TRUE))

hist(df[,1])

log2<-log2(df+1)

write.table(log2,"results_DESeq2_log2_FPKM_plus_1.txt",
            sep="\t",row.names = TRUE,col.names =NA)

################################### ZADANIE 2 ##################################

# Zadanie 2 PODSTAWOWE: Proszę o zrobenie wykresu p-value adjusted (oś y) od 
# log2FoldChange (oś x)
# Proszę zrobić to samo, ale wartość p-value adjusted, przedstawić w formie 
# -log10(padj)
# Na wykresie proszę spróbować zaznaczyć 2 geny które miały co najmniej 
# 5x podwyższoną ekspresje w MUT w skali log2FoldChange (z poprzedniego zadania!)

# Wykres pvalue adjusted vs log2FoldChange
ggplot(as.data.frame(fr1), aes(x = log2FoldChange, y=padj)) + 
  geom_point()

# Wykres -log10 pvalue adjusted vs log2FoldChange
ggplot(as.data.frame(fr1), aes(x = log2FoldChange, y=-log10(padj))) + 
  geom_point()

# Zaznaczenie genów, które miały 5x podwyższoną ekspresję w MUT
ggplot() + 
  geom_point(data = as.data.frame(fr1), 
             aes(x = log2FoldChange, 
                 y=-log10(padj))) +
  geom_point(data = as.data.frame(fr1[w,]),
             aes(x = log2FoldChange, 
                 y=-log10(padj)),
             color = 'red')

# Zadanie 2 DODATKOWE: Proszę policzyć medianę dla  próbek MUT 
# i odobno dla próbek WT
# Następnie proszę zrobić wykres wartości mediany dla MUT i WT
# Na wykresie proszą zaznaczyć geny, które miały co najmniej 5x 
# podwyższoną ekspresje w MUT w skali log2FoldChange (z poprzedniego zadania!)

medWT  <-sapply(1:nrow(log2), function(i) median(as.numeric(log2[i,1:3])))
medMUT <-sapply(1:nrow(log2), function(i) median(as.numeric(log2[i,4:6])))
mediany <- cbind(medWT, medMUT)

ggplot()+ 
  geom_point(data = as.data.frame(mediany),
             aes(x = medWT, 
                 y = medMUT),
             color = 'grey') +
  geom_point(data = as.data.frame(mediany[w,]),
             aes(x = medWT, 
                 y = medMUT), 
             col='red') +
  theme_minimal() +
  labs(x = 'WT', y='MUT', title = 'log2FoldChange plot MUT/WT')

##Dlaczego 1 krpoka jest poniżej 5? skoro różnica miała być >5 log2FoldChange ?
##W DESeq inaczej liczona pewnie była log2FoldChange proszę spojrzeć,że
median(as.numeric(df[w[1],1:3]))
##równa jest 0, a log2(0)=-Inf
mean(as.numeric(df[w[1],1:3]))

################################# PCA ANALYSIS #################################

WT<-c(1:3)
MUT<-c(4:6)

# Rysujemy PCA z naszych danych - wersja dla leniwych - DESeq2
# 1 krok - normalizacja log ze stabilizacją wariancji oraz 
# w odniesieniu do library size
tmp<-rlog(dds1)

# Rysowanie samego wykresu PCA
plotPCA(tmp, intgroup="condition")

# Analiza PCA w wersji bardziej klasycznej
## Analiza PCA zapisana do obiektu pca_matrix
pca_matrix<-prcomp(t(log2), center=TRUE)

### Zobaczmy sobie podsumowanie obiektu pca_matrix
summary(pca_matrix)

###Jaka jest struktura pca_matrix ?
str(pca_matrix)

### Macierz wyników PCA znajduje się wmacierzy pca_matrix$x
pca_matrix$x[1:5,]

### Zobaczmy co jest w środku tej analizy

### Rysujemy barplot z ważnością (importance) naszych głównych składowych
barplot(summary(pca_matrix)$importance[2,],ylab="% of variance",cex.lab=2)

### Rysujemy wykres próbek w przestrzeni 1 i 2 z głównych składowych
par(mfrow=c(1,1),mar=c(5,6,3,2))
plot(pca_matrix$x[,c(1,2)], cex.lab=3,cex.axis=2,bty="l")

### Przypisujemy kolory do naszych próbek
points(pca_matrix$x[WT,1],pca_matrix$x[WT,2], pch = 19, col = "green", cex = 2)
points(pca_matrix$x[MUT,1],pca_matrix$x[MUT,2], pch = 19, col = "red", cex = 2)

### Przypisujemy legendę do naszego wykresu
legend("topleft", legend=c("WT","MUT"),
       col=c("green","red"),pch=19,cex=1.5,bty="n")

### Przypisujemy nazwy próbek do wykresu
text(pca_matrix$x[,1],pca_matrix$x[,2],
     labels=c("WT1","WT2","WT3","MUT1","MUT2","MUT3"),pos=4)

### Mamy problem z opisami, naprawmy
par(mfrow=c(1,1),mar=c(5,6,3,2))
plot(pca_matrix$x[,c(1,2)], cex.lab=3,cex.axis=2,bty="l")
points(pca_matrix$x[WT,1],pca_matrix$x[WT,2], pch = 19, col = "green", cex = 2)
points(pca_matrix$x[MUT,1],pca_matrix$x[MUT,2], pch = 19, col = "red", cex = 2)
legend("topleft", legend=c("WT","MUT"),
       col=c("green","red"),pch=19,cex=1.5,bty="n")
text(pca_matrix$x[,1],pca_matrix$x[,2],
     labels=c("WT1","WT2","WT3","","",""),pos=4)
text(pca_matrix$x[,1],pca_matrix$x[,2],
     labels=c("","","","MUT1","MUT2","MUT3"),pos=2)

##Analizę PCA możemy też zrobić dla genów, a nie dla próbek
pca_matrix<-prcomp(log2, center=TRUE)
plot(pca_matrix$x[,c(1,2)], cex.lab=3,cex.axis=2,bty="l")

##################################### ZADANIE 3 ################################ 
# Narysujmy wykres PCA dla 3 i 4 głównej składowej zarówno dla genów, 
# jak i dla próbek

# Wykres dla próbek:
pca_matrix<-prcomp(t(log2), center=TRUE)
par(mfrow=c(1,1),mar=c(5,6,3,2))
par(mfrow=c(1,1),mar=c(5,6,3,2))
plot(pca_matrix$x[,c(3,4)], cex.lab=3,cex.axis=2,bty="l")
points(pca_matrix$x[WT,3],pca_matrix$x[WT,4], pch = 19, col = "green", cex = 2)
points(pca_matrix$x[MUT,3],pca_matrix$x[MUT,4], pch = 19, col = "red", cex = 2)
legend("topleft", legend=c("WT","MUT"),
       col=c("green","red"),pch=19,cex=1.5,bty="n")
text(pca_matrix$x[,3],pca_matrix$x[,4],
     labels=c("","WT2","WT3","","MUT2","MUT3"),pos=1)
text(pca_matrix$x[,3],pca_matrix$x[,4],
     labels=c("WT1","","","MUT1","",""),pos=4)

### Wykres dla genów
pca_matrix<-prcomp(log2, center=TRUE)
plot(pca_matrix$x[,c(3,4)], cex.lab=3,cex.axis=2,bty="l")


################################# VOLCANO PLOTS ################################

library(EnhancedVolcano)

EnhancedVolcano(fr1,
                lab=rownames(fr1),
                x="log2FoldChange",
                y="padj")

################################ ZADANIE 4 #####################################

# Zróbmy analogiczny wykres ale z nazwami NCBI (hugo names) zamiast 
# id z ENSEMBL-a
m<-match(rownames(fr1),gen[,1])
res3<-fr1[which(!is.na(m)),]
rownames(res3)<-gen[na.omit(m),2]
EnhancedVolcano(res3,
                lab=rownames(res3),
                x="log2FoldChange",
                y="padj")

# wykres -log10(p) < 40
EnhancedVolcano(res3,
                lab=rownames(res3),
                x="log2FoldChange",
                y="padj",
                ylim=c(0,40))

########################### ANALIZA ŚCIEŻEK GO/KEGG ############################

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(enrichplot)

# Zawężenie genów:
w<-which(fr1$padj<0.05)

# Analiza wzbogaceń w ścieżki sygnałowe Gene Ontology (GO)
# dla wszystkich kategorii GO (ALL): BP, MF i CC

kk1 <- enrichGO(gene = res$entrezgene[w],
                OrgDb="org.Hs.eg.db",
                pvalueCutoff = 0.05,
                pAdjustMethod="fdr",
                ont="ALL")

# Tworzenie wykresu dla wzbogaconych ścieżek 
barplot(kk1, showCategory=10) + ggtitle("GO MF barplot MUT/WT")

# Analiza ścieżek KEGG
# Zawężamy geny tylko do tych istotnych na poziomie istotności 0.05
w<-which(fr1$padj<0.05)

# Analiza wzbogaceń w ścieżki sygnałowe KEGG
kk2 <- enrichKEGG(gene         = res$entrezgene[w],
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)

# Tworzenie wykresu dla wzbogaconych ścieżek 
barplot(kk2, showCategory=30) + ggtitle("KEGG barplot MUT/WT")

################################ ZADANIE 5 #####################################

# Zadanie 5 PODSTAWOWE: Zróbmy analizę ścieżek GO BP wybierając geny, 
# które są istotnie podwyższone w MUT

w<-which(fr1$padj<0.05 & fr1$log2FoldChange>0)
# Analiza wzbogacę w ścieżki sygnałowe Gene Ontology (GO) 
# dla wszystkich kategorii GO (ALL): BP, MF i CC
# BP - biological process
# MF - molecular funcrtion
# CC - cellular process

kk.bp <- enrichGO(gene         = res$entrezgene[w],OrgDb="org.Hs.eg.db",
               pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP")

kk.mf <- enrichGO(gene         = res$entrezgene[w],OrgDb="org.Hs.eg.db",
                  pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="MF")

kk.cc <- enrichGO(gene         = res$entrezgene[w],OrgDb="org.Hs.eg.db",
                  pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="CC")

# Tworzenie wykresu dla wzbogaconych ścieżek 
barplot(kk.bp, showCategory=30) + 
  ggtitle("GO BP barplot Upregulated in MUT vs WT")

barplot(kk.mf, showCategory=30) + 
  ggtitle("GO MF barplot Upregulated in MUT vs WT")

barplot(kk.cc, showCategory=30) + 
  ggtitle("GO CC barplot Upregulated in MUT vs WT")

# Zadanie 5A DODATKOWE: Zróbmy analizę ścieżek REACTOME dla takich samych kryteriów
# jak dla analizy GO powyżej, ale wykres zróbmy kropokowy (dotplot)
# Podpowiedź - użyjmy pakietu ReactomePA

w<-which(fr1$padj<0.05 & fr1$log2FoldChange>0)

library("ReactomePA")
kk.reactome <- enrichPathway(gene=res$entrezgene[w],pvalueCutoff=0.05, readable=F)

dotplot(kk.reactome, showCategory=30) + 
  ggtitle("REACTOME barplot MUT/WT")

# Zadanie 5B DODATKOWE: Proszę z poprzednich analiz GO i KEGG oraz z analizy 
# REACTOME spróbować wygenerować tabelkę z wynikami, która będzie zawierać
# także przypisanie ścieżek do genów w formacie NCBI

# Reactome
kk.reactome <- enrichPathway(gene=res$entrezgene[w],
                             pvalueCutoff=0.05, 
                             readable=T)

o<-order(kk.reactome[,6],
         decreasing = FALSE)

Wyniki.reactome <-as.data.frame(kk.reactome)[o,]

write.table(Wyniki.reactome,
            "REACTOME_MUT_vs_WT.txt",
            sep="\t",
            row.names = TRUE,
            col.names = NA)

# GO

o.GO<-order(kk1[,6],
         decreasing = FALSE)

Wyniki.GO <- as.data.frame(kk1)[o.GO,]

write.table(Wyniki.GO,
            "GO_MUT_vs_WT.txt",
            sep="\t",
            row.names = TRUE,
            col.names = NA)

# KEGG
o.KEGG <- order(kk2[,6],
                decreasing = FALSE)

Wyniki.KEGG <- as.data.frame(kk2)[o.GO,]

write.table(Wyniki.KEGG,
            "KEGG_MUT_vs_WT.txt",
            sep="\t",
            row.names = TRUE,
            col.names = NA)

# Inna ciekawa forma prezentacji wyników - cnetplot
# Powtorzmy operacje z poprzenich linijek:
library(ggnewscale)

w<-which(fr1$padj<0.05)
# Analiza wzbogaceń w ścieżki sygnałowe Gene Ontology (GO) dla wszystkich 
# kategorii GO (ALL): BP, MF i CC

kk.ggnewscale <- enrichGO(gene         = res$entrezgene[w],OrgDb="org.Hs.eg.db",
                pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="ALL")
cnetplot(kk.ggnewscale)

# Nie bardzo czytelne, spróbujemy przypisać nazwy genów
edox<-setReadable(kk.ggnewscale,'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox)

# To też nie do końca to, spróbujmy przypisać wartości log2FoldChange 
# do naszych genów
GENELIST<-res$log2FoldChange_MUT_WT
names(GENELIST)<-res$entrezgene

w<-which(fr1$padj<0.05)
kk.ggnewscale2 <- enrichGO(gene = res$entrezgene[w],
                OrgDb="org.Hs.eg.db",
                pvalueCutoff = 0.05,
                pAdjustMethod="fdr",
                ont="ALL")

# Przypiszmy nazwy hugo names zamiast ENTREZID
edox <- setReadable(kk.ggnewscale2, 'org.Hs.eg.db', 'ENTREZID')

# Wersja "gałązka"
cnetplot(edox, foldChange=GENELIST)

# Wersja kolista
cnetplot(edox, foldChange=GENELIST, circular = TRUE, colorEdge = TRUE)

############### tight junction #########################

  # zmiana z amyloid beta clearance na 'tight junction'.
  # Pierwszy szlak nie znajdował się w wynikach analizy

wszystkie <- kk1@result$Description

DFG<- 'tight junction'
tmp<-unlist(strsplit(as.data.frame(kk1)[which(as.data.frame(kk1)[,3]==DFG),9],"/"))
eg1<-bitr(tmp, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
unique.eg1<- eg1[!duplicated(eg1$ENTREZID), ]
m<-match(unique.eg1[,2],rownames(log2))
lim<-log2[m,]
rownames(lim)<-gen[match(rownames(lim),gen[,1]),2]
colnames(lim)<-c("WT1","WT2","WT3","MUT1","MUT2","MUT3")
library(viridis)
library(gplots)

############################## ZADANIE 6 #######################################
lim <- as.matrix(lim)
heatmap.2(lim, col=inferno(10), scale="none",
          key=TRUE, keysize=0.75,symkey=FALSE, density.info="none",
          trace="none", cexRow=1.2,cexCol=2.5,margins=c(0,0),main=DFG)

# Zadanie 6 DODATKOWE: Po naprawieniu poprzedniego błedu proszę 
# pominąć klastrowanie po wierszach i kolumnach
# Proszę ustawić geny oraz próbki tak jak w oryginalnej data.frame - lim
# Proszę też zmienić skalę kolorów na inną i zwiększyć rozpiętość skali kolorów

heatmap.2(as.matrix(lim), col=viridis(50), scale="none",Rowv=FALSE,Colv=FALSE,
          key=TRUE, keysize=0.75,symkey=FALSE, density.info="none",
          trace="none", cexRow=1.2,cexCol=2.5,margins=c(0,0),main=DFG)

############################## Z SCORE #########################################

zscore<-sapply(1:nrow(lim), 
               function(i) scale(as.matrix(lim)[i,],
                                 center=TRUE, scale=TRUE))

zscore
zscore2<-t(zscore)

rownames(zscore2)<-rownames(lim)
colnames(zscore2)<-colnames(lim)

# Uporządkowana heatmapa
heatmap.2(zscore2, col=viridis(90), scale="none",Rowv=TRUE,Colv=TRUE,
          key=TRUE, keysize=0.75,symkey=FALSE, density.info="none",
          trace="none", cexRow=1.2,cexCol=2.5,margins=c(8,8),main=DFG)

############################## ZADANIE 7 #######################################
w<-which(fr1$padj<0.05)
kk1 <- enrichGO(gene         = res$entrezgene[w],OrgDb="org.Hs.eg.db",
                pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="ALL")
sciezki <- kk1@result$Description
sciezki[c(1,2,3,8)]

co.nas.interesuje <- c(1,2,3,8)

tmp <- c()  
for (i in co.nas.interesuje) {
  nowe <- unlist(strsplit(kk1@result$geneID[i], "/"))
  tmp <- c(tmp, nowe)
}

m<-sapply(1:length(tmp), function(i) which(res$entrezgene==tmp[i]))

############################## ZADANIE 8 #######################################
# Zadanie 8 PODSTAWOWE: Proszę zmapować nazwy z tmp w inny sposób - 
# za pomocą funkcji grep i zapiszcie wyniki do zmiennej jkl

jkl<-sapply(1:length(tmp), function(i) grep(tmp[i],res$entrezgene))

M<-sapply(1:length(jkl), function(i) jkl[[i]][1])


lim<-as.matrix(df[m,])
nrow(res) == nrow(df)
dim(lim)

length(rownames(lim))
length(unique(rownames(lim)))

colnames(lim)<-NULL
rownames(lim)<-NULL

library(ggplot2)
library(viridis)
library(gplots)
library(RColorBrewer)

zscore <-sapply(1:nrow(lim), function(i) scale(lim[i,],center=TRUE, scale=TRUE))
rownames(zscore)<-colnames(lim)
colnames(zscore)<-rownames(lim)
zscore2<-t(zscore)

POP=brewer.pal(9,"Set1")

display.brewer.pal(9,"Set1")

names<-sciezki[c(1,2,3,8)]

par(mar=c(1,1,1,1), xpd=TRUE,cex.main=1.2)
asd<-heatmap.2(zscore2, main = 'GO terms DEGs MUT/WT',
               #dendrogram = "row", # no dendrogram for columns
               #Rowv = "NA", # * use self-made dendrogram
               #Colv = "NA",  # make sure the columns follow data's order
               col = inferno(90),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,0.1),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               RowSideColors = c(rep(POP[5],length(unlist(strsplit(kk1@result$geneID[1],"/")))),rep(POP[6],length(unlist(strsplit(kk1@result$geneID[2],"/")))),rep(POP[7],length(unlist(strsplit(kk1@result$geneID[3],"/")))),rep(POP[8],length(unlist(strsplit(kk1@result$geneID[8],"/"))))), #to add nice colored strips
               ColSideColors = c(rep(POP[1],3),rep(POP[2],3)),
               colsep=c(3,6,9,12,15,18,21,24,27,30,33)
)
legend("left", legend=c("WT","MUT",names), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.8)

par(mar=c(5.1, 6.1, 4.1, 4.1), xpd=TRUE,cex.main=1.2)
asd<-heatmap.2(zscore2, main = 'GO terms DEGs MUT/WT',
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(90),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,0.1),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               RowSideColors = c(rep(POP[5],length(unlist(strsplit(kk1@result$geneID[1],"/")))),rep(POP[6],length(unlist(strsplit(kk1@result$geneID[2],"/")))),rep(POP[7],length(unlist(strsplit(kk1@result$geneID[3],"/")))),rep(POP[8],length(unlist(strsplit(kk1@result$geneID[8],"/"))))), #to add nice colored strips
               ColSideColors = c(rep(POP[1],3),rep(POP[2],3)),
               colsep=c(3,6,9,12,15,18,21,24,27,30,33)
)
legend("left", legend=c("WT","MUT",
                        names), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.7,inset=c(-0.01,0),y.intersp=2)


############################## ZADANIE 9 #######################################

# Zadanie 9 PODSTAWOWE: Zróbcie proszę heatmapę bez klastrowania wierszy ani kolumn, 
# czyli tak jak w oryginalnym zscore
# Zróbcie też tak, żeby nazwy wierszy (liczby) były widoczne
# Jeśli jest problem z za długimi nazwami śieżek, popraecie je 
# analogicznie do "amyloid-beta clearance"

par(mar=c(5.1, 6.1, 4.1, 4.1), xpd=TRUE,cex.main=1.2)
asd<-heatmap.2(zscore2, main = 'GO terms DEGs MUT/WT',
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(90),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,7),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               RowSideColors = c(rep(POP[5],length(unlist(strsplit(kk1@result$geneID[1],"/")))),rep(POP[6],length(unlist(strsplit(kk1@result$geneID[2],"/")))),rep(POP[7],length(unlist(strsplit(kk1@result$geneID[3],"/")))),rep(POP[8],length(unlist(strsplit(kk1@result$geneID[8],"/"))))), #to add nice colored strips
               ColSideColors = c(rep(POP[1],3),rep(POP[2],3)),
               colsep=c(3,6,9,12,15,18,21,24,27,30,33)
)
#legend("left", legend=c("WT","MUT",names), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.8,inset=c(-0.3,0))
legend("left", legend=c("WT","MUT", names), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.7,inset=c(-0.1,0),y.intersp=2)



###Zadanie 9 DODATKOWE: Zróbcie heatmapę jak wyżej ale nadajcie jej nazwy genów - Hugo names
###Posortujcie też geny w obrębie ścieżek po wartości log2FoldChange
###Zwróćcie uwagę na unikalność nazw wierszy...

rownames(lim)<-res$hgnc_symbol[m]
length(rownames(lim))
length(unique(rownames(lim)))
###Wartości niezgodne :)
length(rownames(lim))
length(unique(rownames(lim)))
###Wartości niezgodne
which(duplicated(rownames(lim)))
###Możemy usunąć 2gie wystąpienie genu (zalecane), lub zostawić
lim2<-lim[which(!duplicated(rownames(lim))),]
length(rownames(lim2))
zscore<-sapply(1:nrow(lim2), function(i) scale(lim2[i,],center=TRUE, scale=TRUE))
rownames(zscore)<-colnames(lim2)
colnames(zscore)<-rownames(lim2)
zscore2<-t(zscore)
names<-sciezki[c(1,2,3,8)]

###Zmapujmy nasze geny do obiektu res
fg<-res$hgnc_symbol[m][which(!duplicated(rownames(lim)))]
m1<-match(unlist(strsplit(kk1@result$geneID[1],"/")),res$entrezgene)
m2<-match(unlist(strsplit(kk1@result$geneID[2],"/")),res$entrezgene)
m3<-match(unlist(strsplit(kk1@result$geneID[3],"/")),res$entrezgene)
m4<-match(unlist(strsplit(kk1@result$geneID[8],"/")),res$entrezgene)


###Stwórzmy sobie plik w którym mapujemy geny z log2FoldChane i ze ścieżkami
map<-cbind(c(res$hgnc_symbol[m1],res$hgnc_symbol[m2],res$hgnc_symbol[m3],res$hgnc_symbol[m4]),
           c(res$log2FoldChange_MUT_WT[m1],res$log2FoldChange_MUT_WT[m2],res$log2FoldChange_MUT_WT[m3],res$log2FoldChange_MUT_WT[m4]),
           c(rep(names[1],length(m1)),rep(names[2],length(m2)),rep(names[3],length(m3)),rep(names[4],length(m4))))
###Ograniczmy nasz obiekt map do obuektu map2 pozbawionego zduplikowanych nazw
map2<-map[which(!duplicated(rownames(lim))),]
###Teraz zróbmy sobie wektor porządkujący geny po log2Foldchange

###Tworzenie wektorów numerycznych z pozycjami naszych genów z naszych ścieżek
w1<-which(map2[,3]==names[1])
w2<-which(map2[,3]==names[2])
w3<-which(map2[,3]==names[3])
w4<-which(map2[,3]==names[4])

o1<-order(as.numeric(map2[w1,2]))
o2<-order(as.numeric(map2[w2,2]))
o3<-order(as.numeric(map2[w3,2]))
o4<-order(as.numeric(map2[w4,2]))

lim1<-zscore2[w1[o1],]
lim2<-zscore2[w2[o2],]
lim3<-zscore2[w3[o3],]
lim4<-zscore2[w4[o4],]
###Zwróćcie uwagę, że jeśli któryś z tych obiektów lim powyżej miałby tylko 1 gen, to zgubiłby nazwę wiersza, bo 
###zamiast matrix, będzie zwykłym wektorem numerycznym

ZSCORE<-rbind(lim1,lim2,lim3,lim4)


par(mar=c(5, 4, 4, 2)+.1, xpd=TRUE,cex.main=1.2)
asd<-heatmap.2(ZSCORE, main = 'GO terms DEGs MUT/WT',
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(90),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(1,5),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(15,0,3,1)),
               RowSideColors = c(rep(POP[5],length(w1)),rep(POP[6],length(w2)),rep(POP[7],length(w3)),rep(POP[8],length(w4))), #to add nice colored strips
               ColSideColors = c(rep(POP[1],3),rep(POP[2],3)),
               colsep=c(3,6,9,12,15,18,21,24,27,30,33))

legend("left", 
       legend=c("WT","MUT", names),
       col=c(POP[c(1,2,5,6,7,8)]),
       pch=15,cex=0.7,
       inset=c(-0.1,0),y.intersp=2)


###Wersja gdy chcecie miec spojne wymiary ryciny
pdf(file="test.pdf", height=10, width=10)
par(mar=c(5.1, 8.1, 4.1, 4.1), xpd=TRUE,cex.main=1.2)
asd<-heatmap.2(ZSCORE, main = 'GO terms DEGs MUT/WT',
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's ordĪĪĪer
               col = inferno(90),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(2,7),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               RowSideColors = c(rep(POP[5],length(w1)),rep(POP[6],length(w2)),rep(POP[7],length(w3)),rep(POP[8],length(w4))), #to add nice colored strips
               ColSideColors = c(rep(POP[1],3),rep(POP[2],3)),
               colsep=c(3,6,9,12,15,18,21,24,27,30,33)
)
#legend("left", legend=c("WT","MUT",names), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.8,inset=c(-0.3,0))
legend("left", legend=c("WT","MUT",paste("amyloid-beta","clearance",sep="\n"),
                        names[2:4]), col=c(POP[c(1,2,5,6,7,8)]),pch=15,cex=0.7,inset=c(-0.18,0),y.intersp=2)

dev.off()


################################################################
###Każdy kod do analizy "prawdziwych" danych polecam dokładnie przejrzeć
###Zapuścić od początku do najmniej 2 razy, żeby zobaczyć, że nie pomineliscie zapisanie jakiegos kroku
###Jeśli macie wątpliwości, podziellcie się kodem z kimś, kto go puści i potencjalnie zdebuguje

