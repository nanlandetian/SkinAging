library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
#install.packages('hdf5r')
library(hdf5r)
#install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library('SeuratDisk')

Convert('GSE149590_RAW/GSM4505404_tabula-muris-senis-droplet-official-raw-obj.h5ad','h5seurat',overwrite = T,assay = 'RNA')
scRNA_all<-LoadH5Seurat(file= 'GSE149590_RAW/GSM4505404_tabula-muris-senis-droplet-official-raw-obj.h5seurat')


scRNA_all <- NormalizeData(scRNA_all, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

metadata<-read.csv('GSE149590_RAW/GSM4505404_tabula-muris-senis-droplet-official-raw-obj-metadata.csv',row.names = 1)

row.names(metadata)<-colnames(scRNA_all)
scRNA_all<-AddMetaData(scRNA_all,metadata = metadata)

# scRNA_all@meta.data<-metadata
# table(scRNA_all$cell_ontology_class)
# Idents(scRNA_all)<-'age'
# VlnPlot(scRNA_all,features = 'Kdm6a')

table(scRNA_all$cell_ontology_class)
FOXE1<-rep(FALSE,length(scRNA_all$cell_ontology_class))
scRNA_all$FOXE1<-FOXE1

sce<-subset(scRNA_all,subset = Xpnpep2 > 0)

sce$FOXE1<-'TRUE'
subcell_type_value<-as.vector(sce$FOXE1)
names(subcell_type_value)<-names(sce$FOXE1)
scRNA_all$FOXE1[names(subcell_type_value)]<-subcell_type_value
table(scRNA_all$FOXE1)
table(scRNA_all$cell_ontology_class,scRNA_all$FOXE1)
write.csv(file = 'XPNPEP2_expression_over0_celltype_droplet.csv',table(scRNA_all$cell_ontology_class,scRNA_all$FOXE1))

#aa<-table(scRNA_all$cell_ontology_class,scRNA_all$tissue)
#bb<-as.data.frame.array(table(scRNA_all$cell_ontology_class,scRNA_all$age)) #table数据转宽矩阵

old<-scRNA_all$age
old[old=='1m']<-'young'
old[old=='3m']<-'young'
old[old=='18m']<-'middle'
old[old=='21m']<-'middle'
old[old=='24m']<-'old'
old[old=='30m']<-'old'
table(old)
old<-factor(old,levels = c('young','middle','old'))
scRNA_all$old<-old

Idents(scRNA_all)<-'cell_ontology_class'


scRNAsub = scRNA_all[, scRNA_all@meta.data$cell_ontology_class %in% choosed]
sce = scRNAsub[, scRNAsub@meta.data$old %in% c('old','young')]
table(sce$old)
#sce2 <- subset(sce,subset = Pi16 > 0)


VlnPlot(sce,features = 'Idh2',split.by = 'old',idents = c('fibroblast of cardiac tissue',"fibroblast","fibroblast of lung","pulmonary interstitial fibroblast"))
p<-VlnPlot(scRNA_all,features = 'Lrrc4',idents = c(names(table(scRNA_all$cell_ontology_class))[63:123]))
#p+theme(legend.position = 'none')


choosed<-c("fibroblast of cardiac tissue","fibroblast","fibroblast of lung","pulmonary interstitial fibroblast",
          "mesenchymal stem cell",
          'basal cell',
          "luminal epithelial cell of mammary gland" ,
          "ciliated columnar cell of tracheobronchial tree",
          "club cell of bronchiole" ,
          "hepatocyte",
          "kidney distal convoluted tubule epithelial cell",
          "kidney loop of Henle ascending limb epithelial cell",
          "kidney loop of Henle thick ascending limb epithelial cell",
          "Schwann cell" ,
          "leukocyte" 
          )

genes<-c('Mettl1','Kdm6a','Pi16','Foxe1','Lhx9','Trmt1','Wdr4') #2个基因没有 Plpp3 Epb41l3
p<-VlnPlot(sce,features = genes,split.by = 'old')
p & theme(axis.text.x = element_text(angle = 0,hjust = 0.5))

table(scRNAsub$cell_ontology_class)
TPM_mtx<-GetAssayData(object = scRNAsub, slot = "data")
#intersect(rownames(TPM_mtx),genes) #检查是不是需要提取的所有基因在矩阵中存在
TPM_mtx_filter<-TPM_mtx[genes,] 
TPM_mtx_filter_asMtx<-as.matrix(TPM_mtx_filter)
library(reshape2)
long_mtx<-melt(TPM_mtx_filter_asMtx)
colnames(long_mtx)<-c('gene','cell','value')
meta.data<-scRNAsub@meta.data
meta.data

age<-meta.data$age
names(age)<-rownames(meta.data)
long_mtx$age<-age[long_mtx$cell]

cell_ontology_class<-meta.data$cell_ontology_class
names(cell_ontology_class)<-rownames(meta.data)
long_mtx$cell_ontology_class<-cell_ontology_class[long_mtx$cell]

long_mtx$old[long_mtx$age=='1m']<-'Young'
long_mtx$old[long_mtx$age=='3m']<-'Young'
long_mtx$old[long_mtx$age=='18m']<-'other'
long_mtx$old[long_mtx$age=='21m']<-'Old'
long_mtx$old[long_mtx$age=='24m']<-'Old'
long_mtx$old[long_mtx$age=='30m']<-'Old'

long_mtx$old<-factor(long_mtx$old,levels = c('Young','Old','other'))

table(scRNAsub$cell_ontology_class)

library(ggpubr)
library(dplyr)
long_mtx %>% filter((long_mtx$cell_ontology_class=='hepatocyte') & 
                      (long_mtx$old == 'Young' | long_mtx$old == 'Old')  & long_mtx$gene == 'Kdm6a') %>%        #Kdm6a
  ggplot(aes(x=old,y=value,fill=old)) +  #x=gene
  #geom_violin(scale = 'width',alpha=0.8,colour='White',draw_quantiles = c(0.25, 0.75),linetype = "dashed") + 
  geom_violin(scale = 'width',alpha=0.8,colour='White',draw_quantiles = c(0.25,0.5,0.75)) +
  #geom_violin(scale = 'width',fill="transparent",draw_quantiles = 0.5)
  scale_fill_manual(values = c('#1A00FE','#FB0006'))+
  stat_compare_means(method = 't.test',aes(label= paste0("p = ", ..p.format..)))+
  theme_bw()+
  ylab('Expression')+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(colour = 'black'), #,size = 1
    #axis.line = element_line(colour = 'black'),
    axis.ticks.length = unit(.25,'cm'),
    text = element_text(size = 20,colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  )


long_mtx %>% filter((long_mtx$cell_ontology_class=='hepatocyte') & 
                      (long_mtx$old == 'Young' | long_mtx$old == 'Old')  & long_mtx$gene == 'Mettl1' 
                    ) %>%  
  ggplot(aes(x=old,y=value))+ 
  geom_jitter(position=position_jitter(0.2),color="#196D6E") +
  stat_summary(fun.data=data_summary) +
  #geom_violin(scale = 'width',fill="transparent",draw_quantiles = 0.5)
  scale_fill_manual(values = c('#1A00FE','#FB0006'))+
  stat_compare_means(method = 't.test',aes(label= paste0("p= ", ..p.format..)))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.length = unit(.25,'cm')
  )




  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }




long_mtx %>% filter((long_mtx$cell_ontology_class=='mesenchymal stem cell') & long_mtx$old == 'Young' & long_mtx$gene == 'Mettl1' ) %>%  
  ggplot(aes(gene,value,fill=old)) + 
  geom_violin(scale = 'width',alpha=0.5,draw_quantiles = c(0.25, 0.5, 0.75))

quantile(aa$value,probs=c(0.25,0.5,0.75))
quantile(aa$value,0.75) + 1.5*IQR(aa$value)
?IQR()


df <- data.frame(x = 1, y = c(0, 0.25, 0.5, 0.75, 5))
p<-ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_violin(draw_quantiles = c(0.25,0.5,0.75))
quantile(df$y,probs=c(0.25,0.5,0.75))


create_quantile_segment_frame <- function(data, draw_quantiles) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y, ties = "ordered")
  ys <- ecdf(draw_quantiles) # these are all the y-values for quantiles
  
  # Get the violin bounds for the requested quantiles.
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  
  # We have two rows per segment drawn. Each segment gets its own group.
  ggplot2:::new_data_frame(list(
    x = ggplot2:::interleave(violin.xminvs, violin.xmaxvs),
    y = rep(ys, each = 2),
    group = rep(ys, each = 2)
  ))
}
assignInNamespace("create_quantile_segment_frame", create_quantile_segment_frame, "ggplot2")
df <- data.frame(x = 1, y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10))
#<<<<<<

################某个基因在各个亚群中表达大于0的细胞数统计 #########
# scRNA_all@meta.data<-metadata
# table(scRNA_all$cell_ontology_class)
# Idents(scRNA_all)<-'age'
# VlnPlot(scRNA_all,features = 'Kdm6a')

table(scRNA_all$cell_ontology_class)
FOXE1<-rep(FALSE,length(scRNA_all$cell_ontology_class))
scRNA_all$FOXE1<-FOXE1

sce<-subset(scRNA_all,subset = Xpnpep2 > 0)

sce$FOXE1<-'TRUE'
subcell_type_value<-as.vector(sce$FOXE1)
names(subcell_type_value)<-names(sce$FOXE1)
scRNA_all$FOXE1[names(subcell_type_value)]<-subcell_type_value
table(scRNA_all$FOXE1)
table(scRNA_all$cell_ontology_class,scRNA_all$FOXE1)
write.csv(file = 'XPNPEP2_expression_over0_celltype_droplet.csv',table(scRNA_all$cell_ontology_class,scRNA_all$FOXE1))

################# 某个基因在部分亚群（有细胞表达）中的表达 #########
choosed<-c("keratinocyte",
           "basal cell of epidermis",
           "mesenchymal stem cell",
           "kidney proximal convoluted tubule epithelial cell",
           'basal cell',
           "stromal cell" ,
           "bronchial smooth muscle cell",
           "mesenchymal stem cell of adipose" ,
           "bladder cell",
           "bladder urothelial cell",
           "fibroblast of cardiac tissue",
           "epithelial cell of proximal tubule",
           "mesenchymal progenitor cell" ,
           "enterocyte of epithelium of large intestine",
           "adventitial cell",
           'epidermal cell',
           "fibroblast of lung",
           "large intestine goblet cell",
           "fibroblast"
)

scRNAsub = scRNA_all[, scRNA_all@meta.data$cell_ontology_class %in% choosed]
table(scRNAsub$cell_ontology_class)
Idents(scRNAsub)<-'cell_ontology_class'
VlnPlot(scRNAsub,features = 'Xpnpep2')

###
choosed<-c(
           "fibroblast of cardiac tissue",
           "fibroblast of lung",
           "fibroblast",
           "stromal cell" ,
           "mesenchymal stem cell of adipose" ,
           "mesenchymal progenitor cell" ,
           "mesenchymal stem cell"
           
)
scRNAsub = scRNA_all[, scRNA_all@meta.data$cell_ontology_class %in% choosed]
table(scRNAsub$cell_ontology_class)
Idents(scRNAsub)<-'cell_ontology_class'
VlnPlot(scRNAsub,features = 'Xpnpep2',split.by = 'old')

