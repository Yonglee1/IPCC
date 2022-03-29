# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' ipcc_cpdb,using cellphonedb to analysis
#'
#' @param input your scRNA
#' @param species mouse or human
#'
#' @return figures and txt analyzed by cellphonedb
#' @export
#'
#' @examples ipcc_cpdb(input= scRNA,species = "mouse")
ipcc_cpdb <- function(input,species) {

  library(Seurat)
  scRNA <- input
  if(species=="mouse"){

    ###基因转换,mouse
    library(biomaRt)
    human = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useEnsembl(biomart="ensembl", dataset = "mmusculus_gene_ensembl")
    count_data <- as.matrix(scRNA@assays$RNA@data)
    ms.gene <- rownames(count_data)
    mtoh <- getLDS(     #使用getLDS()将基因进行转换
      values = ms.gene,mart = mouse,
      attributes = "mgi_symbol",filters = "mgi_symbol",
      martL=human,
      attributesL = c("hgnc_symbol","chromosome_name"))
    count_data <- count_data[rownames(count_data)%in%c(mtoh$MGI.symbol),]#留下转换后的基因的矩阵
    mtoh <- mtoh[mtoh[,2]!="",]#去除人类基因的空字符
    mtoh <- mtoh[c(match(rownames(count_data),mtoh$MGI.symbol)),]#match相同的基因，输出
    rownames(count_data) <- mtoh$HGNC.symbol#替换基因名
    count_data <- as.data.frame(count_data)
    write.table(as.matrix(count_data), 'cellphonedb_count.txt', sep='\t', quote=F)

  }

  else{
    #####不用转换，human
    write.table(as.matrix(scRNA@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
  }

  ###cellphonedb
  meta_data <- cbind(rownames(scRNA@meta.data), scRNA@meta.data[,'celltype', drop=F])
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA

  write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

  shell_cmd<-"cellphonedb method statistical_analysis  cellphonedb_meta.txt  cellphonedb_count.txt --counts-data=gene_name"
  grep_out<-system(shell_cmd)
  cat(grep_out)

  mypvals <- read.table("out/pvalues.txt",
                        header = T,sep = "\t",stringsAsFactors = F)
  mymeans <- read.table("out/means.txt",
                        header = T,sep = "\t",stringsAsFactors = F)

  sm = as.data.frame(
    do.call(rbind,
            lapply( 12:ncol(mypvals) , function(i){
              return(c( strsplit(colnames(mypvals)[i],'\\.')[[1]],
                        sum(mypvals[,i] <0.05)))
            }))
  )

  colnames(sm)=c('SOURCE' ,'TARGET' ,'count')
  sm$count = as.numeric( sm$count )

  write.table(sm,file = 'count_network.txt',
              sep = '\t',
              quote = F,row.names = F)

  library(reshape2)
  sm_df =dcast(as.data.frame(sm),SOURCE~TARGET )
  sm_df[is.na(sm_df)]=0
  rownames(sm_df) = sm_df[,1]
  sm_df = sm_df[,-1]
  p1=pheatmap::pheatmap(sm_df,display_numbers = T,cellwidth = 20,cellheight = 20)
  p2=pheatmap::pheatmap(log(sm_df+1),display_numbers = T,,cellwidth = 20,cellheight = 20)
  library(patchwork)
  library(cowplot)
  library(ggplotify)
  as.ggplot(p1) + as.ggplot(p2)

  ggplot2::ggsave('out/cpdb_heatmap.png',width = 10,height = 5)

  kp = grepl(pattern = "B", colnames(mypvals)) | grepl(pattern = "空", colnames(mypvals))
  table(kp)
  pos = (1:ncol(mypvals))[kp]
  choose_pvalues <- mypvals[,c(c(1,5,6,8,9),pos  )]
  choose_means <- mymeans[,c(c(1,5,6,8,9),pos)]

  logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.05, 1, sum)
  # 只保留具有细胞特异性的一些相互作用对
  choose_pvalues <- choose_pvalues[logi>=1,]

  # 去掉空值
  logi1 <- choose_pvalues$gene_a != ""
  logi2 <- choose_pvalues$gene_b != ""
  logi <- logi1 & logi2
  choose_pvalues <- choose_pvalues[logi,]

  # 同样的条件保留choose_means
  choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

  library(tidyverse)
  meansdf <- choose_means %>% reshape2::melt()
  meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                        CC = meansdf$variable,
                        means = meansdf$value)
  pvalsdf <- choose_pvalues %>% reshape2::melt()
  pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                        CC = pvalsdf$variable,
                        pvals = pvalsdf$value)

  # 合并p值和mean文件
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")

  # dotplot可视化
  summary((filter(pldf,means >0))$means)
  head(pldf)
  pcc =  pldf%>% filter(means >0) %>%
    ggplot(aes(CC.x,interacting_pair.x) )+
    geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
    scale_size_continuous(range = c(1,3))+
    scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 1  )+
    theme_bw()+
    # scale_color_manual(values = rainbow(100))+
    theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8))
  ggplot2::ggsave('out/cpdb_dotplot.png',width = 10,height = 5)


}

