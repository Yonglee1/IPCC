#' Using cellchat to analysis.1
#'
#' @param input your scRNA
#' @param species mouse or human
#'
#' @return pathwayname
#' @export
#'
#' @examples ipcc_cellchat1(input = scRNA,specis="mouse")
ipcc_cellchat1 <- function(input,species){

  scRNA <- input
  library(CellChat)
  library(ggalluvial)
  library(svglite)
  library(Seurat)
  library(ggplot2)
  library(CellChat)
  library(tidyverse)
  library(ggalluvial)

  options(stringsAsFactors = FALSE)

  # CellChat要求输入标准化后的表达数据
  data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  identity <- subset(scRNA@meta.data, select = "celltype")


  ##创建cellchat对象
  cellchat <- createCellChat(object = data.input)
  cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  groupSize <- as.numeric(table(cellchat@idents)) # 后面有用

  ##设置参考数据库
  # 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
  if(species=="mouse"){CellChatDB <- CellChatDB.mouse  }
  else{CellChatDB <-CellChatDB.human}
  # 使用"Secreted Signaling"用于细胞通讯分析
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  # 将数据库传递给cellchat对象
  cellchat@DB <- CellChatDB.use

  ##配体-受体分析
  # 提取数据库支持的数据子集
  cellchat <- subsetData(cellchat)#cellchat0 <- cellchat
  # 识别过表达基因,多了var.features
  cellchat <- identifyOverExpressedGenes(cellchat)#cellchat1 <- cellchat
  # 识别配体-受体对,多了LR
  cellchat <- identifyOverExpressedInteractions(cellchat)#cellchat2 <- cellchat
  # 将配体、受体投射到PPI网络，多了dataproject
  cellchat <- projectData(cellchat, PPI.mouse)#cellchat3 <- cellchat

  ##推测细胞通讯网络
  cellchat <- computeCommunProb(cellchat)#cellchat4 <- cellchat,多了net
  cellchat <- computeCommunProbPathway(cellchat)#cellchat5 <- cellchat,多了netp
  cellchat <- aggregateNet(cellchat)#cellchat6 <- cellchatm,多了net的count和weight

  #细胞通讯网络系统分析及可视化
  vertex.receiver = c(3, 6)          #指定靶细胞的索引
  cellchat@netP$pathways             #查看富集到的信号通路***

}


#' Using cellchat to analysis.2
#'
#' @param input your scRNA
#' @param species mouse or human
#' @param path input pathway name from console
#'
#' @return figures
#' @export
#'
#' @examples ipcc_cellchat2(input=scRNA,specis="mouese",path="GALECTIN")
ipcc_cellchat2 <- function(input,species,path){

  scRNA <- input
  library(CellChat)
  library(ggalluvial)
  library(svglite)
  library(Seurat)
  library(ggplot2)
  library(CellChat)
  library(tidyverse)
  library(ggalluvial)

  options(stringsAsFactors = FALSE)

  # CellChat要求输入标准化后的表达数据
  data.input <- GetAssayData(scRNA, assay = "RNA", slot = "data")
  identity <- subset(scRNA@meta.data, select = "celltype")


  ##创建cellchat对象
  cellchat <- createCellChat(object = data.input)
  cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  groupSize <- as.numeric(table(cellchat@idents)) # 后面有用

  ##设置参考数据库
  # 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
  if(species=="mouse"){CellChatDB <- CellChatDB.mouse  }
  else{CellChatDB <-CellChatDB.human}

  # 使用"Secreted Signaling"用于细胞通讯分析
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  # 将数据库传递给cellchat对象
  cellchat@DB <- CellChatDB.use

  ##配体-受体分析
  # 提取数据库支持的数据子集
  cellchat <- subsetData(cellchat)#cellchat0 <- cellchat
  # 识别过表达基因,多了var.features
  cellchat <- identifyOverExpressedGenes(cellchat)#cellchat1 <- cellchat
  # 识别配体-受体对,多了LR
  cellchat <- identifyOverExpressedInteractions(cellchat)#cellchat2 <- cellchat
  # 将配体、受体投射到PPI网络，多了dataproject
  cellchat <- projectData(cellchat, PPI.mouse)#cellchat3 <- cellchat

  ##推测细胞通讯网络
  cellchat <- computeCommunProb(cellchat)#cellchat4 <- cellchat,多了net
  cellchat <- computeCommunProbPathway(cellchat)#cellchat5 <- cellchat,多了netp
  cellchat <- aggregateNet(cellchat)#cellchat6 <- cellchatm,多了net的count和weight

  #细胞通讯网络系统分析及可视化
  vertex.receiver = c(3, 6)          #指定靶细胞的索引

  pathways.show <- path            #指定需要展示的通路**

  #dotplot
  p=netVisual_bubble(cellchat,angle = -45,hjust = -0.1,vjust = 0.8)
  ggsave("out/cellchat_dotplot.png", p, width = 9, height = 6)
  # Hierarchy plot
  png(filename = "out/cellchat_sig_pathway_hierarchy.png", width = 1000, height = 650)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = vertex.receiver, vertex.weight = groupSize)
  dev.off()
  # Circle plot
  png(filename = "out/cellchat_sig_pathway_cricle.png", width = 650, height = 600)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
  dev.off()
  #chord
  png(filename = "out/cellchat_sig_pathway_chord.png", width = 650, height = 600)
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", vertex.receiver = vertex.receiver, vertex.weight = groupSize)
  dev.off()

  # 计算配体-受体对信号网络的贡献度
  png(filename = "out/cellchat_sig_pathway_L-R.png", width = 800, height = 600)
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  dev.off()
  # 分析细胞在信号网络中角色
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  png(filename = "out/cellchat_sig_pathway_role.png", width = 400, height = 300)
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show)
  dev.off()
  #使用散点图在 2D 空间中可视化占主导地位的发射器（源）和接收器（目标）。
  png(filename = "out/cellchat_sig_pathway_scatter.png", width = 800, height = 600)
  netAnalysis_signalingRole_scatter(cellchat)
  dev.off()

  ##细胞通讯模式和信号网络
  nPatterns = 4   #默认为5
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  # river plot
  p = netAnalysis_river(cellchat, pattern = "outgoing")
  ggsave("out/cellchat_com_pattern_outgoing_river.png", p, width = 12, height = 6)
  # dot plot
  p = netAnalysis_dot(cellchat, pattern = "outgoing")
  ggsave("out/cellchat_com_pattern_outgoing_dot.png", p, width = 9, height = 6)

  nPatterns = 4
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  # river plot
  p = netAnalysis_river(cellchat, pattern = "incoming")
  ggsave("out/cellchat_com_pattern_incoming_river.png", p, width = 12, height = 6)
  # dot plot
  p = netAnalysis_dot(cellchat, pattern = "incoming")
  ggsave("out/cellchat_com_pattern_incoming_dot.png", p, width = 9, height = 6)


}
