

#' Using nichenet to analysis
#'
#' @param input you scRNA
#' @param species mouse or human
#' @param control the name of your control group
#' @param case the name of your case group
#'
#' @return figures and txt
#' @export
#'
#' @examples ipcc_nichenet (input=scRNA,species="mouse",control="SS",case="LCMV")
ipcc_nichenet <- function(input,species,control,case){
  scRNA <- input
  library(nichenetr)
  library(Seurat)
  library(tidyverse)

  ##读入nichenet先验数据
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
  # weighted_networks列表包含两个数据框，lr_sig是配体-受体权重信号网络，gr是配体-靶基因权重调控网络

  scRNA <- UpdateSeuratObject(scRNA)

  #aggregate是处理条件，SS相当于control，LCMV相当于case。

  Idents(scRNA) <- "celltype"
  #数nichenet_seuratobj_aggregate，可以一步完成seurat对象的配体调控网络分析。
  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNA,
                                                 top_n_ligands = 20,
                                                 receiver = "CD8 T",
                                                 sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
                                                 condition_colname = "aggregate",
                                                 condition_oi = case,
                                                 condition_reference = control,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 lr_network = lr_network,
                                                 weighted_networks = weighted_networks,
                                                 organism = species)
  # top_n_ligands参数指定用于后续分析的高活性配体的数量

  ## 查看配体活性分析结果
  # 主要参考pearson指标，bona_fide_ligand=True代表有文献报道的配体-受体，
  # bona_fide_ligand=False代表PPI预测未经实验证实的配体-受体。
  x <- nichenet_output$ligand_activities
  write.csv(x, "out/ncn_ligand_activities.csv", row.names = F)
  # 查看top20 ligands在各个细胞亚群中表达情况
  p = DotPlot(scRNA, features = nichenet_output$top_ligands, cols = "RdYlBu") + RotatedAxis()
  ggsave("out/ncn_top20_ligands.png", p, width = 12, height = 6)

  # 按"aggregate"的分类对比配体的表达情况
  p = DotPlot(scRNA, features = nichenet_output$top_ligands, split.by = "aggregate") + RotatedAxis()
  ggsave("out/ncn_top20_ligands_compare.png", p, width = 12, height = 8)
  # 用小提琴图对比配体的表达情况
  p = VlnPlot(scRNA, features = nichenet_output$top_ligands,
              split.by = "aggregate",  pt.size = 0, combine = T)
  ggsave("out/ncn_VlnPlot_ligands_compare.png", p, width = 12, height = 10)

  ## 查看配体调控靶基因
  p = nichenet_output$ligand_target_heatmap
  ggsave("out/ncn_ncn_Heatmap_ligand-target.png", p, width = 12, height = 6)
  # 查看top配体调控的靶基因及其评分
  x <- nichenet_output$ligand_target_matrix
  #x2 <- nichenet_output$ligand_target_df
  write.csv(x, "out/ncn_ligand_target.csv", row.names = T)
  # 查看被配体调控靶基因的表达情况
  p = DotPlot(scRNA %>% subset(idents = "CD8 T"),
              features = nichenet_output$top_targets,
              split.by = "aggregate") + RotatedAxis()
  ggsave("out/ncn_Targets_Expression_dotplot.png", p, width = 12, height = 6)
  p = VlnPlot(scRNA %>% subset(idents = "CD8 T"), features = nichenet_output$top_targets,
              split.by = "aggregate", pt.size = 0, combine = T, ncol = 8)
  ggsave("out/ncn_Targets_Expression_vlnplot.png", p, width = 12, height = 8)

  ## 查看受体情况
  # 查看配体-受体互作
  p = nichenet_output$ligand_receptor_heatmap
  ggsave("out/ncn_Heatmap_ligand-receptor.png", p, width = 12, height = 6)
  x <- nichenet_output$ligand_receptor_matrix
  #x <- nichenet_output$ligand_receptor_df
  write.csv(x, "out/ncn_ligand_receptor.csv", row.names = F)
  # 查看受体表达情况
  p = DotPlot(scRNA %>% subset(idents = "CD8 T"),
              features = nichenet_output$top_receptors,
              split.by = "aggregate") + RotatedAxis()
  ggsave("out/ncn_Receptors_Expression_dotplot.png", p, width = 12, height = 6)
  p = VlnPlot(scRNA %>% subset(idents = "CD8 T"), features = nichenet_output$top_receptors,
              split.by = "aggregate", pt.size = 0, combine = T, ncol = 8)
  ggsave("out/ncn_Receptors_Expression_vlnplot.png", p, width = 12, height = 8)
  # 有文献报道的配体-受体
  # Show ‘bona fide’ ligand-receptor links that are described in the literature and not predicted based on PPI
  p = nichenet_output$ligand_receptor_heatmap_bonafide
  ggsave("out/ncn_Heatmap_ligand-receptor_bonafide.png", p, width = 8, height = 4)
  x <- nichenet_output$ligand_receptor_matrix_bonafide
  #x <- nichenet_output$ligand_receptor_df_bonafide
  write.csv(x, "out/ncn_ligand_receptor_bonafide.csv", row.names = F)

}
