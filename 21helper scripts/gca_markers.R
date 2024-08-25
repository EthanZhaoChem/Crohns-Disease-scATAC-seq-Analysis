# make sure this script is shared across three folders: atac_gca2024, final_GCArna, final_GCAatac
# original file is in final_GCArna, other folders have soft links
library(ggsci)
library(ggpubr)
library(colorspace)
library(paletteer) 
library(wesanderson)
library(ggrepel)
library(ggrastr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(scales)


gca.rna.unionHeatmap.markers <- list('Stem'=c('LGR5'),
                                    'TA'=c('TUBA1B'),
                                    'Early Enterocyte'=c('OLFM4'),
                                    'Enterocyte'=c('FABP6'),
                                    'CD Enterocyte1'=c('APOA4'),
                                    'CD Enterocyte2'=c('NOS2'),
                                    'Early Colonocyte'=c('MECOM'),
                                    'Colonocyte'=c('CA2', 'AQP8'),
                                    'CD Colonocyte1'=c('PIGR'),
                                    'CD Colonocyte2'=c('CEACAM7'),
                                    'Ileum Goblet'=c('MUC2'),
                                    'Colon Goblet'=c('ITLN1'),
                                    'BEST4'=c('BEST4'),
                                    'Tuft'=c('FYB1'),
                                    'EEC'=c('CHGA'),
                                    'Paneth'=c('DEFA5'),
                                    'Proliferating T'=c('TUBB'),
                                    'CD4 Tcm'=c('ICOS'),             
                                    'Treg'=c('CTLA4'),
                                    'CD103- CD4'=c('IL7R'),
                                    'CD103+ CD4'=c('RORA'),
                                    'CD103+ CD8'=c('ITGAE', 'CCL5'),
                                    'KLRG1+ CD8'=c('GZMK', 'KLRG1'), 
                                    'gdT'=c('TRDC'),
                                    'NK'=c('KLRF1'),
                                    'ILCs'=c('AREG', 'AFF3'),
                                    'GC B'=c('BCL6', 'BACH2'),
                                    'Naive B'=c('IGHD'),
                                    'Memory B'=c('CD27', 'MS4A1'),
                                    'IgA plasma'=c('IGHA1', 'XBP1'),         
                                    'IgG plasma'=c('IGHG1'),
                                    'Mast'=c('CPA3'),
                                    'MAIT'=c('SLC4A10'),
                                    'Monocyte'=c('S100A9'),
                                    'Macrophage'=c('C1QA'),
                                    'DC'=c('HLA-DRA'),
                                    'Neutrophil'=c('FCGR3B'),
                                    'Fibroblast1'=c("ADAMDEC1"),
                                    'Fibroblast2'=c("IGFBP5"),
                                    'Fibroblast3'=c('NRG1'),                                                         
                                    'Fibroblast4'=c("DCN"),
                                    'Myofibroblast'=c("ACTA2"), 
                                    'Pericyte'=c("NOTCH3"),
                                    'Contractile Pericyte'=c("MUSTN1"), 
                                    'Arterial'=c("FLT1"), 
                                    'Venous'=c("FAM155A"),
                                    'Lymphatic endothelium'=c("LYVE1"),
                                    'Glial'=c("S100B"))


gca_rna.lineage.markers <- list(
  'epithelial' = c('EPCAM', 'PHGR1','FABP1', 'KRT8'),
  'tcell' = c('CD3D', 'CCL5', 'IL7R', 'STAT4'),
  'bcell' = c('CD79A', 'MS4A1'),
  'plasma cell' = c('XBP1', 'SLAMF7', 'IGHA1', 'IGHG1'),
  'myeloid' = c('C1QA', 'LYZ', 'CPA3', 'TYROBP'),
  'others' = c('IGFBP7', 'COL3A1', 'PECAM1', 'S100B')
)


gca_rna.epithelial.markers <- list(
  "Stem" = c('LGR5', 'OLFM4', 'SMOC2'), 
  "TA" = c('MKI67', 'TOP2A', 'TUBA1B'), 
  "EC1-1" = c('ADH1C', 'SI', 'GSTA2'),
  "EC1-2" = c('SLC15A1', 'APOA4', 'CUBN'),
  "EC2-1" = c('CD24', 'FABP5'),
  "EC2-2" = c('CEACAM5', 'CA2', 'CEACAM7', 'AQP8'),
  "Goblet1" = c('MUC2', 'TFF3', 'CLCA1', 'SPINK4', 'FER1L6'),
  "Goblet2" = c('ITLN1',  'KLK1', 'BEST2'),
  "M-like" = c('CCL20'),
  "BEST4" = c('BEST4', 'CA7', 'NOTCH2', 'SPIB'),
  "Paneth" = c('ITLN2', 'DEFA5', 'REG3A'),
  "EEC" = c('CHGA', 'KCNB2', 'RIMBP2'),
  "Tuft" = c('POU2F3', 'FYB1')
)

gca_rna.tcell.markers <- list(
  "CD4 Tcm" = c('CD3D', 'CD4', 'LEF1', 'TCF7', 'CCR7'), 
  "Treg" = c('CTLA4', 'FOXP3', 'IL2RA'), 
  "CD103- CD4 Trm" = c('IL7R'), 
  "CD103+ CD4 Trm" = c('ITGAE', 'IL17A', 'IL26'), 
  "CD103+ CD8 Trm" = c('CD8A','CD8B'), 
  "KLRG1+ CD8 Trm" = c('KLRG1', 'GZMK','IFNG', 'HLA-DRB1', 'HLA-DRA'), 
  "gdT" = c('TRDC', 'GNLY', 'GZMA', 'ENTPD1'), 
  "MAIT" = c('SLC4A10', 'NCR3'),
  "NK T" = c('FGFBP2', 'GZMB', 'CX3CR1'),
  "NK" = c('NCR1', 'KLRF1'), 
  "ILCs" = c('PRKG1', 'PCDH9', 'AFF3', 'AREG','IL1R1', 'IL23R', 'KIT') 
)


gca_rna.bcell.markers <- list(
  "GC B" = c('BCL6'),
  "Naive B" = c('CD19', 'MS4A1', 'IGHD'),
  "Memory B" = c('CD27', 'CD83'),
  "IgA plasma" = c('SDC1', 'SLAMF7', 'IGHA1'),
  "IgG plasma" = c('IGHG1')
)

gca_rna.myeloid.markers <- list("Monocyte" = c('VCAN', 'FCN1', 'EREG', 'S100A8', 'S100A4'),
                            "Macrophage" = c('CD14', 'CD163', 'MMP12', 'C1QA'),
                            "cDC1" = c('CLEC9A'),
                            "cDC2" = c('CD1C'),
                            "Lymphoid DC" = c('LAMP3'),
                            "Mast" = c('CPA3', 'KIT', 'TPSB2'),
                            "Neutrophil" = c('FCGR3B')
)

gca_rna.stromal.markers <- list("Stromal-1" = c('COL1A1', 'ADAMDEC1', 'CCL11', 'CCL13'), 
                           "Stromal-2" = c('NRG1', 'NPY', 'PTGS1'),  
                           "Stromal-3" = c('SOX6', 'COL4A6'), 
                           "Myofibroblast" = c('ACTA2', 'TAGLN'), 
                           "Arterial" = c('PECAM1', 'HEY1', 'EFNB2'), 
                           "Venous" = c('ACKR1', 'VWF'), 
                           "Pericyte" = c('NOTCH3', 'MCAM', 'RGS5'), 
                           "Contractile pericyte" = c('PLN', 'RERGL', 'KCNAB1'), 
                           "Smooth muscle" = c('DES', 'CNN1'), 
                           "Lymphatic endothelium" = c('PROX1', 'LYVE1', 'CCL21'), 
                           "Glial" = c('S100B', 'NRXN1'))


gca_rna.immune.markers <- unique(c(unlist(gca_rna.tcell.markers,
                                          gca_rna.bcell.markers,
                                          gca_rna.myeloid.markers)))


























