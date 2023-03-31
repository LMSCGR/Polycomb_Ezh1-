# To generate Figure 3d

# Load the libraries
library(ggplot2)
library(tidyverse)

setwd("path to: /Figure3d")

writeLines(capture.output(sessionInfo()),"sessionInfo.txt")

# import the data and clean the data
newdata  <- read.delim("data_figure3c.txt",sep = "\t", head = TRUE)

subMat <- filter(newdata, p.value < 0.05 )
subMat <- subMat[grep(pattern ="^Mir[0-9]|^Gm[0-9]|^Vmn1r[0-9]|Ppip5k1|Mirlet7f-1|[0-9]Rik$|Kcnd3os|Olfr1532-ps1|Ifna5|Speer4d|Cabp7|Spink14",
                      rownames(subMat),invert=TRUE),]

subMat$type <- c(" ")
subMat$type[subMat$type == " "] <- "A"
subMat$type[subMat$log2FC_A_ko_vs_wt > 0.58  & subMat$log2FC_R_ko_vs_wt > 0.38] <- "B"
subMat$type[subMat$log2FC_A_ko_vs_wt > 0.58  & subMat$log2FC_R_ko_vs_wt < -0.38] <- "C"
subMat$type[subMat$log2FC_A_ko_vs_wt < -0.58 & subMat$log2FC_R_ko_vs_wt > 0.38] <- "D"
subMat$type[subMat$log2FC_A_ko_vs_wt < -0.58 & subMat$log2FC_R_ko_vs_wt < -0.38] <- "E"
subMat$type[subMat$log2FC_A_ko_vs_wt < 0.58 & subMat$log2FC_A_ko_vs_wt > -0.58 & subMat$log2FC_R_ko_vs_wt > 0.38] <- "F"
subMat$type[subMat$log2FC_A_ko_vs_wt < 0.58 & subMat$log2FC_A_ko_vs_wt > -0.58 & subMat$log2FC_R_ko_vs_wt < -0.38] <- "G"
subMat$type[subMat$log2FC_R_ko_vs_wt > -0.38 & subMat$log2FC_R_ko_vs_wt < 0.38 & subMat$log2FC_A_ko_vs_wt < -0.58] <- "H"

subMat_g1 <- subMat[subMat$type=="A",]
subMat_g2 <- subMat[subMat$type=="B",]
subMat_g3 <- subMat[subMat$type=="C",]
subMat_g4 <- subMat[subMat$type=="D",]
subMat_g5 <- subMat[subMat$type=="E",]
subMat_g6 <- subMat[subMat$type=="F",]
subMat_g7 <- subMat[subMat$type=="G",]
subMat_g8 <- subMat[subMat$type=="H",]

ggplot() +
  geom_point(data=subMat_g1,aes(log2FC_A_ko_vs_wt,log2FC_R_ko_vs_wt), color="#999999",size = 0.5) +
  geom_point(data=subMat_g2,aes(log2FC_A_ko_vs_wt, log2FC_R_ko_vs_wt), color="#FF4040",size = 1) +
  geom_point(data=subMat_g3,aes(log2FC_A_ko_vs_wt, log2FC_R_ko_vs_wt), color="#0000FF",size = 1) +
  geom_point(data=subMat_g6,aes(log2FC_A_ko_vs_wt, log2FC_R_ko_vs_wt), color="#999999",size = 1) +
  geom_point(data=subMat_g7,aes(log2FC_A_ko_vs_wt, log2FC_R_ko_vs_wt), color="#999999",size = 1) +
  geom_jitter() +
  labs(
    x = paste0("Promoter Accessibility: ATAC-seq (Ezh1KO/WT, log2FC)"),
    y = paste0("Expression:RAN-seq (Ezh1KO/WT, log2FC)")) +
  theme_bw() +
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size= 14),
        legend.background = element_rect(fill = "transparent", colour = "transparent")) +
  geom_hline(yintercept= 0.38, linetype = 2, color="#999999" )+
  geom_hline(yintercept=-0.38, linetype = 2, color="#999999")+
  geom_vline(xintercept = 0.58 , linetype= 2, color="#999999")+
  geom_vline(xintercept =-0.58 , linetype= 2, color="#999999")

