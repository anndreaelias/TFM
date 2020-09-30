#################
#################
#### FIGURES ####
#################
#################

# Required packages
if(!require("xlsx")) (install.packages("xlsx")); library(xlsx)
if(!require("readxl")) (install.packages("readxl")); library(readxl)
if(!require("devtools")) (install.packages("devtools")); library(devtools)
if(!require("ggplot2")) (install.packages("ggplot2")) ; library(ggplot2)  
if(!require("dplyr")) (install.packages("dplyr")) ; library(dplyr)  
if(!require("dendextend")) (install.packages("dendextend")) ; library(dendextend)  
if(!require("RColorBrewer")) (install.packages("RColorBrewer")) ; library(RColorBrewer)  
if(!require("GUniFrac")) (install.packages("GUniFrac")) ; library(GUniFrac)  
if(!require("labdsv")) (install.packages("labdsv")) ; library(labdsv)  
if(!require("corrplot")) (install.packages("corrplot")); library(corrplot)
if(!require("ggsignif")) (install.packages("ggsignif")); library(ggsignif)
if(!require("tidyr")) (install.packages("tidyr")); library(tidyr)
if(!require("wesanderson")) (install.packages("wesanderson")); library(wesanderson)
if(!require("ggbiplot")) (install.packages("ggbiplot")); library(ggbiplot)
if(!require("ggendro")) (install.packages("ggdendro")); library(ggdendro)
if(!require("readlxl")) (install.packages("readxl")); library(readxl)
if(!require("vegan")) (install.packages("vegan")) ; library(vegan)  
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

# Read data
source("read_ddbb.R")

##################################################
# Figure 2.2: Follow-up of IBD samples over time #
##################################################

meta_timepoint <- as.data.frame(read.xlsx("New.metadata.xlsx", sheetIndex = 3))
meta_timepoint <- meta_timepoint[,-1]
# CD plot #
cd <- subset(meta_timepoint, meta_timepoint$Disease=="CD")
cd_rem <- subset(meta_timepoint, meta_timepoint$Disease=="CD" & meta_timepoint$Status=="Remission")
cd_rel <- subset(meta_timepoint, meta_timepoint$Disease=="CD" & meta_timepoint$Status=="Relapse")
#table(cd_rem$Timepoint)
#table(cd_rel$Timepoint)
T1 <- c("TP0", "TP1", "TP2", "TP3", "TP4",
        "TP0", "TP1", "TP2", "TP4",
        "TP0", "TP1", "TP3", "TP4",
        "TP0", "TP2", "TP4",
        "TP0", "TP1", "TP2", "TP3",
        "TP0", "TP1", "TP3",
        "TP0", "TP1", "TP2",
        "TP0", "TP1",
        "TP0")
X1 <- c(0,1,2,3,4,
        0,1,2,4,
        0,1,3,4,
        0,2,4,
        0,1,2,3,
        0,1,3,
        0,1,2,
        0,1,
        0)
Size1 <- c(rep("n.a = 13", 5),
           rep("n.b = 5", 4),
           rep("n.c = 2", 4),
           rep("n.d = 1", 3),
           rep("n.e = 4", 4),
           rep("n.f = 2", 3),
           rep("n.g = 4", 3),
           rep("n.h = 2", 2),
           rep("n.i = 1", 1))
cd.plot <- as.data.frame(cbind(T1, X1, Size1))
cd.plot$X1 <- as.numeric(cd.plot$X1)
#pdf("Figures/CDgrouptimepoints.pdf") 
cd_grouptimepoints <- ggplot(cd.plot, aes(T1, factor(Size1))) + 
  geom_point(colour="grey70", size=2.5) +
  annotate("point", shape=4, x = c(4,3,2,4,3), y = c(2,3,4,4,6), colour = "red", size = 4) +
  labs(x="Time Point", y="ID patient")  +
  theme_bw() +
  geom_segment(aes(x = 1, y = 1, xend = 5, yend = 1), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 2, xend = 5, yend = 2), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 3, xend = 5, yend = 3), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 4, xend = 5, yend = 4), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 5, xend = 4, yend = 5), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 6, xend = 4, yend = 6), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 7, xend = 3, yend = 7), colour="grey50", size=0.3) +
  geom_segment(aes(x = 1, y = 8, xend = 2, yend = 8), colour="grey50", size=0.3) + 
  theme(axis.text = element_text(size=14)) + 
  theme(axis.title.x = element_text(size=16, face = "bold")) +
  theme(axis.title.y = element_text(size=16, face = "bold"))
print(cd_grouptimepoints)
#dev.off()

# UC plot #
uc_rem <- subset(meta_timepoint, meta_timepoint$Disease=="UC" & meta_timepoint$Status=="Remission")
uc_rel <- subset(meta_timepoint, meta_timepoint$Disease=="UC" & meta_timepoint$Status=="Relapse")
table(uc_rem$Timepoint)
table(uc_rel$Timepoint)
T1 <- c("TP0", "TP1", "TP2", "TP3", "TP4",
        "TP0", "TP1", "TP2", "TP4",
        "TP0", "TP2", "TP3", "TP4",
        "TP0", "TP4",
        "TP0", "TP1", "TP2", "TP3",
        "TP0", "TP1", "TP2",
        "TP0", "TP1")
X1 <- c(0,1,2,3,4,
        0,1,2,4,
        0,2,3,4,
        0,4,
        0,1,2,3,
        0,1,2,
        0,1)
Size1 <- c(rep("n.a = 12", 5),
           rep("n.b = 1", 4),
           rep("n.c = 1", 4),
           rep("n.d = 4", 2),
           rep("n.e = 2", 4),
           rep("n.f = 3", 3),
           rep("n.g = 8", 2))
uc.plot <- as.data.frame(cbind(T1, X1, Size1))
uc.plot$X1 <- as.numeric(uc.plot$X1)
#pdf("Figures/UCgrouptimepoints.pdf") 
uc_grouptimepoints <- ggplot(uc.plot, aes(T1, factor(Size1))) + 
  geom_point(colour="grey70", size=2.5) + 
  annotate("point", shape=4, x = c(4,2,2,3,4), y = c(2,3,4,4,4), colour = "red", size = 4) +
  labs(x="Time Point", y="ID patient")  +
  theme_bw() +
  geom_segment(aes(x = 1, y = 1, xend = 5, yend = 1), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 2, xend = 5, yend = 2), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 3, xend = 5, yend = 3), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 4, xend = 5, yend = 4), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 5, xend = 4, yend = 5), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 6, xend = 3, yend = 6), colour="grey50", size=0.3) + 
  geom_segment(aes(x = 1, y = 7, xend = 2, yend = 7), colour="grey50", size=0.3) +
  theme(axis.text = element_text(size=14)) +  
  theme(axis.title.x = element_text(size=16, face = "bold")) +
  theme(axis.title.y = element_text(size=16, face = "bold")) 
print(uc_grouptimepoints)
#dev.off()

##########################################################
# Figure 3.1: Bacterial and fungal loads at initial time #
##########################################################

df <- rbind(clinic_data_cd[,c(1,3,4,5,10,11,12,13,14,15,26)], clinic_data_uc[,c(1,3,4,5,10,11,12,13,14,15,26)], 
            clinic_data_relhc)
df$Disease <- as.factor(df$Disease)
levels(df$Disease) <- c("CD", "HC-CD", "HC-UC", "UC")
#head(df)
## For V4_t0 ##
#pdf("Boxplotgroups_V4.pdf") 
figure3.1.V4 <- ggplot(df, aes(x=Disease, y=V4_t0, fill=Disease)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("CD\n(n=34)", "HC-CD\n(n=29)",
                            "HC-UC\n(n=29)", "UC\n(n=31)")) +
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) +
  geom_signif(comparisons = list(c("UC", "HC-UC")) ,
              map_signif_level=TRUE, 
              test = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Bacterial load (V4)")
print(figure3.1.V4)
#dev.off() 

## For ITS2_t0 ##
#pdf("Boxplotgroups_ITS2.pdf") 
figure3.1.ITS2 <- ggplot(df, aes(x=Disease, y=ITS2_t0, fill=Disease)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("CD\n(n=34)", "HC-CD\n(n=29)",
                            "HC-UC\n(n=29)", "UC\n(n=31)")) +
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) + 
  geom_signif(comparisons = list(c("CD", "UC")),
              map_signif_level=TRUE, 
              test = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Fungal load (ITS2)")
print(figure3.1.ITS2)
#dev.off()

figure_3.1 <- ggarrange(figure3.1.V4, figure3.1.ITS2, 
                                 labels = c("A", "B"),
                                 ncol = 2, nrow = 1)
ggsave(file="Boxplotgroups_V4_ITS2.png", figure_3.1, height=5, width=10, units="in")


################################################################
# Figure 3.2: Calprotectin and hemoglobin levels at final time #
################################################################
# As in CD and UC we have the groups of relapse and remission...
# En este caso quiero separar: CD-RM (tiempo final estado remision), CD-RL (tiempo final estado de recaida)
clinic_data_cd1 <- clinic_data_cd
clinic_data_cd1 <- unite(clinic_data_cd1, newgroup, c(Disease, Group), remove=FALSE)
clinic_data_uc1 <- clinic_data_uc
clinic_data_uc1 <- unite(clinic_data_uc1, newgroup, c(Disease, Group), remove=FALSE)
df2 <- rbind(clinic_data_cd1, clinic_data_uc1)

#pdf("Calprotectin_tfinal.pdf") 
df2_subset <- df2[complete.cases(df2$Calprotectin_tfinal),]
df2_subset$newgroup <- as.factor(df2_subset$newgroup)
figure3.2.calprotectin <- ggplot(df2_subset, aes(x=factor(newgroup, levels=c("CD_Remission", "CD_Relapse", "UC_Remission", "UC_Relapse")), 
                       y=Calprotectin_tfinal, fill=newgroup)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) +
  scale_x_discrete(labels=c("CD-REM\n(n=21)", "CD-RL\n(n=13)",
                            "UC-REM\n(n=13)", "UC-RL\n(n=18)")) +
  geom_signif(comparisons = list(c("CD_Relapse", "CD_Remission"), c("UC_Relapse", "UC_Remission")),
              map_signif_level=TRUE) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Fecal calprotectin levels")
print(figure3.2.calprotectin)
# There are not differences between CD_Remission vs UC_Remission, and neither between UC_Relapse vs CD_Relapse
# There are differences between CD_Relapse and CD_Remission, and also between UC_Relapse and UC_Remission
#dev.off()

#pdf("Boxplot_Hemoglobin_tfinal.pdf") 
df3_subset <- df2[complete.cases(df2$Hemoglobin_tfinal),]
df3_subset$newgroup <- as.factor(df3_subset$newgroup)
figure3.2.hemoglobin <- ggplot(df3_subset, aes(x=factor(newgroup, levels=c("CD_Remission", "CD_Relapse", "UC_Remission", "UC_Relapse")), 
                       y=Hemoglobin_tfinal, fill=newgroup)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) +
  scale_x_discrete(labels=c("CD-REM\n(n=21)", "CD-RL\n(n=13)",
                            "UC-REM\n(n=13)", "UC-RL\n(n=18)")) +
  geom_signif(comparisons = list(c("CD_Relapse", "CD_Remission"), c("UC_Relapse", "UC_Remission")),
              map_signif_level=TRUE) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Hemoglobin levels")
print(figure3.2.hemoglobin)
#dev.off()

figure_3.2 <- ggarrange(figure3.2.calprotectin, figure3.2.hemoglobin, 
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)
ggsave(file="Boxplotgroups_CALPRO_HMG.png", figure_3.2, height=5, width=10, units="in")


########################################################
# Figure 3.3: Bacterial and fungal loads at final time #
########################################################
#pdf("Boxplot_V4_tfinal.pdf") 
df3_subset <- df2[complete.cases(df2$V4_tfinal),]
df3_subset$newgroup <- as.factor(df3_subset$newgroup)
figure3.3.v4 <- ggplot(df3_subset, aes(x=factor(newgroup, levels=c("CD_Remission", "CD_Relapse", "UC_Remission", "UC_Relapse")), 
                       y=V4_tfinal, fill=newgroup)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) +
  scale_x_discrete(labels=c("CD-REM\n(n=21)", "CD-RL\n(n=13)",
                            "UC-REM\n(n=13)", "UC-RL\n(n=18)")) +
  geom_signif(comparisons = list(c("CD_Relapse", "CD_Remission"), c("UC_Relapse", "UC_Remission")),
              map_signif_level=TRUE) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Bacterial load (V4)")
print(figure3.3.v4)
#dev.off()

#pdf("Boxplot_ITS2_tfinal.pdf") 
df4_subset <- df2[complete.cases(df2$ITS2_tfinal),]
df4_subset$newgroup <- as.factor(df4_subset$newgroup)
figure3.3.ITS2 <-ggplot(df4_subset, aes(x=factor(newgroup, levels=c("CD_Remission", "CD_Relapse", "UC_Remission", "UC_Relapse")), 
                       y=ITS2_tfinal, fill=newgroup)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 4)) +
  scale_x_discrete(labels=c("CD-REM\n(n=21)", "CD-RL\n(n=13)",
                            "UC-REM\n(n=13)", "UC-RL\n(n=18)")) +
  geom_signif(comparisons = list(c("CD_Relapse", "CD_Remission"), c("UC_Relapse", "UC_Remission")),
              map_signif_level=TRUE) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Fungal load (ITS2)")
print(figure3.3.ITS2)
#dev.off()

figure_3.3 <- ggarrange(figure3.3.v4, figure3.3.ITS2, 
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)
ggsave(file="Boxplotgroups_V4_ITS2_FINALTIME.png", figure_3.3, height=5, width=10, units="in")


##########################################################################################
# Figure 3.4: PCA to discriminate relapsing and remitting conditions at final time point #
##########################################################################################
# Grups: CD Relapse, CD Remission AT BASELINE
# Variables: ITS2, V4, Calpro, Hemoglobin, CRP, Blood whit cell count
# A continuacio hi ha 4 tipus de grafics iguals per a: CD Initial time, CD Final time,
# UC Initial time, UC Final time. A la memoria nomes he afegit els 2 grafics de Final time.

# Initial time point:
cd_init <- clinic_data_cd
cd_init <- cd_init[,c(1,2,3,4,15,16,18,20,22,24)] 
cd_init <- cd_init[complete.cases(cd_init) , ] 
colnames(cd_init)[c(3,4,6,7,8,9)] <- c("ITS2", "V4", "Calprotectin", "Blood white\ncell count", "CRP", "Hemoglobin")
cd_init$Group <- factor(cd_init$Group, levels=c("Remission", "Relapse"))
pr <- cd_init[,-10]
data.pca <- prcomp(pr[,-c(1,2,5)], center = TRUE,scale. = TRUE)

#pdf("PCA_CD_Initialtime.pdf") 
pca_CD_init <- ggbiplot(data.pca, groups=cd_init$Group, ellipse = TRUE) +
  #geom_point(aes(colour=cd_init$Group, size=cd_init$number_relapses)) + # FOR THE SIZE POINTS
  geom_text(aes(label=cd_init$ID)) +
  theme_bw() + 
  labs(color="Group", size="Number relapses") +
  scale_color_manual(values=wes_palette("Darjeeling2", n = 3)[c(2,3)]) + 
  theme(axis.text = element_text(size=14, face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12)) 
print(pca_CD_init)
#dev.off()

# Final time point:
cd_final <- clinic_data_cd
cd_final <- cd_final[,c(1,2,7,8,15,17,19,21,23,24)] 
cd_final <- cd_final[complete.cases(cd_final) , ] 
colnames(cd_final)[c(3,4,6,7,8,9)] <- c("ITS2", "V4", "Calprotectin", "Blood white\ncell count", "CRP", "Hemoglobin")
cd_final$Group <- factor(cd_final$Group, levels=c("Remission", "Relapse"))
pr <- cd_final[,-10]
data.pca <- prcomp(pr[,-c(1,2,5)], center = TRUE,scale. = TRUE)
df <- as.data.frame(data.pca$x)

#pdf("PCA_CD_Finaltimepoint.pdf", width=8, height=6) 
pca_CD_final <- ggbiplot(data.pca, groups=cd_final$Group, ellipse = TRUE) +
  #geom_point(aes(colour=cd_final$Group, size=cd_final$number_relapses)) + # FOR THE SIZE POINTS
 # geom_text( aes(label= factor(subset(cd_final, cd_final$ID=="TP0.CD.35.0")[,1]))) + # to show all id's
  annotate("text", label=cd_final$ID[10], x=df$PC1[10], y=df$PC2[10]-0.15, size=4.5) +  # add only CD.35
  annotate("text", label=cd_final$ID[8], x=df$PC1[8]-2.3, y=df$PC2[8]-0.25, size=4.5) +  # add only CD.18
  theme_bw() + 
  labs(color="Group", size="Number relapses") +
  ylim(c(-3.5,3)) +
  xlim(c(-2.5,4)) + 
  scale_color_manual(values=wes_palette("Darjeeling2", n = 3)[c(2,3)]) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  theme(title = element_text(size=10, face = "bold")) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12)) +  
  ggtitle("CD patients")

print(pca_CD_final)
#dev.off()


# Grups: UC Relapse, UC Remission
# Variables: ITS2, V4, Calpro, Hemoglobin, CRP, Blood whit cell count
# No es veu tan clara la separacio entre els grups...

# For initial time point:
uc_init_tp <- clinic_data_uc
uc_init_tp <- uc_init_tp[,-c(5,6,7,8,9,10,11,12,13,14,17,19,21,23,25,26)] # Initial time point
uc_init_tp$CRP_t0 <- as.numeric(paste(uc_init_tp$CRP_t0))
uc_init_tp$Hemoglobin_t0 <- as.numeric(paste(uc_init_tp$Hemoglobin_t0))
uc_init_tp <- uc_init_tp[complete.cases(uc_init_tp) , ] # Pasem de 31 pacients a 26
colnames(uc_init_tp)[c(3,4,6,7,8,9)] <- c("ITS2", "V4", "Calprotectin", "Blood white\ncell count", "CRP", "Hemoglobin")
uc_init_tp$Group <- factor(uc_init_tp$Group, levels=c("Remission", "Relapse"))
pr <- uc_init_tp[,-10]
data.pca <- prcomp(pr[,-c(1,2,5)], center = TRUE,scale. = TRUE)

#pdf("PCA_UC_Initialtimepoint.pdf") 
pca_UC_init <- ggbiplot(data.pca, groups=uc_init_tp$Group, ellipse = TRUE) +
  #geom_point(aes(colour=uc_init_tp$Group, size=uc_init_tp$number_relapses)) + # FOR THE SIZE POINTS
  geom_text(aes(label=uc_init_tp$ID)) +
  theme_bw() + 
  labs(color="Group", size="Number relapses") +
  xlim(c(-2, 2.8)) +
  ylim(c(-2.5,2.5)) +   
  scale_color_manual(values=wes_palette("Darjeeling2", n = 3)[c(2,3)]) + 
  theme(axis.text = element_text(size=14, face = "bold"),
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12)) 
print(pca_UC_init)
#dev.off()

# For final time point:
uc_final_tp <- clinic_data_uc
uc_final_tp <- uc_final_tp[,-c(3,4,5,6,9,10,11,12,13,14,16,18,20,22,25,26)] 
uc_final_tp$blood_white_cell_count_tfinal <- as.numeric(paste(uc_final_tp$blood_white_cell_count_tfinal))
uc_final_tp$CRP_tfinal <- as.numeric(paste(uc_final_tp$CRP_tfinal))
uc_final_tp$Hemoglobin_tfinal <- as.numeric(paste(uc_final_tp$Hemoglobin_tfinal))

uc_final_tp <- uc_final_tp[complete.cases(uc_final_tp) , ] # Pasem de 31 pacients a 26
colnames(uc_final_tp)[c(3,4,6,7,8,9)] <- c("ITS2", "V4", "Calprotectin", "Blood white\ncell count", "CRP", "Hemoglobin")
uc_final_tp$Group <- factor(uc_final_tp$Group, levels=c("Remission", "Relapse"))
pr <- uc_final_tp[,-10]
data.pca <- prcomp(pr[,-c(1,2,5)], center = TRUE,scale. = TRUE)
df <- as.data.frame(data.pca$x)

#pdf("PCA_UC_Finaltimepoint.pdf",  width=8, height=6) 
pca_UC_final <- ggbiplot(data.pca, groups=uc_final_tp$Group, ellipse = TRUE) +
#  geom_point(aes(colour=uc_final_tp$Group, size=uc_final_tp$number_relapses)) + # FOR THE SIZE POINTS
 # geom_text(aes(label=uc_final_tp$ID)) + # to show ID's
  annotate("text", label=uc_final_tp$ID[18], x=df$PC1[18]-2.5, y=df$PC2[18]+0.25, size=4.5) +  # add only UC.56
  annotate("text", label=uc_final_tp$ID[11], x=df$PC1[11]-0.1, y=df$PC2[11]+0.45, size=4.5) +  # add only UC.17
  annotate("text", label=uc_final_tp$ID[7], x=df$PC1[7]+0.4, y=df$PC2[7]+0.15, size=4.5) +  # add only UC.53
  theme_bw() + 
  labs(color="Group", size="Number relapses") +
  ylim(c(-3.5,3)) + 
  xlim(c(-2.5,4)) + 
  scale_color_manual(values=wes_palette("Darjeeling2", n = 3)[c(2,3)]) + 
  theme(title = element_text(size=10, face = "bold")) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12)) +  
  ggtitle("UC patients")
    
print(pca_UC_final)
#dev.off()

figure_3.4 <- ggarrange(pca_CD_final, pca_UC_final, 
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)
ggsave(file="PCAS_finaltime.png", figure_3.4, height=5, width=12, units="in")


###################################################
# Figure 3.5: Differences in time for CD patients #
####################################################

## ITS2 in CD Remission 
cd_rem <- subset(clinic_data_cd, clinic_data_cd$Group=="Remission")
head(cd_rem)
var_0 <- cd_rem$ITS2_t0
var_f <- cd_rem$ITS2_tfinal
var <- c(var_0, var_f)
time <- c(rep("Initial", length(var_0)), rep("Final", length(var_f)))
df <- as.data.frame(cbind(var, time))
df$var <- as.numeric(paste(df$var))
df$id <- as.character(rep(cd_rem$ID,2))
df <- df[complete.cases(df),]

#pdf("CD_REMISSION_ITS2.pdf") 
its2_rem_cd <- ggplot(df, aes(factor(time, levels=c("Initial", "Final")),var)) +
  geom_point(aes(color=time)) +
  geom_text(data=df %>% filter(var >= 59272800), aes(label=id), size=4.5) + # to show id's
  geom_line(aes(group = id)) +
  theme_bw() + 
  labs(x="", y="") + 
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(legend.position="none") +
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("ITS2 in CD Remission")

print(its2_rem_cd)
#dev.off()

## Calprotectin in CD Relapse

cd_rel <- subset(clinic_data_cd, clinic_data_cd$Group=="Relapse")
head(cd_rel)
var_0 <- cd_rel$Calprotectin_t0
var_f <- cd_rel$Calprotectin_tfinal
var <- c(var_0, var_f)
time <- c(rep("Initial", length(var_0)), rep("Final", length(var_f)))
df <- as.data.frame(cbind(var, time))
df$id <- rep(cd_rel$ID,2)
df$var <- as.numeric(paste(df$var))
df <- df[complete.cases(df),]

#pdf("CD_RELAPSE_CALPRO.pdf") 
calpro_rel_cd <- ggplot(df, aes(factor(time, levels=c("Initial", "Final")),var)) +
  geom_point(aes(color=time)) +
  geom_text(data=df %>% filter(var >= 2000), aes(label=id), size=4.5) + # to show id's
  geom_line(aes(group = id)) +
  theme_bw() + 
  labs(x="", y="") + 
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(axis.text = element_text(size=11, face = "bold")) +
  theme(legend.position="none") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Calprotectin in CD Relapse")

print(calpro_rel_cd)
#dev.off()

differences_time_cd <- ggarrange(its2_rem_cd, calpro_rel_cd, 
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1)
ggsave(file="CD_ITS2_CALPRO_DIFF.png", differences_time_cd, height=4, width=8, units="in")

###################################################
# Figure 3.6: Differences in time for UC patients #
###################################################

##### ITS2 in UC Remission

uc_rem <- subset(clinic_data_uc, clinic_data_uc$Group=="Remission")
head(uc_rem)
var_0 <- uc_rem$ITS2_t0
var_f <- uc_rem$ITS2_tfinal
var <- c(var_0, var_f)
time <- c(rep("Initial", length(var_0)), rep("Final", length(var_f)))
df <- as.data.frame(cbind(var, time))
df$id <- rep(uc_rem$ID,2)
df$var <- as.numeric(paste(df$var))
df <- df[complete.cases(df),]

#pdf("UC_REMISSION_ITS2.pdf") 
its2_rem_uc <- ggplot(df, aes(factor(time, levels=c("Initial", "Final")),var)) +
  geom_point(aes(color=time)) +
  geom_text(data=df %>% filter(var >= 150000000), aes(label=id), size=4.5) + # to show id's
  geom_line(aes(group = id)) +
  theme_bw() + 
  labs(x="", y="") + 
  theme(legend.position="none") +
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(axis.text = element_text(size=11, face = "bold")) +
  theme(legend.position="none") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("ITS2 in UC Remission")
print(its2_rem_uc)
#dev.off()

##### Hemoglobin in UC Relapse

uc_rel <- subset(clinic_data_uc, clinic_data_uc$Group=="Relapse")
head(uc_rel)
var_0 <- uc_rel$Hemoglobin_t0
var_f <- uc_rel$Hemoglobin_tfinal
var <- c(var_0, var_f)
id <- rep(uc_rel$ID,2)
time <- c(rep("Initial", length(var_0)), rep("Final", length(var_f)))
df <- as.data.frame(cbind(id, var, time))
df <- df[complete.cases(df),]
df$var <- as.numeric(paste(df$var))

#pdf("UC_RELAPSE_HEMOGLOBIN.pdf") 
hemogl_rel_uc <- ggplot(df, aes(factor(time, levels=c("Initial", "Final")),var)) +
  geom_point(aes(color=time)) +
  geom_line(aes(group = id)) +
  theme_bw() + 
  labs(x="", y="") + 
  theme(legend.position="none") +
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(axis.text = element_text(size=11, face = "bold")) +
  theme(legend.position="none") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("Hemoglobin in UC Relapse")
  
print(hemogl_rel_uc)
#dev.off()

##### CRP in UC Relapse

uc_rel <- subset(clinic_data_uc, clinic_data_uc$Group=="Relapse")
head(uc_rel)
var_0 <- uc_rel$CRP_t0
var_f <- uc_rel$CRP_tfinal
var <- c(var_0, var_f)
time <- c(rep("Initial", length(var_0)), rep("Final", length(var_f)))
df <- as.data.frame(cbind(id, var, time))
df$id <- rep(uc_rel$ID,2)
df$var <- as.numeric(paste(df$var))
df <- df[complete.cases(df),]

#pdf("UC_RELAPSE_CRP.pdf") 
crp_rel_uc <- ggplot(df, aes(factor(time, levels=c("Initial", "Final")),var)) +
  geom_point(aes(color=time)) +
  geom_text(data=df %>% filter(var >= 20), aes(label=id)) + # to show id's
  geom_line(aes(group = id)) +
  theme_bw() + 
  labs(x="", y="") + 
  theme(legend.position="none") +
  theme(axis.text = element_text(size=11, face = "bold")) + 
  theme(axis.text = element_text(size=11, face = "bold")) +
  theme(legend.position="none") + 
  theme(title = element_text(size=10, face = "bold")) +
  ggtitle("CRP in UC Relapse")

print(crp_rel_uc)
#dev.off()

differences_time_uc <- ggarrange(its2_rem_uc, hemogl_rel_uc, crp_rel_uc, 
                                 labels = c("A", "B", "C"),
                                 ncol = 3, nrow = 1)
ggsave(file="UC_ITS2_HMG_CRP_DIFF.png", differences_time_uc, height=4, width=9, units="in")


#########################################
# Figure 3.7: Quality plots in each run #
#########################################

# These plots are an output from QIIME2 software.

#############################################
# Figure 3.8: Batch effect controls cluster #
#############################################

test <- as.dist(distance, diag = TRUE)
hc <- hclust(test, method="average")
dendr <- dendro_data(hc, type="rectangle") 

cols <- c(rep("plum",2), rep("purple",2), rep("salmon4",2), rep("chocolate1",2),
          rep("darkgoldenrod1",2), rep("lightblue1",2), rep("lightgoldenrod4",2), rep("greenyellow",2),
          rep("hotpink2",2), rep("plum2",2), rep("saddlebrown",2), rep("forestgreen",2), rep("red",2),
          rep("moccasin",2), rep("yellowgreen",2),
          rep("skyblue2",2), rep("deeppink3",2), rep("olivedrab1",2), rep("steelblue",2),
          rep("gold",2), rep("cyan",2), rep("darkblue",2))

long <- rep(-0.19, nrow(distance))
lat <- 1:44
points <- as.data.frame(cbind(long, lat, cols))

#pdf("Cluster_dendogram.pdf") 
cluster <- ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x=x, y=y, label=label, hjust=0), size=4) +
  geom_point(data = points, aes(x = as.numeric(paste(lat)), y = as.numeric(paste(long))), color=cols, size=4) + 
  scale_color_manual(values=cols) + 
  coord_flip() + scale_y_reverse(expand=c(0.3, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())
print(cluster)
#dev.off()

###############################################
# Figure 3.9:  Batch effect controls boxplot # 
###############################################
#Preparing data
# Run5 vs Run8
df1 <- distance[1:10, c("FI.01.run8","FI.02.run8", "FI.03.run8","FI.04.run8","FI.05.run8","FI.06.run8","FI.07.run8","FI.08.run8", "FI.09.run8","FI.10.run8")]
x1 <- diag(as.matrix(df1[,]))
boxplot(diag(as.matrix(df1[,])), ylim=c(0,0.16))
summary(diag(as.matrix(df1[,])))

# Run5 vs Run7
df2 <- distance[11:14, c("P208.F.run7", "P209.F.run7", "P210.F.run7", "P211.F.run7")]
boxplot(diag(as.matrix(df2[,])), ylim=c(0,0.16))
x2  <- diag(as.matrix(df2[,]))
summary(diag(as.matrix(df2[,])))

# Run7 vs Run9
df3 <- distance[15:18, c("158.12_20.run9", "158.16_20.run9", "210.12_20.run9", "210.6_12.run9")]
boxplot(diag(as.matrix(df3[,])), ylim=c(0,0.16))
x3  <- diag(as.matrix(df3[,]))
summary(diag(as.matrix(df3[,])))

# Run8 vs Run9
df4 <- distance[33:36, c("UC.51-4.run9", "UC.52-0.run9", "UC.52-1.run9","UC.53-0.run9")]
boxplot(diag(as.matrix(df4[,])), ylim=c(0,0.16))
x4  <- diag(as.matrix(df4[,]))
summary(diag(as.matrix(df4[,])))

cntrls <- c(x1,x2,x3,x4)
ntcntrls <- c(df1[lower.tri(df1, diag = FALSE)], df2[lower.tri(df2, diag = FALSE)], df3[lower.tri(df3, diag = FALSE)], df4[lower.tri(df4, diag = FALSE)])
group <- c(rep("Same HC samples", length(cntrls)), rep("Different HC samples", length(ntcntrls)))
value <- c(cntrls, ntcntrls)
df <- as.data.frame(cbind(group, value))
df$value <- as.numeric(paste(df$value))

#pdf("Boxplot_ControlRuns.pdf", width=8, height=6) 
boxplot_controls <- df %>%
  ggplot(aes(x=factor(group, levels=c("Same HC samples", "Different HC samples", "All samples")), y=value, fill=group)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1) +
  theme_bw() +
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 3)) +
  theme(legend.position="none", plot.title = element_text(size=11)) +
  labs(x="", y="") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
print(boxplot_controls)
#dev.off()

####################################################### 
####################################################### 
# For the following Figures: compositional_analysis.R #
#######################################################
#######################################################
