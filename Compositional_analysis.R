######################
######################
# Taxonomic analysis #
######################
######################

##############################
#### 1) Required packages ####
##############################
if(!require("readlxl")) (install.packages("readxl")); library(readxl)
if(!require("xlsx")) (install.packages("xlsx")); library(xlsx)
if(!require("ggplot2")) (install.packages("ggplot2")) ; library(ggplot2)  
if(!require("dplyr")) (install.packages("dplyr")) ; library(dplyr)  
if(!require("vegan")) (install.packages("vegan")) ; library(vegan)  
if(!require("ggsignif")) (install.packages("ggsignif")); library(ggsignif)
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")


########################
#### 2) Source data ####
########################

source("read_ddbb.R") # for the metadata

# Otu table and taxonomy:
otu <- as.data.frame(read_excel("table-L6.xlsx")) # otu table in terms of absolute frequency L6 - Genus level
head(otu)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax <- rbind(unlist(strsplit(rownames(otu)[1], "[;]"))) # to create taxonomy table
for(i in 2:nrow(otu)){
  tax <- rbind(tax, unlist(strsplit(rownames(otu)[i], "[;]")))
}
taxonomy_L6 <- as.data.frame(tax)
colnames(taxonomy_L6) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(otu) <- 1:nrow(otu)
dim(otu)
rownames(otu) <- paste0("otu_0", rownames(otu))
rownames(taxonomy_L6) <- paste0("otu_0", rownames(taxonomy_L6))
matches <- rownames(metadata)
otus <- otu[, as.character(matches)]
rownames(otus) <- rownames(taxonomy_L6)

###########################################
#### 3) Create a phyloseq environment: ####
###########################################

OTU.physeq = otu_table(as.matrix(otus), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(taxonomy_L6))
meta.physeq = sample_data(metadata) # rownames de metadata han de ser el ID de los pacientes
physeq.alpha = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq.alpha

# The following code is to create our palette of colors in all of the following barplots for compositional data.
phylum_colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                   "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", 
                   "#5E738F","#D1A33D", "#8A7C64", "#599861", "cyan", "yellow")

######## The following plots are in level L6: Genus

############################################################ 
# Figure 3.10: Alpha-diversity indexes in different groups #
############################################################

newmeta1 <- subset(metadata, ((metadata$t0_tfinal=="Initial" & (metadata$Group=="Healthy"| metadata$Group=="HR_UC" | metadata$Group=="HR_CD")) |
                                metadata$t0_tfinal=="Final" & (metadata$Group=="CD" | metadata$Group=="CD_RL" | 
                                                                 metadata$Group=="UC" | metadata$Group=="UC_RL"
                                )))

newmeta1$Group <- factor(newmeta1$Group, levels=c("CD", "CD_RL", "HR_CD", "Healthy", "HR_UC", "UC_RL", "UC"))
meta.physeq1 = sample_data(newmeta1) # rownames de metadata han de ser el ID de los pacientes
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### a) Shannon Index

#pdf("ShannonIndex_Allgroups.pdf") 
shannon_allgroups <- plot_richness(physeq.alpha1, x="Group",  measures="Shannon") + 
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                      "#AD6F3B")) + 
  geom_signif(comparisons = list(c("CD", "CD_RL"), c("UC","UC_RL"), c("CD", "HR_CD"), c("CD_RL", "Healthy"), 
                                 c("CD", "Healthy"), c("CD_RL", "UC_RL"), c("UC", "CD")),
              map_signif_level=TRUE, 
              y_position=c(3.9, 4.1, 4.2, 4.6, 4.9, 5.2, 5.5), test = "wilcox.test") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=8)) +
  labs(x="", y="") +
  scale_x_discrete(labels=c("CD" = "CD_REM\n(n=21)", "CD_RL" = "CD_REL\n(n=13)",
                            "Healthy" = "HC\n(n=54)", "HR_CD" = "HC_CD\n(n=29)", "HR_UC"="HC_UC\n(n=29)", "UC"="UC_REM\n(n=13)", "UC_RL"="UC_REL\n(n=18)")) +
  ylim(c(0, 5.6)) + 
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9))
print(shannon_allgroups)
#dev.off()

### b) Chao1 Index

#pdf("Chao1Index_Allgroups.pdf") 
chao1_allgroups <- plot_richness(physeq.alpha1, x="Group",  measures="Chao1") +  
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                      "#AD6F3B")) + 
  geom_signif(comparisons = list(c("CD", "CD_RL"), c("UC","UC_RL"), c("CD", "HR_CD"), c("CD_RL", "Healthy"), 
                                 c("CD", "Healthy"),
                                 c("CD_RL", "UC_RL"), c("UC", "CD")),
              map_signif_level=TRUE, 
              y_position=c(85, 85, 94, 105, 114, 122, 129), test = "wilcox.test") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=8)) +
  labs(x="", y="") +
  scale_x_discrete(labels=c("CD" = "CD_REM\n(n=21)", "CD_RL" = "CD_REL\n(n=13)",
                            "Healthy" = "HC\n(n=54)", "HR_CD" = "HC_CD\n(n=29)", "HR_UC"="HC_UC\n(n=29)", "UC"="UC_REM\n(n=13)", "UC_RL"="UC_REL\n(n=18)")) +
  ylim(c(0,130)) + 
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9))
print(chao1_allgroups)
#dev.off()

allgroups <-ggarrange(shannon_allgroups, chao1_allgroups, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
#ggsave(file="shannon_chao1_allgroups.png", allgroups, height=4, width=8)

####################################################################
# Figure 3.11: Shannon index over time in CD Relapse and remission #
####################################################################

### CD Remission ###
newmeta2 <- subset(metadata, metadata$Group=="CD")
meta.physeq1 = sample_data(newmeta2) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### Shannon 
#pdf("ShannonIndex_CD_REM_timepoints.pdf") 
shannon_cd_rem <- plot_richness(physeq.alpha1, x="Timepoint", measures="Shannon") + 
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=9)) +
  labs(x="", y="") +
  ylim(c(0,4)) + 
  #geom_signif(comparisons = list(c("TP1", "TP4"), c("TP0", "TP1"), c("TP0", "TP2"),
  #                              c("TP0", "TP3"), c("TP1","TP2"), c("TP1", "TP3"), c("TP1", "TP4")),
  #          map_signif_level=TRUE, 
  #         y_position=c(3.4, 3.7, 4.0, 4.2, 4.2, 4.5, 4.8), test = "wilcox.test") + # wilcox.test es el test por defecto
  #y_position=c(70, 73, 75, 77, 79, 80, 80.5), test = "wilcox.test") + # wilcox.test es el test por defecto
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9)) 
# No se han encontrado diferencias significativas entre los diferentes grupos
print(shannon_cd_rem)
#dev.off()


### CD Relapse ###
newmeta2 <- subset(metadata, metadata$Group=="CD_RL")
meta.physeq1 = sample_data(newmeta2) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### Shannon
#pdf("ShannonIndex_CD_REL_timepoints.pdf") 
shannon_cd_rel <- plot_richness(physeq.alpha1, x="Timepoint", measures="Shannon") + 
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=9)) +
  labs(x="", y="") +
  ylim(c(0,4)) + 
  # ylim(c(0,80)) + 
  # geom_signif(comparisons = list(c("TP1", "TP4"), c("TP0", "TP1"), c("TP0", "TP2"),
  #                               c("TP0", "TP3"), c("TP1","TP2"), c("TP1", "TP3"), c("TP1", "TP4")),
  #           map_signif_level=TRUE, 
  #          y_position=c(3.4, 3.7, 4.0, 4.2, 4.2, 4.5, 4.8), test = "wilcox.test") + # wilcox.test es el test por defecto
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9)) 
print(shannon_cd_rel)
# No se han encontrado diferencias significativas entre los diferentes grupos
#dev.off()

shannon_tp_cd <-ggarrange(shannon_cd_rem, shannon_cd_rel, 
                          labels = c("A", "B"),
                          ncol = 2, nrow = 1)
#ggsave(file="shannon_tp_cd.png", shannon_tp_cd, height=4, width=8, units="in")

####################################################################
# Figure 3.12: Shannon index over time in UC Relapse and remission #
####################################################################

### UC Remission ###
newmeta2 <- subset(metadata, metadata$Group=="UC")
meta.physeq1 = sample_data(newmeta2) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### Shannon 
#pdf("ShannonIndex_UC_REM_timepoints.pdf") 
shannon_uc_rem <- plot_richness(physeq.alpha1, x="Timepoint", measures="Shannon") + 
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=9)) +
  labs(x="", y="") +
  ylim(c(0,4)) + 
  #geom_signif(comparisons = list(c("TP1", "TP4"), c("TP0", "TP1"), c("TP0", "TP2"),
  #                              c("TP0", "TP3"), c("TP1","TP2"), c("TP1", "TP3"), c("TP1", "TP4")),
  #          map_signif_level=TRUE, 
  #         y_position=c(3.4, 3.7, 4.0, 4.2, 4.2, 4.5, 4.8), test = "wilcox.test") + # wilcox.test es el test por defecto
  #y_position=c(70, 73, 75, 77, 79, 80, 80.5), test = "wilcox.test") + # wilcox.test es el test por defecto
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9)) 
# No se han encontrado diferencias significativas entre los diferentes grupos
print(shannon_uc_rem)
#dev.off()

### UC Relapse ###
newmeta2 <- subset(metadata, metadata$Group=="UC_RL")
meta.physeq1 = sample_data(newmeta2) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### Shannon
#pdf("ShannonIndex_UC_REL_timepoints.pdf") 
shannon_uc_rel <- plot_richness(physeq.alpha1, x="Timepoint", measures="Shannon") + 
  geom_boxplot(fill=c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578")) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=9)) +
  labs(x="", y="") +
  ylim(c(0,4)) + 
  # ylim(c(0,80)) + 
  # geom_signif(comparisons = list(c("TP1", "TP4"), c("TP0", "TP1"), c("TP0", "TP2"),
  #                               c("TP0", "TP3"), c("TP1","TP2"), c("TP1", "TP3"), c("TP1", "TP4")),
  #           map_signif_level=TRUE, 
  #          y_position=c(3.4, 3.7, 4.0, 4.2, 4.2, 4.5, 4.8), test = "wilcox.test") + # wilcox.test es el test por defecto
  theme(axis.text =  element_text( face = "bold")) + 
  theme(axis.text.y = element_text(size=9)) 
print(shannon_uc_rel)
# No se han encontrado diferencias significativas entre los diferentes grupos
#dev.off()
shannon_tp_uc <-ggarrange(shannon_uc_rem, shannon_uc_rel, 
                          labels = c("A", "B"),
                          ncol = 2, nrow = 1)
#ggsave(file="shannon_tp_uc.png", shannon_tp_uc, height=4, width=8, units="in")

###########################################
# Figure 3.13: PCoA's in different groups #
###########################################

newmeta1 <- subset(metadata, metadata$Timepoint=="TP0")
meta.physeq1 = sample_data(newmeta1) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

### a) For healthy individuals

physeq.sub <- physeq.alpha1 %>%
  subset_samples(Disease_group == "HR_CD" | Disease_group == "HR_UC" | Disease_group == "Healthy")
dist.ord <- ordinate(physeq.sub, "PCoA", "bray") 

#pdf("PCOA_HRCD_HRUC_HC.pdf")
pcoa1 <- plot_ordination(physeq.sub, dist.ord, type="samples", color="Disease_group") + 
  geom_point(size=3) + stat_ellipse( type = "norm") +
  scale_color_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                "#AD6F3B")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  theme(axis.title = element_text(size=20,face = "bold")) +
  theme(legend.text = element_text(size=22, face="bold")) 
print(pcoa1)
#dev.off()

### b) For CD, UC and Healthy individuals

physeq.sub <- physeq.alpha1 %>%
  subset_samples(Disease_group == "CD" | Disease_group == "UC" | Disease_group == "Healthy")
dist.ord <- ordinate(physeq.sub, "PCoA", "bray") 

#pdf("PCOA_CD_UC_HC.pdf") 
pcoa2 <- plot_ordination(physeq.sub, dist.ord, type="samples", color="Disease_group") + 
  geom_point(size=3) + stat_ellipse( type = "norm") +
  scale_color_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                "#AD6F3B")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  theme(axis.title = element_text(size=20,face = "bold")) +
  theme(legend.text = element_text(size=22, face="bold")) 
print(pcoa2)
#dev.off()

### c) For CD and UC in relapsing and remitting conditions

physeq.sub <- physeq.alpha1 %>%
  subset_samples(Group == "CD" | Group == "UC" | Group == "CD_RL"| Group=="UC_RL")
dist.ord <- ordinate(physeq.sub, "PCoA", "bray") 

#pdf("PCOA_CDRelRem_UCRelRem.pdf") 
pcoa3 <- plot_ordination(physeq.sub, dist.ord, type="samples", color="Group") + 
  geom_point(size=3) + stat_ellipse(type = "norm") +
  scale_color_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
                                "#AD6F3B")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.y = element_text(size=20, face="bold"),
        axis.text.x = element_text(size=20, face="bold")) + 
  theme(axis.title = element_text(size=20,face = "bold")) +
  theme(legend.text = element_text(size=22, face="bold")) 
print(pcoa3)
#dev.off()
pcoas <-ggarrange(pcoa1, pcoa2, pcoa3, 
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
ggsave(file="PCOAS.png", pcoas, height=14, width=18, units="in")


############################
## COMPOSITIONAL ANALYSIS ##
############################

# Re-read the data to take only the significant genus obtained from the ANCOM results:

otu <- as.data.frame(read_excel("table-L6.xlsx")) # otu table in terms of absolute frequency L6 - Genus level
head(otu)
rownames(otu) <- otu[,1]
otu <- otu[,-1]

ancom_L6 <- as.data.frame(read_excel("ancom_L6_results.xlsx")) # otu table in terms of absolute frequency L6 - Genus level
head(ancom_L6)
ancom_L6$Reject_null_hypothesis <- as.factor(ancom_L6$Reject_null_hypothesis)
signif <- subset(ancom_L6, ancom_L6$Reject_null_hypothesis=="True")

# Select only the significant otus (results from ANCOM)
otu_signif <- subset(otu, rownames(otu) %in% signif$Otu)

otu <- otu_signif
tax <- rbind(unlist(strsplit(rownames(otu)[1], "[;]"))) # to create taxonomy table
for(i in 2:nrow(otu)){
  tax <- rbind(tax, unlist(strsplit(rownames(otu)[i], "[;]")))
}
taxonomy_L6 <- as.data.frame(tax)
colnames(taxonomy_L6) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(otu) <- 1:nrow(otu)
dim(otu)
rownames(otu) <- paste0("otu_0", rownames(otu))
rownames(taxonomy_L6) <- paste0("otu_0", rownames(taxonomy_L6))
matches <- rownames(metadata)
otus <- otu[, as.character(matches)]
rownames(otus) <- rownames(taxonomy_L6)

OTU.physeq = otu_table(as.matrix(otus), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(taxonomy_L6))
meta.physeq = sample_data(metadata) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq.alpha

newmeta1 <- subset(metadata, metadata$Timepoint=="TP0")
meta.physeq1 = sample_data(newmeta1) # rownames de metadata han de ser l'ID dels pacients
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)


# The subsequent data used (physeq.alpha1) is only for individuals at TP0.
# Subset only bacterias:

newmeta1 <- subset(metadata, metadata$Timepoint=="TP0")
meta.physeq1 = sample_data(newmeta1) # rownames de metadata han de ser el ID de los pacientes
physeq.alpha1 = phyloseq(OTU.physeq, tax.physeq, meta.physeq1)

justbacteria <- physeq.alpha1 %>%
  subset_taxa(
    Kingdom == "k__Bacteria"
  )
justbacteria

#################################################
#### Figures 3.14: Genus level in all groups ####
#################################################

genusabundance <- justbacteria %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

pr <- genusabundance[, c("Genus", "Disease", "Abundance", "Disease_group", "Group", "Group2")]

gen <- c()
disease <- c()
means <- c()

for(i in 1:length(levels(pr$Genus))){
  
  phy <- levels(pr$Genus)[i]
  
  for(j in 1:length(levels(pr$Disease_group))){
    
    dis <- levels(pr$Disease_group)[j]
    newdf <- subset(pr, pr$Genus==phy & pr$Disease_group==dis)
    means <- c(means, mean(newdf$Abundance))
    gen <- c(gen, phy)
    disease <- c(disease, dis)
  }
}
all <- as.data.frame(cbind(gen, disease ,means))
colnames(all) <- c("Genus", "Disease", "Abundance")
all$Abundance <- as.numeric(paste(all$Abundance))
#pdf("Compositional_Genus_Allgroups.pdf") 
ggplot(all) +
  geom_col(mapping = aes(x = Disease, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  xlab(NULL)+
  theme_minimal()+
  ylim(0,1) + 
  labs(x="", y="Relative Abundance (%)") + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete(labels=c("CD" = "CD", "Healthy" = "HC",
                            "HR_CD" = "HC-CD", "HR_UC"="HC-UC","UC"="UC")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=13, face = "bold"),
        axis.text.y = element_text(size=13, face = "bold")) +
  theme(title = element_text(size=14, face = "bold")) +
  ggtitle("Genus")
#dev.off()

#########################################################
#### Figure 3.15: Healthy controls and IBD subgroups #### 
#########################################################

head(otu_signif)
newdf <- as.data.frame(t(otu_signif))
newdf$ID <- rownames(newdf)
library(dplyr)
library(reshape)
df = merge(x=newdf,y=metadata,by="ID")
dim(df)

##################################################
######## A) BOXPLOT TP0 FOR IBD vs HEALTHY ####### 
##################################################

tp0_ibd_healthy <- subset(df, (df$Timepoint=="TP0" & (df$Group=="UC" | df$Group=="CD" | df$Group=="CD_RL" | df$Group=="UC_RL" |
                                                        df$Group=="Healthy")))
tp0_ibd_healthy$IBD_HEALTHY <- ifelse(tp0_ibd_healthy$Group=="Healthy", "Healthy", "IBD")
tp0_ibd_healthy <- tp0_ibd_healthy[, c("IBD_HEALTHY",
                                       "k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus",
                                       "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium",
                                       "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium",
                                       "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                                       "k__Bacteria;p__Actinobacteria;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella",
                                       "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium",
                                       "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Oscillospira",
                                       "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus")]
#pdf("Boxplot_TP0_IBD_HEALTHY.pdf") 
p3 <- ggplot(melt(tp0_ibd_healthy), aes(x = variable, y = value, fill=IBD_HEALTHY)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus" = "Coprobacillus",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium" = "Faecalibacterium",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium" = "Clostridium",
                            "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia" = "Escherichia",
                            "k__Bacteria;p__Actinobacteria;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella" = "Collinsella",
                            "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium" = "Fusobacterium",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Oscillospira" = "Oscillospira",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Ruminococcus" = "Ruminococcus")) +
  scale_fill_manual(values = c("#CBD588", "#DA5724")) +
  theme(axis.text.x = element_text(angle = 20, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  labs(fill = "Group") + 
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13, face="bold"))
print(p3)
#dev.off()

##################################################
######## B) BOXPLOT TP0 FOR CD vs HEALTHY ######## 
##################################################

tp0_cd_healthy <- subset(df, (df$Timepoint=="TP0" & (df$Group=="CD" | df$Group=="CD_RL" | df$Group=="Healthy")))
tp0_cd_healthy$CD_HEALTHY <- ifelse(tp0_cd_healthy$Group=="Healthy", "Healthy", "CD")
tp0_cd_healthy <- tp0_cd_healthy[, c("CD_HEALTHY",
                                     "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium",
                                     "k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus",
                                     "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                                     "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium")]
#pdf("Boxplot_TP0_CD_HEALTHY.pdf") 
p1 <- ggplot(melt(tp0_cd_healthy), aes(x = variable, y = value, fill=CD_HEALTHY)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium" = "Clostridium",
                            "k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus" = "Coprobacillus",
                            "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia" = "Escherichia",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium" = "Faecalibacterium")) +
  scale_fill_manual(values = c("#5F7FC7", "#CBD588")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  labs(fill = "Group") + 
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13, face="bold"))
print(p1)
#dev.off()

##################################################
######## C) BOXPLOT TP0 FOR UC vs HEALTHY ######## 
##################################################

tp0_uc_healthy <- subset(df, (df$Timepoint=="TP0" & (df$Group=="UC" | df$Group=="UC_RL" | df$Group=="Healthy")))
tp0_uc_healthy$UC_HEALTHY <- ifelse(tp0_uc_healthy$Group=="Healthy", "Healthy", "UC")
tp0_uc_healthy <- tp0_uc_healthy[, c("UC_HEALTHY",
                                     "k__Bacteria;p__Actinobacteria;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella")]
#pdf("Boxplot_TP0_UC_HEALTHY.pdf") 
p2 <- ggplot(melt(tp0_uc_healthy), aes(x = variable, y = value, fill=UC_HEALTHY)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("k__Bacteria;p__Actinobacteria;c__Coriobacteriia;o__Coriobacteriales;f__Coriobacteriaceae;g__Collinsella" = "Collinsella")) +
  scale_fill_manual(values = c("#CBD588", "orange")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  labs(fill = "Group") + 
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13, face="bold"))
print(p2)
#dev.off()

############################################
######## D) BOXPLOT TP0 FOR UC vs CD ####### 
############################################

tp0_cd_uc <- subset(df, (df$Timepoint=="TP0" & (df$Group=="UC" | df$Group=="CD" | df$Group=="CD_RL" | df$Group=="UC_RL")))
tp0_cd_uc$Group <- factor(tp0_cd_uc$Group)
levels(tp0_cd_uc$Group) <- c("CD", "CD", "UC", "UC")

tp0_cd_uc <- tp0_cd_uc[, c("Group",
                           "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia",
                           "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium",
                           "k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus",
                           "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium")]
#pdf("Boxplot_TP0_CD_UC.pdf") 
p4 <- ggplot(melt(tp0_cd_uc), aes(x = variable, y = value, fill=Group)) +
  geom_boxplot() + 
  theme_bw() + 
  scale_x_discrete(labels=c("k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia" = "Escherichia",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Clostridium" = "Clostridium",
                            "k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Coprobacillus" = "Coprobacillus",
                            "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__Faecalibacterium" = "Faecalibacterium")) +
  scale_fill_manual(values = c("#5F7FC7", "orange")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5, size=12)) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  labs(x="", y="") + 
  labs(fill = "Group") + 
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=13, face="bold"))
print(p4)
#dev.off()

boxplots <-ggarrange(p3, p1, p2, p4 ,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2)
# ggsave(file="boxplots.png", boxplots, height=12, width=14, units="in")

