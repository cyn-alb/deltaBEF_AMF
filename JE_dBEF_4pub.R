# SSU amplicons were sequenced with Illumina MiSeq and the resulting sequences processed with the dadasnake
# pipeline (Weißbecker et al 2021). Sequences were trimmed and filtered and reads assigned to ASV, which
# were matched against the SILVA database to filter for fungi. 
# However for AMF not ASV but virtual taxa (VT) are more commonly used, therefore fungal ASV were transformed
# to VT by blasting against the Maarjam database. 

#### 0 Prep data and session ####
# library("hiReadsProcessor")
# library(seqRFLP)
# library(mixOmics) 
library(vegan)
library(tidyr)
library(dplyr)
library(phyloseq)
library(ggplot2)
library("readxl")
library(RColorBrewer)
library(paletteer)
library(nlme)
library(stringr)
# library("mia")

set.seed(2)


##### 0.1 prep sequences for blast against Maarjam #####
setwd("C:/Users/albracht/Documents/R-Projekte/JE_dBEF")
library("hiReadsProcessor")

# read in fasta file sequences (fungal sequences)
# this fasta needs to be split into subsets of max 300 sequences, 
# as Maarjam can only process max 300 sequences at once
split.seqTab <- readDNAStringSet("filtered.seqs.fasta")
splitSeqsToFiles(split.seqTab, totalFiles = 30, suffix ="fasta", filename = "JE_dBEF_splitFasta.fa", 
                 outDir = getwd())
# blast files against https://maarjam.botany.ut.ee/
# the .result output files were parsed with a python script to a .tsv with only the best VT hit

##### 0.2 compile VT abundance table #####
# read in .tsv output files and merge them
VT1 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_1.tsv", header = TRUE, sep = "\t")
VT2 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_2.tsv", header = TRUE, sep = "\t")
VT3 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_3.tsv", header = TRUE, sep = "\t")
VT4 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_4.tsv", header = TRUE, sep = "\t")
VT5 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_5.tsv", header = TRUE, sep = "\t")
VT6 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_6.tsv", header = TRUE, sep = "\t")
VT7 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_7.tsv", header = TRUE, sep = "\t")
VT8 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_8.tsv", header = TRUE, sep = "\t")
VT9 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_9.tsv", header = TRUE, sep = "\t")
VT10 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_10.tsv", header = TRUE, sep = "\t")
VT11 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_11.tsv", header = TRUE, sep = "\t")
VT12 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_12.tsv", header = TRUE, sep = "\t")
VT13 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_13.tsv", header = TRUE, sep = "\t")
VT14 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_14.tsv", header = TRUE, sep = "\t")
VT15 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_15.tsv", header = TRUE, sep = "\t")
VT16 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_16.tsv", header = TRUE, sep = "\t")
VT17 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_17.tsv", header = TRUE, sep = "\t")
VT18 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_18.tsv", header = TRUE, sep = "\t")
VT19 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_19.tsv", header = TRUE, sep = "\t")
VT20 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_20.tsv", header = TRUE, sep = "\t")
VT21 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_21.tsv", header = TRUE, sep = "\t")
VT22 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_22.tsv", header = TRUE, sep = "\t")
VT23 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_23.tsv", header = TRUE, sep = "\t")
VT24 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_24.tsv", header = TRUE, sep = "\t")
VT25 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_25.tsv", header = TRUE, sep = "\t")
VT26 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_26.tsv", header = TRUE, sep = "\t")
VT27 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_27.tsv", header = TRUE, sep = "\t")
VT28 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_28.tsv", header = TRUE, sep = "\t")
VT29 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_29.tsv", header = TRUE, sep = "\t")
VT30 <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_VT_parse_30.tsv", header = TRUE, sep = "\t")

VT_dBEF <- bind_rows(VT1, VT2, VT3, VT4, VT5, VT6, VT7, VT8, VT9, VT10,
                    VT11, VT12, VT13, VT14, VT15, VT16, VT17, VT18, VT19, VT20,
                    VT21, VT22, VT23, VT24, VT25, VT26, VT27, VT28, VT29, VT30)

# how many ASV did we have? How many VTs were they assigned to? 
length(unique(VT_dBEF$Contig)) #7798 ASVs 
length(unique(VT_dBEF$VTX)) #108 VT

# the VT HitNames come as long string like this one: 
# |gnl||BL_ORD_ID||45479 ||gi||47275||gb||HF559284|| Glomeraceae Glomus Torrecillas14 Fu2 |
# so we shorten the HitName string for only Species and VT name
VT_dBEF$Hitshort <- VT_dBEF$HitName %>%  stringr::str_remove(pattern = "[|]gnl[|]+BL_ORD_ID[|]+[[:digit:]]+ [|]+gi[|]+[[:digit:]]+[|]+gb[|]+[[:alpha:]][[:alpha:]][[:digit:]]+[|]+")
# some contigs (ASV) had several best hits, aggregate them to remove duplets
VT_dBEF_agg <- aggregate(VT_dBEF$VTX, list(VT_dBEF$Contig),function(x)paste(unique(x[x!=""]),sep=";",collapse=";"))

##### 0.3 unassigned VTs #####
# Unfortunately Maarjam can't assign all contigs to VTs
# (check again if these ASV are actually fungi by e.g. blasting against NCBI)
length(which(VT_dBEF_agg$x=="")) #1581 ASVs not assigned to VT
VT_dBEF_VTX <- VT_dBEF_agg %>% filter(x!="")  %>% rename(OTU = Group.1 , VT=x)
VT_dBEF_nVT <- VT_dBEF_agg %>% filter(x=="") %>%  rename(OTU=Group.1)

# how abundant are the nonVT taxa? can we filter singletons? -> load in ASV abundance table to check
filtered.seqTab <- readRDS("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/filtered.seqTab.RDS")
ASV_dBEF <- filtered.seqTab[,c(1,3:242)] 
colnames(ASV_dBEF) = gsub(patter = "-AMF", replacement = "", x=colnames(ASV_dBEF)) # clean colnames
rownames(ASV_dBEF) <- ASV_dBEF$OTU
ASV_dBEF <- ASV_dBEF %>% select(-OTU) # df of 9097 x 240
VT_dBEF_nVT_abund <- merge(ASV_dBEF, VT_dBEF_nVT, by.x = 0, by.y = "OTU", all.x = FALSE) 
# abundance table of the 1581 non-VT-assigned ASV

rowSums(VT_dBEF_nVT_abund[,2:241]) #lots of singletons
VT_dBEF_nVT_abund = VT_dBEF_nVT_abund %>% filter(rowSums(VT_dBEF_nVT_abund[,2:241])>2) 
# 446 ASV left after filtering
VT_dBEF_nVT_abund$OTU = VT_dBEF_nVT_abund$Row.names


##### 0.4 adding the resident species sequence data #####
# ASV were blast against Maarjam and results parsed with python (EVE)
VT1r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_1.tsv", header = TRUE, sep = "\t")
VT2r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_2.tsv", header = TRUE, sep = "\t")
VT3r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_3.tsv", header = TRUE, sep = "\t")
VT4r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_4.tsv", header = TRUE, sep = "\t")
VT5r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_5.tsv", header = TRUE, sep = "\t")
VT6r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_6.tsv", header = TRUE, sep = "\t")
VT7r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_7.tsv", header = TRUE, sep = "\t")
VT8r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_8.tsv", header = TRUE, sep = "\t")
VT9r <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_residents/VT/new_VTs_9.tsv", header = TRUE, sep = "\t")

VT_res <- bind_rows(VT1r, VT2r, VT3r, VT4r, VT5r, VT6r, VT7r, VT8r, VT9r)

length(unique(VT_res$Contig)) # 1371 ASVs
length(unique(VT_res$VTX)) # 80 VT
#filter HitName for only Species and VT name
VT_res$Hitshort <- VT_res$HitName %>%  stringr::str_remove(pattern = "[|]gnl[|]+BL_ORD_ID[|]+[[:digit:]]+ [|]+gi[|]+[[:digit:]]+[|]+gb[|]+[[:alpha:]][[:alpha:]][[:digit:]]+[|]+")

VT_res_agg <- aggregate(VT_res$VTX, list(VT_res$Contig),function(x)paste(unique(x[x!=""]),sep=";",collapse=";"))
length(which(VT_res_agg$x=="")) # 219 ASVs keinem VT zugeordnet
VT_res_VTX <- VT_res_agg %>% filter(x!="")  %>% rename(OTU = Group.1 , VT=x)  # 1115 VTX
VT_res_nVT <- VT_res_agg %>% filter(x=="") %>%  rename(OTU=Group.1)

# residents abundance table
filtered.seqTab.res <- readRDS("C:/Users/albracht/Documents/R-Projekte/JE_residents/filtered.seqTab.RDS")
ASV_res <- filtered.seqTab.res[,c(1,3:376)] #only abundance part
rownames(ASV_res) <- ASV_res$OTU
ASV_res <- ASV_res %>% select(-OTU) # 1371 x 374

# how abundant are the nonVT taxa? can we filter singletons?
VT_res_nVT_abund <- merge(ASV_res, VT_res_nVT, by.x = 0, by.y = "OTU", all.x = FALSE) #200 ASV
rowSums(VT_res_nVT_abund[,2:375]) #lots of singletons
VT_res_nVT_abund = VT_res_nVT_abund %>% filter(rowSums(VT_res_nVT_abund[,2:375])>2) 
dim(VT_res_nVT_abund) # 192 ASV left to be treed
VT_res_nVT_abund$OTU = VT_res_nVT_abund$Row.names


##### 0.5 treeing of ASV that could not be assigned to VT #####
# create fasta file to cluster remaining ASV to VT
VT_dBEF_nVTseq <- filtered.seqTab  %>% select(Row.names, OTU) %>% rename(Sequence= Row.names)
VT_dBEF_nVTseq <- merge(VT_dBEF_nVT_abund, VT_dBEF_nVTseq, by="OTU") %>% select(OTU, Sequence)
dim(VT_dBEF_nVTseq) # 446 sequence variants from dBEF community sampling

VT_res_nVTseq <- filtered.seqTab.res %>% select(Row.names, OTU) %>% rename(Sequence=Row.names)
VT_res_nVTseq <- merge(VT_res_nVT_abund, VT_res_nVTseq, by="OTU") %>% select(OTU, Sequence)
dim(VT_res_nVTseq) # 192 sequence variants from resident species sampling

# to combine dBEF w/ residents dataset, rename "OTU" as they are from different sequencing runs
VT_dBEF_nVTseq$OTU <- str_replace(VT_dBEF_nVTseq$OTU, "OTU", "OTU_d")
VT_res_nVTseq$OTU <- str_replace(VT_res_nVTseq$OTU, "OTU", "OTU_r")
VT_all_nVTseq <- bind_rows(VT_dBEF_nVTseq, VT_res_nVTseq)
dim(VT_all_nVTseq) # 638 ASV sequences

# create fasta file from nonVT ASVs
library(seqRFLP)
nonVTseq_all.fasta <- dataframe2fas(VT_all_nVTseq, file="JE_dBEF_amf_nonVTseq_all.fasta")
# nonVT ASVs were used to construct a maximum likelihood phylogenetic tree based on a general 
# time-reversible, discrete gamma (GTR+G) model using MAFFT and raxML.Thus, following the methodology 
# of the MaarjAM database, these ASVs were assigned custom virtual taxa (VTC) with cophenetic distances
# below 0.03

# we read the results of treeing back into R
cusVT <- read.table("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/VT/JE_dBEF_amf_nonVTseq_all.tree.cluster.003.tsv", sep="\t", header=TRUE)
cusVT <- cusVT %>% mutate(VT = paste0("VTC", formatC(cluster, width=5, flag ="0")))
#amf_cusVT <- amf_cusVT %>% select(OTU,VT)
cusVT_dBEF = cusVT %>% filter(str_detect(cusVT$OTU, "OTU_d_"))
# remove the string that was added to differentiate comm. from resident species ASVs
cusVT_dBEF$OTU <- str_replace(cusVT_dBEF$OTU, "_d_", "_")
dim(cusVT_dBEF) # 446 ASV

# combine VTX und cusVT and remove double assigned VTs
VT_dBEF_com <- bind_rows(VT_dBEF_VTX , cusVT_dBEF) %>% separate(VT, c("VT","two", "three", "four"), sep=";") %>% 
  select(OTU, VT)
dim(VT_dBEF_com) # 6663
 

# add VT IDs to ASV table and aggregate reads at VT level
VT_dBEF_com <- tibble::column_to_rownames(VT_dBEF_com, var = "OTU")
setdiff(VT_dBEF_com$OTU, ASV_dBEF$OTU) # none
ASV_VT_dBEF <- merge(ASV_dBEF, VT_dBEF_com, by=0) 
dim(ASV_VT_dBEF) # 6663 x 242
length(unique(ASV_VT_dBEF$VT)) # 128 VT

ASV_VT_dBEF$VT <- as.factor(ASV_VT_dBEF$VT)
ASV_VT_dBEF_aggr <- ASV_VT_dBEF %>% select(-Row.names) %>% group_by(VT) %>% summarise(across(.cols = everything(), .fns =  sum))
dim(ASV_VT_dBEF_aggr) # 128 x 241


##### 0.6 adding taxonomic information #####
#prepare tax data
VT_dBEF_VTX <- VT_dBEF_VTX %>% tidyr::separate(VT, c("VT","B"), sep = ';', remove = TRUE) %>% select(OTU, VT) # some ASV had double hits
VT_dBEF_VTX_tax <- merge(VT_dBEF_VTX, VT_dBEF, by.x ="VT", by.y = "VTX") %>% select(VT, Hitshort) %>%  distinct(VT, .keep_all = T)
VT_dBEF_VTX_tax <- VT_dBEF_VTX_tax %>% separate(Hitshort, c("x","Family", "Genus", "Species"), sep=" ") %>% select(-"Species", -"x")
VT_dBEF_VTC <- cusVT_dBEF %>% select(VT) %>% distinct(VT)
VT_dBEF_VTC <- VT_dBEF_VTC %>% mutate(Family = paste0("Cluster", substr(VT, 6,8)), Genus = paste0("Cluster", substr(VT, 6,8))) 
VT_dBEF_tax <- bind_rows(VT_dBEF_VTX_tax, VT_dBEF_VTC) # 128 VTX / VTC
VT_dBEF_tax$VTX = VT_dBEF_tax$VT

write.table(VT_dBEF_tax, "JEdBEF_VTtax.txt", sep = "\t", quote=F, row.names = F)

##### 0.7 adding meta data #####
dBEF_samples <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/Plot_Information_DeltaBEF_2020.xlsx", sheet="Plot_Information_DeltaBEF")
dBEF_samples <- tibble::column_to_rownames(dBEF_samples, var="plot_ID")
#clean sample table 
dBEF_samples <- dBEF_samples %>% 
  mutate(NUMHERB = NUMTHERB+NUMSHERB) %>%  #Tall / short herbs functionally just herbs
  mutate(HERB = ifelse(SHERB == 1 | THERB == 1, 1, 
                       0)) %>% 
  select(-COMMENT, -NUMTHERB, -NUMSHERB, -SHERB, -THERB)
str(dBEF_samples)
dBEF_samples$treatment=as.factor(dBEF_samples$treatment)
dBEF_samples$history_plant=as.factor(dBEF_samples$history_plant)
dBEF_samples$history_soil=as.factor(dBEF_samples$history_soil)
dBEF_samples$LOG_SD = log(dBEF_samples$SOWNDIV)
dBEF_samples$SOWNDIV=as.factor(dBEF_samples$SOWNDIV)
dBEF_samples$FUNCGR=as.factor(dBEF_samples$FUNCGR)
dBEF_samples$GRASS=as.factor(dBEF_samples$GRASS)
dBEF_samples$LEG=as.factor(dBEF_samples$LEG)
dBEF_samples$HERB=as.factor(dBEF_samples$HERB)
dBEF_samples$BLOCK=as.factor(dBEF_samples$BLOCK)

# we don't need the 60 species control plots
dBEF_samples = dBEF_samples %>% filter(SOWNDIV != 60)

##### 0.8 Phyloseq object #####
VT_dBEF_phy <- tibble::column_to_rownames(ASV_VT_dBEF_aggr, var="VT")
VT_dBEF_phy = otu_table(as.matrix(VT_dBEF_phy), taxa_are_rows=TRUE)

Tax_dBEF_phy <- tibble::column_to_rownames(VT_dBEF_tax, var="VT")
Tax_dBEF_phy = tax_table(as.matrix(Tax_dBEF_phy))

Sam_dBEF_phy = sample_data(dBEF_samples)

JE_dBEF <- phyloseq(VT_dBEF_phy, Tax_dBEF_phy, Sam_dBEF_phy)
JE_dBEF <- subset_samples(JE_dBEF, SOWNDIV != "60") # 128 taxa in 228 samples

##### 0.9 What do VTC contain? #####
# how many ASV per VT ? What might be their taxonomy according to Silva?
VTC_ov = cusVT_dBEF %>% inner_join(filtered.seqTab, by="OTU") %>% 
  select(OTU, cluster, VT, taxonomy.mothur.SILVA_138_SSURef_NR99_cut)
VTC_ov = VTC_ov %>% separate(taxonomy.mothur.SILVA_138_SSURef_NR99_cut, 
                             c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s"), 
                             sep = ";")
table(VTC_ov$VT) # nr of ASV per VTC
length(unique(VTC_ov$l[VTC_ov$cluster == 1])) #11 Glomus, uncultured_unclassified, Rhizophagus
# Funneliformis, Septoglomus, Glomeromycetes_unclassified, Diversisporaceae_unclassified, 
# Glomeraceae_unclassified, Glomerales_unclassified, Claroideoglomus, Otospora  
length(unique(VTC_ov$l[VTC_ov$cluster == 2])) #2 Glomeraceae_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 3])) #2 Sclerocystis, Glomeraceae_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 4])) #1 Archaeospora
length(unique(VTC_ov$l[VTC_ov$cluster == 5])) #1 Paraglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 6])) #1 Paraglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 12]))#1 Rhizophagus
length(unique(VTC_ov$l[VTC_ov$cluster == 13]))#1 Septoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 14]))#1 Septoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 15]))#2 Archaeospora; uncultured_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 16]))#1 Glomeromycetes_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 17]))#4 Glomeromycetes_unclassified, Glomeromycetes_unclassified,
# Funneliformis, Septoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 18]))#3 Diversisporales_unclassified, Entrophospora_unclassified, Otospora
length(unique(VTC_ov$l[VTC_ov$cluster == 19]))#1 Paraglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 20]))#7 unclassified Diversisporaceae, Diversispora, unclassified Glomeromycetes, 
# Funneliformis, unclassified Diversisporales, Archaeospora, Claroideoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 21]))#1 Scutellospora
length(unique(VTC_ov$l[VTC_ov$cluster == 22]))#1 Diversispora_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 23]))#1 Diversispora_unclassified
length(unique(VTC_ov$l[VTC_ov$cluster == 24]))#1 Claroideoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 25]))#1 Claroideoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 26]))#1 Claroideoglomus
length(unique(VTC_ov$l[VTC_ov$cluster == 27]))#1 uncultured_unclassified

#### 1 data pre-processing ####
# different versions of abundances created to use for different analyses
# unfiltered:                       Phi-Factor
# unfiltered + rarefied:            alpha Diversity
# filtered (Top 50 VT):             relative abundance, differential abundance, beta Diversity / Permanova
# filtered (Top 50 VT) & rarefied:  turnover analysis
# f+rf + CLR-transformed:           Parafac, PCA

# How many reads at start? 
sum(sample_sums(otu_table(JE_dBEF))) # 3967994
mean(sample_sums(otu_table(JE_dBEF))) # mean 17403.48 +- 5980.216)
sd(sample_sums(otu_table(JE_dBEF)))
min(sample_sums(otu_table(JE_dBEF))) # range 2450 to 32249
max(sample_sums(otu_table(JE_dBEF)))

##### 1.1 filtering #####
#filtering to VT that have at least 1 count in min 2 % of samples (5 of 240 samples)
JE_dBEF_f = filter_taxa(JE_dBEF, function(x) sum(x > 0) > (0.2 * length(x)), TRUE) 
# 50 VT / 228 samples --> 75 VT removed
VT_dBEF_f <- as.matrix(otu_table(JE_dBEF_f))
write.table(t(VT_dBEF_f), file = "JEdBEF_fVT_norf_new.txt", sep = "\t", quote = F, col.names = NA) 


#JE_dBEF_f = filter_taxa(JE_dBEF, function(x) sum(x > 0) > (0.3 * length(x)), TRUE) 
#35 taxa / 228 samples --> 92 VT removed
#JE_dBEF_f = filter_taxa(JE_dBEF, function(x) sum(x > 03 > (0.2 * length(x)), TRUE) 
#27 taxa / 228 samples --> 100 VT removed

hist(rowSums(otu_table(JE_dBEF) >0), breaks = 50)

#how many reads were removed? 
min(sample_sums(JE_dBEF)- sample_sums(JE_dBEF_f)) 
max(sample_sums(JE_dBEF)- sample_sums(JE_dBEF_f)) # range 0 - 6013
sum(otu_table(JE_dBEF)) # 3967994
sum(otu_table(JE_dBEF_f)) # 3871665
((sum(otu_table(JE_dBEF))-sum(otu_table(JE_dBEF_f)))/sum(otu_table(JE_dBEF)))*100 #2.45 % reads removed

#how many zeros left in OTU table?
sum(otu_table(JE_dBEF_f) == 0) # 5726
sum(otu_table(JE_dBEF_f) == 0) / (dim(otu_table(JE_dBEF_f))[2] * dim(otu_table(JE_dBEF_f))[1]) * 100 # 49.24 % are zeros 
sum(otu_table(JE_dBEF) == 0) / (dim(otu_table(JE_dBEF))[2] * dim(otu_table(JE_dBEF))[1]) * 100  # before 76.55 % zeros

#how many reads per samples removed?
fps <- ((colSums(otu_table(JE_dBEF)) - colSums(otu_table(JE_dBEF_f)))/colSums(otu_table(JE_dBEF)))*100
plot(fps) + text(fps, names(fps), cex=0.6, pos=4, col="red")
#B4A22D2 - 74 % removed
#B4A22D3 - 30.5 % removed
#B3A01D3 - 32.7 % 
#B1A01D2 - 26 % 
#B4A08D3 - 20 % 

# which taxa were removed? 
VT_out = anti_join(as.data.frame(otu_table(JE_dBEF)), as.data.frame(otu_table(JE_dBEF_f)))
dim(VT_out) # 77 * 228
VT_out$sum = rowSums(VT_out)
VT_out = VT_out %>% mutate(count = rowSums(. > 0)-1)
VT_out = merge(VT_out, VT_dBEF_tax, by.x= 0, by.y = "VT")
VT_out = VT_out[order(VT_out[,232], decreasing = T),]
hist(table(VT_out['Family']))

# how much lost per treatment / Div level? 
hist(colSums(VT_out[,2:229]))
boxplot(colSums(VT_out[,2:229])~dBEF_samples$treatment)
boxplot(colSums(VT_out[,2:229])/colSums(otu_table(JE_dBEF)) ~ dBEF_samples$treatment)
boxplot(colSums(VT_out[,2:229])/colSums(otu_table(JE_dBEF)) ~ dBEF_samples$SOWNDIV)
boxplot(colSums(VT_out[,2:229])/colSums(otu_table(JE_dBEF)) ~ dBEF_samples$treatment+dBEF_samples$SOWNDIV, las=2) 
# 1 outlier for D3.2, rest rather similar
boxplot(colSums(otu_table(JE_dBEF)) ~ dBEF_samples$BLOCK+dBEF_samples$treatment, las=2)

##### 1.2 rarefaction #####
sample_sums(JE_dBEF_f)[order(sample_sums(JE_dBEF_f), decreasing = T)] # 2085 is lowest
JE_dBEF_frf <- rarefy_even_depth(JE_dBEF_f, rngseed=1, sample.size=2000 , replace=FALSE) # no OTUs removed

sample_sums(JE_dBEF)[order(sample_sums(JE_dBEF), decreasing = T)] # 2450
JE_dBEF_rf <- rarefy_even_depth(JE_dBEF, rngseed=1, sample.size=2400 , replace=FALSE) # 20 OTUs removed

rarcol<-colorRampPalette(brewer.pal(8,"Dark2"))(25) #custom color palette based on brewer Dark2 palette

pdf("Graphics/20230504_dBEF_1_rarecurves.pdf")
rarecurve(t(as.data.frame(otu_table(JE_dBEF))),step=50,cex=0.5,col=rarcol, xlab="Number of Reads",ylab="Number of OTUs", label=F) 
rarecurve(t(as.data.frame(otu_table(JE_dBEF_f))),step=50,cex=0.5,col=rarcol, xlab="Number of Reads",ylab="Number of OTUs", label=F) 
dev.off()

##### 1.3 CLR transformation #####
VT_dBEF_frf = as(otu_table(JE_dBEF_frf), "matrix")
write.table(t(VT_dBEF_frf), file = "JEdBEF_fVT_v1_new.txt", sep = "\t", quote = F, col.names = NA) # to export for PARAFAC

library(robCompositions)
# matrix[matrix==0] <- 1
# matrix_clr <- t(cenLR(t(matrix))$x.clr)
VT_dBEF_frf[VT_dBEF_frf == 0] <- 1
VT_dBEF_frf_clr <- t(cenLR(t(VT_dBEF_frf))$x.clr)

# filtered sum(x>0) > 0.2 #50 taxa

# write.table(t(amf_VT_frf), file = "JEdBEF_fVT_v2.txt", sep = "\t",quote = F, col.names = NA) 
# #filtered sum(x>0) > 0.3 #35 taxa
# write.table(t(amf_VT_frf), file = "JEdBEF_fVT_v3.txt", sep = "\t",quote = F, col.names = NA) 
# #filtered sum(x>3) > 0.2 #27 taxa

#CLR phyloseq object
VT_dBEF_phy_clr = otu_table(VT_dBEF_frf_clr, taxa_are_rows = T)
JE_dBEF_clr <- phyloseq(VT_dBEF_phy_clr, Tax_dBEF_phy, Sam_dBEF_phy) # 51 VT on 228 samples


#### 2 alpha-Diversity ####
##### 2.1 alpha-Diversity indices #####
dBEF_rich <- estimate_richness(JE_dBEF_rf, split=TRUE,
                              measures = c("Observed", "Shannon","Simpson", "Chao1"))
dBEF_rich$Evenness<- dBEF_rich$Shannon/log(dBEF_rich$Observed) #pielous evenness

#add AlphaDiv directly to (phyloseq) sample data
dBEF_samples <- merge(dBEF_samples, dBEF_rich, by.x=0, by.y = 0)
dBEF_samples = dBEF_samples %>% tibble::column_to_rownames(var = "Row.names")

dBEF_rich_phy <- sample_data(dBEF_rich)
JE_dBEF <- merge_phyloseq(JE_dBEF, dBEF_rich_phy)

#mean alpha div per treatment
mean(dBEF_samples$Observed[dBEF_samples$treatment=="D1"]) # 18.80
mean(dBEF_samples$Observed[dBEF_samples$treatment=="D2"]) # 18.16
mean(dBEF_samples$Observed[dBEF_samples$treatment=="D3"]) # 19.22

mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==16 | dBEF_samples$treatment == "D1"]) # 19.34
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==8 | dBEF_samples$treatment == "D1"])  # 19.06
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==4 | dBEF_samples$treatment == "D1"])  # 18.82
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==2 | dBEF_samples$treatment == "D1"])  # 18.56
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==1 | dBEF_samples$treatment == "D1"])  # 18.06

mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==16 | dBEF_samples$treatment == "D2"]) # 18.96
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==8 | dBEF_samples$treatment == "D2"])  # 18.89
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==4 | dBEF_samples$treatment == "D2"])  # 18.34
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==2 | dBEF_samples$treatment == "D2"])  # 18.43
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==1 | dBEF_samples$treatment == "D2"])  # 17.36

mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==16 | dBEF_samples$treatment == "D3"]) # 19.83
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==8 | dBEF_samples$treatment == "D3"])  # 19.23
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==4 | dBEF_samples$treatment == "D3"])  # 18.87
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==2 | dBEF_samples$treatment == "D3"])  # 18.81
mean(dBEF_samples$Observed[dBEF_samples$SOWNDIV==1 | dBEF_samples$treatment == "D3"])  # 18.32

#test effects on alpha div
pairwise.wilcox.test(sample_data(JE_dBEF)$Observed, sample_data(JE_dBEF)$LOG_SD)
#         1       2        4       8       16  
#   2  0.29354 -       -       -      
#   4  0.10081 0.74812 -       -      
#   8  0.00098 0.29354 0.74812 -      
#   16 0.00013 0.07472 0.28528 0.65832           --> monocultures different from high-diverse comm
pairwise.wilcox.test(sample_data(JE_dBEF)$Observed, sample_data(JE_dBEF)$treatment)
#        D1   D2  
#   D2 0.78 -   
#   D3 0.96 0.78                                 --> no difference between treatments
pairwise.wilcox.test(sample_data(JE_dBEF)$Simpson, sample_data(JE_dBEF)$history_plant)
# Observed:   0.53
# Shannon:    0.033
# Simpson:    0.0013
# Evenness:   0.00058
pairwise.wilcox.test(sample_data(JE_dBEF)$Evenness, sample_data(JE_dBEF)$history_soil)
# Observed:   0.54
# Shannon:    0.034
# Simpson:    0.0077
# Evenness:   0.0013

## plot alpha diversities
dBEF_rich_plot<- plot_richness(JE_dBEF_rf, x = "LOG_SD", measures=c("Observed"), nrow=3, color="treatment") +
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f" )) +
  geom_point(size=0.4, alpha=0.9) + theme(axis.text.x = element_blank()) +
  facet_grid(.~treatment) + xlab("log(Plant Species Richness)") + ylab("VT Richness") +
  geom_smooth(method="lm")+ theme_bw()

coef(lm(dBEF_samples[dBEF_samples$treatment == "D1",]$Observed~dBEF_samples[dBEF_samples$treatment =="D1",]$LOG_SD))
coef(lm(dBEF_samples[dBEF_samples$treatment == "D2",]$Observed~dBEF_samples[dBEF_samples$treatment =="D2",]$LOG_SD))
coef(lm(dBEF_samples[dBEF_samples$treatment == "D3",]$Observed~dBEF_samples[dBEF_samples$treatment =="D3",]$LOG_SD))
# D1 slope = 15.87 + 2.11   R^2 = 0.07
# D2 slope = 16.10 + 1.48   R^2 = 0.05
# D3 slope = 16.82 + 1.73   R^2 = 0.09
# summary(lm(dBEF_samples[dBEF_samples$treatment == "D3",]$Observed~dBEF_samples[dBEF_samples$treatment =="D3",]$LOG_SD))$r.squared

dBEF_Shannon_plot<- plot_richness(JE_dBEF_rf, x = "LOG_SD", measures=c("Shannon"), nrow=3, color="treatment") + 
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f" )) + 
  geom_point(size=0.4, alpha=0.9) +theme(axis.text.x = element_blank()) +
  facet_grid(.~treatment) + xlab("log(Plant Species Richness)") + ylab("Shannon Diversity") +
  geom_smooth(method = "lm") + theme_bw()

dBEF_Simpson_plot<- plot_richness(JE_dBEF_rf, x = "LOG_SD", measures=c("Simpson"), nrow=3, color="treatment") + 
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) + 
  geom_point(size=0.4, alpha=0.9) + theme_bw() +
  facet_grid(.~treatment) + geom_smooth(method = "lm")+
  xlab("log(Plant Species Richness)") + ylab("Simpson Diversity") 
  
dBEF_Evenn_plot = ggplot(dBEF_samples, aes(x = LOG_SD, y = Evenness, colour= treatment) ) + 
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) +
  geom_point(size=0.6, alpha=0.9) + 
  theme_bw() + geom_smooth(method = "lm") +
  facet_grid(.~treatment) + xlab("log(Plant Species Richness)") + ylab("Pielous Evenness") 
##add slopes with coef(lm(dBEF_samples$Evenness~dBEF_samples$LOG_SD))[2]
coef(lm(dBEF_samples[dBEF_samples$treatment == "D1",]$Evenness~dBEF_samples[dBEF_samples$treatment =="D1",]$LOG_SD))
coef(lm(dBEF_samples[dBEF_samples$treatment == "D2",]$Evenness~dBEF_samples[dBEF_samples$treatment =="D2",]$LOG_SD))
coef(lm(dBEF_samples[dBEF_samples$treatment == "D3",]$Evenness~dBEF_samples[dBEF_samples$treatment =="D3",]$LOG_SD))
# D1 slope = 0.61 + 0.012   R^2 = 0.009
# D2 slope = 0.64 - 0.002   R^2 = 3.689782e-06
# D3 slope = 0.67 + 0.015   R^2 = 0.017
# R^2: summary(lm(dBEF_samples[dBEF_samples$treatment == "D3",]$Evenness~dBEF_samples[dBEF_samples$treatment =="D3",]$LOG_SD))$r.squared

pdf("Graphics/20230504_dBEF_2_alpha.pdf", width = 10, height = 8)
ggarrange(dBEF_rich_plot+ggtitle("A"),
          dBEF_Evenn_plot+ggtitle("B"), 
          dBEF_Shannon_plot+ggtitle("C"), 
          dBEF_Simpson_plot+ggtitle("D"), 
          nrow = 2, ncol = 2,  common.legend = T)
dev.off()

##### 2.2 alpha LM #####
lm_block <- lm(Observed ~ BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
anova(lm_block) 
# F = 5.83, p = 0.0007 ***
lm_blockE <- lm(Evenness ~ BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
anova(lm_blockE) 
# F = 4.67, p = 0.0035 ***
lm_blockS <- lm(Shannon ~ BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
anova(lm_blockS) 
# F = 6.80, p = 0.0002 ***
lm_blockSi <- lm(Simpson ~ BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
anova(lm_blockSi) 

# aov_alpha <- lme(Observed ~ treatment*SOWNDIV, random=~1|BLOCK/plot, data = as(sample_data(JE_dBEF),"data.frame"))
# anova(aov_alpha)
#                    numDF denDF  F-value p-value
# (Intercept)           1   142 346.1176  <.0001
# treatment             2   142   0.9874  0.3751
# SOWNDIV               4    68   5.3853  0.0008 ***
# treatment:SOWNDIV     8   142   0.5612  0.8082
aov_alpha2 <- lme(Observed ~ treatment*LOG_SD, random=~1|BLOCK/plot, data = as(sample_data(JE_dBEF),"data.frame"))
anova(aov_alpha2)
#                   numDF denDF  F-value p-value
# (Intercept)          1   148 349.0883  <.0001
# treatment            2   148   1.0002  0.3703
# LOG_SD               1    71  21.0200  <.0001 ***
# treatment:LOG_SD     2   148   0.1944  0.8235
aov_alphaD1 <- lme(Observed ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D1",])
anova(aov_alphaD1)
#               numDF denDF  F-value p-value
# (Intercept)     1    71 403.6656  <.0001
# LOG_SD          1    71   8.2521  0.0054
aov_alphaD2 <- lme(Observed ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D2",])
anova(aov_alphaD2)
#             numDF denDF   F-value p-value
# (Intercept)     1    71 121.33115  <.0001
# LOG_SD          1    71   6.47219  0.0131
aov_alphaD3 <- lme(Observed ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D3",])
anova(aov_alphaD3)
# n             umDF denDF   F-value p-value
# (Intercept)     1    71 1060.8079  <.0001
# LOG_SD          1    71   10.1335  0.0022

aov_even <- lme(Evenness ~ treatment*LOG_SD, random=~1|BLOCK/plot, data = as(sample_data(JE_dBEF),"data.frame"))
anova(aov_even)
#                   numDF denDF   F-value p-value
# (Intercept)          1   148 1289.0445  <.0001
# treatment            2   148    5.6990  0.0041 **
# LOG_SD               1    71    0.6357  0.4279
# treatment:LOG_SD     2   148    0.3252  0.7229
aov_evenD1 <- lme(Evenness ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D1",])
anova(aov_evenD1)
aov_evenD2 <- lme(Evenness ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D2",])
anova(aov_evenD2)
aov_evenD3 <- lme(Evenness ~ LOG_SD, random=~1|BLOCK/plot, 
                   data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D3",])
anova(aov_evenD3)

aov_Sh <- lme(Shannon ~ treatment*LOG_SD, random=~1|BLOCK/plot, data = as(sample_data(JE_dBEF),"data.frame"))
anova(aov_Sh)
#                   numDF denDF  F-value p-value
# (Intercept)          1   148 522.8727  <.0001
# treatment            2   148   5.4890  0.0050 **
# LOG_SD               1    71   7.3409  0.0084 **
# treatment:LOG_SD     2   148   0.2612  0.7705
aov_ShD1 <- lme(Shannon ~ LOG_SD, random=~1|BLOCK/plot, 
                  data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D1",])
anova(aov_ShD1)
aov_ShD2 <- lme(Shannon ~ LOG_SD, random=~1|BLOCK/plot, 
                  data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D2",])
anova(aov_ShD2)
aov_ShD3 <- lme(Shannon ~ LOG_SD, random=~1|BLOCK/plot, 
                  data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D3",])
anova(aov_ShD3)


aov_Si <- lme(Simpson ~ treatment*LOG_SD, random=~1|BLOCK/plot, data = as(sample_data(JE_dBEF),"data.frame"))
anova(aov_Si)
#                   numDF denDF   F-value p-value
# (Intercept)          1   148 1506.5138  <.0001
# treatment            2   148    4.7822  0.0097 **
# LOG_SD               1    71    2.6358  0.1089
# treatment:LOG_SD     2   148    0.1375  0.8716
aov_SiD1 <- lme(Simpson ~ LOG_SD, random=~1|BLOCK/plot, 
                data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D1",])
anova(aov_SiD1)
aov_SiD2 <- lme(Simpson ~ LOG_SD, random=~1|BLOCK/plot, 
                data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D2",])
anova(aov_SiD2)
aov_SiD3 <- lme(Simpson ~ LOG_SD, random=~1|BLOCK/plot, 
                data = as(sample_data(JE_dBEF),"data.frame")[as(sample_data(JE_dBEF),"data.frame")$treatment == "D3",])
anova(aov_SiD3)

## Do we see patterns when looking at the +/- of plant and soil history individually?
# aov_evenHP <- lme(Observed ~ history_plant*LOG_SD, random=~1|BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
# anova(aov_evenHP)
#                       numDF denDF F-value p-value
# (Intercept)              1   221 349.0145  <.0001
# history_plant            1   221   0.7688  0.3815
# LOG_SD                   1   221  24.9004  <.0001 ***
# history_plant:LOG_SD     1   221   0.0044  0.9474
# aov_alphaHS <- lme(Observed ~ history_soil*LOG_SD, random=~1|BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
# anova(aov_alphaHS)
#                       numDF denDF F-value p-value
# (Intercept)             1   221 349.0151  <.0001
# history_soil            1   221   0.2067  0.6498
# LOG_SD                  1   221  24.8623  <.0001 ***
# history_soil:LOG_SD     1   221   0.2308  0.6314
# aov_evenHS <- lme(Evenness ~ history_soil*LOG_SD, random=~1|BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
# anova(aov_evenHS)
#                     numDF denDF   F-value p-value
# (Intercept)             1   221 1286.3962  <.0001
# history_soil            1   221    4.7168  0.0309
# LOG_SD                  1   221    0.7311  0.3934
# history_soil:LOG_SD     1   221    0.1562  0.6931
# aov_evenHP <- lme(Evenness ~ history_plant*LOG_SD, random=~1|BLOCK, data = as(sample_data(JE_dBEF),"data.frame"))
# anova(aov_evenHP)
#                       numDF denDF   F-value p-value
# (Intercept)              1   221 1285.9754  <.0001
# history_plant            1   221   10.0434  0.0017
# LOG_SD                   1   221    0.7489  0.3878
# history_plant:LOG_SD     1   221    0.1411  0.7075


aov_alphaSDT <- lme(Observed ~ treatment, random=~1|BLOCK, 
                    data = dBEF_samples[dBEF_samples$SOWNDIV==1,])
anova(aov_alphaSDT) #ns for all plant diversities

aov_evenSDT <- lme(Evenness ~ treatment, random=~1|BLOCK, 
                    data = dBEF_samples[dBEF_samples$SOWNDIV==16,])
anova(aov_evenSDT) # 2-species F = 4.09  p = 0.0237 | ns for all other plant diversities


#paired t-test for Alpha Div
library(PairedData)
alpha = group_by(dBEF_samples, plot) %>%
  filter(n() !=1) %>% #filters any unmatched samples
  ungroup()
alphaD1 <- subset(alpha, treatment== "D1", Observed, drop = TRUE)
alphaD2 <- subset(alpha, treatment== "D2", Observed, drop = TRUE)
alphaD3 <- subset(alpha, treatment== "D3", Observed, drop = TRUE)

alphaPairD12 <- paired(alphaD1, alphaD2)
alphaPairD13 <- paired(alphaD1, alphaD3)
alphaPairD23 <- paired(alphaD2, alphaD3)

library(ggpubr)
ggarrange(plot(alphaPairD12, type="profile") + theme_bw(),
          plot(alphaPairD13, type="profile") + theme_bw(),
          plot(alphaPairD23, type="profile") + theme_bw(),
          nrow = 1, ncol = 3
)

#paired t-test of Alpha diversities
alphaD12 <- alpha %>% subset(treatment!= "D3") #%>% select(treatment, Observed)
alphaD23 <- alpha %>% subset(treatment!= "D1")
alphaD13 <- alpha %>% subset(treatment!= "D2")

t.test(Observed~treatment, data=alphaD12, paired=TRUE) # t = 0.6622, df = 75, p-value = 0.5099
adiff12 <- with(alphaD12, 
                Observed[treatment =="D1"] - Observed[treatment=="D2"])
shapiro.test(adiff12) # W = 0.985, p-value = 0.519

t.test(Observed~treatment, data=alphaD13, paired=TRUE) # t = -0.46615, df = 75, p-value = 0.6425
adiff13 <- with(alphaD13, 
                Observed[treatment =="D1"] - Observed[treatment=="D3"])
shapiro.test(adiff13) # W = 0.98678, p-value = 0.6195

t.test(Observed~treatment, data=alphaD23, paired=TRUE) # t = -1.3267, df = 75, p-value = 0.1886
adiff23 <- with(alphaD23, 
                Observed[treatment =="D2"] - Observed[treatment=="D3"])
shapiro.test(adiff23) # W = 0.96885, p-value = 0.05839

t.test(Evenness ~ treatment,data = dBEF_samples[dBEF_samples$treatment != "D3",], paired = TRUE)
# t = -0.4934, df = 75, p-value = 0.6232
t.test(Evenness ~ treatment,data = dBEF_samples[dBEF_samples$treatment != "D2",], paired = TRUE)
# t = -3.3144, df = 75, p-value = 0.001417 --> sig differences Evenness btw D1 and D3
t.test(Evenness ~ treatment,data = dBEF_samples[dBEF_samples$treatment != "D1",], paired = TRUE)
# t = -2.5697, df = 75, p-value = 0.01216    --> sig differences Evenness btw D2 and D3 --> D3 different from other two

mean(dBEF_samples[dBEF_samples$treatment == "D3",]$Evenness)
# 0.6903773
mean(dBEF_samples[dBEF_samples$treatment == "D1",]$Evenness)
# 0.6308699
mean(dBEF_samples[dBEF_samples$treatment == "D2",]$Evenness)
# 0.6408946

##### 2.3 alpha-Diversity trends #####
##filtering out trends of alpha diversity
dBEF_alpha_trends <- dBEF_samples %>% select(plot, treatment, Observed) %>% 
  pivot_wider(names_from = treatment, values_from = Observed)
dBEF_alpha_trends <- dBEF_alpha_trends %>% rowwise() %>% mutate(mean = mean(c(D1,D2,D3)))
dBEF_alpha_trends <- dBEF_alpha_trends %>% rowwise() %>% mutate(D1_nn = mean-D1,
                                                              D2_nn = mean-D2, 
                                                              D3_nn =mean-D3)

winner <-function(a, b, c){
  if (a>b & a>c) 
  {return("decrease")} 
  else if(b>a & b>c) 
  {return("peak")}
  else if(b<a & b<c)
  {return("valley")}
  else
  {return("increase")}  
} 
dBEF_alpha_trends <- dBEF_alpha_trends %>% rowwise() %>% 
  mutate(trend = winner(D1_nn, D2_nn, D3_nn))

dBEF_alpha_trends_long <- dBEF_alpha_trends %>% select(-D1, -D2, -D3, -mean) %>% 
  pivot_longer(cols = c(D1_nn, D2_nn, D3_nn), names_to = "treat", values_to ="Obs")
dBEF_alpha_trends_long <- inner_join(dBEF_alpha_trends_long, dBEF_samples, by ="plot") 

dBEF_alpha_trend_plot <- ggplot(dBEF_alpha_trends_long, aes(x = treat, y = Obs, colour= trend) ) + 
  geom_point() + geom_line(aes(group = plot)) + scale_color_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("VT richness (mean normalized)") + xlab("Treatment")

alpha_trends_bar <- dBEF_alpha_trends_long %>% distinct(plot, .keep_all= TRUE) 
alpha_trends_bar <-  alpha_trends_bar %>% group_by(treat) %>% 
  summarize(increase = sum(str_count(trend, 'increase')), # 15
            valley = sum(str_count(trend, 'valley')),     # 7
            peak = sum(str_count(trend, 'peak')),         # 30
            decrease = sum(str_count(trend, 'decrease'))) # 24
# alpha_trend_leg <- ggplot(dBEF_alpha_trends_long, aes(x = treat, y = Obs, colour= trend) ) + 
#   geom_point() + geom_line(aes(group = plot)) + scale_color_brewer(palette = "Paired") +
#   facet_grid(.~LEG) +  theme(axis.text.x = element_text(angle = 90)) +
#   ylab("VT richness (mean normalized)") + xlab("Treatment")
# alpha_trend_block <- ggplot(dBEF_alpha_trends_long, aes(x = treat, y = Obs, colour= trend) ) + 
#   geom_point() + geom_line(aes(group = plot)) + scale_color_brewer(palette = "Paired") +
#   facet_grid(.~BLOCK) +  theme(axis.text.x = element_text(angle = 90)) +
#   ylab("VT richness (mean normalized)") + xlab("Treatment")
pdf("Graphics/20230504_dBEF_2_alpha_trends.pdf")
dBEF_alpha_trend_plot + facet_grid(.~SOWNDIV)
#;alpha_trend_leg; alpha_trend_block
dev.off()



#### 3 PARAFAC as explorative analysis #### 
# for explorative PARAFAC analysis we exported the filtered, rarefied and CLR-transformed Vt abundance
# data to MATLAB to use package 'N-way toolbox' with a bootstrap approach of 100 initiations for the best
# model.

##### 3.1 load PARAFAC results #####
# Output of all 100 initations/models (1 component / 2 components) was loaded back into R:
model_folder1 <- "C:/Users/albracht/Documents/MATLAB/PARAFAC/test2_v1_Jena/comp.1"
model_folder2 <- "C:/Users/albracht/Documents/MATLAB/PARAFAC/test2_v1_Jena/comp.2" # 1 component is better model

summarize_parafac <- function(directory){
  input_files <- list.files(path=directory,
                            pattern="_input.csv")
  feature_files <- list.files(path=directory,
                              pattern="_feature_mode.csv")
  treatment_files <- list.files(path=directory,
                                pattern="_time_mode.csv")
  plot_files <- list.files(path=directory,
                           pattern="_plots_mode.csv")
  
  compNo <- unique(gsub("model_comp.([[:digit:]]+)_model.[[:digit:]]+_input.csv",
                        "\\1",input_files))
  if(length(compNo)>1) stop("More than one model in folder.") else 
    print(paste(compNo,"component model."))
  
  compNum <- as.numeric(compNo)
  
  inits <- gsub("model_comp.[[:digit:]]+_model.([[:digit:]]+)_input.csv",
                "\\1",input_files)
  print(paste(length(inits),"initializations found."))
  
  features <- lapply(inits,
                     function(i) read.delim(file.path(directory,
                                                      paste0("model_comp.",
                                                             compNo,"_model.",
                                                             i,"_feature_mode.csv")),
                                            sep=",",header=F))
  if(!all(sapply(2:length(features), 
                 function(i) all(features[[1]][,compNum+1]==features[[i]][,compNum+1])))) 
    stop("features don't match")
  
  treatments <- lapply(inits,
                       function(i) read.delim(file.path(directory,
                                                        paste0("model_comp.",
                                                               compNo,"_model.",
                                                               i,"_time_mode.csv")),
                                              sep=",",header=F))
  if(!all(sapply(2:length(treatments), 
                 function(i) all(dim(treatments[[1]])==dim(treatments[[i]]))))) 
    stop("treatments don't match")
  
  plots <- lapply(inits,
                  function(i) read.delim(file.path(directory,
                                                   paste0("model_comp.",
                                                          compNo,"_model.",
                                                          i,"_plots_mode.csv")),
                                         sep=",",header=F))
  if(!all(sapply(2:length(plots), 
                 function(i) all(plots[[1]][,compNum+1]==plots[[i]][,compNum+1])))) 
    stop("plots don't match")
  
  treatment_summary <- lapply(1:compNum,
                              function(comp) data.frame("treatment"=c("D1","D2","D3"),
                                                        "mean"=rowMeans(sapply(treatments,
                                                                               function(i) i[,comp])),
                                                        "sd"=apply(sapply(treatments,
                                                                          function(i) i[,comp]),
                                                                   1,sd),
                                                        "min"=apply(sapply(treatments,
                                                                           function(i) i[,comp]),
                                                                    1,min),
                                                        "max"=apply(sapply(treatments,
                                                                           function(i) i[,comp]),
                                                                    1,max)))
  treatment_mean <- data.frame("treatment"=c("D1","D2","D3"),
                               sapply(1:compNum,function(comp) treatment_summary[[comp]]$mean))
  colnames(treatment_mean) <- c("treatment",paste0("component.",1:compNum))
  
  feature_summary <- lapply(1:compNum,
                            function(comp) data.frame(features[[1]][,compNum+1:3],
                                                      "mean"=rowMeans(sapply(features,
                                                                             function(i) i[,comp])),
                                                      "sd"=apply(sapply(features,
                                                                        function(i) i[,comp]),
                                                                 1,sd),
                                                      "min"=apply(sapply(features,
                                                                         function(i) i[,comp]),
                                                                  1,min),
                                                      "max"=apply(sapply(features,
                                                                         function(i) i[,comp]),
                                                                  1,max)))
  feature_mean <- data.frame(feature_summary[[1]][,1:3],
                             sapply(1:compNum,function(comp) feature_summary[[comp]]$mean))
  colnames(feature_mean) <- c(colnames(feature_mean)[1:3],
                              paste0("component.",1:compNum))
  
  plot_summary <- lapply(1:compNum,
                         function(comp) data.frame(plots[[1]][,compNum+1:3],
                                                   "mean"=rowMeans(sapply(plots,
                                                                          function(i) i[,comp])),
                                                   "sd"=apply(sapply(plots,
                                                                     function(i) i[,comp]),
                                                              1,sd),
                                                   "min"=apply(sapply(plots,
                                                                      function(i) i[,comp]),
                                                               1,min),
                                                   "max"=apply(sapply(plots,
                                                                      function(i) i[,comp]),
                                                               1,max)))
  plot_mean <- data.frame(plot_summary[[1]][,1:3],
                          sapply(1:compNum,function(comp) plot_summary[[comp]]$mean))
  colnames(plot_mean) <- c(colnames(plot_mean)[1:3],
                           paste0("component.",1:compNum))
  
  return(list("plotMeans"=plot_mean,
              "treatmentMeans"=treatment_mean,
              "featureMeans"=feature_mean,
              "plotSummary"=plot_summary,
              "treatmentSummary"=treatment_summary,
              "featureSummary"=feature_summary))
}

##### 3.2 summarize PARAFAC models #####
# 1 component models
# plot_means: V2 = plot / V3 = SONWDIV / V4 = Block
summarizedModels1 <- summarize_parafac(model_folder1)
plot_means1 <- summarizedModels1$plotMeans
feature_means1 <- summarizedModels1$featureMeans
treatment_means1 <- summarizedModels1$treatmentMeans
boxplot(plot_means1$component.1~plot_means1$V4)
kruskal.test(plot_means1$component.1~plot_means1$V4)
boxplot(plot_means1$component.1~plot_means1$V3) #chi^2 = 5.6942, df = 3, p-value = 0.1275
kruskal.test(plot_means1$component.1~plot_means1$V3) #chi^2 = 35.808 df = 4 p = 3.169e-07
# feature means (VT): V2 = VTX / V3 = Family / V4 = Order
boxplot(feature_means1$component.1~feature_means1$V4,las=2)
boxplot(feature_means1$component.1~feature_means1$V3,las=2)
boxplot(feature_means1$component.1~feature_means1$V2,las=2)
boxplot(feature_means1$component.1~feature_means1$V4,las=2)
plot(treatment_means1$component.1)
#text(treatment_means1$component.1,
#     labels=treatment_means1$treatment)

# # 2 component models
# summarizedModels2 <- summarize_parafac(model_folder2)
# # plot_means: V3 = plot / V4 = SONWDIV / V5 = Block
# plot_means2 <- summarizedModels2$plotMeans
# feature_means2 <- summarizedModels2$featureMeans
# treatment_means2 <- summarizedModels2$treatmentMeans
# 
# boxplot(plot_means2$component.1~plot_means2$V5)
# kruskal.test(plot_means2$component.1~plot_means2$V5) #chi^2 = 5.5803 df = 3 p = 0.1339
# 
# boxplot(plot_means2$component.2~plot_means2$V5)
# kruskal.test(plot_means2$component.2~plot_means2$V5) #chi^2 = 10.786 df = 3 p = 0.01294
# 
# boxplot(plot_means2$component.1~plot_means2$V4)
# kruskal.test(plot_means2$component.1~plot_means2$V4) #chi^2 = 38.771 df = 4 p = 7.768e-08
# boxplot(plot_means2$component.2~plot_means2$V4)
# kruskal.test(plot_means2$component.2~plot_means2$V4) #chi^2 = 7.5498 df = 4 p = 0.1095
# 
# # feature means (VT): V3 = VTX / V4 = Family / V5 = Order
# boxplot(feature_means2$component.1~feature_means2$V4,las=2)
# boxplot(feature_means2$component.1~feature_means2$V5,las=2)
# boxplot(feature_means2$component.2~feature_means2$V4,las=2)
# boxplot(feature_means2$component.2~feature_means2$V5,las=2)
# 
# plot(treatment_means2$component.1,
#      treatment_means2$component.2,type="n")
# text(treatment_means2$component.1,
#      treatment_means2$component.2,
#      labels=treatment_means2$treatment)

##### 3.3 plot PARAFAC #####
parafac1_treat = summarizedModels1$treatmentMeans %>% 
  ggplot(aes(x = treatment, y = component.1, color = treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) +
  theme_minimal() +
  theme(legend.position="none", plot.title = element_text(size=11)) + 
  geom_hline(aes(yintercept = 0, alpha=0.7), linetype = "dashed") + 
  ggtitle("C") + ylim(0,1) +
  xlab("Treatment") + ylab("loading means")

parafac1_plots = summarizedModels1$plotMeans %>%
  ggplot(aes(x=as.factor(V3), y=component.1, fill = as.factor(V3))) +
  geom_boxplot() + 
  #scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) +
  scale_fill_brewer(palette = "Greys") +
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_minimal() +
  theme(legend.position="none", plot.title = element_text(size=11)) + 
  geom_hline(aes(yintercept = 0, alpha=0.7), linetype = "dashed") + 
  ggtitle("B") + ylim(-8,8) +
  xlab("Plant Species Richness") + ylab("loading means")

parafac1_VT = summarizedModels1$featureMeans %>%
  ggplot(aes(x=reorder(V2, -component.1, sum), y=component.1, fill = V4)) +
  geom_col() + 
  #scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) +
  #scale_fill_manual(values= moregrey) +
  scale_fill_paletteer_d("ggsci::springfield_simpsons") +
  theme_minimal() +
  theme(legend.position="bottom", plot.title = element_text(size=11),
        axis.text.x = element_blank()) +
  geom_hline(aes(yintercept = 0, alpha=0.7), linetype = "dashed") + 
  ggtitle("A") + ylim(-0.5, 0.5)+
  xlab("Genus / Cluster") + ylab("loading means")

summarizedModels1$featureMeans %>%
  ggplot(aes(x=reorder(V2, -component.1, sum), y=component.1, fill = V4)) +
  geom_col() + 
  #scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) +
  #scale_fill_manual(values= moregrey) +
  scale_fill_paletteer_d("ggsci::springfield_simpsons") +
  theme_minimal() +
  theme(legend.position="bottom", plot.title = element_text(size=11)
        ) +
  geom_hline(aes(yintercept = 0, alpha=0.7), linetype = "dashed") + 
  ggtitle("A") + ylim(-0.4, 0.4)+
  xlab("VT") + ylab("loading means") + coord_flip()

pdf("Graphics/20230504_dBEF_3_parafac.pdf", height = 6, width = 10)
ggarrange(parafac1_VT, #parafac2_treat, 
          parafac1_plots, #parafac2_plots, 
          parafac1_treat,   #parafac2_VT, 
          nrow = 1, ncol =  3, common.legend = T)
parafac1_VT #for legend only
dev.off()

##### 3.4 PCA #####
# dBEF_pca = prcomp(t(VT_dBEF_frf_clr), scale. = TRUE)
# dBEF_pca_plot <- autoplot(dBEF_pca, data = dBEF_samples, colour = "treatment", label = F, shape = "SOWNDIV",
#          frame = TRUE, frame.type = 'norm')
#          #loadings = T, loadings.label = TRUE, loadings.label.size = 2)

library("mixOmics")
dBEF_pca <- pca(t(VT_dBEF_frf_clr), ncomp = 2, center = TRUE, scale = TRUE)
plot(dBEF_pca) # explained variance per component
# pcaVar <- tune.pca(t(VT_dBEF_frf_clr), ncomp = 10, center = T, scale = F)
# plot(pcaVar) # diff way to plot expl variance
# keepX <- c(seq(5, 25, 5))
# plot(pcatune) #diamonds = optimal number of features er component
# pcatune$choice.keepX #25 10 5
dBEF_samples$treatSD <- paste0(dBEF_samples$treatment,"_", dBEF_samples$SOWNDIV)
dBEF_pca_plot <- plotIndiv(dBEF_pca, comp = c(1,2), ind.names = F,
                       group = dBEF_samples$treatment, legend = T, legend.title = "Treatment",
                       pch = dBEF_samples$SOWNDIV, legend.title.pch = "Plant Species Richness",
                       cex = 2, ellipse = T,
                       col.per.group = c("#01665e", "#91cf60", "#fec44f"))
# pcaVar_plot <- plotVar(dBEF_pca, comp = c(1,1), var.names = T, title = "PCA Var")

# sparse PCA to focus on important VT
spcatune <- tune.spca(t(VT_dBEF_frf_clr), ncomp = 3, nrepeat = 5, folds = 3, test.keepX = keepX)
dBEF_spca <- spca(t(VT_dBEF_frf_clr), ncomp = 3, keepX = spcatune$choice.keepX)
dBEF_spca_plot <- plotIndiv(dBEF_spca, comp = c(1,2), ind.names = TRUE,
                       group = dBEF_samples$treatment, legend = T,
                       ellipse = T, col.per.group = c("#01665e", "#91cf60", "#fec44f"))
dBEF_spca_Varplot <- plotVar(dBEF_spca, comp = c(1,2), var.names = T, title = "sPCA Var") #shows cluster of VT and their contribution to each principal component
plotLoadings(dBEF_spca)

pdf("Graphics/20230504_dBEF_3_pca.pdf", width = 6, height = 4)
plotIndiv(dBEF_pca, comp = c(1,2), ind.names = F,
          group = dBEF_samples$treatment, legend = T, legend.title = "Treatment",
          pch = dBEF_samples$SOWNDIV, legend.title.pch = "Plant Species Richness",
          cex = 2, ellipse = T,
          col.per.group = c("#01665e", "#91cf60", "#fec44f"))
biplot(dBEF_spca, cex = 0.7, ind.names = F, 
       group = dBEF_samples$treatment, legend.title = "Treatment", #expand = 10, 
       col.per.group = c("#01665e", "#91cf60", "#fec44f"), ellipse = T,
       pch = dBEF_samples$SOWNDIV, pch.legend.title ="Plant Species Richness", pch.size = 1,
       var.names = T, var.names.size = 3) #,
       #xlim = c(-7,3), ylim = c(-3,7))  #shows how each feature may explain the positioning of the sample
dev.off()
## --> Glomus in direction of D1, Claroideoglomus drive D3 (and D2?)


#### 4 beta-Diversity and permanova ####
##### 4.1 Aitchison distances: extract pairwise comparisons  #####
exJE_ps <- as.data.frame(otu_table(JE_dBEF_f)) # CLR creates NA, so we can only use the filtered+rarefied abundances
exJE_ps_pa <- decostand(exJE_ps, method = "pa") #for test with Soerenson distance
exJE_ps[exJE_ps==0] <- 0.1
exJE_ps <- as.matrix(exJE_ps)

sI <- dBEF_samples[order(rownames(dBEF_samples)),]
sI$sample <- paste0(sI$plot,sI$treatment)

library(robCompositions)
exJE_dist <- aDist(t(exJE_ps)) 
exJE_distM <- as.matrix(exJE_dist)

# check that my metadata makes sense:
all(colnames(exJE_distM)==sI$sample)
all(rownames(exJE_distM)==sI$sample)

# extract groups of interest (from upper triangle):
iMat <- array(list(),dim=c(15,15), # new data will go into an array with lists of values
              dimnames=list(paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                 each=length(unique(sI$treatment))),sep=""),
                                  rep(sort(unique(sI$treatment)),
                                      times=length(unique(sI$SOWNDIV))),sep="_"),
                            paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                 each=length(unique(sI$treatment))),sep=""),
                                  rep(sort(unique(sI$treatment)),
                                      times=length(unique(sI$SOWNDIV))),sep="_")))
# long table version:
iTab <- data.frame("sample1"="","sample2"="","group1"="","group2"="",
                   "common"="",
                   "distance"=rep(0,length(c(exJE_distM))/2-(nrow(exJE_distM)/2)))
cnt <- 1
for(i in 1:(nrow(exJE_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(exJE_distM)){ #j are columns - only read upper triangle
    comp1 <- paste(paste("SD",sI$SOWNDIV[i],sep=""),sI$treatment[i],sep="_")
    comp2 <- paste(paste("SD",sI$SOWNDIV[j],sep=""),sI$treatment[j],sep="_")
    fields <- sort(c(which(rownames(iMat)==comp1),which(rownames(iMat)==comp2))) # only fill upper triangle
    iMat[fields[1],fields[2]] <- list(append(iMat[fields[1],fields[2]][[1]],
                                             exJE_distM[i,j]))
    iTab$distance[cnt] <- exJE_distM[i,j]
    iTab$sample1[cnt] <- rownames(exJE_distM)[i]
    iTab$sample2[cnt] <- colnames(exJE_distM)[j]
    iTab$group1[cnt] <- comp1
    iTab$group2[cnt] <- comp2
    coms <- ""
    if(sI$SOWNDIV[i]==sI$SOWNDIV[j]) coms <- paste(coms,paste0("SD",sI$SOWNDIV[j]),sep="_")
    if(sI$treatment[i]==sI$treatment[j]) coms <- paste(coms,sI$treatment[j],sep="_")
    if(sI$plot[i]==sI$plot[j]) coms <- paste(coms,sI$plot[j],sep="_")
    coms <- sub("^_","",coms)
    iTab$common[cnt] <- coms
    cnt <- cnt + 1
  }
}


iMatP <- array(list(),dim=c(15,15), # new data will go into an array with lists of values
               dimnames=list(paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                  each=length(unique(sI$treatment))),sep=""),
                                   rep(sort(unique(sI$treatment)),
                                       times=length(unique(sI$SOWNDIV))),sep="_"),
                             paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                  each=length(unique(sI$treatment))),sep=""),
                                   rep(sort(unique(sI$treatment)),
                                       times=length(unique(sI$SOWNDIV))),sep="_")))

for(i in 1:(nrow(exJE_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(exJE_distM)){ #j are columns - only read upper triangle
    # for comparisons of D1-D3, use only samples from the same plot, keep other fields empty:
    if(sI$plot[i]==sI$plot[j]){
      comp1 <- paste(paste("SD",sI$SOWNDIV[i],sep=""),sI$treatment[i],sep="_")
      comp2 <- paste(paste("SD",sI$SOWNDIV[j],sep=""),sI$treatment[j],sep="_")
      fields <- sort(c(which(rownames(iMatP)==comp1),which(rownames(iMatP)==comp2))) # only fill upper triangle
      iMatP[fields[1],fields[2]] <- list(append(iMatP[fields[1],fields[2]][[1]],
                                                exJE_distM[i,j]))
    }
  }
}

##### 4.2 Aitchison distances: test pairwise comparisons #####
# are old monocultures more similar to each other than new ones?
wilcox.test(iMat[3,3][[1]],
            iMat[1,1][[1]],paired = T)
# yes  p = 9.047e-06
wilcox.test(iMat[2,2][[1]],
            iMat[1,1][[1]],paired = T)
# no p = 0.7742
wilcox.test(iMat[3,3][[1]],
            iMat[2,2][[1]],paired = T)
# yes p = 0.0006135

# are old high-diverse mixtures more similar to each other than new ones?
wilcox.test(iMat[15,15][[1]],
            iMat[13,13][[1]],paired = T)
#yes p = 0.03306
wilcox.test(iMat[14,14][[1]],
            iMat[13,13][[1]],paired = T)
#no p = 0.7681
wilcox.test(iMat[15,15][[1]],
            iMat[14,14][[1]],paired = T)
#no p = 0.1311

#in the middle:
# 2-species mixes
wilcox.test(iMat[6,6][[1]],
            iMat[4,4][[1]],paired = T)
# no (means of D3 are higher than of D1): p = 0.3925
wilcox.test(iMat[5,5][[1]],
            iMat[4,4][[1]],paired = T)
# yes p = 0.03722
wilcox.test(iMat[6,6][[1]],
            iMat[5,5][[1]],paired = T)
# yes p = 0.003918

# 4-species mixes
wilcox.test(iMat[9,9][[1]],
            iMat[7,7][[1]],paired = T)
# no p = 0.3983
wilcox.test(iMat[8,8][[1]],
            iMat[7,7][[1]],paired = T)
# no p = 0.1873
wilcox.test(iMat[9,9][[1]],
            iMat[8,8][[1]],paired = T)
# no p = 0.5186

# 8-species mixes
wilcox.test(iMat[12,12][[1]],
            iMat[10,10][[1]],paired = T)
# no p = 0.8948
wilcox.test(iMat[11,11][[1]],
            iMat[10,10][[1]],paired = T)
# no p = 0.469
wilcox.test(iMat[12,12][[1]],
            iMat[11,11][[1]],paired = T)
# no p = 0.3587

# are old monoculture plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[1,3][[1]], #D3 vs D1
            iMatP[2,3][[1]]) #D3 vs D2
# (yes) p = 0.04974

# are no-history monoculture plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[1,2][[1]], #D1 vs D2
            iMatP[1,3][[1]]) #D1 vs D3
# no p = 0.4274

# are old high-diversity plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[13,15][[1]], #D3 vs D1
            iMatP[14,15][[1]]) #D3 vs D2
# yes p = 0.004906

# are no-history high-diversity plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[13,14][[1]], #D1 vs D2
            iMatP[13,15][[1]]) #D1 vs D3
# no p = 0.91

# are old 2-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[4,6][[1]], #D3 vs D1
            iMatP[5,6][[1]]) #D3 vs D2
# (yes) p = 0.01356

# are no-history 2-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[4,5][[1]], #D1 vs D2
            iMatP[4,6][[1]]) #D1 vs D3
# no p = 0.7804

# are old 4-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[7,9][[1]], #D3 vs D1
            iMatP[8,9][[1]]) #D3 vs D2
# yes p = 0.00019

# are no-history 4-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[7,8][[1]], #D1 vs D2
            iMatP[7,9][[1]]) #D1 vs D3
# no p = 0.3045

# are old 8-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[10,12][[1]], #D3 vs D1
            iMatP[11,12][[1]]) #D3 vs D2
# yes p = 0.002593

# are no-history 8-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[10,11][[1]], #D1 vs D2
            iMatP[10,12][[1]]) #D1 vs D3
# no p = 0.8672



# are monocultures more dissimilar from each other than high-diversity plots (old treatment)?
wilcox.test(iMat[3,3][[1]],
            iMat[15,15][[1]])
#no p = 0.704

# are monocultures more dissimilar from each other than high-diversity plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[13,13][[1]])
# (yes) p = 0.03835

# are monocultures more dissimilar from each other than 2-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[4,4][[1]])
# yes p = 0.0002782

# are monocultures more dissimilar from each other than 4-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[7,7][[1]])
# no p = 0.5216

# are monocultures more dissimilar from each other than 8-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[10,10][[1]])
# yes p = 0.0004263

# are 2-species mixes more dissimilar from each other than 4-species plots (new treatment)?
wilcox.test(iMat[7,7][[1]],
            iMat[4,4][[1]])
# yes (2 are more dissimilar) p = 0.01241

# are monocultures more dissimilar from each other than high-diversity plots (treatment D2)?
wilcox.test(iMat[2,2][[1]],
            iMat[14,14][[1]])
# (yes) (boarderline) p = 0.04886



##### 4.3 Aitchison distances: plot pairwise comparisons ##### 
limD <- c(min(unlist(iMat)),max(unlist(iMat)))
pdf("Graphics/20230614_dBEF_4_betabox_withinT.pdf")
boxplot(iMat[1,1][[1]], #within SD1_D1
        iMat[2,2][[1]], #within SD1_D2
        iMat[3,3][[1]], #within SD1_D3
        names=c(paste(rownames(iMat)[1],"\n",colnames(iMat)[1],sep=""), #within SD1_D3
                paste(rownames(iMat)[2],"\n",colnames(iMat)[2],sep=""), #within SD1_D2
                paste(rownames(iMat)[3],"\n",colnames(iMat)[3],sep="") #within SD1_D1
        ),
        main="SD1",ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[1,1][[1]]))), factor = 6), iMat[1,1][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[2,2][[1]]))), factor = 3), iMat[2,2][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[3,3][[1]]))), factor = 2), iMat[3,3][[1]], cex = 0.7)

boxplot(iMat[4,4][[1]], #within SD2_D´1
        iMat[5,5][[1]], #within SD2_D2
        iMat[6,6][[1]], #within SD2_D3
        names=c(paste(rownames(iMat)[4],"\n",colnames(iMat)[4],sep=""), #within SD2_D1
                paste(rownames(iMat)[5],"\n",colnames(iMat)[5],sep=""), #within SD2_D2
                paste(rownames(iMat)[6],"\n",colnames(iMat)[6],sep="") #within SD2_D3
        ),
        main="SD2",ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[4,4][[1]]))), factor = 6), iMat[4,4][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[5,5][[1]]))), factor = 3), iMat[5,5][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[6,6][[1]]))), factor = 2), iMat[6,6][[1]], cex = 0.7)

boxplot(iMat[7,7][[1]], #within SD4_D1
        iMat[8,8][[1]], #within SD4_D2
        iMat[9,9][[1]], #within SD4_D3
        names=c(paste(rownames(iMat)[7],"\n",colnames(iMat)[7],sep=""), #within SD4_D1
                paste(rownames(iMat)[8],"\n",colnames(iMat)[8],sep=""), #within SD4_D2
                paste(rownames(iMat)[9],"\n",colnames(iMat)[9],sep="") #within SD4_D3
        ),
        main="SD4",ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[7,7][[1]]))), factor = 6), iMat[7,7][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[8,8][[1]]))), factor = 3), iMat[8,8][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[9,9][[1]]))), factor = 2), iMat[9,9][[1]], cex = 0.7)

boxplot(iMat[10,10][[1]], #within SD8_D1
        iMat[11,11][[1]], #within SD8_D2
        iMat[12,12][[1]], #within SD8_D3
        names=c(paste(rownames(iMat)[10],"\n",colnames(iMat)[10],sep=""), #within SD8_D1
                paste(rownames(iMat)[11],"\n",colnames(iMat)[11],sep=""), #within SD8_D2
                paste(rownames(iMat)[12],"\n",colnames(iMat)[12],sep="") #within SD8_D3
        ),
        main="SD8",ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[10,10][[1]]))), factor = 6), iMat[10,10][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[11,11][[1]]))), factor = 3), iMat[11,11][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[12,12][[1]]))), factor = 2), iMat[12,12][[1]], cex = 0.7)

boxplot(iMat[13,13][[1]], #within SD16_D1
        iMat[14,14][[1]], #within SD16_D2
        iMat[15,15][[1]], #within SD16_D3
        names=c(paste(rownames(iMat)[13],"\n",colnames(iMat)[13],sep=""), #within SD16_D1
                paste(rownames(iMat)[14],"\n",colnames(iMat)[14],sep=""), #within SD16_D2
                paste(rownames(iMat)[15],"\n",colnames(iMat)[15],sep="") #within SD16_D3
        ),
        main="SD16",ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[13,13][[1]]))), factor = 6), iMat[13,13][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[14,14][[1]]))), factor = 3), iMat[14,14][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[15,15][[1]]))), factor = 2), iMat[15,15][[1]], cex = 0.7)
dev.off()

boxplot(iMat[1,1][[1]], #within SD1_D1
        iMat[2,2][[1]], #within SD1_D2
        iMat[3,3][[1]], #within SD1_D3
        iMat[4,4][[1]], #within SD2_D1
        iMat[5,5][[1]], #within SD2_D2
        iMat[6,6][[1]], #within SD2_D3
        iMat[7,7][[1]], #within SD4_D1
        iMat[8,8][[1]], #within SD4_D2
        iMat[9,9][[1]], #within SD4_D3
        iMat[10,10][[1]], #within SD8_D1
        iMat[11,11][[1]], #within SD8_D2
        iMat[12,12][[1]], #within SD8_D3
        iMat[13,13][[1]], #within SD16_D1
        iMat[14,14][[1]], #within SD16_D2
        iMat[15,15][[1]], #within SD16_D3
        names=c(paste(rownames(iMat)[1],"\n",colnames(iMat)[1],sep=""), #within SD1_D1
                paste(rownames(iMat)[2],"\n",colnames(iMat)[2],sep=""), #within SD1_D2
                paste(rownames(iMat)[3],"\n",colnames(iMat)[3],sep=""), #within SD1_D3
                paste(rownames(iMat)[4],"\n",colnames(iMat)[4],sep=""), #within SD2_D1
                paste(rownames(iMat)[5],"\n",colnames(iMat)[5],sep=""), #within SD2_D2
                paste(rownames(iMat)[6],"\n",colnames(iMat)[6],sep=""), #within SD2_D3
                paste(rownames(iMat)[7],"\n",colnames(iMat)[7],sep=""), #within SD4_D1
                paste(rownames(iMat)[8],"\n",colnames(iMat)[8],sep=""), #within SD4_D2
                paste(rownames(iMat)[9],"\n",colnames(iMat)[9],sep=""), #within SD4_D3
                paste(rownames(iMat)[10],"\n",colnames(iMat)[10],sep=""), #within SD8_D1
                paste(rownames(iMat)[11],"\n",colnames(iMat)[11],sep=""), #within SD8_D2
                paste(rownames(iMat)[12],"\n",colnames(iMat)[12],sep=""), #within SD8_D3
                paste(rownames(iMat)[13],"\n",colnames(iMat)[13],sep=""), #within SD16_D1
                paste(rownames(iMat)[14],"\n",colnames(iMat)[14],sep=""), #within SD16_D2
                paste(rownames(iMat)[15],"\n",colnames(iMat)[15],sep="") #within SD16_D3
        ),
        ylim=limD,
        las=2)

pdf("Graphics/20230614_dBEF_4_betabox_betweenT.pdf")
boxplot(iMatP[1,3][[1]], #between SD1_D1 and SD1_D3
        iMatP[2,3][[1]], #between SD1_D2 and SD1_D3
        iMatP[1,2][[1]], #between SD1_D1 and SD1_D2
        names=c(paste(rownames(iMatP)[1],"\n",colnames(iMatP)[3],sep=""), #between SD1_D1 and SD1_D3
                paste(rownames(iMatP)[2],"\n",colnames(iMatP)[3],sep=""), #between SD1_D2 and SD1_D3
                paste(rownames(iMatP)[1],"\n",colnames(iMatP)[2],sep="") #between SD1_D1 and SD1_D2
        ),
        ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[1,3][[1]]))), factor = 6), iMat[1,3][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[2,3][[1]]))), factor = 3), iMat[2,3][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[1,2][[1]]))), factor = 2), iMat[1,2][[1]], cex = 0.7)

boxplot(iMatP[4,6][[1]], #between SD2_D1 and SD2_D3
        iMatP[5,6][[1]], #between SD2_D2 and SD2_D3
        iMatP[4,5][[1]], #between SD2_D1 and SD2_D2
        names=c(paste(rownames(iMatP)[4],"\n",colnames(iMatP)[6],sep=""), #between SD2_D1 and SD2_D3
                paste(rownames(iMatP)[5],"\n",colnames(iMatP)[6],sep=""), #between SD2_D2 and SD2_D3
                paste(rownames(iMatP)[4],"\n",colnames(iMatP)[5],sep="") #between SD2_D1 and SD2_D2
        ),
        ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[4,6][[1]]))), factor = 6), iMat[4,6][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[5,6][[1]]))), factor = 3), iMat[5,6][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[4,5][[1]]))), factor = 2), iMat[4,5][[1]], cex = 0.7)

boxplot(iMatP[7,9][[1]], #between SD4_D1 and SD4_D3
        iMatP[8,9][[1]], #between SD4_D2 and SD4_D3
        iMatP[7,8][[1]], #between SD4_D1 and SD4_D2
        names=c(paste(rownames(iMatP)[7],"\n",colnames(iMatP)[9],sep=""), #between SD4_D1 and SD4_D3
                paste(rownames(iMatP)[8],"\n",colnames(iMatP)[9],sep=""), #between SD4_D2 and SD4_D3
                paste(rownames(iMatP)[7],"\n",colnames(iMatP)[8],sep="") #between SD4_D1 and SD4_D2
        ),
        ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[7,9][[1]]))), factor = 6), iMat[7,9][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[8,9][[1]]))), factor = 3), iMat[8,9][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[7,8][[1]]))), factor = 2), iMat[7,8][[1]], cex = 0.7)

boxplot(iMatP[10,12][[1]], #between SD8_D1 and SD8_D3
        iMatP[11,12][[1]], #between SD8_D2 and SD8_D3
        iMatP[10,11][[1]], #between SD8_D1 and SD8_D2
        names=c(paste(rownames(iMatP)[10],"\n",colnames(iMatP)[12],sep=""), #between SD8_D1 and SD8_D3
                paste(rownames(iMatP)[11],"\n",colnames(iMatP)[12],sep=""), #between SD8_D2 and SD8_D3
                paste(rownames(iMatP)[10],"\n",colnames(iMatP)[11],sep="") #between SD8_D1 and SD8_D2
        ),
        ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[10,12][[1]]))), factor = 6), iMat[10,12][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[11,12][[1]]))), factor = 3), iMat[11,12][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[10,11][[1]]))), factor = 2), iMat[10,11][[1]], cex = 0.7)

boxplot(iMatP[13,15][[1]], #between SD16_D1 and SD16_D3
        iMatP[14,15][[1]], #between SD16_D2 and SD16_D3
        iMatP[13,14][[1]], #between SD16_D1 and SD16_D2
        names=c(paste(rownames(iMatP)[13],"\n",colnames(iMatP)[15],sep=""), #between SD16_D1 and SD16_D3
                paste(rownames(iMatP)[14],"\n",colnames(iMatP)[15],sep=""), #between SD16_D2 and SD16_D3
                paste(rownames(iMatP)[13],"\n",colnames(iMatP)[14],sep="") #between SD16_D1 and SD16_D2
        ),
        ylim=limD,
        las=2)
points(x = jitter(c(rep(1, length(iMat[13,15][[1]]))), factor = 6), iMat[13,15][[1]], cex = 0.7)
points(x = jitter(c(rep(2, length(iMat[14,15][[1]]))), factor = 3), iMat[14,15][[1]], cex = 0.7)
points(x = jitter(c(rep(3, length(iMat[13,14][[1]]))), factor = 2), iMat[13,14][[1]], cex = 0.7)
dev.off()


boxplot(iMat[3,3][[1]], #within SD1_D3
        iMat[6,6][[1]], #within SD2_D3
        iMat[9,9][[1]], #within SD4_D3
        iMat[12,12][[1]], #within SD8_D3
        iMat[15,15][[1]], #within SD16_D3
        names=c(paste(rownames(iMat)[3],"\n",colnames(iMat)[3],sep=""), #within SD1_D3
                paste(rownames(iMat)[6],"\n",colnames(iMat)[6],sep=""), #within SD2_D3
                paste(rownames(iMat)[9],"\n",colnames(iMat)[9],sep=""), #within SD4_D3
                paste(rownames(iMat)[12],"\n",colnames(iMat)[12],sep=""), #within SD8_D3
                paste(rownames(iMat)[15],"\n",colnames(iMat)[15],sep="") #within SD16_D3
        ),
        main="D3",ylim=limD,
        las=2)


boxplot(iMat[2,2][[1]], #within SD1_D2
        iMat[5,5][[1]], #within SD2_D2
        iMat[8,8][[1]], #within SD4_D2
        iMat[11,11][[1]], #within SD8_D2
        iMat[14,14][[1]], #within SD16_D2
        names=c(paste(rownames(iMat)[2],"\n",colnames(iMat)[2],sep=""), #within SD1_D2
                paste(rownames(iMat)[5],"\n",colnames(iMat)[5],sep=""), #within SD2_D2
                paste(rownames(iMat)[8],"\n",colnames(iMat)[8],sep=""), #within SD4_D2
                paste(rownames(iMat)[11],"\n",colnames(iMat)[11],sep=""), #within SD8_D2
                paste(rownames(iMat)[14],"\n",colnames(iMat)[14],sep="") #within SD16_D2
        ),
        main="D2",ylim=limD,
        las=2)

boxplot(iMat[1,1][[1]], #within SD1_D1
        iMat[4,4][[1]], #within SD2_D1
        iMat[7,7][[1]], #within SD4_D1
        iMat[10,10][[1]], #within SD8_D1
        iMat[13,13][[1]], #within SD16_D1
        names=c(paste(rownames(iMat)[1],"\n",colnames(iMat)[1],sep=""), #within SD1_D1
                paste(rownames(iMat)[4],"\n",colnames(iMat)[4],sep=""), #within SD2_D1
                paste(rownames(iMat)[7],"\n",colnames(iMat)[7],sep=""), #within SD4_D1
                paste(rownames(iMat)[10],"\n",colnames(iMat)[10],sep=""), #within SD8_D1
                paste(rownames(iMat)[13],"\n",colnames(iMat)[13],sep="") #within SD16_D1
        ),
        main="D1",ylim=limD,
        las=2)


boxplot(iMat[3,6][[1]], #between SD1_D3 and SD2_D3
        iMat[3,9][[1]], #between SD1_D3 and SD4_D3
        iMat[3,12][[1]], #between SD1_D3 and SD8_D3
        iMat[3,15][[1]], #between SD1_D3 and SD16_D3
        names=c(paste(rownames(iMat)[3],"\n",colnames(iMat)[6],sep=""), #between SD1_D3 and SD2_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[9],sep=""), #between SD1_D3 and SD4_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[12],sep=""), #between SD1_D3 and SD8_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[15],sep="") #between SD1_D3 and SD16_D3
        ),
        main="D3: monocultures against others",ylim=limD,
        las=2)

boxplot(iMat[3,6][[1]], #between SD1_D3 and SD2_D3
        iMat[3,9][[1]], #between SD1_D3 and SD4_D3
        iMat[3,12][[1]], #between SD1_D3 and SD8_D3
        iMat[3,15][[1]], #between SD1_D3 and SD16_D3
        iMat[6,9][[1]], #between SD2_D3 and SD4_D3
        iMat[6,12][[1]], #between SD2_D3 and SD8_D3
        iMat[6,15][[1]], #between SD2_D3 and SD16_D3
        iMat[9,12][[1]], #between SD4_D3 and SD8_D3
        iMat[9,15][[1]], #between SD4_D3 and SD16_D3
        iMat[12,15][[1]], #between SD8_D3 and SD16_D3
        names=c(paste(rownames(iMat)[3],"\n",colnames(iMat)[6],sep=""), #between SD1_D3 and SD2_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[9],sep=""), #between SD1_D3 and SD4_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[12],sep=""), #between SD1_D3 and SD8_D3
                paste(rownames(iMat)[3],"\n",colnames(iMat)[15],sep=""), #between SD1_D3 and SD16_D3
                paste(rownames(iMat)[6],"\n",colnames(iMat)[9],sep=""), #between SD2_D3 and SD4_D3
                paste(rownames(iMat)[6],"\n",colnames(iMat)[12],sep=""), #between SD2_D3 and SD8_D3
                paste(rownames(iMat)[6],"\n",colnames(iMat)[15],sep=""), #between SD2_D3 and SD16_D3
                paste(rownames(iMat)[9],"\n",colnames(iMat)[12],sep=""), #between SD4_D3 and SD8_D3
                paste(rownames(iMat)[9],"\n",colnames(iMat)[15],sep=""), #between SD4_D3 and SD16_D3
                paste(rownames(iMat)[12],"\n",colnames(iMat)[15],sep="") #between SD8_D3 and SD16_D3
        ),
        main="D3: all levels against all",ylim=limD,
        las=2)


boxplot(iMat[2,5][[1]], #between SD1_D2 and SD2_D2
        iMat[2,8][[1]], #between SD1_D2 and SD4_D2
        iMat[2,11][[1]], #between SD1_D2 and SD8_D2
        iMat[2,14][[1]], #between SD1_D2 and SD16_D2
        names=c(paste(rownames(iMat)[2],"\n",colnames(iMat)[5],sep=""), #between SD1_D2 and SD2_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[8],sep=""), #between SD1_D2 and SD4_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[11],sep=""), #between SD1_D2 and SD8_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[14],sep="") #between SD1_D2 and SD16_D2
        ),
        main="D2: monocultures against others",ylim=limD,
        las=2)

boxplot(iMat[2,5][[1]], #between SD1_D2 and SD2_D2
        iMat[2,8][[1]], #between SD1_D2 and SD4_D2
        iMat[2,11][[1]], #between SD1_D2 and SD8_D2
        iMat[2,14][[1]], #between SD1_D2 and SD16_D2
        iMat[5,8][[1]], #between SD2_D2 and SD4_D2
        iMat[5,11][[1]], #between SD2_D2 and SD8_D2
        iMat[5,14][[1]], #between SD2_D2 and SD16_D2
        iMat[8,11][[1]], #between SD4_D2 and SD8_D2
        iMat[8,14][[1]], #between SD4_D2 and SD16_D2
        iMat[11,14][[1]], #between SD8_D2 and SD16_D2
        names=c(paste(rownames(iMat)[2],"\n",colnames(iMat)[5],sep=""), #between SD1_D2 and SD2_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[8],sep=""), #between SD1_D2 and SD4_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[11],sep=""), #between SD1_D2 and SD8_D2
                paste(rownames(iMat)[2],"\n",colnames(iMat)[14],sep=""), #between SD1_D2 and SD16_D2
                paste(rownames(iMat)[5],"\n",colnames(iMat)[8],sep=""), #between SD2_D2 and SD4_D2
                paste(rownames(iMat)[5],"\n",colnames(iMat)[11],sep=""), #between SD2_D2 and SD8_D2
                paste(rownames(iMat)[5],"\n",colnames(iMat)[14],sep=""), #between SD2_D2 and SD16_D2
                paste(rownames(iMat)[8],"\n",colnames(iMat)[11],sep=""), #between SD4_D2 and SD8_D2
                paste(rownames(iMat)[8],"\n",colnames(iMat)[14],sep=""), #between SD4_D2 and SD16_D2
                paste(rownames(iMat)[11],"\n",colnames(iMat)[14],sep="") #between SD8_D2 and SD16_D2
        ),
        main="D2: all levels against all",ylim=limD,
        las=2)

boxplot(iMat[1,4][[1]], #between SD1_D1 and SD2_D1
        iMat[1,7][[1]], #between SD1_D1 and SD4_D1
        iMat[1,10][[1]], #between SD1_D1 and SD8_D1
        iMat[1,13][[1]], #between SD1_D1 and SD16_D1
        names=c(paste(rownames(iMat)[1],"\n",colnames(iMat)[4],sep=""), #between SD1_D1 and SD2_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[7],sep=""), #between SD1_D1 and SD4_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[10],sep=""), #between SD1_D1 and SD8_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[13],sep="") #between SD1_D1 and SD16_D1
        ),
        main="D1: monocultures against others",ylim=limD,
        las=2)

boxplot(iMat[1,4][[1]], #between SD1_D1 and SD2_D1
        iMat[1,7][[1]], #between SD1_D1 and SD4_D1
        iMat[1,10][[1]], #between SD1_D1 and SD8_D1
        iMat[1,13][[1]], #between SD1_D1 and SD16_D1
        iMat[4,7][[1]], #between SD2_D1 and SD4_D1
        iMat[4,10][[1]], #between SD2_D1 and SD8_D1
        iMat[4,13][[1]], #between SD2_D1 and SD16_D1
        iMat[7,10][[1]], #between SD4_D1 and SD8_D1
        iMat[7,13][[1]], #between SD4_D1 and SD16_D1
        iMat[10,13][[1]], #between SD8_D1 and SD16_D1
        names=c(paste(rownames(iMat)[1],"\n",colnames(iMat)[4],sep=""), #between SD1_D1 and SD2_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[7],sep=""), #between SD1_D1 and SD4_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[10],sep=""), #between SD1_D1 and SD8_D1
                paste(rownames(iMat)[1],"\n",colnames(iMat)[13],sep=""), #between SD1_D1 and SD16_D1
                paste(rownames(iMat)[4],"\n",colnames(iMat)[7],sep=""), #between SD2_D1 and SD4_D1
                paste(rownames(iMat)[4],"\n",colnames(iMat)[10],sep=""), #between SD2_D1 and SD8_D1
                paste(rownames(iMat)[4],"\n",colnames(iMat)[13],sep=""), #between SD2_D1 and SD16_D1
                paste(rownames(iMat)[7],"\n",colnames(iMat)[10],sep=""), #between SD4_D1 and SD8_D1
                paste(rownames(iMat)[7],"\n",colnames(iMat)[13],sep=""), #between SD4_D1 and SD16_D1
                paste(rownames(iMat)[10],"\n",colnames(iMat)[13],sep="") #between SD8_D1 and SD16_D1
        ),
        main="D1: all levels against all",ylim=limD,
        las=2)
dev.off()

##### 4.4 Aitchison distances: plot within plot comparisons #####
beta_inplot = iTab %>% 
  mutate(sam1 = substr(sample1, 1,5)) %>% mutate(sam2 = substr(sample2, 1, 5)) %>% 
  filter(sam1 == sam2) %>% select(-sam1, -sam2) %>% 
  mutate(comparison = paste0(substr(sample1, 6,7),"_",substr(sample2, 6,7)))
beta_SD1 = beta_inplot %>% filter(str_detect(group1, pattern = "SD1_"))
beta_SD2 = beta_inplot %>% filter(str_detect(group1, pattern = "SD2_"))
beta_SD4 = beta_inplot %>% filter(str_detect(group1, pattern = "SD4_"))
beta_SD8 = beta_inplot %>% filter(str_detect(group1, pattern = "SD8_"))
beta_SD16 = beta_inplot %>% filter(str_detect(group1, pattern = "SD16_"))

beta_inplot_plants <- merge(beta_inplot, pI, by.x ="sample1" , by.y = 'row.names')

beta_SD1p <- ggplot(beta_SD1) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw() + ylim(15,35)
beta_SD2p <- ggplot(beta_SD2) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(15,35)
beta_SD4p <- ggplot(beta_SD4) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(15,35)
beta_SD8p <- ggplot(beta_SD8) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(15,35)
beta_SD16p <- ggplot(beta_SD16) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(15,35)
beta_inplotp <- ggplot(beta_inplot) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw() + ylim(15,35)

pdf("Graphics/20230922_dBEF_4b_betabox_withinPlot.pdf")
ggarrange(beta_SD1p, beta_SD2p, beta_SD4p, beta_SD8p, beta_SD16p, beta_inplotp, ncol = 3, nrow = 2)
dev.off()

iTab_both <- merge(iTab, iTab_sor, by = c("sample1", "sample2"))
plot(iTab_both$distance.y, iTab_both$distance.x, ylab = "Aitchison distance", xlab = "Soerensen")


# plot beta diversity per plant species
beta_plot_plant <- ggplot(beta_inplot_plants) +
  geom_jitter(aes(x = c(,32:90), y = distance, color = comparison)) +
  theme_bw() + facet_grid(SOWNDIV~.)

# tests if the histories are significantly different
wilcox_test(distance ~ comparison, data = beta_SD1, paired = TRUE)
wilcox_test(distance ~ comparison, data = beta_SD2, paired = TRUE)
wilcox_test(distance ~ comparison, data = beta_SD4, paired = TRUE)
wilcox_test(distance ~ comparison, data = beta_SD8, paired = TRUE)
wilcox_test(distance ~ comparison, data = beta_SD16, paired = TRUE)
wilcox_test(distance ~ comparison, data = beta_inplot, paired = TRUE)
                                                                      
##### 4.5 perMANOVA #####

perm = how(nperm=9999)
setBlocks(perm) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)$BLOCK)

adonis2(exJE_dist ~ treatment*LOG_SD, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
#                     Df SumOfSqs      R2       F Pr(>F)    
# treatment          2     7053 0.09071 11.4954  0.0001 ***
# LOG_SD             1     1606 0.02065  5.2345  0.0001 ***
# treatment:LOG_SD   2      987 0.01270  1.6091  0.0066 **  
# Residual         222    68100 0.87593                   
# Total            227    77746 1.00000 

adonis2(exJE_dist ~ HERB + GRASS + LEG, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
#           Df SumOfSqs      R2      F Pr(>F)    
# HERB       1      503 0.00647 1.4890  0.0737 .  
# GRASS      1      868 0.01116 2.5700  0.0004 ***
# LEG        1      723 0.00930 2.1403  0.0028 ** 
# Residual 224    75652 0.97307                  
# Total    227    77746 1.00000  
adonis2(exJE_dist ~ treatment + LOG_SD + FUNCGR + LEG + GRASS + HERB, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
#                     Df SumOfSqs      R2       F Pr(>F)    
# treatment   2     7053 0.09071 11.6173  0.0001 ***
# LOG_SD      1     1606 0.02065  5.2900  0.0001 ***
# FUNCGR      3      924 0.01189  1.0149  0.3712    
# LEG         1      905 0.01164  2.9804  0.0001 ***
# GRASS       1      534 0.00687  1.7595  0.0385 *  
# HERB        1      553 0.00711  1.8221  0.0110 *  
# Residual  218    66171 0.85113                   
# Total     227    77746 1.00000

# Is history effect NH vs SH > SH vs PSH (and NH vs PSH)?
adonis2(exJE_dist ~ treatment*LOG_SD, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)


# Does functional group p/a interact with treatment?
adonis2(exJE_dist ~ treatment*FUNCGR, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
# no interaction, but FUNCGR sig --> is hidden by plant div in complex model
adonis2(exJE_dist ~ treatment*NUMLEG, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
# no interaction
adonis2(exJE_dist ~ treatment*GRASS, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
# no interaction
adonis2(exJE_dist ~ treatment*HERB, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
# no interaction

adonis2(exJE_dist ~ history_soil+history_plant, data=as(sample_data(JE_dBEF_f), "data.frame"),permutations = perm)
#                 Df SumOfSqs      R2       F Pr(>F)    
# history_soil    1     6252 0.08041 19.8971  0.001 ***
# history_plant   1      801 0.01030  2.5496  0.001 ***
# Residual      225    70693 0.90929                   
# Total         227    77746 1.00000   

##### 4.6 perMANOVA per treatment #####
#calculate distances per treatment
exJE_ps_D1 <- exJE_ps[,grep("D1$",colnames(exJE_ps))]
exJE_ps_D2 <- exJE_ps[,grep("D2$",colnames(exJE_ps))]
exJE_ps_D3 <- exJE_ps[,grep("D3$",colnames(exJE_ps))]
exJE_dist_D1 <- aDist(t(exJE_ps_D1)) 
exJE_dist_D2 <- aDist(t(exJE_ps_D2)) 
exJE_dist_D3 <- aDist(t(exJE_ps_D3)) 

permD1 = how(nperm=9999)
permD2 = how(nperm=9999)
permD3 = how(nperm=9999)
setBlocks(permD1) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D1",]$BLOCK)
setBlocks(permD2) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D2",]$BLOCK)
setBlocks(permD3) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D3",]$BLOCK)

adonis2(exJE_dist_D1 ~ LOG_SD, 
         data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D1",],
         permutations = permD1)
#           Df SumOfSqs      R2     F Pr(>F)   
# LOG_SD    1    614.6 0.02602 1.977 0.0022 **
adonis2(exJE_dist_D2 ~ LOG_SD, 
        data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D2",],
        permutations = permD2)
#           Df SumOfSqs      R2      F Pr(>F)   
# LOG_SD    1    544.9 0.02345 1.7772 0.0048 **
adonis2(exJE_dist_D3 ~ LOG_SD, 
        data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D3",],
        permutations = permD3)
#           Df SumOfSqs      R2      F Pr(>F)    
# LOG_SD    1   1433.5 0.06013 4.7341  1e-04 ***

# Does effect of functional groups increase with age? 
adonis2(exJE_dist_D1 ~ LOG_SD + FUNCGR + LEG + GRASS + HERB, 
        data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D1",],
        permutations = permD1)
#           Df SumOfSqs      R2     F Pr(>F)   
# LOG_SD    1    614.6 0.02602 2.0057 0.0021 **
# FUNCGR    3    787.2 0.03333 0.8563 0.7362   
# LEG       1    483.0 0.02045 1.5765 0.0242 * 
# GRASS     1    534.1 0.02261 1.7432 0.0192 * 
# HERB      1    364.2 0.01542 1.1887 0.1882   
# Residual 68  20835.8 0.88216                 
# Total    75  23618.9 1.00000
adonis2(exJE_dist_D2 ~ LOG_SD + FUNCGR + LEG + GRASS + HERB, 
        data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D2",],
        permutations = permD2)
#           Df SumOfSqs      R2      F Pr(>F)   
# LOG_SD    1    544.9 0.02345 1.7992 0.0039 ** 
# FUNCGR    3    886.8 0.03817 0.9761 0.3974    
# LEG       1    694.3 0.02988 2.2926 0.0001 ***
# GRASS     1    284.1 0.01223 0.9381 0.5279    
# HERB      1    229.7 0.00989 0.7585 0.7844    
# Residual 68  20594.3 0.88638                  
# Total    75  23234.2 1.00000
adonis2(exJE_dist_D3 ~ LOG_SD + FUNCGR + LEG + GRASS + HERB, 
        data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment == "D3",],
        permutations = permD3)
#           Df SumOfSqs      R2      F Pr(>F)    
# LOG_SD    1   1433.5 0.06013 4.7368 0.0001 ***
# FUNCGR    3    729.1 0.03058 0.8031 0.8760    
# LEG       1    475.0 0.01992 1.5695 0.0395 *  
# GRASS     1    269.1 0.01129 0.8893 0.6118    
# HERB      1    355.1 0.01490 1.1735 0.2219    
# Residual 68  20578.3 0.86318                  
# Total    75  23840.1 1.00000    

# FG      R^2 = 0.033 > 0.038 > 0.031 (all ns)
# LEG     R^2 = 0.020 > 0.029 > 0.019 (all sig)
# GRASS   R^2 = 0.023 > 0.012 > 0.011 (only D1 sig)
# HERB    R^2 = 0.015 > 0.009 > 0.015 (all ns)

# Is history effect NH vs SH > SH vs PSH (and NH vs PSH)?
exJE_ps_D12 <- exJE_ps[,!grepl("D3$",colnames(exJE_ps))]
exJE_ps_D23 <- exJE_ps[,!grepl("D1$",colnames(exJE_ps))]
exJE_ps_D13 <- exJE_ps[,!grepl("D2$",colnames(exJE_ps))]
exJE_dist_D12 <- aDist(t(exJE_ps_D12)) 
exJE_dist_D23 <- aDist(t(exJE_ps_D23)) 
exJE_dist_D13 <- aDist(t(exJE_ps_D13)) 

permD12 = how(nperm=9999)
permD23 = how(nperm=9999)
permD13 = how(nperm=9999)
# setBlocks(permD12) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D3",]$BLOCK)
# setBlocks(permD23) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D1",]$BLOCK)
# setBlocks(permD13) <- with(sample_data(JE_dBEF_f),sample_data(JE_dBEF_f)[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D2",]$BLOCK)

adonis2(exJE_dist_D12 ~ plot+treatment, data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D3",],
        permutations = permD12)
# history R2 = 0.096 F = 16.966 p < 0.001
adonis2(exJE_dist_D23 ~ plot+treatment, data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D1",],
        permutations = permD23)
# history R2 = 0.017 F = 3.057 p = 1e-04
adonis2(exJE_dist_D13 ~ plot+treatment, data=as(sample_data(JE_dBEF_f), "data.frame")[as(sample_data(JE_dBEF_f), "data.frame")$treatment != "D2",],
        permutations = permD13)
# history R2 = 0.092 F = 16.54 p = 1e-04


##### 4.7 Soerenson distances: extract pairwise comparisons #####
exJE_ps_pa <- as.matrix(exJE_ps_pa)

exJE_sor_dist <- vegdist(t(exJE_ps_pa), method = "bray", binary = F) 
exJE_sor_distM <- as.matrix(exJE_sor_dist)

# check that my metadata makes sense:
all(colnames(exJE_sor_distM)==sI$sample)
all(rownames(exJE_sor_distM)==sI$sample)

# extract groups of interest (from upper triangle):
iMat_sor <- array(list(),dim=c(15,15), # new data will go into an array with lists of values
              dimnames=list(paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                 each=length(unique(sI$treatment))),sep=""),
                                  rep(sort(unique(sI$treatment)),
                                      times=length(unique(sI$SOWNDIV))),sep="_"),
                            paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                 each=length(unique(sI$treatment))),sep=""),
                                  rep(sort(unique(sI$treatment)),
                                      times=length(unique(sI$SOWNDIV))),sep="_")))
# long table version:
iTab_sor <- data.frame("sample1"="","sample2"="","group1"="","group2"="",
                   "common"="",
                   "distance"=rep(0,length(c(exJE_sor_distM))/2-(nrow(exJE_sor_distM)/2)))
cnt <- 1
for(i in 1:(nrow(exJE_sor_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(exJE_sor_distM)){ #j are columns - only read upper triangle
    comp1 <- paste(paste("SD",sI$SOWNDIV[i],sep=""),sI$treatment[i],sep="_")
    comp2 <- paste(paste("SD",sI$SOWNDIV[j],sep=""),sI$treatment[j],sep="_")
    fields <- sort(c(which(rownames(iMat_sor)==comp1),which(rownames(iMat_sor)==comp2))) # only fill upper triangle
    iMat_sor[fields[1],fields[2]] <- list(append(iMat_sor[fields[1],fields[2]][[1]],
                                             exJE_sor_distM[i,j]))
    iTab_sor$distance[cnt] <- exJE_sor_distM[i,j]
    iTab_sor$sample1[cnt] <- rownames(exJE_sor_distM)[i]
    iTab_sor$sample2[cnt] <- colnames(exJE_sor_distM)[j]
    iTab_sor$group1[cnt] <- comp1
    iTab_sor$group2[cnt] <- comp2
    coms <- ""
    if(sI$SOWNDIV[i]==sI$SOWNDIV[j]) coms <- paste(coms,paste0("SD",sI$SOWNDIV[j]),sep="_")
    if(sI$treatment[i]==sI$treatment[j]) coms <- paste(coms,sI$treatment[j],sep="_")
    if(sI$plot[i]==sI$plot[j]) coms <- paste(coms,sI$plot[j],sep="_")
    coms <- sub("^_","",coms)
    iTab_sor$common[cnt] <- coms
    cnt <- cnt + 1
  }
}


iMatP_sor <- array(list(),dim=c(15,15), # new data will go into an array with lists of values
               dimnames=list(paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                  each=length(unique(sI$treatment))),sep=""),
                                   rep(sort(unique(sI$treatment)),
                                       times=length(unique(sI$SOWNDIV))),sep="_"),
                             paste(paste("SD",rep(sort(unique(sI$SOWNDIV)),
                                                  each=length(unique(sI$treatment))),sep=""),
                                   rep(sort(unique(sI$treatment)),
                                       times=length(unique(sI$SOWNDIV))),sep="_")))

for(i in 1:(nrow(exJE_sor_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(exJE_sor_distM)){ #j are columns - only read upper triangle
    # for comparisons of D1-D3, use only samples from the same plot, keep other fields empty:
    if(sI$plot[i]==sI$plot[j]){
      comp1 <- paste(paste("SD",sI$SOWNDIV[i],sep=""),sI$treatment[i],sep="_")
      comp2 <- paste(paste("SD",sI$SOWNDIV[j],sep=""),sI$treatment[j],sep="_")
      fields <- sort(c(which(rownames(iMatP_sor)==comp1),which(rownames(iMatP_sor)==comp2))) # only fill upper triangle
      iMatP_sor[fields[1],fields[2]] <- list(append(iMatP_sor[fields[1],fields[2]][[1]],
                                                exJE_sor_distM[i,j]))
    }
  }
}

##### 4.8 Soerenson distances: test pairwise comparisons #####
# are old monocultures more similar to each other than new ones?
wilcox.test(iMat_sor[3,3][[1]],
            iMat_sor[1,1][[1]],paired = T)
# no  p = 0.28 // same as Aitch? no
wilcox.test(iMat_sor[2,2][[1]],
            iMat_sor[1,1][[1]],paired = T)
# no p = 0.7742 // same as Aitch? yes
wilcox.test(iMat_sor[3,3][[1]],
            iMat_sor[2,2][[1]],paired = T)
# no p = 0.0006135 // same as Aitch? no

# are old high-diverse mixtures more similar to each other than new ones?
wilcox.test(iMat_sor[15,15][[1]],
            iMat_sor[13,13][[1]],paired = T)
#no p = 0.3054 // same as Aitch? no
wilcox.test(iMat_sor[14,14][[1]],
            iMat_sor[13,13][[1]],paired = T)
#yes p = 0.00027 // same as Aitch? no
wilcox.test(iMat_sor[15,15][[1]],
            iMat_sor[14,14][[1]],paired = T)
#yes p = 2.606e-06 // same as Aitch? no

#in the middle:
# 2-species mixes
wilcox.test(iMat_sor[6,6][[1]],
            iMat_sor[4,4][[1]],paired = T)
# yes (means of D3 are higher than of D1): p = 0.01181 // same as Aitch? no
wilcox.test(iMat_sor[5,5][[1]],
            iMat_sor[4,4][[1]],paired = T)
# no p = 0.4116 // same as Aitch? no
wilcox.test(iMat_sor[6,6][[1]],
            iMat_sor[5,5][[1]],paired = T)
# yes p = 0.007784 // same as Aitch? yes

# 4-species mixes
wilcox.test(iMat_sor[9,9][[1]],
            iMat_sor[7,7][[1]],paired = T)
# no p = 0.2625 // same as Aitch? yes
wilcox.test(iMat_sor[8,8][[1]],
            iMat_sor[7,7][[1]],paired = T)
# yes p = 0.01306 // same as Aitch? no
wilcox.test(iMat_sor[9,9][[1]],
            iMat_sor[8,8][[1]],paired = T)
# yes p = 0.0001078 // same as Aitch? no

# 8-species mixes
wilcox.test(iMat_sor[12,12][[1]],
            iMat_sor[10,10][[1]],paired = T)
# yes p = 0.00158 // same as Aitch? no
wilcox.test(iMat_sor[11,11][[1]],
            iMat_sor[10,10][[1]],paired = T)
# yes p = 0.01559 // same as Aitch? no
wilcox.test(iMat_sor[12,12][[1]],
            iMat_sor[11,11][[1]],paired = T)
# no p = 0.9155 // same as Aitch? yes

# are old monoculture plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP_sor[1,3][[1]], #D3 vs D1
            iMatP_sor[2,3][[1]]) #D3 vs D2
# no p = 0.945 // same as Aitch? no

# are no-history monoculture plots more different from those with soil history and old or new plants?
wilcox.test(iMatP_sor[1,2][[1]], #D1 vs D2
            iMatP_sor[1,3][[1]]) #D1 vs D3
# no p = 0.3945 // same as Aitch? yes

# are old high-diversity plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP_sor[13,15][[1]], #D3 vs D1
            iMatP_sor[14,15][[1]]) #D3 vs D2
# no p = 0.3819

# are no-history high-diversity plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[13,14][[1]], #D1 vs D2
            iMatP[13,15][[1]]) #D1 vs D3
# no p = 0.91

# are old 2-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[4,6][[1]], #D3 vs D1
            iMatP[5,6][[1]]) #D3 vs D2
# (yes) p = 0.01356

# are no-history 2-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[4,5][[1]], #D1 vs D2
            iMatP[4,6][[1]]) #D1 vs D3
# no p = 0.7804

# are old 4-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[7,9][[1]], #D3 vs D1
            iMatP[8,9][[1]]) #D3 vs D2
# yes p = 0.00019

# are no-history 4-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[7,8][[1]], #D1 vs D2
            iMatP[7,9][[1]]) #D1 vs D3
# no p = 0.3045

# are old 8-spec plots more different from those with new plants and with or without soil history?
wilcox.test(iMatP[10,12][[1]], #D3 vs D1
            iMatP[11,12][[1]]) #D3 vs D2
# yes p = 0.002593

# are no-history 8-spec plots more different from those with soil history and old or new plants?
wilcox.test(iMatP[10,11][[1]], #D1 vs D2
            iMatP[10,12][[1]]) #D1 vs D3
# no p = 0.8672



# are monocultures more dissimilar from each other than high-diversity plots (old treatment)?
wilcox.test(iMat[3,3][[1]],
            iMat[15,15][[1]])
#no p = 0.704

# are monocultures more dissimilar from each other than high-diversity plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[13,13][[1]])
# (yes) p = 0.03835

# are monocultures more dissimilar from each other than 2-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[4,4][[1]])
# yes p = 0.0002782

# are monocultures more dissimilar from each other than 4-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[7,7][[1]])
# no p = 0.5216

# are monocultures more dissimilar from each other than 8-species plots (new treatment)?
wilcox.test(iMat[1,1][[1]],
            iMat[10,10][[1]])
# yes p = 0.0004263

# are 2-species mixes more dissimilar from each other than 4-species plots (new treatment)?
wilcox.test(iMat[7,7][[1]],
            iMat[4,4][[1]])
# yes (2 are more dissimilar) p = 0.01241

# are monocultures more dissimilar from each other than high-diversity plots (treatment D2)?
wilcox.test(iMat[2,2][[1]],
            iMat[14,14][[1]])
# (yes) (boarderline) p = 0.04886




##### 4.9 Soerenson distances: plot within plot comparisons #####
beta_sor_inplot = iTab_sor %>% 
  mutate(sam1 = substr(sample1, 1,5)) %>% mutate(sam2 = substr(sample2, 1, 5)) %>% 
  filter(sam1 == sam2) %>% select(-sam1, -sam2) %>% 
  mutate(comparison = paste0(substr(sample1, 6,7),"_",substr(sample2, 6,7)))
beta_sor_SD1 = beta_sor_inplot %>% filter(str_detect(group1, pattern = "SD1_"))
beta_sor_SD2 = beta_sor_inplot %>% filter(str_detect(group1, pattern = "SD2_"))
beta_sor_SD4 = beta_sor_inplot %>% filter(str_detect(group1, pattern = "SD4_"))
beta_sor_SD8 = beta_sor_inplot %>% filter(str_detect(group1, pattern = "SD8_"))
beta_sor_SD16 = beta_sor_inplot %>% filter(str_detect(group1, pattern = "SD16_"))


beta_sor_SD1p <- ggplot(beta_sor_SD1) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw() + ylim(0,1)
beta_sor_SD2p <- ggplot(beta_sor_SD2) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(0,1)
beta_sor_SD4p <- ggplot(beta_sor_SD4) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(0,1)
beta_sor_SD8p <- ggplot(beta_sor_SD8) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(0,1)
beta_sor_SD16p <- ggplot(beta_sor_SD16) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw()+ ylim(0,1)
beta_sor_inplotp <- ggplot(beta_sor_inplot) +
  geom_boxplot(aes(comparison, distance)) + geom_jitter(aes(comparison, distance), width = 0.15) +
  theme_bw() + ylim(0,1)

pdf("Graphics/20230922_dBEF_4b_betabox_withinPlot.pdf")
ggarrange(beta_sor_SD1p, beta_sor_SD2p, beta_sor_SD4p, beta_sor_SD8p, beta_sor_SD16p, beta_sor_inplotp, ncol = 3, nrow = 2)
dev.off()

beta_inplot_t <- t_test(beta_inplot, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                            c("D1_D3","D2_D3")))
beta_SD1_t <- t_test(beta_SD1, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                      c("D1_D3","D2_D3")))
beta_SD2_t <- t_test(beta_SD2, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                      c("D1_D3","D2_D3")))
beta_SD4_t <- t_test(beta_SD4, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                      c("D1_D3","D2_D3")))
beta_SD8_t <- t_test(beta_SD8, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                      c("D1_D3","D2_D3")))
beta_SD16_t <- t_test(beta_SD16, distance~comparison, comparison = list(c("D1_D2","D1_D3"), c("D1_D2","D2_D3"),
                                                                        c("D1_D3","D2_D3")))
##### 4.10 disect Beta Div with Plants ##### 
#### plant data ####
pInfo <- read.delim("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/plotInfo.txt", row.names = 1)
plMat <- pInfo[,grep("[[:upper:]]{3}_[[:upper:]]{3}",colnames(pInfo))]
plMat$CRE_BIE <- pInfo$SOWNDIV-rowSums(pInfo[,grepl("^[[:upper:]]{3}_[[:upper:]]{3}$",
                                                    colnames(pInfo))])
#rownames(plMat) <- pInfo$plot
plMat <- plMat[pInfo$SOWNDIV<60,]

library(vegan)
pJE_dist <- vegdist(plMat,method = "bray",binary = T) # 
pJE_distM <- as.matrix(pJE_dist)

# check that my metadata makes sense:
all(colnames(pJE_distM)==sI$plot[sI$treatment=="D3"])

# extract groups of interest (from upper triangle):

# long table version:
iTab$distancePlants <- 0

for(i in 1:(nrow(pJE_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(pJE_distM)){ #j are columns - only read upper triangle
    s1pattern <- paste0("^",rownames(pJE_distM)[i])
    s2pattern <- paste0("^",colnames(pJE_distM)[j])
    iTab$distancePlants[grepl(s1pattern,iTab$sample1) &
                          grepl(s2pattern,iTab$sample2)] <- pJE_distM[i,j]
  }
}

#subset
beta_plantplot <- unique(iTab[grep("SD.+_D.",iTab$common),])
beta_plantplot$treatment <- gsub(".*_","",beta_plantplot$common)
beta_plantplot$SD <- gsub("_.*","",beta_plantplot$common)

bp <- ggplot(beta_plantplot, 
             mapping =aes(x=distancePlants,y=distance,color=treatment)) +
  geom_point() + 
  geom_smooth(method="lm")
bp + facet_wrap(~SD, ncol=3)

bp1a <- ggplot(beta_plantplot, 
               mapping =aes(x=distancePlants,y=distance)) +
  geom_smooth(method="lm",color="black") +
  geom_point(aes(color=treatment)) + 
  facet_wrap(~SD, ncol=3)

bp2 <- ggplot(beta_plantplot, 
              mapping =aes(x=distancePlants,y=distance,color=SD)) +
  geom_point() + 
  geom_smooth(method="lm")
bp2 + facet_wrap(~treatment, ncol=3)

bp2a <- ggplot(beta_plantplot, 
               mapping =aes(x=distancePlants,y=distance)) +
  geom_smooth(method="lm",color="black") +
  geom_point(aes(color=SD)) + 
  facet_wrap(~treatment, ncol=3)


#### real plant data ####
bm <- read.delim("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/JE_dBEF_PlantBM_2021.txt",
                 dec=",")
bm$sample <- apply(bm[,c("plotcode","treatment")],1,function(x) paste0(x,collapse=""))
bmm <- bm[bm$season=="May",]
matBM <- tapply(bmm$biomass,list(bmm$sample,bmm$species_other),mean)
matBM[is.na(matBM)]<-0
matBMD <- matBM[rownames(matBM) %in% rownames(exJE),
                toupper(gsub(".","_",colnames(matBM),fixed=T)) %in% colnames(plMat)]

p2JE_dist <- vegdist(matBMD,method = "bray",binary = F) # 
p2JE_distM <- as.matrix(p2JE_dist)


# check that my metadata makes sense:
all(colnames(p2JE_distM)==sI$sample)
all(rownames(p2JE_distM)==sI$sample)

# extract groups of interest (from upper triangle):

# long table version:
iTab$distanceRealPlants <- 0
for(i in 1:(nrow(p2JE_distM)-1)){ #i are rows 
  for(j in (i+1):ncol(p2JE_distM)){ #j are columns - only read upper triangle
    iTab$distanceRealPlants[iTab$sample1==rownames(p2JE_distM)[i] &
                              iTab$sample2==colnames(p2JE_distM)[j]] <- p2JE_distM[i,j]
  }
}

#subset
beta_plantplot2 <- iTab[grep("SD.+_D.",iTab$common),]
beta_plantplot2$treatment <- gsub(".*_","",beta_plantplot2$common)
beta_plantplot2$SD <- gsub("_.*","",beta_plantplot2$common)

bpr <- ggplot(beta_plantplot2, 
              mapping =aes(x=distanceRealPlants,y=distance,color=treatment)) +
  geom_point() + 
  geom_smooth(method="lm")
bpr + facet_wrap(~SD, ncol=3)

bpr1a <- ggplot(beta_plantplot2, 
                mapping =aes(x=distanceRealPlants,y=distance)) +
  geom_smooth(method="lm",color="black") +
  geom_point(aes(color=treatment)) + 
  facet_wrap(~SD, ncol=3)

bpr2 <- ggplot(beta_plantplot2, 
               mapping =aes(x=distanceRealPlants,y=distance,color=SD)) +
  geom_point() + 
  geom_smooth(method="lm")
bpr2 + facet_wrap(~treatment, ncol=3)

bpr2a <- ggplot(beta_plantplot2, 
                mapping =aes(x=distanceRealPlants,y=distance)) +
  geom_smooth(method="lm",color="black") +
  geom_point(aes(color=SD)) + 
  facet_wrap(~treatment, ncol=3)


#### 5 relative and differential abundances ####
##### 5.1 relative abundances #####
dBEF_relabu <- JE_dBEF_f %>% 
  tax_glom(taxrank="Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt() %>% 
  #filter(Abundance>0.01) %>% 
  arrange(Genus) 

lm_relabuSD <- lme(Abundance ~ LOG_SD, random = ~ 1|OTU, data=dBEF_relabu)
summary(lm_relabuSD) # ns
lm_relabuT <- lme(Abundance ~ treatment, random = ~ 1|OTU, data=dBEF_relabu)
summary(lm_relabuT) # ns
lm_relabuSDT <- lme(Abundance ~ treatment*LOG_SD, random = ~ 1|OTU, data=dBEF_relabu)
summary(lm_relabuSDT) # ns

unique(dBEF_relabu$Genus) # 14 genera
table(dBEF_relabu$Genus)

# plot relative abundances
dBEF_relabu_plot <- ggplot(dBEF_relabu, aes(x=as.factor(SOWNDIV), y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill", show.legend=TRUE) +
  #scale_fill_manual(values=mycols, name="Genus") +
  scale_fill_paletteer_d("ggthemes::Classic_20", name ="Genus") +
  #scale_fill_paletteer_d("colorBlindness::paletteMartin", name = "Genus") +
  ylab("Relativ Abundance") + xlab("treatment") +
  theme(legend.position="right", legend.text=element_text(size=6), legend.title=element_text(size=10), legend.key.size=unit(0.5, 'cm')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("relative abundance of AMF genera") 

pdf("Graphics/20230504_dBEF_5_relabu.pdf")
dBEF_relabu_plot + coord_flip() + facet_grid(treatment~.) + theme_pubclean()
dev.off()

##### 5.2 Maaslin #####
# Does treatment or plant diversity create significant differential abundances?
library(Maaslin2)
# Maaslin uses 2 input files: asv/otu table and meta data with features as columns and 
# samples as rows for both  
maas_asv <- t(as.data.frame(otu_table(JE_dBEF_f)))
maas_meta <- data.frame(dBEF_samples)
#maas_meta <- tibble::column_to_rownames(maas_meta, var = "Row.names")

maas_meta = maas_meta %>% mutate(treat_SD = paste(treatment, SOWNDIV, sep = "_"))

# you can specifiy different GLMs/normalizations/transforms. We used LM on CLR-transformed data
# and ZINB model on TSS normalized data
fit_data_treat <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_output_T",
  analysis_method = "LM",
  #transform = "AST",
  fixed_effects = c("treatment"), #,"SOWNDIV", "treat_SD"), #last for interaction of treatment*SOWNDIV
  reference = c("treatment,D1"),#,"treat_SD,D1_1"),  
  random_effects = "plot",
  normalization = "CLR",
  standardize = FALSE,
  min_prevalence = 0 # prev filtering already done
)
fit_data_SD <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_output_SD",
  analysis_method = "LM",
  #transform = "AST",
  fixed_effects = c("SOWNDIV"), #,"SOWNDIV", "treat_SD"), #last two for interaction of treatment*SOWNDIV
  reference = c("SOWNDIV,1"),#,"treat_SD,D1_1"),  
  random_effects = c("BLOCK","treatment"),
  normalization = "CLR",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
fit_data_inter <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_output_TSD",
  analysis_method = "LM",
  #transform = "AST",
  fixed_effects = c("treatment", "SOWNDIV", "treat_SD"), #last two for interaction of treatment*SOWNDIV
  reference = c("treatment,D1","treat_SD,D1_1"),  
  random_effects = "BLOCK",
  normalization = "CLR",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

fit_data_treat_ZINB <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_ZINB_T",
  analysis_method = "ZINB",
  transform = "NONE",
  fixed_effects = c("treatment"), #,"SOWNDIV", "treat_SD"), #last for interaction of treatment*SOWNDIV
  reference = c("treatment,D1"),#,"treat_SD,D1_1"),  
  random_effects = "plot",
  normalization = "TMM",
  standardize = FALSE,
  min_prevalence = 0 # prev filtering already done
)
fit_data_SD_ZINB <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_ZINB_SD_v2",
  analysis_method = "ZINB",
  transform = "NONE",
  fixed_effects = c("LOG_SD"), #,"SOWNDIV", "treat_SD"), #last two for interaction of treatment*SOWNDIV
  reference = c("LOG_SD,1"),#,"treat_SD,D1_1"),  
  random_effects = c("BLOCK","treatment"),
  normalization = "TMM",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)

fit_data_inter_ZINB <- Maaslin2(
  maas_asv,
  maas_meta,
  output = "Maaslin_ZINB_TSD_v2",
  analysis_method = "ZINB", 
  transform = "NONE",
  fixed_effects = c("treatment","SOWNDIV","treat_SD"), 
  reference = c("treatment,D1","treat_SD,D1_1"),  
  random_effects = c("BLOCK","plot"),
  normalization = "TMM", #or none
  standardize = FALSE,
  min_prevalence = 0 # prev filtering already done
)

##### 5.3 relabu & diffabu plot #####
maas_VT <- data.frame(VT=as.character(fit_data_treat_ZINB$results$feature),
                 coef=as.numeric(fit_data_treat_ZINB$results$coef))
maas_VT = maas_VT %>% distinct(VT) 
maas_VT = maas_VT[1:39,]

maas_VT_SD <- data.frame(VT=as.character(fit_data_SD_ZINB$results$feature),
                      coef=as.numeric(fit_data_SD_ZINB$results$coef))
maas_VT_SD = maas_VT_SD %>% distinct(VT) 
maas_VT_SD = maas_VT_SD[1:28,]

#tax proportions
dBEF_relabu_prop = taxa_proportions(JE_dBEF_f, "Genus", treatment = c("treatment", "SOWNDIV"))
dBEF_relabu_propT = taxa_proportions(JE_dBEF_f, "Genus", treatment = "treatment")
dBEF_relabu_VTX = taxa_proportions(JE_dBEF_f, 'VTX', treatment = "treatment")
unique_taxa(JE_dBEF_f, treatment = c("treatment", "SOWNDIV")) # none

dBEF_relabu_VTX = dBEF_relabu_VTX %>% 
  filter(treatment == "D1") %>% 
  left_join(VT_dBEF_tax, by = c("VTX" = "VT"))
dBEF_relabu_VTX <- dBEF_relabu_VTX %>% mutate(VTXgen = paste0(VTX," ", Genus))


# only the VTX that show up in Maaslin
dBEF_relabu_VTX_plot <- dBEF_relabu_VTX%>% 
  filter(VTX %in% maas_VT) %>% 
  mutate(VTX = factor(VTX, levels = rev(maas_VT))) %>% 
  mutate(Proportion = log10(Proportion+0.00001)+5.1) %>% 
  ggplot(aes(y = VTX, x = Proportion, fill = Genus)) + geom_col() + #xlim(0,0.25) +
  scale_fill_paletteer_d("ggthemes::Classic_20", name ="Genus") + theme_classic() + 
  theme(legend.position="bottom") + ylab("") + xlab("Rel. Abundance") 


dBEF_relabu_VTX_SD_plot <- dBEF_relabu_VTX %>% 
  filter(VTX %in% maas_VT_SD) %>% 
  mutate(VTX = factor(VTX, levels = rev(maas_VT_SD))) %>% 
  mutate(Proportion = log10(Proportion+0.00001)+5.1) %>% 
  ggplot(aes(y = VTX, x = Proportion, fill = Genus)) + geom_col() + #xlim(0,0.15) +
  scale_fill_manual(values = c("#C7C7C7","#1F77B4","#AEC7E8","#7F7F7F","#2CA02C","#98DF8A","#BCBD22","#FF9896","#9467BD"), name ="Genus") + 
  theme_classic() + 
  theme(legend.position="bottom") + ylab("") + xlab("Rel. Abundance") 

pdf("Graphics/20230706_dBEF_5_relabu.pdf", width = 3, height = 5)
dBEF_relabu_VTX_plot
dBEF_relabu_VTX_SD_plot 
dev.off()

#### 6 species turnover ####
# Are there VT that appear / dosappear over time? 
#otu + meta into long format
dBEF_VT_frf_long = t(as.data.frame(otu_table(JE_dBEF_frf))) 
dBEF_VT_frf_long = as.data.frame(dBEF_VT_frf_long) %>% 
  tibble::rownames_to_column(var="plot") %>% 
  pivot_longer(!plot, names_to = "VTX", values_to = "relabu")
dBEF_VT_frf_long = dBEF_VT_frf_long %>% 
  mutate(treat = ifelse(str_ends(plot, "D1"),"1", 
                        ifelse(str_ends(plot, "D2"),"2","3"))) %>% 
  mutate(treat = as.numeric(treat)) 
dBEF_VT_frf_long$plot = substr(dBEF_VT_frf_long$plot,1,5) #reduce 'plot' to just plot code w/o treatment
dBEF_VT_frf_long$relabu = log(dBEF_VT_frf_long$relabu)

library(codyn)
dBEF_spec_turnover <- turnover(dBEF_VT_frf_long, 
                              time.var = "treat",
                              species.var = "VTX",
                              abundance.var = "relabu",
                              replicate.var = "plot",
                              metric = "total")
dBEF_spec_disappearance <- turnover(dBEF_VT_frf_long, 
                                   time.var = "treat",
                                   species.var = "VTX",
                                   abundance.var = "relabu",
                                   replicate.var = "plot",
                                   metric = "disappearance")
dBEF_spec_appearance <- turnover(dBEF_VT_frf_long, 
                                time.var = "treat",
                                species.var = "VTX",
                                abundance.var = "relabu",
                                replicate.var = "plot",
                                metric = "appearance")
#Format a compiled data frame
dBEF_spec_turnover$metric<-"total"
names(dBEF_spec_turnover)[1]="turnover"

dBEF_spec_appearance$metric<-"appearance"
names(dBEF_spec_appearance)[1]="turnover"

dBEF_spec_disappearance$metric<-"disappearance"
names(dBEF_spec_disappearance)[1]="turnover"

dBEF_spec_allturnover<-bind_rows(dBEF_spec_turnover, dBEF_spec_appearance, dBEF_spec_disappearance) 
#add meta data
dBEF_spec_allturnover<-dBEF_spec_allturnover %>% merge(dBEF_samples[dBEF_samples$treatment == "D1",], by = "plot") %>% 
  select(plot, turnover, treat, metric, SOWNDIV, BLOCK, FUNCGR, GRASS, LEG, HERB)

# What's the mean appearance D1 to D2?
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'appearance' | dBEF_spec_allturnover$treat == 2])
# 0.4076
# and D2 to D3 ? 
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'appearance' | dBEF_spec_allturnover$treat == 3])
# 0.3543
# What's the mean disappearance D1 to D2?
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'disappearance' | dBEF_spec_allturnover$treat == 2])
# 0.3925
# and D2 to D3 ? 
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'disappearance' | dBEF_spec_allturnover$treat == 3])
# 0.3651
# What's the mean total turnover D1 to D2?
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'total' | dBEF_spec_allturnover$treat == 2])
# 0.4694
# and D2 to D3 ? 
mean(dBEF_spec_allturnover$turnover[dBEF_spec_allturnover$metric == 'total' | dBEF_spec_allturnover$treat == 3])
# 0.4425

#Create the graph
pdf("Graphics/20230617_dBEF_6_turnover.pdf")
#turnover_graph <- 
ggplot(dBEF_spec_allturnover, aes(x=treat, y=turnover, color=SOWNDIV)) + #, group = plot)) + 
  #geom_line() + geom_point() + 
  geom_boxplot(aes(x=treat, y=turnover, color=SOWNDIV)) + geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_grid(rows = vars(metric), cols = vars(SOWNDIV)) + 
  theme_classic() + 
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
dev.off()

ggplot(dBEF_spec_allturnover, aes(x=treat, y=turnover, color=metric)) + 
  geom_line(size = 1) + 
  facet_wrap(~SOWNDIV) + 
  theme_bw() + 
  theme(legend.position="bottom")
ggplot(dBEF_spec_allturnover, aes(x=treat, y=turnover, color=metric)) + 
  geom_line(size = 1) + 
  facet_wrap(~metric) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

# testing significant differences
summary(aov(turnover ~ SOWNDIV*treat,data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total",]))
# SD: F = 3.102 p = 0.018 | treat: F = 17.849 p < 0.001 | SD*T F = 0.137 p = 0.968
summary(aov(turnover ~ SOWNDIV*treat,data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "appearance",]))
# SD: F = 0.820 p = 0.515 | treat: F = 0.003 p = 0.975 | SD*T F = 0.877 p = 0.480
summary(aov(turnover ~ SOWNDIV*treat,data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance",]))
# SD: F = 0.744 p = 0.563 | treat: F = 13.804 p < 0.001 | SD*T F = 1.037 p = 0.390

t.test(turnover ~ treat, data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total",], paired = TRUE)
# t = 5.1815, df = 75, p-value = 1.797e-06
t.test(turnover ~ treat, data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "appearance",], paired = TRUE)
#t = 0.046436, df = 75, p-value = 0.9631
t.test(turnover ~ treat,data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance",], paired = TRUE)
# t = 3.4272, df = 75, p-value = 0.0009931

# Is turnover stronger in low-diverse grasslands? 
anova(lme(turnover ~ SOWNDIV + treat, random = ~1|BLOCK, data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total",]))

t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total" |dBEF_spec_allturnover$SOWNDIV==1 ,], paired = TRUE)
# SD1: t = 3.7722, df = 103, p-value = 0.0002704
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total" |dBEF_spec_allturnover$SOWNDIV==2 ,], paired = TRUE)
# SD2: t = 3.9515, df = 107, p-value = 0.0001396
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total" |dBEF_spec_allturnover$SOWNDIV==4 ,], paired = TRUE)
# SD4: t = 3.4175, df = 107, p-value = 0.0008947
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total" |dBEF_spec_allturnover$SOWNDIV==8 ,], paired = TRUE)
# SD8: t = 3.8104, df = 107, p-value = 0.000232
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "total" |dBEF_spec_allturnover$SOWNDIV==16 ,], paired = TRUE)
# SD16: t = 4.2522, df = 103, p-value = 4.667e-05

t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance" |dBEF_spec_allturnover$SOWNDIV==1 ,], paired = TRUE)
# SD1: t = 3.7667, df = 103, p-value = 0.0002757
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance" |dBEF_spec_allturnover$SOWNDIV==2 ,], paired = TRUE)
# SD2: t = 3.207, df = 107, p-value = 0.001769
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance" |dBEF_spec_allturnover$SOWNDIV==4 ,], paired = TRUE)
# SD4: t = 2.8628, df = 107, p-value = 0.005053
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance" |dBEF_spec_allturnover$SOWNDIV==8 ,], paired = TRUE)
# SD8: t = 2.9197, df = 107, p-value = 0.004272
t.test(turnover ~ treat, 
       data = dBEF_spec_allturnover[dBEF_spec_allturnover$metric == "disappearance" |dBEF_spec_allturnover$SOWNDIV==16 ,], paired = TRUE)
# SD16: t = 3.4784, df = 103, p-value = 0.0007404

#rank shift 
# = relative changes in species rank abundances -> degree of species reordering between two time points
dBEF_spec_rankshift <- rank_shift(dBEF_VT_frf_long, 
                                 time.var = "treat",
                                 species.var = "VTX",
                                 abundance.var = "relabu",
                                 replicate.var = "plot"
                                 )
dBEF_spec_rankshift$treat <- as.character(substr(dBEF_spec_rankshift$year_pair, 3,3))
dBEF_spec_rankshift <- dBEF_spec_rankshift %>% left_join(dBEF_samples[,c(1,5:7)], by = "plot") %>% distinct(plot, treat, .keep_all = T)
ggplot(dBEF_spec_rankshift, aes(treat, MRS)) + 
  geom_line(size= 2) + theme_bw() + facet_grid(.~SOWNDIV)

dBEF_spec_ratechange <- rate_change(dBEF_VT_frf_long, 
                                  time.var = "treat",
                                  species.var = "VTX",
                                  abundance.var = "relabu",
                                  replicate.var = "plot")

ggplot(dBEF_spec_ratechange, aes(interval, distance, color = replicate)) + facet_wrap(~replicate) + 
  geom_point() + theme_bw() + stat_smooth(method = "lm", se = F, size = 2)



#### 7 generalists / specialist (Phi-Factor) ####
# Do some VT only occur in presence of specific plant species?
phi.coeff <- function(x,y){
  # x is the vector of absolute frequencies of the species in a subset, for which 
  #    the degree of specialisation is calculated
  # y is the vector of absolute frequencies of species in total
  sx <- sum(x)
  sy <- sum(y)
  if (sx!=0 & sy!=0){
    phi <- unlist(sapply(1:length(x),function(i){
      a <- as.double(x[i])
      b <- as.double(sx-x[i])
      c <- as.double(y[i]-x[i])
      d <- as.double(sy-y[i]-(sx-x[i]))
      (a*d-b*c)/sqrt((a+c)*(b+d)*(a+b)*(c+d))
    }))
    return(phi)
  } else {
    return(rep(0,length(x))) 
  }
}

phiOnPlotMatrix <- function(OTU_table_wide,sampleInfo_plantComm){ #first argument is the OTU table in wide format (OTU names are colnames), 
  #second argument is a table with the plant community information - samples should be in the same order as the first, every column should represent one plant
  if(nrow(sampleInfo_plantComm)!=nrow(OTU_table_wide)){
    print("warning: the dimensions of the input table don't match")
    return(NA)
  }
  if(any(!OTU_table_wide %in% c(0,1))){
    require(vegan)
    otuTab2 <- decostand(OTU_table_wide,"pa",1)
  }else{
    otuTab2 <- as.matrix(OTU_table_wide)
  }
  if(any(!sampleInfo_plantComm %in% c(0,1))){
    require(vegan)
    sampleInfo_plantComm2 <- decostand(sampleInfo_plantComm,"pa",1)
  }else{
    sampleInfo_plantComm2 <- as.matrix(sampleInfo_plantComm)
  }
  result <- matrix(0,nrow=ncol(sampleInfo_plantComm2),
                   ncol=ncol(otuTab2),
                   dimnames=list(colnames(sampleInfo_plantComm2),
                                 colnames(otuTab2)))
  for(j in 1:ncol(sampleInfo_plantComm2)){
    sampleInfo_plant <- sampleInfo_plantComm2[,j]
    otu_x1 <- aggregate(otuTab2,list(sampleInfo_plant),sum)
    rownames(otu_x1) <- otu_x1[,1]
    otu_x1 <- otu_x1[,-1]
    otu_y1 <- table(sampleInfo_plant)
    result[j,] <- apply(otu_x1,2,function(x) phi.coeff(x,otu_y1)[2])
  }
  return(result)
}

#example
pI <- read.delim("Plot_Information_DeltaBEF_2020.txt",row.names = 1)
pIp <- pI[,25:ncol(pI)] 
# I seem to have only 59 plant species, although I thought there are 60?
# exJE <- read.delim("JEdBEF_fVT_v1.txt",row.names = 1) --> Phi on unfiltered data
exJE <- as.data.frame(t(otu_table(JE_dBEF_f)))

hist(rowSums(exJE),breaks=25)
length(which(rowSums(exJE)<6000))
exJE <- exJE[rowSums(exJE)>=6000,] # why ? 

# rarefaction of all samples with more than 6000 reads, remove the other 5
exJE <- rrarefy(exJE, min(rowSums(exJE)))

pIp <- pIp[sapply(rownames(exJE),function(x) which(rownames(pIp)==x)),]
all(rownames(pIp)==rownames(exJE))

# subset treatments
pIp1 <- pIp[grep("D1$",rownames(pIp)),]
pIp2 <- pIp[grep("D2$",rownames(pIp)),]
pIp3 <- pIp[grep("D3$",rownames(pIp)),]

exJE1 <- exJE[grep("D1$",rownames(exJE)),]
exJE2 <- exJE[grep("D2$",rownames(exJE)),]
exJE3 <- exJE[grep("D3$",rownames(exJE)),]

# phi coefficient per treatment
phi1 <- phiOnPlotMatrix(exJE1,pIp1)
phi2 <- phiOnPlotMatrix(exJE2,pIp2)
phi3 <- phiOnPlotMatrix(exJE3,pIp3)

# change community matrix to presence absence in order to be able to perform swapping
exJE1_pa <- decostand(exJE1,"pa",1)
exJE2_pa <- decostand(exJE2,"pa",1)
exJE3_pa <- decostand(exJE3,"pa",1)

# plot phi-matrices
# pdf("Graphics/20230510_dBEF_7_phi1heat.pdf")
# heatmap(phi1, Rowv = NA, col = cm.colors(256))
# heatmap(phi2, Rowv = NA, revC = TRUE, col = cm.colors(256))
# heatmap(phi3, Rowv = NA,revC = TRUE, col = cm.colors(256))
# dev.off()
pdf("Graphics/20230913_dBEF_7_phihist.pdf")
phi2_hist <- hist(phi2, col = rgb(145, 207, 96, alpha = 80, maxColorValue = 255), 
                  xlim = c(-0.5, 0.5), ylim = c(0,650)) # "#91cf60"
phi1_hist <- hist(phi1, col = rgb(1, 102, 94, alpha = 80, maxColorValue = 255),
                  xlim = c(-0.5, 0.5), ylim = c(0,650), add = TRUE)  #"#01665e"
phi3_hist <- hist(phi3, col = rgb(254, 196, 79, alpha = 80, maxColorValue = 255),
                  xlim = c(-0.5, 0.5), ylim = c(0,650), add = TRUE ) #"#fec44f"

plot(density(phi1), col = rgb(1, 102, 94, maxColorValue = 255), lwd = 3)
lines(density(phi2), col = rgb(145, 207, 96, maxColorValue = 255),lwd = 3) # "#91cf60"
lines(density(phi3), col = rgb(254, 196, 79, maxColorValue = 255),lwd = 3)
legend("topright", c("NH","SH", "PSH"),lty = c(1,1,1), lwd = c(2,2,2),col = c("#01665e","#91cf60","#fec44f"),
      bty = "n", y.intersp = 0.6)
dev.off()
# plot 2 using facets to separate as an alternative
ggplot(b, aes(x = population, y = value)) + geom_histogram(stat = "identity") + facet_grid(. ~ lambda)

phi1_long = as.data.frame(phi1) %>% tibble::rownames_to_column(var = "plantspec") %>% 
  pivot_longer(!plantspec, names_to = "VTX", values_to = "phi") %>% 
  left_join(VT_dBEF_tax, by = c("VTX" = "VT"))
phi2_long = as.data.frame(phi2) %>% tibble::rownames_to_column(var = "plantspec") %>% 
  pivot_longer(!plantspec, names_to = "VTX", values_to = "phi") %>% 
  left_join(VT_dBEF_tax, by = c("VTX" = "VT"))
phi3_long = as.data.frame(phi3) %>% tibble::rownames_to_column(var = "plantspec") %>% 
  pivot_longer(!plantspec, names_to = "VTX", values_to = "phi") %>% 
  left_join(VT_dBEF_tax, by = c("VTX" = "VT"))

ggplot(phi1_long, aes(x = Genus, y = phi)) + geom_jitter(aes(color=plantspec)) + geom_boxplot()
ggplot(phi2_long, aes(x = Genus, y = phi)) + geom_jitter(aes(color=plantspec)) + geom_boxplot()
ggplot(phi3_long, aes(x = Genus, y = phi)) + geom_jitter(aes(color=plantspec)) + geom_boxplot()



##### 7.1 nullmodel 1 (c0) #####
# number of nullmodels
nnull <- 1000 # 100 is okay, 1000 would be better

# nullmodel 1 ("c0") keeps frequency of VTs intact, but richness may vary
exJE1_c0 <- simulate(nullmodel(exJE1_pa,"c0"),nnull,seed=101) #nullmodel generation
phi1_c0 <- apply(exJE1_c0,3, function(x) phiOnPlotMatrix(x,pIp1)) #calculate phis for all nullmodel instances
p1_c0 <- matrix(c(1-sapply(1:nrow(phi1_c0),
                           #count how often phi is larger in real data than in nullmodel -> p-value
                           function(i) length(which(c(phi1)[i]-phi1_c0[i,]>0)))/nnull),
                nrow=nrow(phi1),ncol=ncol(phi1)) 
phi1[p.adjust(p1_c0,"fdr")<0.5] #phi values with significant p-value
plot(c(phi1),col=c("black","red")[1+as.numeric(p.adjust(p1_c0,"fdr")<0.5)]) #plot all phis highlighting the significant ones
phi1_c0Sig <- matrix(p.adjust(p1_c0,"fdr")<0.5,nrow=nrow(phi1),ncol=ncol(phi1),
                     dimnames=dimnames(phi1)) #matrix with T/F for each combi
# VTC00020 w/ Bro hor   || VTX00156 w/ Gal mol 
# VTX00130 w/ Aju rep   || VTX00214 w/ Alo pra
# VTX00340 w/ Phl pra, Pla med, Tar off


# same for treatment 2:
exJE2_c0 <- simulate(nullmodel(exJE2_pa,"c0"),nnull,seed=101)
phi2_c0 <- apply(exJE2_c0,3, function(x) phiOnPlotMatrix(x,pIp2))
p2_c0 <- matrix(c(1-sapply(1:nrow(phi2_c0),
                           function(i) length(which(c(phi2)[i]-phi2_c0[i,]>0)))/nnull),
                nrow=nrow(phi2),ncol=ncol(phi2))
phi2[p.adjust(p2_c0,"fdr")<0.5]
plot(c(phi2),col=c("black","red")[1+as.numeric(p.adjust(p2_c0,"fdr")<0.5)])
phi2_c0Sig <- matrix(p.adjust(p2_c0,"fdr")<0.5,nrow=nrow(phi2),ncol=ncol(phi2),
                     dimnames=dimnames(phi2))
# VTC00004 w/ Tri cam  || VTX00166 w/ Poa pra
# VTX00245 w/ Poa pra  || VTX00281 w/ poa pra
# VTX00335 w/ Poa pra, Ant odo

# same for treatment 3:
exJE3_c0 <- simulate(nullmodel(exJE3_pa,"c0"),nnull,seed=101)
phi3_c0 <- apply(exJE3_c0,3, function(x) phiOnPlotMatrix(x,pIp3))
p3_c0 <- matrix(c(1-sapply(1:nrow(phi3_c0),
                           function(i) length(which(c(phi3)[i]-phi3_c0[i,]>0)))/nnull),
                nrow=nrow(phi3),ncol=ncol(phi3))
phi3[p.adjust(p3_c0,"fdr")<0.5] # no sig
plot(c(phi3),col=c("black","red")[1+as.numeric(p.adjust(p3_c0,"fdr")<0.5)])
phi3_c0Sig <- matrix(p.adjust(p3_c0,"fdr")<0.5,nrow=nrow(phi3),ncol=ncol(phi3),
                     dimnames=dimnames(phi3))

# there are less specializations in D2 and D3 than D1
which(phi1_c0Sig & phi2_c0Sig) # 0
which(phi1_c0Sig & phi3_c0Sig) # 0
which(phi2_c0Sig & phi3_c0Sig) # 0



##### 7.2 nullmodel 2 (greedyswap) #####
# nullmodel 2 ("greedyqswap") keeps frequency and richness intact - everything else as above
exJE1_gs <- simulate(nullmodel(exJE1_pa,"greedyqswap"),nnull,101,thin=100)
phi1_gs <- apply(exJE1_gs,3, function(x) phiOnPlotMatrix(x,pIp1))
p1_gs <- matrix(c(1-sapply(1:nrow(phi1_gs),
                           function(i) length(which(c(phi1)[i]-phi1_gs[i,]>0)))/nnull),
                nrow=nrow(phi1),
                ncol=ncol(phi1))
phi1[p.adjust(p1_gs,"fdr")<0.05] # just 1
plot(c(phi1),col=c("black","red")[1+as.numeric(p.adjust(p1_gs,"fdr")<0.05)])
phi1_gsSig <- matrix(p.adjust(p1_gs,"fdr")<0.05,nrow=nrow(phi1),ncol=ncol(phi1),
                     dimnames=dimnames(phi1))
# VTX00156 w/ Gal mol

# agreement between the two nullmodels:
which(phi1_c0Sig & phi1_gsSig) # VTX00156 w/ Gal mol
plot(p1_c0,p1_gs,cex=0.3)

# nullmodel 2 on treatment 2
exJE2_gs <- simulate(nullmodel(exJE2_pa,"greedyqswap"),nnull,101,thin=100)
phi2_gs <- apply(exJE2_gs,3, function(x) phiOnPlotMatrix(x,pIp2))
p2_gs <- matrix(c(1-sapply(1:nrow(phi2_gs),
                           function(i) length(which(c(phi2)[i]-phi2_gs[i,]>0)))/nnull),
                nrow=nrow(phi2),
                ncol=ncol(phi2))
phi2[p.adjust(p2_gs,"fdr")<0.05] # 2 left
plot(c(phi2),col=c("black","red")[1+as.numeric(p.adjust(p2_gs,"fdr")<0.05)])
phi2_gsSig <- matrix(p.adjust(p2_gs,"fdr")<0.05,nrow=nrow(phi2),ncol=ncol(phi2),
                     dimnames=dimnames(phi2))
# VTX 00245 w/ Poa pra || VTX00153 w/ Pim maj

which(phi2_c0Sig & phi2_gsSig) # VTX 00245 w/ Poa pra
plot(p2_c0,p2_gs,cex=0.3)

# nullmodel 2 on treatment 3
exJE3_gs <- simulate(nullmodel(exJE3_pa,"greedyqswap"),nnull,101,thin=100)
phi3_gs <- apply(exJE3_gs,3, function(x) phiOnPlotMatrix(x,pIp3))
p3_gs <- matrix(c(1-sapply(1:nrow(phi3_gs),
                           function(i) length(which(c(phi3)[i]-phi3_gs[i,]>0)))/nnull),
                nrow=nrow(phi3),
                ncol=ncol(phi3))
phi3[p.adjust(p3_gs,"fdr")<0.05] # 1 left
plot(c(phi3),col=c("black","red")[1+as.numeric(p.adjust(p3_gs,"fdr")<0.05)])
phi3_gsSig <- matrix(p.adjust(p3_gs,"fdr")<0.05,nrow=nrow(phi3),ncol=ncol(phi3),
                     dimnames=dimnames(phi3))
# VTX00166 w/ Bro ere

which(phi3_c0Sig & phi3_gsSig) # none
plot(p3_c0,p3_gs,cex=0.3)

which(phi1_gsSig & phi2_gsSig & phi3_gsSig) # 1487 1488 1496 1503 1515 1531 in all 3 treatments
                                            # VTX00061 w/ Luz cam, Phl pra, Leo his, Ver cha, Kna Arv, Tri hyb
which(phi1_gsSig & phi3_gsSig)
which(phi2_gsSig & phi3_gsSig)


# and with the 2nd 0-model, treatment 1 has less specialisations

pdf("Graphics/20230510_dBEF_7_nullmodel.pdf")
plot(p1_c0,p1_gs,cex=0.3)
plot(p2_c0,p2_gs,cex=0.3)
plot(p3_c0,p3_gs,cex=0.3)
dev.off()

## --> only two VT show specificity: extract values that go into calculations of phi
# VTX00156 w/ Gal mol in D1
# VTX00245 w/ Poa pra in D2
# in how many plots are the plants present?
table(pI$GAL_MOL[pI$treatment == "D1" & pI$SOWNDIV != 60]) # 6 / 76
table(pI$POA_PRA[pI$treatment == "D2" & pI$SOWNDIV != 60]) # 8 / 76

# in how many plots are the VT present(in that treatment) ?
colSums(exJE1 != 0) # VTX00156 25 / 76
colSums(exJE2 != 0) # VTX00245 16 / 76

Galmol = pI %>% filter(GAL_MOL == 1) %>% filter(treatment == "D1", SOWNDIV != 60)
Galmol = merge(as.data.frame(exJE), Galmol, by = 0)

Poapra = pI %>% filter(POA_PRA == 1) %>% filter(treatment == "D2", SOWNDIV != 60)
Poapra = merge(as.data.frame(exJE), Poapra, by = 0)

#### 8 add resident species ####
# let's add data on plant species-specific AMF communities (resident species)
# the VT has already been added, as the non-assigned ASV were treed together with dBEF
# let's finish up the data prep for residents species data
# nonVT assigned ASV were treed together with dBEF ASV
cusVT_res = cusVT %>% filter(str_detect(cusVT$OTU, "OTU_r_"))
# remove the string that was added to differentiate comm. from resident species ASVs
cusVT_res$OTU <- str_replace(cusVT_res$OTU, "_r_", "_") 
length(setdiff(VT_res_nVT$OTU, cusVT_res$OTU)) # 20 (singletons that wer filtered out vefore treeing)

# combine VTX und cusVT and remove double assigned VTs
VT_res_com <- bind_rows(VT_res_VTX, cusVT_res) %>% separate(VT, c("VT","two", "three", "four"), sep=";") %>% 
  select(OTU, VT)
dim(VT_res_com) # 1295 x 2

# add VT IDs to ASV table and aggregate reads at VT level
length(setdiff(row.names(ASV_res),VT_res_com$OTU)) # 76 --> 20 filtered singletons. What are the other 56 ?? 
VT_res_com = tibble::column_to_rownames(VT_res_com, var = "OTU")
ASV_VT_res <- merge(ASV_res, VT_res_com, by = 0) 
dim(ASV_VT_res) # 1295 x 376  
length(unique(ASV_VT_res$VT)) #87 VT
ASV_VT_res$VT <- as.factor(ASV_VT_res$VT)
ASV_VT_res_aggr <- ASV_VT_res %>% select(-Row.names) %>% group_by(VT) %>% summarise(across(.cols = everything(), .fns =  sum))
dim(ASV_VT_res_aggr) # 87 x 375

#prepare tax data
VT_res_VTX <- VT_res_VTX %>% tidyr::separate(VT, c("VT","B"), sep = ';', remove = TRUE) %>% select(OTU, VT) # some ASV had double hits
VT_res_VTX_tax <- merge(VT_res_VTX, VT_res, by.x = "VT", by.y = "VTX") %>% select(VT, Hitshort) %>%  distinct(VT, .keep_all = T)
VT_res_VTX_tax <- VT_res_VTX_tax %>% separate(Hitshort, c("x","Family", "Genus", "Species"), sep=" ") %>% select(-"Species", -"x")
VT_res_VTC <- cusVT_res %>% select(VT) %>% distinct(VT)
VT_res_VTC <- VT_res_VTC %>% mutate(Family = paste0("Cluster", substr(VT, 6,8)), Genus = paste0("Cluster", substr(VT, 6,8))) 
VT_res_tax <- bind_rows(VT_res_VTX_tax, VT_res_VTC) # 87 VT


# add resident plant info
res_samples <- as.data.frame(read_excel("C:/Users/albracht/Documents/R-Projekte/JE_residents/JE_res_SampleList.xlsx"))
# resident data contains some non-Main Exp. mono plots
pIm <- as.data.frame(read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/Plot_Information_Monos_2020.xlsx", sheet = "Plot_Monos_2020" ))
pIm <- dplyr::bind_rows(pIm, pI[pI$treatment == "D3",])
pIm <- pIm %>% select(-COMMENT, -history_plant, -history_soil)
res_samples <- res_samples %>% inner_join(pIm, by = "plot")
res_samples <- left_join(res_samples, res_mono_samples, by = "plot")
res_samples <- tibble::column_to_rownames(res_samples, var="Sample_ID")

# prepare phyloseq
ASV_VT_res_aggr <- ASV_VT_res_aggr %>% tibble::column_to_rownames(var = "VT")
VT_res_phy = otu_table(as.matrix(ASV_VT_res_aggr), taxa_are_rows=TRUE)
Tax_res_phy <- tibble::column_to_rownames(VT_res_tax,var = "VT")
Tax_res_phy =tax_table(as.matrix(Tax_res_phy))
Sam_res_phy =sample_data(res_samples)
#JE_res_phy@sam_data$sampletype = factor(JE_res_phy@sam_data$sampletyp, levels = c("BS","RS","RO"))

JE_res<-phyloseq(VT_res_phy, Tax_res_phy, Sam_res_phy)
JE_res@sam_data$sampletype = factor(JE_res@sam_data$sampletyp, levels = c("BS","RS","RO"))


#### 9 resident species: Bulk Soil ####
# how similar are dBEF and Res Bulk Soil communities?
setdiff(rownames(VT_res_phy), rownames(VT_dBEF_phy))

ASV_VT_dBEFres = dplyr::bind_rows(as.data.frame(t(ASV_VT_dBEF_aggr)), as.data.frame(t(ASV_VT_res_aggr))) #combined ASV table
ASV_VT_dBEFres[is.na(ASV_VT_dBEFres)] <- 0

pIall <- dplyr::bind_rows(res_samples, pI)
pIall$sampletype[is.na(pIall$sampletype)] <- "D3"
pIall$LOG_SD = log(pIall$SOWNDIV)
pIall$SOWNDIV = as.factor(pIall$SOWNDIV)

VT_dBEFres_tax = dplyr::union(VT_dBEF_tax, VT_res_tax)
VT_dBEFres_tax = distinct(VT_dBEFres_tax, VT, .keep_all = T) # 1 VT duplicated
VT_dBEFres_tax = tibble::column_to_rownames(VT_dBEFres_tax, var = "VT")


VT_dBEFres_phy = otu_table(as.matrix(ASV_VT_dBEFres), taxa_are_rows = F)
Sam_dBEFres_phy = sample_data(pIall)
Tax_dBEFres_phy = tax_table(as.matrix(VT_dBEFres_tax))

JE_dBEFres = phyloseq(VT_dBEFres_phy, Tax_dBEFres_phy, Sam_dBEFres_phy)
JE_dBEFres = subset_samples(JE_dBEFres, SOWNDIV !="60")
JE_dBEFres_BS = subset_samples(JE_dBEFres, sampletype == "BS" | sampletype == "D3")


dist_BS <- ordinate(JE_dBEFres_BS, method = "NMDS", type="samples")
BS_ordi <- plot_ordination(JE_dBEFres_BS, dist_BS, "samples", color = "sampletype", shape = "SOWNDIV") 

BS_relabu <- JE_dBEFres_BS %>% 
  tax_glom(taxrank="Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt() %>% 
  filter(Abundance>0.01) %>% 
  arrange(Genus) 

unique(BS_relabu$Genus) # 12 genera

# plot relative abundances
BS_relabu_plot <- ggplot(BS_relabu, aes(x=as.factor(SOWNDIV), y=Abundance, fill=Genus)) +
  geom_bar(stat="identity", position="fill", show.legend=TRUE) +
  #scale_fill_manual(values=mycols, name="Genus") +
  scale_fill_paletteer_d("ggthemes::Classic_20", name ="Genus") +
  #scale_fill_paletteer_d("colorBlindness::paletteMartin", name = "Genus") +
  ylab("Relativ Abundance") + xlab("Plant Species Richness") +
  theme(legend.position="right", legend.text=element_text(size=6), legend.title=element_text(size=10), legend.key.size=unit(0.5, 'cm')) 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

pdf("Graphics/20230424_dBEF_9_BSres.pdf")
ggarrange(BS_ordi, BS_relabu_plot + coord_flip() + facet_grid(sampletype~.))
dev.off()



#### 10 resident species: Rhizosphere ####
# Do we find same species / patterns in Rhizosphere?

##### 10.1 Specialization of AMF in Rhizosphere / Roots #####
#pI <- read.delim("Plot_Information_DeltaBEF_2020.txt",row.names = 1)
pIres <- res_samples[,16:74] 
# I seem to have only 59 plant species, although I thought there are 60?
# exJE <- read.delim("JEdBEF_fVT_v1.txt",row.names = 1) --> Phi on unfiltered data
exJEres <- as.data.frame(t(otu_table(JE_res)))

pIres <- pIres[sapply(rownames(exJEres),function(x) which(rownames(pIres)==x)),]
all(rownames(pIres)==rownames(exJEres))

pIresBS <- pIp[grep("BS$",rownames(pIres)),]
pIresRS <- pIp[grep("RS$",rownames(pIres)),]
pIresRO <- pIp[grep("RO$",rownames(pIres)),]

exJEresBS <- exJEres[grep("BS$",rownames(exJEres)),]
exJEresRS <- exJEres[grep("RS$",rownames(exJEres)),]
exJEresRO <- exJEres[grep("RO$",rownames(exJEres)),]

phiBS <- phiOnPlotMatrix(exJEresBS,pIresBS)
phiRS <- phiOnPlotMatrix(exJEresRS,pIresRS)
phiRO <- phiOnPlotMatrix(exJEresRO,pIresRO)

sort(c(phiBS), decreasing = T) # 27 > 0.34
specBS <- phiBS > 0.34
table(apply(specBS, 2, function(x) length(which(x)))) # per VT how many times > 0.37 ?
#  0  1  2  3  5  9 
# 76  6  2  1  1  1
sort(apply(specBS, 2, function(x) length(which(x))), decreasing = T) # which VT are these?
# VTX00156 VTC00004 VTX00052 VTX00114 VTX00338 
# 9        5        3        2        2 
# VTC00014 VTC00019 VTX00005 VTX00135 VTX00354 
# 1        1        1        1        1 
# VTX00356  
# 1        
sort(c(phiRS), decreasing = T) # 31 > 0.30
specRS <- phiRS > 0.30
table(apply(specRS, 2, function(x) length(which(x)))) # per VT how many times > 0.37 ?
#  0  1  2  3  4  5 
# 73  7  1  3  2  1 
sort(apply(specRS, 2, function(x) length(which(x))), decreasing = T) 
# VTX00356 VTX00052 VTX00338 VTC00004 VTX00281    --> VTX00052 w/ Gle hed, Arr ela, Pas sat and Ant syl
# 5        4        4        3        3 
# VTX00357 VTX00212 VTX00005 VTX00115 VTX00130    --> VTX00357 w/ Cyn cri, Pri ver, Med lup
# 3        2        1        1        1 
# VTX00140 VTX00155 VTX00222 VTX00392             --> VTX00140 w/ Vic cra
# 1        1        1        1 

sort(c(phiRO), decreasing = T) # 37 > 0.36
specRO <- phiRO > 0.36
table(apply(specRO, 2, function(x) length(which(x)))) # per VT how many times > 0.37 ?
# 0  1  2  3  5  6 
# 70  9  1  5  1  1 
sort(apply(specRO, 2, function(x) length(which(x))), decreasing = T) 
# VTX00247 VTX00223 VTC00002 VTC00004 VTX00005 
# 6        5        3        3        3 
# VTX00049 VTX00052 VTC00006 VTC00009 VTX00129   --> VTX00052 w/ Gle hed, Arr ela, Pas sat
# 3        3        2        1        1 
# VTX00137 VTX00140 VTX00155 VTX00301 VTX00304   --> VTX00140 w/ Dac glo
# 1        1        1        1        1 
# VTX00357 VTX00444                              --> VTX00357 w/ Med lup
# 1        1 

##### 10.2 turnover of species BS > RS > RO #####
JE_res_f = filter_taxa(JE_res, function(x) sum(x > 0) > (0.2 * length(x)), TRUE)
sample_sums(JE_res_f)[order(sample_sums(JE_res_f), decreasing = T)] # 2279
JE_res_frf <- rarefy_even_depth(JE_res_f, rngseed=1, sample.size=2200 , replace=FALSE) # no OTUs removed

res_VT_frf_long = t(as.data.frame(otu_table(JE_res_frf))) 
res_VT_frf_long = as.data.frame(res_VT_frf_long) %>% 
  tibble::rownames_to_column(var="sample") %>% 
  pivot_longer(!sample, names_to = "VTX", values_to = "relabu")
res_VT_frf_long = res_VT_frf_long %>% 
  mutate(comp = ifelse(str_ends(sample, "BS"),"1", 
                        ifelse(str_ends(sample, "RS"),"2","3"))) %>% 
  mutate(comp = as.numeric(comp)) 
res_VT_frf_long$sample = substr(res_VT_frf_long$sample,1,13) #reduce 'plot' to just plot code w/o treatment
res_VT_frf_long$relabu = log(res_VT_frf_long$relabu)

library(codyn)
res_spec_turnover <- turnover(res_VT_frf_long, 
                               time.var = "comp",
                               species.var = "VTX",
                               abundance.var = "relabu",
                               replicate.var = "sample",
                               metric = "total")
res_spec_disappearance <- turnover(res_VT_frf_long, 
                                    time.var = "comp",
                                    species.var = "VTX",
                                    abundance.var = "relabu",
                                    replicate.var = "sample",
                                    metric = "disappearance")
res_spec_appearance <- turnover(res_VT_frf_long, 
                                 time.var = "comp",
                                 species.var = "VTX",
                                 abundance.var = "relabu",
                                 replicate.var = "sample",
                                 metric = "appearance")
#Format a compiled data frame
res_spec_turnover$metric<-"total"
names(res_spec_turnover)[1]="turnover"

res_spec_appearance$metric<-"appearance"
names(res_spec_appearance)[1]="turnover"

res_spec_disappearance$metric<-"disappearance"
names(res_spec_disappearance)[1]="turnover"

res_spec_allturnover<-rbind(res_spec_turnover, res_spec_appearance, res_spec_disappearance) 
#add meta data
res_spec_allturnover<-res_spec_allturnover %>% merge(res_samples, by = "plot") %>% 
  select(plot, turnover, treat, metric, SOWNDIV, BLOCK, FUNCGR, GRASS, LEG, HERB)

#Create the graph
pdf("Graphics/20230424_res_6_turnover.pdf")
#turnover_graph <- 
ggplot(res_spec_allturnover, aes(x=treat, y=turnover, color=SOWNDIV)) + 
  geom_line(size = 1) + 
  facet_wrap(~metric) + 
  theme_bw() + 
  theme(legend.position="bottom")
ggplot(res_spec_allturnover, aes(x=treat, y=turnover, color=metric)) + 
  geom_line(size = 1) + 
  facet_wrap(~SOWNDIV) + 
  theme_bw() + 
  theme(legend.position="bottom")
ggplot(res_spec_allturnover, aes(x=treat, y=turnover, color=metric)) + 
  geom_line(size = 1) + 
  facet_wrap(~metric) + 
  theme_bw() + 
  theme(legend.position="bottom")
dev.off()

#### 11. Environemntal Data ####

##### 11.1 Plant Above-/Belowground Biomass #####
pI_pbm21 <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/envData/JE_dBEF_PlantBiomass_2021.xlsx", sheet="Data")
pI_rbm21 <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/envData/JE_dBEF_RootBiomass_2021.xlsx", sheet="Data")

plant_spec <- read.csv("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/SpeciesList.csv", sep=";")


pI_pbm21 <- left_join(pI_pbm21, plant_spec, by="species_other") %>% 
  mutate(plot_treat = paste0(plotcode, treatment)) %>% mutate(FG = replace_na(FG,"weed")) %>% mutate(full_name = replace_na(full_name,"Weed"))
pI_pbm21$biomass = as.numeric(pI_pbm21$biomass)

grp <- c('plot_treat', 'FG')
pI_pbm21_sum <- pI_pbm21 %>% filter(season =="May" | subsample == "1")  %>% #what is the subsample? 
  group_by(across(all_of(grp))) %>% summarise(bm_sum = sum(biomass))
pI_pbm21_sum <- merge(pI_pbm21_sum, dBEF_samples, by.x="plot_treat", by.y = 0)


pI_env <- pI_pbm21_sum  %>% pivot_wider(names_from = FG, values_from=bm_sum) %>% 
  replace(is.na(.),0) %>% 
  rename(BM_dead = dead, BM_grass = grass, BM_leg = leg, BM_sherb = sherb,BM_therb = therb, BM_weed = weed) %>% 
  mutate(BM_total = BM_dead + BM_grass + BM_leg + BM_sherb + BM_therb + BM_weed) %>% 
  mutate(BM_target = BM_grass + BM_leg + BM_sherb + BM_therb)

pdf("Graphics/20231002_dBEF_11_pbm.pdf")
ggplot(pI_env, aes(x = treatment.x , y = BM_target, color = SOWNDIV))  + 
  geom_violin()  + facet_grid(.~SOWNDIV) + theme_classic() + ylab("Aboveground Biomass of Target Species [g]") +
  xlab("Treatment") +
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

ggplot(pI_env, aes(x = LOG_SD, y = BM_target, color = treatment.x))  + 
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ylab("Aboveground Biomass of Target Species [g]") +
  xlab("Treatment") +
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f"))
  
pI_env %>% select(plot.x, treatment.x, SOWNDIV, BM_target, BM_grass, BM_leg, BM_herb, BM_weed, BM_dead, root.bm) %>% 
  pivot_longer(cols = -c(plot.x,treatment.x,SOWNDIV), names_to = "FG", values_to = "BM") %>% 
  mutate(FG = factor(FG, levels = c("BM_target", "BM_grass", "BM_leg", "BM_herb", "BM_dead", "root.bm"))) %>% 
  ggplot(aes(x = treatment.x, y = BM, color = SOWNDIV))  + 
  geom_violin()+ ylim(0, 1500) + facet_grid(FG~SOWNDIV) + theme_classic() + 
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))

pI_env %>% select(plot.x, treatment.x, LOG_SD, BM_target, BM_grass, BM_leg, BM_herb) %>% #, BM_weed, BM_dead, root.bm) %>% 
  pivot_longer(cols = -c(plot.x,treatment.x,LOG_SD), names_to = "FG", values_to = "BM") %>% 
  mutate(FG = factor(FG, levels = c("BM_target", "BM_grass", "BM_leg", "BM_herb"))) %>% 
  ggplot(aes(x = LOG_SD, y = BM, color = treatment.x))  + 
  geom_point()+geom_smooth(method = "lm")+ ylim(0, 1500) + theme_classic() + facet_grid(.~FG) +
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f")) 

dev.off()

#standardized biomass
pI_env_std <- pI_env %>% select(plot.x, treatment.x, LOG_SD, BM_target, BM_grass, BM_leg, BM_herb) %>% 
  group_by(treatment.x) %>% mutate(BM_target_std = as.numeric(scale(BM_target))) %>% 
  mutate(BM_leg_std = as.numeric(scale(BM_leg))) %>% mutate(BM_grass_std = as.numeric(scale(BM_grass))) %>% 
  mutate(BM_herb_std = as.numeric(scale(BM_herb)))

pdf("Graphics/20231005_dBEF_11_pbm_std.pdf")
ggplot(pI_env_std, aes(x = LOG_SD, y = BM_target_std, color = treatment.x))  + 
  geom_point() + geom_smooth(method = "lm") + theme_classic() + ylab("Standardized Aboveground Biomass of Target Species [g]") +
  xlab("Treatment") +
  theme(legend.position="bottom", panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  scale_color_manual(values=c("#01665e", "#91cf60", "#fec44f"))
dev.off()

# Is Plant Biomass influenced by Block, treatment or Plant Diversity? 
summary(lm(BM_target ~ BLOCK, data= pI_pbm21_wide)) # ns
anova(lme(BM_target ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env)) # Div is significant p = 4.55e-07
lm(BM_target ~ LOG_SD, data = pI_env[pI_env$treatment.x == "D1",]) # Intercept = 223.4  Log_SD = 187.6
lm(BM_target ~ LOG_SD, data = pI_env[pI_env$treatment.x == "D2",]) # Intercept = 173.8  Log_SD = 128.6
lm(BM_target ~ LOG_SD, data = pI_env[pI_env$treatment.x == "D3",]) # Intercept = 223.4  Log_SD = 189.9

lm(BM_target_std ~ LOG_SD, data = pI_env_std[pI_env_std$treatment.x == "D1",]) # Intercept = -0.745  Log_SD = 0.537
lm(BM_target_std ~ LOG_SD, data = pI_env_std[pI_env_std$treatment.x == "D2",]) # Intercept = -0.51  Log_SD = 0.369
lm(BM_target_std ~ LOG_SD, data = pI_env_std[pI_env_std$treatment.x == "D3",]) # Intercept = -0.81  Log_SD = 0.586


mean(pI_env[pI_env$treatment.x == "D3",]$BM_target, na.rm = T)
# 398.4436
mean(pI_env[pI_env$treatment.x == "D2",]$BM_target)
# 352.0966
mean(pI_env[pI_env$treatment.x == "D1",]$BM_target)
# 483.5287

ggplot(pI_env, aes(treatment.x, BM_target)) + geom_boxplot()
ggplot(pI_env, aes(SOWNDIV, BM_target)) + geom_boxplot()

summary(lm(BM_target ~ FUNCGR, data = pI_env )) # bm sig higher with 3-4 FUNCGR
summary(lm(BM_target ~ LEG, data = pI_env )) # bm sig higher with legumes present
summary(lm(BM_target ~ HERB, data = pI_env ))# ns
summary(lm(BM_target ~ GRASS, data = pI_env )) # bm sig higher with grasses present

anova(lme(BM_target ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env )) # total includes weeds | Div is sig
anova(lme(BM_dead ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_weed ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_grass ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_herb ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_leg ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env ))


anova(lme(BM_total ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB, random = ~1|BLOCK/plot.x , data = pI_env )) # total includes weeds | Div is sig
anova(lme(BM_grass ~ treatment.x + LOG_SD + FUNCGR + LEG + HERB, random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_herb ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS , random = ~1|BLOCK/plot.x, data = pI_env ))
anova(lme(BM_leg ~ treatment.x + LOG_SD + FUNCGR + GRASS + HERB, random = ~1|BLOCK/plot.x, data = pI_env ))

ggplot(pI_env[pI_env$LEG == 1,], aes(treatment.x, BM_leg)) + geom_boxplot()
ggplot(pI_env[pI_env$LEG == 1,], aes(SOWNDIV, BM_leg)) + geom_boxplot()
ggplot(pI_env[pI_env$GRASS == 1,], aes(treatment.x, BM_grass)) + geom_boxplot()
ggplot(pI_env[pI_env$GRASS == 1,], aes(SOWNDIV, BM_grass)) + geom_boxplot()
ggplot(pI_env[pI_env$HERB == 1,], aes(treatment.x, BM_therb+BM_sherb)) + geom_boxplot()
ggplot(pI_env[pI_env$HERB == 1,], aes(SOWNDIV, BM_therb+BM_sherb)) + geom_boxplot()

# Is root biomass affected by treatment or plant diversity? 
pI_rbm21 <- pI_rbm21 %>%  mutate(plot_treat = paste0(plot, treatment)) %>% 
   filter(month == "6") %>%  filter(depth.min == 0.00)
pI_env = left_join(pI_env, pI_rbm21, by = "plot_treat")

anova(lm(root.bm ~ BLOCK, data= pI_env)) # ns
anova(lm(root.bm ~ month, data= pI_env)) #ns
anova(lme(root.bm ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit ))
anova(lme(root.bm ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB , random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit )) 


anova(lm(fine_root.bm ~ BLOCK, data= pI_env)) # ns
anova(lme(fine_root.bm ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit ))
anova(lme(fine_root.bm ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB , random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit )) 


anova(lme(coarse_root.bm ~ BLOCK, data= pI_env)) # ns
anova(lme(coarse_root.bm ~ month, data= pI_env)) # 
anova(lme(coarse_root.bm ~ treatment.x * LOG_SD , random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit )) 
anova(lme(coarse_root.bm ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB , random = ~1|BLOCK/plot.x, data = pI_env, na.action = na.omit )) 


##### 11.2 Soil Data #####
pI_resp21 <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/envData/JE_dBEF_MicrobRespiration_2021.xlsx", sheet="Data")
pI_sP21 <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/envData/JE_dBEF_SoilP_2021.xlsx", sheet="Data")
pI_sSWC21 <- read_excel("C:/Users/albracht/Documents/R-Projekte/JE_dBEF/envData/JE_dBEF_SWC_2021.xlsx", sheet="Data")
# SWC data contains only D1 and D2 treatmnents

pI_resp21 <- pI_resp21 %>% rename(plot_treat = plot, SWC_res = soil_water_content) %>% select(-year) 
pI_sP21 <- pI_sP21 %>% mutate(plot_treat = paste0(plot, treatment)) %>% select(plot_treat, SWC, PO4_P) %>% rename(SWC_P = SWC)
pI_sSWC21 <- pI_sSWC21 %>% filter(str_detect(pI_sSWC21$date, "-05-30")) %>%  #|-06-")) %>%  #data for every week
  mutate(plot_treat = paste0(plot, treatment))  %>% 
  select(plot_treat,Vol10, Vol20) %>% rename(SWC_10 = Vol10, SWC_20 = Vol20) %>% mutate_at(c('SWC_10', 'SWC_20'), as.numeric)
pI_sSWC21 <- pI_sSWC21 %>% group_by(plot_treat) %>% summarise(SWC_10 = mean(SWC_10), SWC_20 = mean(SWC_20))

pI_env = left_join(pI_env, pI_sP21, by = "plot_treat")
pI_env = left_join(pI_env, pI_sSWC21, by = "plot_treat")
pI_env = left_join(pI_env, pI_resp21, by = "plot_treat")
pI_env$BM_herb = pI_env$BM_sherb + pI_env$BM_therb

pdf("Graphics/20230510_dBEF_11_soil.pdf")
ggarrange(
ggplot(pI_env, aes(x = LOG_SD , y = SWC_10, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("SWC in 0-10 cm depth"),
ggplot(pI_env, aes(x = LOG_SD , y = SWC_20, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("SWC in 10-20 cm depth"),
ggplot(pI_env, aes(x = LOG_SD , y = SWC_P, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("SWC (gravimetric)"),
ggplot(pI_env, aes(x = LOG_SD , y = SWC_res, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("SWC (respiration)"),

ggplot(pI_env, aes(x = LOG_SD , y = PO4_P, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("Soil P"),

ggplot(pI_env, aes(x = LOG_SD , y = basal_respiration, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("Basal Microbial Soil Respiration"),
ggplot(pI_env, aes(x = LOG_SD , y = soil_microbial_biomass_C, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("Microbial Biomass C"),
ggplot(pI_env, aes(x = LOG_SD , y = respiratory_quotient, group = treatment))  + 
  geom_point()  + facet_grid(.~treatment) + geom_smooth(method = lm) + ggtitle("Respiratory Quotient"),

nrow = 4, ncol = 2, common.legend = T)
dev.off()

# Is microbial respiration different across treatments and/or Plant Diversity? 
summary(lm(basal_respiration ~ BLOCK, data= pI_env)) # sig
anova(lme(basal_respiration ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x,data= pI_env)) 
anova(lme(basal_respiration ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB, random = ~1|BLOCK/plot.x, data = pI_env ))

summary(lm(soil_microbial_biomass_C ~ BLOCK, data= pI_env)) # B2/B3 higher than B1, B4 not
anova(lme(soil_microbial_biomass_C ~ treatment.x * LOG_SD, random = ~1|BLOCK/plot.x, data= pI_env)) # D2, treat/Div interaction sig
anova(lme(soil_microbial_biomass_C ~ treatment.x + LOG_SD + FUNCGR + LEG + GRASS + HERB, random = ~1|BLOCK/plot.x, data = pI_env ))


summary(lm(respiratory_quotient ~ BLOCK, data= pI_env)) # B3 / B4 sig higher
summary(lm(respiratory_quotient ~ treatment * LOG_SD, data= pI_env)) # ns

summary(lm(respiratory_quotient ~ LEG, data= pI_env)) # ns
summary(lm(respiratory_quotient ~ GRASS, data= pI_env)) # sig higher with grasses present
summary(lm(respiratory_quotient ~ HERB, data= pI_env)) # ns


# Is soil influenced by Block, treatment or Plant Diversity? 
summary(lm(PO4_P ~ BLOCK, data= pI_env)) # Block 4 sig lower
summary(lm(PO4_P ~ BLOCK*treatment.x, data= pI_env))
summary(lm(SWC_P ~ BLOCK, data= pI_env)) # Block 3 and 4 significantly wetter
summary(lm(PO4_P~SWC_P, data = pI_env))


summary(lm(PO4_P ~ treatment.x * LOG_SD, data= pI_env)) # P sig higher in D1
summary(lm(SWC ~ treatment * LOG_SD, data= pI_env)) # SWC sig lower in D2 and D3, in D3 sig interaction of treat*Div

summary(lm(PO4_P ~ LEG, data= pI_env)) # P sig lower when LEG present
summary(lm(SWC ~ LEG, data= pI_env)) # ns
summary(lm(PO4_P ~ GRASS, data= pI_env)) # P sig higher when GRASS present
summary(lm(SWC_P ~ GRASS, data= pI_env)) # SWC sig higher when GRASS present
summary(lm(PO4_P ~ HERB, data= pI_env)) # ns
summary(lm(SWC_P ~ HERB, data= pI_env)) # ns

# different SWC correlated?
summary(lm(SWC_P ~ SWC_res, data = pI_env)) #yes
summary(lm(SWC_P ~ SWC_10, data = pI_env)) #yes
summary(lm(SWC_P ~ SWC_20, data = pI_env)) #yes

summary(lm(SWC_res ~ SWC_10, data = pI_env)) #yes
summary(lm(SWC_res ~ SWC_20, data = pI_env)) #no
summary(lm(SWC_10 ~ SWC_20, data = pI_env)) #yes

# different SWC against Experiment
summary(lm(SWC_res ~ BLOCK, data= pI_env)) # all sig higher than B1
summary(lm(SWC_P ~ BLOCK, data= pI_env)) # B3/B4 sig higher, B2 lower than B1
summary(lm(SWC_10 ~ BLOCK, data= pI_env)) # all sig higher than B1
summary(lm(SWC_20 ~ BLOCK, data= pI_env)) # ns

summary(lm(SWC_res ~ treatment * LOG_SD, data= pI_env)) # ns
summary(lm(SWC_P ~ treatment * LOG_SD, data= pI_env)) # treat sig, in D3 interaction treat/Div sig
summary(lm(SWC_10 ~ treatment * LOG_SD, data= pI_env)) # ns
summary(lm(SWC_20 ~ treatment * LOG_SD, data= pI_env)) # ns

summary(lm(SWC_res ~ LEG, data= pI_env)) # ns
summary(lm(SWC_res ~ GRASS, data= pI_env)) # ns
summary(lm(SWC_res ~ HERB, data= pI_env)) # ns

summary(lm(SWC_P ~ LEG, data= pI_env)) # ns
summary(lm(SWC_P ~ GRASS, data= pI_env)) # sig higher with grasses present
summary(lm(SWC_P ~ HERB, data= pI_env)) # ns

summary(lm(SWC_10 ~ LEG, data= pI_env)) # ns
summary(lm(SWC_10 ~ GRASS, data= pI_env)) # sig lower with grasses present
summary(lm(SWC_10 ~ HERB, data= pI_env)) # ns

summary(lm(SWC_20 ~ LEG, data= pI_env)) # ns
summary(lm(SWC_20 ~ GRASS, data= pI_env)) # sig lower with grasses present
summary(lm(SWC_20 ~ HERB, data= pI_env)) # ns

##### 11.5 Are Env Parameter correlated? #####
# is there a correlation between the different environmental data?
# Do coares and fine root biomasses correlate?
summary(lm(coarse_root.bm ~ fine_root.bm , data = pI_env)) #no

# Do plant biomasses of functional groups correlate?
summary(lm(BM_target ~ BM_leg, data = pI_env)) #yes
summary(lm(BM_target ~ BM_grass, data = pI_env)) #yes
summary(lm(BM_target ~ BM_sherb, data = pI_env)) #yes
summary(lm(BM_target ~ BM_therb, data = pI_env)) #yes
summary(lm(BM_sherb ~ BM_therb, data = pI_env)) #yes
summary(lm(BM_sherb ~ BM_leg, data = pI_env)) #np
summary(lm(BM_sherb ~BM_grass, data = pI_env)) #yes
summary(lm(BM_therb ~ BM_leg, data = pI_env)) #no
summary(lm(BM_therb ~ BM_grass, data = pI_env)) #yes
summary(lm(BM_grass ~ BM_leg, data = pI_env)) #no  --> BM leg does not correlate w/ other FG

# Do root and aboveground biomasses correlate?
summary(lm(root.bm ~ BM_target, data = pI_env)) #yes
summary(lm(coarse_root.bm ~ BM_target, data = pI_env)) # no
summary(lm(fine_root.bm ~ BM_target, data = pI_env)) #yes

summary(lm(root.bm ~ BM_leg, pI_env)) #no
summary(lm(coarse_root.bm ~ BM_leg, pI_env)) #no
summary(lm(fine_root.bm ~ BM_leg, pI_env)) #no

summary(lm(root.bm ~ BM_grass, pI_env)) #yes
summary(lm(coarse_root.bm ~ BM_grass, pI_env)) #no
summary(lm(fine_root.bm ~ BM_grass, pI_env)) #yes

summary(lm(root.bm ~ BM_sherb, pI_env)) #no
summary(lm(coarse_root.bm ~ BM_sherb, pI_env)) #no
summary(lm(fine_root.bm ~ BM_sherb, pI_env)) #no

summary(lm(root.bm ~ BM_therb, pI_env)) #no
summary(lm(coarse_root.bm ~ BM_therb, pI_env)) #yes
summary(lm(fine_root.bm ~ BM_therb, pI_env)) #no
# --> coarse root biomass driven by therbs, fine and total root bm driven by grasses

# Do soil parameteres correlate? 
summary(lm(PO4_P ~ SWC_P, pI_env)) #no
summary(lm(PO4_P ~ SWC_res, pI_env))#no
summary(lm(PO4_P ~ SWC_10, pI_env))#yes
summary(lm(PO4_P ~ SWC_20, pI_env))#no

summary(lm(basal_respiration ~ SWC_P, pI_env)) #yes
summary(lm(basal_respiration ~ SWC_res, pI_env)) #yes
summary(lm(basal_respiration ~ SWC_10, pI_env)) #yes
summary(lm(basal_respiration ~ SWC_20, pI_env)) #no

summary(lm(respiratory_quotient ~ SWC_P, pI_env)) #yes
summary(lm(respiratory_quotient ~ SWC_res, pI_env)) #yes
summary(lm(respiratory_quotient ~ SWC_10, pI_env)) #yes
summary(lm(respiratory_quotient ~ SWC_20, pI_env)) #no

summary(lm(soil_microbial_biomass_C ~ SWC_P, pI_env)) #yes
summary(lm(soil_microbial_biomass_C ~ SWC_res, pI_env)) #yes
summary(lm(soil_microbial_biomass_C ~ SWC_10, pI_env)) #no
summary(lm(soil_microbial_biomass_C ~ SWC_20, pI_env)) #no

summary(lm(PO4_P ~ basal_respiration, pI_env)) #no
summary(lm(PO4_P ~ respiratory_quotient, pI_env))#yes
summary(lm(PO4_P ~ soil_microbial_biomass_C, pI_env))#yes

# is Plant Biomass affected by Soil ?
summary(lm(BM_target ~ PO4_P, data = pI_env)) #marginal yes
summary(lm(BM_target ~ SWC_P, data = pI_env)) #yes

summary(lm(root.bm ~ PO4_P, data = pI_env)) #yes
summary(lm(root.bm ~ SWC_P, data = pI_env)) #yes




# Which soil factor affects AMF alpha-Div?
summary(lm(Observed ~ PO4_P, pI_env))#yes
summary(lm(Shannon ~ PO4_P, pI_env))#yes
summary(lm(Evenness ~ PO4_P, pI_env))#no

summary(lm(Observed ~ SWC_res, pI_env))#no
summary(lm(Shannon ~ SWC_res, pI_env))#no
summary(lm(Evenness ~ SWC_res, pI_env))#no

# Does plant biomass affect AMF alpha-Div?
summary(lm(Observed ~ BM_target, pI_env))#no
summary(lm(Shannon ~ BM_target, pI_env))#no
summary(lm(Evenness ~ BM_target, pI_env))#no

summary(lm(Observed ~ BM_target, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Shannon ~ BM_target, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Evenness ~ BM_target, pI_env[pI_env$treatment.x == "D1",]))#no

summary(lm(Observed ~ BM_target, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Shannon ~ BM_target, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Evenness ~ BM_target, pI_env[pI_env$treatment.x == "D2",]))#no

summary(lm(Observed ~ BM_target, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Shannon ~ BM_target, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Evenness ~ BM_target, pI_env[pI_env$treatment.x == "D3",]))#no

summary(lm(Observed ~ BM_leg, pI_env))#yes
summary(lm(Observed ~ BM_leg, pI_env[pI_env$treatment.x == "D1",]))#yes
summary(lm(Observed ~ BM_leg, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Observed ~ BM_leg, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Shannon ~ BM_leg, pI_env))#yes
summary(lm(Shannon ~ BM_leg, pI_env[pI_env$treatment.x == "D1",]))#yes
summary(lm(Shannon ~ BM_leg, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Shannon ~ BM_leg, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Evenness ~ BM_leg, pI_env))#no
summary(lm(Evenness ~ BM_leg, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Evenness ~ BM_leg, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Evenness ~ BM_leg, pI_env[pI_env$treatment.x == "D3",]))#no

summary(lm(Observed ~ BM_grass, pI_env))#yes
summary(lm(Observed ~ BM_grass, pI_env[pI_env$treatment.x == "D1",]))#yes
summary(lm(Observed ~ BM_grass, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Observed ~ BM_grass, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Shannon ~ BM_grass, pI_env))#no
summary(lm(Shannon ~ BM_grass, pI_env[pI_env$treatment.x == "D1",]))#yes
summary(lm(Shannon ~ BM_grass, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Shannon ~ BM_grass, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Evenness ~ BM_grass, pI_env))#no
summary(lm(Evenness ~ BM_grass, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Evenness ~ BM_grass, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Evenness ~ BM_grass, pI_env[pI_env$treatment.x == "D3",]))#no

summary(lm(Observed ~ BM_sherb, pI_env))#yes
summary(lm(Observed ~ BM_sherb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Observed ~ BM_sherb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Observed ~ BM_sherb, pI_env[pI_env$treatment.x == "D3",]))#yes
summary(lm(Shannon ~ BM_sherb, pI_env))#yes
summary(lm(Shannon ~ BM_sherb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Shannon ~ BM_sherb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Shannon ~ BM_sherb, pI_env[pI_env$treatment.x == "D3",]))#yes
summary(lm(Evenness ~ BM_sherb, pI_env))#no
summary(lm(Evenness ~ BM_sherb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Evenness ~ BM_sherb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Evenness ~ BM_sherb, pI_env[pI_env$treatment.x == "D3",]))#yes

summary(lm(Observed ~ BM_therb, pI_env))#no
summary(lm(Observed ~ BM_therb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Observed ~ BM_therb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Observed ~ BM_therb, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Shannon ~ BM_therb, pI_env))#no
summary(lm(Shannon ~ BM_therb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Shannon ~ BM_therb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Shannon ~ BM_therb, pI_env[pI_env$treatment.x == "D3",]))#no
summary(lm(Evenness ~ BM_therb, pI_env))#no
summary(lm(Evenness ~ BM_therb, pI_env[pI_env$treatment.x == "D1",]))#no
summary(lm(Evenness ~ BM_therb, pI_env[pI_env$treatment.x == "D2",]))#no
summary(lm(Evenness ~ BM_therb, pI_env[pI_env$treatment.x == "D3",]))#no

summary(lm(Observed ~ BM_dead, pI_env))#no
summary(lm(Shannon ~ BM_dead, pI_env))#no
summary(lm(Evenness ~ BM_dead, pI_env))#no

summary(lm(Observed ~ BM_weed, pI_env))#yes
summary(lm(Shannon ~ BM_weed, pI_env))#yes
summary(lm(Evenness ~ BM_weed, pI_env))#yes
# --> richness driven by leg, grass, sherb and weeds
# --> Shannon driven by leg, sherb, weed
# --> Evenness driven by weed
# biomass of plant tendentially decrease with age (ns)


#### 12 Mediation ####
install.packages("LDM_5.0.tar.gz", repos=NULL)
library(LDM)
# build mediation formula
# otu.table | (confounders) ~ (first set of covariates) + (second set) + (last set)
# otu.table | (confounders) ~ (factors: treatment+Div) + (have influence on e.g BM)
# otu.table | (treatment + PLant Div) ~ (plant biomass + root biomass) +
#                                       (Soil P + SWC) +
#                                       (Micr Respiration + Sil Micr BM C + Resp Q)
pI_env = pI_env %>% tibble::column_to_rownames(var = "plot_treat")

med_P <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD) + (PO4_P), 
                      data = pI_env, test.mediation = TRUE, 
                      dist.method = "bray", seed = 4321)
med_P$R.squared
# cov1  (exposure)      cov2 (outcome)
#        0.101519820    0.008422292
med_P$F.statistics
# cov1     8.4784447    |  cov2   0.7033891
med_P$p.permanova
# cov1       0.00019996   | cov2  0.02040000
med_P$med.p.permanova # 0.0204

med_SWC <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD) + (SWC_P), 
                     data = pI_env, test.mediation = TRUE, 
                     dist.method = "bray", seed = 4321)
med_SWC$R.squared
# cov1  (exposure)      cov2 (outcome)
#       0.10151982      0.02685448  
med_SWC$F.statistics
# cov1     8.657738   |  cov2   2.290184
med_SWC$p.permanova
# cov1      0.00019996    | cov2  0.00019996
med_SWC$med.p.permanova # 0.00019996

## PERMANOVAFL
# how much of plant biomass is mediated by AMF ? 
# [treatment+Diversity (+ Soil-P + SWC)] -> [AMF] -> [Plant Biomass]
med_bm <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                                (BM_target), #+BM_leg+BM_grass+BM_sherb+BM_therb),  
                                 #+ (coarse_root.bm+fine_root.bm) +  
                                #(basal_respiration+soil_microbial_biomass_C),
                              data = pI_env, test.mediation = TRUE, 
                              dist.method = "bray", seed = 4321)
med_bm$R.squared
# cov1  (exposure)      cov2 (outcome)
# 0.136031954           0.004965351 
med_bm$F.statistics
# cov1  6.9995268     |  cov2 0.2554922  
med_bm$p.permanova
# cov1   0.00019996     | cov2 0.24120000 
med_bm$med.p.permanova # 0.2412

med_bm2 <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD) +
                        (BM_target), #+BM_leg+BM_grass+BM_sherb+BM_therb),  
                      #+ (coarse_root.bm+fine_root.bm) +  
                      #(basal_respiration+soil_microbial_biomass_C),
                      data = pI_env, test.mediation = TRUE, 
                      dist.method = "bray", seed = 4321)
med_bm2$R.squared
# cov1  (exposure)      cov2 (outcome)
# 0.101519820            0.004538899  
med_bm2$F.statistics
# cov1  8.4416133      |  cov2 0.3774202   
med_bm2$p.permanova
# cov1   0.00019996      | cov2 0.34040000  
med_bm2$med.p.permanova # 0.3404

med_bm3 <- permanovaFL(formula = maas_asv | confound ~ (PO4_P+SWC_P) +
                        (BM_target), #+BM_leg+BM_grass+BM_sherb+BM_therb),  
                      #+ (coarse_root.bm+fine_root.bm) +  
                      #(basal_respiration+soil_microbial_biomass_C),
                      data = pI_env, test.mediation = TRUE, 
                      dist.method = "bray", seed = 4321)
med_bm3$R.squared
# cov1  (exposure)      cov2 (outcome)
# 0.063638481             0.007994527   
med_bm3$F.statistics
# cov1  7.6774702       |  cov2 0.9644753    
med_bm3$p.permanova
# cov1   0.00019996      | cov2 0.03520000   
med_bm3$med.p.permanova # 0.0352


# per plant functional group
# legumes
maas_asv_l <- maas_asv[pI_env$LEG == "1",]
med_bml <- permanovaFL(formula = maas_asv_l | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                        (BM_leg), 
                      data = pI_env[pI_env$LEG == "1",], test.mediation = TRUE, 
                      dist.method = "bray", seed = 4321)
med_bml$R.squared
# cov1 (exp)   0.16349713      | cov2 (out) 0.01455302    
med_bml$F.statistics
# cov1  4.3761025       |  cov2 0.3895205    
med_bml$p.permanova
# cov1   0.00019996      | cov2 0.03800000   
med_bml$med.p.permanova # 0.038

med_bml2 <- permanovaFL(formula = maas_asv_l | confound ~ (treatment.x+LOG_SD) +
                         (BM_leg), 
                       data = pI_env[pI_env$LEG == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bml2$R.squared
# cov1 (treat)   0.12369338      | cov2 (div) 0.01475931    
med_bml2$F.statistics
# cov1  5.3599914       |  cov2 0.6395637    
med_bml2$p.permanova
# cov1   0.00019996      | cov2 0.04260000  
med_bml2$med.p.permanova # 0.0426

med_bml3 <- permanovaFL(formula = maas_asv_l | confound ~ (PO4_P+SWC_P) +
                          (BM_leg), 
                        data = pI_env[pI_env$LEG == "1",], test.mediation = TRUE, 
                        dist.method = "bray", seed = 4321)
med_bml3$R.squared
# cov1 (treat)   0.07510037      | cov2 (div) 0.01620358    
med_bml3$F.statistics
# cov1  4.669516       |  cov2 1.007490    
med_bml3$p.permanova
# cov1   0.00019996      | cov2 0.032800000 
med_bml3$med.p.permanova # 0.0328

# grasses
maas_asv_g <- maas_asv[pI_env$GRASS == "1",]
med_bmg <- permanovaFL(formula = maas_asv_g | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                         (BM_grass), 
                       data = pI_env[pI_env$GRASS == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmg$R.squared
# cov1 (exp)   0.148292670      | cov2 (out) 0.009417177     
med_bmg$F.statistics
# cov1  3.8732956        |  cov2 0.2459697     
med_bmg$p.permanova
# cov1   0.00019996     | cov2 0.25960000   
med_bmg$med.p.permanova # 0.2596

med_bmg2 <- permanovaFL(formula = maas_asv_g | confound ~ (treatment.x+LOG_SD) +
                         (BM_grass), 
                       data = pI_env[pI_env$GRASS == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmg2$R.squared
# cov1 (exp)   0.11679152      | cov2 (out) 0.01200791  
med_bmg2$F.statistics
# cov1  5.004837          |  cov2 0.514572    
med_bmg2$p.permanova
# cov1   0.00019996     | cov2 0.11620000  
med_bmg2$med.p.permanova # 0.1162

med_bmg3 <- permanovaFL(formula = maas_asv_g | confound ~ (PO4_P+SWC_P) +
                         (BM_grass), 
                       data = pI_env[pI_env$GRASS == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmg3$R.squared
# cov1 (exp)   0.09290180     | cov2 (out) 0.01002962    
med_bmg3$F.statistics
# cov1  3.7972564       |  cov2 0.4099496    
med_bmg3$p.permanova
# cov1   0.00019996     | cov2 0.25960000  
med_bmg3$med.p.permanova

# herbs
maas_asv_h <- maas_asv[pI_env$HERB == "1",]
med_bmh <- permanovaFL(formula = maas_asv_h | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                         (BM_herb), 
                       data = pI_env[pI_env$HERB == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmh$R.squared
# cov1 (exp)   0.146757940       | cov2 (out) 0.003558743  
med_bmh$F.statistics
# cov1  5.9761379        |  cov2 0.1449158   
med_bmh$p.permanova
# cov1   0.00019996     | cov2 0.69760000     
med_bmh$med.p.permanova #0.6976

med_bmh2 <- permanovaFL(formula = maas_asv_h | confound ~ (treatment.x+LOG_SD) +
                         (BM_herb), 
                       data = pI_env[pI_env$HERB == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmh2$R.squared
# cov1 (exp)   0.106850898       | cov2 (out) 0.004015466  
med_bmh2$F.statistics
# cov1  7.0101600         |  cov2 0.2634424     
med_bmh2$p.permanova
# cov1   0.00019996     | cov2 0.63740000     
med_bmh2$med.p.permanova #0.6374

med_bmh3 <- permanovaFL(formula = maas_asv_h | confound ~ (PO4_P+SWC_P) +
                         (BM_herb), 
                       data = pI_env[pI_env$HERB == "1",], test.mediation = TRUE, 
                       dist.method = "bray", seed = 4321)
med_bmh3$R.squared
# cov1 (exp)   0.067683380       | cov2 (out) 0.003791136  
med_bmh3$F.statistics
# cov1  6.4146193         |  cov2 0.3593008     
med_bmh3$p.permanova
# cov1   0.00019996     | cov2 0.71060000      
med_bmh3$med.p.permanova #0.7106


##### 12.2 Mediation of Root Biomass #####
# all roots
med_rbm <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                        (root.bm), test.mediation = TRUE, 
                      data = pI_env,
                      dist.method = "bray", seed = 4321)
med_rbm$R.squared
# cov1 (exp)   0.137384587      | cov2 (out) 0.006973734  
med_rbm$F.statistics
# cov1  7.000557     |  cov2 0.355353   
med_rbm$p.permanova
# cov1   0.00019996     | cov2 0.05300000  
med_rbm$med.p.permanova # 0.053

med_rbm2 <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD) +
                         (root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_rbm2$R.squared
# cov1 (exp)   0.10255655      | cov2 (out) 0.00744013  
med_rbm2$F.statistics
# cov1  8.4503211     |  cov2 0.6130421   
med_rbm2$p.permanova
# cov1   0.00019996     | cov2 0.04760000   
med_rbm2$med.p.permanova  # 0.0476

med_rbm3 <- permanovaFL(formula = maas_asv | confound ~ (PO4_P+SWC_P) +
                         (root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_rbm3$R.squared
# cov1 (exp)   0.06479026      | cov2 (out)  0.01187117 
med_rbm3$F.statistics
# cov1   7.753737    |  cov2   1.420676 
med_rbm3$p.permanova
# cov1   0.00019996      | cov2   0.00239952 
med_rbm3$med.p.permanova # 0.002
 
# coarse roots
med_cbm <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                         (coarse_root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_cbm$R.squared
# cov1 (expo)    0.137384587  | cov2 (out)   0.002810349 
med_cbm$F.statistics
# cov1    6.9666582     |  cov2   0.1425104
med_cbm$p.permanova
# cov1    0.00019996     | cov2    0.73260000 
med_cbm$med.p.permanova # 0.7326

med_cbm2 <- permanovaFL(formula = maas_asv | confound ~ (treatment.x+LOG_SD) +
                         (coarse_root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_cbm2$R.squared
# cov1 (expo)   0.10255655  | cov2 (out)   0.00323464 
med_cbm2$F.statistics
# cov1   8.4105790      |  cov2   0.2652702
med_cbm2$p.permanova
# cov1   0.00019996     | cov2  0.64220000   
med_cbm2$med.p.permanova # 0.6422

med_cbm3 <- permanovaFL(formula = maas_asv | confound ~ (PO4_P+SWC_P) +
                         (coarse_root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_cbm3$R.squared
# cov1 (expo)   0.064790265  | cov2 (out)   0.003331523
med_cbm3$F.statistics
# cov1    7.6826823   |  cov2   0.3950444 
med_cbm3$p.permanova
# cov1    0.00019996   | cov2    0.65720000
med_cbm3$med.p.permanova # 0.6571

# fine roots
med_frbm <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                         (fine_root.bm), test.mediation = TRUE, 
                       data = pI_env,
                       dist.method = "bray", seed = 4321)
med_frbm$R.squared
# cov1 (exp)   0.13738459    | cov2 (out)   0.01022537 
med_frbm$F.statistics
# cov1   7.0272618      |  cov2    0.5230306
med_frbm$p.permanova
# cov1   0.00019996    | cov2    0.00419916 
med_frbm$med.p.permanova # 0.0042

med_frbm2 <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD) +
                          (fine_root.bm), test.mediation = TRUE, 
                        data = pI_env,
                        dist.method = "bray", seed = 4321)
med_frbm2$R.squared
# cov1 (exp)   0.102556552       | cov2 (out)    0.009953592
med_frbm2$F.statistics
# cov1    8.4742533     |  cov2    0.8224658
med_frbm2$p.permanova
# cov1    0.00019996     | cov2    0.00499900
med_frbm2$med.p.permanova # 0.0049

med_frbm3 <- permanovaFL(formula = maas_asv |  confound ~ (PO4_P+SWC_P) +
                          (fine_root.bm), test.mediation = TRUE, 
                        data = pI_env,
                        dist.method = "bray", seed = 4321)
med_frbm3$R.squared
# cov1 (exp)   0.06479026    | cov2 (out)    0.01751472
med_frbm3$F.statistics
# cov1     7.801420    |  cov2    2.108954
med_frbm3$p.permanova
# cov1     0.00019996     | cov2    0.00019996 
med_frbm3$med.p.permanova # 0.00019

##### 12.3 Mediation of Respiration / Micr. Biomass #####
# microbial biomass
med_mb <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                          (soil_microbial_biomass_C), test.mediation = TRUE, 
                        data = pI_env,
                        dist.method = "bray", seed = 4321)
med_mb$R.squared
# cov1 (exp)    0.136031954       | cov2 (out)   0.007629341 
med_mb$F.statistics
# cov1    7.0213016     |  cov2    0.3937891
med_mb$p.permanova
# cov1    0.00019996       | cov2    0.03580000
med_mb$med.p.permanova # 0.0358

med_mb2 <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD) +
                        (soil_microbial_biomass_C), test.mediation = TRUE, 
                      data = pI_env,
                      dist.method = "bray", seed = 4321)
med_mb2$R.squared
# cov1 (exp)     0.101519820     | cov2 (out)    0.008077648
med_mb2$F.statistics
# cov1       8.475163  |  cov2   0.674345 
med_mb2$p.permanova
# cov1     0.00019996      | cov2    0.02880000
med_mb2$med.p.permanova # 0,0288

med_mb3 <- permanovaFL(formula = maas_asv |  confound ~ (PO4_P+SWC_P) +
                        (soil_microbial_biomass_C), test.mediation = TRUE, 
                      data = pI_env,
                      dist.method = "bray", seed = 4321)
med_mb3$R.squared
# cov1 (exp)    0.06363848      | cov2 (out)    0.01376660
med_mb3$F.statistics
# cov1    7.725503     |  cov2    1.671220
med_mb3$p.permanova
# cov1    0.00019996       | cov2    0.00059988
med_mb3$med.p.permanova # 0.00059988

# microbial respiration
med_mres <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD+PO4_P+SWC_P) +
                        (basal_respiration), test.mediation = TRUE, 
                      data = pI_env,
                      dist.method = "bray", seed = 4321)
med_mres$R.squared
# cov1 (exp)    0.136031954   |  cov2 (out)    0.005966609
med_mres$F.statistics
# cov1    7.0076950     |  cov2    0.3073703
med_mres$p.permanova
# cov1     0.00019996      | cov2    0.11140000
med_mres$med.p.permanova # 0.1114

med_mres2 <- permanovaFL(formula = maas_asv |  confound ~ (treatment.x+LOG_SD) +
                          (basal_respiration), test.mediation = TRUE, 
                        data = pI_env,
                        dist.method = "bray", seed = 4321)
med_mres2$R.squared
# cov1 (exp)     0.1015198      | cov2 (out)   0.0166281 
med_mres2$F.statistics
# cov1     8.557338   |  cov2    1.401620
med_mres2$p.permanova
# cov1     0.00019996      | cov2    0.00019996
med_mres2$med.p.permanova # 0.000199996

med_mres3 <- permanovaFL(formula = maas_asv |  confound ~ (PO4_P+SWC_P) +
                          (basal_respiration), test.mediation = TRUE, 
                        data = pI_env,
                        dist.method = "bray", seed = 4321)
med_mres3$R.squared
# cov1 (exp)    0.063638481      | cov2 (out)   0.008604342 
med_mres3$F.statistics
# cov1     7.682517    |  cov2   1.038727 
med_mres3$p.permanova
# cov1     0.00019996      | cov2   0.02100000 
med_mres3$med.p.permanova # 0.021

