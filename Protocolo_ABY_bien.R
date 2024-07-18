library(ampvis2)         #para análisis de prevalencia y PCAs
library(BiocManager)
library(ape)             #para exportar arboles en nwk
library(ggplot2)         #para barplots
library(dplyr)
library(phyloseq)        #para análisis de composicion y estructura
library(ranacapa)
library(readxl)          #para exportar en csv
library(tidyverse)
library(vegan)

# CREAR OBJETO PHYLOSEQ-------------------------------------------------------------------------------------------------------------

#Importar archivos 
otu_mat = read.delim("D:/Postdoc/Otros/Abby/R_analysis/R_chido/ABY_table.txt",row.names=1)
tax_mat = read.delim("D:/Postdoc/Otros/Abby/R_analysis/R_chido/ABY_clasification.txt",row.names=1)
#tax_mat = read.delim("D:/Postdoc/Otros/Abby/R_analysis/R_chido/ABY_clasification_pls.txt",row.names=1)
samples_df = read.delim("D:/Postdoc/Otros/Abby/R_analysis/R_chido/ABY_metadata.txt",row.names=1)

#Transformar en matrices otu y tax tables
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transformar a objetos phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
#root_tree <- read_tree("C:/Users/ST/Desktop/Ber/Postdoctorado/tree.nwk")
OTS <- phyloseq (OTU, TAX, samples)
OTS <- prune_samples(sample_names(OTS) != "Mock", OTS) # Remover potenciales muestras sint?ticas
OTS

#EXPLORACION DE COMPOSICION BARPLOTS-------------------------------------------------------------------------------------------------------------------------------------------------

#Definir colores
phylum_colors <-c("burlywood3","#00FECF","chartreuse","blue3","#5E738F", "red4","996600","gold","pink","orange","#009E73","#619CFF","chocolate3","purple","red","black", "#B05B6F",
                          "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1","#1E6E00","99FFCC", "#320033", "#66E1D3", "#CFCDAC", "#4FC601", "#4A3B53", "#FF2F80",
                          "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", 
                          "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                          "#3B5DFF", "#549E79","#7B4F4B", "#A1C299", "purple", "#999999", "#E69F00",  "#009E73","darkorange1","darkgrey","#FC4E07","#D14285",  
                          "#652926", "red4" , "#CC79A7",  "#009E73","#00BA38" , "#252A52" , "brown",    "#D55E00" , "cyan1", "royalblue4","#CBD588", "#5F7FC7", 
                          "#CC79A7","orange","#56B4E9","#DA5724", "#508578","#C84248", "darkorchid", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                          "darkseagreen", "#5E738F","#D1A33D","#8A7C64", "#599861", "yellow","darkgoldenrod1", "#56B4E9","darkolivegreen1","#F0E442", "#0072B2", "#D55E00")

#Con PRIMARY LIFESTYLE--------------------------------------------------------------
#Cortar datos menores a 0.5%
xP_01 <- OTS%>%
	tax_glom(taxrank = "primary_lifestyle") %>%                                  # agglomerate at family level
	transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
	psmelt() %>%                                                      # Melt to long format
	filter(Abundance > 0.0005) %>%                                      # Filter out low abundance taxa
	arrange(primary_lifestyle)     
#Graficar
zP_01<- ggplot(xP_01, aes(x = Name, y = Abundance, fill = primary_lifestyle)) +
	geom_bar(stat = "identity",position = "Fill") + #Fill, Stack,
	scale_fill_manual(values = phylum_colors) +
	theme(axis.title.x = element_blank()) +
	guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
	ylab("Relative Abundance (Phyla > 0.005%) \n") +
	ggtitle("primary_lifestyle_abundance")
zP_01 <- zP_01+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zP_01
zP_01grid <- zP_01 + facet_grid(~Location, scales = "free", space='free')
zP_01grid

#Con Phylum--------------------------------------------------------------
#Cortar datos menores a 1%
xP_01 <- OTS%>%
  tax_glom(taxrank = "Phylum") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.01) %>%                                      # Filter out low abundance taxa
  arrange(Phylum)     
#Graficar
zP_01<- ggplot(xP_01, aes(x = Name, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",position = "Fill") + #Fill, Stack,
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 1%) \n") +
  ggtitle("Phylum_abundance")
zP_01 <- zP_01+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zP_01
zP_01grid <- zP_01 + facet_grid(~Location, scales = "free", space='free')
zP_01grid
#Con Order ---------------------------------------------------------------
#Cortar datos menores a 5%
xO_05 <- OTS%>%
  tax_glom(taxrank = "Order") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.03) %>%                                      # Filter out low abundance taxa
  arrange(Order)
#
zO_05<- ggplot(xO_05, aes(x = Name, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 3%) \n") +
  ggtitle("Order_abundance")
zO_05 <- zO_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zO_05
zO_05grid <- zO_05 + facet_grid(~Location, scales = "free", space='free')
zO_05grid

#Con Family ---------------------------------------------------------------
#Cortar datos menores a 5%
xF_05 <- OTS%>%
  tax_glom(taxrank = "Family") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.03) %>%                                      # Filter out low abundance taxa
  arrange(Family)
#
zF_05<- ggplot(xF_05, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 3%) \n") +
  ggtitle("Family_abundance")
zF_05 <- zF_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zF_05
zF_05grid <- zF_05 + facet_grid(~Location, scales = "free", space='free')
zF_05grid


#Cortar datos menores a 5%
xF_05 <- OTS%>%
	tax_glom(taxrank = "Genus") %>%                                  # agglomerate at family level
	transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
	psmelt() %>%                                                      # Melt to long format
	filter(Abundance > 0.03) %>%                                      # Filter out low abundance taxa
	arrange(Genus)
#
zF_05<- ggplot(xF_05, aes(x = Sample, y = Abundance, fill = Genus)) +
	geom_bar(stat = "identity",position = "Fill") + 
	scale_fill_manual(values = phylum_colors) +
	theme(axis.title.x = element_blank()) +
	guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
	ylab("Relative Abundance (Phyla > 3%) \n") +
	ggtitle("Genus_abundance")
zF_05 <- zF_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zF_05
zF_05grid <- zF_05 + facet_grid(~Location, scales = "free", space='free')
zF_05grid


#EXPORTAR A KRONA2 y STAMP sin corte de %-------------------------------------------------------------------

glom_o <- tax_glom(OTS, taxrank = "Order")
write.csv(glom_o@tax_table, file = "C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/Order_tax.csv")
write.csv(glom_o@otu_table, file = "C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/Order_table.csv")
#
glom_f <- tax_glom(OTS, taxrank = "Family")
write.csv(glom_f@tax_table, file = "C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/Family_tax.csv")
write.csv(glom_f@otu_table, file = "C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/Family_table.csv")
#

#ALFA DIVERSIDAD--------------------------------------------------------------------------------------------------------------------

#Otu_table debe ser numeros enteros no decimales
#Graficar varios indices
alpha_meas = c("Observed","Chao1", "Shannon", "Simpson")
alf_varios <- (p <- plot_richness(OTS, "Name", color = "Location", measures=alpha_meas))
alf_varios
#Graficar solo riqueza "observed"
alpha_meas = c("Observed")
alf_obs <- (p <- plot_richness(OTS, "Name", color = "Location", measures=alpha_meas))
alf_obs
alf_obs_grid <- alf_obs + facet_grid(~Location, scales = "free", space='free') + 
  ggtitle("Richness")
alf_obs_grid

#Graficar solo riqueza "observed"
alpha_meas = c("Shannon")
alf_obs <- (p <- plot_richness(OTS, "Name", color = "Location", measures=alpha_meas))
alf_obs
alf_obs_grid <- alf_obs + facet_grid(~Location, scales = "free", space='free') + 
	ggtitle("Richness")
alf_obs_grid
#
#Generar tabla de riqueza con todos los indices
rich = estimate_richness(OTS)
rich
write.csv(rich, file = "D:/Postdoc/Otros/Abby/R_analysis/R_chido/alpha_diversity.csv")

#BETA DIVERSIDAD-------------------------------------------------------------
 
bet_bray = distance(OTS, method='bray')                          #Bray curtis toma en cuenta la abundancia absoluta
tree1= plot(hclust(bet_bray, method="ward"), xlab="bray_curtis")

bet_jacc = distance(OTS, method='jaccard')                      #Jaccard toma en cuenta ausencia preserencia
tree1= plot(hclust(bet_jacc, method="ward"), xlab="jaccard")

#graficar beta_heatmap

heatmap(as.matrix(bet_bray), xlab="bray_curtis")
heatmap(as.matrix(bet_jacc), xlab="jaccard")

#exportar dendograma a formato nwk

hc1 = hclust(bet_bray, method="ward")
hc2 = hclust(bet_jacc, method="ward")
class(hc1)
class(hc2)
my_tree1 <- as.phylo(hc1) 
my_tree2 <- as.phylo(hc2) 
write.tree(phy=my_tree1, file="tree_bray.newick")
write.tree(phy=my_tree2, file="tree_jacc.newick")

#tree <- read.tree ("tree_jacc.newick")
#p <- ggtree(tree)



#AMPVIS2 Prevalencia y PCA--------------------------------------------------------------------------------------------------------------

#Convertir de phyloseq a ampvis
Totu_table =t(otu_table(OTS))
otu_table(OTS)=Totu_table
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(OTS)@.Data)),
                           t(phyloseq::otu_table(OTS)@.Data),
                           phyloseq::tax_table(OTS)@.Data,
                           check.names = F)
#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(OTS), 
                           check.names = F)
av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)
#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)
amp_rankabundance(av2_obj, group_by = "Name")
sqrt
log10

#PCA--------------------------------------------------------------------------------
amp_ordinate(
  data = av2_obj,
  type = "pca",
  sample_color_by = "Location",
  sample_colorframe = FALSE,
  sample_colorframe_label = "Name",
  species_plotly = TRUE)    

#Preavalencia Heatmap--------------------------------------------------------------------------
amp_heatmap(av2_obj,
            group_by = "Name",
            facet_by = "Location",
            tax_aggregate = "Phylum",
            #tax_add = "Family",
            tax_show = 13,
            color_vector = c("white", "blue"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="left")
#

#NMDS----------------------------------------------------------------------------------------------
#crear una table con top30
#importar tabla y metadatos promediados
tab = read.delim("C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/ABY_TOP30_table.txt",row.names=1)
env = read.delim("C:/Users/ST/Desktop/Ber/Postdoctorado/ABY/ABY_metadata.txt",row.names=1)
tab
env
#definir columnas tabla y metadata
com = tab[,1:8]
com
envi = env[,3:11]
envi
#convert com to a matrix
m_com = as.matrix(com)
##Perform the NMDS ordination,nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
#Graficar
plot(nmds, type = "t") 
plot(envi)


#CAP CON VECTORES FISICOQUIMICOS---------------------------------------------------------------------------------------------------------------------------------

erie <- OTS
erie
bray_not_na <- phyloseq::distance(physeq = erie, method = "bray")

#CAP ordinate
cap_ord <- ordinate(
  physeq = erie, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~  COS+MOS+pH +Zn+CT+Humidity+Temp+CE+CaCo3)

# CAP plot
cap_plot <- plot_ordination(
  physeq = erie, 
  ordination = cap_ord, 
  color = "Bray_Group", 
  axes = c(1,2)
) + 
  #aes(shape = biome) + 
  geom_text(aes(label = Name),hjust = 1, 
            nudge_y = 0.1, colour= "black") + # label es la etiqueta que quieres que aparezca en el punto
  geom_point(aes(colour = Name),  size = 3) + 
  geom_point(colour = "white", size = 1) + 
  scale_color_manual(values = c( 
    "cyan2",  "blue3", "yellow", "red", "purple", "green", "orange", "pink",
    "#70dc85", "#78ac83", "#006602","pink", "gold","#8A7C64", "#599861", "navy", "#5F7FC7" , "tomato", "#673770",  
    "#008080", "#2F4F4F", "#FAEBD7", "#ff1493", "#5e738f","#808000", "#D14285", "#ffa500", "cbd588", "wheat", 
    "#d2b48c", "cyan2","black",  "#BC8F8F", "#800000","#008B8B",  "#BC8F8F", "red", "BF5650",
    "#0089A3","#66796D","orange","#3A2465" , "green","purple","black","deeppink" , "blue",
    "cyan1", "gold", "coral")) 

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.05* CAP1, 
                 y = 1.15 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
p <- cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = 0.5, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 5,
    color= "black",
    data = arrowdf, 
    show.legend = TRUE,
  ) + theme(text=element_text(size=5, family="Times New Roman"))+ theme_light()

p


#plot_redes de coocurrecnia----------------------------------------------------------------------------------------------------------
#
p <- plot_newet(OTS, maxdist=0.4, color="Name", shape="Location")

plot_net(OTS, type = "taxa", point_label = "Genus", point_size = 10, point_alpha = 0.5, maxdist = 0.5, color = "Phylum", distance = "bray", laymeth = "auto") 

#
p1 <- ggviolin(OTS, x = "species", y = "Shannon",
               add = "boxplot", fill = "species", palette = c("#a6cee3", "#b2df8a", "#fdbf6f"))  
print(p1)

#
amp_venn(av2_obj, group_by = "Location", cut_a = 0, cut_f = 50, text_size = 3)


glom_o <- tax_glom(MX55, taxrank = "order")
glom_o

glom_c <- tax_glom(MX55, taxrank = "class")
glom_c

glom_f <- tax_glom(MX55, taxrank = "family")
glom_f

glom_g <- tax_glom(MX55, taxrank = "genus")
glom_g

#Script_correlation
#cor script
library(Hmisc)
library(corrplot)
library(corrtable)
library(readxl)
install.packages("corrtable")

#setwd
setwd("D:/Postdoc/Otros/Abby/R_analysis")

data<- read_excel("D:/Postdoc/Otros/Abby/R_analysis/Abby_fisicoquimicos.xlsx")

#create data frame
names(data)
Index<-data.frame(data$COS,data$MOS,data$pH,data$Zn,data$CT,data$Humidity,data$Temp,data$CE,data$CaCo3,data$Richness)



corr<-cor(Index, method = c("spearman"))
corr

#plot
corrplot(corr, type = "upper", order = "original",
									tl.col = "black", tl.srt = 45)

#P values

Index<-as.matrix(Index)
rcx=rcorr(Index, type = c("spearman"))

df.rcx.r=data.frame(rcx$r)
df.rcx.p=data.frame(rcx$P)

#write csv

write.csv(df.rcx.r,'correlationmatrix_index.csv')
write.csv(df.rcx.p,'correlationmatrix_index_p.csv')



###########################################################################################################
cap_ord <- ordinate(
	physeq = erie, 
	method = "CAP",
	distance = bray_not_na,
	formula = ~  COS+MOS+pH +Zn+CT+Humidity+Temp+CE+CaCo3)



#Analisis p value para fisicoquimicos (ANOVA)#############################################################################################################################
anova(cap_ord)


# Calculate bray curtis distance matrix
MCA_bray <- phyloseq::distance(erie, method = "bray")
# make a data frame from the sample_data

sampledf <- data.frame(sample_data(erie))

# Adonis test
ad_table <- adonis2(MCA_bray ~ COS+MOS+pH +Zn+CT+Humidity+Temp+CE+CaCo3, data = sampledf) 

write.csv(ad_table, "table_CAP.csv")



