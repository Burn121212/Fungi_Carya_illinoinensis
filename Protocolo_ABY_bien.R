# Análisis ecologógicos bioestadísticos Micobioma de Carya ilionensis
## Cargar librerias
```
library(ampvis2)         #para análisis de prevalencia y PCAs
library(BiocManager)     #para instalar phyloseq
library(ape)             #para exportar arboles en nwk
library(ggplot2)         #para barplots
library(dplyr)           #para filtrar
library(phyloseq)        #para análisis de composicion y estructura
library(ranacapa)        #para ver curvas rarefacción
library(readxl)          #para exportar en csv
library(tidyverse)       #analisis de diversidad
library(vegan)           #analisis de diversidad
library(corrplot)        #analisis de correlación 
```
---
## Crear archivo Phyloseq
```	
#Importar archivos 
otu_mat = read.delim("ABY_table.txt",row.names=1)
tax_mat = read.delim("ABY_clasification.txt",row.names=1)
pls_mat = read.delim("ABY_clasification_pls.txt",row.names=1)
samples_df = read.delim("ABY_metadata.txt",row.names=1)
#Transformar en matrices otu y tax tables
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
pls_mat <- as.matrix(pls_mat)
#Transformar a objetos phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
PLS = tax_table(pls_mat)
samples = sample_data(samples_df)
#root_tree <- read_tree("C:/Users/ST/Desktop/Ber/Postdoctorado/tree.nwk")
OTS <- phyloseq (OTU, TAX, samples)
OPS <- phyloseq (OTU, PLS, samples)
#OTS <- prune_samples(sample_names(OTS) != "Mock", OTS) # Remover potenciales muestras sint?ticas
OTS
#Conteo de tax
tax_data <- tax_table(OTS) %>% as.data.frame()
# Calcular el número de elementos únicos por nivel taxonómico
taxonomy_summary <- sapply(c("Phylum", "Class", "Order", "Family", "Genus","Species"), function(rank) {
	length(unique(na.omit(tax_data[[rank]])))
})
taxonomy_summary
```
## Rarefacción
```
#Checar sample depth rare
summary(sample_sums(OTS))
hist(sample_sums(OTS))
#Hacer rarefaccion
set.seed(123)
OTS_rar <- rarefy_even_depth(OTS, sample.size = 10000, rngseed = 123, replace = FALSE, verbose = TRUE)
OTS_rar
hist(sample_sums(OTS_rar))
# Conteo tax
tax_data_rar <- tax_table(OTS_rar) %>% as.data.frame()
# Calcular el número de elementos únicos por nivel taxonómico
taxonomy_summary_rar <- sapply(c("Phylum", "Class", "Order", "Family","Genus","Species"), function(rank) {
	length(unique(na.omit(tax_data_rar[[rank]])))
})
taxonomy_summary_rar
```
---
## Composición con Barplots
```
#Definir colores

phylum_colors <-c("purple","#FFAA92","#009E73","green", "#5E738F", "red4","blue","#DA5724","#F0E442","#619CFF","maroon4","#575329", "#00FECF", "#B05B6F",
                          "#8CD0FF", "darkolivegreen1", "#999999", "#C8A1A1","#1E6E00", "#320033", "#66E1D3", "#CFCDAC", "#C84248", "#4A3B53", "#FF2F80",
                           "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#575329", "#252A52", "#B05B6F", "#8CD0FF", 
                          "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                          "#3B5DFF", "#549E79","#7B4F4B", "#A1C299", "purple", "#999999", "#E69F00",  "#009E73","darkorange1","darkgrey","#FC4E07","#D14285",  
                          "#652926", "red4" , "#CC79A7",  "#009E73","#00BA38" , "#252A52" , "brown",    "#D55E00" , "cyan1", "royalblue4","#CBD588", "#5F7FC7", 
                          "#CC79A7","orange","#56B4E9","#DA5724", "#508578","#C84248", "darkorchid", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                           "darkseagreen", "#5E738F", "#D1A33D","#8A7C64", "#599861", "yellow","darkgoldenrod1", "#56B4E9","darkolivegreen1","#F0E442","#0072B2","#D55E00")
#Aglomerar Phylum
xP_01 <- OTS_rar%>%
  tax_glom(taxrank = "Phylum") %>%                                  # agglomerate 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.01) %>%                                      # Filter out low abundance taxa
  arrange(Phylum)     
#Graficar Phylum
zP_01<- ggplot(xP_01, aes(x = Name, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity",position = "Fill") + #Fill, Stack,
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 1%) \n") +
  ggtitle("Phylum_abundance")
zP_01 <- zP_01+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zP_01grid <- zP_01 + facet_grid(~Location, scales = "free", space='free')
zP_01grid
#Aglomerar Class
xC_02 <- OTS_rar%>%
	tax_glom(taxrank = "Class") %>%                                  # agglomerate at family level
	transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
	psmelt() %>%                                                      # Melt to long format
	filter(Abundance > 0.02) %>%                                      # Filter out low abundance taxa
	arrange(Class)
#Grafica Class
zC_05<- ggplot(xC_02, aes(x = Name, y = Abundance, fill = Class)) +
	geom_bar(stat = "identity",position = "Fill") + 
	scale_fill_manual(values = phylum_colors) +
	theme(axis.title.x = element_blank()) +
	guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
	ylab("Relative Abundance (Phyla > 2%) \n") +
	ggtitle("Class_abundance")
zC_05 <- zC_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zC_05grid <- zC_05 + facet_grid(~Location, scales = "free", space='free')
zC_05grid
#Con Order
xO_02 <- OTS_rar%>%
  tax_glom(taxrank = "Order") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.02) %>%                                      # Filter out low abundance taxa
  arrange(Order)
#
zO_05<- ggplot(xO_02, aes(x = Name, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Order_abundance")
zO_05 <- zO_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zO_05grid <- zO_05 + facet_grid(~Location, scales = "free", space='free')
zO_05grid
#Aglomerar Family
xF_05 <- OTS_rar%>%
  tax_glom(taxrank = "Family") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.03) %>%                                      # Filter out low abundance taxa
  arrange(Family)
#Graficar Family
zF_05<- ggplot(xF_05, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 5%) \n") +
  ggtitle("Family_abundance")
zF_05 <- zF_05+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zF_05grid <- zF_05 + facet_grid(~Location, scales = "free", space='free')
zF_05grid
#Aglomerar Genus 
xG_04 <- OTS_rar%>%
  tax_glom(taxrank = "Genus") %>%                                  # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%              # Transform to rel. abundance
  psmelt() %>%                                                      # Melt to long format
  filter(Abundance > 0.04) %>%                                      # Filter out low abundance taxa
  arrange(Genus)
#Graficar Genus
zG_04<- ggplot(xG_04, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 4%) \n") +
  ggtitle("Genus_abundance")
zG_04 <- zG_04+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
zG_04grid <- zG_04 + facet_grid(~Location, scales = "free", space='free')
zG_04grid
```
---
## Exportar KRONA2 sin corte de abundancia relativa
```
glom_o <- tax_glom(OTS, taxrank = "Order")
write.csv(glom_o@tax_table, file = "Order_tax.csv")
write.csv(glom_o@otu_table, file = "Order_table.csv")
#
glom_f <- tax_glom(OTS, taxrank = "Family")
write.csv(glom_f@tax_table, file = "Family_tax.csv")
write.csv(glom_f@otu_table, file = "Family_table.csv")
#
glom_g <- tax_glom(OTS, taxrank = "Genus")
write.csv(glom_g@tax_table, file = "Genus_tax.csv")
write.csv(glom_g@otu_table, file = "Genus_table.csv")
```
---
## Diversidad Beta
```
d_bray = distance(OTS_rar, method='bray')
plot(hclust(d_bray, method="ward"), xlab="bray_curtis")
d_jacc = distance(OTS_rar, method='jaccard')
plot(hclust(d_jacc, method="ward"), xlab="jaccard")
hc <- hclust(d, method="ward")
my_tree <- as.phylo(hc) 
write.tree(phy=my_tree, file="bray_tree.newick")
heatmap(as.matrix(d), xlab="bray_curtis")
bet_bray = distance(OTS_rar, method='bray')                     #Bray curtis toma en cuenta la abundancia absoluta
bet_bray
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
```
---	
## Diversidad Alfa
```
#Otu_table debe ser numeros enteros no decimales
#Graficar varios indices
alpha_meas = c("Observed","Chao1", "Shannon", "Simpson")
alf_varios <- (p <- plot_richness(OTS_rar, "Name", color = "Location", measures=alpha_meas))
alf_varios
alf_varios_grid <- alf_varios + facet_grid(~Location, scales = "free", space='free')alf_varios
alf_varios_grid
#Graficar solo riqueza "observed"
alpha_meas = c("Observed")
alf_obs <- (p <- plot_richness(OTS, "Name", color = "Location", measures=alpha_meas))
alf_obs_grid <- alf_obs + facet_grid(~Location, scales = "free", space='free')
alf_obs_grid
alf_obs_grid <- alf_obs + facet_grid(~Location, scales = "free", space='free') + 
  ggtitle("Richness")
alf_obs_grid
#
#Generar tabla de riqueza con todos los indices
rich = estimate_richness(OTS_rar)
rich
write.csv(rich, file = "alpha_diversity_b.csv")
```
---
## RDA analysis ajustado con Hellinger 
```
otu_table <- read.csv("otu_table.csv", row.names = 1) #la tabla de otus debe estar en forma de transpocicion, filas son sitios
env_data  <- read.csv("env_data.csv", row.names = 1)
# Asegúrate que las muestras estén en el mismo orden
otu_table <- otu_table[rownames(env_data), ]
#Tranformacion de Hellinger
otu_hell <- decostand(otu_table, method = "hellinger")
#Relizar RDA 
rda_model <- rda(otu_hell ~ COS + MOS + pH + Zn + CT + Humidity + Temp + CE + CaCo3, data = env_data)
#Ver resumen estadistico
summary(rda_model)               # Resumen del modelo
anova(rda_model)                 # Prueba de significancia global (permutación)
anova(rda_model, by = "axis")    # Significancia por eje
anova(rda_model, by = "terms")   # Significancia por variable ambiental
#Visualizacion
plot(rda_model, scaling = 2, main = "RDA with Hellinger-transformed data")
```
---
## Permanova entre dos grupos 
```
# Asegurar que Zone sea factor (centro y periferia)
env_data$Zone <- as.factor(env_data$Zone)
# Crear matriz de distancias
dist_matrix <- vegdist(otu_hell, method = "euclidean")
# Ejecutar PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Zone, data = env_data, permutations = 999)
# Mostrar resultados
print(permanova_result)
# Exportar tabla
permanova_df <- as.data.frame(permanova_result) # Convertir los resultados a data.frame
write.csv(permanova_df, "permanova_result.csv")
```
---
## AMPVIS2 Prevalencia y PCA
```
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
#PCA
amp_ordinate(
  data = av2_obj,
  type = "pca",
  sample_color_by = "Location",
  sample_colorframe = FALSE,
  sample_colorframe_label = "Name",
  species_plotly = TRUE)    
#Preavalencia Heatmap
amp_heatmap(av2_obj,
            group_by = "Name",
            facet_by = "Location",
            tax_aggregate = "Species",
            #tax_add = "Family",
            tax_show = 30,
            color_vector = c("white", "cyan"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="left")
#
```
---
## NMDS
```
#crear una table con top30
#importar tabla y metadatos promediados
tab = read.delim("C:/Users/ST/Desktop/Ber/Proyecto_Abby/ABY/ABY_TOP25_Table.txt",row.names=1)
env = read.delim("C:/Users/ST/Desktop/Ber/Proyecto_Abby/ABY/ABY_metadata.txt",row.names=1)
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
```
---
## CAP con vectores
```
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
  geom_point(aes(colour = Location),  size = 3) + 
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
```
---
## ANOVA fisicoquímicos
```
anova(cap_ord)
#Analisis correlacion
#setwd
setwd("D:/Postdoc/Otros/Abby/R_analysis")
data<- read_excel("C:/Users/ST/Desktop/Ber/Proyecto Abby/ABY/aby_richnes_correlation.xlsx")
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
rcx= rcorr(Index, type = c("spearman"))
df.rcx.r=data.frame(rcx$r)
df.rcx.p=data.frame(rcx$P)
#write csv
write.csv(corr,'correlatio_spearman.csv')
write.csv(df.rcx.p,'correlationmatrix_index_p.csv')
```
---
## Redes de coocurrecnia
```
p <- plot_net(OTS, distance ="bray", type= "taxa", maxdist=0.4, color= "NULL", shape= "NULL")
p
#

plot_net(OTS, type = "taxa", point_label = "Genus", point_size = 10, point_alpha = 0.5, maxdist = 0.5, color = "Phylum", distance = "bray", laymeth = "auto") 
co_occurrence_network(xP_01, treatment = "Location", 
                      classification = 'Phylum')
plot_net(xP_01, treatment = "Location", 
                      classification = 'Phylum')
# otro tipo de plot
plot_net(BC, maxdist=0.4, color="Site", shape="Biome")
#
p1 <- ggviolin(OTS, x = "species", y = "Shannon",
               add = "boxplot", fill = "species", palette = c("#a6cee3", "#b2df8a", "#fdbf6f"))  
print(p1)
```
---
## Diagrama de Venn
```
amp_venn(av2_obj, group_by = "Location", cut_a = 0, cut_f = 50, text_size = 3)
```





