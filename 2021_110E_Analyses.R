
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# load packages ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library("mvabund")
library("knitr")
library("BiocStyle")
library("phyloseq")
library("gridExtra")
library("ggplot2")
library("dada2")
library("DECIPHER")
library("phangorn")
library("dplyr")
library("data.table")
library("BBmisc")
library("vegan")
library("usdm")
library(reshape2)
library(RColorBrewer)
library(dendextend)
# library(metagMisc)

load(file = "dada2_110E.RData")
metadata_mat <- read.csv("data/metadata_110.csv", sep = ",", colClasses=c("ID"="character"))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Explore metadata ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# determine multicolinearity by correlation
# use PCA to generate new
# determine multicolinearity on new table by correlation
# AIC will be used to refine this in new models

colnames(metadata_mat)
summary(metadata_mat)
meta <- normalize(metadata_mat[ , c(17:26)], method="standardize", range =c(0,1), margin=2L)
meta$latitude <- metadata_mat$Latitude
summary(meta)
pairs(meta,
      panel = panel.smooth,
      main = "Bivariate plots of abiotic variables"
)

summary(meta)

# explore by PCA
pc.met <- prcomp(meta[1:10], center = FALSE, scale. = FALSE) 
pc.met_sub <- prcomp(meta[ -c(1,2,3), -c(9,10,11)] , center = FALSE, scale. = FALSE)

library(ggfortify)
autoplot(pc.met, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 2, label = TRUE)
autoplot(pc.met_sub, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 2, label = TRUE)
autoplot(pc.met,  
         x = 1,    # PC3
         y = 3,    # PC4
         loadings = TRUE, loadings.label = TRUE,loadings.label.size = 2, label = TRUE)
autoplot(pc.met_sub, 
         x = 1,    # PC3
         y = 3,    # PC4
         loadings = TRUE, loadings.label = TRUE,loadings.label.size = 2, label = TRUE)

plot(meta$NOx~meta$latitude)
plot(meta$Phosphate_.µM.~meta$latitude)
plot(meta$Ammonia_.uM.~meta$latitude)
plot(meta$Nitrite_.µM.~meta$latitude)

plot(meta$SALINITY_.PSU.~meta$latitude)
plot(meta$Silicate_.µM.~meta$latitude)
plot(meta$OXYGEN_SATURATION_...~meta$latitude)

plot(meta$Oxygen_.µM.~meta$latitude)
plot(meta$SEA_SURFACE_TEMP_.C..~meta$latitude)

# Apparent collinearity
# Pearson r linear correlation among environmental variables 
met.pearson <- cor(meta) # default method = "pearson" 
pearson.metadata <- round(met.pearson, 2)
pearson.metadata
write.csv(pearson.metadata, file = "pearson.metadata.csv")
# # THERE IS:
# apparent multicolinearity between variables
# oxygen saturation, which is a function of temperature and oxygen concentration, is over 100%. This indicates that either the water is rapidly warming or photosynthesis (correlation with fluorescence >0.7)
# oxygen and temperature highly correlated/indistinguishable (-.99)
# fluorescence and nutrients are also highly correlated

pc.nutrients <-prcomp(meta[ , c(6,7,9,10)], center = FALSE, scale. = FALSE) # PC1 98.43%
pc.temp_oxyn <-prcomp(meta[ , c(4,5)], center = FALSE, scale. = FALSE) # PC1 99.58%

cor(cbind.data.frame(meta$SALINITY_.PSU., meta$Silicate_.µM., meta$FLUORESCENCE_.Units., pc.nutrients$x[,1], meta$OXYGEN_SATURATION_..., pc.temp_oxyn$x[,1]))
# Fluorescence is too problematic, correlated with too many variables at this sampling time
cor(cbind.data.frame(meta$SALINITY_.PSU., meta$Silicate_.µM., pc.nutrients$x[,1], meta$OXYGEN_SATURATION_..., pc.temp_oxyn$x[,1]))

meta_min <- cbind.data.frame(pc.nutrients$x[,1], pc.temp_oxyn$x[,1], meta$Silicate_.µM., meta$SALINITY_.PSU., meta$OXYGEN_SATURATION_... )
cor(meta_min)
colnames(meta_min) <- c("nutrient", "temp_oxyg", "silicate", "salinity", "oxygen_sat")
meta_min$latitude <- metadata_mat$Latitude
meta_min$sampleID <- rownames(meta_min)

d <- melt(meta_min[1:6], id.vars="latitude")

ggplot(d, aes(latitude,value)) + 
  geom_point() + 
  stat_smooth() +
  facet_wrap(~variable) 
rm(d)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Make phyloseq object ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

seqtabNoChim <- as_tibble(seqtabNoC)
seqtabNoChim <- seqtabNoC[match(rownames(meta), rownames(seqtabNoC)),]
seqtabNoChim <- as_tibble(seqtabNoChim)

# Create phyloseq object:
ps <- phyloseq(otu_table(seqtabNoChim, taxa_are_rows=FALSE), 
               sample_data(meta_min))
# view elements
ps

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Preprocessing data frame ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Remove low abundance ASVs

prevdf = apply(X = otu_table(ps),#vector with prevalence
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
mult = apply(X = otu_table(ps),#vector with product of all non-0 abundances
             MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
             FUN = function(x){prod(x[which(x > 0)])})
prevdf = data.frame(Prevalence = prevdf,#data frame with above information
                    TotalAbundance = taxa_sums(ps),
                    Product = mult)
#Set abundance/prevalence threshold for filter:
prev_abundance_Threshold = 100 # abundance threshold, which is a product

#Execute filter:
keepTaxa = rownames(prevdf)[(prevdf$Product > prev_abundance_Threshold)]
prev_ps = prune_taxa(keepTaxa, ps)
dim(otu_table(ps))
# how many reads is this?
sum(otu_table(prev_ps))
# how many reads trimmed?
sum(otu_table(ps)) - sum(otu_table(prev_ps))

meta_min$count <- rowSums(otu_table(prev_ps)) # make an offset value for read depth
write.csv(meta_min, "meta_min.csv")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# make phylogenetic tree ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

seqs <- getSequences(otu_table(prev_ps)@.Data)
names(seqs) <- seqs # names variant  as its sequence
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
writeXStringSet(alignment, file = "DECIPHER_alignment.fasta")
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
# best molecular evolutionary model
(mT <- modelTest(phangAlign, treeNJ)) # find model with lowest BIC for the alignment object
# choose best model from the table according to BIC
bestmodel <- mT$Model[which.min(mT$BIC)]
bestmodel
# search for a better tree using NNI rearrangements
fit = pml(treeNJ, data=phangAlign, model=bestmodel)
fitGTR <- update(fit, k=4, inv=0.2) # use update to change parameters
# search for a better tree using NNI rearrangements
fitter <- optim.pml(fitGTR, optNni=TRUE)
fitGTR <- optim.pml(fitter, model="GTR", multicore=TRUE, optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# Bootstrap
bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, control = pml.control(trace = 0))
BA <- plotBS(midpoint(fitGTR$tree), bs, p = 50, type="n") # option type=="n" just assigns the bootstrap values and return the tree without plotting it.

#Merge phylogenetic tree to phyloseq object:
prev_ps_t<-merge_phyloseq(prev_ps,BA)
prev_ps_t<-merge_phyloseq(prev_ps_t,alignment)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Do taxonomy ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

table_prev_ps_t <- as.matrix(otu_table(prev_ps_t))
pr2_database <- "edited_ps2.fasta"
PR2_tax_levels <- c("Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species")

set.seed(123)
pr_taxa_out <- assignTaxonomy(table_prev_ps_t,
                              refFasta=pr2_database,
                              taxLevels = PR2_tax_levels,
                              minBoot = 50, # Changed from zero <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,
                              outputBootstraps = TRUE,
                              multithread=TRUE,
                              verbose = TRUE)
View(pr_taxa_out$tax)
View(pr_taxa_out$boot)

as_taxa_out <- as.data.frame(pr_taxa_out$tax)

# Change all columns to characters (otherwise everything becomes NA
for (i in 1:8){ as_taxa_out[,i] <- as.character(as_taxa_out[,i])}
# Change NA to and empty string (changing the script to use is.na() is also an option
as_taxa_out[is.na(as_taxa_out)] <- ""
# Fill missing taxonomy

for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,2] == ""){
    x <- paste(as_taxa_out[i,1], "_unknown_supergroup", sep = "")
    as_taxa_out[i, 2] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,3] == ""){
    x <- paste(as_taxa_out[i,2], "_unknown_division/phyla", sep = "")
    as_taxa_out[i, 3] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,4] == ""){
    x <- paste(as_taxa_out[i,3], "_unknown_class", sep = "")
    as_taxa_out[i, 4] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,5] == ""){
    x <- paste(as_taxa_out[i,4], "_unknown_order", sep = "")
    as_taxa_out[i, 5] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,6] == ""){
    x <- paste(as_taxa_out[i,5], "_unknown_family", sep = "")
    as_taxa_out[i, 6] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,7] == ""){
    x <- paste(as_taxa_out[i,6], "_unknown_genus", sep = "")
    as_taxa_out[i, 7] <- x
  }
}
for (i in 1:nrow(as_taxa_out)){
  if (as_taxa_out[i,8] == ""){
    x <- paste(as_taxa_out[i,7], "_unknown_species", sep = "")
    as_taxa_out[i, 8] <- x
  }
}

as_taxa_out$composite_species <- paste(as_taxa_out$Kingdom, as_taxa_out$Supergroup, as_taxa_out$Division, as_taxa_out$Class, as_taxa_out$Order, as_taxa_out$Family, as_taxa_out$Genus, as_taxa_out$Species, sep = ";")
as_taxa_out$composite_genus <- paste(as_taxa_out$Kingdom, as_taxa_out$Supergroup, as_taxa_out$Division, as_taxa_out$Class, as_taxa_out$Order, as_taxa_out$Family, as_taxa_out$Genus, sep = ";")
# add taxa table to phyloseq object
tax_table(prev_ps_t) <- as.matrix(as_taxa_out)
prev_ps_t

#

rank_names(prev_ps_t) # what are the taxa at each rank?
table(tax_table(prev_ps_t)[, "Supergroup"], exclude = NULL) # identify outgroup at division level
table(tax_table(prev_ps_t)[, "Division"], exclude = NULL) # identify outgroup at division level
table(tax_table(prev_ps_t)[, "Class"], exclude = NULL)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 
#### root phylogeny with Hacrobia ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

colSums(otu_table(subset_taxa(prev_ps_t ,Supergroup == "Hacrobia")))
out_hac <- getSequences(otu_table(subset_taxa(prev_ps_t ,Supergroup == "Hacrobia")))[1]
is.rooted(phy_tree(prev_ps_t))
fitGTR$tree <- root(phy_tree(prev_ps_t),
                    outgroup = out_hac, 
                    resolve.root = TRUE)

#Merge phylogenetic tree to phyloseq object:
phy_tree(prev_ps_t) <- fitGTR$tree
table(tax_table(prev_ps_t)[, "Class"], exclude = NULL)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Rename ASVs ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

taxa_names(prev_ps_t) <-  paste("ASV_", 1:ntaxa(prev_ps_t), sep ="")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Prepare for GLM ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##### remove single occurrences ----

temp <- otu_table(prev_ps_t)
temp[temp > 0] <- 1
keep <- colnames(temp[, colSums(temp) > 2])

x1 = prune_taxa(colnames(otu_table(prev_ps_t)) %in% keep, prev_ps_t) 
is.numeric(otu_table(x1))

#### Agglomerate ----

x1_aglom = tip_glom(x1, h = 0.025)
x1_aglom

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# manyGLM ---- 
# prev_ps_t matrix
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

write.table(as.data.frame(otu_table(x1_aglom)), file ="ps_filt_aglom.csv", sep = ",")
mv_table <- read.table(file="ps_filt_aglom.csv",  header = TRUE, sep = ",", row.names = 1, check.names=FALSE)

mv_table <- mvabund(mv_table)
meanvar.plot(mv_table)

mod1 <- manyglm(mv_table ~
                  nutrient +
                  temp_oxyg +
                  salinity +
                  silicate +
                  oxygen_sat +
                  offset(log(count)),
                data = meta_min,
                family = "negative binomial")
plot(mod1)
# Drop variables to minimise AIC
drop1(mod1)
x <- c(drop1(mod1)[[2]] > min(drop1(mod1)[[2]])) # Minimise AIC
y <- rownames(drop1(mod1))
y[!x]
mod2 <- update(mod1, .~. - silicate); drop1(mod2)
x <- c(drop1(mod2)[[2]] > min(drop1(mod2)[[2]]))
y <- rownames(drop1(mod2))
y[!x]
plot(mod2)

set.seed(123)
mv_aglom_9999 <- anova(mod2 , p.uni="adjusted", nBoot=9999)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# graph significant manyGLM ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

write.csv(t(mv_aglom_9999$uni.p), file = "mv_3_aglom_9999$uni.p.csv")
write.csv(mv_aglom_9999$table, file = "mv_3_aglom_9999$table.csv")
write.csv((as.data.frame(cbind(colSums(x1_aglom@otu_table), colMeans(x1_aglom@otu_table)))), file = "prev_ps_t$mean.csv")
# get significant ASVs

x <- as.data.frame(t(mv_aglom_9999$uni.p))
sig_nut <- rownames(x[x$nutrient < 0.051, ]) # make a vector of significant ASVs
rownames(x[x$nutrient < 0.05, ])
rownames(x[x$nutrient > 0.05 & x$nutrient < 0.1, ])
sig_temp <- rownames(x[x$temp_oxyg < 0.051, ]) # make a vector of significant ASVs
rownames(x[x$temp_oxyg < 0.05, ])
rownames(x[x$temp_oxyg > 0.05 & x$temp_oxyg < 0.1, ])
sig_nut_subset <- as.data.frame(t(subset(t(mv_table), colnames(mv_table) %in% sig_nut)))
paste(colnames(mv_table)[colnames(mv_table) %in% sig_nut], (tax_table(x1_aglom)[,3])[colnames(otu_table(x1_aglom)) %in% sig_nut], (tax_table(x1_aglom)[,7])[colnames(otu_table(x1_aglom)) %in% sig_nut], sep = " ")
sig_temp_subset <- as.data.frame(t(subset(t(mv_table), colnames(mv_table) %in% sig_temp)))
paste(colnames(mv_table)[colnames(mv_table) %in% sig_temp], (tax_table(x1_aglom)[,3])[colnames(otu_table(x1_aglom)) %in% sig_temp], (tax_table(x1_aglom)[,7])[colnames(otu_table(x1_aglom)) %in% sig_temp], sep = " ")
# what are the taxa?
write.csv(subset(prev_ps_t@tax_table, rownames(prev_ps_t@tax_table) %in% colnames(sig_nut_subset)), file = "3_sig_nut_subset_taxa.csv")
write.csv(subset(prev_ps_t@tax_table, rownames(prev_ps_t@tax_table) %in% colnames(sig_temp_subset)), file = "3_sig_temp_subset_taxa.csv")
sig_nut_subset$latitude <- meta_min$latitude
sig_temp_subset$latitude <- meta_min$latitude

d <- melt(sig_nut_subset, id.vars="latitude")
ggplot(d, aes(latitude,value)) +
  geom_point() +
  stat_smooth() +
  scale_x_reverse() +
  facet_wrap(~variable, ncol = 2) +
  scale_y_continuous(trans='log')
rm(d)

d <- melt(sig_temp_subset, id.vars="latitude")
ggplot(d, aes(latitude,value)) +
  geom_point() +
  stat_smooth() +
  scale_x_reverse() +
  facet_wrap(~variable, ncol = 2) +
  scale_y_continuous(trans='log')
rm(d)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Summary stats ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sort(rowSums(otu_table(prev_ps_t)))
sort(colSums(otu_table(prev_ps_t)))
sort(colSums(otu_table(prev_ps.percent)/sum(otu_table(prev_ps.percent))), decreasing = TRUE)[1:20]
sum(sort(colSums(otu_table(prev_ps.percent)/sum(otu_table(prev_ps.percent))), decreasing = TRUE)[1:5])
sum(sort(colSums(otu_table(prev_ps.percent)/sum(otu_table(prev_ps.percent))), decreasing = TRUE)[1:18])

prev_ps.percent = transform_sample_counts(prev_ps_t, function(x) 100 * x/sum(x))

d <- as.data.frame(otu_table(prev_ps.percent))[1:36]
d$latitude <- meta$latitude

d <- melt(d, id.vars="latitude")

ggplot(d, aes(latitude,value)) + 
  geom_point() + 
  stat_smooth() +
  facet_wrap(~variable) +
  scale_y_continuous(trans='log10')
rm(d)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# subsample ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#### subset by taxonomy ----
ps_met = subset_taxa(prev_ps_t, Supergroup == "Opisthokonta")
ps_Hexanauplia = subset_taxa(prev_ps_t, Class == "Hexanauplia")

##### rarefy ----

min_rar <- min(rowSums(otu_table(prev_ps_t)))/2
set.seed(123)
ps_rare = rarefy_even_depth(prev_ps_t, replace=TRUE, sample.size = min_rar)

min_rar <- min(rowSums(otu_table(ps_met)))/2
set.seed(123)
ps_met_rare = rarefy_even_depth(ps_met, replace=TRUE, sample.size = min_rar)

min_rar <- min(rowSums(otu_table(ps_Hexanauplia)))/2
set.seed(123)
ps_Hex_rare = rarefy_even_depth(ps_Hexanauplia, replace=TRUE, sample.size = min_rar)

#### transform ----
prev_ps_t_rel_abund <- transform_sample_counts(ps_rare, function(x) 100*x/sum(x))
prev_ps_t_rel_4th <- transform_sample_counts(prev_ps_t_rel_abund, function(x) x^(1/4))

ps_met_rel_abund <- transform_sample_counts(ps_met_rare, function(x) 100*x/sum(x))
ps_met_rel_4th <- transform_sample_counts(ps_met_rel_abund, function(x) x^(1/4))

ps_Hexanauplia_rel_abund <- transform_sample_counts(ps_Hex_rare, function(x) 100*x/sum(x))
ps_Hexanauplia_rel_4th <- transform_sample_counts(ps_Hexanauplia_rel_abund, function(x) x^(1/4))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Distance metrics ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Phylogenetic diversity and Shannon entropy ----
library(picante)

par(mfrow=c(3,1))

df.pd <- pd(as.data.frame(prev_ps_t_rel_4th@otu_table), prev_ps_t_rel_4th@phy_tree, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
df.pd$H <- diversity(as.data.frame(prev_ps_t_rel_4th@otu_table), "shannon")
df.pd
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(df.pd$PD ~ meta_min$latitude, xlab = "", ylab = "", main = "Phylogenetic a-diversity and entropy - all Eukarya")
par(new = TRUE) # tell R to not forget old plot
plot(df.pd$H ~ meta_min$latitude, axes = FALSE, bty = "n", xlab = "", ylab = "", col="red")
axis(side=4, at = pretty(range(df.pd$H)), col="red", col.axis="red")
mtext("latitude °S", side=1, line=3)
mtext("PD", side=2, line=3)
mtext("H", side=4, line=3, col="red")

df.pd <- pd(as.data.frame(ps_met_rel_4th@otu_table), ps_met_rel_4th@phy_tree, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
df.pd$H <- diversity(as.data.frame(ps_met_rel_4th@otu_table), "shannon")
df.pd
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(df.pd$PD ~ meta_min$latitude, xlab = "", ylab = "", main = "Phylogenetic a-diversity and entropy - all metazoa")
par(new = TRUE) # tell R to not forget old plot
plot(df.pd$H ~ meta_min$latitude, axes = FALSE, bty = "n", xlab = "", ylab = "", col="red")
axis(side=4, at = pretty(range(df.pd$H)), col="red", col.axis="red")
mtext("latitude °S", side=1, line=3)
mtext("PD", side=2, line=3)
mtext("H", side=4, line=3, col="red")


# plot(df.pd$SR ~ meta_min$latitude)


df.pd <- pd(as.data.frame(ps_Hexanauplia_rel_4th@otu_table), ps_Hexanauplia_rel_4th@phy_tree, include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
df.pd$H <- diversity(as.data.frame(ps_Hexanauplia_rel_4th@otu_table), "shannon")
df.pd
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(df.pd$PD ~ meta_min$latitude, xlab = "", ylab = "", main = "Phylogenetic a-diversity and entropy - all copepods")
par(new = TRUE) # tell R to not forget old plot
plot(df.pd$H ~ meta_min$latitude, axes = FALSE, bty = "n", xlab = "", ylab = "", col="red")
axis(side=4, at = pretty(range(df.pd$H)), col="red", col.axis="red")
mtext("latitude °S", side=1, line=3)
mtext("PD", side=2, line=3)
mtext("H", side=4, line=3, col="red")

#### ranked abundance top 50%, top 75%
ps_rare_taxa <- as.data.frame(tax_table(ps_rare)[,1:8]) 
ps_rare_taxa$relabund  <- 100*colSums(otu_table(ps_rare))/sum(otu_table(ps_rare))
ranked_abund_rare <- ps_rare_taxa[order(ps_rare_taxa$relabund, decreasing = TRUE), ]

sum(ps_rare_taxa$relabund[1])
sum(ps_rare_taxa$relabund[1:5])
sum(ps_rare_taxa$relabund[1:18])

#### table all ----
write.csv(ranked_abund_rare, file = "ranked_abund_rare.csv")

#### table  ----
Hexanauplia <- ranked_abund_rare[ranked_abund_rare$Supergroup == "Opisthokonta" & ranked_abund_rare$Class == "Hexanauplia", ]
write.csv(Hexanauplia, file = "Hexanauplia.csv")
Opist_not_Hexanauplia <- ranked_abund_rare[ranked_abund_rare$Supergroup == "Opisthokonta" & ranked_abund_rare$Class != "Hexanauplia", ]
write.csv(Opist_not_Hexanauplia, file = "Opist_not_Hexanauplia.csv")
Ciliophora <- ranked_abund_rare[ranked_abund_rare$Division == "Ciliophora" , ]
write.csv(Ciliophora, file = "Ciliophora")
Dinoflagellata <- ranked_abund_rare[ranked_abund_rare$Division == "Dinoflagellata" , ]
write.csv(Dinoflagellata, file = "Dinoflagellata")

Rhizaria <- ranked_abund_rare[ranked_abund_rare$Supergroup == "Rhizaria" , ]
write.csv(Rhizaria, file = "Rhizaria.csv")
Hacrobia <- ranked_abund_rare[ranked_abund_rare$Supergroup == "Hacrobia" , ]
write.csv(Hacrobia, file = "Hacrobia.csv")

Alveolata <- ranked_abund_rare[ranked_abund_rare$Supergroup == "Alveolata" , ]
write.csv(Alveolata, file = "Alveolata_ALL.csv")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# ordinations ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
par(mfrow=c(1,1))

### Make color vectors for plots
x <- as.data.frame(table(tax_table(prev_ps_t_rel_4th)[, "Supergroup"], exclude = NULL))
Supergroup <- as.character(x$Var1)             #make vector to unify coloration of legend in graphs
Supergroup_cols <- c("#66c2a5", "#ffff99", "#fb8072", "#8da0cb")
x <- as.data.frame(table(tax_table(ps_met)[, "Class"], exclude = NULL))
Class <- as.character(x$Var1)             #make vector to unify coloration of legend in graphs
# Class_cols <- c("#8dd3c7", "#fb8072", "#ffffb3", "#80b1d3", "#fdb462", "#bebada")

library("dichromat")
Class_cols <- palette(brewer.pal(n = 11, name = "BrBG"))
Class_cols <- colorRampPalette(Class_cols)(length(Class)+1)
pie(rep(1, length(Class_cols)), col = Class_cols , main="") 
Class_cols <- Class_cols[-7]
# Zissou palette
as_taxa_out <- as.data.frame(tax_table(ps_Hexanauplia_rel_4th))
as_taxa_out$composite_family <- paste(as_taxa_out$Kingdom, as_taxa_out$Supergroup, as_taxa_out$Division, as_taxa_out$Class, as_taxa_out$Order, as_taxa_out$Family, sep = ";")
# add taxa table to phyloseq object
tax_table(ps_Hexanauplia_rel_4th) <- as.matrix(as_taxa_out)

Family <- sort(unique(tax_table(ps_Hexanauplia_rel_4th)[,  "composite_family"]))

library("wesanderson")
Calanoida_cols <- wes_palette("Zissou1", length(grep("Calanoida", Family)), type = "continuous")
Cyclopoida_cols <- wes_palette("Zissou1", length(grep("Cyclopoida", Family)), type = "continuous")
Harpacticoida_cols <- wes_palette("Zissou1", length(grep("Harpacticoida", Family)), type = "continuous")
Hexanauplia_unknown_order_cols <- "#737373"
Poecilostomatoida_cols <- wes_palette("Zissou1", length(grep("Poecilostomatoida", Family)), type = "continuous")
Family_cols <- c(Calanoida_cols, Cyclopoida_cols, Harpacticoida_cols, Hexanauplia_unknown_order_cols, Poecilostomatoida_cols)
  
# create an "other" colour
other <- "other"
other_cols <- "#737373"

taxa_col <- c(Supergroup_cols, other_cols)
names(taxa_col) <- c(Supergroup,  other)

Metazoa_col <- c(Class_cols, other_cols)
names(Metazoa_col) <- c(Class, other)

names(Family_cols) <- Family

plot(seq_len(length(taxa_col)), rep_len(1, length(taxa_col)),
     col = taxa_col, pch = 16, cex = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
plot(seq_len(length(Metazoa_col)), rep_len(1, length(Metazoa_col)),
     col = Metazoa_col, pch = 16, cex = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
plot(seq_len(length(Family_cols)), rep_len(1, length(Family_cols)),
     col = Family_cols, pch = 16, cex = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')

#### All taxa DPCoA ----

dpcoa_dist_all = phyloseq::distance(prev_ps_t_rel_4th, method = "dpcoa")
out.dpcoa.dpcoa <- ordinate(prev_ps_t_rel_4th, method = "DPCoA", distance = dpcoa_dist_all)
evals <- out.dpcoa.dpcoa$eig
plot_scree(out.dpcoa.dpcoa)

Axes_124_Full <- out.dpcoa.dpcoa$tab[1:4]

prev_ps_sam = plot_ordination(prev_ps_t_rel_4th, out.dpcoa.dpcoa, color = "latitude") +
  labs(col = "latitude Age") +
  theme_bw() + 
  theme(legend.position="none") +
  geom_text(mapping = aes(label = sampleID), size = 1) #+
  #geom_point(size = 4, shape = 1, alpha = 0.6) 
prev_ps_sam$layers <- prev_ps_sam$layers[-1]
prev_ps_tax = plot_ordination(prev_ps_t_rel_4th, color = "Supergroup", out.dpcoa.dpcoa, type = "species") +
  theme_bw() +
  theme(legend.position="left") +
  labs(color = "Supergroup") +
  scale_color_manual(values = taxa_col) +
  geom_point(size = 4, shape = 16, alpha = 0.6) 
prev_ps_tax$layers <- prev_ps_tax$layers[-1]
prev_ps_tax

#### Metazoa DPCoA ----

dpcoa_dist_met = phyloseq::distance(ps_met_rel_4th, method = "dpcoa")
out.dpcoa.dpcoa <- ordinate(ps_met_rel_4th, method = "DPCoA", distance = dpcoa_dist_met)
evals <- out.dpcoa.dpcoa$eig

plot_scree(out.dpcoa.dpcoa)
Axes_124_metazoans <- out.dpcoa.dpcoa$tab[1:4]

ps_met_sam = plot_ordination(ps_met_rel_4th, out.dpcoa.dpcoa, color = "latitude") +
  labs(col = "latitude Age") +
  theme_bw() + 
  theme(legend.position="none") +
  geom_text(mapping = aes(label = sampleID), size = 1) #+
  #geom_point(size = 4, shape = 1, alpha = 0.6) 
ps_met_sam$layers <- ps_met_sam$layers[-1]
ps_met_tax = plot_ordination(ps_met_rel_4th, color = "Class", out.dpcoa.dpcoa, type = "species") +
  theme_bw() +
  theme(legend.position="left") +
  labs(color = "Class") +
  scale_color_manual(values = Metazoa_col) +
  geom_point(size = 4, shape = 16, alpha = 0.6) 
ps_met_tax$layers <- ps_met_tax$layers[-1]
ps_met_tax + theme(legend.position="none")
ps_met_tax 

#### Hexanauplia DPCoA ----

dpcoa_dist_art = phyloseq::distance(ps_Hexanauplia_rel_4th, method = "dpcoa")
out.dpcoa.dpcoa <- ordinate(ps_Hexanauplia_rel_4th, method = "DPCoA", distance = dpcoa_dist_art)
evals <- out.dpcoa.dpcoa$eig

plot_scree(out.dpcoa.dpcoa)
Axes_124_Arthropoda <- out.dpcoa.dpcoa$tab[1:4]

ps_Hexanauplia_sam = plot_ordination(ps_Hexanauplia_rel_4th, out.dpcoa.dpcoa, color = "latitude") +
  labs(col = "latitude Age") +
  theme_bw() + 
  theme(legend.position="none") +
  geom_text(mapping = aes(label = sampleID), size = 4) #+
  #geom_point(size = 4, shape = 1, alpha = 0.6) 
ps_Hexanauplia_sam$layers <- ps_Hexanauplia_sam$layers[-1]
ps_Hexanauplia_tax = plot_ordination(ps_Hexanauplia_rel_4th, color = "composite_family", out.dpcoa.dpcoa, shape = "Order" ,type = "species") +
  theme_bw() +
  labs(color = "Family") +
  scale_shape_manual(values = c(16,15,17,0,18)) +
  scale_color_manual(values = Family_cols) +
  geom_point(size = 4, alpha = 0.6) 
ps_Hexanauplia_tax$layers <- ps_Hexanauplia_tax$layers[-1]
ps_Hexanauplia_tax + theme(legend.position="bottom")

#### group plot ----

library("ggpubr")
ggpubr::ggarrange(prev_ps_sam, prev_ps_tax + theme(legend.position="none"), ps_met_sam, ps_met_tax + theme(legend.position="none"), ps_Hexanauplia_sam, ps_Hexanauplia_tax + theme(legend.position="none"),
          ncol = 2, 
          nrow = 3)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Heat map() ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(phytools)

par(mfrow=c(1,1))

plot(x1_aglom@phy_tree, edge.width = 2, label.offset = 2) + 
  nodelabels(col = "blue", frame = "none")

phylo.heatmap(x1_aglom@phy_tree, t(as.data.frame(x1_aglom@otu_table))^(1/4))

sort(100*colSums(otu_table(x1_aglom))/sum(otu_table(x1_aglom)))

plot_tree(x1_aglom, label.tips = "Class", color="sampleID", plot.margin=0.1, title = "By Height")

detach("package:phytools", unload=TRUE) 

library(ggOceanMaps)
metadata_mat$Latitude <- metadata_mat$Latitude * -1
x <- basemap(limits = c(95, 135, -5, -45), bathymetry = F, legends = F)  
x +  geom_spatial_point(data = metadata_mat, aes(x = Longitude, y = Latitude), col = "blue", size = 4, shape = 21 )


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
## Fuzzy c-means clustering  ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(cluster)

par(mfrow=c(1,1))

fuzzy.2 <- fanny(dpcoa_dist_all, k = 2, memb.exp = 1.5)
two_k <- as.data.frame(cbind(fuzzy.2$membership, 1:20, rep("k=2", length(rownames(fuzzy.2$membership)))))
colnames(two_k) <- c("1", "2", "station", "k")
rownames(two_k) <- paste(rownames(two_k), two_k$k, sep = "_")
fuzzy.3 <- fanny(dpcoa_dist_all, k = 3, memb.exp = 1.5)
three_k <- as.data.frame(cbind(fuzzy.3$membership, 1:20, rep("k=3", length(rownames(fuzzy.3$membership)))))
colnames(three_k) <- c("1", "2", "3", "station", "k")
rownames(three_k) <- paste(rownames(three_k), three_k$k, sep = "_")
fuzzy.4 <- fanny(dpcoa_dist_all, k = 4, memb.exp = 1.5)
four_k <- as.data.frame(cbind(fuzzy.4$membership, 1:20, rep("k=4", length(rownames(fuzzy.4$membership)))))
colnames(four_k) <- c("1", "2", "3", "4", "station", "k")
rownames(four_k) <- paste(rownames(four_k), four_k$k, sep = "_")
# sandwich and melt
x <- as.data.frame(bind_rows(two_k, three_k, four_k))
z <- melt(x, id.vars=c("station", "k"))
z$value <- as.numeric(as.character(z$value)) # value is not numeric, fix this
z$station <- as.numeric(as.character(z$station)) # value is not numeric, fix this

# Stacked + percent
ggplot(data = z, aes(fill=variable, y=value, x=station)) + 
  geom_bar(position="fill", stat="identity") + 
  #scale_x_reverse() +
  ggtitle(label="All taxa fuzzy cluster") +
  coord_flip() + 
  facet_wrap(~k) 

## Fuzzy c-means clustering metazoa ----
fuzzy.2 <- fanny(dpcoa_dist_met, k = 2, memb.exp = 1.5)
two_k <- as.data.frame(cbind(fuzzy.2$membership, 1:20, rep("k=2", length(rownames(fuzzy.2$membership)))))
colnames(two_k) <- c("1", "2", "station", "k")
rownames(two_k) <- paste(rownames(two_k), two_k$k, sep = "_")
fuzzy.3 <- fanny(dpcoa_dist_met, k = 3, memb.exp = 1.5)
three_k <- as.data.frame(cbind(fuzzy.3$membership, 1:20, rep("k=3", length(rownames(fuzzy.3$membership)))))
colnames(three_k) <- c("1", "2", "3", "station", "k")
rownames(three_k) <- paste(rownames(three_k), three_k$k, sep = "_")
fuzzy.4 <- fanny(dpcoa_dist_met, k = 4, memb.exp = 1.5)
four_k <- as.data.frame(cbind(fuzzy.4$membership, 1:20, rep("k=4", length(rownames(fuzzy.4$membership)))))
colnames(four_k) <- c("1", "2", "3", "4", "station", "k")
rownames(four_k) <- paste(rownames(four_k), four_k$k, sep = "_")
# sandwich and melt
x <- as.data.frame(bind_rows(two_k, three_k, four_k))
z <- melt(x, id.vars=c("station", "k"))
z$value <- as.numeric(as.character(z$value)) # value is not numeric, fix this
z$station <- as.numeric(as.character(z$station)) # value is not numeric, fix this

# Stacked + percent
ggplot(data = z, aes(fill=variable, y=value, x=station)) + 
  geom_bar(position="fill", stat="identity") + 
  #scale_x_reverse() +
  ggtitle(label="Animal fuzzy cluster") +
  coord_flip() + 
  facet_wrap(~k) 

## Fuzzy c-means clustering arthropods ----
fuzzy.2 <- fanny(dpcoa_dist_art, k = 2, memb.exp = 1.5)
two_k <- as.data.frame(cbind(fuzzy.2$membership, 1:20, rep("k=2", length(rownames(fuzzy.2$membership)))))
colnames(two_k) <- c("1", "2", "station", "k")
rownames(two_k) <- paste(rownames(two_k), two_k$k, sep = "_")
fuzzy.3 <- fanny(dpcoa_dist_art, k = 3, memb.exp = 1.5)
three_k <- as.data.frame(cbind(fuzzy.3$membership, 1:20, rep("k=3", length(rownames(fuzzy.3$membership)))))
colnames(three_k) <- c("1", "2", "3", "station", "k")
rownames(three_k) <- paste(rownames(three_k), three_k$k, sep = "_")
fuzzy.4 <- fanny(dpcoa_dist_art, k = 4, memb.exp = 1.5, maxit = 1000)
four_k <- as.data.frame(cbind(fuzzy.4$membership, 1:20, rep("k=4", length(rownames(fuzzy.4$membership)))))
colnames(four_k) <- c("1", "2", "3", "4", "station", "k")
rownames(four_k) <- paste(rownames(four_k), four_k$k, sep = "_")
# sandwich and melt
x <- as.data.frame(bind_rows(two_k, three_k, four_k))
z <- melt(x, id.vars=c("station", "k"))
z$value <- as.numeric(as.character(z$value)) # value is not numeric, fix this
z$station <- as.numeric(as.character(z$station)) # value is not numeric, fix this

# Stacked + percent
ggplot(data = z, aes(fill=variable, y=value, x=station)) + 
  geom_bar(position="fill", stat="identity") + 
  #scale_x_reverse() +
  coord_flip() + 
  ggtitle(label="Copepod fuzzy cluster") +
  facet_wrap(~k) 

## Fuzzy c-means clustering arthropods Bray Curtis ----

braycurtis_dist_art = phyloseq::distance(ps_Hexanauplia_rel_4th, method = "bray")

fuzzy.2 <- fanny(braycurtis_dist_art, k = 2, memb.exp = 1.5)
two_k <- as.data.frame(cbind(fuzzy.2$membership, 1:20, rep("k=2", length(rownames(fuzzy.2$membership)))))
colnames(two_k) <- c("1", "2", "station", "k")
rownames(two_k) <- paste(rownames(two_k), two_k$k, sep = "_")
fuzzy.3 <- fanny(braycurtis_dist_art, k = 3, memb.exp = 1.5)
three_k <- as.data.frame(cbind(fuzzy.3$membership, 1:20, rep("k=3", length(rownames(fuzzy.3$membership)))))
colnames(three_k) <- c("1", "2", "3", "station", "k")
rownames(three_k) <- paste(rownames(three_k), three_k$k, sep = "_")
fuzzy.4 <- fanny(braycurtis_dist_art, k = 4, memb.exp = 1.5, maxit = 1000)
four_k <- as.data.frame(cbind(fuzzy.4$membership, 1:20, rep("k=4", length(rownames(fuzzy.4$membership)))))
colnames(four_k) <- c("1", "2", "3", "4", "station", "k")
rownames(four_k) <- paste(rownames(four_k), four_k$k, sep = "_")
# sandwich and melt
x <- as.data.frame(bind_rows(two_k, three_k, four_k))
z <- melt(x, id.vars=c("station", "k"))
z$value <- as.numeric(as.character(z$value)) # value is not numeric, fix this
z$station <- as.numeric(as.character(z$station)) # value is not numeric, fix this

# Stacked + percent
ggplot(data = z, aes(fill=variable, y=value, x=station)) + 
  geom_bar(position="fill", stat="identity") + 
  #scale_x_reverse() +
  ggtitle(label="Copepod Bray Curtis") +
  coord_flip() + 
  facet_wrap(~k) 

## Fuzzy c-means clustering Environment ----
fuzzy.2 <- fanny(meta_min[,1:5], k = 2, memb.exp = 1.5)
two_k <- as.data.frame(cbind(fuzzy.2$membership, 1:20, rep("k=2", length(fuzzy.2$clustering))))
colnames(two_k) <- c("1", "2", "station", "k")
rownames(two_k) <- paste(rownames(two_k), two_k$k, sep = "_")
fuzzy.3 <- fanny(meta_min[,1:5], k = 3, memb.exp = 1.5)
three_k <- as.data.frame(cbind(fuzzy.3$membership, 1:20, rep("k=3", length(fuzzy.2$clustering))))
colnames(three_k) <- c("1", "2", "3", "station", "k")
rownames(three_k) <- paste(rownames(three_k), three_k$k, sep = "_")
fuzzy.4 <- fanny(meta_min[,1:5], k = 4, memb.exp = 1.5, maxit = 1000)
four_k <- as.data.frame(cbind(fuzzy.4$membership, 1:20, rep("k=4", length(fuzzy.2$clustering))))
colnames(four_k) <- c("1", "2", "3", "4", "station", "k")
rownames(four_k) <- paste(rownames(four_k), four_k$k, sep = "_")
# sandwich and melt
x <- as.data.frame(bind_rows(two_k, three_k, four_k))
z <- melt(x, id.vars=c("station", "k"))
z$value <- as.numeric(as.character(z$value)) # value is not numeric, fix this
z$station <- as.numeric(as.character(z$station)) # value is not numeric, fix this

# Stacked + percent
ggplot(data = z, aes(fill=variable, y=value, x=station)) + 
  geom_bar(position="fill", stat="identity") + 
  #scale_x_reverse() +
  coord_flip() + 
  ggtitle(label="Environmental") +
  facet_wrap(~k) 

#@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Vegemite and Copepoda PA map ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@

library(gplots)
art_df_RA <- as.matrix(otu_table(ps_Hexanauplia_rel_4th))
colnames(art_df_RA) <- paste(tax_table(ps_Hexanauplia_rel_4th)[,5], colnames(art_df_RA))
vegemite(x=art_df_RA, scale="Hult") # quick look

art_df_PA <- art_df_RA
art_df_PA[art_df_PA > 0] <- 1 # make binary
art_df_PA <- art_df_PA[nrow(art_df_PA):1, ] # reverse order of rows

set.seed(123)
h.2 <- heatmap.2(art_df_PA,
          dendrogram = "column",
        col = c("white", "black"), # for binary data
        scale = "none", 
        density.info = "none", 
        Rowv = FALSE,
        Colv = TRUE, 
        trace = "none",
        margins = c(8, 4)
        )
col_order <- colnames(art_df_PA)[h.2$colInd]

#

library("indicspecies")
# a: Warmer, low salinity water = black, 
# b: Warm, higher salinity water = purple,
# c: Cool, higher salinity water = orange, 
# d: Cold, high nitrate + nitrite water = yellow
meta_min$groups <- c("d","d","d","c","c","c","b","b","b","a","a","a","a","a","a","a","a","a","a","a")

set.seed(123)
indval_RA = multipatt(art_df_RA, meta_min$groups, control = how(nperm=999))
summary(indval_RA)
set.seed(123)
indval_PA = multipatt(art_df_PA, meta_min$groups, control = how(nperm=999))
summary(indval_PA)

RA_sig <- indval_RA$sign[c(1:4,7)]
RA_sig[5][is.na(RA_sig[5]),] <- 1 # make sure NA results are recorded as not significant
RA_sig[1:4][RA_sig[5] > 0.05, ] <- 0 # Don't graph insignificant (turn to zero)
RA_sig <- as.matrix(RA_sig[1:4][col_order,]) # order significant species so they match the dendrogram

heatmap.2(t(RA_sig),
          col = c("white", "RED"), # for binary data
          scale = "none", 
          density.info = "none", 
          Rowv = FALSE,
          Colv = FALSE, 
          trace = "none",
          margins = c(8, 8)
)

PA_sig <- indval_PA$sign[c(1:4,7)]
PA_sig[5][is.na(PA_sig[5]),] <- 1 # make sure NA results are recorded as not significant
PA_sig[1:4][PA_sig[5] > 0.05, ] <- 0 # Don't graph insignificant (turn to zero)
PA_sig <- as.matrix(PA_sig[1:4][col_order,]) # order significant species so they match the dendrogram

heatmap.2(t(PA_sig),
          col = c("white", "blue"), # for binary data
          scale = "none", 
          density.info = "none", 
          Rowv = FALSE,
          Colv = FALSE, 
          trace = "none",
          margins = c(8, 8)
)

###############################

# piecewise regression ----

###############################

#### choose a breakpoint  ----
#   Take a stab at where the breakpoints is, then choose values either side
#   Create a variable called breaks to hold these breakpoints:
x <- meta_min$latitude
# breaks <- x[which(x >= min(x) & x <= max(x))] 
breaks <- unique(x)
breaks <- (breaks[2:length(breaks)] + breaks[1:length(breaks)-1]) / 2 # shift breakpoints between observed sites

# Now iteratively search these breakpoints for the model that has the lowest residual MSE

##################################################################################

#### temperature/oxygen ----

breaks <- c(0,1,2,3)
r_sq_best <- as.data.frame(breaks)

lm_temp_oxyg <- lm(meta_min$temp_oxyg~x)

par(mfrow=c(4,1), mai=c(.1,.6,.1,.1))
plot(data=meta_min, temp_oxyg ~ x, xlab="latitude", col="blue")
lines(c(min(x),max(x)),c(meta_min$temp_oxyg[20],meta_min$temp_oxyg[1]),col="red")


mse <- numeric(length(breaks))
for(i in 1:length(breaks)){
  piecewise1 <- lm(meta_min$temp_oxyg ~ x*(x<breaks[i]) + x*(x>=breaks[i]))
  mse[i] <- summary(piecewise1)[6]
}
mse <- as.numeric(mse)

mse2 <- numeric(length(combn(breaks,2)[1,]))
for(i in 1:length(mse2)){
  j <- min(combn(breaks,2)[1:2,i])
  k <- max(combn(breaks,2)[1:2,i])
  piecewise2 <- lm(meta_min$temp_oxyg ~ x*(x<j) + x*(j<=x & x<k) + x*(x>=k))
  mse2[i] <- summary(piecewise2)[6]
}
mse2 <- as.numeric(mse2)

mse3 <- numeric(length(combn(breaks,3)[1,]))  # there are 969 combinations of three
for(i in 1:length(mse3)){                     #
  j <- min(combn(breaks,3)[1:3,i])
  k <- median(combn(breaks,3)[1:3,i])
  l <- max(combn(breaks,3)[1:3,i])
  piecewise3 <- lm(meta_min$temp_oxyg ~ x*(x<j) + x*(j<=x & x<k) + x*(k<=x & x<l) + x*(x>=l))
  mse3[i] <- summary(piecewise3)[6]
}
mse3 <- as.numeric(mse3)

# pick the breakpoint with the lowest error:

best <- breaks[which(mse==min(mse))]
best 
best2 <- sort(combn(breaks,2)[1:2,which(mse2==min(mse2))])
best2 
best3 <- sort(combn(breaks,3)[1:3,which(mse3==min(mse3))])
best3 

temp_oxyg <- lm(meta_min$temp_oxyg~x*(x<best)+x*(x>=best))

summary(temp_oxyg)

a1 <- summary(temp_oxyg)[[4]][1]+summary(temp_oxyg)[[4]][3] # intercept 1
a2 <- summary(temp_oxyg)[[4]][1] # intercept 2
b1 <- summary(temp_oxyg)[[4]][2]+summary(temp_oxyg)[[4]][4] # slope 1
b2 <- summary(temp_oxyg)[[4]][2] # slope 2

plot(data=meta_min, temp_oxyg ~ x, xlab="latitude")
rect(ci$normal[2],min(meta_min$temp_oxyg),ci$normal[3],26, col="red", lwd=0.1, density=30, border = "transparent")
lines(c(min(x),best),c(a1+b1*min(x),a1+b1*best),col="red")
lines(c(best,max(x)),c(a2+b2*best,a2+b2*max(x)),col="red")
abline(v=best, lty=3, col="red")

# 2 breakpoints
temp_oxyg2 <- lm(meta_min$temp_oxyg ~ (x*(x<min(best2)) + x*(min(best2)<=x & x<max(best2))) + x*(x>=max(best2)))
summary(temp_oxyg2)

sub <- meta_min[x<=max(best2),]

a1 <- summary(temp_oxyg2)[[4]][1]+summary(temp_oxyg2)[[4]][3] # intercept 1
a2 <- summary(temp_oxyg2)[[4]][1]+summary(temp_oxyg2)[[4]][4]
a3 <- summary(temp_oxyg2)[[4]][1] # intercept 3

b1 <- summary(temp_oxyg2)[[4]][2]+summary(temp_oxyg2)[[4]][5] # slope 1
b2 <- summary(temp_oxyg2)[[4]][2]+summary(temp_oxyg2)[[4]][6]
b3 <- summary(temp_oxyg2)[[4]][2] # slope 3


plot(data=meta_min, temp_oxyg ~ x, xlab="latitude", col="blue")
lines(c(min(x),min(best2)),c(a1+b1*min(x),a1+b1*min(best2)),col="red")
lines(c(min(best2),max(best2)),c(a2+b2*min(best2),a2+b2*max(best2)),col="red")
lines(c(max(best2),max(x)),c(a3+b3*max(best2),a3+b3*max(x)),col="red")
abline(v=min(best2), lty=3, col="red")
abline(v=max(best2), lty=3, col="red")

# 3 breakpoints

temp_oxyg3 <- lm(meta_min$temp_oxyg ~ (x*(x<min(best3)) + x*(min(best3)<=x & x<median(best3)) + x*(median(best3)<=x & x<max(best3)) + x*(x>=max(best3))))
summary(temp_oxyg3)

sub <- meta_min[x<=max(best3),]

a1 <- summary(temp_oxyg3)[[4]][1]+summary(temp_oxyg3)[[4]][3] # intercept 1
a2 <- summary(temp_oxyg3)[[4]][1]+summary(temp_oxyg3)[[4]][4]
a3 <- summary(temp_oxyg3)[[4]][1]+summary(temp_oxyg3)[[4]][5]
a4 <- summary(temp_oxyg3)[[4]][1] # intercept 3

b1 <- summary(temp_oxyg3)[[4]][2]+summary(temp_oxyg3)[[4]][6] # slope 1
b2 <- summary(temp_oxyg3)[[4]][2]+summary(temp_oxyg3)[[4]][7]
b3 <- summary(temp_oxyg3)[[4]][2]+summary(temp_oxyg3)[[4]][8]
b4 <- summary(temp_oxyg3)[[4]][2] # slope 3

plot(data=meta_min, temp_oxyg ~ x, xlab="latitude")
lines(c(min(x),min(best3)),c(a1+b1*min(x),a1+b1*min(best3)),col="green")
lines(c(min(best3),median(best3)),c(a2+b2*min(best3),a2+b2*median(best3)),col="green")
lines(c(median(best3),max(best3)),c(a3+b3*median(best3),a3+b3*max(best3)),col="green")
lines(c(max(best3),max(x)),c(a4+b4*max(best3),a4+b4*max(x)),col="green")
abline(v=min(best3), lty=3, col="green")
abline(v=median(best3), lty=3, col="green")
abline(v=max(best3), lty=3, col="green")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Tree for supplementary materials ----
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

tt <- as.data.frame(tax_table(prev_ps_t))
tt$ASV_genus <- paste(rownames(tt), tt$Species)
tax_table(prev_ps_t) <- as.matrix(tt)
plot_tree(prev_ps_t, label.tips = "ASV_genus", color="latitude", ladderize = "left", plot.margin=0.1, title = "By Height")

tt <- as.data.frame(tax_table(x1_aglom))
tt$ASV_genus <- paste(rownames(tt), tt$Species)
tax_table(x1_aglom) <- as.matrix(tt)
plot_tree(x1_aglom, label.tips = "composite_genus", color="latitude", ladderize = "left", plot.margin=0.1, title = "By Height")
