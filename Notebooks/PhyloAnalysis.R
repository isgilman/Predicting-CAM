# A time stamp to console prompt
h <- taskCallbackManager()
h$add(function(expr, value, ok, visible) {
   options("prompt"=format(Sys.time(), "%H:%M:%S $ "));
   return(TRUE) },
   name = "simpleHandler")

library("ape")
library("geiger")
library("phytools")
library("RPANDA")
library("ggplot2")
library("ggtree")
library("treeio")
library("deeptime")
library("dplyr")
library(nlme)
library(caper)
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

## Load tree
setwd("~/Dropbox/GitHub_repos/Predicting-CAM/Notebooks/")
fullTree = read.beast("../Results/PhyloWeld.timetree.nex")
fullTree@phylo$tip.label = gsub("_", " ", fullTree@phylo$tip.label)
fullTree

## Load anatomy data
morphoData = read.csv("../Data/Portullugo-Anatomy-Data.2023-08-09.csv")
morphoData$tipLabel = gsub("_", " ", morphoData$tipLabel)
morphoData[,"logMA"] = log10(morphoData$MA)
morphoData[,"logLT"] = log10(morphoData$LT)
morphoData[,"logIAS"] = log10(morphoData$IAS)

## Match tree and data
haveData = morphoData[!is.na(morphoData$logMA) | !is.na(morphoData$logIAS) | !is.na(morphoData$logLT),]
rownames(haveData) = haveData$tipLabel
phy = drop.tip(fullTree, setdiff(fullTree@phylo$tip.label, haveData$tipLabel))
phy

# See https://github.com/mlandis/vib_div/blob/master/code/plot/plot_fig3_mcc.R for details on why this is needed
phy@data$age_0.95_HPD = lapply(phy@data$CI_height, function(z) { if (is.na(z)) { return(c(NA,NA)) } else { return(z) } })

################
### Plotting ###
################
# I'm not the greatest using ggtree, so I'll try to get 90% of the way there in R
# and do the rest in Illustrator

## Clades with calibration points
opuntiodeae = getMRCA(phy@phylo, c("Opuntia cochenillifera","Grusonia bradtiana"))
portulaca = getMRCA(phy@phylo, c("Portulaca umbraticola","Portulaca oleracea"))
didieroideae = getMRCA(phy@phylo, c("Alluaudia procera","Alluaudia dumosa"))
cactoideae = getMRCA(phy@phylo, c("Ariocarpus retusus","Copiapoa desertorum"))
coreCactaceae = getMRCA(phy@phylo, c("Ariocarpus retusus","Opuntia cochenillifera"))
cactaceae = getMRCA(phy@phylo, c("Ariocarpus retusus","Pereskia sacharosa 67"))
talinaceae = getMRCA(phy@phylo, c("Talinum paniculatum","Talinum fruticosum"))
acp = getMRCA(phy@phylo, c("Portulaca umbraticola","Tunilla corrugata"))
acpt = getMRCA(phy@phylo, c("Portulaca umbraticola","Talinum portulacifolium"))
montiaceae = getMRCA(phy@phylo,c("Phemeranthus teretifolius 400","Calandrinia flava 79"))
portulacineae = getMRCA(phy@phylo, c("Phemeranthus teretifolius 400","Tacinga lilae"))
molluginaceae = getMRCA(phy@phylo, c("Mollugo verticillata","Mollugo nudicaulis"))
portullugo = getMRCA(phy@phylo, c("Mollugo verticillata","Tacinga lilae"))

calNodes = c(opuntiodeae, portulaca, didieroideae, cactoideae, coreCactaceae, cactaceae, talinaceae, acp, acpt, montiaceae, portulacineae, molluginaceae, portullugo)
node.label.cols <- rep(NA, phy@phylo$Nnode + length(phy@phylo$tip.label))
node.label.cols[calNodes] <- 'Black' 
nodes.text <- rep(NA, phy@phylo$Nnode + length(phy@phylo$tip.label))
nodes.text[calNodes] <- c(1:length(calNodes))

## Get node data for plotting
p = ggtree(phy@phylo)
# p$data = p$data %>% mutate(node = as.character(node))
p$data = p$data %>% mutate(node = as.integer(node))
p$data = dplyr::left_join(p$data, phy@data, by="node")
p$data = p$data %>% mutate(node = as.numeric(node))

# size refers to edge width which is set to 0 here because q will be overwritten by p later 
hpd = p + geom_tree() + geom_rootedge(rootedge = 9) +
   coord_cartesian(xlim = c(-15, max(nodeHeights(phy@phylo)) * 1.2), ylim = c(-3.5, phy@phylo$Nnode * 1.05), expand = FALSE) +
   geom_range("age_0.95_HPD", color='blue', size=0.5, alpha=.3)
hpdDF = hpd$data[, c("CI_date", "CI_height", "x", "y")]
hpdDF$min = NA
hpdDF$max = NA
for(i in 1:nrow(hpdDF)){
   hpdDF$min[i] = max(nodeHeights(phy@phylo)) - hpdDF$CI_height[i][[1]][1]
   hpdDF$max[i] = max(nodeHeights(phy@phylo)) - hpdDF$CI_height[i][[1]][2]
}

## Define major epochs
greyEpochs = p +
   annotate(geom = "rect", xmin = max(nodeHeights(phy@phylo))-epochs$max_age[3], xmax = max(nodeHeights(phy@phylo)), ymin = -4, ymax = phy@phylo$Nnode * 1.05, fill = "gray95") +
   annotate(geom = "rect", xmin = max(nodeHeights(phy@phylo))-epochs$max_age[5], xmax = max(nodeHeights(phy@phylo))-epochs$min_age[5], ymin = -4, ymax = phy@phylo$Nnode * 1.05, fill = "gray95") +
   annotate(geom = "rect", xmin = max(nodeHeights(phy@phylo))-epochs$max_age[7], xmax = max(nodeHeights(phy@phylo))-epochs$min_age[7], ymin = -4, ymax = phy@phylo$Nnode * 1.05, fill = "gray95") +
   annotate(geom = "text",x = max(nodeHeights(phy@phylo))-mean(c(epochs$max_age[3], epochs$min_age[1])),y = -2,label = "Pliocene", hjust = 0.5, size = 3) +
   annotate(geom = "text",x = max(nodeHeights(phy@phylo))-mean(c(epochs$max_age[4], epochs$min_age[4])),y = -2,label = "Miocene",hjust = 0.5, size = 3) +
   annotate(geom = "text",x = max(nodeHeights(phy@phylo))-mean(c(epochs$max_age[5], epochs$min_age[5])),y = -2,label = "Oligocene",hjust = 0.5, size = 3) +
   annotate(geom = "text",x = max(nodeHeights(phy@phylo))-mean(c(epochs$max_age[6], epochs$min_age[6])),y = -2,label = "Eocene",hjust = 0.5, size = 3) +
   annotate(geom = "text",x = max(nodeHeights(phy@phylo))-mean(c(epochs$max_age[7], epochs$min_age[7])),y = -2,label = "Paleocene",hjust = 0.5, size = 3) 
greyEpochs

## Plot annotated tree
labels = seq(70, 0, -10)
breaks = (max(nodeHeights(phy@phylo)) - labels)
annotated = greyEpochs + geom_tree() + 
   geom_text2(aes(label = label, subset = isTip), color="#636363", size = 2.5, hjust=0, nudge_x = 1, angle = 0, family = "Helvetica", fontface = 4) +
   scale_x_continuous(breaks = breaks, labels = labels, ) + 
   geom_rootedge(rootedge = 9.8) + 
   coord_cartesian(xlim = c(-15, max(nodeHeights(phy@phylo))*1.1), ylim = c(-3, length(phy@phylo$tip.label)*1.1), expand = FALSE) +
   geom_segment(data = hpdDF, aes(x = min, y = y, xend = max, yend = y), size = 2, alpha = 0.5, color = "blue") +
   geom_point(data = ggtree(phy@phylo)$data, aes(x, y), color = node.label.cols,size = 5, shape = 21, fill = "#ffffc6", stroke = 0.5) +
   geom_text(data = ggtree(phy@phylo)$data, aes(x, y),label = nodes.text, size = 3, hjust = 0.5, vjust = 0.5) +
   theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
         axis.line.x = element_line(color = "black", size = 0.25),
         axis.ticks.x = element_line(color = "black", size = 0.25),
         axis.text.x = element_text(color = "black", size = 7),
         axis.title.x = element_text(color = "black", size = 8),
         legend.position = "none")
annotated

tipStates = data.frame(haveData$CAMphenoCat)
rownames(tipStates) = haveData$tipLabel

colorMapping = c("nonCAM"="#753071", "mCAM"="#f1a15b", "pCAM"="#7dab68")
p = gheatmap(annotated, tipStates, offset=0, width=0.01, colnames = FALSE, color="black") + 
   scale_fill_gradientn(colors=colorMapping)+
   theme(legend.position = "none")
p

ggsave("../Figures/Tips-With-Data.timetree.pdf", plot=p, device = "pdf", width=8.5, height = 11, units="in")
# The tip labels are cut off but I got around this by releasing masking layers in Illustrator.

##########################
### Threshold analysis ###
##########################
tree = phy@phylo
## MA data
maData = morphoData[!is.na(morphoData[,c("MA")]),c("CAMphenoCat","MA", "logMA")]
rownames(maData) = morphoData[!is.na(morphoData[,c("MA")]),"tipLabel"]
maTree = keep.tip(tree, tip = intersect(rownames(maData), tree$tip.label))
maData = maData[intersect(rownames(maData), tree$tip.label), ]

Xma = maData$logMA
names(Xma) = maTree$tip.label
maANC = fastAnc(maTree, maData$logMA)
ggtree(ladderize(maTree), size=2, aes(colour=c(maData$logMA, maANC)), continuous = 'color',) + geom_tiplab() + 
   scale_colour_continuous(type = "viridis", breaks=seq(from=2.7, to=4.2, length.out=6))

## IAS data
iasData = morphoData[!is.na(morphoData[,c("IAS")]),c("CAMphenoCat","IAS", "logIAS")]
rownames(iasData) = morphoData[!is.na(morphoData[,c("IAS")]),"tipLabel"]
iasTree = keep.tip(tree, tip = intersect(rownames(iasData), tree$tip.label))
iasData = iasData[intersect(rownames(iasData), tree$tip.label), ]

Xias = iasData$logIAS
names(Xias) = iasTree$tip.label
iasANC = fastAnc(iasTree, iasData$logIAS)
ggtree(ladderize(iasTree), size=2, aes(colour=c(iasData$logIAS, iasANC)), continuous = 'color',) + geom_tiplab() + 
   scale_colour_continuous(type = "viridis", breaks=seq(from=-0.4, to=-1.5, length.out=6))

## LT data
ltData = morphoData[!is.na(morphoData[,c("LT")]),c("CAMphenoCat","LT", "logLT")]
rownames(ltData) = morphoData[!is.na(morphoData[,c("LT")]),"tipLabel"]
ltTree = keep.tip(tree, tip = intersect(rownames(ltData), tree$tip.label))
ltData = ltData[intersect(rownames(ltData), tree$tip.label), ]

Xlt = ltData$logLT
names(Xlt) = ltTree$tip.label
ltANC = fastAnc(ltTree, ltData$logLT)
ggtree(ladderize(ltTree), size=2, aes(colour=c(ltData$logLT, ltANC)), continuous = 'color',) + geom_tiplab() + 
   scale_colour_continuous(type = "viridis", breaks=seq(from=2.1, to=3.6, length.out=6))

# Sample parameters
# The full analysis used 10e+06 generations, but that was run on an HPCC; using a smaller number
# for an example here.
sample = 1000
ngen = 5e+04
burnin = 0.2*ngen

## MA threshold
maThresh = threshBayes(maTree, maData[,c("CAMphenoCat","logMA")], ngen = ngen, control = list(sample = sample), plot=TRUE)
# write.csv(maThresh$par, "../Results/maThreshPar.csv")
density(maThresh)
plot(maThresh)
plot(density(maThresh))

## IAS threshold
iasThresh = threshBayes(iasTree, iasData[,c("CAMphenoCat","logIAS")], ngen = ngen, control = list(sample = sample), plot=TRUE)
# write.csv(iasThresh$par, "../Results/iasThreshPar.csv")
density(iasThresh)
plot(iasThresh)
plot(density(iasThresh))

## LT threshold
ltThresh = threshBayes(ltTree, ltData[,c("CAMphenoCat","logLT")], ngen = ngen, control = list(sample = sample, plot=TRUE))
# write.csv(ltThresh$par, "../Data/PhyloData/ltThreshPar.csv")
density(ltThresh)
plot(ltDensity)
plot(density(ltThresh))

######################################
### Ancestral State Reconstruction ###
######################################
bigPhy = read.beast("../Results/PhyloWeld.timetree.nex")
bigPhy = bigPhy@phylo
bigPhyData = read.csv("../Data/253-Taxon-CAM-States.csv")
rownames(bigPhyData) = bigPhy$tip.label

pruned = drop.tip(bigPhy, bigPhyData[bigPhyData$Drop==TRUE,"tip.label"])
prunedData = bigPhyData[bigPhyData$Drop==FALSE,]

## Plot tree circular
ggtree(pruned, layout = 'circular', ) %<+% prunedData[,c('tip.label',"CAMpheno")] + 
   geom_tippoint(aes(color=CAMpheno), show.legend = FALSE) +
   scale_color_manual(values = c("#f1a15b", "#753071", "#7dab68")) +
   geom_cladelab(label="Anacampserotaceae", node = getMRCA(pruned, c("Anacampseros_filamentosa","Talinopsis_frutescens")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Basellaceae", node = getMRCA(pruned, c("Anredera_baselloides_38","Basella_alba")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Cactaceae", node = getMRCA(pruned, c("Ariocarpus_retusus","Pereskia_sacharosa_67")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Didiereaceae", node = getMRCA(pruned, c("Alluaudia_dumosa","Ceraria_fruticulosa_72")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Montiaceae", node = getMRCA(pruned, c("Calandrinia_arenicola_90","Phemeranthus_confertiflorus_2")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Portulacaceae", node = getMRCA(pruned, c("Portulaca_sp_221","Portulaca_cf_digyna_184")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Talinaceae", node = getMRCA(pruned, c("Talinella_pachypoda_73","Talinum_arnotii_4")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1) +
   geom_cladelab(label="Molluginaceae", node = getMRCA(pruned, c("Mollugo_verticillata","Mollugo_pentaphylla")), offset.text=1.5, barsize=0.5, fontsize=1, angle='auto', offset=1)
 
######################################
### Ancestral State Reconstruction ###  
######################################
### Simmaps ###
CAMpheno = prunedData$CAMpheno
names(CAMpheno) = pruned$tip.label

# Adjust zero-edge length branches
bigdstCAM = pruned
bigdstCAM$edge.length[bigdstCAM$edge.length==0] = max(nodeHeights(pruned))*1e-6
bigTH = prunedData$CAMpheno
bigTH = setNames(bigTH, prunedData$tip.label)

## ARD
# Just doing a small number as an example; we used 10,000 simmaps in the full analysis
nsim = 10
ARDmaps = make.simmap(tree = ladderize(pruned), x = CAMpheno, model = "ARD", nsim=nsim, pi=c(0,1,0))
ARDsumm = summary(ARDmaps, plot=FALSE)

cols = setNames(c("#753071","#f1a15b","#7dab68"), c("nonCAM","mCAM","pCAM"))
plotTree(ladderize(bigdstCAM), setEnv = TRUE, fsize=0.1)
tiplabels(pie = to.matrix(CAMpheno, c("nonCAM","mCAM","pCAM"))[ladderize(bigdstCAM)$tip.label,], piecol = cols, cex = 0.1)
nodelabels(pie=ARDsumm$ace, piecol = cols[c("mCAM","nonCAM","pCAM")], cex=0.2)

## ARD, no reversions from pCAM
# I took the fitted Q matrix from the above ARD reconstruction and restricted
# reversions from pCAM
QnoREV = t(matrix(c(c(-0.003610314,0.002497845, 0.001112469),
           c(0.002401805,-0.002401805, 0.000000000),
           c(0,0.000000000,0)), ncol=3))
rownames(QnoREV) = c("mCAM", "nonCAM","pCAM")
colnames(QnoREV) = c("mCAM", "nonCAM","pCAM")
NoREVmaps = make.simmap(tree = ladderize(pruned), x = CAMpheno, Q=QnoREV, nsim=nsim, pi=c(0,1,0))
NoREVsumm = summary(NoREVmaps, plot=FALSE)

plotTree(ladderize(bigdstCAM), setEnv = TRUE, fsize=0.1)
tiplabels(pie = to.matrix(CAMpheno, c("nonCAM","mCAM","pCAM"))[ladderize(bigdstCAM)$tip.label,], piecol = CAMpheno, cex = 0.1)
nodelabels(pie=NoREVsumm$ace, piecol = cols[c("mCAM","nonCAM","pCAM")], cex=0.2)

###########################
### Phylogenetic signal ###
###########################
## MA
maPhyloData = as.matrix(maData$logMA)
rownames(maPhyloData) = rownames(maData)
maK = phylosignal(phy = maTree, x = maPhyloData, reps = 1000)
maLambda = phylosig(maTree, maPhyloData, method = "lambda", test=TRUE)
maK
maLambda
# maK = 1.149206 (p = 0.000999001); maLambda = 0.886041 (p = 3.35657e-09)

## IAS
iasPhyloData = as.matrix(iasData$logIAS)
rownames(iasPhyloData) = rownames(iasData)
iasK = phylosignal(phy = iasTree, x = iasPhyloData, reps=1000)
iasLambda = phylosig(iasTree, iasPhyloData, method = "lambda", test=TRUE)
iasK
iasLambda
# iasK = 0.5823167 (p = 0.02597403); iasLamdba = 0.557936 (p = 0.00548264)

## LT
ltPhyloData = as.matrix(ltData$logLT)
rownames(ltPhyloData) = rownames(ltData)
ltK = phylosignal(phy = ltTree, x = ltPhyloData, reps=1000)
ltLambda = phylosig(ltTree, ltPhyloData, method = "lambda", test=TRUE)
ltK
ltLambda
# ltK = 1.089129 (p = 0.000999001); ltLambda = 0.907159 (p = 0.000634405)

############
### PGLS ###
############
## IAS as a function of MA
IasMaData = morphoData[!is.na(morphoData$logMA) & !is.na(morphoData$logIAS) ,]
IasMaTree = keep.tip(tree, IasMaData$tipLabel)
pglsIasMaBM = gls(model = logIAS ~ logMA, data = IasMaData, correlation = corBrownian(phy = IasMaTree), method = "ML")
pglsIasMaOU = gls(model = logIAS ~ logMA, data = IasMaData, correlation = corPagel(1, phy = IasMaTree, fixed = FALSE), method = "ML")
summary(pglsIasMaBM)
summary(pglsIasMaOU)
plot(IasMaData$logMA, IasMaData$logIAS)
abline(a = coef(pglsIasMaBM)[1], b = coef(pglsIasMaBM)[2])
abline(a = coef(pglsIasMaOU)[1], b = coef(pglsIasMaOU)[2], lty="dashed")

## LT as a function of MA
LtMaData = morphoData[!is.na(morphoData$logMA) & !is.na(morphoData$logLT),]
LtMaTree = keep.tip(tree, LtMaData$tipLabel)
pglsLtMaBM = gls(model = logLT ~ logMA, data = LtMaData, correlation = corBrownian(phy = LtMaTree), method = "REML")
pglsLtMaOU = gls(model = logLT ~ logMA, data = LtMaData, correlation = corPagel(1, phy = LtMaTree, fixed = FALSE), method = "REML")
summary(pglsLtMaBM)
summary(pglsLtMaOU)
plot(LtMaData$logMA, LtMaData$logLT)
abline(a = coef(pglsLtMaBM)[1], b = coef(pglsLtMaBM)[2])
abline(a = coef(pglsLtMaOU)[1], b = coef(pglsLtMaOU)[2], lty="dashed")

## LT as a function of IAS
LtIasData = morphoData[!is.na(morphoData$logIAS) & !is.na(morphoData$logLT),]
rownames(LtIasData) = LtIasData$tipLabel
LtIasTree = keep.tip(tree, LtIasData$tipLabel)
pglsLtIasBM = gls(model = logLT ~ logIAS, data = LtIasData, correlation = corBrownian(phy = LtIasTree), , method = "REML")
pglsLtIasOU = gls(model = logLT ~ logIAS, data = LtIasData, correlation = corPagel(value = 1, phy = LtIasTree, fixed = FALSE), , method = "REML")
summary(pglsLtIasBM)
summary(pglsLtIasOU)
plot(LtIasData$logIAS, LtIasData$logLT)
abline(a = coef(pglsLtIasBM)[1], b = coef(pglsLtIasBM)[2])
abline(a = coef(pglsLtIasOU)[1], b = coef(pglsLtIasOU)[2], lty="dashed")

## MA as a function of CAM phenotype
pglsMaCAMbm = gls(model = logMA ~ CAMphenoCat, data = maData, correlation = corBrownian(phy = maTree), method = "ML")
pglsMaCAMou = gls(model = logMA ~ CAMphenoCat, data = maData, correlation = corPagel(1, phy = maTree, fixed = FALSE), method = "ML")
summary(pglsMaCAMbm)
summary(pglsMaCAMou)
plot(maData$CAMphenoCat, maData$logMA)
abline(a = coef(pglsMaCAMbm)[1], b = coef(pglsMaCAMbm)[2])
abline(a = coef(pglsMaCAMou)[1], b = coef(pglsMaCAMou)[2], lty="dashed")

## LT as a function of CAM phenotype
pglsLtCAMbm = gls(model = logLT ~ CAMphenoCat, data = ltData, correlation = corBrownian(phy = ltTree), method = "ML")
pglsLtCAMou = gls(model = logLT ~ CAMphenoCat, data = ltData, correlation = corPagel(1, phy = ltTree, fixed = FALSE), method = "ML")
summary(pglsLtCAMbm)
summary(pglsLtCAMou)
plot(ltData$CAMphenoCat, ltData$logLT)
abline(a = coef(pglsLtCAMbm)[1], b = coef(pglsLtCAMbm)[2])
abline(a = coef(pglsLtCAMou)[1], b = coef(pglsLtCAMou)[2], lty="dashed")

## IAS as a function of CAM phenotype
pglsIasCAMbm = gls(model = logIAS ~ CAMphenoCat, data = iasData, correlation = corBrownian(phy = iasTree), method = "ML")
pglsIasCAMou = gls(model = logIAS ~ CAMphenoCat, data = iasData, correlation = corPagel(1, phy = iasTree, fixed = FALSE), method = "ML")
summary(pglsIasCAMbm)
summary(pglsIasCAMou)
plot(iasData$CAMphenoCat, iasData$logIAS)
abline(a = coef(pglsIasCAMbm)[1], b = coef(pglsIasCAMbm)[2])
abline(a = coef(pglsIasCAMou)[1], b = coef(pglsIasCAMou)[2], lty="dashed")
