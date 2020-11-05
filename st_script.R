#################################
## Physiological comparison SC ##
#################################

library(ggplot2)
library(tidyr)
library(ggpubr)
library(corrplot)
library(factoextra)
library (vegan)
library (car)

rm(list=ls())

###DATASETS 
##raw
sf = read.table('Data/st_data.txt', sep='\t', header=T)
##results
sf.res = sf[,1:3]

###PREPARATION
## ISOTOPES
sf.res$H15N = sf$H15N
sf.res$H13C = sf$H13C
sf.res$A15N = sf$A15N
sf.res$A13C = sf$A13C
sf.res$H.C.N = sf$H.C.N
sf.res$A.C.N = sf$A.C.N

## DW, AFDW, OM AND IM
# wo:weight of filter
# w1:weight of filter + 1mL tissue
# w2:weight of filter + 1mL tissue (ashed)
# w3:weight of tube
# w4:weight of tube + 1mL tissue (bleached)
dw = sf$w1-sf$w0 # dry weight tissue in g.mL-1
afdw = sf$w1-sf$w2 # ash free dry weight tissue in g.mL-1
dwsc = sf$w4-sf$w3 # dry weight sclerite in g.mL-1
sf.res$perc_om = afdw/dw*100 # Organic Matter (%)
sf.res$perc_im = 100-sf.res$perc_om # Inorganic Matter (%)

## ALGAL CELL DENSITY
# zX_X:cell counting 5 x 1 mm2 repeated five time
# dil:dilution 
z1 = rowMeans(sf[,c('z1_1','z1_2', 'z1_3','z1_4','z1_5')])
z2 = rowMeans(sf[,c('z2_1','z2_2', 'z2_3','z2_4','z2_5')])
z3 = rowMeans(sf[,c('z3_1','z3_2', 'z3_3','z3_4','z3_5')])
z4 = rowMeans(sf[,c('z4_1','z4_2', 'z4_3','z4_4','z4_5')])
z5 = rowMeans(sf[,c('z5_1','z5_2', 'z5_3','z5_4','z5_5')])
count = (z1+z2+z3+z4+z5)*10000/(5*(1/sf$dil)) # Algal cell in 1 mL 
sf.res$cell.afdw = count/afdw # algal cell.gafdw-1 


## CHLOROPHYLLS

Ca = 11.85*(sf$A2-sf$A1)-1.54*(sf$A3-sf$A1)-0.08*(sf$A4-sf$A1)
Cb = -5.43*(sf$A2-sf$A1)+21.03*(sf$A3-sf$A1)-2.66*(sf$A4-sf$A1)
Cc = -1.67*(sf$A2-sf$A1)-7.6*(sf$A3-sf$A1)+24.52*(sf$A4-sf$A1)
sf.res$chla.afdw = Ca*sf$V_extra/afdw  # ug chl / gafdw 
sf.res$chlb.afdw = Cb*sf$V_extra/afdw  # ug chl / gafdw 
sf.res$chlc.afdw = Cc*sf$V_extra/afdw  # ug chl / gafdw
sf.res$ratio = sf.res$chlc.afdw/sf.res$chla.afdw # ratio chl c / chl a

##PROTEINS #proteins 1ML used for host by 5mL used for the zoox!!!!
#
#
abs.prot.zoox1.corr = sf$abs.prot.zoox1-sf$abs.protZ.blank # correction zoox absorbance with blank
abs.prot.zoox2.corr = sf$abs.prot.zoox2-sf$abs.protZ.blank # correction zoox absorbance with blank
abs.prot.zoox3.corr = sf$abs.prot.zoox3-sf$abs.protZ.blank # correction zoox absorbance with blank
abs.prot.host1.corr = sf$abs.prot.host1-sf$abs.protH.blank # correction host absorbance with blank
abs.prot.host2.corr = sf$abs.prot.host2-sf$abs.protH.blank # correction host absorbance with blank
abs.prot.host3.corr = sf$abs.prot.host3-sf$abs.protH.blank # correction host absorbance with blank
prot.zoox.1 = sf$STCZX4*abs.prot.zoox1.corr^4+sf$STCZX3*abs.prot.zoox1.corr^3+sf$STCZX2*abs.prot.zoox1.corr^2+sf$STCZX1*abs.prot.zoox1.corr+sf$STCZE # conc.zoox1.SDS (ug mL-1)
prot.zoox.2 = sf$STCZX4*abs.prot.zoox2.corr^4+sf$STCZX3*abs.prot.zoox2.corr^3+sf$STCZX2*abs.prot.zoox2.corr^2+sf$STCZX1*abs.prot.zoox2.corr+sf$STCZE # conc.zoox2.SDS (ug mL-1)
prot.zoox.3 = sf$STCZX4*abs.prot.zoox3.corr^4+sf$STCZX3*abs.prot.zoox3.corr^3+sf$STCZX2*abs.prot.zoox3.corr^2+sf$STCZX1*abs.prot.zoox3.corr+sf$STCZE # conc.zoox3.SDS (ug mL-1)
prot.host.1 = sf$STCHX4*abs.prot.host1.corr^4+sf$STCHX3*abs.prot.host1.corr^3+sf$STCHX2*abs.prot.host1.corr^2+sf$STCHX1*abs.prot.host1.corr+sf$STCHE # conc.host1 (ug mL-1)
prot.host.2 = sf$STCHX4*abs.prot.host2.corr^4+sf$STCHX3*abs.prot.host2.corr^3+sf$STCHX2*abs.prot.host2.corr^2+sf$STCHX1*abs.prot.host2.corr+sf$STCHE # conc.host2 (ug mL-1)
prot.host.3 = sf$STCHX4*abs.prot.host3.corr^4+sf$STCHX3*abs.prot.host3.corr^3+sf$STCHX2*abs.prot.host3.corr^2+sf$STCHX1*abs.prot.host3.corr+sf$STCHE # conc.host3 (ug mL-1)

sf.res$protZ.afdw = apply(rbind(prot.zoox.1,prot.zoox.2,prot.zoox.3),2,mean,na.rm = TRUE)*sf$SDS1V/sf$V.prot.slurry /afdw # ug prot/g afdw  
# x2 (SDS) because  previous solution is  ug/ml | vol centrifuged (5ml) | afdw (g/ml) 
sf.res$protH.afdw = apply(rbind(prot.host.1,prot.host.2,prot.host.3),2,mean,na.rm = TRUE)/afdw # ug prot / gafdw  

sf.res$prot.afdw = sf.res$protH.afdw + sf.res$protZ.afdw
# sf.res$prot_per = sf.res$prot.afdw * afdw / dw / 1000


##LIPIDS
# w0_l: dish
# wt_l: dish + total sample
# w1_l: dish + subsample
dw.lip10 = sf$w1_l - sf$w0_l
dw.lip100 = sf$wt_l-sf$w0_l
afdw.lip100 = dw.lip100*(sf.res$perc_om/100)
sf.res$lip.cont = (sf$wplate1-sf$wplate0)*(dw.lip100/(dw.lip100-dw.lip10)) # g of lipid in 100% tissue
sf.res$lip.afdw = sf.res$lip.cont/afdw.lip100 # g lipid / gafdw

### PLOT ###

## OM $ IM
#OM
perc_om.plot = sf.res[,c('species','perc_om')]
s = ggplot(data=perc_om.plot, aes(x=species, y=perc_om, fill=species))
perc_om_main_title = expression(paste("Organic Matter"))
perc_om_y_title = expression(paste("organic matter", ' [%]')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= perc_om_main_title, x="Species", y=perc_om_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))
#IM
perc_im.plot=sf.res[,c('species','perc_im')]
s = ggplot(data=perc_im.plot, aes(x=species, y=perc_im, fill=species))
perc_im_main_title = expression(paste("Inorganic Matter"))
perc_im_y_title = expression(paste("inorganic matter", ' [%]')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= perc_im_main_title, x="Species", y=perc_im_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))
#IM vs OM
compo = perc_om.plot
compo$OM=compo$perc_om
compo$perc_om=NULL
compo$IM=perc_im.plot$perc_im
compo=aggregate(.~species, compo, mean)
compo= compo %>% gather(component, percentage, OM:IM)

ggplot(compo, aes(x=species, y=percentage, fill=component))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=c('light grey','black'))+
  theme_bw()+
  labs(title= "Tissue composition", x="Species", y="Percentage (%)")+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))

## CHLOROPHYLL
# chl a
chla.plot=sf.res[,c('species','chla.afdw')]
s=ggplot(data=chla.plot, aes(x=species, y=chla.afdw, fill=species))
chla_main_title = expression(paste("Chlorophyll", italic(" a")))
chla_y_title = expression(paste("chl", italic(' a'),' [',mu,'g ', afdw^-1,']')) 
s + geom_violin(trim = FALSE)+
    geom_jitter(shape=16, position=position_jitter(0.1))+
    labs(title= chla_main_title, x="Species", y=chla_y_title)+
    theme_bw()+
    theme(legend.position='none')+
    theme(axis.text.x= element_text(face="italic"))+
    theme(plot.title= element_text(hjust="0.5"))+
    scale_y_continuous(limits=c(500,3000))

# chl c
chlc.plot=sf.res[,c('species','chlc.afdw')]
s=ggplot(data=chlc.plot, aes(x=species, y=chlc.afdw, fill=species))
chla_main_title = expression(paste("Chlorophyll", italic(" c")))
chla_y_title = expression(paste("chl", italic(' c'),' [',mu,'g ', afdw^-1,']')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= chla_main_title, x="Species", y=chla_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))+
  scale_y_continuous(limits=c(100,1000))

# ratio
ratio.plot=sf.res[,c('species','ratio')]
s=ggplot(data=ratio.plot, aes(x=species, y=ratio, fill=species))
ratio_main_title = expression(paste("Ratio chlorophylls", italic(" c/a")))
ratio_y_title = expression(paste("ratio", italic(' c/a'))) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= ratio_main_title, x="Species", y=ratio_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))+
  scale_y_continuous(limits=c(0.2,0.6))

## ALGAL CELL DENSITY
# Plot
cell.plot=sf.res[,c('species','cell.afdw')]
s=ggplot(data=cell.plot, aes(x=species, y=cell.afdw, fill=species))
cell_main_title = expression(paste("Algal symbioont content"))
cell_y_title = expression(paste("cell", ' [','count ', afdw^-1,']')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= cell_main_title, x="Species", y=cell_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))+
  scale_y_continuous(limits=c(0,1.0e+09))

## Chl a / cell
# Plot
a.cell.plot=sf.res[,c('species','cell.afdw', 'chla.afdw')]
s=ggplot(data=a.cell.plot, aes(x=species, y=chla.afdw/cell.afdw, fill=species))
a.cell_main_title = expression(paste("Chlorophyll", italic(" a"), " per cell"))
a.cell_y_title = expression(paste("chl", italic(' a'),' [',mu,'g ', cell^-1,']')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= a.cell_main_title, x="Species", y=a.cell_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))


## PROTEINS
# Plot
prot.plot=sf.res[,c('species','prot.afdw')]
s=ggplot(data=prot.plot, aes(x=species, y=prot.afdw, fill=species))
prot_main_title = expression(paste("Protein"))
prot_y_title = expression(paste("protein", ' [',mu,'g ' ,afdw^-1,']')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= prot_main_title, x="Species", y=prot_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))

## LIPIDS
# Plot
lip.plot=sf.res[,c('species','lip.afdw')]
s=ggplot(data=lip.plot, aes(x=species, y=lip.afdw, fill=species))
lip_main_title = expression(paste("Lipid"))
lip_y_title = expression(paste("lipid", ' [', 'g ' ,afdw^-1,']')) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= lip_main_title, x="Species", y=lip_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))
## ISOTOPES
Sys.setlocale("LC_ALL","English") # to successfully display "delta" symbol

isoH=sf[,c('species','H15N','H13C')]
isoA=sf[,c('species','A15N','A13C')]

carbon_title=expression(paste(delta^{13}, "C (\u2030)"))
nitrogen_title=expression(paste(delta^{15}, "N (\u2030)"))

ggscatterhist(
  isoH, x = "H13C", y = "H15N",
  color = "species", size = 3, alpha = 0.6,
  palette = c("#F1796F", "#06BCBD"),
  margin.plot = "boxplot",
  margin.params = list(fill = "species", color = "black", size = 0.2),
  title = 'Animal fraction', xlab = carbon_title, ylab = nitrogen_title,
  legend = "bottom",
  ylim=c(4.5,6.5),
  xlim=c(-23,-17),
  ggtheme = theme_bw()
  )

ggscatterhist(
  isoA, x = "A13C", y = "A15N",
  color = "species", size = 3, alpha = 0.6,
  palette = c("#F1796F", "#06BCBD"),
  margin.plot = "boxplot",
  margin.params = list(fill = "species", color = "black", size = 0.2),
  title = 'Algal fraction', xlab = carbon_title, ylab = nitrogen_title,
  legend = "bottom",
  ylim=c(4.5,6.5),
  xlim=c(-23,-17),
  ggtheme = theme_bw()
)

## Host C/N ratio
# Plot
H.CN.plot=sf[,c('species','H.C.N')]
s=ggplot(data=H.CN.plot, aes(x=species, y=H.C.N, fill=species))
H.CN_main_title = expression(paste("Host C/N ratio"))
H.CN_y_title = expression(paste("ratio")) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= H.CN_main_title, x="Species", y=H.CN_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))+
  scale_y_continuous(limits=c(3,12))

## Algal C/N ratio
# Plot
A.CN.plot=sf[,c('species','A.C.N')]
s=ggplot(data=A.CN.plot, aes(x=species, y=A.C.N, fill=species))
A.CN_main_title = expression(paste("Algal C/N ratio"))
A.CN_y_title = expression(paste("ratio")) 
s + geom_violin(trim = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.1))+
  labs(title= A.CN_main_title, x="Species", y=A.CN_y_title)+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.text.x= element_text(face="italic"))+
  theme(plot.title= element_text(hjust="0.5"))+
  scale_y_continuous(limits=c(2, 8))


## MULTIVARIATES
#selection variables
comp=sf.res[, -which(names(sf.res) %in% c("species","depth"))]
row.names(comp)=comp$sample_no
comp$sample_no=NULL
comp$lip.cont=NULL
comp$perc_im=NULL
comp$protH.afdw=NULL
comp$protZ.afdw=NULL
comp$chlb.afdw=NULL


colnames(comp)[which(names(comp) == "A13C")] = "algae d13C"
colnames(comp)[which(names(comp) == "A15N")] = "algae d15N"
colnames(comp)[which(names(comp) == "H13C")] = "host d13C"
colnames(comp)[which(names(comp) == "H15N")] = "host d15N"
colnames(comp)[which(names(comp) == "ratio")] = "ratio chl c / chl a"
colnames(comp)[which(names(comp) == "perc_om")] = "OM"
colnames(comp)[which(names(comp) == "H.C.N")] = "host C:N"
colnames(comp)[which(names(comp) == "A.C.N")] = "algae C:N"
colnames(comp)[which(names(comp) == "cell.afdw")] = "Algal symbioont content"
colnames(comp)[which(names(comp) == "chla.afdw")] = "chl a"
colnames(comp)[which(names(comp) == "chlc.afdw")] = "chl c"
colnames(comp)[which(names(comp) == "prot.afdw")] = "total lipid"
colnames(comp)[which(names(comp) == "lip.afdw")] = "total protein"


#correlation
comp.full = comp[complete.cases(comp),]
mat = cor(comp.full)
res1 = cor.mtest(comp.full, conf.level = .95) # significance test

p_value = res1$p
row.names(p_value) = row.names(mat)
colnames(p_value) = colnames(mat)
p_value[which(p_value > 0.05)] = "N.S."

#visualisation
col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mat, method = "color", col = col(200),
         type = "upper",tl.cex=0.6, order = "hclust",
         addCoef.col = "black",
         tl.col = "black", tl.srt = 90,
         p.mat=res1$p,sig.level = 0.05,
         diag = FALSE, number.cex=0.6)
#pca
var1=dim(comp.full)[2]
pca=prcomp(comp.full, scale=T)
# (scale=TRUE) : scales the variables so that all have unit variance . This is necessary if the data has different units

fviz_pca_biplot(pca,
                label='var',
                geom.ind = "point",
                repel = TRUE,
                pointshape = 21,
                pointsize = 3,
                fill.ind = na.omit(sf.res)$species,
                col.ind = "black",
                #habillage=na.omit(sf.res)$species,
                addEllipses=TRUE, ellipse.level=0.8,
                legend.title = "Species")
# PERMANOVA
Genus = rep(c("Stereonephthya", "Litophyton"), c(8, 5))
adonis(pca[["x"]] ~ Genus, permutations = 999, method = "euclidean") # p < 0.01
# dispersion test
mod = betadisper(d = dist(pca[["x"]], method = "euclidean"), group = Genus, type = "centroid")
boxplot(mod)
permutest(mod, pairwise = TRUE, permutations = 999) # ANOVA like permutation test # p < 0.001

###################
## Data Analysis ##
###################

## Mean ##
library(dplyr)
sf.res_tb = as_tibble(sf.res[, -2])
Trait_mean = sf.res_tb %>% group_by(species) %>% summarise_all (mean, na.rm = T)

## SD ##
Trait_sd = sf.res_tb %>% group_by(species) %>% summarise_all (sd, na.rm = T)

## OM *** ##
leveneTest(perc_om ~ species, sf.res) # p = 0.621
t.test(perc_om ~ species, sf.res, na.rm = T) # p = 1.373e-06

## Cell ##
leveneTest(cell.afdw ~ species, sf.res) # p = 0.2161
t.test(cell.afdw ~ species, sf.res, na.rm = T) # p = 0.07903

## Chlorophyll a ** ##
leveneTest(chla.afdw ~ species, sf.res) # p = 0.02045
t.test(chla.afdw ~ species, sf.res, na.rm = T) # p = 0.003353
 
## Chlorophyll c ##
leveneTest(chlc.afdw ~ species, sf.res) # p = 0.1271
t.test(chlc.afdw ~ species, sf.res, na.rm = T) # p = 0.5672

## Chlorophyll c/a *** ##
leveneTest(ratio ~ species, sf.res) # p = 0.2757
t.test(ratio ~ species, sf.res, na.rm = T) # p = 5.224e-08

## Chlorophyll a/ cell *** ##
leveneTest(sf.res$chla.afdw/sf.res$cell.afdw, sf.res$species) # p = 0.07025
t.test(sf.res[1:10, ]$chla.afdw/sf.res[1:10, ]$cell.afdw, sf.res[11:20, ]$chla.afdw/sf.res[11:20, ]$cell.afdw, na.rm = T)
# p = 0.0001801

## Protein ##
leveneTest(prot.afdw ~ species, sf.res) # p = 0.1995
t.test(prot.afdw ~ species, sf.res, na.rm = T) # p = 0.2195

## Lipid *** ##
leveneTest(lip.afdw ~ species, sf.res) # p = 0.3825
t.test(lip.afdw ~ species, sf.res, na.rm = T) # p = 1.371e-07

# Host d15N
leveneTest(H15N ~ species, sf.res) # p = 0.1438
t.test(H15N ~ species, sf.res, na.rm = T) # p = 2.909e-06 ***

# Host d13C
leveneTest(H13C ~ species, sf.res) # p = 0.3572
t.test(H13C ~ species, sf.res, na.rm = T) # p = 0.01995 *

# Algae d15N
leveneTest(A15N ~ species, sf.res) # p = 0.8417
t.test(A15N ~ species, sf.res, na.rm = T) # p = 0.005131 **

# Algae d13C
leveneTest(A13C ~ species, sf.res) # p = 0.002009
t.test(A13C ~ species, sf.res, na.rm = T) # p = 0.3507

# Host C/N ratio
leveneTest(H.C.N ~ species, sf.res) # p = 0.02613 *
t.test(H.C.N ~ species, sf.res, na.rm = T) # p = 0.07908

# Algae C/N ratio
leveneTest(A.C.N ~ species, sf.res) # p = 0.3416
t.test(A.C.N ~ species, sf.res, na.rm = T) # p = 0.1765
