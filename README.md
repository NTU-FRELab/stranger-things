# stranger-things

## Data and script for data analysis of Hsu et al. Stranger things: physiological performance of two octocorals associated with singular Symbioniaceae in a marginal coral community from northern Taiwan - [FRElab](https://www.dipintothereef.com/) 

### Run all the packages needed
```{r}
 
library(ggplot2)
library(tidyr)
library(ggpubr)
library(corrplot)
library(factoextra)
library(vegan)
library(car)
 
```

### Use `rm(list=ls())` to clear the environment

### Datasets
```{r}
# raw data
sf = read.table('Data/st_data.txt', sep='\t', header=T)

# creating a dataframe for the coming results
sf.res = sf[,1:3]
```

### Traits

#### Isotope
```{r}
# H: coral host
# A: algae 
# C.N: carbon to nitrogen ratio
sf.res$H15N = sf$H15N
sf.res$H13C = sf$H13C
sf.res$A15N = sf$A15N
sf.res$A13C = sf$A13C
sf.res$H.C.N = sf$H.C.N
sf.res$A.C.N = sf$A.C.N
```

#### DW, AFDW, OM and IM

```{r}
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
```
#### Alagal cell density

```{r}
# zX_X:cell counting 5 x 1 mm2 repeated five time
# dil:dilution
z1 = rowMeans(sf[,c('z1_1','z1_2', 'z1_3','z1_4','z1_5')])
z2 = rowMeans(sf[,c('z2_1','z2_2', 'z2_3','z2_4','z2_5')])
z3 = rowMeans(sf[,c('z3_1','z3_2', 'z3_3','z3_4','z3_5')])
z4 = rowMeans(sf[,c('z4_1','z4_2', 'z4_3','z4_4','z4_5')])
z5 = rowMeans(sf[,c('z5_1','z5_2', 'z5_3','z5_4','z5_5')])
count = (z1+z2+z3+z4+z5)*10000/(5*(1/sf$dil)) # Algal cell in 1 mL 
sf.res$cell.afdw = count/afdw # algal cell.gafdw-1 
```
#### Chlorophyll

```{r}
Ca = 11.85*(sf$A2-sf$A1)-1.54*(sf$A3-sf$A1)-0.08*(sf$A4-sf$A1) # chl a
Cb = -5.43*(sf$A2-sf$A1)+21.03*(sf$A3-sf$A1)-2.66*(sf$A4-sf$A1) # chl b
Cc = -1.67*(sf$A2-sf$A1)-7.6*(sf$A3-sf$A1)+24.52*(sf$A4-sf$A1) # chl c
# V_extra: volume of extraction solution (mL)
sf.res$chla.afdw = Ca*sf$V_extra/afdw  # ug chl / gafdw 
sf.res$chlb.afdw = Cb*sf$V_extra/afdw  # ug chl / gafdw 
sf.res$chlc.afdw = Cc*sf$V_extra/afdw  # ug chl / gafdw
```

#### Protein

```{r}
# correction zoox (algae) absorbance with blank
abs.prot.zoox1.corr = sf$abs.prot.zoox1-sf$abs.protZ.blank 
abs.prot.zoox2.corr = sf$abs.prot.zoox2-sf$abs.protZ.blank
abs.prot.zoox3.corr = sf$abs.prot.zoox3-sf$abs.protZ.blank
abs.prot.host1.corr = sf$abs.prot.host1-sf$abs.protH.blank
abs.prot.host2.corr = sf$abs.prot.host2-sf$abs.protH.blank
abs.prot.host3.corr = sf$abs.prot.host3-sf$abs.protH.blank

# substitute absorbance into trichromatic equations (Jeffrey and Humphrey, 1975)
prot.zoox.1 = sf$STCZX4*abs.prot.zoox1.corr^4+sf$STCZX3*abs.prot.zoox1.corr^3+sf$STCZX2*abs.prot.zoox1.corr^2+sf$STCZX1*abs.prot.zoox1.corr+sf$STCZE # conc.zoox1.SDS (ug mL-1)
prot.zoox.2 = sf$STCZX4*abs.prot.zoox2.corr^4+sf$STCZX3*abs.prot.zoox2.corr^3+sf$STCZX2*abs.prot.zoox2.corr^2+sf$STCZX1*abs.prot.zoox2.corr+sf$STCZE # conc.zoox2.SDS (ug mL-1)
prot.zoox.3 = sf$STCZX4*abs.prot.zoox3.corr^4+sf$STCZX3*abs.prot.zoox3.corr^3+sf$STCZX2*abs.prot.zoox3.corr^2+sf$STCZX1*abs.prot.zoox3.corr+sf$STCZE # conc.zoox3.SDS (ug mL-1)
prot.host.1 = sf$STCHX4*abs.prot.host1.corr^4+sf$STCHX3*abs.prot.host1.corr^3+sf$STCHX2*abs.prot.host1.corr^2+sf$STCHX1*abs.prot.host1.corr+sf$STCHE # conc.host1 (ug mL-1)
prot.host.2 = sf$STCHX4*abs.prot.host2.corr^4+sf$STCHX3*abs.prot.host2.corr^3+sf$STCHX2*abs.prot.host2.corr^2+sf$STCHX1*abs.prot.host2.corr+sf$STCHE # conc.host2 (ug mL-1)
prot.host.3 = sf$STCHX4*abs.prot.host3.corr^4+sf$STCHX3*abs.prot.host3.corr^3+sf$STCHX2*abs.prot.host3.corr^2+sf$STCHX1*abs.prot.host3.corr+sf$STCHE # conc.host3 (ug mL-1)

# initial concentration of algal protein
sf.res$protZ.afdw = apply(rbind(prot.zoox.1,prot.zoox.2,prot.zoox.3),2,mean,na.rm = TRUE)*sf$SDS1V/sf$V.prot.slurry /afdw # ug prot/g afdw  
# SDS1V = 2 mL SDS (added to release of algal protein)
# V.prot.slurry = 5 mL (original volume being centrifuged) 

# initial concentration of host protein
sf.res$protH.afdw = apply(rbind(prot.host.1,prot.host.2,prot.host.3),2,mean,na.rm = TRUE)/afdw # ug prot / gafdw
# no dilution

# total protein
sf.res$prot.afdw = sf.res$protH.afdw + sf.res$protZ.afdw
```

#### Lipid

```{r}
# w0_l: dish
# wt_l: dish + total sample
# w1_l: dish + subsample
# wplate0: dish
# wplate1: dish + lipid
dw.lip10 = sf$w1_l - sf$w0_l
dw.lip100 = sf$wt_l-sf$w0_l
afdw.lip100 = dw.lip100*(sf.res$perc_om/100) # using OM percentage to convert dw into afdw

# total lipid amounts
sf.res$lip.cont = (sf$wplate1-sf$wplate0)*(dw.lip100/(dw.lip100-dw.lip10)) # g of lipid in 100% tissue

# total lipid concentration
sf.res$lip.afdw = sf.res$lip.cont/afdw.lip100 # g lipid / gafdw
```

### Plot

#### OM
```{r}
perc_om.plot = sf.res[,c('species','perc_om')]
ggplot(data = perc_om.plot, aes(x= species, y= perc_om, fill=species))
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
```
#### IM
```{r}
perc_im.plot=sf.res[,c('species','perc_im')]
ggplot(data = perc_im.plot, aes(x= species, y= perc_im, fill=species))
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
```