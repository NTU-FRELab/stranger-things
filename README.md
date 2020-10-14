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
# sf.res$chlac.afdw = sf.res$chla.afdw + sf.res$chlc.afdw
# sf.res$chla.cell = sf.res$chla.afdw/sf.res$cell.afdw
```

#### Protein

```{r}

```


