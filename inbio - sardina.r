# inbio: estimation of biological parameters
# (Software developed by P Sampedro, V Trujillo, M Sainza) 
# J Cebrian, G Costas, E Velasco
# May 2018
#
# Versions used R 3.0.0
# -------------

# Preparing the appropiate file csv
# Requirements: International Metric System (decimals with dots, not separation in thousands)
# Files must have the following columns: 
# tal: length in cm or mm
# pes: weight in g
# ed: age (0, 1, 2, ...) >> this column is optional
# mad: maturity stage (0=inmature, 1=mature)
# sex: male (1), female (2), indeterminate (3)

# Each raw in the file represent a fish. 
# Empty cell for missing values (not -1, not 0)
# The name of the file will be used during the inbio importation process
# Your file must be separated by comma (,)

library (dplyr)

# Setting working directory
setwd("C:/Users/Eva/Dropbox/Misión Turkey/inbio")

# read data and prepare data set
data <- read.csv2("C:/Users/Eva/Dropbox/Misión Turkey/inbio/biology_market.csv"); head(data) 
biol_sardine <- select(data,TALLA, PESO, MADURACIÓN, SEXO, EDAD); head (biol_sardine)

str(biol_sardine)
table(biol_sardine$SEXO)
table(biol_sardine$MADURACIÓN)
table(biol_sardine$EDAD) # in this case, we have few age data, so we'll calculate by length

biol_sardine<-biol_sardine%>%rename(tal=TALLA, pes=PESO, mad=MADURACIÓN, sex=SEXO, ed=EDAD)%>%
  mutate(mad=as.factor(as.character(mad)), pes=as.numeric(as.character(pes)), sex = recode(sex, "H"="2", "M"="1", "U"="3"),mad = recode(mad, "1"="0", "2"="1", "3"="1","4"="1","5"="1","6"="1"), tal=tal/10)

summary(biol_sardine)

write.table(biol_sardine, "biol_sardine.csv", row.names = FALSE, sep=",")


# Results will be exported in files 'txt' and 'pdf'

# 4 routines Inbio:
# Peso = weight
# Madurez = maturity
# Crecimiento = growth
# Sex ratio

# Install package inbio 2.1 from local zip
library(inbio)

#######################################################################
# Routine Weight: calculate length-weight relationship (nor linear estimation, minimum squares fitting)
# Method bootstrap for calculate variation coefficients
# peso(especie = "spp name", cl = 0.5, unid="cm", a=5e-4, b=3, sex = F, n = 1000)

# Arguments: 
# especie: file's name with data (without extension, without directory)
# cl: interval of Length (0.5-1). Default=1
# unid: cm
# a: start value of the fitting (should be closed to the expected value)
# b: start value of the fitting (should be closed to the expected value)
# sex: if sex=F (by default), proccess analyze males and females together,
# if sex= T, analyze separately. 
# n: number of bootstrap replications. Recommended 1000

# Results: 'txt' file and 'pdf' file in the working directory 

# Examples:

# calculate L-W relationship for males and females all together
peso(especie = "biol_sardine", cl = 0.5, unid="cm", a=5e-4, b=3, sex = F, n = 100)

# calculate L-W relationship for males and females separately
# peso(especie="biol_sardine", sex=T) 
# Calculate L-W relationship of a specie measured nearest 0.5 cm, giving starting values
# peso(especie="biol_sardine", cl=0.5, a=0.04, b=2.8)

# When action is done, it appears the next message 
# "EL PROGRAMA HA TERMINADO.PUEDES VER LOS RESULTADOS EN EL DIRECTORIO DE TRABAJO
# That means the process has end and you can see results saved in the chosen working directory
# The proccess takes much or less time depending on your PC's memmory RAM, the volume of data and the self routine.

# You can run another routine or change data file. Pay attention because running two routines, overwrite the results! 

#######################################################################
# Routine Maturity: calculate size and age of sexual maturity (glm)
# Method bootstrap

biol_sardine_ogive<-filter(biol_sardine,!is.na(mad)) # if there are problems with NAs
write.table(biol_sardine_ogive, "biol_sardine.csv" ,row.names = FALSE,sep=",")

# madurez(especie = "spp name", cl = 0.5, unid="cm", edad=T, sex = F, n = 100)

# Arguments: 
# especie: file's name with data (without extension, without directory)
# cl: interval of Length (0.5-1). Default=1
# unid: cm
# edad: age. If edad=T (by default), age data is required, age of maturity is calculated. 
# If edad=F, age of maturity is not calculated. 
# sex: if sex=F (by default), proccess analyze males and females together,
# if sex= T, analyze separately. 
# n: number of bootstrap replications. Recommended 1000

# Results: 'txt' file and 'pdf' file in the working directory 

# Examples: 

# Calculate size and age at sexual maturity and its variation coefficients for both sexes (age required).
# madurez(especie="biol_sardine")

# If you don't have age data:
# Calculate size at sexual maturity and its variation coefficient for males and females separately
madurez(especie="biol_sardine", edad=F, sex=T, n=100)

#######################################################################
# Routine Growth: estimation of weights and lengths fitting the von Bertalanffy growth curve
# Lt=Linf(1-e^(-k(t-t0))
# Wt=Winf(1-e^(-k(t-t0)))^b

# crecimiento(especie = "spp name", cl = 0.5, unid="cm", sex = F, b = 3, Li = 28, Lfija=F, Ki = 0.7, T0i = -1, Wi = 2500, Wfijo=F, Kwi = 0.4, T0wi = -1, n = 1000)

# Arguments: 
# especie: file's name with data (without extension, without directory)
# cl: interval of Length (0.5-1). Default=1
# unid: cm
# sex: if sex=F (by default), proccess analyze males and females together,
# if sex= T, analyze separately. 
# b: L-W relationship (default 3)
# Li: start value of Linf
# Lfija: if Lfija=T, Li=Linf. If Lfija=F (default) Linf is calculated and its coefficient of variation
# Ki: start value of Ki
# T0i: start value of T0i (when fish Length = 0)
# Wi: start value of Winf
# Wfijo: if Wfijo=T, Wi=Winf. If Wfija=F (default) Winf is calculated and its coefficient of variation
# Kwi: start value of Kwi
# T0wi: start value of T0wi (when fish weight = 0)
# n: number of bootstrap replications. Recommended 1000

# Results: 'txt' file and 'pdf' file in the working directory 
# Examples: 
# Calculate the parameters of the von Bertalanffy growth curve for both sexes together
crecimiento(especie="biol_sardine", Li = 28, Ki = 0.27, T0i = -1.9, Wi = 2500, Kwi = 0.4, T0wi = -1)

# Calculate the parameters of the von Bertalanffy growth curve for both sexes separately. 
# Fix Linf=20 cm and Winf=100g
crecimiento(especie="biol_sardine", sex=T, Li = 28, Lfija=T, Ki = 0.27, T0i = -1.9, Wi = 100, Wfijo=T, Kwi = 0.4, T0wi = -1)

#######################################################################
# Routine Sex Ratio: calculate the porcentage of males and females by size and age
# indeterminated are not considered.

# sexratio(especie = "spp_name", cl = 0.5, unid="cm", edad=T, n = 1000)
# Arguments: 
# especie: file's name with data (without extension, without directory)
# cl: interval of Length (0.5-1). Default=1
# unid: cm
# edad: age. If edad=T (by default), age data is required, age of maturity is calculated. 
# n: number of bootstrap replications. Recommended 1000

# Results: 'txt' file and 'pdf' file in the working directory 
# Examples: 
# Calculate sexratio by lenth interval of 1cm and sex ratio by age (Age required).
# sexratio(especie="biol_sardine")

# Calculate sexratio by lenth interval of 0.5cm, not by age.
sexratio(especie="biol_sardine", cl=0.5, edad=F)
