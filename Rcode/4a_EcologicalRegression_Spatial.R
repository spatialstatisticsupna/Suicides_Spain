rm(list=ls())
library(bigDM)
library(dplyr)
library(INLA)
library(tidyr)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


######################################################################################
######################################################################################
##########        Ecological regression using spatial random effects        ##########
######################################################################################
######################################################################################

#####################################
## Load data and cartography files ##
#####################################
load("../Data/Suicides_Spain.Rdata")
print(Suicides)

## Set the observed and population data for the province of Madrid to NA for the years 2010–2012 ##
Suicides[Suicides$PROV=="28" & Suicides$Year %in% 2010:2012, c("O","Pop")] <- NA

## Delete from our analysis the age group "0-9" ##
Suicides <- Suicides |> filter(Age!="0-9")


#############################################################
## Compute aggregated data by sex, province and age-groups ##
#############################################################
Data <- Suicides |> 
  group_by(Sex, PROV, Age) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)

n.area <- length(unique(Data$PROV))
n.age <- length(unique(Data$Age))


##################################################################################
## Compute the percentage of population living in rural areas for each province ##
##################################################################################

## Classification of Spanish municipalities by degree of urbanisation (EUROSTAT)
## Urban=cities, and towns or suburbs / Rural=rural areas
load("../Data/DGURBA_Spain.Rdata")

RuralAreas.PROV <- DGURBA_Spain |> 
  mutate(PROV=substr(CODMUNI,1,2)) |> 
  filter(PROV %in% Carto_SpainPROV$PROV) |> 
  group_by(PROV, Tipo) |> 
  summarise(POB=sum(POB), .groups="drop") |> 
  group_by(PROV) |> 
  mutate(Rural_percent=round(100*POB/sum(POB),1)) |> 
  ungroup() |> 
  filter(Tipo=="Rural") |>
  select(PROV, Rural_percent)

carto <- merge(Carto_SpainPROV, RuralAreas.PROV)

Fig9a <- tm_shape(carto) +
  tm_polygons(fill="Rural_percent",
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_layout(panel.labels=c("Rural population (%)")) + 
  tm_options(component.autoscale = FALSE)

Data <- Data |> 
  left_join(RuralAreas.PROV, by="PROV")


###########################################
## Unemployment rate by sex and province ##
###########################################
load("../Data/Unemployment.Rdata")

carto <- Carto_SpainPROV |>
  left_join(Unemployment.PROV.S |>
              pivot_wider(names_from=Sex,
                          values_from=U.Rate,
                          names_prefix="U.Rate."),
            by="PROV")

Fig9b <- tm_shape(carto) +
  tm_polygons(fill="U.Rate.Males",
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_layout(panel.labels=c("Unemployment rate (males)")) + 
  tm_options(component.autoscale = FALSE)

Fig9c <- tm_shape(carto) +
  tm_polygons(fill="U.Rate.Females",
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_layout(panel.labels=c("Unemployment rate (females)")) + 
  tm_options(component.autoscale = FALSE)

Data <- Data |> 
  left_join(Unemployment.PROV.S, by=c("Sex","PROV"))

Fig9 <- tmap_arrange(Fig9a, Fig9b, Fig9c, nrow=1, ncol=3)
tmap_save(Fig9, filename="./Figures/Figure9.pdf", width=12, height=3.5)


################################################
## Fit age-space interaction (type III) model ##
################################################

## Hyperprior distributions ##
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


## Precision matrices ##
aux <- bigDM::connect_subgraphs(Carto_SpainPROV)
W <- aux$W
Rs <- as(Diagonal(n.area,colSums(W))-W,"Matrix")
Rs <- INLA::inla.scale.model(Rs, constr=list(A=matrix(1,1,n.area), e=0))

D <- diff(diag(n.age), differences=1)
Ra <- as(t(D)%*%D, "Matrix")
Ra <- INLA::inla.scale.model(Ra, constr=list(A=matrix(1,1,n.age), e=0))


## Data for male/female population  ##
data.M <- Data |> 
  filter(Sex=="Males") |> 
  mutate(Rural_percent=as.numeric(scale(Rural_percent)), 
         U.Rate=as.numeric(scale(U.Rate)), 
         ID.area=as.numeric(as.factor(PROV)),
         ID.age=as.numeric(as.factor(Age)),
         ID.area.age=seq(1,n.area*n.age))

data.F <- Data |> 
  filter(Sex=="Females") |> 
  mutate(Rural_percent=as.numeric(scale(Rural_percent)), 
         U.Rate=as.numeric(scale(U.Rate)),
         ID.area=as.numeric(as.factor(PROV)),
         ID.age=as.numeric(as.factor(Age)),
         ID.area.age=seq(1,n.area*n.age))


## Model without random effects ##
f.null <- O ~ 1 + Rural_percent + U.Rate

Null.M <- inla(f.null, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

Null.F <- inla(f.null, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Type III interaction model ##
R <- kronecker(Diagonal(n.area),Ra)
A.constr <- kronecker(Diagonal(n.area),matrix(1,1,n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ 1 + Rural_percent + U.Rate + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIII.M <- inla(f.TypeIII, family="poisson", data=data.M, E=Pop,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

TypeIII.F <- inla(f.TypeIII, family="poisson", data=data.F, E=Pop,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Type III interaction model (spatial+) ##
Rs.eigen <- eigen(Rs)

k <- 5
for(x in c("Rural_percent","U.Rate")){
  
  ## Male population ##
  X <- data.M |> select(PROV, x) |> distinct() |> pull(x)
  Z <- Rs.eigen$vectors[,seq(1,n.area-1-k)] %*% solve(Rs.eigen$vectors, X)[seq(1,n.area-1-k)]
  
  data.M <- data.M |> 
    left_join(tibble(PROV=unique(data.M$PROV), Z=as.numeric(Z)) |> 
                rename(!!paste0(x,".UC") := Z),
              by="PROV")
  
  ## Female population ##
  X <- data.F |> select(PROV, x) |> distinct() |> pull(x)
  Z <- Rs.eigen$vectors[,seq(1,n.area-1-k)] %*% solve(Rs.eigen$vectors, X)[seq(1,n.area-1-k)]
  
  data.F <- data.F |> 
    left_join(tibble(PROV=unique(data.M$PROV), Z=as.numeric(Z)) |> 
                rename(!!paste0(x,".UC") := Z),
              by="PROV")
}

f.spatialplus <- O ~ 1 + Rural_percent.UC + U.Rate.UC + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))
 
SpatialPlus.M <- inla(f.spatialplus, family="poisson", data=data.M, E=Pop,
                      control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

SpatialPlus.F <- inla(f.spatialplus, family="poisson", data=data.F, E=Pop,
                      control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Comparison of estimates ##
covariates <- c("Rural_percent","U.Rate")
names(covariates) <- covariates

lapply(covariates, function(x){
  
  rbind(Null=Null.M$summary.fixed[x,1:5],
        TypeIII=TypeIII.M$summary.fixed[x,1:5],
        SpatialPlus=SpatialPlus.M$summary.fixed[paste0(x,".UC"),1:5])
})

lapply(covariates, function(x){
  
  rbind(Null=Null.F$summary.fixed[x,1:5],
        TypeIII=TypeIII.F$summary.fixed[x,1:5],
        SpatialPlus=SpatialPlus.F$summary.fixed[paste0(x,".UC"),1:5])
})


###########################################################################################
## Table 4: Parameter estimates for models fitted under the simplified spatial+ approach ##
###########################################################################################
Table4 <- list(Males=round(SpatialPlus.M$summary.fixed[-1,3:5],3),
               Females=round(SpatialPlus.F$summary.fixed[-1,3:5],3))
Table4


###############################################################################
## Interpretation of fixed effect estimates as rate ratios in original scale ##
###############################################################################
aux <- Data |>
  group_by(Sex) |>
  summarise(SD.Rural_percent=sd(Rural_percent),
            SD.URate=sd(U.Rate))
print(aux)

## Auxiliary function ##
interpret_beta <- function(beta.marginal=NULL, c=NULL, sd=NULL){
  beta.transform <- inla.tmarginal(function(x) exp(c*x/sd), beta.marginal)
  res <- inla.qmarginal(c(0.5,0.025,0.975), beta.transform)
  names(res) <- c("mean","0.025quant","0.975quant")
  
  return(res)
}


## Null model 
##############
Null <- list(Males=rbind(interpret_beta(Null.M$marginals.fixed$Rural_percent, c=10, sd=aux[aux$Sex=="Males",]$SD.Rural_percent),
                         interpret_beta(Null.M$marginals.fixed$U.Rate, c=5, sd=aux[aux$Sex=="Males",]$SD.URate)),
             Females=rbind(interpret_beta(Null.F$marginals.fixed$Rural_percent, c=10, sd=aux[aux$Sex=="Females",]$SD.Rural_percent),
                           interpret_beta(Null.F$marginals.fixed$U.Rate, c=5, sd=aux[aux$Sex=="Females",]$SD.URate)))

lapply(Null, function(x){
  rownames(x) <- c("Rural_percent","U.Rate")
  round(x,3)
})


## Spatial+ model 
##################
SpatialPlus <- list(Males=rbind(interpret_beta(SpatialPlus.M$marginals.fixed$Rural_percent, c=10, sd=aux[aux$Sex=="Males",]$SD.Rural_percent),
                                interpret_beta(SpatialPlus.M$marginals.fixed$U.Rate, c=5, sd=aux[aux$Sex=="Males",]$SD.URate)),
                    Females=rbind(interpret_beta(SpatialPlus.F$marginals.fixed$Rural_percent, c=10, sd=aux[aux$Sex=="Females",]$SD.Rural_percent),
                                  interpret_beta(SpatialPlus.F$marginals.fixed$U.Rate, c=5, sd=aux[aux$Sex=="Females",]$SD.URate)))

lapply(SpatialPlus, function(x){
  rownames(x) <- c("Rural_percent","U.Rate")
  round(x,3)
})


# For males, a 10 percentage point increase in the proportion of the population living in rural areas 
# is associated with an estimated 5.4% increase in suicide mortality
# rate ratio: 1.054; 95% credible interval: [1.022, 1.086]

# Among females, no meaningful association is detected, as the rate ratio remains close to one with wide uncertainty
# rate ratio: 0.997; 95% credible interval: [0.955, 1.040]
