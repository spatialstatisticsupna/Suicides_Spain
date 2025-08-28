rm(list=ls())
library(bigDM)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(INLA)
library(tidyr)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#######################################################################################
#######################################################################################
##########        Ecological regression using temporal random effects        ##########
#######################################################################################
#######################################################################################

#####################################
## Load data and cartography files ##
#####################################
load("../Data/Suicides_Spain.Rdata")
print(Suicides)

## Set the observed and population data for the province of Madrid to NA for the years 2010–2012 ##
Suicides[Suicides$PROV=="28" & Suicides$Year %in% 2010:2012, c("O","Pop")] <- NA

## Delete from our analysis the age group "0-9" ##
Suicides <- Suicides |> filter(Age!="0-9")


#########################################################
## Compute aggregated data by sex, year and age-groups ##
#########################################################
Data <- Suicides |> 
  group_by(Sex, Year, Age) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)

n.year <- length(unique(Data$Year))
n.age <- length(unique(Data$Age))


#######################################
## Unemployment rate by sex and year ##
#######################################
load("../Data/Unemployment.Rdata")

aux <- Unemployment.PROV.ST |> 
  group_by(Sex, Year) |> 
  summarise(U.Rate=mean(U.Rate), .groups="drop")

Fig10a <- ggplot(aux, aes(x=Year, y=U.Rate, color=Sex, group=Sex)) +
  geom_line(linewidth=0.8) +
  xlab("Year") + 
  ylab("Unemployment rate (%)") +
  scale_color_manual(values=c("Males"="#4f81bd", "Females"="#d35400"), name=NULL) +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="top",
        legend.text=element_text(size=14),
        axis.title.y.right=element_text(angle=90))

Data <- left_join(Data, aux, by=c("Sex","Year"))


##################
## Poverty rate ##
##################
load("../Data/Poverty.Rdata")

aux <- Poverty.TA |> 
  mutate(Year=as.character(Year)) |> 
  group_by(Sex, Year) |> 
  summarise(P.Rate=mean(Rate), .groups="drop")

Fig10b <- ggplot(aux, aes(x=Year, y=P.Rate, color=Sex, group=Sex)) +
  geom_line(linewidth=0.8) +
  xlab("Year") + 
  ylab("At risk of poverty rate (%)") +
  scale_color_manual(values=c("Males"="#4f81bd", "Females"="#d35400"), name=NULL) +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="top",
        legend.text=element_text(size=14),
        axis.title.y.right=element_text(angle=90))

Data <- left_join(Data, aux, by=c("Sex","Year"))

Fig10 <- ggarrange(Fig10a, Fig10b, nrow=1, ncol=2)
ggsave("./Figures/Figure10.pdf", Fig10, width=12, height=6)


#################################################
## Fit age-time interaction (type II/IV) model ##
#################################################

## Hyperprior distributions ##
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


## Precision matrices ##
Dt1 <- diff(diag(n.year), differences=1)
Rt1 <- as(t(Dt1)%*%Dt1, "Matrix")
Rt1 <- INLA::inla.scale.model(Rt1, constr=list(A=matrix(1,1,n.year), e=0))

Dt2 <- diff(diag(n.year), differences=2)
Rt2 <- as(t(Dt2)%*%Dt2, "Matrix")
Rt2 <- INLA::inla.scale.model(Rt2, constr=list(A=rbind(matrix(1,1,n.year),matrix(1:n.year,1,n.year)), e=c(0,0)))

Da <- diff(diag(n.age), differences=1)
Ra <- as(t(Da)%*%Da, "Matrix")
Ra <- INLA::inla.scale.model(Ra, constr=list(A=matrix(1,1,n.age), e=0))


## Data for male/female population  ##
data.M <- Data |> 
  filter(Sex=="Males") |> 
  mutate(U.Rate=as.numeric(scale(U.Rate)), 
         P.Rate=as.numeric(scale(P.Rate)), 
         ID.year=as.numeric(as.factor(Year)),
         ID.age=as.numeric(as.factor(Age)),
         ID.year.age=seq(1,n.year*n.age))

data.F <- Data |> 
  filter(Sex=="Females") |> 
  mutate(U.Rate=as.numeric(scale(U.Rate)), 
         P.Rate=as.numeric(scale(P.Rate)), 
         ID.year=as.numeric(as.factor(Year)),
         ID.age=as.numeric(as.factor(Age)),
         ID.year.age=seq(1,n.year*n.age))


## Model without random effects ##
f.null <- O ~ 1 + U.Rate + P.Rate

Null.M <- inla(f.null, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

Null.F <- inla(f.null, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Type II/IV interaction model ##
R <- kronecker(Rt1,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeII <- O ~ 1 + U.Rate + P.Rate + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeII.M <- inla(f.TypeII, family="poisson", data=data.M, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

R <- kronecker(Rt2,Ra)
A1 <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ 1 + U.Rate + P.Rate + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year+2*n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIV.F <- inla(f.TypeIV, family="poisson", data=data.F, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Type II/IV interaction model (spatial+) ##
Rt1.eigen <- eigen(Rt1)
Rt2.eigen <- eigen(Rt2)

k <- 1
for(x in c("U.Rate","P.Rate")){
  
  ## Male population ##
  X <- data.M |> select(Year, x) |> distinct() |> pull(x)
  Z <- Rt1.eigen$vectors[,seq(1,n.year-1-k)] %*% solve(Rt1.eigen$vectors, X)[seq(1,n.year-1-k)]
  
  data.M <- data.M |> 
    left_join(tibble(Year=unique(data.M$Year), Z=as.numeric(Z)) |> 
                rename(!!paste0(x,".UC") := Z),
              by="Year")
  
  ## Female population ##
  X <- data.F |> select(Year, x) |> distinct() |> pull(x)
  Z <- Rt2.eigen$vectors[,seq(1,n.year-2-k)] %*% solve(Rt2.eigen$vectors, X)[seq(1,n.year-2-k)]

  data.F <- data.F |> 
    left_join(tibble(Year=unique(data.F$Year), Z=as.numeric(Z)) |> 
                rename(!!paste0(x,".UC") := Z),
              by="Year")
}

R <- kronecker(Rt1,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.spatialplus <- O ~ 1 + U.Rate.UC + P.Rate.UC + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

SpatialPlus.M <- inla(f.spatialplus, family="poisson", data=data.M, E=Pop,
                      control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

R <- kronecker(Rt2,Ra)
A1 <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.spatialplus <- O ~ 1 + U.Rate.UC + P.Rate.UC + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year+2*n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

SpatialPlus.F <- inla(f.spatialplus, family="poisson", data=data.F, E=Pop,
                      control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))


## Comparison of estimates ##
covariates <- c("U.Rate","P.Rate")
names(covariates) <- covariates

lapply(covariates, function(x){
  
  rbind(Null=Null.M$summary.fixed[x,1:5],
        TypeII=TypeII.M$summary.fixed[x,1:5],
        SpatialPlus=SpatialPlus.M$summary.fixed[paste0(x,".UC"),1:5])
})

lapply(covariates, function(x){
  
  rbind(Null=Null.F$summary.fixed[x,1:5],
        TypeIV=TypeIV.F$summary.fixed[x,1:5],
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
  summarise(SD.URate=sd(U.Rate),
            SD.PRate=sd(P.Rate))

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
Null <- list(Males=rbind(interpret_beta(Null.M$marginals.fixed$U.Rate, c=1, sd=aux[aux$Sex=="Males",]$SD.URate),
                         interpret_beta(Null.M$marginals.fixed$P.Rate, c=1, sd=aux[aux$Sex=="Males",]$SD.PRate)),
             Females=rbind(interpret_beta(Null.F$marginals.fixed$U.Rate, c=1, sd=aux[aux$Sex=="Females",]$SD.URate),
                           interpret_beta(Null.F$marginals.fixed$P.Rate, c=1, sd=aux[aux$Sex=="Females",]$SD.PRate)))

lapply(Null, function(x){
  rownames(x) <- c("U.Rate","P.Rate")
  round(x,3)
})


## Spatial+ model 
##################
SpatialPlus <- list(Males=rbind(interpret_beta(SpatialPlus.M$marginals.fixed$U.Rate, c=1, sd=aux[aux$Sex=="Males",]$SD.URate),
                                interpret_beta(SpatialPlus.M$marginals.fixed$P.Rate, c=1, sd=aux[aux$Sex=="Males",]$SD.PRate)),
                    Females=rbind(interpret_beta(SpatialPlus.F$marginals.fixed$U.Rate, c=1, sd=aux[aux$Sex=="Females",]$SD.URate),
                                  interpret_beta(SpatialPlus.F$marginals.fixed$P.Rate, c=1, sd=aux[aux$Sex=="Females",]$SD.PRate)))

lapply(SpatialPlus, function(x){
  rownames(x) <- c("U.Rate","P.Rate")
  round(x,3)
})


SpatialPlus <- list(Males=rbind(exp(1/aux[aux$Sex=="Males",]$SD.URate*SpatialPlus.M$summary.fixed["U.Rate",c("mean","0.025quant","0.975quant")]),
                                exp(1/aux[aux$Sex=="Males",]$SD.PRate*SpatialPlus.M$summary.fixed["P.Rate",c("mean","0.025quant","0.975quant")])),
                    Females=rbind(exp(1/aux[aux$Sex=="Females",]$SD.URate*SpatialPlus.F$summary.fixed["U.Rate",c("mean","0.025quant","0.975quant")]),
                                  exp(1/aux[aux$Sex=="Females",]$SD.PRate*SpatialPlus.F$summary.fixed["P.Rate",c("mean","0.025quant","0.975quant")])))

lapply(SpatialPlus, function(x) round(x,3))

# A one percentage point increase in annual unemployment rate is associated with a 2.4% increase in female suicide mortality
# rate ratio: 1.024; 95% credible interval: [1.004, 1.045]

# A weaker and more uncertain effect is estimated for males
# rate ratio: 1.015; 95% credible interval: [1.000, 1.031]
