rm(list=ls())
library(bigDM)
library(dplyr)
library(INLA)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


################################################################################
################################################################################
##########                Age-space interaction models                 #########
################################################################################
################################################################################

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


##############################
## Hyperprior distributions ##
##############################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


###############################
## Define precision matrices ##
###############################
aux <- bigDM::connect_subgraphs(Carto_SpainPROV)
W <- aux$W
Rs <- as(Diagonal(n.area,colSums(W))-W,"Matrix")
Rs <- INLA::inla.scale.model(Rs, constr=list(A=matrix(1,1,n.area), e=0))

D <- diff(diag(n.age), differences=1)
Ra <- as(t(D)%*%D, "Matrix")
Ra <- INLA::inla.scale.model(Ra, constr=list(A=matrix(1,1,n.age), e=0))


###################################################################
## Compute posterior estimates of age-group and spatial patterns ##
###################################################################
age.pattern <- inla.make.lincombs('(Intercept)'=rep(1,n.age), ID.age=diag(n.age))
names(age.pattern) <- paste("Age",as.character(1:n.age),sep=".")

spatial.pattern <- inla.make.lincombs('(Intercept)'=rep(1,n.area), ID.area=diag(n.area))
names(spatial.pattern) <- paste("Spatial",as.character(1:n.area),sep=".")


################################################################
## INLA models for MALE population                            ##
## Age (RW1) + Space (iCAR) + Age-Space interaction (4 types) ##
################################################################
data.M <- Data |> 
  filter(Sex=="Males") |> 
  mutate(ID.area=as.numeric(as.factor(PROV)),
         ID.age=as.numeric(as.factor(Age)),
         ID.area.age=seq(1,n.area*n.age))


## Additive model ##
f.Additive <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif)))
  
Additive <- inla(f.Additive, family="poisson", data=data.M, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                 lincomb=c(age.pattern, spatial.pattern))
  
## Type I interaction ##
f.TypeI <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
TypeI <- inla(f.TypeI, family="poisson", data=data.M, E=Pop,
              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
              lincomb=c(age.pattern, spatial.pattern))
  
## Type II interaction ##
R <- kronecker(Rs,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.area),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")
  
f.TypeII <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeII <- inla(f.TypeII, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, spatial.pattern))

## Type III interaction ##
R <- kronecker(Diagonal(n.area),Ra)
A.constr <- kronecker(Diagonal(n.area),matrix(1,1,n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIII <- inla(f.TypeIII, family="poisson", data=data.M, E=Pop,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                lincomb=c(age.pattern, spatial.pattern))

## Type IV interaction ##
R <- kronecker(Rs,Ra)
A1 <- kronecker(matrix(1,1,n.area),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.area),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area+n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIV <- inla(f.TypeIV, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, spatial.pattern))


MODELS.M <- list('Additive'=Additive, 'TypeI'=TypeI, 'TypeII'=TypeII, 'TypeIII'=TypeIII, 'TypeIV'=TypeIV)


################################################################
## INLA models for FEMALE population                          ##
## Age (RW1) + Space (iCAR) + Age-Space interaction (4 types) ##
################################################################
data.F <- Data |> 
  filter(Sex=="Females") |> 
  mutate(ID.area=as.numeric(as.factor(PROV)),
         ID.age=as.numeric(as.factor(Age)),
         ID.area.age=seq(1,n.area*n.age))


## Additive model ##
f.Additive <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif)))

Additive <- inla(f.Additive, family="poisson", data=data.F, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                 lincomb=c(age.pattern, spatial.pattern))

## Type I interaction ##
f.TypeI <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))

TypeI <- inla(f.TypeI, family="poisson", data=data.F, E=Pop,
              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
              lincomb=c(age.pattern, spatial.pattern))

## Type II interaction ##
R <- kronecker(Rs,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.area),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeII <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeII <- inla(f.TypeII, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, spatial.pattern))

## Type III interaction ##
R <- kronecker(Diagonal(n.area),Ra)
A.constr <- kronecker(Diagonal(n.area),matrix(1,1,n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIII <- inla(f.TypeIII, family="poisson", data=data.F, E=Pop,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                lincomb=c(age.pattern, spatial.pattern))

## Type IV interaction ##
R <- kronecker(Rs,Ra)
A1 <- kronecker(matrix(1,1,n.area),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.area),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ 1 + 
  f(ID.area, model='besag', graph=Rs, scale.model=TRUE, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.age, model='generic0', Cmatrix=R, rankdef=n.area+n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIV <- inla(f.TypeIV, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, spatial.pattern))


MODELS.F <- list('Additive'=Additive, 'TypeI'=TypeI, 'TypeII'=TypeII, 'TypeIII'=TypeIII, 'TypeIV'=TypeIV)


##################
## Save results ##
##################
save(list=c("MODELS.M","MODELS.F"), file="INLAmodels_AgeSpace.Rdata")
