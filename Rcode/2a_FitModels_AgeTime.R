rm(list=ls())
library(dplyr)
library(INLA)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


################################################################################
################################################################################
##########                Age-time interaction models                 ########## 
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


##############################
## Hyperprior distributions ##
##############################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


###############################
## Define precision matrices ##
###############################
Dt1 <- diff(diag(n.year), differences=1)
Rt1 <- as(t(Dt1)%*%Dt1, "Matrix")
Rt1 <- INLA::inla.scale.model(Rt1, constr=list(A=matrix(1,1,n.year), e=0))

Dt2 <- diff(diag(n.year), differences=2)
Rt2 <- as(t(Dt2)%*%Dt2, "Matrix")
Rt2 <- INLA::inla.scale.model(Rt2, constr=list(A=rbind(matrix(1,1,n.year),matrix(1:n.year,1,n.year)), e=c(0,0)))

Da <- diff(diag(n.age), differences=1)
Ra <- as(t(Da)%*%Da, "Matrix")
Ra <- INLA::inla.scale.model(Ra, constr=list(A=matrix(1,1,n.age), e=0))


####################################################################
## Compute posterior estimates of age-group and temporal patterns ##
####################################################################
age.pattern <- inla.make.lincombs('(Intercept)'=rep(1,n.age), ID.age=diag(n.age))
names(age.pattern) <- paste("Age",as.character(1:n.age),sep=".")

time.pattern <- inla.make.lincombs('(Intercept)'=rep(1,n.year), ID.year=diag(n.year))
names(time.pattern) <- paste("Time",as.character(1:n.year),sep=".")


#############################################################
## INLA models for MALE population                         ##
## Age (RW1) + Time (RW1) + Age-Time interaction (4 types) ##
#############################################################
data.M <- Data |> 
  filter(Sex=="Males") |> 
  mutate(ID.year=as.numeric(as.factor(Year)),
         ID.age=as.numeric(as.factor(Age)),
         ID.year.age=seq(1,n.year*n.age))


## Additive model ##
f.Additive <- O ~ 1 + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif)))
  
Additive <- inla(f.Additive, family="poisson", data=data.M, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                 lincomb=c(age.pattern, time.pattern))

## Type I interaction ##
f.TypeI <- O ~ 1 + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
TypeI <- inla(f.TypeI, family="poisson", data=data.M, E=Pop,
              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
              lincomb=c(age.pattern, time.pattern))

## Type II interaction ##
R <- kronecker(Rt1,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")
  
f.TypeII <- O ~ 1 + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))
  
TypeII <- inla(f.TypeII, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, time.pattern))

## Type III interaction ##
R <- kronecker(Diagonal(n.year),Ra)
A.constr <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(A.constr[-1,],"matrix")
  
f.TypeIII <- O ~ 1 + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))
  
TypeIII <- inla(f.TypeIII, family="poisson", data=data.M, E=Pop,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                lincomb=c(age.pattern, time.pattern))

## Type IV interaction ##
R <- kronecker(Rt1,Ra)
A1 <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")
  
f.TypeIV <- O ~ 1 + 
  f(ID.year, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year+n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIV <- inla(f.TypeIV, family="poisson", data=data.M, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, time.pattern))


MODELS.M <- list('Additive'=Additive, 'TypeI'=TypeI, 'TypeII'=TypeII, 'TypeIII'=TypeIII, 'TypeIV'=TypeIV)


#############################################################
## INLA models for FEMALE population                       ##
## Age (RW1) + Time (RW2) + Age-Time interaction (4 types) ##
#############################################################
data.F <- Data |> 
  filter(Sex=="Females") |> 
  mutate(ID.year=as.numeric(as.factor(Year)),
         ID.age=as.numeric(as.factor(Age)),
         ID.year.age=seq(1,n.year*n.age))

## Additive model ##
f.Additive <- O ~ 1 + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif)))

Additive <- inla(f.Additive, family="poisson", data=data.F, E=Pop,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                 lincomb=c(age.pattern, time.pattern))

## Type I interaction ##
f.TypeI <- O ~ 1 + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=matrix(rep(1:n.age,each=T),1,n.year*n.age),e=0))

TypeI <- inla(f.TypeI, family="poisson", data=data.F, E=Pop,
              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
              lincomb=c(age.pattern, time.pattern))

## Type II interaction ##
R <- kronecker(Rt2,Diagonal(n.age))
A.constr <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeII <- O ~ 1 + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=2*n.age, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeII <- inla(f.TypeII, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, time.pattern))

## Type III interaction ##
R <- kronecker(Diagonal(n.year),Ra)
A.constr <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ 1 + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIII <- inla(f.TypeIII, family="poisson", data=data.F, E=Pop,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
                lincomb=c(age.pattern, time.pattern))

## Type IV interaction ##
R <- kronecker(Rt2,Ra)
A1 <- kronecker(matrix(1,1,n.year),Diagonal(n.age))
A2 <- kronecker(Diagonal(n.year),matrix(1,1,n.age))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ 1 + 
  f(ID.year, model="rw2", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.age, model="rw1", constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year.age, model='generic0', Cmatrix=R, rankdef=n.year+2*n.age-1, constr=TRUE,
    hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

TypeIV <- inla(f.TypeIV, family="poisson", data=data.F, E=Pop,
               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE),
               lincomb=c(age.pattern, time.pattern))


MODELS.F <- list('Additive'=Additive, 'TypeI'=TypeI, 'TypeII'=TypeII, 'TypeIII'=TypeIII, 'TypeIV'=TypeIV)


##################
## Save results ##
##################
save(list=c("MODELS.M","MODELS.F"), file="INLAmodels_AgeTime.Rdata")
