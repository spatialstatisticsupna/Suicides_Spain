rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(INLA)
library(tidyr)
library(tmap)

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


########################################
## Load previously fitted INLA models ##
########################################
load("INLAmodels_AgeSpace.Rdata")

## Model comparison ##
aux1 <- lapply(MODELS.M, function(x){
  data.frame(DIC=x$dic$dic,
             WAIC=x$waic$waic)
})

aux2 <- lapply(MODELS.F, function(x){
  data.frame(DIC=x$dic$dic,
             WAIC=x$waic$waic)
})

Table3 <- list(Males=do.call(rbind,aux1),
               Females=do.call(rbind,aux2))
print(Table3, digits=5)


## Select final model ##
FinalModel.M <- MODELS.M$TypeIII
FinalModel.F <- MODELS.F$TypeIII


## Percentage of explained variability (MALES) ##
var.age <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.age`)
var.area <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.area`)
var.inter <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.area.age`)
round(100*c(var.age,var.area,var.inter)/(var.age+var.area+var.inter),1)

## Percentage of explained variability (FEMALES) ##
var.age <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.age`)
var.area <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.area`)
var.inter <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.area.age`)
round(100*c(var.age,var.area,var.inter)/(var.age+var.area+var.inter),1)


########################################################################
## Figure 5: Sex-specific spatial patterns of suicide mortality rates ##
########################################################################
patterns.M <- lapply(FinalModel.M$marginals.lincomb.derived, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
patterns.M.post <- lapply(patterns.M, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

patterns.F <- lapply(FinalModel.F$marginals.lincomb.derived, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
patterns.F.post <- lapply(patterns.F, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

AGE.patterns <- bind_rows(do.call(rbind,patterns.M.post[grep("^Age\\.", names(patterns.M.post))]) |> 
                            as_tibble() |> 
                            rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                            mutate(Sex="Males", Age=unique(Suicides$Age)),
                          do.call(rbind,patterns.F.post[grep("^Age\\.", names(patterns.F.post))]) |> 
                            as_tibble() |> 
                            rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                            mutate(Sex="Females", Age=unique(Suicides$Age)))

SPATIAL.patterns <- bind_rows(do.call(rbind,patterns.M.post[grep("^Spatial\\.", names(patterns.M.post))]) |> 
                                as_tibble() |> 
                                rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                                mutate(Sex="Males", PROV=unique(Suicides$PROV)),
                              do.call(rbind,patterns.F.post[grep("^Spatial\\.", names(patterns.F.post))]) |> 
                                as_tibble() |> 
                                rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                                mutate(Sex="Females", PROV=unique(Suicides$PROV)))


## Plot of spatial patterns ##
aux.PROV <- SPATIAL.patterns |>
  select(PROV, Sex, `0.5quant`) |> 
  pivot_wider(names_from=Sex,
              values_from=`0.5quant`,
              names_prefix="Est.Rates.") |> 
  mutate(Prob.Males=unlist(lapply(FinalModel.M$marginals.random$ID.area, function(x) 1-inla.pmarginal(0,x))),
         Prob.Females=unlist(lapply(FinalModel.F$marginals.random$ID.area, function(x) 1-inla.pmarginal(0,x))))

carto <- merge(Carto_SpainPROV, aux.PROV)

colors1 <- c("#ffffd9", "#c7e9b4", "#7fcdbb", "#225ea8", "#0c2c84")
colors2 <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")

Fig5a <- tm_shape(carto) +
  tm_polygons(fill="Est.Rates.Males",
              fill.scale=tm_scale_intervals(values=colors1,
                                            breaks=c(-Inf,9,11,12,15,Inf)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(a)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig5b <- tm_shape(carto) +
  tm_polygons(fill="Prob.Males",
              fill.scale=tm_scale_intervals(values=colors2,
                                            breaks=c(0,0.1,0.2,0.8,0.9,1)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(b)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig5c <- tm_shape(carto) +
  tm_polygons(fill="Est.Rates.Females",
              fill.scale=tm_scale_intervals(values=colors1,
                                            breaks=c(-Inf,3,3.5,4.5,5,Inf)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(c)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig5d <- tm_shape(carto) +
  tm_polygons(fill="Prob.Females",
              fill.scale=tm_scale_intervals(values=colors2,
                                            breaks=c(0,0.1,0.2,0.8,0.9,1)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(d)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig5 <- tmap_arrange(Fig5a, Fig5b, Fig5c, Fig5d, nrow=2, ncol=2)
tmap_save(Fig5, filename="./Figures/Figure5.pdf", width=12, height=8)


###########################################################################
## Figure 7/8: Geographical distribution of mortality rates by age group ##
###########################################################################
rates.M <- lapply(FinalModel.M$marginals.linear.predictor, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
rates.M.post <- lapply(rates.M, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

rates.F <- lapply(FinalModel.F$marginals.linear.predictor, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
rates.F.post <- lapply(rates.F, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

EST.rates <- bind_rows(do.call(rbind,rates.M.post) |> 
                         as_tibble() |> 
                         rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                         mutate(Sex="Males",
                                PROV=FinalModel.M$.args$data$PROV,
                                Age=as.factor(FinalModel.M$.args$data$Age)),
                       do.call(rbind,rates.F.post) |> 
                         as_tibble() |> 
                         rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                         mutate(Sex="Females",
                                PROV=FinalModel.F$.args$data$PROV,
                                Age=as.factor(FinalModel.F$.args$data$Age)))

colors <- c("#ffffd9", "#c7e9b4", "#7fcdbb", "#225ea8", "#0c2c84")


## Maps of posterior median estimates for MALE population ##
Maps <- vector("list",8)
names(Maps) <- unique(EST.rates$Age)

titles <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)")
names(titles) <- unique(EST.rates$Age)
  
for(i in names(Maps)){
  carto <- Carto_SpainPROV
  carto$Est.Rates <- EST.rates |> filter(Sex=="Males", Age==i) |> pull(`0.5quant`)
  
  breaks <- quantile(carto$Est.Rates,(0:5)/5)
  labels <- sprintf("%.1f - %.1f", head(breaks,-1), tail(breaks,-1))
  
  Maps[[i]] <- tm_shape(carto) +
    tm_polygons(fill="Est.Rates",
                fill.scale=tm_scale_intervals(values=colors, breaks=breaks, labels=labels),
                fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                      position=tm_pos_out("right","center"))) +
    tm_title(text=titles[i], size=1.5) + 
    tm_options(component.autoscale = FALSE)
}

Fig7 <- tmap_arrange(Maps[[1]],Maps[[2]],Maps[[3]],Maps[[4]],
                     Maps[[5]],Maps[[6]],Maps[[7]],Maps[[8]],
                     nrow=3, ncol=3)

tmap_save(Fig7, filename="./Figures/Figure7.pdf", width=12, height=8)


## Maps of posterior median estimates for FEMALE population ##
Maps <- vector("list",8)
names(Maps) <- unique(EST.rates$Age)

titles <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)")
names(titles) <- unique(EST.rates$Age)

for(i in names(Maps)){
  carto <- Carto_SpainPROV
  carto$Est.Rates <- EST.rates |> filter(Sex=="Females", Age==i) |> pull(`0.5quant`)
  
  breaks <- quantile(carto$Est.Rates,(0:5)/5)
  labels <- sprintf("%.1f - %.1f", head(breaks,-1), tail(breaks,-1))
  
  Maps[[i]] <- tm_shape(carto) +
    tm_polygons(fill="Est.Rates",
                fill.scale=tm_scale_intervals(values=colors, breaks=breaks, labels=labels),
                fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                      position=tm_pos_out("right","center"))) +
    tm_title(text=titles[i], size=1.5) + 
    tm_options(component.autoscale = FALSE)
}

Fig8 <- tmap_arrange(Maps[[1]],Maps[[2]],Maps[[3]],Maps[[4]],
                     Maps[[5]],Maps[[6]],Maps[[7]],Maps[[8]],
                     nrow=3, ncol=3)

tmap_save(Fig8, filename="./Figures/Figure8.pdf", width=12, height=8)


###############################################################
## Figures S2/S3: Maps of posterior exceedance probabilities ##
###############################################################

## Maps of posterior exceedance probabilities for MALE population ##
ages <- FinalModel.M$summary.lincomb.derived$mean[1:8]
names(ages) <- unique(FinalModel.M$.args$data$Age)

Maps <- vector("list",8)
names(Maps) <- names(ages)

colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")

for(i in names(ages)){
  
  loc <- which(FinalModel.M$.args$data$Age==i)
  
  carto <- Carto_SpainPROV
  carto$Prob <- unlist(lapply(FinalModel.M$marginals.linear.predictor[loc], function(x) 1-inla.pmarginal(ages[i], x)))
  
  Maps[[i]] <- tm_shape(carto) +
    tm_polygons(fill="Prob",
                fill.scale=tm_scale_intervals(values=colors,
                                              breaks=c(0,0.1,0.2,0.8,0.9,1)),
                fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                      position=tm_pos_out("right","center"))) +
    tm_title(text=paste("Males: ",i)) + 
    tm_options(component.autoscale = FALSE)
}

FigS2 <- tmap_arrange(Maps[[1]],Maps[[2]],Maps[[3]],Maps[[4]],
                      Maps[[5]],Maps[[6]],Maps[[7]],Maps[[8]],
                      nrow=3, ncol=3)

tmap_save(FigS2, filename="./Figures/FigureS2.pdf", width=12, height=8)


## Maps of posterior exceedance probabilities for FEMALE population ##
ages <- FinalModel.F$summary.lincomb.derived$mean[1:8]
names(ages) <- unique(FinalModel.F$.args$data$Age)

Maps <- vector("list",8)
names(Maps) <- names(ages)

colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")

for(i in names(ages)){
  
  loc <- which(FinalModel.F$.args$data$Age==i)
  
  carto <- Carto_SpainPROV
  carto$Prob <- unlist(lapply(FinalModel.F$marginals.linear.predictor[loc], function(x) 1-inla.pmarginal(ages[i], x)))
  
  Maps[[i]] <- tm_shape(carto) +
    tm_polygons(fill="Prob",
                fill.scale=tm_scale_intervals(values=colors,
                                              breaks=c(0,0.1,0.2,0.8,0.9,1)),
                fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                      position=tm_pos_out("right","center"))) +
    tm_title(text=paste("Females: ",i)) + 
    tm_options(component.autoscale = FALSE)
}

FigS3 <- tmap_arrange(Maps[[1]],Maps[[2]],Maps[[3]],Maps[[4]],
                      Maps[[5]],Maps[[6]],Maps[[7]],Maps[[8]],
                      nrow=3, ncol=3)

tmap_save(FigS3, filename="./Figures/FigureS3.pdf", width=12, height=8)
