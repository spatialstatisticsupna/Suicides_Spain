rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
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


########################################
## Load previously fitted INLA models ##
########################################
load("INLAmodels_AgeTime.Rdata")

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
print(Table3, digits=4)


## Select final model ##
FinalModel.M <- MODELS.M$TypeII
FinalModel.F <- MODELS.F$TypeIV


## Percentage of explained variability (MALES) ##
var.age <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.age`)
var.time <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.year`)
var.inter <- inla.emarginal(function(x) 1/x, FinalModel.M$marginals.hyperpar$`Precision for ID.year.age`)
round(100*c(var.age,var.time,var.inter)/(var.age+var.time+var.inter),1)
round(100*c(var.time,var.inter)/(var.time+var.inter),1)

## Percentage of explained variability (FEMALES) ##
var.age <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.age`)
var.time <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.year`)
var.inter <- inla.emarginal(function(x) 1/x, FinalModel.F$marginals.hyperpar$`Precision for ID.year.age`)
round(100*c(var.age,var.time,var.inter)/(var.age+var.time+var.inter),1)
round(100*c(var.time,var.inter)/(var.time+var.inter),1)


#############################################################################
## Figure 4: Sex-specific age and time patterns of suicide mortality rates ##
#############################################################################
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

TIME.patterns <- bind_rows(do.call(rbind,patterns.M.post[grep("^Time\\.", names(patterns.M.post))]) |> 
                            as_tibble() |> 
                            rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                            mutate(Sex="Males", Year=unique(Suicides$Year)),
                          do.call(rbind,patterns.F.post[grep("^Time\\.", names(patterns.F.post))]) |> 
                            as_tibble() |> 
                            rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                            mutate(Sex="Females", Year=unique(Suicides$Year)))

## Plot of age group patterns ##
aux.AGE <- Suicides |> 
  group_by(Sex, Age) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5) |> 
  left_join(AGE.patterns, by=c("Sex","Age"))

Fig4a <- ggplot(aux.AGE, aes(x = Age)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`, group=Sex, fill=Sex), alpha=0.4) +
  geom_line(aes(y=`0.5quant`, color=Sex, group=Sex), linewidth=0.8) + 
  scale_y_continuous(
    name = "Estimated mortality rates (per 100,000 inhabitants)",
    expand = c(0,0),
    limits = c(0,41)
  ) +
  scale_x_discrete(expand = c(0.02,0.02)) + 
  xlab("Age groups") + 
  scale_color_manual(values = c("Males" = "#4f81bd", "Females" = "#d35400"), name = NULL) +
  scale_fill_manual(values = c("Males" = "#4f81bd", "Females" = "#d35400"), name = NULL) +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="top",
        legend.text=element_text(size=14),
        axis.title.y.right=element_text(angle=90))


## Plot of time patterns ##
aux.YEAR <- Suicides |> 
  group_by(Sex, Year) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5) |> 
  left_join(TIME.patterns, by=c("Sex","Year"))

Fig4b <- ggplot(aux.YEAR, aes(x = Year)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`, group=Sex, fill=Sex), alpha=0.4) +
  geom_line(aes(y=`0.5quant`, color=Sex, group=Sex), linewidth=0.8) + 
  scale_y_continuous(
    name = "Estimated mortality rates (per 100,000 inhabitants)",
    expand = c(0,0),
    limits = c(0,15)
  ) +
  scale_x_discrete(expand = c(0.02,0.02)) + 
  xlab("Years") + 
  scale_color_manual(values = c("Males" = "#4f81bd", "Females" = "#d35400"), name = NULL) +
  scale_fill_manual(values = c("Males" = "#4f81bd", "Females" = "#d35400"), name = NULL) +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="top",
        legend.text=element_text(size=14),
        axis.title.y.right=element_text(angle=90))

Fig4 <- ggarrange(Fig4a, Fig4b, nrow=1, ncol=2)
ggsave("./Figures/Figure4.pdf", Fig4, width=12, height=6)


##################################################################
## Figure 6: Temporal evolution of mortality rates by age group ##
##################################################################
rates.M <- lapply(FinalModel.M$marginals.linear.predictor, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
rates.M.post <- lapply(rates.M, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

rates.F <- lapply(FinalModel.F$marginals.linear.predictor, function(x) inla.tmarginal(function(y) 1e+5*exp(y), x))
rates.F.post <- lapply(rates.F, function(x) inla.qmarginal(c(0.025,0.5,0.975), x))

EST.rates <- bind_rows(do.call(rbind,rates.M.post) |> 
                         as_tibble() |> 
                         rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                         mutate(Sex="Males",
                                Year=FinalModel.M$.args$data$Year,
                                Age=as.factor(FinalModel.M$.args$data$Age)),
                       do.call(rbind,rates.F.post) |> 
                         as_tibble() |> 
                         rename('0.025quant'=V1, '0.5quant'=V2, '0.975quant'= V3) |> 
                         mutate(Sex="Females",
                                Year=FinalModel.F$.args$data$Year,
                                Age=as.factor(FinalModel.F$.args$data$Age)))

Fig6a <- ggplot(EST.rates |> filter(Sex=="Males"), aes(x=Year)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`, group=Age, fill=Age),
              alpha=0.4, show.legend=FALSE) +
  geom_line(aes(y=`0.5quant`, color=Age, group=Age), linewidth=0.8) + 
  xlab("Year") + 
  ylab("Estimated mortality rates (per 100,000 inhabitants)") +
  ggtitle("Male population") + 
  scale_colour_discrete(name=NULL) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_continuous(expand=c(0,0.1), limits=c(0,50)) + 
  guides(colour=guide_legend(reverse=TRUE)) + 
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="right",
        legend.text=element_text(size=10))

Fig6b <- ggplot(EST.rates |> filter(Sex=="Females"), aes(x=Year)) +
  geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`, group=Age, fill=Age),
              alpha=0.4, show.legend=FALSE) +
  geom_line(aes(y=`0.5quant`, color=Age, group=Age), linewidth=0.8) + 
  xlab("Year") + 
  ylab("Estimated mortality rates (per 100,000 inhabitants)") +
  ggtitle("Female population") + 
  scale_colour_discrete(name=NULL) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_continuous(expand=c(0,0.1), limits=c(0,8)) + 
  guides(colour=guide_legend(reverse=TRUE)) + 
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="right",
        legend.text=element_text(size=10))

Fig6 <- ggarrange(Fig6a, Fig6b, nrow=1, ncol=2)
ggsave("./Figures/Figure6.pdf", Fig6, width=12, height=6)
