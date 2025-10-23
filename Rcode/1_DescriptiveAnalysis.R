rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(tmap)
library(tidyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


####################################
## Load data and cartography file ##
####################################
load("../Data/Suicides_Spain.Rdata")
print(Suicides)

## Set the observed and population data for the province of Madrid to NA for the years 2010–2012 ##
Suicides[Suicides$PROV=="28" & Suicides$Year %in% 2010:2012, c("O","Pop")] <- NA

## Delete from our analysis the age group "0-9" ##
Suicides <- Suicides |> filter(Age!="0-9")

n.year <- length(unique(Suicides$Year))
n.area <- length(unique(Suicides$PROV))
n.sex <- length(unique(Suicides$Sex))
n.age <- length(unique(Suicides$Age))
nrow(Suicides)==n.year*n.area*n.sex*n.age

## Average mortality rate over the the period 2010-2022 ##
Suicides |> 
  group_by(Sex) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)


###########################################################################
## Figure S1: Map with the administrative division of Spain by provinces ##
###########################################################################
FigS1 <- tm_shape(Carto_SpainPROV) + 
  tm_polygons(fill="lavender") + 
  tm_text(text="Name", size=0.8, fontface=2)

print(FigS1)
tmap_save(FigS1, filename="./Figures/FigureS1.pdf", width=12, height=8)


#################################################################################
## Figure 1: Average annual number of suicide deaths and crude mortality rates ##
##           by age group and sex in continental Spain, 2010-2022              ##
#################################################################################
aux <- Suicides |> 
  group_by(Sex, Age) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)

aux$RateLabel <- ifelse(aux$Sex=="Females", "Female mortality rates", "Male mortality rates")
aux$DeathLabel <- ifelse(aux$Sex=="Females", "Female deaths", "Male deaths")

scale_factor <- max(aux$Rate)/max(aux$O)

Fig1 <- ggplot(aux, aes(x=Age)) +
  geom_bar(aes(y=O*scale_factor, fill=DeathLabel),
           stat="identity", position="dodge", alpha=0.7, width=0.8) +
  geom_line(aes(y=Rate, color=RateLabel, group=RateLabel), linewidth=1.2) +
  scale_y_continuous(
    name = "Suicide mortality rate per 100,000 inhabitants",
    expand = c(0,0),
    sec.axis = sec_axis(~./scale_factor, name="Number of deaths by year"),
    limits = c(0,40)
  ) +
  scale_x_discrete(expand = c(0,0)) + 
  xlab("Age groups") + 
  scale_fill_manual(values=c("Female deaths"="#f4a582", "Male deaths"="#92c5de"), name=NULL) +
  scale_color_manual(values=c("Female mortality rates"="#d35400", "Male mortality rates"="#4f81bd"), name=NULL) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linetype = 1, shape = NA)),
    fill = guide_legend(order = 2)
  ) +
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="top",
        legend.text=element_text(size=14),
        axis.title.y.right=element_text(angle=90))

plot(Fig1)
ggsave("./Figures/Figure1.pdf", Fig1, width=12, height=8)


################################################################################
## Figure 2: Maps of suicide mortality rates by province in continental Spain ##
################################################################################
aux <- Suicides |> 
  group_by(Sex, PROV) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)

carto <- merge(Carto_SpainPROV, aux, by="PROV")

colors <- c("#ffffd9", "#c7e9b4", "#7fcdbb", "#1d91c0", "#225ea8", "#0c2c84")

Fig2a <- tm_shape(carto |> filter(Sex=="Males")) +
  tm_polygons(fill="Rate",
              fill.scale=tm_scale_intervals(values=colors,
                                            breaks=c(-Inf,10,13,14,15,18,Inf)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(a)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig2b <- tm_shape(carto |> filter(Sex=="Females")) +
  tm_polygons(fill="Rate",
              fill.scale=tm_scale_intervals(values=colors,
                                            breaks=c(-Inf,3,4,5,6,Inf)),
              fill.legend=tm_legend(title="", reverse=TRUE, frame=FALSE,
                                    position=tm_pos_out("right","center"))) +
  tm_title(text="(b)", size=2) + 
  tm_options(component.autoscale = FALSE)

Fig2 <- tmap_arrange(Fig2a, Fig2b, nrow=1, ncol=2)
tmap_save(Fig2, filename="./Figures/Figure2.pdf", width=12, height=5)


##########################################################################
## Figure 3: Annual suicide mortality rates by sex and age groups.      ##
##           The dashed black line indicates the overall mortality rate ##
##########################################################################
aux.Males <- Suicides |> 
  filter(Sex=="Males") |> 
  group_by(Year) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5, Age=NA)

aux.Females <- Suicides |> 
  filter(Sex=="Females") |> 
  group_by(Year) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5, Age=NA)

## Average annual increase ##
100*((aux.Males$Rate[n.year]/aux.Males$Rate[1])^(1/(n.year))-1)
100*((aux.Females$Rate[n.year]/aux.Females$Rate[1])^(1/(n.year))-1)

aux <- Suicides |> 
  group_by(Sex, Year, Age) |> 
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=O/Pop*1e+5)

Fig3a <- ggplot(aux |> filter(Sex=="Males"),
       aes(x=Year, y=Rate, colour=Age, group=Age)) + 
  geom_line(linewidth=0.8) +
  geom_line(data=aux.Males, aes(x=Year, y=Rate),
            color="black", linewidth=1.2, linetype = "dashed",
            show.legend=FALSE) + 
  xlab("Year") + 
  ylab("Suicide mortality rates per 100,000 inhabitants") +
  ggtitle("(a)") + 
  scale_colour_discrete(name=NULL) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_continuous(expand=c(0,0.1), limits=c(0,50)) + 
  guides(colour=guide_legend(reverse=TRUE)) + 
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="right",
        legend.text=element_text(size=10),
        plot.title=element_text(size=rel(2)))

Fig3b <- ggplot(aux |> filter(Sex=="Females"),
                aes(x=Year, y=Rate, colour=Age, group=Age)) + 
  geom_line(linewidth=0.8) +
  geom_line(data=aux.Females, aes(x=Year, y=Rate),
            color="black", linewidth=1.2, linetype = "dashed",
            show.legend=FALSE) + 
  xlab("Year") + 
  ylab("Suicide mortality rates per 100,000 inhabitants") +
  ggtitle("(b)") + 
  scale_colour_discrete(name=NULL) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_continuous(expand=c(0,0.1), limits=c(0,8)) + 
  guides(colour=guide_legend(reverse=TRUE)) + 
  theme_minimal(base_size=13) +
  theme(axis.text.x=element_text(angle=45, vjust=0.7),
        legend.position="right",
        legend.text=element_text(size=10),
        plot.title=element_text(size=rel(2)))

Fig3 <- ggarrange(Fig3a, Fig3b, nrow=1, ncol=2)
ggsave("./Figures/Figure3.pdf", Fig3, width=12, height=6)


#############################################################################################
## Table 1: Suicide mortality rates (per 100,000 inhabitants), with 95% credible intervals ##
##          in brackets, by sex, age groups, and rurality level in Spanish provinces       ##
#############################################################################################

## Classification of Spanish municipalities by degree of urbanisation (EUROSTAT)
## Urban=cities, and towns or suburbs / Rural=rural areas
load("../Data/DGURBA_Spain.Rdata")

## Calculate the percentage of the population living in rural areas at province level ##
RuralAreas.PROV <- DGURBA_Spain |> 
  mutate(PROV=substr(CODMUNI,1,2)) |> 
  group_by(PROV, Tipo) |> 
  summarise(POB=sum(POB), .groups="drop") |> 
  group_by(PROV) |> 
  mutate(Rural_percent=round(POB/sum(POB),3)) |> 
  ungroup() |> 
  filter(Tipo=="Rural") |> 
  select(-Tipo, -POB)

## Define the variable 'Level of rurality' ##
RuralAreas.PROV <- RuralAreas.PROV |> 
  mutate(Rural_level=case_when(Rural_percent > 0.4 ~ "High",
                               Rural_percent > 0.2 ~ "Intermediate",
                               TRUE               ~ "Low"),
         Rural_level=factor(Rural_level, levels=c("Low","Intermediate","High")))

## Compute the table ##
aux <- Suicides |> 
  left_join(RuralAreas.PROV, by="PROV") |> 
  mutate(Age2=case_when(Age=="10-19" ~ "40-",
                        Age=="20-29" ~ "40-",
                        Age=="30-39" ~ "40-",
                        Age=="40-49" ~ "40-70",
                        Age=="50-59" ~ "40-70",
                        Age=="60-69" ~ "40-70",
                        Age=="70-79" ~ "70+",
                        Age=="80+" ~ "70+")) |> 
  group_by(Rural_level, Age2, Sex) |>
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=1e+5*O/Pop,
         IC.inf=Rate-1.96*sqrt(O)/(Pop/1e+5),
         IC.sup=Rate+1.96*sqrt(O)/(Pop/1e+5))

Table1 <- lapply(aux |> group_split(Age2), function(x){
  
  x |> 
    select(Sex, Rural_level, Rate, IC.inf, IC.sup) |> 
    pivot_longer(cols=c(Rate,IC.inf,IC.sup),
                 names_to="Var",
                 values_to="Value") |> 
    mutate(Value=round(Value,1)) |> 
    pivot_wider(names_from=Sex,
                values_from=Value)
  
})
names(Table1) <- unique(aux$Age2)

Table1$Total <- aux |> 
  group_by(Rural_level, Sex) |>
  summarise(O=sum(O, na.rm=TRUE),
            Pop=sum(Pop, na.rm=TRUE),
            .groups = "drop") |> 
  mutate(Rate=1e+5*O/Pop,
         IC.inf=Rate-1.96*sqrt(O)/(Pop/1e+5),
         IC.sup=Rate+1.96*sqrt(O)/(Pop/1e+5)) |> 
  select(Sex, Rural_level, Rate, IC.inf, IC.sup) |> 
  pivot_longer(cols=c(Rate,IC.inf,IC.sup),
               names_to="Var",
               values_to="Value") |> 
  mutate(Value=round(Value,1)) |> 
  pivot_wider(names_from=Sex,
              values_from=Value)

print(Table1)
