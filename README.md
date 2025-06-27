# Suicides_Spain

This repository contains the data and R code to reproduce the analyses of the paper entitled *"Spatio-temporal analysis of suicide mortality rates in Spain: exploring sex, age, and socio-demographic disparities using Bayesian disease mapping models"* (Adin et al., 2025)

## Table of contents

-   [Data](#data)
-   [R code](#r-code)
-   [Supplementary Material](#supplementary-material)
-   [Acknowledgements](#Acknowledgements)
-   [References](#references)

# Data

Sex- and age-specific suicide mortality counts (International Classification of Diseases-10 codes X60-X84), along with corresponding population data, stratified by year and province within continental Spain.

-   [**Suicides_Spain.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/raw/master/Data/Suicides_Spain.Rdata)

    This .Rdata file contains the following objects:

    -   `Carto_SpainPROV`: `sf` object containing polygon geometries for the provinces of continental Spain

    -   `Suicides`: `tibble` object with 10.998 rows and 6 columns

        -   `Year`: character vector indicating the year (2010-2022)
        -   `PROV`: character vector representing province identifiers
        -   `Sex`: factor with 2 levels ("Males" and "Females")
        -   `Age`: character vector indicating age groups ("0-9", "10-19", ..., "70-79", "80+")
        -   `O`: observed number of suicide deaths
        -   `Pop`: population at risk

    **Data source: INE (Spanish Statistical Office)**

-   [**DGURBA_Spain.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/data/DGURBA_Spain.Rdata)

    This .Rdata file contains a `tibble` object with 8109 rows and 4 columns

    -   `CODMUNI`: character vector representing municipality identifiers
    -   `POB`: total population
    -   `DGURBA`: classification of Spanish municipalities by degree of urbanisation (1=cities, 2=towns or suburbs, 3=rural areas)
    -   `Tipo`: character vector indicating urban/rural areas

    **Data source: [GISCO (Geographic Information System of the European Commission)](https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/degree-urbanisation)**

# R code

R code to reproduce all analyses presented in this paper, including the fitting of age–time and age–space interaction models using INLA (<http://www.r-inla.org/>), as well as the code to generate all figures and tables.

-   [**1_DescriptiveAnalysis.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/1_DescriptiveAnalysis.R)

    Include description

-   [**2a_FitModels_AgeTime.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/2a_FitModels_AgeTime.R)

    Include description

-   [**2b_FitModels_AgeSpace.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/2b_FitModels_AgeSpace.R)

    Include description

-   [**3a_Results_AgeTime.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/3a_Results_AgeTime.R)

    Include description

-   [**3b_Results_AgeSpace.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/3b_Results_AgeSpace.R)

    Include description

-   [**4_EcologicalRegression.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/R/4_EcologicalRegression.R)

    PENDIENTE

# Supplementary Material

PENDIENTE

# Acknowledgements

This work has been supported by projects PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033 (Spanish Ministry of Science and Innovation) and UNEDPAM/PI/PR24/01A (Centro Asociado UNED - Pamplona).

<p float="left">
  <img src="https://github.com/spatialstatisticsupna/Suicides_Spain/blob/main/micin-aei.jpg" width="49%" />
  <img src="https://github.com/spatialstatisticsupna/Suicides_Spain/blob/main/UNED_Pamplona.jpg" width="40%" />
</p>


# References

PENDIENTE
