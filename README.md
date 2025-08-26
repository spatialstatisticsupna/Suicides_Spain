# Suicides_Spain

This repository contains the data and R code to reproduce the analyses of the paper entitled *"Spatio-temporal analysis of suicide mortality rates in Spain: exploring sex, age, and socio-demographic disparities using Bayesian disease mapping models"* (Adin et al., 2025)

## Table of contents

-   [Data](#data)
-   [R code](#r-code)
-   [Acknowledgements](#Acknowledgements)


# Data

Sex- and age-specific suicide mortality counts (International Classification of Diseases-10 codes X60-X84), along with corresponding population data, stratified by year and province within continental Spain.

-   [**Suicides_Spain.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Data/Suicides_Spain.Rdata)

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

-   [**DGURBA_Spain.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Data/DGURBA_Spain.Rdata)

    This .Rdata file contains a `tibble` object with 8109 rows and 4 columns

    -   `CODMUNI`: character vector representing municipality identifiers
    -   `POB`: total population
    -   `DGURBA`: classification of Spanish municipalities by degree of urbanisation (1=cities, 2=towns or suburbs, 3=rural areas)
    -   `Tipo`: character vector indicating urban/rural areas

    **Data source: [GISCO (Geographic Information System of the European Commission)](https://ec.europa.eu/eurostat/web/gisco/geodata/population-distribution/degree-urbanisation)**

-   [**Unemployment.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Data/Unemployment.Rdata)

    The object `Unemployment.PROV.S` has 94 rows and 3 columns:
    
    -   `Sex`: factor with 2 levels ("Males" and "Females")
    -   `PROV`: character vector representing province identifiers
    -   `U.Rate`: average unemployment rate (%) between 2010 and 2022
    
    The object `Unemployment.PROV.ST` has 1222 rows and 4 columns:
    
    -   `Sex`: factor with 2 levels ("Males" and "Females")
    -   `PROV`: character vector representing province identifiers
    -   `Year`: character vector indicating the year (2010-2022)
    -   `U.Rate`: annual unemployment rate (%)
    
    **Data source: INE (Spanish Statistical Office)**
    
-   [**Poverty.Rdata**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Data/Poverty.Rdata)

    The object `Poverty.TA` has 104 rows and 4 columns:

    -   `Sex`: factor with 2 levels ("Males" and "Females")
    -   `Age`: character vector indicating age groups ("16-29", "30-44", "45-64", and "65+")
    -   `Year`: character vector indicating the year (2010-2022)
    -   `Rate`: annual at-risk-of-poverty rate (%)
    
    **Data source: INE (Spanish Statistical Office)**
    
    
# R code

R code to reproduce all analyses presented in this paper, including the fitting of age–time and age–space interaction models using INLA (<http://www.r-inla.org/>), as well as the code to generate all figures and tables.

-   [**1_DescriptiveAnalysis.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/1_DescriptiveAnalysis.R)

    Performs the descriptive analyses outlined in Section 2, including the generation of Figures 1–3 and Table 1.

-   [**2a_FitModels_AgeTime.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/2a_FitModels_AgeTime.R)

    Fits Bayesian hierarchical models that incorporate age–time interaction effects, as detailed in Section 3.1.

-   [**2b_FitModels_AgeSpace.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/2b_FitModels_AgeSpace.R)

    Fits Bayesian hierarchical models that incorporate age-space interaction effects, as detailed in Section 3.2.

-   [**3a_Results_AgeTime.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/3a_Results_AgeTime.R)

    Performs model comparison and analyzes the estimated rates for age–time interaction models (Section 4.1), including the generation of Table 3, Figure 4 and Figure 6.

-   [**3b_Results_AgeSpace.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/3b_Results_AgeSpace.R)

    Performs model comparison and analyzes the estimated rates for age-space interaction models (Section 4.1), including the generation of Table 3, Figure 5, Figures 7-8 and Figures S2-S3.

-   [**4a_EcologicalRegression_Spatial.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/4a_EcologicalRegression_Spatial.R)

    Fits ecological regression model under the simplified spatial+ appproach to address confounding for spatial covariates.

-   [**4b_EcologicalRegression_Temporal.R**](https://github.com/spatialstatisticsupna/Suicides_Spain/blob/master/Rcode/4b_EcologicalRegression_Temporal.R)

    Fits ecological regression model under the simplified spatial+ appproach to address confounding for temporal covariates.
    

# Acknowledgements

This work has been supported by projects PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033 (Spanish Ministry of Science and Innovation) and UNEDPAM/PI/PR24/01A (Centro Asociado UNED - Pamplona).

<p float="left">
  <img src="https://github.com/spatialstatisticsupna/Suicides_Spain/blob/main/micin-aei.jpg" width="49%" />
  <img src="https://github.com/spatialstatisticsupna/Suicides_Spain/blob/main/UNED_Pamplona.jpg" width="40%" />
</p>

