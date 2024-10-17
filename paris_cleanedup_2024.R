
# Things to run before any session:

# miscellaneous
library(gdata)
library(plyr)
library(countrycode)
library(doBy)
library(reshape2)

# Population prediction libraries
library(wpp2022)
library(bayesTFR)
library(bayesLife)
library(bayesPop)
library(popReconstruct)

# For model fitting with JAGS
library(rjags)
library(coda)

# for plotting, displaying
library(xtable)
library(ggplot2)
library(gridExtra)
library(scales)
library(MASS)
library(rworldmap)

library(car)

data.location <- "NatureData/"
sims.location <- paste0(data.location, "Simulations/")
plot.location <- "NatureData/Plots/"


year.start <- 2023
year.end <- 2100
n.trajectories <- 1000


load(file=paste0(data.location, paste0("poppreds_formatted_", year.start, "_", year.end, ".rda")))
load(file=paste0(data.location, "model_results_ar1const_2023.rda"))
load(file=paste0(data.location, "proj_evals_ar1const_2023.rda"))
load(file=paste0(data.location, "data_medium_2024.Rda"))
paris.objective <- read.csv(paste0(data.location, 'paris_agreement_new.csv')) # load 'paris_agreement_new_ndc.csv' for 2022 updated NDCs
paris.objective.category <- read.csv(paste0(data.location, 'paris_agreement_summary.csv'))
data.medium = data_medium_new
load(file=paste0(data.location, "cleaned_all_data_2024.Rda"))
cleaned_all_data <- cleaned_all_data_2024
load(file=paste0(data.location,'real_hadcrut5.Rda'))
# load("~/Documents/Research/Climate/Data/proj.evals.2023.adjusted.new.Rda")
real_data = hadcrut5
## Paris Agreement Objectives
paris.objective <- paris.objective[, c(1, 3:8)]
paris.objective <- paris.objective[!is.na(paris.objective$size), ]

EU.countries <- c('Belgium', 'Bulgaria', 'Croatia',
                  'Czech Republic', 'Denmark', 'Germany', 'Estonia', 'Ireland',
                  'Greece', 'Spain', 'France', 'Italy', 'Cyprus', 'Latvia', 'Lithuania',
                  'Luxembourg', 'Hungary', 'Malta', 'Netherlands', 'Austria',
                  'Poland', 'Portugal', 'Romania', 'Slovenia', 'Slovakia', 'Finland',
                  'Sweden', 'United Kingdom')

eu.objective <- data.frame(Country = EU.countries,
                           country_code = countrycode(EU.countries, 'country.name', 'iso3c'),
                           size = 40,
                           year = 2030,
                           old_year = 1990,
                           Type = 'Emission', ratio = NA)

paris.objective <- rbind(paris.objective, eu.objective)

paris.objective.category <- paris.objective.category[, c('Country', 'country_code', 'Type')]
eu.category <- data.frame(Country = EU.countries,
                          country_code = countrycode(EU.countries, 'country.name', 'iso3c'), Type = 'Emission')
paris.objective.category <- rbind(paris.objective.category, eu.category)
paris.objective.category.included <- subset(paris.objective.category, country_code %in% unique(data.medium$Isocode))

# China, USA, India, Japan, Germany, Russia, Canada, Korea, Brazil, UK,
# debug this fucntion
convert.traj <- function(proj.evals, paris.objective, policy.continue=FALSE)
{
  ipat.components.bycountry <- proj.evals$ipat.components.bycountry
  param.indeces <- round(seq(1, floor(1e5 / 20), length.out=(length(ipat.components.bycountry) / 5)))
  corr_level <- as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[1]][param.indeces, 'rho'])
  for (j in 2:5)
  {
    corr_level <- c(corr_level, as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[j]][param.indeces, 'rho']))
  }
  adjusted_countries <- c()
  for (l in 1:nrow(paris.objective))
  {
    # if(paris.objective$Country[l] == 'Gambia'){browser()}
    country_code <- as.character(paris.objective$country_code[l])
    type <- paris.objective$Type[l]
    if (!(country_code %in% as.character(unique(data.medium$Isocode))))
    {
      next
    }
    cat('country ', as.character(paris.objective$Country[l]), 'Started.\n')
    if (type == 'Emission' || (type == 'BAU' && !is.na(paris.objective$ratio[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- ifelse(type == 'Emission', paris.objective$size[l], size <- paris.objective$ratio[l])
      year_ind <- ceiling((paris.objective$year[l] - 2023) / 1) # change to 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1) } # change to 77 means starting at 2023 (2023+77 = 2100)
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023 (2022+78=2100)
      ref_level <- data.medium$CO2.total[(data.medium$Isocode == country_code) & (data.medium$Year == old_year)] * (1 - size / 100) * 11/3 * 1e6
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$CO2[year_ind] < ref_level)
        {
          next
        }
        ratio <- ref_level / ipat.components.bycountry[[j]][[country_code]]$CO2[year_ind]

        tech_ratio <- ratio^(1/(1 - corr_level[j]))
        tech_ratio <- tech_ratio^exponential_change
        if (country_code != 'USA')
        {
          gdp_ratio <- ratio^(corr_level[j]/(1 - corr_level[j]))
          gdp_ratio <- gdp_ratio^exponential_change
          ipat.components.bycountry[[j]][[country_code]]$FrontierGap <-
            ipat.components.bycountry[[j]][[country_code]]$FrontierGap + log(gdp_ratio)
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
            proj.evals$ipat.components.bycountry[[j]]$USA$GDPpercapita /
            exp(ipat.components.bycountry[[j]][[country_code]]$FrontierGap)
        }
        else
        {
          # browser()
          gdp_ratio <- ratio^(-corr_level[j]/(1 - corr_level[j]))
          gdp_ratio <- gdp_ratio^exponential_change
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
            ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * gdp_ratio
        }

        ipat.components.bycountry[[j]][[country_code]]$Tech <-
          ipat.components.bycountry[[j]][[country_code]]$Tech * tech_ratio
        ipat.components.bycountry[[j]][[country_code]]$GDP <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Pop * 1e3
        ipat.components.bycountry[[j]][[country_code]]$CO2 <-
          ipat.components.bycountry[[j]][[country_code]]$GDP * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
        ipat.components.bycountry[[j]][[country_code]]$CO2percapita <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
      }
    }
    else if (type == 'Intensity'|| (type == 'BAU' && !is.na(paris.objective$old_year[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- paris.objective$size[l]
      year_ind <- ceiling((paris.objective$year[l] - 2023) / 1) # change to 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1) } # change to 77 means starting at 2023 (2023+77 = 2100)
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023 (2022+78=2100)
      ref_level <- data.medium$Tech[(data.medium$Isocode == country_code) & (data.medium$Year == old_year)] * (1 - size / 100)
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$Tech[year_ind] < ref_level)
        {
          next
        }
        ratio <- ref_level / ipat.components.bycountry[[j]][[country_code]]$Tech[year_ind]

        tech_ratio <- ratio
        tech_ratio <- tech_ratio^exponential_change
        gdp_ratio <- ratio^(corr_level[j])
        gdp_ratio <- gdp_ratio^exponential_change
        ipat.components.bycountry[[j]][[country_code]]$FrontierGap <-
          ipat.components.bycountry[[j]][[country_code]]$FrontierGap + log(gdp_ratio)
        ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
          proj.evals$ipat.components.bycountry[[j]]$USA$GDPpercapita /
          exp(ipat.components.bycountry[[j]][[country_code]]$FrontierGap)

        ipat.components.bycountry[[j]][[country_code]]$Tech <-
          ipat.components.bycountry[[j]][[country_code]]$Tech * tech_ratio
        ipat.components.bycountry[[j]][[country_code]]$GDP <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Pop * 1e3
        ipat.components.bycountry[[j]][[country_code]]$CO2 <-
          ipat.components.bycountry[[j]][[country_code]]$GDP * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
        ipat.components.bycountry[[j]][[country_code]]$CO2percapita <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
      }
    }
    else
    {
      next
    }
    adjusted_countries <- c(adjusted_countries, as.character(paris.objective$country_code[l]))
    cat('country ', as.character(paris.objective$Country[l]), 'Finished.\n')
  }

  return (list(ipat.components.bycountry = ipat.components.bycountry, adjusted_countries = adjusted_countries))
}
# not useful, can be ingored
compute.ref.level <- function(proj.evals, paris.objective, policy.continue=FALSE)
{
  ipat.components.bycountry <- proj.evals$ipat.components.bycountry
  param.indeces <- round(seq(1, floor(1e5 / 20), length.out=(length(ipat.components.bycountry) / 5)))
  corr_level <- as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[1]][param.indeces, 'rho'])
  for (j in 2:5)
  {
    corr_level <- c(corr_level, as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[j]][param.indeces, 'rho']))
  }
  adjusted_countries <- c()
  ref_levels <- c()
  types <- c()
  for (l in 1:nrow(paris.objective))
  {
    # if(paris.objective$Country[l] == 'Gambia'){browser()}
    country_code <- as.character(paris.objective$country_code[l])
    type <- paris.objective$Type[l]
    if (!(country_code %in% as.character(unique(data.medium$Isocode))))
    {
      next
    }
    cat('country ', as.character(paris.objective$Country[l]), 'Started.\n')
    if (type == 'Emission' || (type == 'BAU' && !is.na(paris.objective$ratio[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- ifelse(type == 'Emission', paris.objective$size[l], size <- paris.objective$ratio[l])
      year_ind <- ceiling((paris.objective$year[l] - 2023) / 1) # change to 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1) } # change to 77 means starting at 2023 (2023+77=2100)
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023 (2022+78=2100)
      ref_level <- data.medium$CO2.total[(data.medium$Isocode == country_code) & (data.medium$Year == old_year)] * (1 - size / 100) * 11/3 * 1e6
      type <- 'Emission'
      ref_levels <- c(ref_levels, ref_level)
      types <- c(types, type)
    }
    else if (type == 'Intensity'|| (type == 'BAU' && is.na(paris.objective$old_year[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- paris.objective$size[l]
      year_ind <- ceiling((paris.objective$year[l] - 2023) / 1) # change to 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1) } # change to 77 means starting at 2023 (2023+77=2100)
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023 (2022+78=2100)
      ref_level <- data.medium$Tech[(data.medium$Isocode == country_code) & (data.medium$Year == 2023)] * (1 - size / 100)
      type <- 'Intensity'
      ref_levels <- c(ref_levels, ref_level)
      types <- c(types, type)
    }
    else
    {
      next
    }
    adjusted_countries <- c(adjusted_countries, as.character(paris.objective$country_code[l]))
    cat('country ', as.character(paris.objective$Country[l]), 'Finished.\n')
  }

  return (list(adjusted_countries = adjusted_countries, types = types, ref_levels = ref_levels))
}
# not useful, can be ingored
#res <- compute.ref.level(proj.evals.2023.ar1.const, paris.objective)
# output_tbl <- data.frame(Country = paris.objective$Country, country_code = paris.objective$country_code,
#                          year = paris.objective$year, reference_year = paris.objective$old_year,
#                          type = paris.objective$Type)
# res_tbl <- data.frame(country_code = res$adjusted_countries, types = res$types, ref_levels = res$ref_levels)
# res_tbl <- base::merge(output_tbl, res_tbl, by.x = "country_code", by.y = "country_code")

# debug this fucntion
convert.traj.annual <- function(proj.evals, paris.objective, policy.continue=FALSE)
{
  ipat.components.bycountry <- proj.evals$ipat.annual.components.bycountry
  param.indeces <- round(seq(1, floor(1e5 / 20), length.out=(length(ipat.components.bycountry) / 5)))
  corr_level <- as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[1]][param.indeces, 'rho'])
  for (j in 2:5)
  {
    corr_level <- c(corr_level, as.numeric(model.ar1.const.results.2023$corr.model.out[[2]][[j]][param.indeces, 'rho']))
  }
  adjusted_countries <- c()
  for (l in 1:nrow(paris.objective))
  {
    country_code <- as.character(paris.objective$country_code[l])
    type <- paris.objective$Type[l]
    if (!(country_code %in% as.character(unique(data.medium$Isocode))))
    {
      next
    }
    cat('country ', as.character(paris.objective$Country[l]), 'Started.\n')
    if (type == 'Emission' || (type == 'BAU' && !is.na(paris.objective$ratio[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- ifelse(type == 'Emission', paris.objective$size[l], size <- paris.objective$ratio[l])
      year_ind <- paris.objective$year[l] - 2023 # change to 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1)} # change to 77 means starting at 2023
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023(2022)
      ref_level <- data.medium$CO2.total[(data.medium$Isocode == country_code) &
                                           (data.medium$Year == old_year)] * (1 - size / 100) * 11/3 * 1e6
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$CO2[year_ind] < ref_level)
        {
          next
        }
        ratio <- ref_level / ipat.components.bycountry[[j]][[country_code]]$CO2[year_ind]

        tech_ratio <- ratio^(1/(1 - corr_level[j]))
        tech_ratio <- tech_ratio^exponential_change
        if (country_code != 'USA')
        {
          gdp_ratio <- ratio^(corr_level[j]/(1 - corr_level[j]))
          gdp_ratio <- gdp_ratio^exponential_change
          ipat.components.bycountry[[j]][[country_code]]$FrontierGap <-
            ipat.components.bycountry[[j]][[country_code]]$FrontierGap + log(gdp_ratio)
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
            ipat.components.bycountry[[j]]$USA$GDPpercapita /
            exp(ipat.components.bycountry[[j]][[country_code]]$FrontierGap)
        }
        else
        {
          # browser()
          gdp_ratio <- ratio^(-corr_level[j]/(1 - corr_level[j]))
          gdp_ratio <- gdp_ratio^exponential_change
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
            ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * gdp_ratio
        }

        ipat.components.bycountry[[j]][[country_code]]$Tech <-
          ipat.components.bycountry[[j]][[country_code]]$Tech * tech_ratio
        ipat.components.bycountry[[j]][[country_code]]$GDP <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Pop * 1e3
        ipat.components.bycountry[[j]][[country_code]]$CO2 <-
          ipat.components.bycountry[[j]][[country_code]]$GDP * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
        ipat.components.bycountry[[j]][[country_code]]$CO2percapita <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
      }
    }
    else if (type == 'Intensity'|| (type == 'BAU' && !is.na(paris.objective$old_year[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- paris.objective$size[l]
      year_ind <- paris.objective$year[l] - 2023 # change 2023
      if (policy.continue) {exponential_change <- (0:77) / (year_ind - 1)} # change to 77 means starting at 2023
      else {exponential_change <- c((0:(year_ind - 1))/(year_ind - 1), rep(1, 78-year_ind))} # change to 78 means starting at 2023
      ref_level <- data.medium$Tech[(data.medium$Isocode == country_code) & (data.medium$Year == old_year)] * (1 - size / 100)
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$Tech[year_ind] < ref_level)
        {
          next
        }
        ratio <- ref_level / ipat.components.bycountry[[j]][[country_code]]$Tech[year_ind]

        tech_ratio <- ratio
        tech_ratio <- tech_ratio^exponential_change
        gdp_ratio <- ratio^(corr_level[j])
        gdp_ratio <- gdp_ratio^exponential_change
        ipat.components.bycountry[[j]][[country_code]]$FrontierGap <-
          ipat.components.bycountry[[j]][[country_code]]$FrontierGap + log(gdp_ratio)
        ipat.components.bycountry[[j]][[country_code]]$GDPpercapita <-
          ipat.components.bycountry[[j]]$USA$GDPpercapita /
          exp(ipat.components.bycountry[[j]][[country_code]]$FrontierGap)

        ipat.components.bycountry[[j]][[country_code]]$Tech <-
          ipat.components.bycountry[[j]][[country_code]]$Tech * tech_ratio
        ipat.components.bycountry[[j]][[country_code]]$GDP <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Pop * 1e3
        ipat.components.bycountry[[j]][[country_code]]$CO2 <-
          ipat.components.bycountry[[j]][[country_code]]$GDP * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
        ipat.components.bycountry[[j]][[country_code]]$CO2percapita <-
          ipat.components.bycountry[[j]][[country_code]]$GDPpercapita * ipat.components.bycountry[[j]][[country_code]]$Tech / 1e4
      }
    }
    adjusted_countries <- c(adjusted_countries, as.character(paris.objective$country_code[l]))
    cat('country ', as.character(paris.objective$Country[l]), 'Finished.\n')
  }

  return (list(ipat.components.bycountry = ipat.components.bycountry, adjusted_countries = adjusted_countries))
}
# debug this fucntion
calculate.prob <- function(proj.evals, paris.objective)
{
  ipat.components.bycountry <- proj.evals$ipat.annual.components.bycountry
  param.indeces <- round(seq(1, floor(1e5 / 20), length.out=(length(ipat.components.bycountry) / 5)))
  countries <- names(proj.evals.2023.ar1.const$ipat.annual.components.bycountry[[1]])
  probs <- data.frame(country = countries, probability = NA)
  for (l in 1:nrow(paris.objective))
  {
    country_code <- as.character(paris.objective$country_code[l])
    type <- paris.objective$Type[l]
    if (!(country_code %in% as.character(unique(data.medium$Isocode))))
    {
      next
    }
    cat('country ', as.character(paris.objective$Country[l]), 'Started.\n')
    if (type == 'Emission' || (type == 'BAU' && !is.na(paris.objective$ratio[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- ifelse(type == 'Emission', paris.objective$size[l], size <- paris.objective$ratio[l])
      year_ind <- paris.objective$year[l] - 2022 # change 2014 here, forget to change last time
      ref_level <- data.medium$CO2.total[(data.medium$Isocode == country_code) &
                                           (data.medium$Year == old_year)] * (1 - size / 100) * 11/3 * 1e6
      count <- 0
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$CO2[year_ind] < ref_level) count <- count + 1
      }
      count <- count / j
      probs$probability[probs$country == country_code] <- count

    }
    else if (type == 'Intensity'|| (type == 'BAU' && !is.na(paris.objective$old_year[l])))
    {
      old_year <- paris.objective$old_year[l]
      size <- paris.objective$size[l]
      year_ind <- paris.objective$year[l] - 2022 # change 2014 here, forget to change last time
      ref_level <- data.medium$Tech[(data.medium$Isocode == country_code) & (data.medium$Year == old_year)] * (1 - size / 100)
      count <- 0
      for (j in 1:length(ipat.components.bycountry))
      {
        if (ipat.components.bycountry[[j]][[country_code]]$Tech[year_ind] < ref_level) count <- count + 1
      }
      count <- count / j
      probs$probability[probs$country == country_code] <- count
    }
    # adjusted_countries <- c(adjusted_countries, as.character(paris.objective$country_code[l]))
    cat('country ', as.character(paris.objective$Country[l]), 'Finished.\n')
  }

  return (probs)
}


converted_trajs <- convert.traj(proj.evals.2023.ar1.const, paris.objective)
converted_annual_trajs <- convert.traj.annual(proj.evals.2023.ar1.const, paris.objective)
converted_annual_trajs_cont <- convert.traj.annual(proj.evals.2023.ar1.const, paris.objective, TRUE)

proj.evals.2023.adjusted <- list(ipat.components.bycountry = converted_trajs$ipat.components.bycountry,
                                 ipat.annual.components.bycountry = converted_annual_trajs$ipat.components.bycountry,
                                 adjusted_countries = converted_annual_trajs$adjusted_countries)

proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont <-
  converted_annual_trajs_cont$ipat.components.bycountry
## Worldwide Trajectories
get.ipat.components.total <- function(ipat.components.bycountry.traj,
                                      year.sequence) {
  # get.ipat.components.total gets IPAT components for the world for our
  # projections.

  # US is special

  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  data.world <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.world) <- names.tmp
  data.world[, "Year"] <- year.sequence
  index <- 0
  pop.projections <- numeric(length(year.sequence))
  gdp.projections <- numeric(length(year.sequence))
  co2.projections <- numeric(length(year.sequence))
  for (country in names(ipat.components.bycountry.traj))
  {
    pop.projections <- pop.projections + ipat.components.bycountry.traj[[country]]$Pop * 1000
    gdp.projections <- gdp.projections + ipat.components.bycountry.traj[[country]]$GDP
    co2.projections <- co2.projections + ipat.components.bycountry.traj[[country]]$CO2
  }
  data.world$Pop <- pop.projections
  data.world$GDP <- gdp.projections
  data.world$CO2 <- co2.projections

  data.world$GDPpercapita <- data.world$GDP / data.world$Pop
  # Get carbon per capita. Note: this is what is CO2 in data.medium.
  data.world$CO2percapita <- data.world$CO2 / data.world$Pop
  # Get tech
  data.world$Tech <- data.world$CO2percapita / data.world$GDPpercapita * 10^4

  data.world
}


trajs.worldwide <- lapply(1:length(proj.evals.2023.adjusted$ipat.components.bycountry), function(i) {
  x <- proj.evals.2023.adjusted$ipat.components.bycountry[[i]]
  get.ipat.components.total(x, year.sequence=seq(year.start, year.end,1))
})

trajs.annual.worldwide <- lapply(1:length(proj.evals.2023.adjusted$ipat.annual.components.bycountry), function(i) {
  x <- proj.evals.2023.adjusted$ipat.annual.components.bycountry[[i]]
  get.ipat.components.total(x, year.sequence=year.start:year.end)
})

trajs.annual.worldwide.cont <- lapply(1:length(proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont), function(i) {
  x <- proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont[[i]]
  get.ipat.components.total(x, year.sequence=year.start:year.end)
})

trajs.annual.worldwide.usa <- lapply(1:length(proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont), function(i) {
  x <- proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont[[i]]
  x$USA <- proj.evals.2023.ar1.const$ipat.annual.components.bycountry[[i]]$USA
  get.ipat.components.total(x, year.sequence=year.start:year.end)
})

proj.evals.2023.adjusted$trajs.worldwide <- trajs.worldwide
proj.evals.2023.adjusted$trajs.annual.worldwide <- trajs.annual.worldwide
proj.evals.2023.adjusted$trajs.annual.worldwide.cont <- trajs.annual.worldwide.cont
proj.evals.2023.adjusted$trajs.annual.worldwide.usa <- trajs.annual.worldwide.usa

proj.evals.2023.adjusted$probability <- calculate.prob(proj.evals.2023.ar1.const, paris.objective)
## World Quantiles:
get.trajs.quants <- function(trajs.worldwide, quantiles=c(0.1,0.5,0.9)) {
  # get.trajs.quants gets quantiles for worldwide CO2 emissions
  year.seq <- trajs.worldwide[[1]]$Year
  co2_trajs <- matrix(0, nrow = length(trajs.worldwide), ncol = length(year.seq))
  for (j in 1:length(trajs.worldwide))
  {
    co2_trajs[j, ] <- trajs.worldwide[[j]]$CO2
  }
  trajs.quants <- apply(co2_trajs, 2, quantile, probs = quantiles) / 1e9
  trajs.quants <- cbind(data.frame(Quantile=quantiles, trajs.quants))
  names(trajs.quants)[-1] <- paste0("CO2", year.seq)

  trajs.quants
}


proj.evals.2023.adjusted$trajs.quants <- get.trajs.quants(proj.evals.2023.adjusted$trajs.worldwide,
                                                          quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
proj.evals.2023.adjusted$trajs.annual.quants <-
  get.trajs.quants(proj.evals.2023.adjusted$trajs.annual.worldwide,
                   quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
proj.evals.2023.adjusted$trajs.annual.quants.cont <-
  get.trajs.quants(proj.evals.2023.adjusted$trajs.annual.worldwide.cont,
                   quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))
proj.evals.2023.adjusted$trajs.annual.quants.usa <-
  get.trajs.quants(proj.evals.2023.adjusted$trajs.annual.worldwide.usa,
                   quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975))

# rcp.list <- list()
# for (rcp.num in c("RCP26", "RCP45", "RCP60", "RCP85")) {
#   rcp.list[[rcp.num]] <- read.csv(paste0(data.location, "IPCC_", rcp.num, "_data.csv"))
# }
rcp.colors <- c("green", "red", "black", "purple")
library(readxl)
rcp.carbon.yearly.tmp <- read_excel(paste0(data.location, "ssp_db_carbon_emissions.xlsx"),
                                    n_max = 4)
rcp.carbon.yearly <- rcp.carbon.yearly.tmp[, c("Scenario", "Unit",
                                               c(2000, 2005, seq(2010, 2100, by=10)))]
rcp.carbon.yearly[, -c(1,2)] <- (1/1000) * rcp.carbon.yearly[, -c(1,2)] # 1/1000 means Mt to Pg
rcp.carbon.yearly$Unit <- rep("PgCO2/yr", 4)
names(rcp.carbon.yearly)[3:length(names(rcp.carbon.yearly))] <- paste0("Carbon", names(rcp.carbon.yearly)[3:length(names(rcp.carbon.yearly))])

rcp.carbon.yearly$Scenario <- c("SSP3-7.0", "SSP2-4.5", "SSP1-2.6", "SSP5-8.5")

rcp.carbon.cum <- rcp.carbon.yearly[, c("Scenario", "Unit",
                                        paste0("Carbon", seq(2010, 2100, by=10)))]

rcp.carbon.yearly.complete <- rcp.carbon.yearly[, 1:2]
for (year in 2005:2100)
{
  var.name <- paste0('Carbon', year)
  if (year == 2005 || (year %% 10 == 0 && year <= 2100))
  {
    rcp.carbon.yearly.complete[, var.name] <- rcp.carbon.yearly[, var.name]
  }
  else if (year < 2010)
  {
    rcp.carbon.yearly.complete[, var.name] <- (2010 - year)/5 * rcp.carbon.yearly[, 'Carbon2005'] +
      (year - 2005)/5 * rcp.carbon.yearly[, 'Carbon2010']
  }
  else if (year < 2100)
  {
    rcp.carbon.yearly.complete[, var.name] <- (10 - year %% 10) / 10 * rcp.carbon.yearly[, paste0('Carbon', year - (year %% 10))] +
      (year %% 10) / 10 * rcp.carbon.yearly[, paste0('Carbon', year - (year %% 10) + 10)]
  }
  else
  {
    rcp.carbon.yearly.complete[, var.name] <- rcp.carbon.yearly[, 'Carbon2100']
  }
}

rcp.carbon.cum.complete <- rcp.carbon.yearly.complete
# for (i in 1:4)
# {
#   rcp.carbon.cum.complete[i, 4:98] <- cumsum(as.numeric(rcp.carbon.yearly.complete[i, 4:98]))
# }

for (i in 1:4) {
  cumulative_values <- cumsum(as.numeric(rcp.carbon.yearly.complete[i, 4:98]))
  rcp.carbon.cum.complete[i, 4:98] <- matrix(cumulative_values, nrow = 1)
}

rcp.carbon.cum.complete <- rcp.carbon.cum.complete[,-3]

rcp.carbon.cum.complete.adjusted <- t(rcp.carbon.cum.complete[, -c(1,2)])
rcp.carbon.cum.complete.adjusted <- as.data.frame(rcp.carbon.cum.complete.adjusted)
rcp.carbon.cum.complete.adjusted <- rcp.carbon.cum.complete.adjusted[, c(3,2,1,4)]
rcp.carbon.cum.complete.adjusted$year <- 2006:2100
names(rcp.carbon.cum.complete.adjusted)[1:4] <- paste0('ssp', c(126, 245, 370, 585))

rcp.carbon.cum <- rcp.carbon.yearly[, c("Scenario", "Unit",
                                        paste0("Carbon", seq(2010, 2100, by=10)))]
carbon.cum <- rcp.carbon.yearly[, paste0('Carbon', 2005)] * 2.5 + rcp.carbon.yearly[, paste0('Carbon', 2010)] * 2.5
for (year in seq(2010, 2100, by=10)) {
  var.name.tmp <- paste0("Carbon", year)
  if (year == 2010) {
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
  } else if (year == 2100) {
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
  } else {
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
    rcp.carbon.cum[, var.name.tmp] <- carbon.cum
    carbon.cum <- carbon.cum + 5*rcp.carbon.yearly[, var.name.tmp]
  }
}


all_models_split <- unlist(strsplit(names(cleaned_all_data$Average), '_s_'))[-1]
all_models <- unique(all_models_split[seq(2, 5282, 4)])
selected_indices <- c()
for (model in all_models)
{
  all_trajectories <- all_models_split[which(all_models_split == model) + 2]
  all_scenarios <- all_models_split[which(all_models_split == model) + 1]
  indices <- (which(all_models_split == model) + 2) / 4
  if ('ave.dat' %in% all_trajectories)
  {
    if (!is.na(which(all_trajectories == 'ave.dat')[1]))
    {
      selected_indices <- c(selected_indices, indices[which(all_trajectories == 'ave.dat')])
    }
    else
    {
      selected_indices <- c(selected_indices, indices[which(all_trajectories == 'ave.dat')])
    }
  }
  else
  {
    selected_indices <- c(selected_indices, indices)
  }
}


all_df <- data.frame()
for (index in (selected_indices + 1))
{
  temp_df <- cleaned_all_data$Average[, c(1, index)]
  names(temp_df) <- c('year', 'temperature')
  pre_industrial_avg <- mean(subset(temp_df, year >= 1981 & year <= 2005)$temperature)
  temp_df$temperature <- temp_df$temperature - pre_industrial_avg
  temp_df <- subset(temp_df, year >= 2006 & year <= 2100)
  val <- strsplit(names(cleaned_all_data$Average)[index], split = '_s_')[[1]]
  temp_df$model <- val[2]
  temp_df$scenario <- val[3]
  temp_df$co2 <- rcp.carbon.cum.complete.adjusted[, val[3]]
  all_df <- rbind(all_df, temp_df)
}

r2 <- c()
slope <- c()
for (model_x in unique(all_df$model))
{
  temp_df <- subset(all_df, model == model_x)
  m1 <- lm(temperature ~ co2, data = temp_df) # with intercept
  r2 <- c(r2, summary(m1)$r.squared)
  slope <- c(slope, summary(m1)$coef[2,1])
}

all_hist_df <- data.frame()
model <- 'initial'
for (index in (selected_indices + 1))
{
  temp_df <- cleaned_all_data$Average[, c(1, index)]
  val <- strsplit(names(cleaned_all_data$Average)[index], split = '_s_')[[1]]
  model_x <- val[2]
  if (model_x == model)
  {
    next
  }
  else
  {
    model <- model_x
  }
  names(temp_df) <- c('year', 'temperature')
  pre_industrial_avg <- mean(subset(temp_df, year >= 1981 & year <= 2005)$temperature)
  temp_df$temperature <- temp_df$temperature - pre_industrial_avg
  temp_df <- subset(temp_df, year >= 1861 & year <= 2005)
  temp_df$model <- val[2]
  temp_df$real_data <- real_data[12:156, 'V6'] - mean(real_data[132:156, 'V6'])
  all_hist_df <- rbind(all_hist_df, temp_df)
}

make_bugs_data <- function(all_hist_df, model_x)
{
  year <- 1861:2005
  temp_df <- subset(all_hist_df, model==model_x)
  model_forecast <- temp_df$temperature
  hadCrut_observation <- temp_df$real_data
  sd_hadCrut <- real_data$V3[real_data$V1 %in% year]
  model_forecast_anomaly <- model_forecast
  observed_anomaly <- hadCrut_observation
  n_years <- length(year)

  return (list(n_years = n_years,
               sd_delta = sd_hadCrut,
               model_forecast = model_forecast_anomaly,
               observed_anomaly = observed_anomaly))
}

fit_ar1_model <- function(all_hist_df, model_x,
                          n.iterations = 11000, n.adapt = 1000, n.chains = 3, thin = 1,
                          model_name = 'bayes_ar1.bug',
                          var.list = c("true_anomaly", "rho", "sd_w"))
{
  library(rjags)
  library(coda)
  library(reshape2)
  print("Calling make.bugs.data")
  jags.input <- make_bugs_data(all_hist_df, model_x)
  # browser()
  print("Successfully called make.bugs.data")
  jags.out <- jags.model(model_name,
                         data=jags.input,
                         n.chains=n.chains, n.adapt=n.adapt)
  jags.output <- coda.samples(jags.out,
                              var.list,
                              n.iter=n.iterations,
                              thin=thin)
  return (list(input = jags.input, output = jags.output))

}

####
# rcp.carbon.yearly: Yearly Emission from RCP estimate for every five years
# rcp.carbon.cum: Cumulative Emission from RCP estimate for every five years
# rcp.carbon.yearly.complete: Yearly Emission from RCP estimate for every year (2006-2100)
# rcp.carbon.cum.complete: Cumulative Emission from RCP estimate for every year (2006-2100)

# project_xt <- function(carbon_trajectory, rcp_co2, rcp_temp_data, year.range = c(2005, 2110))
# {
#   co2.trajectory <- carbon_trajectory$CO2 / 1e9
#   carbon.cum <- 0
#   proj_length <- length(min(2005, year.range[1]):year.range[2]) - 1
#   co2.trajectory.combined <- co2.trajectory[-1]
#   for (year in 2023:2006)
#   {
#     co2.trajectory.combined <- c(sum(data.medium$CO2.total[data.medium$Year == year] * 11/3)/1000, co2.trajectory.combined)
#   }
#   year_together <- 2006:year.range[2]
#   cum_co2 <- cumsum(co2.trajectory.combined)
#   cum_co2 <- cum_co2
#   # cum_co2 <- co2.trajectory[1] * 5 + cum_co2
#   rcp_co2_adjusted <- rcp_co2
#   rcp_co2_adjusted <- rcp_co2_adjusted[c(3,2,1,4), ]
#   rcp_co2_adjusted <- t(as.matrix(rcp_co2_adjusted[, -(1:2)]))
#   rcp_co2_adjusted <- as.data.frame(rcp_co2_adjusted)
#   names(rcp_co2_adjusted) <- paste0('rcp', c(26, 45, 60, 85))
#   rcp_co2_adjusted$year <- 2006:year.range[2]
#   rcp_temp_adjusted <- rcp_temp_data
#   names(rcp_temp_adjusted) <- c('year', paste0('rcp', c(26, 45, 60, 85)))
#   # projected_temp <- mean(as.numeric(rcp_temp_adjusted[1, -1]))
#   projected_temp <- c()
#   for (i in 1:proj_length)
#   {
#     year <- min(year.range[1], 2005) + i
#     co2_combined <- c(cum_co2[i], as.numeric(rcp_co2_adjusted[rcp_co2_adjusted$year == year, -5]))
#     index_traj <- which(sort(co2_combined, index.return = TRUE)$ix == 1)
#     if (index_traj == 1)
#     {
#       # browser()
#       all_rcp <- as.numeric(rcp_co2_adjusted[rcp_co2_adjusted$year == year, -5])
#       all_rcp_temp <- as.numeric(rcp_temp_adjusted[rcp_temp_adjusted$year == year, c('rcp26', 'rcp45', 'rcp60', 'rcp85')])
#       all_rcp <- sort(all_rcp)
#       all_rcp_temp <- sort(all_rcp_temp)
#       weight1 <- all_rcp[1] - cum_co2[i]
#       weight2 <- all_rcp[2] - cum_co2[i]
#       projected_temp_i <- (weight2 * all_rcp_temp[1] - weight1 * all_rcp_temp[2]) / (weight2 - weight1)
#       projected_temp <- c(projected_temp, projected_temp_i)
#       # projected_temp <- c(projected_temp, rcp_temp_adjusted[rcp_temp_adjusted$year == year, 'rcp26'])
#     }
#     else if (index_traj == 5)
#     {
#       # browser()
#       projected_temp <- c(projected_temp, rcp_temp_adjusted[rcp_temp_adjusted$year == year, 'rcp85'])
#     }
#     else
#     {
#       index_near <- sort(co2_combined, index.return = TRUE)$ix[which(sort(co2_combined, index.return = TRUE)$ix == 1) + c(-1, 1)]
#       #if(i == 19) {browser()}
#       distance_to_left_right <- abs(co2_combined[1] - co2_combined[index_near])
#       weights <- distance_to_left_right[2:1] / sum(distance_to_left_right)
#       projected_temp <- c(projected_temp, weighted.mean(rcp_temp_adjusted[rcp_temp_adjusted$year == year, index_near], weights))
#     }
#   }
#   return (projected_temp)
# }
#
#
# project_zt <- function(fit_result_chain, index)
# {
#   zt2005 <- fit_result_chain[index, 'true_anomaly[145]']
#   rho <- fit_result_chain[index, 'rho']
#   sd_w <- fit_result_chain[index, 'sd_w']
#   zt2005 <- zt2005 - fit_result$input$model_forecast[145]
#   zt_proj <- numeric(96)
#   zt_proj[1] <- zt2005
#   for (i in 2:96)
#   {
#     zt_proj[i] <- zt_proj[i - 1] * rho + rnorm(1, 0, sd_w)
#   }
#   return (zt_proj[-1])
# }
#
# project_temperature <- function(fit_result_chain, carbon.projection.trajectories, rcp_co2, cmip_temp, n.traj = 1000,
#                                 year.range = c(2005, 2110), ratio_co2 = 1)
# {
#   num_carbon_trajectories <- length(carbon.projection.trajectories)
#   num_zt_trajectories <- nrow(fit_result_chain)
#   carbon_trajs <- sample(num_carbon_trajectories, n.traj, replace = TRUE)
#   zt_trajs <- sample(num_zt_trajectories, n.traj, replace = TRUE)
#   length_proj <- (year.range[2] - min(year.range[1], 2005))
#   proj_temps <- matrix(nrow = length_proj, ncol = n.traj + 1)
#   proj_xts <- matrix(nrow = length_proj, ncol = n.traj + 1)
#   proj_zts <- matrix(nrow = length_proj, ncol = n.traj + 1)
#   proj_temps[, 1] <- proj_zts[, 1] <- proj_xts[, 1] <- (min(2005, year.range[1]) + 1):year.range[2]
#     # seq(year.range[1], year.range[2], 5)
#   mean1861_1880 <- numeric(n.traj)
#   n_years <- nrow(carbon.projection.trajectories[[carbon_trajs[i]]])
#   for (i in 1:n.traj)
#   {
#     test_func <- function(ratio, x, obj_ratio)
#     {
#       return ((sum(x * ratio^((1:length(x) - 1))) - sum(x * obj_ratio))^2)
#     }
#
#     ratio_annual <- optimize(test_func, c(0,1), tol = 0.0001, x = carbon.projection.trajectories[[carbon_trajs[i]]]$CO2 / 1e9,
#                           obj_ratio = ratio_co2)
#     #browser()
#     ratio_annual <- ratio_annual$minimum
#     # proj_xt <- project_xt(carbon.projection.trajectories[[carbon_trajs[i]]] * ratio_co2,
#     #                       rcp_co2, cmip_temp, year.range)
#     proj_xt <- project_xt(carbon.projection.trajectories[[carbon_trajs[i]]] * ratio_annual^((1:n_years - 1)),
#                           rcp_co2, cmip_temp, year.range)
#
#     proj_zt <- project_zt(fit_result_chain, zt_trajs[i])
#     # browser()
#     # proj_zt_cutting <- proj_zt[5 * (1:length_proj)]
#     # proj_zt_cutting <- proj_zt[seq(year.range[1], year.range[2], 5) - 2005]
#     proj_temp <- proj_xt + proj_zt#_cutting
#     proj_temps[, i + 1] <- proj_temp
#     proj_xts[, i + 1] <- proj_xt
#     proj_zts[, i + 1] <- proj_zt# _cutting
#   }
#
#   return(list(projection = proj_temps, xt_projection = proj_xts, zt_projection = proj_zts))
# }
### project temp without intercept
# project_xt <- function(carbon_trajectory, rcp_co2, model_x, year.range = c(2005, 2110), start.from.2023=TRUE)
# {
#   co2.trajectory <- carbon_trajectory$CO2 / 1e9
#   carbon.cum <- 0
#   proj_length <- length(year.range[1]:year.range[2]) - 1
#   co2.trajectory.combined <- co2.trajectory[-1]
#   if (start.from.2023)
#   {
#     for (year in 2023:2006)
#     {
#       co2.trajectory.combined <- c(sum(data.medium$CO2.total[data.medium$Year == year] * 11/3)/1000, co2.trajectory.combined)
#     }
#   }
#   cum_co2 <- cumsum(co2.trajectory.combined)
#   # cum_co2 <- co2.trajectory[1] * 5 + cum_co2
#   rcp_co2_adjusted <- rcp_co2
#   rcp_co2_adjusted <- rcp_co2_adjusted[c(3,2,1,4), ]
#   rcp_co2_adjusted <- t(as.matrix(rcp_co2_adjusted[, -(1:2)]))
#   rcp_co2_adjusted <- rcp_co2_adjusted[1:(length(year.range[1]:year.range[2]) - 1), ]
#   rcp_co2_adjusted <- as.data.frame(rcp_co2_adjusted)
#   names(rcp_co2_adjusted) <- paste0('rcp', c(26, 45, 60, 85))
#   rcp_co2_adjusted$year <- (year.range[1]+1):year.range[2]
#   rcp_temp_adjusted <- subset(all_df, model == model_x)
#   temp1 <- rcp_temp_adjusted$temperature[1]
#   rcp_temp_adjusted$temperature <- rcp_temp_adjusted$temperature - temp1
#   #temp_model <- lm(temperature ~ co2 + 0, data = rcp_temp_adjusted)
#
#   # projected_temp <- mean(as.numeric(rcp_temp_adjusted[1, -1]))
#   projected_temp <- as.numeric(predict(temp_model, newdata = data.frame(co2=cum_co2))) + temp1
#   return (projected_temp)
# }


project_xt <- function(carbon_trajectory, year.range = c(2005, 2110), start.from.2023=TRUE)
{
  co2.trajectory <- carbon_trajectory$CO2 / 1e9
  carbon.cum <- 0
  proj_length <- length(year.range[1]:year.range[2]) - 1
  co2.trajectory.combined <- co2.trajectory
  cum_co2 <- cumsum(co2.trajectory.combined)

  y1 <- real_data_adjusted$V6[nrow(real_data_adjusted)] # temp at 2023
  x1 <- co2.trajectory[1] # cum co2 at 2023
  slope <- rnorm(1, mean = 0.45, sd = 0.18)/1000 # AR6 distr
  intercept <- y1 - slope * x1

  projected_temp <- as.numeric(intercept + slope * cum_co2)
  return (projected_temp)
}



real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[132:156])

project_zt <- function(fit_result, index, year.end = 2100, year.start = 2005, proj_xt)
{
  fit_result_chain <- fit_result$output[[1]]
  zt2005 <- fit_result_chain[index, 'true_anomaly[145]']
  rho <- fit_result_chain[index, 'rho']
  sd_w <- fit_result_chain[index, 'sd_w']
  # D <- fit_result_chain[index, 'D']
  D <- 0
  zt_proj <- numeric(year.end - year.start + 1)
  if (year.start == 2005)
  {
    zt2005 <- zt2005 - fit_result$input$model_forecast[145]
    zt_proj[1] <- zt2005
  }
  else
  {
    idx <- which(real_data_adjusted$V1 == year.start)
    zt_start <- rnorm(1, real_data_adjusted$V6[idx], real_data_adjusted$V3[idx])
    zt_proj[1] <- zt_start - proj_xt[year.start - 2005]
  }
  for (i in 2:(year.end - year.start + 1))
  {
    zt_proj[i] <- (zt_proj[i - 1] - D) * rho + D + rnorm(1, 0, sd_w)
  }
  return (zt_proj[-1])
}

project_temperature <- function(fit_result_chain, carbon.projection.trajectories, rcp_co2, model_x, n.traj = 1000,
                                year.range = c(2005, 2110), ratio_co2 = 1, start.from.2023=TRUE)
{
  num_carbon_trajectories <- length(carbon.projection.trajectories)
  num_zt_trajectories <- nrow(fit_result_chain$output[[1]])
  carbon_trajs <- sample(num_carbon_trajectories, n.traj, replace = TRUE)
  zt_trajs <- sample(num_zt_trajectories, n.traj, replace = TRUE)
  length_proj <- year.range[2] - year.range[1]
  proj_temps <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_xts <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_zts <- matrix(nrow = length_proj, ncol = n.traj + 1)
  proj_temps[, 1] <- proj_zts[, 1] <- proj_xts[, 1] <- (year.range[1]+1):year.range[2]
  mean1861_1880 <- numeric(n.traj)
  for (i in 1:n.traj)
  {
    proj_xt <- project_xt(carbon.projection.trajectories[[carbon_trajs[i]]] * ratio_co2,
                          year.range, start.from.2023=start.from.2023)
    proj_zt <- project_zt(fit_result_chain, zt_trajs[i], year.end = year.range[2], year.start = year.range[1], proj_xt=proj_xt)
    proj_xt <- proj_xt[((year.range[1]+1):year.range[2])-2023] # change 2005 to 2023
    # proj_zt_cutting <- proj_zt[5 * (1:length_proj)]
    # proj_zt_cutting <- proj_zt[seq(year.range[1], year.range[2], 5) - 2005]
    proj_temp <- proj_xt + proj_zt # proj_zt_cutting
    proj_temps[, i + 1] <- proj_temp
    proj_xts[, i + 1] <- proj_xt
    proj_zts[, i + 1] <- proj_zt
  }

  return(list(projection = proj_temps, xt_projection = proj_xts, zt_projection = proj_zts))
}


rcp_temp_data <- cleaned_all_data$Average[, c(1, which(names(cleaned_all_data$Average) %in%
                                                         c('_s_one_s_ssp126_s_192',
                                                           '_s_one_s_ssp245_s_192',
                                                           '_s_one_s_ssp370_s_192',
                                                           '_s_one_s_ssp585_s_192')))]
rcp_temp_data <- subset(rcp_temp_data, year %in% 2006:2100)
names(rcp_temp_data) <- c('year', 'rcp26', 'rcp45', 'rcp60', 'rcp85') # in fact: 126 245 370 585
##########################################################################
all_fits <- list()
all_projections <- list()
proj_list <- list()
proj_list$projection <- c()
proj_list$xt_projection <- c()
proj_list$zt_projection <- c()

for (i in 1:54)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000) # original = 10000

  proj_list_model <- project_temperature(fit_result, proj.evals.2023.ar1.const$trajs.annual.worldwide, rcp.carbon.cum.complete,
                                   model_x, year.range = c(2023, 2100), n.traj = 1000)
  all_fits[[model_x]] <- fit_result
  all_projections[[model_x]] <- proj_list_model
  proj_list$projection <- cbind(proj_list$projection, all_projections[[model_x]]$projection)
  proj_list$xt_projection <- cbind(proj_list$projection, all_projections[[model_x]]$xt_projection)
  proj_list$zt_projection <- cbind(proj_list$projection, all_projections[[model_x]]$zt_projection)
}

all_fits <- list()
all_projections <- list()
proj_list_adjusted <- list()
proj_list_adjusted$projection <- c()
proj_list_adjusted$xt_projection <- c()
proj_list_adjusted$zt_projection <- c()

for (i in 1:54)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000)

  proj_list_model <- project_temperature(fit_result, proj.evals.2023.adjusted$trajs.annual.worldwide, rcp.carbon.cum.complete,
                                         model_x, year.range = c(2023, 2100), n.traj = 1000)
  all_fits[[model_x]] <- fit_result
  all_projections[[model_x]] <- proj_list_model
  proj_list_adjusted$projection <- cbind(proj_list_adjusted$projection, all_projections[[model_x]]$projection)
  proj_list_adjusted$xt_projection <- cbind(proj_list_adjusted$projection, all_projections[[model_x]]$xt_projection)
  proj_list_adjusted$zt_projection <- cbind(proj_list_adjusted$projection, all_projections[[model_x]]$zt_projection)
}

all_fits <- list()
all_projections <- list()
proj_list_adjusted_cont <- list()
proj_list_adjusted_cont$projection <- c()
proj_list_adjusted_cont$xt_projection <- c()
proj_list_adjusted_cont$zt_projection <- c()

for (i in 1:54)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000)

  proj_list_model <- project_temperature(fit_result, proj.evals.2023.adjusted$trajs.annual.worldwide.cont,
                                         rcp.carbon.cum.complete, model_x, year.range = c(2023, 2100), n.traj = 1000)
  all_fits[[model_x]] <- fit_result
  all_projections[[model_x]] <- proj_list_model
  proj_list_adjusted_cont$projection <- cbind(proj_list_adjusted_cont$projection, all_projections[[model_x]]$projection)
  proj_list_adjusted_cont$xt_projection <- cbind(proj_list_adjusted_cont$projection, all_projections[[model_x]]$xt_projection)
  proj_list_adjusted_cont$zt_projection <- cbind(proj_list_adjusted_cont$projection, all_projections[[model_x]]$zt_projection)
}

all_fits <- list()
all_projections <- list()
proj_list_adjusted_usa <- list()
proj_list_adjusted_usa$projection <- c()
proj_list_adjusted_usa$xt_projection <- c()
proj_list_adjusted_usa$zt_projection <- c()

for (i in 1:54)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000)

  proj_list_model <- project_temperature(fit_result, proj.evals.2023.adjusted$trajs.annual.worldwide.usa,
                                         rcp.carbon.cum.complete, model_x, year.range = c(2023, 2100), n.traj = 1000)
  all_fits[[model_x]] <- fit_result
  all_projections[[model_x]] <- proj_list_model
  proj_list_adjusted_usa$projection <- cbind(proj_list_adjusted_usa$projection, all_projections[[model_x]]$projection)
  proj_list_adjusted_usa$xt_projection <- cbind(proj_list_adjusted_usa$projection, all_projections[[model_x]]$xt_projection)
  proj_list_adjusted_usa$zt_projection <- cbind(proj_list_adjusted_usa$projection, all_projections[[model_x]]$zt_projection)
}


### Try plot map
map <- getMap(resolution = "coarse")
sPDF <- map[-which(map$ADMIN == "Antarctica"), ]
sPDF$UN <- sPDF$ISO_N3
sPDF$UN[sPDF$ISO3 == "CYN"] <- 196
sPDF$UN[sPDF$ISO3 == "KOS"] <- 688
sPDF$UN[sPDF$ISO3 == "ESH"] <- 732
prob <- rep(NA, length(sPDF$UN))
probability <- data.frame(iso = proj.evals.2023.adjusted$probability$country, un = NA,
                          prob = proj.evals.2023.adjusted$probability$probability)
for (i in 1:nrow(probability))
{
  probability$un[i] <- sPDF$UN[which(sPDF$ISO3 == as.character(probability$iso[i]))]
}

valididx <- which(is.element(sPDF$UN, probability$un))
prob[valididx] <- probability$prob[sapply(sPDF$UN[valididx], function(x, y) which(y==x), probability$un)]
sPDF$prob <- prob
main <- 'Probability'
mapParams <- mapCountryData(sPDF, nameColumnToPlot = 'prob', colourPalette = c('red', 'green'), numCats = 10,
                            ylim = c(-90,90), aspect = 'variable',catMethod = seq(0,1, 0.1), lwd = 1, borderCol = 'black',
                            addLegend = F, mapTitle = 'Probability', missingCountryCol = 'grey')
do.call(addMapLegend, c(mapParams, legendLabels = 'all', legendWidth=0.5, legendMar = 2))

mapParams <- mapCountryData(sPDF, nameColumnToPlot = 'prob', colourPalette = c('red', 'green'), numCats = 10,
                            mapRegion = "europe", aspect = 'variable',catMethod = seq(0,1, 0.1), lwd = 1, borderCol = 'black',
                            addLegend = F, mapTitle = 'Probability', missingCountryCol = 'grey')
do.call(addMapLegend, c(mapParams, legendLabels = 'all', legendWidth=0.5, legendMar = 2))


### Plot GMT Forecast
get_summary_temp_proj <- function(proj_list)
{
  proj_temp <- proj_list$projection
  summary_proj <- data.frame(year = proj_temp[,1])
  summary_proj$q025 <- apply(proj_temp[,-1], 1, quantile, probs = 0.025) - mean(real_data_adjusted$V6[12:31])
  summary_proj$q050 <- apply(proj_temp[,-1], 1, quantile, probs = 0.05) - mean(real_data_adjusted$V6[12:31])
  summary_proj$q100 <- apply(proj_temp[,-1], 1, quantile, probs = 0.1) - mean(real_data_adjusted$V6[12:31])
  summary_proj$q900 <- apply(proj_temp[,-1], 1, quantile, probs = 0.9) - mean(real_data_adjusted$V6[12:31])
  summary_proj$q950 <- apply(proj_temp[,-1], 1, quantile, probs = 0.95) - mean(real_data_adjusted$V6[12:31])
  summary_proj$q975 <- apply(proj_temp[,-1], 1, quantile, probs = 0.975) - mean(real_data_adjusted$V6[12:31])
  summary_proj$median <- apply(proj_temp[,-1], 1, median) - mean(real_data_adjusted$V6[12:31])
  summary_proj
}

summary_proj <- get_summary_temp_proj(proj_list)
summary_proj_adjusted <- get_summary_temp_proj(proj_list_adjusted)
summary_proj_adjusted_cont <- get_summary_temp_proj(proj_list_adjusted_cont)
summary_proj_adjusted_usa <- get_summary_temp_proj(proj_list_adjusted_usa)

true_anomaly_names <- paste('true_anomaly[', 1:145, ']', sep = '')
posterior_anomaly <- as.matrix(fit_result[[2]][[1]][, true_anomaly_names])
summary_anomaly <- data.frame(year = 1861:2005)
summary_anomaly$median <- apply(posterior_anomaly, 2, median)
summary_anomaly$q025 <- apply(posterior_anomaly, 2, quantile, probs = 0.025)
summary_anomaly$q050 <- apply(posterior_anomaly, 2, quantile, probs = 0.05)
summary_anomaly$q100 <- apply(posterior_anomaly, 2, quantile, probs = 0.1)
summary_anomaly$q900 <- apply(posterior_anomaly, 2, quantile, probs = 0.9)
summary_anomaly$q950 <- apply(posterior_anomaly, 2, quantile, probs = 0.95)
summary_anomaly$q975 <- apply(posterior_anomaly, 2, quantile, probs = 0.975)
summary_anomaly$observed_anomaly <- fit_result$input$observed_anomaly #?
summary_anomaly[,2:9] <- summary_anomaly[, 2:9] - mean(real_data_adjusted$V6[12:31])



plot_temp <- function(summary_proj, summary_anomaly, title = 'Anomaly Forecast')
{
  # summary_proj <- rbind(summary_anomaly[145, 1:8], summary_proj)
  other_estimates <- data.frame(year = 2006:2023, observed_anomaly = real_data$V6[157:174] - mean(real_data$V6[12:31])) # change to 157:174(start=2006,end=2023)
  other_quantiles <- other_estimates$observed_anomaly +
    matrix(rep(qnorm(c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)), each = 18), nrow= 18)# change to 18 because 2006 to 2023
  other_quantiles <- data.frame(other_quantiles)
  names(other_quantiles) <- c("q025", "q050", "q100", "median", "q900", "q950", "q975")
  other_estimates <- cbind(other_estimates, other_quantiles)
  summary_anomaly <- rbind(summary_anomaly, other_estimates)
  p1 <- ggplot(data = summary_anomaly, aes(x = year))  +
    theme(text=element_text(size=17),
          panel.background=element_rect(fill="white"),
          panel.grid.major=element_line(color="white") ,
          panel.grid.minor=element_line(color="white")) +
    geom_hline(yintercept = seq(-0.5,4.5,0.5), color = 'grey') +
    geom_vline(xintercept = c(2010, 2030, 2050, 2070, 2090, 2100), color = 'grey') +
    geom_line(aes(y = median, color = 'True'), size = 2) +
    geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    scale_color_manual(values = c('HadCrut5' = 'black', 'True' = 'red', 'Projection' = 'blue')) +
    geom_line(aes(y = observed_anomaly, color = 'HadCrut5')) +
    ylim(-0.5, 4.5)
  p1 <- p1 +
    #geom_ribbon(aes(ymin = q025, ymax = q975), fill = 'red', alpha = 0.2) +
    geom_ribbon(aes(ymin = q050, ymax = q950), fill = 'red', alpha = 0.3) +
    #geom_ribbon(aes(ymin = q025, ymax = q975), data = summary_proj,  fill = 'blue', alpha = 0.2) +
    geom_ribbon(aes(ymin = q050, ymax = q950, x = year), data = summary_proj, fill = 'blue', alpha = 0.3) +
    xlab('year') + ylab('Anomaly') + ggtitle(title)
  return (p1)

}

p1 <- plot_temp(summary_proj, summary_anomaly, 'Anomaly Unadjusted')
p2 <- plot_temp(summary_proj_adjusted, summary_anomaly, 'Anomaly Adjusted')
p3 <- plot_temp(summary_proj_adjusted_cont, summary_anomaly, 'Anomaly Adjusted and Policy Continued')
p4 <- plot_temp(summary_proj_adjusted_usa, summary_anomaly, 'Anomaly Adjusted without USA')


plot.paristarget <- function(Iso, objective, indeces, intensity = FALSE, newdata = NULL,
                             title, ipat.quantiles.bycountry, yaxis = TRUE)
{

  data.historical <- subset(data.medium, Isocode == Iso)
  color.us <- "#e41a1c"
  # browser()
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=15)) +
    scale_x_continuous(breaks=seq(2010,2100,by=20))

  if (!yaxis)
  {
    plot.obj <- plot.obj + theme(axis.title.y=element_blank(),
                                 axis.text.y=element_blank(),
                                 axis.ticks.y=element_blank())
  }
  # browser()
  data.plot <- data.frame(Year = ipat.quantiles.bycountry[[Iso]][["0.025"]]$Year[indeces])
  quantiles.need <- c(0.025, 0.05, 0.5, 0.95, 0.975)
  obj.data <- data.frame(Year = objective[1], obj = objective[2])
  if(intensity)
  {
    range_min <- min(data.medium$Tech[data.medium$Isocode == Iso],
                     proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry[[Iso]]$`0.025`$Tech,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry[[Iso]]$`0.025`$Tech,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont[[Iso]]$`0.025`$Tech)
    range_max <- max(data.medium$Tech[data.medium$Isocode == Iso],
                     proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry[[Iso]]$`0.975`$Tech,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry[[Iso]]$`0.975`$Tech,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont[[Iso]]$`0.975`$Tech)
    for(quant in quantiles.need)
    {
      data.plot[[paste0('quantiles', quant)]] <- ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$Tech[indeces]
    }
    plot.obj <- plot.obj +
      ylab('tonnes of CO2 per $10,000') +
      geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                  alpha=0.3, fill=color.us) +
      geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
                  alpha=0.2, fill=color.us) +
      geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                    color="Projections"),
                alpha=1, size=1.3, color=color.us) +
      ylim(range_min, range_max)
    plot.obj <- plot.obj + geom_line(data = data.historical, aes(x=Year, y = Tech), size = 1.2) +
      geom_point(aes(x = Year, y = obj), data = obj.data, color="Blue", size = 4)
  }
  else
  {
    # browser()
    range_min <- min(data.medium$CO2.total[data.medium$Isocode == Iso] * 11 / 3e3,
                     proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry[[Iso]]$`0.025`$CO2 / 1e9,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry[[Iso]]$`0.025`$CO2 / 1e9,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont[[Iso]]$`0.025`$CO2 / 1e9)
    range_max <- max(data.medium$CO2.total[data.medium$Isocode == Iso] * 11 / 3e3,
                     proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry[[Iso]]$`0.975`$CO2 / 1e9,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry[[Iso]]$`0.975`$CO2 / 1e9,
                     proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont[[Iso]]$`0.975`$CO2 / 1e9)

    data.historical$CO2Total <- data.historical$CO2 * (1000 * data.historical$Pop)/1e9
    # data.2023 <- data.frame(Year = 2023, CO2 = newdata)
    for(quant in quantiles.need)
    {
      data.plot[[paste0('quantiles', quant)]] <- ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$CO2[indeces]/1e9
    }
    plot.obj <- plot.obj +
      ylab('Total CO2 Emissions (gt)') +
      geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                  alpha=0.3, fill=color.us) +
      geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
                  alpha=0.2, fill=color.us) +
      geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                    color="Projections"),
                alpha=1, size=1.3, color=color.us) +
      ylim(range_min, range_max)
    plot.obj <- plot.obj + geom_line(data = data.historical, aes(x=Year, y = CO2Total), size = 1.2) +
      geom_point(aes(x = Year, y = obj), data = obj.data, color="Blue", size = 4) # +
    # geom_point(aes(x = Year, y = CO2), data = data.2023,  color="black", size = 4)
  }
  plot.obj <- plot.obj + ggtitle(title)
  return (plot.obj)
}

get.ipat.quantiles.bycountry <- function(ipat.components.bycountry, quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                                         na.rm=FALSE) {
  # get.ipat.quantiles.bycountry gives IPAT quantiles on a country level.
  # It's a variation on get.ipat.medians, despite the name change from "medians" to "quantiles"

  ipat.quants.bycountry <- list()
  for (iso in names(ipat.components.bycountry[[1]])) {
    ipat.quants.bycountry[[iso]] <- list()
    for (q in quantiles) {
      ipat.quants.bycountry[[iso]][[as.character(q)]] <- ipat.components.bycountry[[1]][[iso]]
      ipat.quants.bycountry[[iso]][[as.character(q)]][, -1] <- NA
    }

    row.index <- 0
    for (year in ipat.components.bycountry[[1]][[iso]]$Year) {
      row.index <- row.index + 1
      for (var.name in names(ipat.components.bycountry[[1]][[iso]])[-1]) {
        if (ipat.components.bycountry[[1]][[iso]]$Year[1] >= 1990) {
          quants.tmp <- quantile(sapply(1:length(ipat.components.bycountry), function(i)
            ipat.components.bycountry[[i]][[iso]][row.index, var.name]), probs=quantiles,
            na.rm=na.rm)
        } else {
          # Need to remove missing because of missing data for some countries
          quants.tmp <- quantile(sapply(1:length(ipat.components.bycountry), function(i)
            ipat.components.bycountry[[i]][[iso]][row.index, var.name]), probs=quantiles,
            na.rm=T)
        }
        for (q_ind in 1:length(quantiles)) {
          q <- quantiles[q_ind]
          ipat.quants.bycountry[[iso]][[as.character(q)]][row.index, var.name] <- quants.tmp[q_ind]
        }
      }
    }
  }

  ipat.quants.bycountry
}

get.ipat.cum.quantiles.bycountry <- function(ipat.components.bycountry,
                                             quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                                             na.rm=FALSE) {
  # get.ipat.quantiles.bycountry gives IPAT quantiles on a country level.
  # It's a variation on get.ipat.medians, despite the name change from "medians" to "quantiles"

  ipat.quants.bycountry <- list()
  for (iso in names(ipat.components.bycountry[[1]])) {
    # if (iso == 'USA') {browser()}
    ipat.quants.bycountry[[iso]] <- data.frame(matrix(0, nrow = nrow(ipat.components.bycountry[[1]][[iso]]),
                                               ncol = length(quantiles) + 1))
    names(ipat.quants.bycountry[[iso]]) <- c('Year', as.character(quantiles))
    ipat.quants.bycountry[[iso]][, 'Year'] <- ipat.components.bycountry[[1]][[iso]][, 'Year']

    var.name <- 'CO2'
    cumsum_co2 <- sapply(1:length(ipat.components.bycountry), function(i)
      cumsum(ipat.components.bycountry[[i]][[iso]][, var.name])) / 2 +
      rbind(0, sapply(1:length(ipat.components.bycountry), function(i)
        cumsum(ipat.components.bycountry[[i]][[iso]][-78, var.name]))) / 2 # change to 78 (2023 to 2100)
    ipat.quants.bycountry[[iso]][, 2:6] <- t(apply(cumsum_co2, 1, quantile, probs = quantiles))
  }

  ipat.quants.bycountry
}

get.ipat.cum.quantiles <- function(trajs.worldwide,
                                   quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975),
                                   na.rm=FALSE) {
  # get.ipat.quantiles.bycountry gives IPAT quantiles on a country level.
  # It's a variation on get.ipat.medians, despite the name change from "medians" to "quantiles"

  ipat.quants <- data.frame(matrix(0, nrow = nrow(trajs.worldwide[[1]]),
                                                  ncol = length(quantiles) + 1))
  names(ipat.quants) <- c('Year', as.character(quantiles))
  ipat.quants[, 'Year'] <- trajs.worldwide[[1]][, 'Year']

  var.name <- 'CO2'
  cumsum_co2 <- sapply(1:length(trajs.worldwide), function(i)
    cumsum(trajs.worldwide[[i]][, var.name])) / 2 +
    rbind(0, sapply(1:length(trajs.worldwide), function(i)
      cumsum(trajs.worldwide[[i]][-78, var.name]))) / 2 # change to 78
  ipat.quants[, 2:6] <- t(apply(cumsum_co2, 1, quantile, probs = quantiles)) / 1e9

  ipat.quants
}


ipat.quantiles.bycountry <- get.ipat.quantiles.bycountry(proj.evals.2023.adjusted$ipat.components.bycountry)
ipat.annual.quantiles.bycountry <-
  get.ipat.quantiles.bycountry(proj.evals.2023.adjusted$ipat.annual.components.bycountry)
ipat.annual.quantiles.bycountry.cont <-
  get.ipat.quantiles.bycountry(proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont)

proj.evals.2023.adjusted$ipat.quantiles.bycountry <- ipat.quantiles.bycountry
proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry <- ipat.annual.quantiles.bycountry
proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont <- ipat.annual.quantiles.bycountry.cont

proj.evals.2023.ar1.const$ipat.quantiles.cum.bycountry <-
  get.ipat.cum.quantiles.bycountry(proj.evals.2023.ar1.const$ipat.annual.components.bycountry)
proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry <-
  get.ipat.cum.quantiles.bycountry(proj.evals.2023.adjusted$ipat.annual.components.bycountry)
proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry.cont <-
  get.ipat.cum.quantiles.bycountry(proj.evals.2023.adjusted$ipat.annual.components.bycountry.cont)

save(proj.evals.2023.adjusted, file = 'proj.evals.2023.adjusted.new.Rda')
font <- 'Times'
plot.co2.projections <- function(data.medium,
                                 trajs.quants,
                                 ylabel="Yearly CO2 emissions (gt CO2)",
                                 plot.rcp=TRUE,
                                 ybreaks=seq(0,120,by=20)) {
  # plot.co2.projections plots projections of yearly CO2 emissions worldwide

  library(doBy)
  library(ggplot2)
  data.restriction <- data.medium[!is.na(data.medium$CO2),]

  co2.sums.na.rm <- summaryBy(CO2.total~ Year, data.restriction, FUN=function(x) sum(x, na.rm=T))
  names(co2.sums.na.rm)[2] <- "Carbon"
  co2.sums.na.rm$Carbon <- co2.sums.na.rm$Carbon / 1e3 * 11/3

  year.start <- 2010
  year.end <- 2100
  x.vals <- seq(year.start, year.end, by=5)

  # Create plot
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=22)) +
    scale_y_continuous(breaks=ybreaks) +
    scale_x_continuous(breaks=seq(1960,2100,by=20))

  trajs.quants.tmp <- t(trajs.quants)
  trajs.quants.transp <- data.frame(trajs.quants.tmp[-1, ])
  names(trajs.quants.transp) <- paste0("Quant", trajs.quants.tmp[1,])
  # browser()
  trajs.quants.transp$Year <- 2023:2100

  # Plot our projections
  color.us <- "#e41a1c"
  plot.obj <- plot.obj +
    ylab(ylabel) +
    geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=trajs.quants.transp,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=Quant0.025, ymax=Quant0.975), data=trajs.quants.transp,
                alpha=0.2, fill=color.us) +
    geom_line(data=trajs.quants.transp, aes(x=Year, y=Quant0.5,
                                            color="Projections"),
              alpha=1, size=1.3, color=color.us)

  if (plot.rcp) {
    # Plot RCP numbers
    years.tmp <- seq(2010, 2100, by=10)
    rcp.yearly.tmp <- data.frame(sapply(data.frame(t(rcp.carbon.yearly)[-c(1:4),],
                                                   stringsAsFactors=F),
                                        as.numeric))
    names(rcp.yearly.tmp) <- rcp.carbon.yearly$Scenario
    rcp.yearly.tmp$Year <- years.tmp
    library(reshape2)
    rcp.yearly.long <- melt(rcp.yearly.tmp, id.vars=c("Year"))
    plot.obj <- plot.obj + geom_line(data=rcp.yearly.long,
                                     aes(x=Year, y=value,
                                         color=variable), size=0.7, alpha=0.9,
                                     linetype=2) +
      scale_color_manual(limits=c("SSP5-8.5", "SSP3-7.0", "SSP2-4.5", "SSP1-2.6"),
                         values=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00"))
  }

  # Plot historical data
  plot.obj <- plot.obj + geom_line(data=co2.sums.na.rm,
                                   aes(x=Year,y=Carbon))

  plot.obj
}

# plot.co2.projections(data.medium, proj.evals.2023.adjusted$trajs.annual.quants)
plot.co2.projections(data.medium, proj.evals.2023.ar1.const$trajs.annual.quants) # change
proj.evals.2023.adjusted$ipat.quantiles.bycountry <- ipat.quantiles.bycountry
proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry <- ipat.annual.quantiles.bycountry
proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont <- ipat.annual.quantiles.bycountry.cont

## Plot Yearly Emission by Country
plot.location <- 'Paris/'
font <- 'Times'
for (country_name in proj.evals.2023.adjusted$adjusted_countries)
{
  if (country_name == 'DMA') {next}
  paris.obj.country <- paris.objective[paris.objective$country_code == country_name, ]
  if (paris.obj.country$Type == 'Intensity')
  {
    obj <- c(paris.obj.country$year, (1 - paris.obj.country$size / 100) * subset(data.medium, Year == paris.obj.country$old_year & Isocode == country_name)$Tech)
    p <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = T,  # change to 78
                          title = paste0("Adjusted"),
                          ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry, yaxis = F)
    p1 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = T, # change to 78
                           title = paste0("Continued"),
                           ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont, yaxis = F)
    p2 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = T,# change to 78
                           title = paste0("None"),
                           ipat.quantiles.bycountry = proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry)
  }
  else if (paris.obj.country$Type == 'Emission')
  {
    base <- subset(data.medium, Year == paris.obj.country$old_year & Isocode == country_name)$CO2 *
      (1000 * subset(data.medium, Year == paris.obj.country$old_year & Isocode == country_name)$Pop)
    obj <- c(paris.obj.country$year, (1 - paris.obj.country$size / 100) * base / 1e9)
    p <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                          title = paste0('Adjusted'),
                          ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry, yaxis = F)
    p1 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                           title = paste0("Continued"),
                           ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont, yaxis = F)
    p2 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                           title = paste0("None"),
                           ipat.quantiles.bycountry = proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry)
  }
  else
  {
    if (!is.na(paris.obj.country$ratio))
    {
      base <- subset(data.medium, Year == paris.obj.country$old_year & Isocode == country_name)$CO2 *
        (1000 * subset(data.medium, Year == paris.obj.country$old_year & Isocode == country_name)$Pop)
      obj <- c(paris.obj.country$year, (1 - paris.obj.country$ratio / 100) * base / 1e9)
      p <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                            title = paste0(paris.obj.country$Country, "'s Emission"),
                            ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry, yaxis = F)
      p1 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                             title = paste0(paris.obj.country$Country, "'s Emission"),
                             ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont, yaxis = F)
      p2 <- plot.paristarget(country_name, obj, indeces = 1:78, intensity = F, # change
                             title = paste0(paris.obj.country$Country, "'s Emission"),
                             ipat.quantiles.bycountry = proj.evals.2023.ar1.const$ipat.annual.quantiles.bycountry)
    }
  }
  pdf(file = paste0('Paris/Comparison/',
                    paris.obj.country$Country, ".pdf"), width=10, height=5)
  grid.arrange(p2, p, p1, ncol=3)
  dev.off()
}

plot.cum.by.country <- function(Iso, indeces = 1:78, title, ipat.quantiles.bycountry) # change
{
  color.us <- "#e41a1c"
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=15)) +
    scale_x_continuous(breaks=seq(2010,2100,by=20))

  # browser()
  data.plot <- data.frame(Year = ipat.quantiles.bycountry[[Iso]]$Year[indeces])
  quantiles.need <- c(0.025, 0.05, 0.5, 0.95, 0.975)
  range_min <- min(proj.evals.2023.ar1.const$ipat.quantiles.cum.bycountry[[Iso]]$`0.025` / 1e9,
                   proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry[[Iso]]$`0.025` / 1e9,
                   proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry.cont[[Iso]]$`0.025` / 1e9)
  range_max <- max(proj.evals.2023.ar1.const$ipat.quantiles.cum.bycountry[[Iso]]$`0.975` / 1e9,
                   proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry[[Iso]]$`0.975` / 1e9,
                   proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry.cont[[Iso]]$`0.975` / 1e9)

  # data.2023 <- data.frame(Year = 2023, CO2 = newdata)
  for(quant in quantiles.need)
  {
    data.plot[[paste0('quantiles', quant)]] <-
      ipat.quantiles.bycountry[[Iso]][[as.character(quant)]][indeces]/1e9
  }
  plot.obj <- plot.obj +
    ylab('Cumulative CO2 Emissions (gt)') +
    geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                alpha=0.3, fill=color.us) +
    geom_ribbon(aes(x=Year, ymin=quantiles0.025, ymax=quantiles0.975), data=data.plot,
                alpha=0.2, fill=color.us) +
    geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                  color="Projections"),
              alpha=1, size=1.3, color=color.us) +
    ylim(range_min, range_max)
  plot.obj <- plot.obj + ggtitle(title)
  return (plot.obj)
}

for (country_name in proj.evals.2023.adjusted$adjusted_countries)
{
  if (country_name == 'DMA') {next}
  p <- plot.cum.by.country(country_name,
                           title = paste0(paris.objective$Country[paris.objective$country_code == country_name], "'s Emission"),
                           ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry)
  p1 <- plot.cum.by.country(country_name,
                            title = paste0(paris.objective$Country[paris.objective$country_code == country_name], "'s Emission"),
                           ipat.quantiles.bycountry = proj.evals.2023.adjusted$ipat.quantiles.cum.bycountry.cont)
  p2 <- plot.cum.by.country(country_name,
                            title = paste0(paris.objective$Country[paris.objective$country_code == country_name], "'s Emission"),
                           ipat.quantiles.bycountry = proj.evals.2023.ar1.const$ipat.quantiles.cum.bycountry)

  pdf(file = paste0(plot.location, 'Cumulative/',
                    paris.objective$Country[paris.objective$country_code == country_name], ".pdf"),
      width=10, height=5)
  grid.arrange(p2, p, p1, ncol=3)
  cat(country_name, '\n')
  dev.off()
}


proj.evals.2023.adjusted$trajs.cum.quants <- get.ipat.cum.quantiles(proj.evals.2023.adjusted$trajs.annual.worldwide)
proj.evals.2023.adjusted$trajs.cum.quants.cont <- get.ipat.cum.quantiles(proj.evals.2023.adjusted$trajs.annual.worldwide.cont)
proj.evals.2023.ar1.const$trajs.cum.quants <- get.ipat.cum.quantiles(proj.evals.2023.ar1.const$trajs.annual.worldwide)

plot.cum.emission <- function(trajs.cum.quants, title)
{
  names(trajs.cum.quants)[2:6] <- paste0('q', names(trajs.cum.quants)[2:6])
  p1 <- ggplot(data = trajs.cum.quants, aes(x = Year))  +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=17)) +
    geom_line(aes(y = q0.5, color = 'red'), size = 2) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = 'red', alpha = 0.2) +
    geom_ribbon(aes(ymin = q0.05, ymax = q0.95), fill = 'red', alpha = 0.3) +
    xlab('year') + ylab('Cumulative Emission') +
    ggtitle(title)

  p1
}
p1 <- plot.cum.emission(proj.evals.2023.ar1.const$trajs.cum.quants, 'Cumulative Emission')
p2 <- plot.cum.emission(proj.evals.2023.adjusted$trajs.cum.quants, 'Cumulative Emission with Paris Agreement')
p3 <- plot.cum.emission(proj.evals.2023.adjusted$trajs.cum.quants.cont, 'Cumulative Emission with Policy Continued')


save(proj.evals.2023.adjusted, file = 'proj_evals_2023_adjusted.Rda')

#### What's needed?
# bi_section <- function(objective = 2, prob = 0.5)
# {
#   obj <- 2
#   ratio_co2_r <- ratio_co2 <- 1
#   ratio_co2_l <- 0
#   count <- 1
#   while (abs(obj) > 1e-3)
#   {
#     ratio_co2 <- ratio_co2_l / 2 + ratio_co2_r / 2
#     proj_list <- project_temperature(fit_result, proj.evals.2023.adjusted$trajs.annual.worldwide.cont, rcp.carbon.cum.complete, rcp_temp_data,
#                                      year.range = c(2023, 2100), n.traj = 1000, ratio_co2 = ratio_co2)
#     proj_temp <- proj_list$projection
#     summary_proj <- data.frame(year = proj_temp[,1])
#     # browser()
#     summary_proj$quantile <- apply(proj_temp[,-1], 1, quantile, prob = prob) - (mean(cleaned_all_data$Average_Diff$modmean[11:30]) + 273.15)
#     global_warming <- summary_proj$quantile[summary_proj$year == 2100]
#     obj <- global_warming - objective
#     if (obj > 0)
#     {
#       ratio_co2_r <- ratio_co2
#     }
#     else
#     {
#       ratio_co2_l <- ratio_co2
#     }
#     cat('Iteration: ', count, ' Objective: ', obj, 'Ratio: ', ratio_co2, '\n')
#     count <- count + 1
#   }
#   return (ratio_co2)
# }
# ratio_2 <- bi_section()
# ratio_2_95 <- bi_section(prob = 0.95)
# ratio_1_5 <- bi_section(objective = 1.5)
# obj_co2 <- ((-(mean(cleaned_all_data$Average_Diff$modmean[11:30]) + 273.15) - 1.5 + rev(rcp_temp_data$rcp45)[1]) * rcp.carbon.cum.complete$Carbon2100[3] +
#               ((mean(cleaned_all_data$Average_Diff$modmean[11:30]) + 273.15) + 1.5 - rev(rcp_temp_data$rcp26)[1]) * rcp.carbon.cum.complete$Carbon2100[2]) /
#   ((-(mean(cleaned_all_data$Average_Diff$modmean[11:30]) + 273.15) - 1.5 + rev(rcp_temp_data$rcp45)[1]) +
#      ((mean(cleaned_all_data$Average_Diff$modmean[11:30]) + 273.15) + 1.5 - rev(rcp_temp_data$rcp26)[1]))
#
#
# cumulative_co2 <- sapply(1:1000, function(i) {return (sum(proj.evals.2023.adjusted$trajs.annual.worldwide.cont[[i]]$CO2/1e9))})
# obj_co2 / median(cumulative_co2)
##########################

plot.intensity <- function(Iso, version = 2010)
{
  color.us <- "#e41a1c"
  plot.obj <- ggplot() +
    labs(color="Scenario") +
    theme(panel.background=element_rect(fill="white"),
          panel.grid.minor=element_line(color="white"),
          panel.grid.major=element_line(color="grey90"),
          text=element_text(family=font, size=15)) +
    scale_x_continuous(breaks=seq(2010,2100,by=20))

  # browser()
  quantiles.need <- c(0.05, 0.5, 0.95)
  range_min <- min(proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]]$`0.05`$Tech,
                   proj.evals.2023.ar1.const$ipat.quantiles.bycountry[[Iso]]$`0.05`$Tech,
                   data.medium$Tech[data.medium$Isocode == Iso],
                   data.medium.new$Tech[data.medium.new$Isocode == Iso])
  range_max <- max(proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]]$`0.95`$Tech,
                   proj.evals.2023.ar1.const$ipat.quantiles.bycountry[[Iso]]$`0.95`$Tech,
                   data.medium$Tech[data.medium$Isocode == Iso],
                   data.medium.new$Tech[data.medium.new$Isocode == Iso])

  if (version == 2010)
  {
    data.plot <- data.frame(Year = seq(2010, 2100, 5))
  }
  else
  {
    data.plot <- data.frame(Year = seq(2023, 2100, 1))
  }


  # data.2023 <- data.frame(Year = 2023, CO2 = newdata)
  for(quant in quantiles.need)
  {
    if (version == 2010)
    {
      data.plot[[paste0('quantiles', quant)]] <-
        proj.evals.2010.ar1.const$ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$Tech
    }
    else
    {
      data.plot[[paste0('quantiles', quant)]] <-
        proj.evals.2023.ar1.const$ipat.quantiles.bycountry[[Iso]][[as.character(quant)]]$Tech
    }
  }
  plot.obj <- plot.obj +
    ylab('Intensity') +
    geom_ribbon(aes(x=Year, ymin=quantiles0.05, ymax=quantiles0.95), data=data.plot,
                alpha=0.3, fill=color.us) +
    geom_line(data=data.plot, aes(x=Year, y=quantiles0.5,
                                  color="Projections"),
              alpha=1, size=1.3, color=color.us)
  if (version == 2010)
  {
    plot.obj <- plot.obj +
      geom_line(data = data.medium[data.medium$Isocode == Iso,], aes(x = Year, y = Tech), lwd = 2) +
      ylim(range_min, range_max)
  }
  else
  {
    plot.obj <- plot.obj +
      geom_line(data = data.medium.new[data.medium.new$Isocode == Iso,], aes(x = Year, y = Tech), lwd = 2) +
      ylim(range_min, range_max)
  }

  plot.obj <- plot.obj + ggtitle(paste0(countrycode(Iso, origin = 'iso3c', destination = 'country.name'),
                                        "'s Intensity Version ", version))
  return (plot.obj)
}
## can not find 'proj.evals.2010.ar1.const'
# for (Iso in unique(data.medium$Isocode))
# {
#   p <- plot.intensity(Iso, version = 2010)
#   p1 <- plot.intensity(Iso, version = 2023)
#
#   pdf(file = paste0(plot.location, 'Intensity/',
#                     countrycode(Iso, origin = 'iso3c', destination = 'country.name'), ".pdf"),
#       width=10, height=5)
#   grid.arrange(p, p1, ncol=2)
#   cat(countrycode(Iso, origin = 'iso3c', destination = 'country.name'), '\n')
#   dev.off()
# }

ratio_emission <- data.frame(Iso = names(proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry))
ratio_emission$country_name <- countrycode(ratio_emission$Iso, origin = 'iso3c', destination = 'country.name')
ratio_emission$proportion <- NA
total <- proj.evals.2023.adjusted$trajs.annual.quants.cont$CO22100[3]
for (i in 1:nrow(ratio_emission))
{
  iso <- as.character(ratio_emission$Iso[i])
  proportion <- rev(proj.evals.2023.adjusted$ipat.annual.quantiles.bycountry.cont[[iso]]$`0.5`$CO2)[1]/1e9
  ratio_emission$proportion[i] <- proportion / total
}

ratio_emission1 <- ratio_emission
ratio_emission2 <- ratio_emission

plot(c(0, 7000), c(287, 292), type = 'n', xlab = 'cumulative emission', ylab = 'CMIP6 Temperature',
     main = 'Temperature and Cumulative Emission from SSP Scenarios')
points(rcp.carbon.cum.complete.adjusted$ssp126, rcp_temp_data$rcp26, pch = 16, col = 'red', cex = 0.3)
points(rcp.carbon.cum.complete.adjusted$ssp245, rcp_temp_data$rcp45, pch = 16, col = 'blue', cex = 0.3)
points(rcp.carbon.cum.complete.adjusted$ssp370, rcp_temp_data$rcp60, pch = 16, col = 'green', cex = 0.3)
points(rcp.carbon.cum.complete.adjusted$ssp585, rcp_temp_data$rcp85, pch = 16, col = 'purple', cex = 0.3)
abline(lm(c(rcp_temp_data$rcp26, rcp_temp_data$rcp45, rcp_temp_data$rcp60, rcp_temp_data$rcp85) ~
            c(rcp.carbon.cum.complete.adjusted$ssp126, rcp.carbon.cum.complete.adjusted$ssp245,
              rcp.carbon.cum.complete.adjusted$ssp370, rcp.carbon.cum.complete.adjusted$ssp585)))
legend('topleft', legend = c('SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5'),
       pch = 16, col = c('red', 'blue', 'green', 'purple'))

all_fits <- list()
all_projections <- list()

real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[132:156])


for (i in 1:54)
{
  model_x <- all_models[i]
  fit_result <- fit_ar1_model(all_hist_df, model_x, n.iterations = 10000) # fit_result

  proj_list <- project_temperature(fit_result, proj.evals.2023.ar1.const$trajs.annual.worldwide, rcp.carbon.cum.complete,
                                   model_x, year.range = c(2023, 2100), n.traj = 1000)
  all_fits[[model_x]] <- fit_result
  all_projections[[model_x]] <- proj_list


}

combined_projections <- c()
for(i in 1:54)
{
  model_x <- all_models[i]
  combined_projections <- cbind(combined_projections, all_projections[[model_x]]$xt_projection)
}

combined.results <- t(apply(combined_projections, 1, quantile, probs=c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)))
combined.results <- as.data.frame(combined.results)
combined.results$Year <- 2024:2100 # change

names(combined.results)[1:7] <- paste0('Quant', c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

combined.results[, 1:7] <- combined.results[, 1:7] + mean(real_data$V6[132:156]) - mean(real_data$V6[12:31])

combined.results.paris <- combined.results
combined.results.cont <- combined.results

## 4.2
plots <- list()

font <- 'Times'
plot.obj <- ggplot() +
  labs(color="Scenario") +
  theme(panel.background=element_rect(fill="white"),
        panel.grid.minor=element_line(color="white"),
        panel.grid.major=element_line(color="grey90"),
        text=element_text(family=font, size=12))

color.us <- "#e41a1c"
plot.obj <- plot.obj + ylab('Anomaly') +
  geom_ribbon(aes(x=Year, ymin=Quant0.05, ymax=Quant0.95), data=combined.results,
              alpha=0.2, fill=color.us) +
  geom_ribbon(aes(x=Year, ymin=Quant0.1, ymax=Quant0.9), data=combined.results,
              alpha=0.3, fill=color.us) +
  geom_line(data=combined.results, aes(x=Year, y=Quant0.5),
            alpha=1, size=1.3, color=color.us) +
  # scale_x_continuous(breaks=c(2010, seq(2025, year.end, 25)), expand = c(0,0)) +
  scale_y_continuous(breaks = c(-0.5,0:4, 4.5), limits = c(-0.5, 4.5))

real_data_adjusted <- real_data
real_data_adjusted$V6 <- real_data_adjusted$V6 - mean(real_data_adjusted$V6[11:32])

plot.obj <- plot.obj +
  geom_line(aes(x = V1, y = V6), data = real_data_adjusted)

plot.obj <- plot.obj + ggtitle('Current Forecast')
plots[[1]] <- plot.obj

#grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=3) # only have plots[[1]]
#plots[[1]]
combined <- summary_proj
combined$adjusted <- 'None'
test <- summary_proj_adjusted
test$adjusted <- 'Adjusted'
combined <- rbind(combined, test)
test <- summary_proj_adjusted_cont
test$adjusted <- 'Continued'
combined <- rbind(combined, test)
test <- summary_proj_adjusted_usa
test$adjusted <- 'USA Excluded'
combined <- rbind(combined, test)

## plots for 4 scenarios
plot_temp_new <- function(summary_proj, summary_anomaly, title = 'Anomaly Forecast')
{
  # summary_proj <- rbind(summary_anomaly[145, 1:8], summary_proj)
  other_estimates <- data.frame(year = 2006:2023, observed_anomaly = real_data$V6[157:174] - mean(real_data$V6[12:31])) # change to 174 due to 2023
  other_quantiles <- other_estimates$observed_anomaly +
    matrix(rep(qnorm(c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)), each = 18), nrow=18) # change to 18
  other_quantiles <- data.frame(other_quantiles)
  names(other_quantiles) <- c("q025", "q050", "q100", "median", "q900", "q950", "q975")
  other_estimates <- cbind(other_estimates, other_quantiles)
  summary_anomaly <- rbind(summary_anomaly, other_estimates)
  p1 <- ggplot(data = summary_proj, aes(x = year))  +
    theme(text=element_text(size=17),
          panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(color="white") ,
          panel.grid.minor = element_line(color="white"),
          axis.text.x = element_text(angle = 45),
          legend.title = element_blank()) +
    geom_hline(yintercept = seq(-0.5,4.5,0.5), color = 'grey') +
    geom_vline(xintercept = c(1900, 2000, 2050, 2100), color = 'grey') +
    geom_line(data = summary_anomaly, aes(x = year, y = median, color = 'Estimate'), size = 2) +
    geom_line(mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    # geom_line(data = summary_proj, mapping = aes(x = year, y = median, color = 'Projection')) +
    scale_color_manual(values = c('HadCrut5' = 'black', 'Estimate' = 'red', 'Projection' = 'blue')) +
    geom_line(data = summary_anomaly, aes(x = year, y = observed_anomaly, color = 'HadCrut5')) +
    ylim(-0.5, 4.5)
  p1 <- p1 +
    #geom_ribbon(aes(ymin = q025, ymax = q975), fill = 'red', alpha = 0.2) +
    geom_ribbon(aes(ymin = q050, ymax = q950), data = summary_anomaly[1:145,], fill = 'red', alpha = 0.3) +
    #geom_ribbon(aes(ymin = q025, ymax = q975), data = summary_proj,  fill = 'blue', alpha = 0.2) +
    geom_ribbon(aes(ymin = q050, ymax = q950, x = year), fill = 'blue', alpha = 0.3) +
    xlab('year') + ylab('Anomaly') + ggtitle(title)
  p1 <- p1 + facet_wrap(~ adjusted, ncol = 2)
  return (p1)

}

plot_temp_new(combined, summary_anomaly) # interval problem
test_df <- data.frame(year = rep(1:10, 3), val = rnorm(30), group = rep(1:3, each=10))
test_df$group <- factor(test_df$group)
head(test_df)

ggplot(data= test_df) + geom_line(aes(x=year, y =val, col = group))

#### Extra percentage
calculate_a <- function(X)
{
  c_func <-function(a)
  {
    return ((1-exp(-a * 85))/a - X/36.6)
  }
  library(pracma)
  root <- bisect(fun = c_func, 0.00001,0.3)$root
  return (root)
}

X <- c(2083, 1579, 839, 471)
rates <- numeric(4)
for(i in 1:4) {rates[i] <- calculate_a(X[i])}

extras <- matrix(nrow=10, ncol=3)
count <- 1
for (country in c('USA', 'JPN', 'DEU', 'RUS', 'CAN', 'KOR', 'BRA', 'GBR'))
{
  l <- which(paris.objective$country_code == country)
  old_year <- paris.objective$old_year[l]
  type <- paris.objective$Type[l]
  size <- ifelse(type == 'Emission', paris.objective$size[l], size <- paris.objective$ratio[l])
  year_ind <- paris.objective$year[l] - - 2022 # change 2014 here, forget to change last time
  ref_level <- data.medium$CO2.total[(data.medium$Isocode == country) &
                                       (data.medium$Year == old_year)] * (1 - size / 100)
  rate <- ref_level / data.medium$CO2.total[(data.medium$Isocode == country) &
                                              (data.medium$Year == 2023)]
  rate <- 1 - rate^(1/(paris.objective$year[l] - 2023))
  for(j in 1:3)
  {
    obj <- data.medium$CO2.total[(data.medium$Isocode == country) &
                                   (data.medium$Year == 2023)] * (1 - rate * rates[j+1]/rates[1])^(paris.objective$year[l] - 2023)
    extras[count, j] <- (1-obj/data.medium$CO2.total[(data.medium$Isocode == country) &
                                                       (data.medium$Year == old_year)])*100/size
  }
  count <- count + 1

}
extras[9,1] <- (1-(1-((1-(308.1 * 0.4/141.7)^(1/15)) * 1.802)) ^15*141.7/308.1)/0.6
extras[9,2] <- (1-(1-((1-(308.1 * 0.4/141.7)^(1/15)) * rates[3]/rates[1])) ^15*141.7/308.1)/0.6
extras[9,3] <- (1-(1-((1-(308.1 * 0.4/141.7)^(1/15)) * rates[4]/rates[1])) ^15*141.7/308.1)/0.6
extras[10,1] <- (1-(1-((1-((333.4/59) * 0.67/(633.3/113.7))^(1/15)) * 1.802)) ^15*(633.3/113.7)/(333.4/59))/0.33
extras[10,2] <- (1-(1-((1-((333.4/59) * 0.67/(633.3/113.7))^(1/15)) * rates[3]/rates[1])) ^15*(633.3/113.7)/(333.4/59))/0.33
extras[10,3] <- (1-(1-((1-((333.4/59) * 0.67/(633.3/113.7))^(1/15)) * rates[4]/rates[1])) ^15*(633.3/113.7)/(333.4/59))/0.33
numeric_columns <- sapply(extras, is.numeric)
extras[numeric_columns] <- (extras[numeric_columns] - 1) * 100
extras <- as.data.frame(extras)
extras$country <- c(c('USA', 'JPN', 'DEU', 'RUS', 'CAN', 'KOR', 'BRA', 'GBR'), 'CHN', 'IND')
xtable(extras)
