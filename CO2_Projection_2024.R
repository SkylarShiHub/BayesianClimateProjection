
# Things to run before any session:

# miscellaneous
library(gdata)
library(plyr)
library(countrycode)
library(doBy)
library(reshape2)

# Population prediction libraries
library(wpp2017)
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

library(car)

#=========================================================
## "peak within last 5 year"
find_maxes <- function(data.medium, plot.countries=FALSE) {
  # find_maxes finds when each country peaks in carbon intensity.
  # If the peak for a country is within the last 5 years of the data,
  # we list that country under rejects.late and don't report a peak
  # for it.
  # If there are fewer than 20 years of carbon intensity (Tech) data
  # for a country, we list that country under rejects.insuf and don't
  # report a peak for it.
  max.vals <- list()
  max.vals[["Iso"]] <- c()
  max.vals[["GDP"]] <- c()
  max.vals[["Tech"]] <- c()
  max.vals[["CO2"]] <- c()
  max.vals[["Year"]] <- c()
  rejects.early <- c()
  rejects.late <- c()
  rejects.insuf <- c()
  isolist <- unique(data.medium$Isocode)

  if (plot.countries) i <- 1
  for (iso in isolist) {
    data.tmp <- data.medium[data.medium$Isocode == iso,]
    data.tmp <- data.tmp[!is.na(data.tmp$Tech),]
    if (sum(!is.na(data.tmp$Tech)) < 20) {
      rejects.insuf <- c(rejects.insuf, iso)
      next
    }
    model.smoothed <- loess(Tech ~ Year, data=data.tmp, span=0.25)
    predicted <- predict(model.smoothed)
    max.ind <- which(predicted == max(predicted))
    # stopifnot(length(max.ind) == 1)
    max.ind <- min(max.ind)
    # Record unsmoothed observations at this point, if we're not near
    # the start or end of the period
    year.last <- max(data.medium$Year) # Allow for out-of-sample validation
    year.last.cutoff <- min(2003, year.last - 5)
    if (1965 <= data.tmp$Year[max.ind] &
        data.tmp$Year[max.ind] <= year.last.cutoff) {
      max.vals[["Iso"]] <- c(max.vals[["Iso"]], iso)
      max.vals[["GDP"]] <- c(max.vals[["GDP"]], data.tmp$GDP[max.ind])
      max.vals[["Tech"]] <- c(max.vals[["Tech"]], data.tmp$Tech[max.ind])
      max.vals[["CO2"]] <- c(max.vals[["CO2"]], data.tmp$CO2[max.ind])
      max.vals[["Year"]] <- c(max.vals[["Year"]], data.tmp$Year[max.ind])
    } else if (data.tmp$Year[max.ind] <= 1964) {
      rejects.early <- c(rejects.early, iso)
    } else {
      rejects.late <- c(rejects.late, iso)
    }

    # Create and save plot
    if (plot.countries) {
      n.perpage <- 8
      if (i %% n.perpage == 1) {
        if (i != 1) dev.off()
        png(paste0(plot.location, "max_plots", floor(i/n.perpage), ".png"),
            width=800, height=1200)
        par(mfrow=c(4,2))
      }
      plot(data.tmp$Year, data.tmp$Tech,
           main=paste0(getCountry(iso), " (", iso, ")"))
      points(data.tmp$Year[!is.na(data.tmp$Tech)], predicted, type="l")
      i <- i + 1
    }
  }
  if (plot.countries) dev.off()

  print("Rejects because peak was before 1965:")
  print(rejects.early)
  print(length(rejects.early))
  print(paste(getCountry(rejects.early), collapse=", "))
  print(paste0("Rejects because peak was after ", min(2003, year.last - 5)))
  print(rejects.late)
  print(length(rejects.late))
  print(paste(getCountry(rejects.late), collapse=", "))
  print("Rejects because of insufficient data:")
  print(rejects.insuf)
  print(length(rejects.insuf))
  if (length(rejects.insuf) > 0)
  {
    print(paste(getCountry(rejects.insuf), collapse=", "))
  }
  else
  {
    print('None.\n')
  }


  list(max.vals, rejects.early, rejects.late, rejects.insuf)
}

#==========================================


predict.population <- function(year.present=2023, year.end=2100,
                               n.iter=10000, n.burnin=1000,
                               make.new=TRUE) {
  # predict.population uses the packages bayesPop, bayesTFR,
  # and bayesLife to create population projections by country,
  # saving results into sims.location.
  # Set make.new=FALSE if you've already run predict.population
  # for the years of interest to load previous results.
  library(bayesPop)
  if (make.new) {
    library(bayesTFR)
    library(bayesLife)
    sim.dir.tfr <- paste0(sims.location, "Sims_TFR_", year.present, "_", year.end)
    # m.tfr <- run.tfr.mcmc(present.year=year.present,
    #                       iter=n.iter,
    #                       burnin=n.burnin,
    #                       output.dir=sim.dir.tfr, verbose=FALSE,
    #                       replace.output=TRUE)
    # if (year.present != 1980) {
    #   m3.tfr <- run.tfr3.mcmc(present.year=year.present,
    #                           sim.dir=sim.dir.tfr,
    #                           iter=n.iter,
    #                           thin=1,
    #                           burnin=n.burnin,
    #                           verbose=TRUE, replace.output=TRUE)
    #   pred.tfr <- tfr.predict(end.year=year.end,
    #                           sim.dir=sim.dir.tfr, m.tfr,
    #                           burnin=n.burnin,
    #                           use.tfr3=FALSE,
    #                           verbose=TRUE,
    #                           replace.output=TRUE)
    # } else {
    #   pred.tfr <- tfr.predict(end.year=year.end,
    #                           sim.dir=sim.dir.tfr, m.tfr,
    #                           burnin=n.burnin, burnin3=n.burnin,
    #                           verbose=TRUE,
    #                           replace.output=TRUE)
    # }

    sim.dir.e0 <- paste0(sims.location, "Sims_E0_", year.present, "_", year.end)
    # m.e0 <- run.e0.mcmc(present.year=year.present,
    #                     my.e0.file=paste0(sims.location, "aids_include.txt"),
    #                     iter=n.iter,
    #                     thin=1,
    #                     verbose=FALSE,
    #                     output.dir=sim.dir.e0,
    #                     replace.output=TRUE)

    # if (year.present < 2000) {
    #   pred.e0 <- e0.predict(end.year=year.end,
    #                         max.e0.eq1 = 78, max.e0.eq1.pred=78,
    #                         sim.dir=sim.dir.e0, m.e0, burnin=n.burnin, verbose=FALSE)
    # } else {
    #   pred.e0 <- e0.predict(end.year=year.end,
    #                         sim.dir=sim.dir.e0, m.e0, burnin=n.burnin, verbose=FALSE)
    # }
    # summary(pred.e0)

    sim.dir.pop <- paste0(sims.location, "Sims_Pop_", year.present, "_", year.end)
    # make changes
    pred.pop <- pop.predict(output.dir=sim.dir.pop, nr.traj = 1000,
                            keep.vital.events = FALSE,
                            end.year = year.end, start.year = 1950, wpp.year = 2022,
                            present.year = year.present,
                            verbose=TRUE, replace.output=TRUE,
                            inputs = list(tfr.sim.dir=sim.dir.tfr, e0F.sim.dir=sim.dir.e0,
                                          e0M.sim.dir='joint_'),
                            annual = TRUE
    )

    # pred.pop <- pop.predict(present.year=year.present, end.year=year.end,
    #                         output.dir=sim.dir.pop, verbose=FALSE,
    #                         inputs = list(tfr.sim.dir=sim.dir.tfr,
    #                                       e0F.sim.dir=sim.dir.e0, e0M.sim.dir='joint_'),
    #                         replace.output=TRUE)
    summary(pred.pop)
  } else {
    # Loading an already-made file
    sim.dir.pop <- paste0(sims.location, "Sims_Pop_", year.present, "_", year.end)
    pred.pop <- get.pop.prediction(sim.dir.pop)
  }

  # 9 for pred.pop$quantiles[,9,] is the 0.5 quantile
  preds.countries <- data.frame(cbind(as.character(pred.pop$countries[,2]),
                                      pred.pop$quantiles[,9,]),
                                stringsAsFactors=F)
  names(preds.countries) <- c("Country", paste0("Pop", seq(year.present, year.end, by=1))) #by=5
  # delete the year 2021 and 2022
  library(dplyr)
  preds.countries <- select(preds.countries, -c(Pop2021, Pop2022)) # may be changed for felxibale later


  library(countrycode)
  preds.countries$Isocode <- getIso(preds.countries$Country)

  list(preds.countries=preds.countries, pred.pop=pred.pop)
}



find.pop.trajectories <- function(data.medium, pred.pop, n.trajs, isos,
                                  n.iter.predict.population=1000,
                                  year.start=2010, year.end=2100,
                                  make.new=TRUE) {
  # find.pop.trajectories converts the format of the results from
  # predict.population to a format we use in our later analyses.
  trajectories.seq <- seq(from=1, to=n.iter.predict.population, length.out=n.trajs)
  if (make.new) {
    countries.codes <- pred.pop$countries$code
    preds.countries <- list()
    index <- 0
    for (traj in trajectories.seq) {
      index <- index + 1
      preds.countries[[index]] <- data.frame(matrix(ncol=(2+(year.end-year.start)/1), #(2+(year.end-year.start)/5)
                                                    nrow=length(isos)))
      # Picking arbitrary country (the 356 for India) for pulling variable names
      year_range <- names(bayesPop:::get.pop.trajectories(pred.pop, 356)$trajectories[,1])
      year_range <- year_range[-c(1, 2)] # delete year 2021 and 2022
      names(preds.countries[[index]]) <- c("Isocode", paste0("Pop",year_range))
    }
    row.ind <- 0
    for (iso in isos) {
      row.ind <- row.ind + 1
      country.code <- countrycode(iso, "iso3c", "un")
      if (iso == "TWN") {
        country.code <- 158
      }
      if (iso == 'DMA') {
        country.code <- 214
      }
      country.trajs <- bayesPop:::get.pop.trajectories(pred.pop, country.code)$trajectories[-c(1, 2), ] # delete year 2021 and 2022
      index <- 0
      for (traj in trajectories.seq) {
        index <- index + 1
        preds.countries[[index]][row.ind, ] <- c(iso, country.trajs[, traj])
      }
    }
    index <- 0
    for (traj in trajectories.seq) {
      index <- index + 1
      # Convert back to numeric. Silly R!
      preds.countries[[index]][,-1] <- sapply(preds.countries[[index]][,-1], as.numeric)
    }

    if (!dir.exists(paste0(sims.location, 'Sims_Pop')))
      dir.create(paste0(sims.location, 'Sims_Pop'))
    if (!dir.exists(paste0(sims.location, 'Sims_Pop/predictions')))
      dir.create(paste0(sims.location, 'Sims_Pop/predictions'))
    save(preds.countries,
         file=paste0(sims.location, "Sims_Pop/predictions/converted_pop_",
                     year.start, "_", year.end, ".RDA"))
  } else {
    load(paste0(sims.location, "Sims_Pop/predictions/converted_pop_",
                year.start, "_", year.end, ".RDA"))
  }

  preds.countries
}

convert.projections <- function(tech.projections.log,
                                gdp.projections.log, frontier.projections.log,
                                tech.USA.projections.log) {
  # convert.projections converts intensity (tech) and GDP projections
  # out of the log scale.
  tech.projections <- cbind(tech.projections.log[,1], exp(tech.projections.log[,-1]))
  names(tech.projections) <- gsub("Log", "", names(tech.projections))
  names(tech.projections)[1] <- "Isocode"
  gdp.projections <- cbind(gdp.projections.log[,1], exp(gdp.projections.log[,-1]))
  names(gdp.projections) <- gsub("Log", "", names(gdp.projections))
  names(gdp.projections)[1] <- "Isocode"
  frontier.projections <- cbind(frontier.projections.log[,1],
                                exp(frontier.projections.log[,-1]))
  names(frontier.projections) <- gsub("Log", "", names(frontier.projections))
  names(frontier.projections)[1] <- "Isocode"
  tech.USA.projections <- cbind(tech.USA.projections.log[,1],
                                exp(tech.USA.projections.log[,-1]))
  names(tech.USA.projections) <- gsub("Log", "", names(tech.USA.projections))
  names(tech.USA.projections)[1] <- "Isocode"

  list(GDPFrontier=frontier.projections, GDPData=gdp.projections,
       TechData=tech.projections, TechUSA=tech.USA.projections)
}


estimate.co2.projections <- function(tech.projections, gdp.projections,
                                     frontier.projections, tech.USA.projections,
                                     pop.projections,
                                     year.range) {
  # estimate.co2.projections combines GDP, CO2 intensity (tech), and
  # population projections to get CO2 emissions.
  # The United States must be handled separately due to its unique
  # role in the GDP model.
  stopifnot(year.range[1] %% 1 == 0 & year.range[2] %% 1 == 0)# %% 5

  # Deal with the USA specially
  pop.projections.USA <- pop.projections[pop.projections$Isocode == "USA", ]
  stopifnot(dim(pop.projections.USA)[1] == 1)

  isos.tmp <- tech.projections$Isocode
  # Population projections match up by country with tech and GDP projections
  pop.projections.sorted <- pop.projections[-1, ] # Remove the USA row
  # print(pop.projections.sorted$Isocode)
  # print(isos.tmp)
  # print(pop.projections.sorted$Isocode == isos.tmp)
  stopifnot(all(pop.projections.sorted$Isocode == isos.tmp))

  stopifnot(all(gdp.projections$Isocode == isos.tmp))

  year.vals <- seq(year.range[1], year.range[2], by=1) # by=5
  emissions.projections <- data.frame(matrix(nrow=length(isos.tmp), ncol=(length(year.vals)+1)))
  names(emissions.projections)[1] <- "Isocode"
  emissions.projections$Isocode <- isos.tmp
  names(emissions.projections)[-1] <- paste0("CO2", year.vals)
  # Add in the US specially
  emissions.projections <- rbind(c("USA", rep(NA, length(year.vals))),
                                 emissions.projections)
  n.countries <- dim(emissions.projections)[1]
  for (year in year.vals) {
    # dividing by 10 comes from intensity which has $10,000 in the denominator,
    # population which is divided by a factor of 1,000, and CO2 is in megatonnes
    var.name <- paste0("CO2", year)
    emissions.projections[-1,var.name] <- gdp.projections[, paste0("GDP", year)] *
      tech.projections[, paste0("Tech", year)] *
      pop.projections.sorted[, paste0("Pop", year)] / 10
    emissions.projections[1, var.name] <- frontier.projections[, paste0("GDP", year)] *
      tech.USA.projections[, paste0("Tech", year)] *
      pop.projections.USA[, paste0("Pop", year)] / 10
  }

  # Fuck R's conversion away from numeric type
  emissions.projections[,-1] <- sapply(emissions.projections[,-1], as.numeric)
  emissions.projections
}
#estimate.annual.co2.projections is now same as estimate.co2.projections
estimate.annual.co2.projections <- function(tech.projections, gdp.projections,
                                            frontier.projections, tech.USA.projections,
                                            pop.projections,
                                            year.range) {
  # estimate.co2.projections combines GDP, CO2 intensity (tech), and
  # population projections to get CO2 emissions.
  # The United States must be handled separately due to its unique
  # role in the GDP model.
  stopifnot(year.range[1] %% 1 == 0 & year.range[2] %% 1 == 0) #%% 5

  # Deal with the USA specially
  pop.projections.USA <- pop.projections[pop.projections$Isocode == "USA", ]
  stopifnot(dim(pop.projections.USA)[1] == 1)

  isos.tmp <- tech.projections$Isocode
  # Population projections match up by country with tech and GDP projections
  pop.projections.sorted <- pop.projections[-1, ] # Remove the USA row
  # print(pop.projections.sorted$Isocode)
  # print(isos.tmp)
  # print(pop.projections.sorted$Isocode == isos.tmp)
  stopifnot(all(pop.projections.sorted$Isocode == isos.tmp))

  stopifnot(all(gdp.projections$Isocode == isos.tmp))

  year.vals <- year.range[1]:year.range[2]
  emissions.projections <- data.frame(matrix(nrow=length(isos.tmp), ncol=(length(year.vals)+1)))
  names(emissions.projections)[1] <- "Isocode"
  emissions.projections$Isocode <- isos.tmp
  names(emissions.projections)[-1] <- paste0("CO2", year.vals)
  # Add in the US specially
  emissions.projections <- rbind(c("USA", rep(NA, length(year.vals))),
                                 emissions.projections)
  n.countries <- dim(emissions.projections)[1]
  for (year in year.vals) {
    # dividing by 10 comes from intensity which has $10,000 in the denominator,
    # population which is divided by a factor of 1,000, and CO2 is in megatonnes
    var.name <- paste0("CO2", year)
    if (year %% 1 == 0) #%% 5
    {
      emissions.projections[-1,var.name] <- gdp.projections[, paste0("GDP", year)] *
        tech.projections[, paste0("Tech", year)] *
        pop.projections.sorted[, paste0("Pop", year)] / 10
      emissions.projections[1, var.name] <- frontier.projections[, paste0("GDP", year)] *
        tech.USA.projections[, paste0("Tech", year)] *
        pop.projections.USA[, paste0("Pop", year)] / 10
    }
    else
    {
      year.left <- year - year %% 5
      year.right <- year.left + 5
      pop.projections.sorted.temp <- pop.projections.sorted[, paste0("Pop", year.right)] * ((year %% 5) / 5) +
        pop.projections.sorted[, paste0("Pop", year.left)] * (1 - (year %% 5) / 5)
      pop.projections.USA.temp <- pop.projections.USA[, paste0("Pop", year.right)] * ((year %% 5) / 5) +
        pop.projections.USA[, paste0("Pop", year.left)] * (1 - (year %% 5) / 5)
      emissions.projections[-1,var.name] <- gdp.projections[, paste0("GDP", year)] *
        tech.projections[, paste0("Tech", year)] *
        pop.projections.sorted.temp / 10
      emissions.projections[1, var.name] <- frontier.projections[, paste0("GDP", year)] *
        tech.USA.projections[, paste0("Tech", year)] *
        pop.projections.USA.temp / 10
    }
  }

  # Fuck R's conversion away from numeric type
  emissions.projections[,-1] <- sapply(emissions.projections[,-1], as.numeric)
  emissions.projections
}

# ignore this function, no use
get.ipat.pastdata <- function(data.medium, year.sequence, names.countries,
                              na.remove=F) {
  # get.ipat.pastdata gets IPAT variables worldwide for the past data.
  # If na.remove is FALSE, there will be missing data for GDP prior to
  # 1990 (before the end of the Soviet Union) and in 2010 (some GDP data
  # missing), and Tech will also be missing for these countries.
  # If na.remove is TRUE, there will be full observations for GDP and Tech,
  # though this will be misleading, as GDP numbers will be too low for
  # some years, which can also bias intensity upwards.
  # The results of get.ipat.pastdata are used for plotting,
  # and not model fitting.
  library(doBy)
  # browser()

  data.tmp1 <- subset(data.medium, Year %% 1 == 0) #Year %% 5
  data.tmp <- subset(data.tmp1, Isocode %in% names.countries)
  data.tmp$Pop <- data.tmp$PopTotal * 1000
  pop.data <- summaryBy(Pop ~ Year, data=data.tmp, FUN=sum)
  data.tmp$GDPTotal <- data.tmp$GDP * data.tmp$Pop
  gdp.data <- summaryBy(GDPTotal ~ Year, data=data.tmp, FUN=sum, na.rm=na.remove)
  data.tmp$CO2Total <- data.tmp$CO2 * data.tmp$Pop
  co2.data <- summaryBy(CO2Total ~ Year, data=data.tmp, FUN=sum, na.rm=na.remove)
  gdp.percapita <- gdp.data$GDPTotal.sum / pop.data$Pop.sum
  co2.percapita <- co2.data$CO2Total.sum / pop.data$Pop.sum
  tech.data <- co2.percapita / gdp.percapita * 10^4

  data.frame(Year=year.sequence, Pop=pop.data$Pop.sum, GDP=gdp.data$GDPTotal.sum,
             GDPpercapita=gdp.percapita,
             CO2=co2.data$CO2Total.sum, CO2percapita=co2.percapita,
             Tech=tech.data)
}



get.ipat.components.total <- function(tech.projections, gdp.projections,
                                      frontier.projections, tech.USA.projections,
                                      pop.projections, co2.projections,
                                      year.sequence) {
  # get.ipat.components.total gets IPAT components for the world for our
  # projections.

  # US is special
  pop.USA.projections <- pop.projections[1,]
  stopifnot(pop.USA.projections$Isocode == "USA")
  pop.projections.mod <- pop.projections[-1,]

  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  data.world <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.world) <- names.tmp
  data.world[, "Year"] <- year.sequence
  index <- 0
  for (year in year.sequence) {
    index <- index + 1
    var.name.co2 <- paste0("CO2", year)
    var.name.gdp <- paste0("GDP", year)

    # Sum over population
    if (year %% 1 == 0) ## year %% 5
    {
      var.name.pop <- paste0("Pop", year)
      data.world[index, "Pop"] <- sum(pop.projections[, var.name.pop]) * 1000
      # Sum over GDP
      stopifnot(all(gdp.projections$Isocode == pop.projections.mod$Isocode))
      gdp.sum.tmp1 <- sum(gdp.projections[, var.name.gdp] *
                            pop.projections.mod[, var.name.pop]) * 1000
      gdp.sum.tmp2 <- frontier.projections[, var.name.gdp] *
        pop.USA.projections[, var.name.pop] * 1000
      data.world[index, "GDP"] <- gdp.sum.tmp1 + gdp.sum.tmp2
    }
    else
    {
      year.left <- year - year %% 5
      year.right <- year.left + 5
      var.name.pop.left <- paste0("Pop", year.left)
      var.name.pop.right <- paste0("Pop", year.right)
      data.world[index, "Pop"] <- sum(pop.projections[, var.name.pop.left]) * 1000 * (1 - (year %% 5) / 5) +
        sum(pop.projections[, var.name.pop.right]) * 1000 * ((year %% 5) / 5)
      stopifnot(all(gdp.projections$Isocode == pop.projections.mod$Isocode))
      gdp.sum.tmp1 <- sum(gdp.projections[, var.name.gdp] *
                            (pop.projections.mod[, var.name.pop.left] * (1 - (year %% 5) / 5) +
                               pop.projections.mod[, var.name.pop.right] * ((year %% 5) / 5))) * 1000
      gdp.sum.tmp2 <- frontier.projections[, var.name.gdp] *
        (pop.USA.projections[, var.name.pop.left] * (1 - (year %% 5) / 5) +
           pop.USA.projections[, var.name.pop.right] * ((year %% 5) / 5)) * 1000
      data.world[index, "GDP"] <- gdp.sum.tmp1 + gdp.sum.tmp2
    }

    # Sum over emisssions
    data.world[index, "CO2"] <- sum(co2.projections[, var.name.co2])
  }
  # Get GDP per capita
  data.world$GDPpercapita <- data.world$GDP / data.world$Pop
  # Get carbon per capita. Note: this is what is CO2 in data.medium.
  data.world$CO2percapita <- data.world$CO2 / data.world$Pop
  # Get tech
  data.world$Tech <- data.world$CO2percapita / data.world$GDPpercapita * 10^4

  data.world
}



get.ipat.components.bycountry <- function(tech.projections, gdp.projections,
                                          frontier.projections, tech.USA.projections,
                                          pop.projections, co2.projections,
                                          year.sequence) {
  # get.ipat.components.bycountry reformats projections of IPAT components
  # at a country level.

  # US is special
  pop.USA.projections <- pop.projections[1,]
  stopifnot(pop.USA.projections$Isocode == "USA")
  pop.projections.mod <- pop.projections[-1,]
  stopifnot(all(pop.projections.mod$Isocode == gdp.projections$Isocode))

  #  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "FrontierGap", "CO2", "CO2percapita", "Tech")
  data.ipat.bycountry <- list()
  # US first:
  data.ipat.bycountry[["USA"]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.ipat.bycountry[["USA"]]) <- names.tmp
  data.ipat.bycountry[["USA"]]$Year <- year.sequence
  rowindex <- 0
  co2.rowindex.USA <- which(co2.projections$Isocode == "USA")
  tech.rowindex.USA <- which(tech.projections$Isocode == "USA")
  for (year in year.sequence) {
    rowindex <- rowindex + 1
    var.name.pop <- paste0("Pop", year)
    var.name.gdp <- paste0("GDP", year)
    var.name.co2 <- paste0("CO2", year)
    var.name.tech <- paste0("Tech", year)
    data.ipat.bycountry[["USA"]][rowindex, "Pop"] <- pop.USA.projections[var.name.pop]
    data.ipat.bycountry[["USA"]][rowindex, "GDPpercapita"] <- frontier.projections[var.name.gdp]
    data.ipat.bycountry[["USA"]][rowindex, "CO2"] <- co2.projections[co2.rowindex.USA, var.name.co2]
    data.ipat.bycountry[["USA"]][rowindex, "Tech"] <- tech.USA.projections[var.name.tech]
  }
  data.ipat.bycountry[["USA"]][, "FrontierGap"] <- rep(0, length(year.sequence))
  data.ipat.bycountry[["USA"]]$CO2percapita <- data.ipat.bycountry[["USA"]]$CO2 /
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)
  data.ipat.bycountry[["USA"]]$GDP <- data.ipat.bycountry[["USA"]]$GDPpercapita *
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)

  for (iso in gdp.projections$Isocode) {
    data.ipat.bycountry[[iso]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
    names(data.ipat.bycountry[[iso]]) <- names.tmp
    data.ipat.bycountry[[iso]]$Year <- year.sequence
    pop.rowindex <- which(pop.projections$Isocode == iso)
    gdp.rowindex <- which(gdp.projections$Isocode == iso)
    co2.rowindex <- which(co2.projections$Isocode == iso)
    tech.rowindex <- which(tech.projections$Isocode == iso)
    rowindex <- 0
    for (year in year.sequence) {
      rowindex <- rowindex + 1
      var.name.pop <- paste0("Pop", year)
      var.name.gdp <- paste0("GDP", year)
      var.name.co2 <- paste0("CO2", year)
      var.name.tech <- paste0("Tech", year)
      data.ipat.bycountry[[iso]][rowindex, "Pop"] <- pop.projections[pop.rowindex, var.name.pop]
      # Confusingly, for hte projections GDP is GDP per capita while CO2 is total co2 emissions
      data.ipat.bycountry[[iso]][rowindex, "GDPpercapita"] <- gdp.projections[gdp.rowindex, var.name.gdp]
      data.ipat.bycountry[[iso]][rowindex, "CO2"] <- co2.projections[co2.rowindex, var.name.co2]
      data.ipat.bycountry[[iso]][rowindex, "Tech"] <- tech.projections[tech.rowindex, var.name.tech]
    }
    data.ipat.bycountry[[iso]]$FrontierGap <- log(data.ipat.bycountry[["USA"]]$GDPpercapita /
                                                    data.ipat.bycountry[[iso]]$GDPpercapita)
    data.ipat.bycountry[[iso]]$CO2percapita <- data.ipat.bycountry[[iso]]$CO2 /
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
    data.ipat.bycountry[[iso]]$GDP <- data.ipat.bycountry[[iso]]$GDPpercapita *
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
  }

  data.ipat.bycountry
}
#get.ipat.annual.components.bycountry is now the same as get.ipat.components.bycountry
get.ipat.annual.components.bycountry <- function(tech.projections, gdp.projections,
                                                 frontier.projections, tech.USA.projections,
                                                 pop.projections, co2.projections,
                                                 year.sequence) {
  # get.ipat.components.bycountry reformats projections of IPAT components
  # at a country level.

  # US is special
  pop.USA.projections <- pop.projections[1,]
  stopifnot(pop.USA.projections$Isocode == "USA")
  pop.projections.mod <- pop.projections[-1,]
  stopifnot(all(pop.projections.mod$Isocode == gdp.projections$Isocode))
  # browser()
  #  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "FrontierGap", "CO2", "CO2percapita", "Tech")
  data.ipat.bycountry <- list()
  # US first:
  data.ipat.bycountry[["USA"]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.ipat.bycountry[["USA"]]) <- names.tmp
  data.ipat.bycountry[["USA"]]$Year <- year.sequence
  rowindex <- 0
  co2.rowindex.USA <- which(co2.projections$Isocode == "USA")
  tech.rowindex.USA <- which(tech.projections$Isocode == "USA")
  for (year in year.sequence) {
    rowindex <- rowindex + 1
    var.name.gdp <- paste0("GDP", year)
    var.name.co2 <- paste0("CO2", year)
    var.name.tech <- paste0("Tech", year)
    data.ipat.bycountry[["USA"]][rowindex, "GDPpercapita"] <- frontier.projections[var.name.gdp]
    data.ipat.bycountry[["USA"]][rowindex, "CO2"] <- co2.projections[co2.rowindex.USA, var.name.co2]
    data.ipat.bycountry[["USA"]][rowindex, "Tech"] <- tech.USA.projections[var.name.tech]
    if (year %% 1 == 0)
    {
      var.name.pop <- paste0("Pop", year)
      data.ipat.bycountry[["USA"]][rowindex, "Pop"] <- pop.USA.projections[var.name.pop]

    }
    else
    {
      year.left <- year - year %% 5
      year.right <- year.left + 5
      var.name.pop.left <- paste0("Pop", year.left)
      var.name.pop.right <- paste0("Pop", year.right)
      data.ipat.bycountry[["USA"]][rowindex, "Pop"] <-
        pop.USA.projections[var.name.pop.right] * ((year %% 5) / 5) +
        pop.USA.projections[var.name.pop.left] * (1 - (year %% 5) / 5)

    }
  }
  data.ipat.bycountry[["USA"]][, "FrontierGap"] <- rep(0, length(year.sequence))
  data.ipat.bycountry[["USA"]]$CO2percapita <- data.ipat.bycountry[["USA"]]$CO2 /
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)
  data.ipat.bycountry[["USA"]]$GDP <- data.ipat.bycountry[["USA"]]$GDPpercapita *
    (10^3 * data.ipat.bycountry[["USA"]]$Pop)

  for (iso in gdp.projections$Isocode) {
    data.ipat.bycountry[[iso]] <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
    names(data.ipat.bycountry[[iso]]) <- names.tmp
    data.ipat.bycountry[[iso]]$Year <- year.sequence
    pop.rowindex <- which(pop.projections$Isocode == iso)
    gdp.rowindex <- which(gdp.projections$Isocode == iso)
    co2.rowindex <- which(co2.projections$Isocode == iso)
    tech.rowindex <- which(tech.projections$Isocode == iso)
    rowindex <- 0
    for (year in year.sequence) {
      rowindex <- rowindex + 1
      var.name.gdp <- paste0("GDP", year)
      var.name.co2 <- paste0("CO2", year)
      var.name.tech <- paste0("Tech", year)
      # Confusingly, for hte projections GDP is GDP per capita while CO2 is total co2 emissions
      data.ipat.bycountry[[iso]][rowindex, "GDPpercapita"] <- gdp.projections[gdp.rowindex, var.name.gdp]
      data.ipat.bycountry[[iso]][rowindex, "CO2"] <- co2.projections[co2.rowindex, var.name.co2]
      data.ipat.bycountry[[iso]][rowindex, "Tech"] <- tech.projections[tech.rowindex, var.name.tech]
      if (year %% 1 == 0)
      {
        var.name.pop <- paste0("Pop", year)
        data.ipat.bycountry[[iso]][rowindex, "Pop"] <- pop.projections[pop.rowindex, var.name.pop]
      }
      else
      {
        year.left <- year - year %% 5
        year.right <- year.left + 5
        var.name.pop.left <- paste0("Pop", year.left)
        var.name.pop.right <- paste0("Pop", year.right)
        data.ipat.bycountry[[iso]][rowindex, "Pop"] <-
          pop.projections[pop.rowindex, var.name.pop.right] * ((year %% 5) / 5) +
          pop.projections[pop.rowindex, var.name.pop.left] * (1 - (year %% 5) / 5)
      }
    }
    data.ipat.bycountry[[iso]]$FrontierGap <- log(data.ipat.bycountry[["USA"]]$GDPpercapita /
                                                    data.ipat.bycountry[[iso]]$GDPpercapita)
    data.ipat.bycountry[[iso]]$CO2percapita <- data.ipat.bycountry[[iso]]$CO2 /
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
    data.ipat.bycountry[[iso]]$GDP <- data.ipat.bycountry[[iso]]$GDPpercapita *
      (10^3 * data.ipat.bycountry[[iso]]$Pop)
  }

  data.ipat.bycountry
}



get.ipat.components.ssa <- function(tech.projections, gdp.projections,
                                    pop.projections, co2.projections,
                                    year.sequence) {
  # get.ipat.components.total gets IPAT components for Sub-Saharan Africa
  # for our projections. Similar to get.ipat.components.total

  # Subset to Sub-Saharan Africa (SSA)
  tech.projections <- subset(tech.projections, Isocode %in% ssa.isos)
  gdp.projections <- subset(gdp.projections, Isocode %in% ssa.isos)
  pop.projections <- subset(pop.projections, Isocode %in% ssa.isos)
  co2.projections <- subset(co2.projections, Isocode %in% ssa.isos)

  names.tmp <- c("Year", "Pop", "GDP", "GDPpercapita", "CO2", "CO2percapita", "Tech")
  data.ssa <- data.frame(matrix(nrow=length(year.sequence), ncol=length(names.tmp)))
  names(data.ssa) <- names.tmp
  data.ssa[, "Year"] <- year.sequence
  index <- 0
  for (year in year.sequence) {
    index <- index + 1
    var.name.co2 <- paste0("CO2", year)
    var.name.gdp <- paste0("GDP", year)
    var.name.pop <- paste0("Pop", year)
    # Sum over population
    data.ssa[index, "Pop"] <- sum(pop.projections[, var.name.pop]) * 1000
    # Sum over GDP
    stopifnot(all(gdp.projections$Isocode == pop.projections$Isocode))
    data.ssa[index, "GDP"] <- sum(gdp.projections[, var.name.gdp] *
                                    pop.projections[, var.name.pop]) * 1000
    # Sum over emisssions
    data.ssa[index, "CO2"] <- sum(co2.projections[, var.name.co2])
  }
  # Get GDP per capita
  data.ssa$GDPpercapita <- data.ssa$GDP / data.ssa$Pop
  # Get carbon per capita. Note: this is what is CO2 in data.medium.
  data.ssa$CO2percapita <- data.ssa$CO2 / data.ssa$Pop
  # Get tech
  data.ssa$Tech <- data.ssa$CO2percapita / data.ssa$GDPpercapita * 10^4

  data.ssa
}




get.trajs.quants <- function(co2.projections, year.seq, quantiles=c(0.1,0.5,0.9)) {
  # get.trajs.quants gets quantiles for worldwide CO2 emissions
  trajs.quants <- data.frame(Quantile=quantiles)
  for (year in year.seq) {
    var.name <- paste0("CO2", year)
    trajs.quants <- cbind(trajs.quants,
                          quantile(sapply(1:length(co2.projections), function(i)
                            #  median(sapply(1:length(co2.projections), function(i)
                            sum(co2.projections[[i]][, var.name])), probs=quantiles) / 10^9)
  }
  names(trajs.quants)[-1] <- paste0("CO2", year.seq)

  trajs.quants
}



get.trajs.quants.bycountry <- function(co2.projections, year.seq,
                                       quantiles=c(0.1,0.5,0.9)) {
  # get.trajs.quants.bycountry gets quantiles for CO2 emissions by country.
  # Compare with get.trajs.quants
  isos <- co2.projections[[1]]$Isocode
  trajs.quants <- list()
  for (year in year.seq) {
    trajs.quants[[paste0("CO2", year)]] <- data.frame(Isocode=isos)
    for (q in quantiles) {
      trajs.quants[[paste0("CO2", year)]][, paste0("CO2_", q)] <- rep(NA, length(isos))
    }
  }
  for (indx in 1:length(isos)) {
    iso <- isos[indx]
    for (year in year.seq) {
      var.name <- paste0("CO2", year)
      trajs.quants[[var.name]][indx, -1] <-  t(quantile(sapply(1:length(co2.projections), function(i)
        co2.projections[[i]][indx, var.name]), probs=quantiles, na.rm=T) / 10^9)
    }
  }

  trajs.quants
}


get.ipat.medians <- function(trajs.worldwide, quantiles=c(0.1, 0.5, 0.9)) {
  # get.ipat.medians gets worldwide quantiles (not just medians) of IPAT variables,
  # based on projections trajs.worldwide
  trajs.quants <- list()
  for (q in quantiles) {
    trajs.quants[[as.character(q)]] <- trajs.worldwide[[1]]
    trajs.quants[[as.character(q)]][, -1] <- NA
  }

  row.index <- 0
  for (year in trajs.quants[[1]]$Year) {
    row.index <- row.index + 1
    for (var.name in names(trajs.worldwide[[1]])[-1]) {
      if (trajs.quants[[1]]$Year[1] >= 1990) {
        quants.tmp <- quantile(sapply(1:length(trajs.worldwide), function(i)
          trajs.worldwide[[i]][row.index, var.name]), probs=quantiles)
      } else {
        quants.tmp <- quantile(sapply(1:length(trajs.worldwide), function(i)
          trajs.worldwide[[i]][row.index, var.name]), probs=quantiles, na.rm=T)
      }
      for (q_ind in 1:length(quantiles)) {
        q <- quantiles[q_ind]
        trajs.quants[[as.character(q)]][row.index, var.name] <- quants.tmp[q_ind]
      }
    }
  }

  trajs.quants
}


get.ipat.quantiles.bycountry <- function(ipat.components.bycountry, quantiles=c(0.1, 0.5, 0.9),
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


# =========================================



# Modelling Carbon Intensity:
make.bugs.data.ar1 <- function(data.medium, max.vals, remove.before.peak=TRUE,
                               isos.remove=c()) {
  # make.bugs.data.ar1 creates the variables that will be used by JAGS in fitting
  # the GDP and intensity models given in corr_model_ar1trend_const.bug

  # Only removing indeces for the Intensity data, not for GDP data
  indeces.remove <- data.medium$Isocode %in% isos.remove
  indeces.remove.gdp <- data.medium$Isocode == "USA"
  data.tmp <- data.medium[!indeces.remove, c("Year","Isocode","Tech")]
  data.tmp$Tech <- log(data.tmp$Tech)
  data.wide.tech <- dcast(data.tmp[!duplicated(data.tmp[,c("Year","Isocode")]),],
                          Year ~ Isocode, value.var="Tech")
  data.tmp.gdp <- data.medium[!indeces.remove.gdp, c("Year","Isocode","GDP")]
  data.tmp.gdp$GDP <- log(data.tmp.gdp$GDP)
  data.wide.gdp <- dcast(data.tmp.gdp[!duplicated(data.tmp.gdp[,c("Year","Isocode")]),],
                         Year ~ Isocode, value.var="GDP")

  #browser()
  tau <- data.wide.tech
  n.years <- dim(tau)[1]
  n.countries <- dim(tau)[2]
  n.countries.gdp <- dim(data.wide.gdp)[2]

  # Get isos for GDP and intensity. 1st element is "Year", which we remove
  isos.gdp <- names(data.wide.gdp)[-1]
  isos.tau <- names(tau)[-1]

  # Make the frontier:
  frontier.gdp <- log(data.medium$GDP[data.medium$Isocode == "USA"])
  frontier.gap <- frontier.gdp - data.wide.gdp # first column is year and is meaningless

  frontier.gap.old.tmp <- frontier.gap[1:(n.years-1), 2:n.countries.gdp]
  frontier.gap.new.tmp <- frontier.gap[2:n.years, 2:n.countries.gdp]
  tau.delta.tmp.gdp <- frontier.gap.new.tmp - frontier.gap.old.tmp
  tau.old.matrix <- tau[1:(n.years-1), 2:n.countries]
  tau.new.matrix <- tau[2:n.years, 2:n.countries]
  n.years <- n.years - 1
  n.countries <- n.countries - 1
  n.countries.gdp <- n.countries.gdp - 1

  # Adding in the USA separately for intensity:
  us.tech <- log(data.medium$Tech[data.medium$Isocode == "USA"])
  tau.new.USA <- us.tech[-1]
  tau.old.USA <- us.tech[-length(us.tech)]

  tau.delta.tmp <- tau.new.matrix - tau.old.matrix
  # This block removes years before the peak
  if (remove.before.peak) {
    print("Removing data before the peak...")
    c <- 0
    # Assume they peaked at the start by default for now.
    # Using max.vals from find_max.R
    countries.peak <- rep(0, n.countries)
    for (iso.name in names(tau.delta.tmp)) {
      c <- c + 1
      if (sum(max.vals$Iso == iso.name) == 1) {
        # we have a max value for that country
        countries.peak[c] <- max.vals$Year[which(max.vals$Iso == iso.name)]
      }
    }
    smallest.year <- min(data.medium$Year)
    for (c in 1:n.countries) {
      for (t in 1:n.years) {
        if (t <= (countries.peak[c] - smallest.year)) {
          tau.delta.tmp[t,c] <- NA
        }
        # if (t == (countries.peak[c] - smallest.year) && !is.na(tau.delta.tmp[t,c]))
        # {
        #   if (tau.delta.tmp[t,c] < 0) {tau.delta.tmp[t,c] <- NA}
        # }
      }
    }
  } else {
    print("Not removing data before the peak.")
  }



  n.GoodObs <- sum(!is.na(tau.delta.tmp.gdp))
  n.GoodObs.c <- rep(0, n.countries.gdp)
  frontier.gap.old <- frontier.gap.new <- rep(NA, n.GoodObs)
  tau.new <- tau.old <- years.intensity <- year.inds <- country.inds <- country.gdp.inds <- rep(NA, n.GoodObs)
  tau.new.funny.index <- c()
  intensity.Obs <- c()
  z <- 0
  z.intensity <- 0
  for (c in 1:n.countries.gdp){
    for (t in 1:n.years){
      if (!is.na(tau.delta.tmp.gdp[t,c])) {
        z <- z + 1
        n.GoodObs.c[c] <- n.GoodObs.c[c] + 1
        if (isos.gdp[c] %in% isos.tau) {
          c.tau <- which(isos.tau == isos.gdp[c])
          stopifnot(length(c.tau) == 1)
          if ((!is.na(tau.delta.tmp[t,c.tau])) & (!is.na(tau.new.matrix[t, c.tau]))
              & (!is.na(tau.old.matrix[t, c.tau]))) {
            z.intensity <- z.intensity + 1
            tau.new[z] <- tau.new.matrix[t, c.tau]
            tau.new.funny.index <- c(tau.new.funny.index, tau.new.matrix[t, c.tau])
            tau.old[z] <- tau.old.matrix[t, c.tau]
            years.intensity[z] <- t + 1960 - 1
            intensity.Obs <- c(intensity.Obs, z)
          }
          country.inds[z] <- c.tau
        }
        frontier.gap.old[z] <- frontier.gap.old.tmp[t, c]
        frontier.gap.new[z] <- frontier.gap.new.tmp[t, c]

        year.inds[z] <- t # not used now
        country.gdp.inds[z] <- c
      }
    }
    # cat(c, ' : ', n.GoodObs.c[c], '\n')
  }
  n.GoodObsIntensity <- z.intensity
  # browser()
  list(
    tau.new=tau.new, tau.old=tau.old,
    tau.new.funny.index=tau.new.funny.index,
    years.intensity=years.intensity,
    frontier.gap.old=frontier.gap.old,
    frontier.gap.new=frontier.gap.new, frontier.gdp=frontier.gdp,
    year.inds=year.inds,
    country.gdp.inds=country.gdp.inds,
    country.inds=country.inds,
    n.countries=n.countries, n.countries.gdp=n.countries.gdp,
    n.years=n.years, n.GoodObs=n.GoodObs,
    names.countries=isos.tau,
    names.countries.gdp=isos.gdp,
    n.GoodObs.c=n.GoodObs.c,
    tau.new.USA=tau.new.USA,
    tau.old.USA=tau.old.USA,
    n.GoodObsIntensity=n.GoodObsIntensity,
    intensity.Obs=intensity.Obs)
}


fit.corr.model.ar1 <- function(data.medium, isos.remove=c("USA"), max.vals,
                               n.iterations=1000, n.adapt=300, n.chains=5, thin=20,
                               model.name="corr_model_ar1trend.bug",
                               var.list=c('delta', 'sig.eps', 'mu.delta', 'sig.delta',
                                          'phi', 'sig.gap', 'gamma', 'gamma.pre1973',
                                          'sig.frontier', 'mu.phi', 'sig.phi',
                                          'mu.ln.gap', 'sig.ln.gap', 'sig.eps.mean','sig.eps.sd',
                                          'rho',
                                          'eta', 'mu.eta', 'sig.eta',
                                          'beta', 'mu.beta', 'sig.beta',
                                          'delta.USA', 'sig.eps.USA')) {
  # fit.corr.model.ar1 uses JAGS to fit the model described in
  # corr_model_ar1trend_const.bug for GDP and intensity.
  # It calls make.bugs.data.ar1() to get the variables
  # that will be used by JAGS, when calling jags.model() and coda.samples()

  library(rjags)
  library(coda)
  library(reshape2)
  print("Calling make.bugs.data")
  jags.input <- make.bugs.data.ar1(data.medium, max.vals,
                                   isos.remove=isos.remove, remove.before.peak=T)
  # browser()
  print("Successfully called make.bugs.data")
  tau.new <- jags.input$tau.new
  tau.new.funny.index <- jags.input$tau.new.funny.index
  tau.old <- jags.input$tau.old
  years.intensity <- jags.input$years.intensity
  frontier.gap.old <- jags.input$frontier.gap.old
  frontier.gap <- jags.input$frontier.gap.new
  frontier.gdp <- jags.input$frontier.gdp
  year.inds <- jags.input$year.inds

  country.gdp.inds <- jags.input$country.gdp.inds # todo
  country.inds <- jags.input$country.inds
  n.countries <- jags.input$n.countries
  n.countries.gdp <- jags.input$n.countries.gdp
  n.years <- jags.input$n.years
  n.GoodObs <- jags.input$n.GoodObs
  names.countries <- jags.input$names.countries
  names.countries.gdp <- jags.input$names.countries.gdp
  n.GoodObs.c <- jags.input$n.GoodObs.c
  tau.new.USA <- jags.input$tau.new.USA
  tau.old.USA <- jags.input$tau.old.USA
  # Treat GDP and Intensity separately
  n.GoodObsIntensity <- jags.input$n.GoodObsIntensity
  intensity.Obs <- jags.input$intensity.Obs

  frontier.change <- frontier.gdp[2:(n.years+1)] - frontier.gdp[1:n.years]
  jags.input$frontier.change <- frontier.change

  print("here")
  print(n.countries)
  print("here")
  #jags.out <- jags.model(model.name, n.chains=n.chains, n.adapt=n.adapt)
  jags.out <- jags.model(model.name,
                         data=jags.input,
                         n.chains=n.chains, n.adapt=n.adapt)
  jags.output <- coda.samples(jags.out,
                              var.list,
                              n.iter=n.iterations,
                              thin=thin)
  list(jags.input, jags.output)
}



project.corr.ar1.const <- function(data.medium, input.data, jags.corr.output, year.range,
                                   param.indeces=seq(400, 500, by=10)) {
  # project.corr.ar1.const uses the fitted model for GDP and intensity
  # as fit by fit.corr.model.ar1 and generations projections forward

  jags.input <- input.data

  tau.change <- jags.input$tau.delta
  frontier.gap.old <- jags.input$frontier.gap.old
  frontier.gap <- jags.input$frontier.gap.new
  frontier.gdp <- jags.input$frontier.gdp
  year.inds <- jags.input$year.inds
  country.gdp.inds <- jags.input$country.gdp.inds
  country.inds <- jags.input$country.inds
  n.countries <- jags.input$n.countries
  n.countries.gdp <- jags.input$n.countries.gdp
  n.years <- jags.input$n.years
  n.GoodObs <- jags.input$n.GoodObs
  names.countries <- jags.input$names.countries
  names.countries.gdp <- jags.input$names.countries.gdp
  n.GoodObs.c <- jags.input$n.GoodObs.c
  tau.change.USA <- jags.input$tau.change.USA
  # Treat GDP and Intensity separately
  n.GoodObsIntensity <- jags.input$n.GoodObsIntensity
  intensity.Obs <- jags.input$intensity.Obs

  #CHANGE
  year.first <- year.range[1]
  # year.first <- min(year.range[1], 2008)
  year.last <- year.range[2]
  print("Time to create the base datasets")
  # Create base GDP data:
  gdp.data.tmp1 <- data.medium[data.medium$Year == year.first, c("Isocode", "GDP")]
  gdp.data.tmp2 <- gdp.data.tmp1[gdp.data.tmp1$Isocode %in% names.countries.gdp,
                                 c("Isocode", "GDP")]
  gdp.data.tmp <- cbind(gdp.data.tmp2, t(rep(NA, (year.last - year.first))))
  names(gdp.data.tmp) <- c("Isocode", paste0("LogGDP", year.first:year.last))
  gdp.data.tmp[, 2] <- log(gdp.data.tmp[, 2])

  # Create base GDP frontier data:
  stopifnot(all(gdp.data.tmp$Isocode == names.countries.gdp))
  stopifnot(! "USA" %in% names.countries.gdp)
  gdp.frontier <- data.frame(Isocode="USA",
                             GDP=data.medium$GDP[data.medium$Isocode == "USA" &
                                                   data.medium$Year == year.first],
                             t(rep(NA, (year.last - year.first))))
  names(gdp.frontier) <- c("Isocode", paste0("LogGDP", year.first:year.last))
  gdp.frontier[, 2] <- log(gdp.frontier[, 2])
  # Create base Intensity data:
  tech.data.tmp1 <- data.medium[data.medium$Year == year.first, c("Isocode", "Tech")]
  # Use names.countries.gdp b/c we'll project ahead for all countries.
  tech.data.tmp2 <- tech.data.tmp1[tech.data.tmp1$Isocode %in% names.countries.gdp,
                                   c("Isocode", "Tech")]
  tech.data.tmp <- cbind(tech.data.tmp2, t(rep(NA, (year.last - year.first))))
  names(tech.data.tmp) <- c("Isocode", paste0("LogTech", year.first:year.last))
  tech.data.tmp[, 2] <- log(tech.data.tmp[, 2])
  # Create base USA Intensity data
  tech.USA.tmp <- data.frame(Isocode="USA",
                             Tech=data.medium$Tech[data.medium$Isocode == "USA" &
                                                     data.medium$Year == year.first],
                             t(rep(NA, (year.last - year.first))))
  names(tech.USA.tmp) <- c("Isocode", paste0("LogTech", year.first:year.last))
  tech.USA.tmp[, 2] <- log(tech.USA.tmp[, 2])

  # Project ahead for most countries
  print("Projecting ahead for countries")
  names.excluded <- setdiff(names.countries.gdp,
                            c("USA", names.countries))
  n.excluded <- length(names.excluded)
  n.chains <- length(jags.corr.output)
  list.proj <- list()
  list.index <- 0
  for (chain.num in 1:n.chains) {
    # browser()
    for (index in param.indeces) {
      list.index <- list.index + 1
      nrow.corr <- nrow(jags.corr.output[[chain.num]])
      while(1)
      {
        vars.tmp1 <- jags.corr.output[[chain.num]][index,]
        vars.tmp <- data.frame(names(vars.tmp1), vars.tmp1)

        eta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "mu.eta"]
        if (eta.tmp < 0)
        {
          break
        }
        else
        {
          index <- index + 1
          if (index > nrow.corr)
          {
            index <- 1
          }
        }
      }


      vars.tmp <- data.frame(names(vars.tmp1), vars.tmp1)
      eta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "mu.eta"]
      beta.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 7) == "mu.beta"]
      deltas.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 6) == "delta["]
      delta.USA <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 9) == "delta.USA"]
      phis <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 3) == "phi"]
      rho <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 3) == "rho"]
      sig.eps.vals.tmp <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 8) == "sig.eps["]
      sig.eps.USA <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 11) == "sig.eps.USA"]
      sig.gap <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 7) == "sig.gap"]
      sig.frontier <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 12) == "sig.frontier"]
      gamma <- vars.tmp$vars.tmp[substr(vars.tmp$names.vars.tmp1, 1, 5) == "gamma"]
      stopifnot(length(gamma) == 1) # make sure not using gamma.pre1973

      # Used for prepeak countries. Pull out hyperparameters, then generate country parameters.
      deltas <- rep(NA, n.countries.gdp)
      sig.eps.vals <- rep(NA, n.countries.gdp)
      stopifnot(length(deltas) == length(phis))
      # Fill out deltas with postpeak countries
      gdp.intensity.inds <- names.countries.gdp %in% names.countries
      gdp.excluded.inds <- ! names.countries.gdp %in% names.countries
      stopifnot(names.countries.gdp[gdp.intensity.inds] == names.countries)
      deltas[gdp.intensity.inds] <- deltas.tmp
      sig.eps.vals[gdp.intensity.inds] <- sig.eps.vals.tmp
      # Get hyperdistributions for prepeak countries
      var.names <- varnames(jags.corr.output[[chain.num]])
      mu.delta.tmp <- vars.tmp1[which(var.names == "mu.delta")]
      sig.delta.tmp <- vars.tmp1[which(var.names == "sig.delta")]
      sig.eps.mean.tmp <- vars.tmp1[which(var.names == "sig.eps.mean")]
      sig.eps.sd.tmp <- vars.tmp1[which(var.names == "sig.eps.sd")]
      # Now generate parameters for prepeak countries
      delta.c <- rep(NA, n.excluded)
      sigma.c <- rep(NA, n.excluded)
      for (c.ind in 1:n.excluded) {
        delta.c[c.ind] <- max(0, rnorm(1, mu.delta.tmp, sig.delta.tmp))
        sigma.c[c.ind] <- rlnorm(1, sig.eps.mean.tmp, sig.eps.sd.tmp)
      }
      deltas[gdp.excluded.inds] <- delta.c
      sig.eps.vals[gdp.excluded.inds] <- sigma.c

      # Create dataframe for GDP. First year is last year in the data
      # Create dataframe for Intensity
      gdp.frontier.copy <- gdp.frontier
      gdp.data.copy <- gdp.data.tmp
      tech.data.copy <- tech.data.tmp
      tech.USA.copy <- tech.USA.tmp
      #CHANGE
      # if (year.range[1] == 2010) {
      #   # We have missing data for some countries for 2009 and 2010, so we have
      #   # this very awkward kludge to deal with this inconsistency.
      #   # In this case we started with a base in 2008, so we now have to fill in
      #   # 2009 and 2010 data, for
      #   # gdp.frontier.copy, gdp.data.copy, tech.data.copy, tech.USA.copy
      #   for (year in 2009:2010) {
      #     print("Doing special 2009, 2010 stuff for USA")
      #     gdp.frontier.copy[, paste0("LogGDP", year)] <- log(subset(data.medium, Isocode == "USA" &
      #                                                                 Year == year)$GDP)
      #     tech.USA.copy[, paste0("LogTech", year)] <- log(subset(data.medium, Isocode == "USA" &
      #                                                              Year == year)$Tech)
      #   }
      #   for (year in 2009:2010) {
      #     print("Doing special 2009, 2010 stuff for non-USA countries")
      #     # Project non-USA countries as if starting from 2008. We'll replace
      #     # simulated values with true values we know after this for loop
      #     year.old <- year - 1
      #     gdp.name.now <- paste0("LogGDP", year)
      #     gdp.name.old <- paste0("LogGDP", year.old)
      #     tech.name.now <- paste0("LogTech", year)
      #     tech.name.old <- paste0("LogTech", year.old)
      #     counter <- 0
      #     while (counter == 0 ||
      #            (max(tech.data.copy[, tech.name.now], na.rm=T) > log(50))) {
      #       counter <- 1
      #       # Project other countries GDP
      #       gdp.epsilon <- rnorm(n.countries.gdp, 0, sig.gap)
      #       gdp.data.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.now] - phis *
      #         (gdp.frontier.copy[, gdp.name.old] - gdp.data.copy[, gdp.name.old]) -
      #         gdp.epsilon
      #       # Project Intensity
      #       tech.epsilon <- rnorm(n.countries.gdp,
      #                             rho*(sig.eps.vals/sig.gap)*
      #                               gdp.epsilon,
      #                             sqrt(1-rho)*sig.eps.vals)
      #       tech.data.copy[, tech.name.now] <- eta.tmp * (year - 1990) +
      #         beta.tmp * tech.data.copy[, tech.name.old] -
      #         deltas + tech.epsilon
      #     }
      #   }
      #   for (year in 2009:2010) {
      #     # Now replace simulated values with true values we know for 2009 and 2010
      #     gdp.name.now <- paste0("LogGDP", year)
      #     tech.name.now <- paste0("LogTech", year)
      #
      #     # GDP
      #     log.gdp.data.tmp <- log(subset(data.medium, Year == year &
      #                                      Isocode %in% names.countries.gdp)[, "GDP"])
      #     which.have.gdp.data <- ! is.na(log.gdp.data.tmp)
      #     # Replace simulated values with observed values
      #     gdp.data.copy[which.have.gdp.data, gdp.name.now] <- log.gdp.data.tmp[which.have.gdp.data]
      #
      #     # Intensity
      #     log.tech.data.tmp <- log(subset(data.medium, Year == year &
      #                                       Isocode %in% names.countries.gdp)[, "Tech"])
      #     which.have.tech.data <- ! is.na(log.tech.data.tmp)
      #     # Replace simulated values with observed values
      #     tech.data.copy[which.have.tech.data, tech.name.now] <- log.tech.data.tmp[which.have.tech.data]
      #   }
      # }
      for (year in (year.range[1]+1):year.range[2]) {
        # year is the year we're projecting towards
        year.old <- year - 1
        gdp.name.now <- paste0("LogGDP", year)
        gdp.name.old <- paste0("LogGDP", year.old)
        tech.name.now <- paste0("LogTech", year)
        tech.name.old <- paste0("LogTech", year.old)
        counter <- 0
        while (counter == 0 ||
               (max(tech.data.copy[, tech.name.now], na.rm=T) > log(50))) {
          #counter <- 1
          counter <- counter + 1
          if (counter > 1) {
            tech.data.copy[, tech.name.now] <- sapply(tech.data.copy[, tech.name.now],
                                                      function(x) min(x, log(50))
            )
            break
          }
          stopifnot(counter < 10)
          # Project frontier GDP
          gdp.epsilon.frontier <- rnorm(1, 0, sig.frontier)
          gdp.frontier.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.old] +
            gamma + gdp.epsilon.frontier
          # Project other countries GDP
          gdp.epsilon <- rnorm(n.countries.gdp, 0, sig.gap)
          gdp.data.copy[, gdp.name.now] <- gdp.frontier.copy[, gdp.name.now] - phis *
            (gdp.frontier.copy[, gdp.name.old] - gdp.data.copy[, gdp.name.old]) -
            gdp.epsilon
          # Project Intensity
          tech.epsilon <- rnorm(n.countries.gdp,
                                rho*(sig.eps.vals/sig.gap)*
                                  gdp.epsilon,
                                sqrt(1-rho)*sig.eps.vals)
          tech.data.copy[, tech.name.now] <- eta.tmp * (year - 1990) +
            beta.tmp * tech.data.copy[, tech.name.old] -
            deltas + tech.epsilon
        }
        # Project USA Intensity
        tech.epsilon.USA <- rnorm(1,
                                  -rho*(sig.eps.USA/sig.frontier)*gdp.epsilon.frontier,
                                  sqrt(1-rho)*sig.eps.USA)
        tech.USA.copy[, tech.name.now] <- eta.tmp * (year - 1990) +
          beta.tmp * tech.USA.copy[, tech.name.old] -
          delta.USA + tech.epsilon.USA
      }
      list.proj[[list.index]] <- list(GDPFrontier=gdp.frontier.copy,
                                      GDPData=gdp.data.copy, TechData=tech.data.copy,
                                      TechUSA=tech.USA.copy)
    }
  }
  print("returning list.proj")

  list.proj
}





# =========================================

get.cumulative.emissions <- function(co2.projection) {
  # get.cumulative.emissions is for getting cumulative emissions from 2010 to 2100

  # Warning: only works if starting from 2010
  stopifnot(names(co2.projection)[1] == "Isocode")
  emit.yearly <- colSums(co2.projection[,-1])
  cum.emissions <- emit.yearly
  cum.emissions.total <- 0
  count <- 0
  for (var.name in names(emit.yearly)) {
    if (count != 0) {
      cum.emissions.total <- cum.emissions.total + 2.5*emit.yearly[var.name]
    }
    cum.emissions[var.name] <- cum.emissions.total
    count <- count + 1
    if (count != length(names(emit.yearly))) {
      cum.emissions.total <- cum.emissions.total + 2.5*emit.yearly[var.name]
    }
  }
  cum.emissions.total <- unname(cum.emissions.total)

  return(list(cum.emissions.total=cum.emissions.total,
              cum.emissions=cum.emissions,
              emit.yearly=emit.yearly))
}


is.carbon.below.limit <- function(co2.projection,
                                  co2.cutoff=11000,
                                  year.list=seq(2023,2100,by=1)) {
  # is.carbon.below.limit checks whether a projection from 2010 to 2100
  # has cumulative CO2 emissions below the cutoff of 11000 gt CO2 (check 11000)
  #browser()
  stopifnot(names(co2.projection)[1] == "Isocode")
  name.list <- paste0("CO2", year.list)
  carbon.cum.tmp <- 0
  # It's as if we're starting from 2023, as in the paper we cite talking about
  # remaining CO2
  for (name in name.list) {
    if (name == "CO22023" || name == "CO22100") {
      carbon.cum.tmp <- carbon.cum.tmp + 1*sum(co2.projection[, name]) # change 2.5 to 1 due to one year step
    } else {
      carbon.cum.tmp <- carbon.cum.tmp + 1*sum(co2.projection[, name]) # change 5 to 1 due to one year step
    }
  }
  carbon.cum <- carbon.cum.tmp / 10^9 # convert to gigatonnes
  stopifnot(!is.na(carbon.cum))
  print(carbon.cum)

  if (carbon.cum < co2.cutoff) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


fit.project.model <- function(data.medium, year.start=2010, year.end=2100,
                              fit.ar1.const=F, n.trajectories=1000) {
  # fit.project.model fits the models and generates forward projections.
  # It uses population predictions already calculated above, calls fit.corr.model.ar1
  # to fit the model for GDP and intensity, and calls project.corr.ar1.const
  # to get projections for GDP and intensity.
  # It then has to combine the population, GDP, and intensity projections.

  library(rjags)
  library(coda)
  library(reshape2)
  data.medium.tmp <- subset(data.medium, Year <= year.start)

  # Calling find_maxes() from find_max.R
  print("======================================")
  print("Find max intensities using find_max.R")
  print("")
  maxes.countries <- find_maxes(data.medium.tmp)
  max.vals <- maxes.countries[[1]]
  rejects.late <- maxes.countries[[3]]
  rejects.insuf <- maxes.countries[[4]]

  # These functions are from estimate_model_corr.R
  print("======================================")
  print("Fit the model! Using fit.corr.model")
  print("")
  n.chains <- 5
  n.iter <- 100000
  # n.iter <- 5000
  n.burnin <- 5000
  # n.burnin <- 1000
  thin <- 20
  #thin <- 1
  stopifnot(n.trajectories %% n.chains == 0)
  param.indeces <- round(seq(1, floor(n.iter / thin), length.out=(n.trajectories / n.chains)))
  #CHANGE
  #year.start.tmp <- min(2008, year.start)
  year.start.tmp <- year.start
  if (! fit.ar1.const) {
    corr.model.out <- fit.corr.model(data.medium.tmp,
                                     isos.remove=c("USA",rejects.late, rejects.insuf), max.vals,
                                     n.iterations=n.iter, n.adapt=n.burnin, n.chains=5,
                                     thin=thin)
    names.countries <- corr.model.out[[1]]$names.countries.gdp # Now we're including all countries
    names.countries.intensity <- corr.model.out[[1]]$names.countries # For only post-peak countries

    print("======================================")
    print("Projecting intensity and GDP using project.corr")
    print("")
    list.proj <- project.corr(data.medium.tmp, input.data=corr.model.out[[1]],
                              jags.corr.output=corr.model.out[[2]], year.range=c(year.start.tmp,year.end),
                              param.indeces=param.indeces)
  } else {
    # can't get gamma.pre1973, b/c when pattern matching later we'd pull it out along
    # with gamma
    #CHANGE
    #corr.model.out <- fit.corr.model.ar1(data.medium,
    corr.model.out <- fit.corr.model.ar1(data.medium,
                                         isos.remove=c("USA", rejects.late, rejects.insuf), max.vals,
                                         n.iterations=n.iter, n.adapt=n.burnin, n.chains=5,
                                         thin=thin,
                                         model.name="corr_model_ar1trend_const.bug",
                                         var.list=c('delta', 'sig.eps', 'mu.delta', 'sig.delta',
                                                    'phi', 'sig.gap', 'gamma',
                                                    'sig.frontier', 'mu.phi', 'sig.phi',
                                                    'mu.ln.gap', 'sig.ln.gap', 'sig.eps.mean','sig.eps.sd',
                                                    'rho',
                                                    'mu.eta',
                                                    'mu.beta',
                                                    'delta.USA', 'sig.eps.USA',
                                                    'deviance'))
    names.countries <- corr.model.out[[1]]$names.countries.gdp # Now we're including all countries
    names.countries.intensity <- corr.model.out[[1]]$names.countries # For only post-peak countries

    print("======================================")
    print("Projecting intensity and GDP using project.corr")
    print("")
    list.proj <- project.corr.ar1.const(data.medium.tmp, input.data=corr.model.out[[1]],
                                        jags.corr.output=corr.model.out[[2]], year.range=c(year.start.tmp,year.end),
                                        param.indeces=param.indeces)
  }

  # From predict_pop.R
  # Getting population trajectories are slow
  print("======================================")
  print("Fitting the population model and getting population trajectories!")
  print("")
  if (year.start == 2010 & year.end == 2100) {
    # gives preds.countries.trajs
    load(paste0(data.location, "poppreds_formatted_2010_2100.rda"))
  } else if (year.start == 2000 & year.end == 2010) {
    load(paste0(data.location, "poppreds_formatted_2000_2010.rda"))
  } else {
    load(paste0(data.location, "poppreds_formatted_",
                year.start, "_", year.end, ".rda"))
    # preds.countries.tmp <- predict.population(year.present=year.start, year.end=year.end,
    #                                           make.new=FALSE)
    # #                                            n.iter=500, n.burnin=100)
    # print("Have grabbed population model fit.")
    # preds.countries <- preds.countries.tmp$preds.countries
    # pred.pop <- preds.countries.tmp$pred.pop
    # print(paste0("Converting format of population trajectories, selecting first ",
    #              n.trajectories, " trajectories."))
    # preds.countries.trajs <- find.pop.trajectories(data.medium.tmp, pred.pop, 1:n.trajectories,
    #                                                c("USA", names.countries), year.start=year.start, year.end=year.end)
  }

  # These functions are from estimate_model_corr.R
  print("======================================")
  print("Combining IPAT projections using list.proj and estimate.co2.projections")
  print("")
  data.proj <- lapply(list.proj, function(x) convert.projections(x$TechData,
                                                                 x$GDPData, x$GDPFrontier, x$TechUSA))
  data.proj <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    tmp <- list()
    tmp$GDPFrontier <- x$GDPFrontier
    tmp$GDPData <- subset(x$GDPData, Isocode %in% preds.countries.trajs[[i]]$Isocode)
    tmp$TechData <- subset(x$TechData, Isocode %in% preds.countries.trajs[[i]]$Isocode)
    tmp$TechUSA <- x$TechUSA
    return (tmp)
  })
  co2.projections.tmp <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    preds.countries <- preds.countries.trajs[[i]]
    estimate.co2.projections(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                             preds.countries, year.range=c(year.start, year.end))
  })

  co2.annual.projections.tmp <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    preds.countries <- preds.countries.trajs[[i]]
    estimate.annual.co2.projections(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                                    preds.countries, year.range=c(year.start, year.end))
  })

  if (year.end >= 2100) {
    co2.annual.projections.restrictions <- sapply(co2.annual.projections.tmp, is.carbon.below.limit)
    # co2.projections.restrictions <- sapply(co2.projections.tmp, is.carbon.below.limit)
    co2.projections <- co2.projections.tmp[co2.annual.projections.restrictions]
    co2.annual.projections <- co2.annual.projections.tmp[co2.annual.projections.restrictions]

    data.proj.restrict <- data.proj[co2.annual.projections.restrictions]
    preds.countries.trajs.restrict <- preds.countries.trajs[co2.annual.projections.restrictions]
  } else {
    co2.projections.restrictions <- 1:n.trajectories
    co2.projections <- co2.projections.tmp
    co2.annual.projections.restrictions <- 1:n.trajectories
    co2.annual.projections <- co2.annual.projections.tmp
    data.proj.restrict <- data.proj
    preds.countries.trajs.restrict <- preds.countries.trajs
  }

  list(corr.model.out=corr.model.out,
       data.proj=data.proj.restrict,
       preds.countries.trajs=preds.countries.trajs.restrict,
       co2.projections=co2.projections,
       co2.annual.projections=co2.annual.projections,
       co2.projections.restrictions=co2.annual.projections.restrictions,
       names.countries.intensity=names.countries.intensity)
}


evaluate.projections <- function(data.medium.full, model.results,
                                 year.start=2010, year.end=2100,
                                 quantiles=c(0.025, 0.1, 0.5, 0.9, 0.975),
                                 outofsample.validate=FALSE,
                                 restrict.isos=FALSE,
                                 isos.restricted=c()) {
  # evaluate.projections uses the results of fit.project.model(), and
  # creates various summaries of the model and its projections by calling functions
  # like get.ipat.components.total, get.ipat.components.bycountry.
  # If year.start is 1980, 1990, or 2000, we would typically also pass
  # outofsample.validate=TRUE, and checks of the out of sample validation are
  # performed.

  if (year.start == 2010 || year.start == 2023) {
    data.proj <- model.results$data.proj
    preds.countries.trajs <- model.results$preds.countries.trajs
    co2.projections <- model.results$co2.projections
    co2.annual.projections <- model.results$co2.annual.projections
  } else {
    co2.proj.tmp <- model.results$co2.projections
    isos.include <- co2.proj.tmp[[1]][
      !is.na(co2.proj.tmp[[1]][, paste0("CO2", year.start)]), "Isocode"]
    data.proj <- lapply(model.results$data.proj, function(x) {
      lapply(x, function(y) {
        subset(y, Isocode %in% isos.include)
      })
    })
    preds.countries.trajs <- lapply(model.results$preds.countries.trajs,
                                    function(x)
                                      subset(x, Isocode %in% isos.include))
    co2.projections <- lapply(model.results$co2.projections, function(x)
      subset(x, Isocode %in% isos.include))
  }
  year.start.tmp <- year.start
  # Get worldwide data
  trajs.worldwide <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.total(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                              preds.countries.trajs[[i]], co2.projections[[i]],
                              year.sequence=seq(year.start, year.end,1))
  })
  trajs.annual.worldwide <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.total(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                              preds.countries.trajs[[i]], co2.annual.projections[[i]],
                              year.sequence=year.start:year.end)
  })
  trajs.quants <- get.trajs.quants(co2.projections,
                                   seq(year.start, year.end, by=1),
                                   quantiles=quantiles)
  trajs.annual.quants <- get.trajs.quants(co2.annual.projections,
                                          year.start:year.end,
                                          quantiles=quantiles)

  trajs.quants.bycountry <- get.trajs.quants.bycountry(co2.projections,
                                                       seq(year.start, year.end, 1),
                                                       quantiles=quantiles)
  trajs.annual.quants.bycountry <- get.trajs.quants.bycountry(co2.annual.projections,
                                                              year.start:year.end,
                                                              quantiles=quantiles)

  # Get data for Sub-Saharan Africa (SSA)
  trajs.ssa <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.ssa(x$TechData, x$GDPData,
                            preds.countries.trajs[[i]], co2.projections[[i]],
                            year.sequence=seq(year.start, year.end,1)) # seq(year.start, year.end,5)
  })

  ipat.quantiles <- get.ipat.medians(trajs.worldwide, quantiles=quantiles)
  ipat.annual.quantiles <- get.ipat.medians(trajs.annual.worldwide, quantiles=quantiles)

  ipat.ssa.quantiles <- get.ipat.medians(trajs.ssa, quantiles=quantiles)

  if (! outofsample.validate) {
    # stopifnot(year.start == 2010)
    # Worldwide cumulative emissions for each scenario,
    # and quantiles of cumulative emissions
    n.projections <- length(co2.projections)
    cum.emissions.tmp <- lapply(1:n.projections,
                                function(i) get.cumulative.emissions(co2.projections[[i]])$cum.emissions)
    cum.emissions <- data.frame(matrix(nrow=n.projections, ncol=length(cum.emissions.tmp[[1]])))
    cum.emissions.quantile <- data.frame(matrix(nrow=length(cum.emissions.tmp[[1]]), ncol=length(quantiles)))
    names(cum.emissions.quantile) <- paste0("Quant", quantiles)
    for (i in 1:n.projections) {
      cum.emissions[i,] <- cum.emissions.tmp[[i]]
    }
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.quantile[year.tmp,] <- quantile(cum.emissions[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.quantile$Year <- seq(year.start.tmp,year.end,1) #seq(year.start.tmp,year.end,5)
  }

  # Carbon emissions for Sub-Saharan Africa (SSA)
  if (year.start == 2010 || year.start == 2005 || year.start == 2023) {
    co2.ssa <- as.data.frame(matrix(nrow=length(co2.projections),
                                    ncol=length(names(co2.projections[[1]])[-1]),
                                    dimnames=list(c(),names(co2.projections[[1]])[-1])))
    which.ssa.rows <- co2.projections[[1]]$Isocode %in% ssa.isos
    for (i in 1:length(co2.projections)) {
      co2.ssa[i,] <- colSums(co2.projections[[i]][which.ssa.rows, -1])
    }
    # Carbon quantiles for SSA
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))
    co2.ssa.quants$Quantiles <- quantiles
    for (col.name in names(co2.ssa)) {
      co2.ssa.quants[, col.name] <- quantile(co2.ssa[, col.name], probs=quantiles) / 10^9
    }

    #GDP per capita for Sub-Saharan Africa (SSA)
    gdp.ssa <- pop.ssa <- as.data.frame(matrix(nrow=length(co2.projections),
                                               ncol=length(seq(year.start, year.end, 1)),
                                               dimnames=list(c(), paste0("GDP", seq(year.start, year.end, 1)))))
    which.ssa.rows.gdp <- data.proj[[1]][["GDPData"]]$Isocode %in% ssa.isos
    which.ssa.rows.pop <- preds.countries.trajs[[1]]$Isocode %in% ssa.isos

    gdp.names.tmp <- paste0("GDP", seq(year.start, year.end, 1))
    for (i in 1:length(co2.projections)) {
      gdp.ssa[i,] <- colSums(data.proj[[i]][["GDPData"]][which.ssa.rows.gdp, gdp.names.tmp] *
                               preds.countries.trajs[[i]][which.ssa.rows.pop, -1] * 1000)
    }

    # Carbon quantiles for SSA
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))
    co2.ssa.quants$Quantiles <- quantiles
    for (col.name in names(co2.ssa)) {
      co2.ssa.quants[, col.name] <- quantile(co2.ssa[, col.name], probs=quantiles) / 10^9
    }

    # GDP per capita quantiles for SSA
    gdp.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(gdp.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(gdp.ssa)))))
    gdp.ssa.quants$Quantiles <- quantiles
    for (col.name in names(gdp.ssa)) {
      gdp.ssa.quants[, col.name] <- quantile(gdp.ssa[, col.name], probs=quantiles)
    }
  } else {
    # Fill out dummy categories
    co2.ssa <- as.data.frame(matrix(nrow=length(co2.projections),
                                    ncol=length(names(co2.projections[[1]])[-1]),
                                    dimnames=list(c(),names(co2.projections[[1]])[-1])))
    co2.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(co2.ssa)))))

    gdp.ssa <- as.data.frame(matrix(nrow=length(co2.projections),
                                    ncol=length(seq(year.start, year.end, 1)),
                                    dimnames=list(c(), paste0("GDP", seq(year.start, year.end, 1)))))
    gdp.ssa.quants <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(gdp.ssa)[2]),
             dimnames=list(c(), c("Quantiles", names(gdp.ssa)))))
  }

  # Cumulative SSA emissions
  if (! outofsample.validate) {
    co2.cum.ssa <- co2.ssa # just to size the array
    for (i in 1:dim(co2.ssa)[1]) {
      co2.cum.ssa[i,] <- get.cumulative.emissions(cbind(Isocode="SSA", co2.ssa[i,]))$cum.emissions
    }
    cum.emissions.ssa.quantile <- data.frame(matrix(nrow=length(seq(year.start, year.end,1)), ncol=length(quantiles)))
    names(cum.emissions.ssa.quantile) <- paste0("Quant", quantiles)
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.ssa.quantile[year.tmp,] <- quantile(co2.cum.ssa[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.ssa.quantile$Year <- seq(year.start.tmp,year.end,1)
  }

  # Cumulative emissions, excluding SSA
  if (! outofsample.validate) {
    # stopifnot(year.start == 2010)
    # Worldwide cumulative emissions for each scenario,
    # and quantiles of cumulative emissions
    n.projections <- length(co2.projections)
    row.inds.notssa <- which(! co2.projections[[1]]$Isocode %in% ssa.isos)
    cum.emissions.tmp.notssa <- lapply(1:n.projections,
                                       function(i)
                                         get.cumulative.emissions(co2.projections[[i]][row.inds.notssa,])$cum.emissions)
    cum.emissions.notssa <- data.frame(matrix(nrow=n.projections, ncol=length(cum.emissions.tmp.notssa[[1]])))
    cum.emissions.notssa.quantile <- data.frame(matrix(nrow=length(cum.emissions.tmp.notssa[[1]]), ncol=length(quantiles)))
    names(cum.emissions.notssa.quantile) <- paste0("Quant", quantiles)
    for (i in 1:n.projections) {
      cum.emissions.notssa[i,] <- cum.emissions.tmp.notssa[[i]]
    }
    for (year.tmp in 1:dim(cum.emissions)[2]) {
      cum.emissions.notssa.quantile[year.tmp,] <- quantile(cum.emissions.notssa[, year.tmp], probs=quantiles) / 10^9
    }
    cum.emissions.notssa.quantile$Year <- seq(year.start.tmp,year.end,1)
  }


  # GDP split by OECD/non-OECD countries. Find quantiles
  if (! outofsample.validate) {
    data.proj <- model.results$data.proj
    pop.preds <- model.results$preds.countries.trajs
    which.oecd.not.us <- data.proj[[1]][["GDPData"]][, "Isocode"] %in% oecd.15.iso
    which.world.notoecd <- data.proj[[1]][["GDPData"]][, "Isocode"] %in% oecd.world.minusoecd15
    n.projections <- length(data.proj)
    stopifnot(pop.preds[[1]][1, "Isocode"] == "USA")
    # bc GDP is projected yearly but Pop isn't
    gdp.columns <- paste0("GDP", seq(year.start.tmp,year.end,1))
    gdp.by.oecd <- gdppercapita.by.oecd <- list()
    gdp.tot.names <- paste0("GDPTotal", seq(year.start.tmp,year.end,1))
    gdp.by.oecd[["OECD"]] <- gdp.by.oecd[["NotOECD"]] <-
      gdp.by.oecd[["WorldNotOECD"]] <-
      gdppercapita.by.oecd[["OECD"]] <- gdppercapita.by.oecd[["NotOECD"]] <- as.data.frame(
        matrix(nrow=n.projections, ncol=length(gdp.columns)))
    names(gdp.by.oecd[["OECD"]]) <- names(gdp.by.oecd[["NotOECD"]]) <-
      names(gdp.by.oecd[["WorldNotOECD"]]) <- gdp.tot.names
    names(gdppercapita.by.oecd[["OECD"]]) <- names(gdppercapita.by.oecd[["NotOECD"]]) <- gdp.tot.names
    for (i in 1:n.projections) {
      popdata.nonus <- pop.preds[[i]][-1,]
      # -1 to remove Isocode column
      gdp.oecd <- colSums(data.proj[[i]][["GDPData"]][which.oecd.not.us, gdp.columns] *
                            popdata.nonus[which.oecd.not.us, -1] * 1000) +
        data.proj[[1]][["GDPFrontier"]][, gdp.columns] * pop.preds[[i]][1,-1] * 1000
      gdp.worldnotoecd <- colSums(data.proj[[i]][["GDPData"]][which.world.notoecd, gdp.columns] *
                                    popdata.nonus[which.world.notoecd, -1] * 1000)
      pop.oecd <- colSums(popdata.nonus[which.oecd.not.us, -1] * 1000) +
        pop.preds[[i]][1,-1] * 1000
      gdp.not.oecd <- colSums(data.proj[[i]][["GDPData"]][! which.oecd.not.us, gdp.columns] *
                                popdata.nonus[! which.oecd.not.us, -1] * 1000)
      pop.not.oecd <- colSums(popdata.nonus[! which.oecd.not.us, -1] * 1000)
      gdp.by.oecd[["OECD"]][i, ] <- gdp.oecd
      gdp.by.oecd[["WorldNotOECD"]][i, ] <- gdp.worldnotoecd
      gdp.by.oecd[["NotOECD"]][i, ] <- gdp.not.oecd
      gdppercapita.by.oecd[["OECD"]][i, ] <- gdp.oecd / pop.oecd
      gdppercapita.by.oecd[["NotOECD"]][i, ] <- gdp.not.oecd / pop.not.oecd
    }
    gdp.by.oecd.quant <- gdppercapita.by.oecd.quant <- list()
    gdp.by.oecd.quant[["OECD"]] <- gdp.by.oecd.quant[["NotOECD"]] <-
      gdp.by.oecd.quant[["WorldNotOECD"]] <- as.data.frame(
        matrix(nrow=length(quantiles), ncol=(1 + length(gdp.columns))))
    gdppercapita.by.oecd.quant[["OECD"]] <- gdppercapita.by.oecd.quant[["NotOECD"]] <- as.data.frame(
      matrix(nrow=length(quantiles), ncol=(1 + length(gdp.columns))))
    names(gdp.by.oecd.quant[["OECD"]]) <- names(gdp.by.oecd.quant[["WorldNotOECD"]]) <-
      names(gdp.by.oecd.quant[["NotOECD"]]) <- c("Quantile", gdp.tot.names)
    names(gdppercapita.by.oecd.quant[["OECD"]]) <- names(gdppercapita.by.oecd.quant[["NotOECD"]]) <- c("Quantile", gdp.tot.names)
    gdp.by.oecd.quant[["OECD"]][, "Quantile"] <-
      gdp.by.oecd.quant[["WorldNotOECD"]][, "Quantile"] <-
      gdp.by.oecd.quant[["NotOECD"]][, "Quantile"] <- quantiles
    gdppercapita.by.oecd.quant[["OECD"]][, "Quantile"] <-
      gdppercapita.by.oecd.quant[["NotOECD"]][, "Quantile"] <- quantiles
    for (var.name in gdp.tot.names) {
      gdp.by.oecd.quant[["OECD"]][, var.name] <- quantile(gdp.by.oecd[["OECD"]][, var.name],
                                                          probs=quantiles)
      gdp.by.oecd.quant[["WorldNotOECD"]][, var.name] <- quantile(gdp.by.oecd[["WorldNotOECD"]][, var.name],
                                                                  probs=quantiles)
      gdp.by.oecd.quant[["NotOECD"]][, var.name] <- quantile(gdp.by.oecd[["NotOECD"]][, var.name],
                                                             probs=quantiles)
      gdppercapita.by.oecd.quant[["OECD"]][, var.name] <- quantile(gdppercapita.by.oecd[["OECD"]][, var.name],
                                                                   probs=quantiles)
      gdppercapita.by.oecd.quant[["NotOECD"]][, var.name] <- quantile(gdppercapita.by.oecd[["NotOECD"]][, var.name],
                                                                      probs=quantiles)
    }
  }

  # Carbon emissions by 5 IPCC regions
  co2.5regions <- list()
  # Uses global variables rcp.5region.names and reg.5.dataframe
  #rcp.5region.names <- c("ASIA", "LAM", "MAF", "OECD", "REF") # adding
  for (reg in rcp.5region.names) {
    if (reg == "REF" & year.start == 1980) next
    co2.5regions[[reg]] <- as.data.frame(matrix(nrow=length(co2.projections),
                                                ncol=length(names(co2.projections[[1]])[-1]),
                                                dimnames=list(c(),names(co2.projections[[1]])[-1])))
    reg.isos <- subset(reg.5.dataframe, reg == Reg5)$Isocode
    which.reg.rows <- co2.projections[[1]]$Isocode %in% reg.isos
    for (i in 1:length(co2.projections)) {
      # Remove the isocode column and sum up over each year
      if (year.start < 1990) {
        # Namibia has missing CO2 data from 1980 to 1989. Just exclude
        # them from the count
        co2.5regions[[reg]][i,] <- colSums(co2.projections[[i]][which.reg.rows, -1], na.rm=T)
      } else {
        co2.5regions[[reg]][i,] <- colSums(co2.projections[[i]][which.reg.rows, -1])
      }
    }
  }
  # Carbon quantiles by 5 IPCC regions
  co2.5regions.quants <- list()
  for (reg in rcp.5region.names) {
    if (reg == "REF" & year.start == 1980) next
    print(reg)
    co2.5regions.quants[[reg]] <- as.data.frame(
      matrix(nrow=length(quantiles),
             ncol=(1 + dim(co2.5regions[[reg]])[2]),
             dimnames=list(c(), c("Quantiles", names(co2.5regions[[reg]])))))
    co2.5regions.quants[[reg]]$Quantiles <- quantiles
    for (col.name in names(co2.5regions[[reg]])) {
      if (year.start < 1990) {
        # Namibia has missing CO2 data from 1980 to 1989. Just exclude
        # them from the count
        co2.5regions.quants[[reg]][, col.name] <-
          quantile(co2.5regions[[reg]][, col.name], probs=quantiles, na.rm=T)
      } else {
        co2.5regions.quants[[reg]][, col.name] <-
          quantile(co2.5regions[[reg]][, col.name], probs=quantiles)
      }
    }
  }

  # IPAT quantiles by country
  ipat.components.bycountry <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.components.bycountry(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                                  preds.countries.trajs[[i]], co2.projections[[i]],
                                  year.sequence=seq(year.start,year.end,1))
  })
  ipat.quantiles.bycountry <- get.ipat.quantiles.bycountry(ipat.components.bycountry,
                                                           quantiles=quantiles)

  ipat.annual.components.bycountry <- lapply(1:length(data.proj), function(i) {
    x <- data.proj[[i]]
    get.ipat.annual.components.bycountry(x$TechData, x$GDPData, x$GDPFrontier, x$TechUSA,
                                         preds.countries.trajs[[i]], co2.annual.projections[[i]],
                                         year.sequence=year.start:year.end)
  })

  ipat.annual.quantiles.bycountry <- get.ipat.quantiles.bycountry(ipat.annual.components.bycountry,
                                                                  quantiles=quantiles)

  # Look at gap from the frontier.
  frontier.gaps.pred <- list()
  for (iso in names(ipat.components.bycountry[[1]])) {
    if (iso != "USA") {
      frontier.gaps.pred[[iso]] <- list()
    }
  }
  for (i in 1:length(ipat.components.bycountry)) {
    us.gdp <- ipat.components.bycountry[[i]]$USA$GDPpercapita
    for (iso in names(ipat.components.bycountry[[1]])) {
      if (iso != "USA") {
        frontier.gaps.pred[[iso]][[i]] <-
          log(us.gdp / ipat.components.bycountry[[i]][[iso]]$GDPpercapita)
      }
    }
  }

  if (outofsample.validate) {
    # Warning: Assuming that the quantiles are symmetric!!!!
    stopifnot(all(abs((quantiles + rev(quantiles)) - 1) < 0.0001))
    stopifnot(year.end < 2023)
    n.intervals <- floor(length(quantiles)/2)
    tests.belowmedian <- prop.belowmedian <- list()
    tests.true <- prop.true.tests <- list()
    tests.true.ipat <- prop.true.tests.ipat <- list()
    n.tests.ipat <- list()
    for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
      tests.belowmedian[[ipat.var]] <- prop.belowmedian[[ipat.var]] <- list()
      tests.true.ipat[[ipat.var]] <- prop.true.tests.ipat[[ipat.var]] <- list()
      n.tests.ipat[[ipat.var]] <- list()
    }
    # Have to ignore the first year b/c it's deterministic
    year.index <- 1
    for (year in seq(year.start + 1, year.end, 1)) {
      year.index <- year.index + 1
      year.str <- as.character(year)
      n.tests <- 0
      tests.true[[year.str]] <- rep(0, n.intervals)
      for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
        tests.belowmedian[[ipat.var]][[year.str]] <- 0
        tests.true.ipat[[ipat.var]][[year.str]] <- rep(0, n.intervals)
        n.tests.ipat[[ipat.var]][[year.str]] <- 0
      }

      trajs.quants.bycountry.tmp <- trajs.quants.bycountry[[paste0("CO2", year)]]
      data.true.tmp <- subset(data.medium.full, Year == year)
      for (i in 1:dim(trajs.quants.bycountry.tmp)[1]) {
        iso <- trajs.quants.bycountry.tmp[i, "Isocode"]
        if ((! restrict.isos) | iso %in% isos.restricted) {
          co2.pred.tmp <- subset(trajs.quants.bycountry.tmp, Isocode == iso)
          data.tmp1 <- subset(data.true.tmp, Isocode == iso)
          co2.true <- data.tmp1$CO2 * (data.tmp1$PopTotal*10^3) / 10^9
          for (q_ind in 1:n.intervals) {
            co2.lower.name <- paste0("CO2_", quantiles[q_ind])
            co2.upper.name <- paste0("CO2_", quantiles[length(quantiles)+1-q_ind])
            co2.lower <- co2.pred.tmp[, co2.lower.name]
            co2.upper <- co2.pred.tmp[, co2.upper.name]
            if ((!is.na(co2.lower)) & (!is.na(co2.upper))) {
              # Need to check is.na b/c of 1980 validation
              stopifnot(co2.lower < co2.upper)
              tests.true[[year.str]][q_ind] <- tests.true[[year.str]][q_ind] +
                (co2.lower < co2.true & co2.true < co2.upper)
            }
          }
          n.tests <- n.tests + 1
        }
      }
      prop.true.tests[[year.str]] <- tests.true[[year.str]] / n.tests

      for (ipat.var in names(ipat.quantiles.bycountry[["USA"]][[1]])[-1]) {
        for (i in 1:dim(trajs.quants.bycountry.tmp)[1]) {
          iso <- as.character(trajs.quants.bycountry.tmp[i, "Isocode"])
          if ((! restrict.isos) | iso %in% isos.restricted) {
            data.tmp1 <- subset(data.true.tmp, Isocode == iso & Year == year)

            if (ipat.var == "GDPpercapita") {
              var.true <- data.tmp1["GDP"]
            } else if (ipat.var == "GDP") {
              var.true <- data.tmp1["GDP"] * (10^3 * data.tmp1["PopTotal"])
            } else if (ipat.var == "CO2") {
              var.true <- data.tmp1["CO2"] * (10^3 * data.tmp1["PopTotal"])
            } else if (ipat.var == "CO2percapita") {
              var.true <- data.tmp1["CO2"]
            } else if (ipat.var == "Pop") {
              var.true <- data.tmp1["PopTotal"]
            } else if (ipat.var == "FrontierGap") {
              next
            } else { # Tech
              var.true <- data.tmp1[ipat.var]
            }
            if (!is.na(var.true)) {
              q_med <- n.intervals+1
              var.median <- ipat.quantiles.bycountry[[iso]][[q_med]][year.index, ipat.var]
              if ((!is.na(var.true)) & (!is.na(var.median))) {
                # Have to have this check because of 1980 out of sample validation
                tests.belowmedian[[ipat.var]][[year.str]] <- tests.belowmedian[[ipat.var]][[year.str]] +
                  (var.true < var.median)
              }
              for (q_ind in 1:n.intervals) {
                var.lower.name <- paste0(ipat.var, "_", quantiles[q_ind])
                var.upper.name <- paste0(ipat.var, "_", quantiles[length(quantiles)+1-q_ind])
                var.lower <- ipat.quantiles.bycountry[[iso]][[q_ind]][year.index, ipat.var]
                var.upper <- ipat.quantiles.bycountry[[iso]][[length(quantiles)+1-q_ind]][year.index, ipat.var]
                if ((!is.na(var.lower)) & (!is.na(var.upper))) {
                  # Need to check is.na b/c of 1980 validation
                  stopifnot(var.lower < var.upper)
                  tests.true.ipat[[ipat.var]][[year.str]][q_ind] <- tests.true.ipat[[ipat.var]][[year.str]][q_ind] +
                    (var.lower < var.true & var.true < var.upper)
                }
              }
              if ((!is.na(var.true)) & (!is.na(var.median))) {
                n.tests.ipat[[ipat.var]][[year.str]] <- n.tests.ipat[[ipat.var]][[year.str]] + 1
              }
            }
          }
        }
        prop.true.tests.ipat[[ipat.var]][[year.str]] <- tests.true.ipat[[ipat.var]][[year.str]] /
          n.tests.ipat[[ipat.var]][[year.str]]
        prop.belowmedian[[ipat.var]][[year.str]] <- tests.belowmedian[[ipat.var]][[year.str]] /
          n.tests.ipat[[ipat.var]][[year.str]]
      }
    }
  }

  if (outofsample.validate) {
    return(list(trajs.worldwide=trajs.worldwide, trajs.quants=trajs.quants, trajs.annual.quants=trajs.annual.quants,
                trajs.quants.bycountry=trajs.quants.bycountry,
                trajs.ssa=trajs.ssa,
                ipat.ssa.quantiles=ipat.ssa.quantiles,
                co2.ssa=co2.ssa,
                co2.ssa.quants=co2.ssa.quants,
                gdp.ssa=gdp.ssa,
                gdp.ssa.quants=gdp.ssa.quants,
                co2.5regions=co2.5regions,
                co2.5regions.quants=co2.5regions.quants,
                ipat.quantiles=ipat.quantiles,
                ipat.components.bycountry=ipat.components.bycountry,
                frontier.gaps.pred=frontier.gaps.pred,
                ipat.quantiles.bycountry=ipat.quantiles.bycountry,
                quantiles=quantiles,
                prop.true.tests=prop.true.tests,
                prop.true.tests.ipat=prop.true.tests.ipat,
                prop.belowmedian=prop.belowmedian))
  } else {
    return(list(trajs.worldwide=trajs.worldwide, trajs.quants=trajs.quants, trajs.annual.quants=trajs.annual.quants,
                trajs.annual.worldwide=trajs.annual.worldwide,
                trajs.quants.bycountry=trajs.quants.bycountry,
                trajs.ssa=trajs.ssa,
                ipat.ssa.quantiles=ipat.ssa.quantiles,
                cum.emissions=cum.emissions, cum.emissions.quantile=cum.emissions.quantile,
                co2.ssa=co2.ssa,
                co2.ssa.quants=co2.ssa.quants,
                gdp.ssa=gdp.ssa,
                gdp.ssa.quants=gdp.ssa.quants,
                co2.cum.ssa=co2.cum.ssa,
                cum.emissions.ssa.quantile=cum.emissions.ssa.quantile,
                cum.emissions.notssa=cum.emissions.notssa,
                cum.emissions.notssa.quantile=cum.emissions.notssa.quantile,
                gdp.by.oecd=gdp.by.oecd,
                gdp.by.oecd.quant=gdp.by.oecd.quant,
                gdppercapita.by.oecd=gdppercapita.by.oecd,
                gdppercapita.by.oecd.quant=gdppercapita.by.oecd.quant,
                co2.5regions=co2.5regions,
                co2.5regions.quants=co2.5regions.quants,
                ipat.quantiles=ipat.quantiles,
                ipat.components.bycountry=ipat.components.bycountry,
                ipat.annual.quantiles=ipat.annual.quantiles,
                ipat.annual.components.bycountry=ipat.annual.components.bycountry,
                frontier.gaps.pred=frontier.gaps.pred,
                ipat.quantiles.bycountry=ipat.quantiles.bycountry,
                ipat.annual.quantiles.bycountry=ipat.annual.quantiles.bycountry,
                quantiles=quantiles))
  }
}

# ==========================================
library(readxl)
data.location <- "NatureData/"
sims.location <- paste0(data.location, "Simulations/")
plot.location <- "NatureData/Plots/"

# Get RCP Data
# rcp.list <- list()
# for (rcp.num in c("RCP26", "RCP45", "RCP60", "RCP85")) {
#   rcp.list[[rcp.num]] <- read.csv(paste0(data.location, "IPCC_", rcp.num, "_data.csv"))
# }
rcp.colors <- c("green", "red", "black", "purple")

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
carbon.cum <- rep(0, 4)
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
print(xtable(rcp.carbon.cum), include.rownames=F)
print(xtable(rcp.carbon.yearly), include.rownames=F)

# RCP 8.5, yearly worldwide emissions of CO2 (fossil fuels and industry).
# From https://tntcat.iiasa.ac.at/AR5DB/dsd?Action=htmlpage&page=regions
# Converting to gigatonnes of carbon
# rcp.8.5.yearly <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
#                              Carbon=c(29227.000, 32727.200, 42304.533, 50743.000, 61550.867, 74083.533, 86519.400, 95194.733, 100489.033, 103901.233, 105380.367)
#                              / ((11/3) * 10^3))
# # Starting RCP 8.5 from 2010
# rcp.8.5.cumulative <- data.frame(Year=c(2010, seq(2023, 2090, by=10), 2100),
#                                  Carbon=NA)
# cum.emissions.tmp <- 0
# rcp.8.5.cumulative$Carbon[1] <- 0
# cum.emissions.tmp <- 5*rcp.8.5.yearly$Carbon[2] # 2010 emissions until 2023
# rcp.8.5.cumulative$Carbon[2] <- cum.emissions.tmp
# # add on 5 years for 2010 start
# for (i in 3:dim(rcp.8.5.cumulative)[1]) {
#   if (rcp.8.5.cumulative$Year[i] == 2100) {
#     cum.emissions.tmp <- cum.emissions.tmp + 5*rcp.8.5.yearly$Carbon[i]
#   } else {
#     cum.emissions.tmp <- cum.emissions.tmp + 10*rcp.8.5.yearly$Carbon[i]
#   }
#   rcp.8.5.cumulative$Carbon[i] <- cum.emissions.tmp
# }
#
# # RCP 8.5 GDP data, in billions of 2005USD per year
# rcp.8.5.yearly.gdp <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
#                                  GDP=c(42987.656, 48302.922, 63588.657, 83849.204, 110376.903, 139114.084, 165405.548, 190946.622, 215806.785, 239972.143, 262928.537))
# # RCP 8.5 Tech data, tonnes of carbon per $10,000
# rcp.8.5.yearly.tech <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
#                                   Tech=(rcp.8.5.yearly$Carbon / rcp.8.5.yearly.gdp$GDP * 10^4))
# # RCP 8.5 Population data, in millions of people in a year
# rcp.8.5.pop <- data.frame(Year=c(2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100),
#                           PopTotal=c(6469.985, 6884.830, 7814.060, 8712.230, 9548.040, 10245.410, 10864.550, 11392.090, 11847.140, 12202.210, 12386.320))

# OECD world GDP forecasts, taken from https://data.oecd.org/gdp/gdp-long-term-forecast.htm
# In millions of USD per year (2010USD)
# oecd.gdp.tmp <- read.csv(paste0(data.location, "oecd_worldgdp_forecast.csv"), # not be used in the following
#                          sep=",", stringsAsFactors=FALSE)
# oecd.gdp <- oecd.gdp.tmp[, c("TIME", "Value")]
# names(oecd.gdp) <- c("Year", "GDP")
# oecd.gdp$GDP <- as.numeric(oecd.gdp$GDP)

# Load dataset:
load(file=paste0(data.location, "data_medium_2024.Rda"))
data.medium <- data_medium_new
# For help in calling the plot function for different countries
countrylist <- as.character(unique(data.medium$Isocode))
isolist <- countrycode(sourcevar=countrylist, origin="iso3c", destination="country.name")
cbind(countrylist, isolist)

getIso <- function(country) {
  library(countrycode)
  # take the first element in case we have a vector of strings
  if (nchar(country[1]) == 3) return(country)
  else return(countrycode(sourcevar=country, origin="country.name", destination="iso3c"))
}

getCountry <- function(iso) {
  library(countrycode)
  # take the first element in case we have a vector of strings
  if (nchar(iso[1]) > 3) return(iso)
  else return(countrycode(sourcevar=iso, origin="iso3c", destination="country.name"))
}

countries.isos <- unique(data.medium$Isocode)


countries.eu <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus",
                  "Czech Republic", "Denmark", "Estonia", "Finland", "France",
                  "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia",
                  "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland",
                  "Portugal", "Romania", "Slovakia", "Slovenia", "Spain",
                  "Sweden", "United Kingdom")
isos.eu <- sapply(countries.eu, getIso)

rcp.5region.tmp <- read_xlsx(paste0(data.location, 'ssp_db_5regions.xlsx'),
                            n_max=20)
rcp.5region <- rcp.5region.tmp[, c("Region", "Scenario", "Unit",
                                   c(2000, 2005, seq(2010, 2100, by=10)))]
names(rcp.5region)[4:length(names(rcp.5region))] <- paste0("Carbon", names(rcp.5region)[4:length(names(rcp.5region))])


rcp.5region[,-c(1:3)] <- (1/1000) * rcp.5region[,-c(1:3)]
rcp.5region$Unit <- "PgCO2/yr"
rcp.5region.names <- c("ASIA", "LAM", "MAF", "OECD", "REF")

# World population data:
library(wpp2022)
data(UNlocations)
UNlocations$Isocode <- sapply(as.character(UNlocations$name), getIso)


#=========================================================


oecd.world.iso <- c("AUS", "AUT", "BEL", "BRA", "CAN", "CHE", "CHL", "CHN",
                    "CZE", "DEU", "DNK", "ESP", "EST", "FIN", "FRA",
                    "GBR", "GRC", "HUN", "IDN", "IND", "IRL", "ISL", "ISR",
                    "ITA", "JPN", "KOR", "LUX", "MEX", "NLD", "NOR", "NZL",
                    "POL", "PRT", "RUS", "SVK", "SVN", "SWE", "TUR",
                    "USA", "ZAF")
oecd.world.minusoecd15 <- c("BRA", "CHN", "IDN", "IND", "RUS", "ZAF")
oecd.15.iso <- oecd.world.iso[! oecd.world.iso %in% oecd.world.minusoecd15]

# ========================================================

# Sub-Saharan Africa
ssa.countries <- c("Angola", "Benin", "Botswana",
                   "Burkina Faso", "Burundi", "Cameroon", "Cape Verde", "Central African Republic",
                   "Chad", "Comoros", "Congo, Republic of", "Congo, the Democratic Republic of the", "C么te d'Ivoire",
                   "Djibouti", "Equatorial Guinea", "Eritrea", "Ethiopia", "Gabon", "The Gambia", "Ghana",
                   "Guinea", "Guinea-Bissau", "Kenya", "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali",
                   "Mauritania", "Mauritius", "Mozambique", "Namibia", "Niger", "Nigeria", "Réunion", "Rwanda",
                   "Sao Tome and Principe", "Senegal", "Seychelles", "Sierra Leone", "Somalia", "South Africa", "Sudan",
                   "Swaziland", "Tanzania", "Togo", "Uganda", "Western Sahara", "Zambia", "Zimbabwe")
ssa.isos <- getIso(ssa.countries)


get.rcp.region.isos <- function(countries.str) {
  # get.rcp.region.isos splits a string of countries, gets isocodes
  library(countrycode)
  return(sapply(strsplit(countries.str, ", "), getIso))
}


# IPCC RCP regions
#OECD90 = Includes the OECD 90 countries, therefore encompassing the countries included in the regions Western Europe (Austria, Belgium, Denmark, Finland, France, Germany, Greece, Iceland, Ireland, Italy, Luxembourg, Netherlands, Norway, Portugal, Spain, Sweden, Switzerland, Turkey, United Kingdom), Northern America (Canada, United States of America) and Pacific OECD (Australia, Fiji, French Polynesia, Guam, Japan, New Caledonia, New Zealand, Samoa, Solomon Islands, Vanuatu) .
oecd.90.countries.str <- "Austria, Belgium, Denmark, Finland, France, Germany, Greece, Iceland, Ireland, Italy, Luxembourg, Netherlands, Norway, Portugal, Spain, Sweden, Switzerland, Turkey, United Kingdom, Canada, United States of America, Australia, Fiji, French Polynesia, Guam, Japan, New Caledonia, New Zealand, Samoa, Solomon Islands, Vanuatu"
oecd.90.iso <- get.rcp.region.isos(oecd.90.countries.str)
#REF = Countries from the Reforming Ecomonies region (Albania, Armenia, Azerbaijan, Belarus, Bosnia and Herzegovina, Bulgaria, Croatia, Cyprus, Czech Republic, Estonia, Georgia, Hungary, Kazakhstan, Kyrgyzstan, Latvia, Lithuania, Malta, Poland, Republic of Moldova, Romania, Russian Federation, Slovakia, Slovenia, Tajikistan, TFYR Macedonia, Turkmenistan, Ukraine, Uzbekistan, Yugoslavia).
ref.countries.str <- "Albania, Armenia, Azerbaijan, Belarus, Bosnia and Herzegovina, Bulgaria, Croatia, Cyprus, Czech Republic, Estonia, Georgia, Hungary, Kazakhstan, Kyrgyzstan, Latvia, Lithuania, Malta, Poland, Republic of Moldova, Romania, Russian Federation, Slovakia, Slovenia, Tajikistan, TFYR Macedonia, Turkmenistan, Ukraine, Uzbekistan, Yugoslavia"
ref.iso <- get.rcp.region.isos(ref.countries.str)
#ASIA = The countries included in the regions China + (China, China Hong Kong SAR, China Macao SAR, Mongolia, Taiwan) , India + (Afghanistan, Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka) and Rest of Asia (Brunei Darussalam, Cambodia, Democratic People's Republic of Korea, East Timor, Indonesia, Lao People's Democratic Republic, Malaysia, Myanmar, Papua New Guinea, Philippines, Republic of Korea, Singapore, Thailand, Viet Nam) are aggregated into this region.
asia.countries.str <- "China, China Hong Kong SAR, China Macao SAR, Mongolia, Taiwan, India, Afghanistan, Bangladesh, Bhutan, India, Maldives, Nepal, Pakistan, Sri Lanka, Brunei Darussalam, Cambodia, Democratic People's Republic of Korea, East Timor, Indonesia, Lao People's Democratic Republic, Malaysia, Myanmar, Papua New Guinea, Philippines, Republic of Korea, Singapore, Thailand, Viet Nam"
asia.iso <- get.rcp.region.isos(asia.countries.str)
#MAF = This region includes the Middle East (Bahrain, Iran (Islamic Republic of), Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Qatar, Saudi Arabia, Syrian Arab Republic, United Arab Emirates, Yemen) and African (Algeria, Angola, Benin, Botswana, Burkina Faso, Burundi, Cote d'Ivoire, Cameroon, Cape Verde, Central African Republic, Chad, Comoros, Congo, Democratic Republic of the Congo, Djibouti, Egypt, Equatorial Guinea, Eritrea, Ethiopia, Gabon, Gambia, Ghana, Guinea, Guinea-Bissau, Kenya, Lesotho, Liberia, Libyan Arab Jamahiriya, Madagascar, Malawi, Mali, Mauritania, Mauritius, Morocco, Mozambique, Namibia, Niger, Nigeria, Reunion, Rwanda, Senegal, Sierra Leone, Somalia, South Africa, Sudan, Swaziland, Togo, Tunisia, Uganda, United Republic of Tanzania, Western Sahara, Zambia, Zimbabwe) countries.
maf.countries.str <- "Bahrain, Iran (Islamic Republic of), Iraq, Israel, Jordan, Kuwait, Lebanon, Oman, Qatar, Saudi Arabia, Syrian Arab Republic, United Arab Emirates, Yemen, Algeria, Angola, Benin, Botswana, Burkina Faso, Burundi, Cote d'Ivoire, Cameroon, Cape Verde, Central African Republic, Chad, Comoros, Congo, Democratic Republic of the Congo, Djibouti, Egypt, Equatorial Guinea, Eritrea, Ethiopia, Gabon, Gambia, Ghana, Guinea, Guinea-Bissau, Kenya, Lesotho, Liberia, Libyan Arab Jamahiriya, Madagascar, Malawi, Mali, Mauritania, Mauritius, Morocco, Mozambique, Namibia, Niger, Nigeria, Reunion, Rwanda, Senegal, Sierra Leone, Somalia, South Africa, Sudan, Swaziland, Togo, Tunisia, Uganda, United Republic of Tanzania, Western Sahara, Zambia, Zimbabwe"
maf.iso <- get.rcp.region.isos(maf.countries.str)
#LAM = This region includes the Latin American countries (Argentina, Bahamas, Barbados, Belize, Bolivia, Brazil, Chile, Colombia, Costa Rica, Cuba, Dominican Republic, Ecuador, El Salvador, Guadeloupe, Guatemala, Guyana, Haiti, Honduras, Jamaica, Martinique, Mexico, Netherlands Antilles, Nicaragua, Panama, Paraguay, Peru, Puerto Rico, Suriname, Trinidad and Tobago, Uruguay, Venezuela).
lam.countries.str <- "Argentina, Bahamas, Barbados, Belize, Bolivia, Brazil, Chile, Colombia, Costa Rica, Cuba, Dominican Republic, Ecuador, El Salvador, Guadeloupe, Guatemala, Guyana, Haiti, Honduras, Jamaica, Martinique, Mexico, Netherlands Antilles, Nicaragua, Panama, Paraguay, Peru, Puerto Rico, Suriname, Trinidad and Tobago, Uruguay, Venezuela"
lam.iso <- get.rcp.region.isos(lam.countries.str)


oecd.world.minusoecd90 <- oecd.world.iso[! oecd.world.iso %in% oecd.90.iso]

reg.5.dataframe <- data.frame(Isocode=c(oecd.90.iso, ref.iso, asia.iso, maf.iso, lam.iso),
                              Reg5=c(rep("OECD", length(oecd.90.iso)),
                                     rep("REF", length(ref.iso)),
                                     rep("ASIA", length(asia.iso)),
                                     rep("MAF", length(maf.iso)),
                                     rep("LAM", length(lam.iso))))

rcp.5region.names <- c("ASIA", "LAM", "MAF", "OECD", "REF")
rcp.5region.names.altorder <- c("OECD", "ASIA", "LAM", "MAF", "REF")


# Note: we're just ignoring missing data here, which is just a few countries.
olddata.5region <- data.frame(Region=rcp.5region.names)
for (year in 1960:2010) {
  var.name <- paste0("Carbon", year)
  olddata.5region[, var.name] <- rep(NA, 5)
  i.row <- 0
  for (reg in rcp.5region.names) {
    i.row <- i.row + 1
    reg.isos <- subset(reg.5.dataframe, Reg5 == reg)$Isocode
    data.tmp <- subset(data.medium, Isocode %in% reg.isos & Year == year)
    olddata.5region[i.row, var.name] <- sum(data.tmp$CO2 * data.tmp$PopTotal * 1000,
                                            na.rm=T)
  }
}

#========================================================

year.start <- 2023
year.end <- 2100
n.trajectories <- 1000

#preds.countries.tmp <- predict.population(year.present=year.start, year.end=year.end, make.new=TRUE)
preds.countries.tmp <- predict.population(year.present=2021, year.end=year.end,
                                          make.new=TRUE) #2021 and 2022 were deleted under the function

pred.pop <- preds.countries.tmp$pred.pop
names.countries.tmp <- as.character(unique(subset(data.medium, Isocode != "USA")$Isocode))
preds.countries.trajs <- find.pop.trajectories(data.medium, pred.pop, n.trajectories,
                                               c("USA", names.countries.tmp), year.start=year.start, year.end=year.end)
save(preds.countries.trajs, file=paste0(data.location, paste0("poppreds_formatted_", year.start, "_", year.end, ".rda")))
load(file=paste0(data.location, paste0("poppreds_formatted_", year.start, "_", year.end, ".rda")))

# =======================================

model.ar1.const.results.2023 <- fit.project.model(data.medium, year.start = 2023, year.end = 2100, fit.ar1.const = T)

save(model.ar1.const.results.2023, file=paste0(data.location, "model_results_ar1const_2023.rda"))
load(file=paste0(data.location, "model_results_ar1const_2023.rda"))

proj.evals.2023.ar1.const <- evaluate.projections(data.medium, model.ar1.const.results.2023,
                                                  year.start=2023, year.end=2100,
                                                  outofsample.validate=F,
                                                  quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975))
save(proj.evals.2023.ar1.const, file=paste0(data.location, "proj_evals_ar1const_2023.rda"))
load(file=paste0(data.location, "proj_evals_ar1const_2023.rda"))

