# libraries
library(tidyverse)
library(sf)
library(taxize)
#---get movebank study data----
#The data from all available studies is downloaded as a csv with individual
#information for each study. We import them and analyze each one

# load information about movebank studies. These were downloaded previously.
# They are all the individual information from all studies that are freely
# available
flist <- list.files("../moveSTIR tests project/movebank/individual_info/", full.names = TRUE)
studies_info <- read.csv("data/movebank_studies_info.csv")

# list all the species in any study. Nearly 400 different species
all_species <- sapply(flist, function(x) {unique(read.csv(x)[,'taxon_canonical_name'])}) |> unlist()|> unique()

# There are some studies that do not list the species
# sapply(templist, function(x) unique(x[,'taxon_canonical_name']))

# create database list with taxonomic information for all
taxa_info <- classification(all_species, db = 'ncbi')
# Transform the list to a database with only class, order, genus, family and sp.
taxa_info <- do.call(rbind, taxa_info) %>% rownames_to_column(var = "taxon_canonical_name") %>% 
  mutate(taxon_canonical_name = gsub(pattern = "[[:punct:]]+[[:digit:]]+",replacement = "", taxon_canonical_name)) %>%
  filter(rank %in% c("class", "order", "family", "genus", "species")) %>% 
  select(1:3) %>% pivot_wider(values_from = name, names_from = rank)
rm(all_species)
# Export file with information about species
write.csv(taxa_info, "data/taxonomic_info.csv", quote = FALSE, row.names = FALSE)
taxa_info <- read.csv("data/taxonomic_info.csv")

##---Fill missing species (mammals only for now)----

# Check studies for which I have access that do not have the species, and read
# their description to fill in the species. I have to see if these match those
# studies that don't have the taxon in the individual info file either.
studies_info %>% left_join(taxa_info, by = c("taxon_ids" = "taxon_canonical_name")) %>% 
  filter(i_have_download_access == "true", taxon_ids=='', is_test == 'false') %>% 
  select(id, acknowledgements, study_objective, sensor_type_ids) %>% view()
# There are only a few studies with mammals. I think it might be safe to assume
# that those that have no description at all will not be useful
data.frame(id =                  c(1268085753, 2200106206,         1707353581,1397833814, 1736857677, 1863036795, 1863036795),
           taxon_canonical_ids = c('Petrogale','Vulpes rueppellii','Tayassu pecari', 'Cebus capucinus', 'Felis catus','Panthera leo','Crocuta crocuta'))
# Select only studies I can download, of mammals, of GPS/Argos
studies_info_filtered <- studies_info %>% left_join(taxa_info, by = c("taxon_ids" = "taxon_canonical_name")) %>% 
  filter(i_have_download_access == "true", 
         grepl("GPS", sensor_type_ids) | grepl("Argos", sensor_type_ids)) %>% 
  filter(number_of_individuals>=10, # studies with enough individuals
         class == "Mammalia", # mammals only
         order != "Chiroptera", # remove bats
         !family %in% c("Phocidae", "Balaenopteridae", "Physeteridae", "Odobenidae", "Balaenidae", "Delphinidae"))# remove marine mammals
# This gives 50 studies. 
# names to import files
flist <- paste0("../moveSTIR tests project/movebank/individual_info/",studies_info_filtered$id, ".csv")

# Create pdf to export graphs with tracking periods.The plots show the period for every
# individual, to find periods of overlap. I am selecting a three-month period
# for each study
pdf(file = "study_periods.pdf")
# Function to create plots for all files. 
plot_times <- function(x) {
  dat <- read.csv(x)
  dat$timestamp_start <- as.POSIXct(strptime(dat$timestamp_start, format = "%Y-%m-%d %H:%M:%S"))
  dat$timestamp_end <- as.POSIXct(strptime(dat$timestamp_end, format = "%Y-%m-%d %H:%M:%S"))
  id <- strsplit(basename(x),split = "\\.")[[1]][1]
  sp <- studies_info_filtered[studies_info_filtered$id==id,"taxon_ids"]
  ggplot(dat, aes(xmin = timestamp_start, xmax = timestamp_end, y = factor(id)))+
    geom_linerange()+
    labs(title = paste(id,sp))+
    theme(axis.text.x = element_text(angle = -45, hjust = 0))
}
lapply(flist, plot_times)
dev.off()

# Based on manual examination of the pdf figures, I excluded studies with little
# overlap, with short periods, and studies of humans. I also selected a 3-month
# period for the remaining ones. I save this information in a csv file.
# studies_info_filtered <- studies_info_filtered %>% filter(!id %in% c(1294989952, 195130114, 64283289, 194854525,10857031, # short term domestic cat studies
#                                                                      334874110,194626601,193984609, 1666470874, 47449884, 300222057, # short term domestic cat studies
#                                                                      412206724, 7023252, 53460105, # cow, babboon (x2), short term
#                                                                      416289710, # red deer, apparent test
#                                                                      193984609, 243688044, # leopard, rabbit, little overlap
#                                                                      399246220, 406811869, 958426706 # human studies
#                                                                      )
#                                                           )
# Import the csv file
movestir_periods <- read.csv("data/movestir_study_periods.csv")
# These dates have to be transformed to the format yyyyMMddHHmmssSSS to include
# in movebank queries
movestir_periods$start <- format(strptime(movestir_periods$start, "%m/%d/%Y", tz = "GMT"), format = "%Y%m%d%H%M%S000")
movestir_periods$end <- format(strptime(movestir_periods$end, "%m/%d/%Y", tz = "GMT"), format = "%Y%m%d%H%M%S000")
head(movestir_periods)

# With this information I can create queries to download data from specific
# studies in the given periods.
### This needs to be adapted to account for licences.###
for (i in 1:nrow(movestir_periods)) {
  id = movestir_periods[i,"id"]
  ts = movestir_periods[i,"start"]
  te = movestir_periods[i,"end"]
  curl_req <- paste0("https://www.movebank.org/movebank/service/direct-read?entity_type=event&study_id=",id,
                     "&timestamp_start=",ts,
                     "&timestamp_end=",te,
                     "&api-token=74aec176-be5f-4941-85f4-840071b15022")
  # Check the header to see if it's an html file describing the licence
  # agreement
  is_html <- grepl("html",readLines(curl_req,1))
  # accept the license if required and download the data.
  if (is_html) { 
    curl_download(curl_req, destfile = "license_terms.txt")
    checksum <- tools::md5sum("license_terms.txt")
    curl_req <- paste0("https://www.movebank.org/movebank/service/direct-read?entity_type=event&study_id=",id,"&api-token=74aec176-be5f-4941-85f4-840071b15022&license-md5=",checksum)
    curl_download(curl_req, handle = h, destfile = paste0("data/movebank_tracks/",i,".csv"))
  } else { # if it's not a license file
    curl::curl_download(curl_req, destfile = paste0("data/movebank_tracks/",id,".csv"))
  }
}

#### Check files ####
dat <- read.csv("data/movebank_tracks/1241071371.csv")
head(dat)

# 
# ## Summary of all available studies
# study_summary <- indiv_info_db %>% mutate(start = lubridate::ymd_hms(timestamp_start),
#                                           end = lubridate::ymd_hms(timestamp_end),
#                                           trk_len = as.numeric(difftime(end,start,units = "days"))) %>% 
#   group_by(study_id, taxon_canonical_name) %>% summarise(n_indiv = n(), 
#                                                          min = min(trk_len),
#                                                          max = max(trk_len), 
#                                                          median_len = median(trk_len)) %>% 
#   filter(n_indiv>=10) %>% # filter to only studies with >= 10 individuals of a single species
#   left_join(taxa_class_db)
# 
# dat <- read.csv(flist[1])
# head(dat)
# # The species in the study
# (spp <- unique(dat$taxon_canonical_name))
# # There are more than 600 studies. Let's create a database with key information about the studies
# for (f in flist) {
#   
# }
# dat <- read.csv()
# movebank_db <- read.csv("../../moveSTIR tests project/movebank/studies_summary.csv")
# head(movebank_db)


# Check if the study is already downloaded, download the data if not
downloaded_file_ids <- list.files("../../moveSTIR tests project/movebank/position_data/", pattern = "[[:digit:]]+.csv") |>strsplit(split = "\\.")|> sapply(\(x) x[1])
studies_to_download <- studies[!studies %in% downloaded_file_ids]


# use the previous commands to see if all studies are now available
downloaded_file_ids <- list.files("../../moveSTIR tests project/movebank/position_data/", pattern = "[[:digit:]]+.csv") |>strsplit(split = "\\.")|> sapply(\(x) x[1])
studies[!studies %in% downloaded_file_ids]
# there are still three left, but it should be good for now


# get median fix rate
dat %>% group_by(individual_ID) %>% mutate(datetime = strptime(date_time, format = "%Y-%m-%d %H:%M:%S"))
head(dat)


# Plot trajectories
colrs <- hcl.colors(length(uid))
plot(dat$location_long, dat$location_lat,type = 'n')
for (i in seq_along(uid)) {
  id = uid[i]
  d = dat[dat$tag_id == id,]
  lines(d$location_long, d$location_lat, col = colrs[i])
}

# Function to transform data in movebank format to the expected input for
# movestir. 
convert_to_movestir <- function(x) {
  dat <- read.csv(x)
  # remove rows with no coordinates
  dat <- dat[complete.cases(dat),]
  study_id <- sub(".csv","",basename(x))
  # get utm zone for transforming to planar coordinates
  utmzone <- median(floor((180+dat$location_long)/6))
  newcrs <- paste0("+proj=utm +zone=",utmzone," +datum=WGS84")
  
  sfg <- st_multipoint(cbind(dat$location_long, dat$location_lat))
  sfc <- st_sfc(sfg, crs = 4326) %>% st_transform(crs = newcrs)
  dat$UTMx <- st_coordinates(sfc)[,"X"]
  dat$UTMy <- st_coordinates(sfc)[,"Y"]
  # rename
  dat <- dat[,c("tag_id", "timestamp", "UTMx", "UTMy")]
  names(dat)[c(1,2)] <- c("individual_ID", "date_time")
  outfile <- paste0("data/",study_id,"_movements.csv")
  write.csv(dat, outfile)
}

# Apply function to files from movebank
for (i in list.files("data/movebank_tracks/", full.names = TRUE)) {
  convert_to_movestir(i)
}
sapply(list.files("data/movebank_tracks/", full.names = TRUE), convert_to_movestir)





rm(list = ls())
