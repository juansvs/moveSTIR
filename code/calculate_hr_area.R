library(data.table)
library(adehabitatHR)
library(raster)
library(sp)

## Script calculates the area of the 95% utilization distribution for the pigs used in 
## the MoveSTIR study.
##
## We then use the Average FOI (converted to a raster) generated from MoveSTIR to estimate
## how average FOI on the landscape is distributed within the home range of pigs.
##
## Script also calculates and home range overlap between pigs to compare the implications of
## ignoring fine-scale heterogeneity in space us on epidemilogical dynamics.

dat = fread("../data/pig_movements.csv")
dat[, datetime:=as.POSIXct(strptime(date_time, format="%Y-%m-%d %H:%M:%S", tz="GMT"))]
dat[, t:=as.numeric(datetime)]
unq_collars = sort(unique(dat$individual_ID))

# Get home ranges
home_ranges = array(NA, dim=length(unq_collars))
my_hrs = list()
for(i in 1:length(unq_collars)){


	id = unq_collars[i]
	tdat = dat[id == individual_ID, .(UTMx, UTMy)]
	sp_dat = SpatialPoints(tdat)
	crs(sp_dat) = "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

	tkern = kernelUD(sp_dat, h="href")
	hr = kernel.area(tkern, percent=95, unin="m", unout="m")
	my_hrs[[id]] = getverticeshr(tkern)
	home_ranges[i] = hr

}

hr_dat = as.data.table(data.frame(pig=unq_collars, hr_m2=home_ranges))
fwrite(hr_dat, "../data/home_ranges_for_pigs.csv")

# Load the FOI raster. 
# NOTE: You first need to run moveSTIR_pig_movements.ipynb to generate foi_raster.tif

if(file.exists("../data/foi_raster.tif")){

	foi = raster("../data/foi_raster.tif")

	# Average FOI in the full area 
	ratios = list()
	for(id in unq_collars){

		mras = mask(foi, my_hrs[[id]])
		vals = values(mras)
		inarea = vals[!is.na(vals)]
		lorenz = cumsum(sort(inarea, decreasing=T)) / sum(inarea)
		inds = which(lorenz < 0.8)
		num = inds[length(inds)]
		area_of_high_foi = (num * 15.02 * 15.01)

		# 80% of the cumulative average FOI occurs in what percent of the home range area?
		print(id)
		ratios[[id]] = area_of_high_foi / hr_dat[pig == id]$hr_m2
	}
}


### Calculate home range overlap ###


# The first metric of home range overlap we compute is based on the calculation
# given in Appendix S9 such that the metric is directly comparable to 
# the direct + indirect MoveSTIR output.  Specifically, we use the area of
# overlap to directly calculate an FOI between two individuals.

# Get the area of overlap for each pair
collars = sort(unq_collars[unq_collars != 12])
area_of_overlap_m2 = array(NA, dim=c(length(collars), length(collars)))
home_range_mult = array(NA, dim=c(length(collars), length(collars)))

for(i in 1:length(collars)){
	ci = collars[i]
	ai = area(my_hrs[[ci]])

	for(j in 1:length(collars)){
		cj = collars[j]
		aj = area(my_hrs[[cj]])

		int_shp = intersect(my_hrs[[ci]], my_hrs[[cj]])

		# Check if there is any area of overlap
		if(is.null(int_shp)){
			int_area = 0
		} else {
			int_area = area(int_shp)
		}
		area_of_overlap_m2[i, j] = int_area
		home_range_mult[i, j] = ai*aj

	}
}

pathogen_decay = (1 / (24 * 60 * 5)) # 5 days on the minute scale

# You have to scale foi by pathogen decay. See Appendix S9
foi_scaled = data.table((area_of_overlap_m2 / home_range_mult)) / (pathogen_decay)
diag(foi_scaled) = 0 # You can't infect yourself, set diagonal to 0
fwrite(foi_scaled, paste0("../data/home_range_overlaps_", "area", ".csv"))

# Calculate different metrics of overlap.  These are not directly comparable to MoveSTIR
# because their units are different (they are not hazard rates). 
# However, we can still ask questions about what they
# predict regarding how individuals contribute to pathogen spread on the network

home_ranges_trun = home_ranges[unq_collars != 12]
dat[, Name:=individual_ID]
trun_dat = dat[individual_ID != 12][, .(Name, UTMx, UTMy)]
spdat = SpatialPointsDataFrame(trun_dat[, .(UTMx, UTMy)], data=trun_dat[, .(Name)])
crs(spdat) = "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Get utilization distributions on a grid
kud = kernelUD(spdat[, 1], grid=500, same4all=TRUE)

# Methods are defined in the documentation for kerneloverlaphr
methods = c("HR", "BA", "VI")

for(method in methods){
	print(paste("Working on", method))
	method = method
	koverlap = kerneloverlaphr(kud, meth=method, conditional=TRUE)
	diag(koverlap) = 0 # Don't worry about individual's home range overlap with itself
	koverlap_normed = koverlap 
	dt_kover = data.table(koverlap_normed)

	fwrite(dt_kover, paste0("../data/home_range_overlaps_", method, ".csv"))
}



