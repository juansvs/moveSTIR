library(data.table)
library(adehabitatHR)
library(raster)
library(sp)

## Script calculates the area of the 95% utilization distribution for the pigs used in 
## the MoveSTIR study.
##
## We then use the Average FOI (converted to a raster) generated from MoveSTIR to estimate
## how average FOI on the landscape is distributed within the home range of pigs.

dat = fread("../data/pig_movements.csv")
dat[, datetime:=as.POSIXct(strptime(date_time, format="%Y-%m-%d %H:%M:%S", tz="GMT"))]
dat[, t:=as.numeric(datetime)]
unq_collars = unique(dat$individual_ID)

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



