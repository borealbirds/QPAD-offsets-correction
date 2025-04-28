#1. SET UP#########

#load libraries
library(raster)
library(sp)
library(suntools)
library(lutz)
library(suncalc)
library(tidyverse)

#get toy data
dat <- read.csv("TSSRExample.csv")

#get the raster
rtz <- raster("utcoffset.tif")

#2. MAKE X CODE#########

#adapted out of the original function

#format to dtm
date <- dat$datetime
dtm <- strptime(date, format="%Y-%m-%d %H:%M:%S", tz="America/Edmonton")

#get day and hour
day <- as.integer(dtm$yday)
hour <- as.numeric(round(dtm$hour + dtm$min/60, 2))

#get coords
lon <- dat$lon
lat <- dat$lat

#project for extraction
xy <- data.frame(x=lon, y=lat)
xy$x[is.na(xy$x)] <- mean(xy$x, na.rm=TRUE)
xy$y[is.na(xy$y)] <- mean(xy$y, na.rm=TRUE)
coordinates(xy) <- ~ x + y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- invisible(spTransform(xy, crs(rtz)))

#get offset
ltz <- raster::extract(rtz, xy) + 7

#get sr
sr <- sunriset(cbind("X"=lon, "Y"=lat),
               as.POSIXct(dtm, tz="America/Edmonton"),
               direction="sunrise", POSIXct.out=FALSE) * 24

#get tssr using original approach
TSSR <- round(unname(hour - sr + ltz), 4)

#get tssr using alternate approach
TSSR2 <- round(unname(hour - sr - ltz), 4)

#put some things on a dataframe
out1 <- dat |>
    cbind(data.frame(TSSR = TSSR,
                     TSSR2 = TSSR2,
                     ltz = ltz,
                     sr = sr))

#3. ALTERNATE APPROACH

#look up timezones with lutz package
out1$tz <- tz_lookup_coords(lat, lon, method="accurate")

#calculate sunrise for each timezone
tzs <- unique(out1$tz)
out2 <- data.frame()
for(i in 1:length(tzs)){

    #get just the data for that timezone
    dat.i <- dplyr::filter(out1, tz==tzs[i]) |>
        mutate(datetime = ymd_hms(datetime, tz=tzs[i]),
               date = as.Date(datetime, tz=tzs[i]))

    #get sunrise
    dat.i$sr2 <- getSunlightTimes(data=dat.i, keep="sunrise", tz=tzs[i])$sunrise

    #output
    out2 <- rbind(out2, dat.i)

}

#calculate tssr
out2$TSSR3 <- as.numeric(difftime(out2$datetime, out2$sr), units="hours")

#4. COMPARE ----

#get differences
out3 <- out2 |>
    mutate(DIFF1 = TSSR - TSSR2,
           DIFF2 = TSSR - TSSR3,
           DIFF3 = TSSR2 - TSSR3)

#Difference between original make_x() with alternate make_x() ----
ggplot(out3) +
    geom_point(aes(x=lon, y=DIFF1))

#Difference between original make_x() with tz loop approach ----
ggplot(out3) +
    geom_point(aes(x=lon, y=DIFF2))

#Difference between original make_x() alternate with tz loop approach ----
ggplot(out3) +
    geom_point(aes(x=lon, y=DIFF3))
