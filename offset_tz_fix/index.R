# --- packages: some are not on CRAN

if (!requireNamespace("QPAD")) {
  if (!requireNamespace("remotes"))
    install.packages("remotes")
  remotes::install_github("psolymos/QPAD")
}
if (!requireNamespace("sp"))
  install.packages("sp")
if (!requireNamespace("maptools"))
  devtools::install_url('https://cran.r-project.org/src/contrib/Archive/maptools/maptools_1.1-8.tar.gz')
if (!requireNamespace("raster"))
  install.packages("raster")
if (!requireNamespace("intrval"))
  install.packages("intrval")

library(QPAD)
library(maptools)
library(intrval)
# library(raster)
library(terra)

# source the functions
source("https://raw.githubusercontent.com/borealbirds/qpad-offsets/refs/heads/main/functions.R")

# update the make_x() function with the TSSR fix
make_x_2 <- function(dt, tm, lon, lat, dur, dis, ..., revert_fix=FALSE, check_xy=FALSE) {
  ## checking lengths
  nn <- c(dt=length(dt), tm=length(tm), lon=length(lon), lat=length(lat), dur=length(dur), dis=length(dis))
  n1 <- nn[nn == 1L]
  n2 <- nn[nn > 1L]
  if (!all(n2 == n2[1L]))
    stop("input lengths must be equal or 1")
  n <- unname(if (length(n2)) n2[1L] else n1[1L])
  if (length(dt) == 1L)
    dt <- rep(dt, n)
  if (length(tm) == 1L)
    tm <- rep(tm, n)
  if (length(lon) == 1L)
    lon <- rep(lon, n)
  if (length(lat) == 1L)
    lat <- rep(lat, n)
  if (length(dur) == 1L)
    dur <- rep(dur, n)
  if (length(dis) == 1L)
    dis <- rep(dis, n)
  ## types
  lat <- as.numeric(lat)
  lon <- as.numeric(lon)
  dur <- as.numeric(dur)
  dis <- as.numeric(dis)
  ## parse date+time into POSIXlt
  dt <- as.character(dt)
  tm <- as.character(tm)
  dtm <- strptime(paste0(dt, " ", tm, ":00"),
    format="%Y-%m-%d %H:%M:%S", tz="America/Edmonton")
  day <- as.integer(dtm$yday)
  hour <- as.numeric(round(dtm$hour + dtm$min/60, 2))
  ## checks
  checkfun <- function(x, name="", range=c(-Inf, Inf)) {
    if (any(x[!is.na(x) & !is.infinite(x)] %)(% range))
      stop(sprintf("Parameter %s is out of range [%.0f, %.0f]", name, range[1], range[2]))
    invisible(NULL)
  }
  ## BCR 4:14 included
  ## crs: WGS84 (EPSG: 4326)
  ## "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #         min       max
  #x -163.89547 -52.66936
  #y   39.66214  68.98741
  if (check_xy) {
    checkfun(lon, "lon", c(-164, -52))
    checkfun(lat, "lat", c(39, 69))
  }
  checkfun(day, "day", c(0, 360))
  checkfun(hour, "hour", c(0, 24))
  checkfun(dur, "dur", c(0, Inf))
  checkfun(dis, "dis", c(0, Inf))
  if (any(is.infinite(lon)))
    stop("Parameter lon must be finite")
  if (any(is.infinite(lat)))
    stop("Parameter lat must be finite")
  ## handling missing values
  ok_xy <- !is.na(lon) & !is.na(lat)

  ## intersect here
  xy <- data.frame(x=lon, y=lat)
  xy$x[is.na(xy$x)] <- mean(xy$x, na.rm=TRUE)
  xy$y[is.na(xy$y)] <- mean(xy$y, na.rm=TRUE)
  xy <- vect(xy, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  xy <- project(xy, crs)

  ## LCC4 and LCC2
  vlcc <- extract(rlcc, xy)$lcc
  # 0: No data (NA/NA)
  # 1: Temperate or sub-polar needleleaf forest (Conif/Forest)
  # 2: Sub-polar taiga needleleaf forest (Conif/Forest)
  # 5: Temperate or sub-polar broadleaf deciduous (DecidMixed/Forest)
  # 6:  Mixed Forest (DecidMixed/Forest)
  # 8: Temperate or sub-polar shrubland (Open/OpenWet)
  # 10: Temperate or sub-polar grassland (Open/OpenWet)
  # 11: Sub-polar or polar shrubland-lichen-moss (Open/OpenWet)
  # 12: Sub-polar or polar grassland-lichen-moss (Open/OpenWet)
  # 13: Sub-polar or polar barren-lichen-moss (Open/OpenWet)
  # 14: Wetland (Wet/OpenWet)
  # 15: Cropland (Open/OpenWet)
  # 16: Barren Lands (Open/OpenWet)
  # 17: Urban and Built-up (Open/OpenWet)
  # 18: Water (NA/NA)
  # 19: Snow and Ice (NA/NA)
  lcclevs <- c("0"="", "1"="Conif", "2"="Conif", "3"="", "4"="",
               "5"="DecidMixed", "6"="DecidMixed", "7"="", "8"="Open", "9"="",
               "10"="Open", "11"="Open", "12"="Open", "13"="Open", "14"="Wet",
               "15"="Open", "16"="Open", "17"="Open", "18"="", "19"="")
  lcc4 <- factor(lcclevs[vlcc+1], c("DecidMixed", "Conif", "Open", "Wet"))
  lcc2 <- lcc4
  levels(lcc2) <- c("Forest", "Forest", "OpenWet", "OpenWet")

  ## TREE
  vtree <- extract(rtree, xy)$tree
  TREE <- vtree / 100
  TREE[TREE %)(% c(0, 1)] <- 0

  ## extract seedgrow value (this is rounded)
  d1 <- extract(rd1, xy)$seedgrow
  ## UTC offset + 7 makes Alberta 0 (MDT offset)
  tz <- raster::extract(rtz, xy)$utcoffset + 7

  ## transform the rest
  JDAY <- round(day / 365, 4) # 0-365
  TREE <- round(vtree / 100, 4)
  MAXDIS <- round(dis / 100, 4)
  MAXDUR <- round(dur, 4)

  ## sunrise time adjusted by offset
  ok_dt <- !is.na(dtm)
  dtm[is.na(dtm)] <- mean(dtm, na.rm=TRUE)
  sr <- suntools::sunriset(cbind("X"=lon, "Y"=lat),
                           as.POSIXct(dtm, tz="America/Edmonton"),
                           direction="sunrise", POSIXct.out=FALSE) * 24

  if (revert_fix) {
    TSSR <- round(unname((hour - (sr - tz)) / 24), 4) # fix here
  } else {
    TSSR <- round(unname((hour - sr - tz) / 24), 4)
  }

  ## days since local spring
  DSLS <- (day - d1) / 365

  out <- data.frame(
    TSSR=TSSR,
    JDAY=JDAY,
    DSLS=DSLS,
    LCC2=lcc2,
    LCC4=lcc4,
    TREE=TREE,
    MAXDUR=MAXDUR,
    MAXDIS=MAXDIS,
    ...)
  out$TSSR[!ok_xy | !ok_dt] <- NA
  out$DSLS[!ok_xy] <- NA
  out$LCC2[!ok_xy] <- NA
  out$LCC4[!ok_xy] <- NA
  out$TREE[!ok_xy] <- NA
  out
}

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

# rlcc <- raster("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/lcc.tif")
# rtree <- raster("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/tree.tif")
# rtz <- raster("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/utcoffset.tif")
# rd1 <- raster("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/seedgrow.tif")
# crs <- proj4string(rtree)
rlcc <- rast("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/lcc.tif")
rtree <- rast("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/tree.tif")
rtz <- rast("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/utcoffset.tif")
rd1 <- rast("https://github.com/borealbirds/qpad-offsets/raw/refs/heads/main/data/seedgrow.tif")
crs <- crs(rtree)

if (F) {
spp <- "ALFL"
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)
dur <- 10 # mins
dis <- 100 # meters

lon <- -113.4938
lat <- 53.5461
x1 <- make_x_2(dt, tm, lon, lat, dur, dis, revert_fix = TRUE)
o1 <- make_off(spp, x1)
x2 <- make_x_2(dt, tm, lon, lat, dur, dis)
o2 <- make_off(spp, x1)
rbind(x1, x2)
rbind(o1, o2)

lon <- -56.26
lat <- 48.58
x1 <- make_x_2(dt, tm, lon, lat, dur, dis, revert_fix = TRUE)
o1 <- make_off(spp, x1)
x2 <- make_x_2(dt, tm, lon, lat, dur, dis)
o2 <- make_off(spp, x1)
rbind(x1, x2)
rbind(o1, o2)

mi <- bestmodelBAMspecies(spp, type="BIC")
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)

}

all_mods <- getBAMmodellist()
all_spp <- getBAMspecieslist()
mods_with_tssr <- all_mods$sra[grepl("TSSR", all_mods$sra)]

spp_with_tssr <- sapply(structure(all_spp, names=all_spp), \(spp) {
    bestmodelBAMspecies(spp, type="BIC")$sra %in% names(mods_with_tssr)
})

spp_to_fix <- names(spp_with_tssr)[spp_with_tssr]

e <- new.env()
load("~/Dropbox/collab/bam-gnm-v4/BAMdb-GNMsubset-2020-01-08.RData", envir=e)
names(e)


dd <- e$dd
tmp <- strsplit(as.character(dd$DATI), " ")
dd$Date <- sapply(tmp, "[[", 1)
dd$Time <- substr(sapply(tmp, "[[", 2), 1, 5)

library(Matrix)
off0 <- as.matrix(e$off)
off0 <- off0[,spp_to_fix]

# see if we can replicate this

X1 <- make_x_2(
    dt = dd$Date, 
    tm = dd$Time, 
    lon = dd$X, 
    lat = dd$Y, 
    dur = dd$MAXDUR, 
    dis = dd$MAXDIS,
    check_xy = FALSE,
    revert_fix = TRUE)
X2 <- make_x_2(
    dt = dd$Date, 
    tm = dd$Time, 
    lon = dd$X, 
    lat = dd$Y, 
    dur = dd$MAXDUR, 
    dis = dd$MAXDIS,
    check_xy = FALSE)
summary(X1$TSSR)
summary(X2$TSSR)

library(ggplot2)

data.frame(TSSR1=X1$TSSR, TSSR2=X2$TSSR) |>
    ggplot(aes(x=TSSR1, y=TSSR2)) +
    geom_bin2d() +
    geom_abline() +
    theme_light()

all(rownames(dd) == rownames(off0))
off1 <- off0
off1[] <- 0
off2 <- off1

for (spp in spp_to_fix) {
    message(spp)
    o1 <- make_off(spp, X1)
    o2 <- make_off(spp, X2)
    off1[,spp] <- o1$offset
    off2[,spp] <- o2$offset
}

save(off0, off1, off2, file="~/Dropbox/collab/bam-gnm-v4/offset_fix_2025-04-03.RData")


# organize results by BCR subunit

detbcr <- e$detbcr

bcr_can <- c("BCR_4", "BCR_5", "BCR_60", 
    "BCR_61", "BCR_70", "BCR_71", "BCR_80", "BCR_81", "BCR_82", "BCR_83", 
    "BCR_9", "BCR_10", "BCR_11", "BCR_12", "BCR_13", "BCR_14")

x <- NULL
for (spp in spp_to_fix) {
    for (bcr in bcr_can) {
        r <- dd[[bcr]] > 0
        z <- data.frame(
            spp = spp, 
            bcr = bcr,
            detbcr = detbcr[spp, bcr],
            mean_off0 = mean(off0[r, spp]),
            mean_off1 = mean(off1[r, spp]),
            mean_off2 = mean(off2[r, spp]))
        x <- rbind(x, z)
    }
}
x$adj <- x$mean_off1 - x$mean_off2
x$alpha <- exp(x$adj)

summary(x[x$detbcr>0,])

boxplot(adj ~ gsub("BCR_", "", bcr), x, range=0, ylim=c(-2,2))
abline(h=0, col=4)

boxplot(alpha ~ gsub("BCR_", "", bcr), x, range=0, ylim=c(0,2))
abline(h=1, col=4)

write.csv(x, row.names=F, file="offset_tz_fix/offset-adjustments-2025-04-04.csv")
