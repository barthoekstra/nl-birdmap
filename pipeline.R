library(tidyverse)
library(bioRad)
library(imager)
library(stars)
library(pracma)
# library(ncdf4)
library(pbmcapply)
library(sf)
library(raster)
# library(rdrop2)

visual_filter <- function(pvolfile, vp) {
  rl <- 160000  # Was 160000
  # 0. Load pvol
  print("Loading pvol")
  pvol <- read_pvolfile(pvolfile, param = c("DBZH", "DBZV", "RHOHV"))

  # 1. Deal with EM interference
  ## Identify
  print("Identify and interpolate EM interference")
  eminterference <- identify_em_interference(pvol)

  ## Interpolate beams where it was found
  pvol$scans <- mapply(function(scan, beams) {
    if (length(beams) > 0) {
      interpolate_em(scan, beams)
    } else {
      scan
    }
  }, pvol$scans, eminterference, SIMPLIFY = FALSE)

  # 2. Remove dual-prf 0.3 degree scan
  print("Remove dual-prf 0.3 degree scan")
  lowest_dualprf_scan <- which(sapply(pvol$scans, function(x) {
    x$attributes$how$highprf != 0 & x$attributes$how$lowprf != 0 & round(x$geo$elangle, 1) == 0.3
  }))

  pvol$scans[[lowest_dualprf_scan]] <- NULL

  # 3. Classify rain
  ## Calculate DPR
  print("Calculating DPR")
  pvol <- suppressWarnings(
    calculate_param(pvol,
                    ZDRL = 10 ** ((DBZH - DBZV) /10),
                    DPR = 10 * log10((ZDRL + 1 - 2 * ZDRL^0.5 * RHOHV) / (ZDRL + 1 + 2 * ZDRL^ 0.5 * RHOHV)))
  )

  ## Remove pixels for which no DPR can be calculated
  pvol$scans <- lapply(pvol$scans, function(x) {
    x$params[["DBZH"]][is.na(x$params[["DPR"]])] <- NA
    return(x)
  })

  ## Create rainmasks
  print("Creating rainmasks")

  # vp <- calculate_vp(pvolfile)
  lon <- vp$attributes$where$lon
  lat <- vp$attributes$where$lat
  # e <- ((s <- st_point(c(lon, lat)) %>% st_sfc() %>% st_set_crs(4326) %>% st_transform(3857) %>% as_Spatial()) %>% extent()) + 160000 * 2
  # ras <- raster::raster(e, resolution = 500, crs = CRS(proj4string(s)))
  # proj4 <- CRS(paste("+proj=aeqd +lat_0=", lat, " +lon_0=", lon, " +units=m", sep = ""))
  # e_rain <- projectExtent(ras, proj4)
  # ras_rain <- raster::raster(e_rain)
  # proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # newext <- extent(projectExtent(ras,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # ras_rain <- raster::raster(newext, resolution = 500, crs = CRS(proj4))

  masks <- lapply(pvol$scans, function(x) {
    if (x$geo$elangle < 90) {
      ppi <- project_as_ppi(get_param(x, "DPR"), grid_size = 500, range_max = rl)
      # ppi <- project_as_ppi(get_param(x, "DPR"), grid_size = (rl*2) / 640, range_max = rl)
      # ppi <- project_as_ppi(get_param(x, "DPR"), grid_size = ras)
      # ppi <- project_as_ppi(get_param(x, "DPR"), grid_size = 500, range_max = round(mean(abs(as.matrix(ras_rain@extent))), -3))
      dpr <- as.cimg(as.matrix(ppi$data))
      (dpr <= -12 & !is.na(dpr)) %>%
        clean(4) %>%  # Remove speckles by shrinking then growing using a 4px radius
        fill(10) %>%  # Fill small gaps within contiguous areas
        identity() -> filled
      if (sum(filled) > 0) {
        filled %>%
          split_connected() %>%  # Split image in contiguous areas classified as rain
          purrr::keep(~ sum(.) > 100) %>%  # Only keep contiguous rain areas if area is > 100 pixels
          parany() -> contiguous  # Merge to 1 image again
      }

      if (exists("contiguous") && !is.null(contiguous)) {  # Only buffer if any rain areas of > 50 pixels are retained
        contiguous %>%
          distance_transform(1, 2) %>%  # Calculate Euclidean distance (2nd argument) to pixels classified as 1
          threshold(5) -> dpr_mask

        dpr_mask %>%
          distance_transform(1, 2) -> dist_mask

        dpr_mask <- -dpr_mask

        dpr_mask <- resize(dpr_mask, size_x = 640, size_y = 640)
        dist_mask <- resize(dist_mask, size_x = 640, size_y = 640)

        m <- as.matrix(dpr_mask)
        m[m == 0] <- 1
        m[m == -1] <- NA

        d <- as.matrix(dist_mask)
        d[d == 0] <- NA
        return(list(m, d))
      } else {
        o <- matrix(nrow = 640, ncol = 640)
        return(list(o, o))
      }
    } else {
      o <- matrix(nrow = 640, ncol = 640)
      return(list(o, o))
    }
  })
  rain_elevations <- simplify2array(lapply(masks, function(x) x[[1]]))
  # rain_elevations_low <- rain_elevations
  rain_elevations <- apply(rain_elevations, c(1, 2), sum, na.rm = TRUE)
  # low_elevations <- get_elevation_angles(pvol) < 5
  # rain_elevations_low <- apply(rain_elevations_low[, , low_elevations], c(1, 2), sum, na.rm = TRUE)
  # rain_elevations[rain_elevations <= 1] <- NA

  # 4. Apply RBC
  print("Apply RBC")

  bioRad::sd_vvp_threshold(vp) <- 2
  # rbc <- integrate_to_ppi(pvol, vp, xlim = c(-160000, 160000), ylim = c(-160000, 160000), res = 500, param = "DBZH")
  # rbc <- integrate_to_ppi(pvol, vp, xlim = c(-rl, rl), ylim = c(-rl, rl), res = 500, param = "DBZH")
  # rbc <- integrate_to_ppi(pvol, vp, xlim = c(-rl, rl), ylim = c(-rl, rl), res = (rl*2) / 640, param = "DBZH")
  rbc <- integrate_to_ppi(pvol, vp, xlim = c(-rl, rl), ylim = c(-rl, rl), res = 500, param = "DBZH")
  # pvol_low <- pvol
  # pvol_low$scans <- pvol_low$scans[low_elevations]
  # rbc_low <- integrate_to_ppi(pvol_low, vp, raster = ras, param = "DBZH")

  ## Add rain mask to RBC
  rbc$data$rain <- as.vector(rain_elevations)
  # rbc$data$rain_low <- as.vector(rain_elevations_low)
  # rbc$data$VIR_low <- as.vector(rbc_low$data$VIR)

  # 5. Save RBC
  print("Save RBC")
  # saveRDS(rbc, paste0("data/rbc/", tools::file_path_sans_ext(basename(pvolfile)), ".png"))
  rbc
  # t <- do.call("c", lapply(pvol$scans, function(x) {
  #   s <- strptime(paste(x$attributes$what[c("startdate", "starttime")], collapse = " "), "%Y%m%d %H%M%S", tz = "UTC")
  #   e <- strptime(paste(x$attributes$what[c("enddate", "endtime")], collapse = " "), "%Y%m%d %H%M%S", tz = "UTC")
  #   return(s + (e - s) / 2)
  # }))
  #
  # rbc_dset <- stars::st_as_stars(rbc$data)
  # attr(rbc_dset, "bbox") <- rbc$geo$bbox
  # attr(rbc_dset, "radar_lat") <- pvol$geo$lat
  # attr(rbc_dset, "radar_lon") <- pvol$geo$lon
  # attr(rbc_dset, "radar") <- pvol$radar
  # attr(rbc_dset, "time") <- pvol$datetime
  #
  # saveRDS(rbc, file = paste0("data/RBC/RDS/", str_replace(rbc_file_name(rbc_dset), ".nc", ".RDS")))
  # save_rbc_nc(rbc_dset)
  #
  # if (pvol$radar == "NL52") {
  #   drop_upload(paste0("data/RBC/", rbc_file_name(rbc_dset)),
  #               path = "UvA/ClearSkies/NLHRW20181019/", mode = "add", verbose = FALSE, dtoken = token)
  #   print(paste0("Uploaded NLHRW data to Dropbox: ", rbc_file_name(rbc_dset)))
  # }
  # if (pvol$radar == "NL51") {
  #   drop_upload(paste0("data/RBC/", rbc_file_name(rbc_dset)),
  #               path = "UvA/ClearSkies/NLDHL20181019/", mode = "add", verbose = FALSE, dtoken = token)
  #   print(paste0("Uploaded NLDHL data to Dropbox: ", rbc_file_name(rbc_dset)))
  # }
}

range_coverage <- function(scan) {
  s <- as.matrix(scan$params$DBZH)
  class(s) <- "matrix"

  hasvalue <- s
  hasvalue[!is.na(s)] <- 1
  coverage <- colSums(hasvalue, na.rm = TRUE)
  coverage / dim(s)[1]
}

calc_linearity <- function(scan) {
  s <- as.matrix(scan$params$DBZH)
  class(s) <- "matrix"

  apply(s, 2, function(d) {
    nearestbin <- round(50000 / scan$geo$rscale)
    r <- nearestbin:dim(s)[1]
    dbzh <- d[r]
    data <- data.frame(r = r, dbzh = dbzh) %>% drop_na()
    if (nrow(data) > 20) {
      m <- lm(dbzh ~ r, data = data)
      if (coef(m)[2] > 0) {  # Only return non-NA if slope is positive
        summary(m)$r.squared
      } else {
        NA
      }
    } else {
      NA
    }
  })
}

interpolate_em <- function(scan, beams) {
  s <- as.matrix(scan$params$DBZH)
  class(s) <- "matrix"

  consecutive_beams <- split(beams, cumsum(c(1, diff(beams) != 1)))

  for (cb in consecutive_beams) {
    extract_beams <- c(min(cb) - 1, cb, max(cb) + 1)
    if (any(extract_beams < 1)) {
      extract_beams <- ((extract_beams - 1) %% 360) + 1
    }
    if (any(extract_beams > 360)) {
      extract_beams[which(extract_beams > 360)] <- extract_beams[which(extract_beams > 360)] - 360
    }
    m <- s[, extract_beams]
    m[is.na(m)] <- -9999
    m[, 2:(length(extract_beams) - 1)] <- NA
    x <- 1:dim(m)[1]  # Ranges
    y <- c(1, length(extract_beams))  # Azimuths
    z <- t(m[x, y])
    xp <- x
    yp <- 2:(length(extract_beams) - 1)
    ip <- expand.grid(x, yp)
    mi <- matrix(interp2(x, y, z, ip[, 1], ip[, 2], method = "nearest"), nrow = length(x))
    mi[mi == -9999] <- NA
    # Now that we have interpolated NA values using nearest-neighbor, we can interpolate reflectivity
    ip2 <- which(!is.na(mi), arr.ind = TRUE)
    mp <- interp2(x, y, z, ip2[, 1], ip2[, 2], method = "linear")
    mi[cbind(ip2[, 1], ip2[, 2])] <- mp
    s[, cb] <- mi
  }

  scan$params$DBZH <- s
  attr(scan$params$DBZH, "class") <- c("param", "matrix", "array")
  return(scan)
}

identify_em_interference <- function(pvol) {
  beams <- lapply(pvol$scans, function(x) {
    classified_em <- which(calc_linearity(x) > 0.75 & range_coverage(x) > 0.75)
    classified_em <- sort(unique(c(classified_em - c(1), classified_em, classified_em + c(1))))
    classified_em
  })

  elevs <- round(get_elevation_angles(pvol), 1)

  beams <- lapply(elevs, function(x) {
    identical <- which(elevs == x)
    if (length(identical) > 0) {
      b <- unique(unlist(beams[identical]))
    } else {
      NULL
    }
  })
  beams
}

# save_rbc_nc <- function(rbc_dset) {
#   time_stamp_raw <- attr(rbc_dset, "time")
#   tz_str <- attributes(time_stamp_raw)$tz
#   time_str <- paste(as.character(time_stamp_raw),tz_str,sep=' ')
#
#   fpath_vir <- paste0("data/RBC/", tools::file_path_sans_ext(rbc_file_name(rbc_dset)), ".nc")
#   # fpath_vir_low <- paste0("data/RBC/", tools::file_path_sans_ext(rbc_file_name(rbc_dset)), "_low.nc")
#   fpath_rain <- paste0("data/RBC/", tools::file_path_sans_ext(rbc_file_name(rbc_dset)), "_rain.nc")
#   # fpath_rain_low <- paste0("data/RBC/", tools::file_path_sans_ext(rbc_file_name(rbc_dset)), "_rain_low.nc")
#
#   write_stars(rbc_dset, dsn = fpath_vir, layer = "VIR")
#   # write_stars(rbc_dset, dsn = fpath_vir_low, layer = "VIR_low")
#   write_stars(rbc_dset, dsn = fpath_rain, layer = "rain")
#   # write_stars(rbc_dset, dsn = fpath_rain_low, layer = "rain_low")
#
#   nv_vir <- ncdf4::nc_open(filename = fpath_vir, write = TRUE)
#   # nv_vir_low <- ncdf4::nc_open(filename = fpath_vir_low, write = TRUE)
#   nv_rain <- ncdf4::nc_open(filename = fpath_rain, write = TRUE)
#   # nv_rain_low <- ncdf4::nc_open(filename = fpath_rain_low, write = TRUE)
#   nv_vir <- ncvar_rename(nv_vir, "Band1", "VIR")
#   # nv_vir_low <- ncvar_rename(nv_vir_low, "Band1", "VIR_low")
#   nv_rain <- ncvar_rename(nv_rain, "Band1", "rain")
#   # nv_rain_low <- ncvar_rename(nv_rain_low, "Band1", "rain_low")
#
#   # Write VIR_low to nv_vir
#   # xdim <- nv_vir_low$dim[["x"]]
#   # ydim <- nv_vir_low$dim[["y"]]
#   # mv <- NA
#   # VIR_low <- ncvar_def("VIR_low", "VIR", list(xdim, ydim), mv)
#   # nv_vir <- ncvar_add(nv_vir, VIR_low)
#   # ncvar_put(nv_vir, "VIR_low", ncvar_get(nv_vir_low, "VIR_low"))
#
#   # Write rain to nv_vir
#   xdim <- nv_rain$dim[["x"]]
#   ydim <- nv_rain$dim[["y"]]
#   mv <- NA
#   rain <- ncvar_def("rain", "nr", list(xdim, ydim), mv)
#   nv_vir <- ncvar_add(nv_vir, rain)
#   ncvar_put(nv_vir, "rain", ncvar_get(nv_rain, "rain"))
#
#   # Write rain low to nv_vir
#   # xdim <- nv_rain_low$dim[["x"]]
#   # ydim <- nv_rain$dim[["y"]]
#   # mv <- NA
#   # rain_low <- ncvar_def("rain_low", "nr", list(xdim, ydim), mv)
#   # nv_vir <- ncvar_add(nv_vir, rain_low)
#   # ncvar_put(nv_vir, "rain_low", ncvar_get(nv_rain_low, "rain_low"))
#
#   # Time
#   ncatt_put(nc = nv_vir, varid = 0, attname = "time", attval = time_str,
#             prec = typeof(time_str), definemode = FALSE)
#   # Radar Lat
#   ncatt_put(nc = nv_vir, varid =  0, attname = 'radar_lat',
#             attval = attributes(rbc_dset)$radar_lat,
#             prec = typeof(attributes(rbc_dset)$radar_lat), definemode = FALSE)
#   # Radar Lon
#   ncatt_put(nc = nv_vir, varid =  0, attname = 'radar_lon',
#             attval = attributes(rbc_dset)$radar_lon,
#             prec = typeof(attributes(rbc_dset)$radar_lon), definemode = FALSE)
#   # BioRad version
#   ncatt_put(nc = nv_vir, varid =  0, attname = 'BioRadVersion',
#             attval = getNamespaceVersion('bioRad')[[1]],
#             prec = typeof(getNamespaceVersion('bioRad')[[1]]),
#             definemode = FALSE)
#   # Bounding box
#   ncatt_put(nc = nv_vir, varid = 0, attname = 'bbox',
#             attval = attributes(rbc_dset)$bbox,
#             prec = typeof(attributes(rbc_dset)$bbox), definemode = FALSE)
#
#   nc_close(nc = nv_vir)
#   # nc_close(nc = nv_vir_low)
#   nc_close(nc = nv_rain)
#   # nc_close(nc = nv_rain_low)
#   file.remove(fpath_rain)
#   # file.remove(fpath_rain_low)
#   # file.remove(fpath_vir_low)
#   message("RBC succesfully saved to disk")
# }
#
# rbc_file_name <- function(rbc_dset){
#   # Function Naming convention
#   # In RBC DSET
#   # OUT fname
#
#   time_raw <- attr(rbc_dset, "time")
#   time_str <- format(x = time_raw,format = '%Y%m%dT%H%M' ,tz = 'UTC') # YYYmmddTHHMM
#   radar <- attr(rbc_dset, "radar")
#   if (radar == "NL52") {
#     radar_str <- "NLHRW"
#   } else if (radar == "NL51") {
#     radar_str <- "NLDHL"
#   } else {
#     radar_str <- "NULL"
#   }
#   dtype = 'RBC' # Is not known from the format, adding it explicitly now.
#   fname = paste0(paste(radar_str,dtype,time_str, sep = '_'),".nc")
#   return(fname)
# }

# vp_file_name <- function(pvolfile) {
#   filename <- strsplit(tools::file_path_sans_ext(basename(pvolfile)), "_")[[1]]
#   fname <- paste0(dirname(pvolfile), "/vp/", paste(filename[1], "vp", filename[3], filename[4], "v0-3-20.h5", sep = "_"))
#   return(fname)
# }

# token <- readRDS("token.RDS")

# cores <- 7
# # files <- list.files(path = "data/20201001", full.names = TRUE)
# # processing1 <- pbmclapply(files, visual_filter, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
# # saveRDS(processing1, file = paste0("data/20201001_processing_3.RDS"))
#
# # files2 <- list.files(path = "data/20201002", full.names = TRUE)
# files <- list.files(path = "data/20181019/", full.names = TRUE)
# processing2 <- pbmclapply(files, visual_filter, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
# saveRDS(processing2, file = paste0("data/20201002_processing_5.RDS"))

# reprocess <- sapply(processing1, function(x) !is.null(x))

# visual_filter(files[1])
# visual_filter(files[410])
# visual_filter("data/20201001/NLHRW_pvol_20201001T1740_6356.h5")

# visual_filter("data/20201002/NLHRW_pvol_20201002T1205_6356.h5")

# files_rbc <- list.files(path = "data/RBC", full.names = TRUE)
# files_dropbox <- paste0("UvA/ClearSkies/Data/", basename(files_rbc))


# mapply(function(x, y) drop_upload(x, path = "UvA/ClearSkies/Data_bioRad_Defaults/", mode = "add", verbose = FALSE, dtoken = token), files_rbc, files_dropbox)
# lapply(files_rbc, function(x) drop_upload(x, path = "UvA/ClearSkies/Data_bioRad_Defaults", mode = "add", verbose = FALSE, dtoken = token))