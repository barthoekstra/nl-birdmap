library(tidyverse)
library(bioRad)
library(imager)
library(stars)
library(pracma)
library(pbmcapply)
library(sf)
library(raster)
library(patchwork)

visual_filter <- function(pvolfile, overwrite = FALSE) {
  rl <- 160000

  rbc_filename <- paste0("data/rbc/", tools::file_path_sans_ext(basename(pvolfile)), ".RDS")
  rbc_plot_filename <- paste0("data/rbc/", tools::file_path_sans_ext(basename(pvolfile)), ".png")
  if (file.exists(rbc_filename) & file.exists(rbc_plot_filename) & !overwrite) {
    cat("RBC already exists, skipping...\n")
    return("exists")
  }

  # 0. Load pvol
  cat(paste0("Loading pvol: ", pvolfile, "\n"))
  pvol <- read_pvolfile(pvolfile, param = c("DBZH", "DBZV", "RHOHV", "VRADH"))

  if (!"RHOHV" %in% names(pvol$scans[[1]]$params)) {
    cat("No dual-pol products found, moving on\n")
    return("no dual-pol")
  }
  # 0.5 Store PPI from raw scan
  ppi <- project_as_ppi(pvol$scans[[1]], grid_size = 500, range_max = 160000)

  # 1. Deal with EM interference
  ## Identify
  cat("Identify and interpolate EM interference\n")
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
  cat("Remove dual-prf 0.3 degree scan\n")
  lowest_dualprf_scan <- which(sapply(pvol$scans, function(x) {
    x$attributes$how$highprf != 0 & x$attributes$how$lowprf != 0 & round(x$geo$elangle, 1) == 0.3
  }))

  pvol$scans[[lowest_dualprf_scan]] <- NULL

  # 3. Remove azimuth-effect
  cat("Removing azimuth-effect\n")
  pvol <- remove_azimuth_effect(pvol)

  # 4. Classify rain
  ## Calculate DPR
  cat("Calculating DPR\n")
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
  cat("Creating rainmasks\n")

  masks <- lapply(pvol$scans, function(x) {
    if (x$geo$elangle < 90) {
      ppi <- project_as_ppi(get_param(x, "DPR"), grid_size = 500, range_max = rl)

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
  rain_elevations <- apply(rain_elevations, c(1, 2), sum, na.rm = TRUE)

  # 5. Apply RBC
  cat("Apply RBC\n")

  vp_filename <- paste0("data/vp/", tools::file_path_sans_ext(basename(pvolfile)), "_vp.h5")
  if (!file.exists(vp_filename)) {
    cat("Calculating vp\n")
    vp <- calculate_vp(pvolfile, vpfile = vp_filename)
  } else {
    cat("Loading vp\n")
    vp <- read_vpfiles(vp_filename)
  }
  bioRad::sd_vvp_threshold(vp) <- 2

  rbc <- suppressWarnings(integrate_to_ppi(pvol, vp, xlim = c(-rl, rl), ylim = c(-rl, rl), res = 500, param = "DBZH"))

  ## Add rain mask to RBC
  rbc$data$rain <- as.vector(rain_elevations)

  # 6. Save RBC
  cat("Save RBC\n")
  saveRDS(rbc, file = rbc_filename)

  # 7. Save Plot
  cat("Save plot\n")
  (plot(ppi, param = "DBZH") + plot(ppi, param = "VRADH")) / (plot(rbc, param = "VID") + plot(rbc, param = "rain", zlim = c(0, 10))) +
    plot_annotation(title = paste0(toupper(pvol$radar), " ", pvol$datetime)) -> rbc_plot
  ggsave(rbc_plot_filename, plot = rbc_plot, width = 10, height = 10)
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

remove_azimuth_effect <- function(pvol) {
  #@TODO: Consider storing summaries of the effects somewhere

  avg_range <- c(5000, 35000) # min-max meters
  avg_altitude <- c(200, 5000) # min-max meters

  for (i in seq_along(pvol$scans)) {
    if (pvol$scans[[i]]$attributes$where$elangle == 90) next
    avg_start_rangegate_range <- round((avg_range[1] / pvol$scans[[i]]$geo$rscale) + pvol$scans[[i]]$geo$rstart)
    avg_end_rangegate_range <- round((avg_range[2] / pvol$scans[[i]]$geo$rscale) + pvol$scans[[i]]$geo$rstart)

    rangegate_altitudes <- suppressWarnings(
      pvol$geo$height + beam_height(range = seq(from = avg_range[1], to = avg_range[2], by = pvol$scans[[i]]$geo$rscale),
                                    pvol$scans[[i]]$attributes$where$elangle)
    )
    avg_start_rangegate_altitude <- which(rangegate_altitudes > avg_altitude[1])
    if (length(avg_start_rangegate_altitude) > 0) {
      avg_start_rangegate_altitude <- avg_start_rangegate_altitude[[1]]
    } else {
      # Apparently all rangegates within the range are are below the minimum altitude, so skip scan
      next
    }
    avg_end_rangegate_altitude <- which(rangegate_altitudes > avg_altitude[2])
    if (length(avg_end_rangegate_altitude) > 0) {
      avg_end_rangegate_altitude <- avg_end_rangegate_altitude[[1]] - 1
    } else {
      # Apparently all rangegates within the range are below the maximum altitude, so the final rangegate we select is the
      # final rangegate within the range
      avg_end_rangegate_altitude <- length(rangegate_altitudes)
    }

    avg_start_rangegate <- max(avg_start_rangegate_range, avg_start_rangegate_altitude)
    avg_end_rangegate <- min(avg_end_rangegate_range, avg_end_rangegate_altitude)

    DBZH_orig <- pvol$scans[[i]]$params$DBZH

    DBZH <- pvol$scans[[i]]$params$DBZH[avg_start_rangegate:avg_end_rangegate,]
    DBZH[DBZH > 20] <- NA

    # Remove rain
    RHOHV <- pvol$scans[[i]]$params$RHOHV[avg_start_rangegate:avg_end_rangegate,]
    DBZH[RHOHV > 0.90] <- NA

    DBZH_avg <- apply(DBZH, MARGIN = 2, FUN = mean, na.rm = TRUE)

    df <- data.frame("elangle" = rep(pvol$scans[[i]]$attributes$where$elangle, 360),
                     "azimuth" = 1:360,
                     "DBZH_avg" = DBZH_avg,
                     "DBZH_avg_corrected" = NA,
                     "datetime" = rep(pvol$datetime, 360))

    # Sine-fitting
    # Following https://stats.stackexchange.com/a/77865
    fit.lm <- lm(DBZH_avg ~ sin((2 * pi * azimuth) / 180) + cos((2 * pi * azimuth) / 180), data = df)
    b0 <- coef(fit.lm)[1]
    alpha <- coef(fit.lm)[2]
    beta <- coef(fit.lm)[3]

    amplitude <- sqrt(alpha^2 + beta^2)
    phase <- atan2(beta, alpha)

    # Correct DBZH
    azimuths <- 1:360
    azimuths_rad <- (2 * pi * azimuths)/180
    sine <- amplitude * sin(azimuths_rad + phase)
    df$sine <- sine

    DBZH_corrected <- sweep(pvol$scans[[i]]$params$DBZH, 2, sine, "-")

    pvol$scans[[i]]$params$DBZH <- DBZH_corrected
    class(pvol$scans[[i]]$params$DBZH) <- c("param", "matrix", "array")
  }
  return(pvol)
}

cores <- 12
files <- list.files(path = "data/pvol", full.names = TRUE)
files <- files[!files %in% c("data/pvol/dhl_files.sh", "data/pvol/hrw_files.sh")]
remaining_files <- str_replace(files, ".h5", ".RDS")
remaining_files <- str_replace(remaining_files, "/pvol/", "/rbc/")
remaining_files <- files[!file.exists(remaining_files)]
processing <- pbmclapply(remaining_files, visual_filter, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
saveRDS(processing, paste0("data/logs/processing_", format(Sys.time(), "%Y%m%dT%H%M"), ".RDS"))