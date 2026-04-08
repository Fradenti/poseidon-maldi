rescale_mz <- function(mz, min, max){
  mz01 <- (mz - min(round(mz)))/(max(round(mz))-min(round(mz)))
  round((max-min)*mz01 + min,4)
}

simulateImage_upd <- function (pixelData, featureData, preset, from = 0.9 * min(mz),
    to = 1.1 * max(mz), by = 400, sdrun = 1, sdpixel = 1, spcorr = 0.3,
    units = c("ppm", "mz"), representation = c("profile", "centroid"),
    nchunks = getCardinalNChunks(), verbose = getCardinalVerbose(),
    BPPARAM = getCardinalBPPARAM(), ...)
{
    # if (!missing(preset) || !is.null(preset)) {
    if (!missing(preset) && !is.null(preset)) {
        preset <- presetImageDef(preset, ...)
        featureData <- preset$featureData
        pixelData <- preset$pixelData
    }
    nms <- intersect(names(pixelData), names(featureData))
    if (length(nms) == 0L)
        stop("no shared column names between pixelData and featureData")
    pData <- pixelData[nms]
    fData <- featureData[nms]
    mz <- mz(featureData)
    if ((!missing(from) || !missing(to)) && !missing(preset)) {
        mz <- (mz - min(mz))/max(mz - min(mz))
        mz <- (from + 0.1 * (to - from)) + (0.8 * (to - from)) *
            mz
        mz(featureData) <- mz
    }
    units <- match.arg(units)
    representation <- match.arg(representation)
    domain <- mz(from = from, to = to, by = by, units = units)
    # seeds <- RNGStreams(nrun(pixelData))
    seeds <- list(trunc(1e+8 * runif(nrun(pixelData))))
    FUN <- function(irun) {
        group <- as.matrix(pData[irun == run(pixelData), , drop = FALSE])
        intensity <- as.matrix(fData)
        runerr <- rnorm(nrow(fData), sd = sdrun)
        pixelerr <- rnorm(nrow(group), sd = sdpixel)
        if (spcorr > 0) {
            pDataRun <- pixelData[irun == run(pixelData), , drop = FALSE]
            W <- as.matrix(spatialWeights(pDataRun, r = 1, matrix = TRUE))
            IrW <- as(diag(nrow(W)) - spcorr * W, "sparseMatrix")
            SARcov <- as(Matrix::solve(t(IrW) %*% IrW), "denseMatrix")
            SARcovL <- Matrix::chol((SARcov + t(SARcov))/2)
            pixelerr <- as.numeric(t(SARcovL) %*% pixelerr)
            pixelerr <- sdpixel * ((pixelerr - mean(pixelerr))/sd(pixelerr))
        }
        spectra <- lapply(seq_len(nrow(group)), function(i) {
            j <- group[i, ]
            if (any(j)) {
                x <- rowSums(intensity[, j, drop = FALSE])
            }
            else {
                x <- rep(0, nrow(fData))
            }
            nz <- x != 0
            x[nz] <- pmax(0, x[nz] + runerr[nz] + pixelerr[i])
            simulateSpectra(1L, mz = mz, intensity = x, from = from,
                to = to, by = by, units = units,
                representation = representation,
                ...)$intensity
        })
        spectra <- do.call(cbind, spectra)
        if (representation == "profile") {
            MSImagingExperiment(spectra, featureData = MassDataFrame(mz = domain),
                pixelData = pixelData[irun == run(pixelData),
                  , drop = FALSE], centroided = FALSE)
        }
        else {
            MSImagingExperiment(spectra, featureData = MassDataFrame(mz = mz),
                pixelData = pixelData[irun == run(pixelData),
                  , drop = FALSE], centroided = TRUE)
        }
    }
    ans <- matter::chunkMapply(FUN, runNames(pixelData), nchunks = nchunks,
        verbose = verbose, seeds = seeds, BPPARAM = BPPARAM)
    ans <- do.call("cbind", ans)
    if (representation == "centroid") {
        nms <- setdiff(names(featureData), names(featureData(ans)))
        featureData(ans)[nms] <- featureData[nms]
        centroided(ans) <- TRUE
    }
    design <- list(pixelData = pixelData, featureData = featureData)
    metadata(ans)[["design"]] <- design
    ans
}


get_peaks_mean <- function(out){
  design <- out@metadata$design
  pixelData <- design$pixelData
  featureData <- design$featureData
  nms <- intersect(names(pixelData), names(featureData))
  pData <- pixelData[nms]
  fData <- featureData[nms]
  mz <- mz(featureData)
  irun = runNames(pixelData)
  group <- as.matrix(pData[irun == run(pixelData), , drop = FALSE])
  intensity <- as.matrix(fData)

  peaks_mean_values <- lapply(seq_len(nrow(group)), function(i) {
    j <- group[i, ]
    if (any(j)) {
      x <- rowSums(intensity[, j, drop = FALSE])
    }
    else {
      x <- rep(0, nrow(fData))
    }
    x
  })
  peaks_mean_values <- do.call(cbind, peaks_mean_values)
  peaks_mean_values
}

get_clustering <- function(out){
  design <- out@metadata$design

  n_pix <- nrow(design$pixelData)
  unik_val_pix <- unique(design$pixelData[,4:6])
  K_pix <- nrow(unik_val_pix)
  cl_pix <- numeric(n_pix)
  for(k in 1:K_pix){
    cl_pix[ which(design$pixelData[,4:6] == unik_val_pix[k,]) ] <- k
  }
  if(K_pix != 6){
    cl_pix2 = as.factor(cl_pix)
  } else {
    labels <- c("sq1","no","sq1+c","c","sq2+c","sq2")
    cl_pix2 <-  labels[cl_pix]
    cl_pix2 <- factor(cl_pix2, levels = labels[c(2,1,3,4,5,6)])
  }
  return(cl_pix2)
}




