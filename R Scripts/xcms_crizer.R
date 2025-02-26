library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(pheatmap)
library(MsExperiment)
library(BiocManager)
library(MetaboAnnotation)
library(Spectra)
library(msdata)
ms1_features <- read.table(system.file("extdata", "MS1_example.txt",
                                       package = "MetaboAnnotation"),
                           header = TRUE, sep = "\t")
target_df <- read.table(system.file("extdata", "LipidMaps_CompDB.txt",
                                    package = "MetaboAnnotation"),
                        header = TRUE, sep = "\t")
# parm <- Mass2MzParam(adducts = MetaboCoreUtils::adductNames(polarity = "positive"),
#                      tolerance = 0.005, ppm = 0)
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
                       tolerance = 0.005, ppm = 0)

matched_features <- matchValues(ms1_features[1:100, ], target_df, parm)
matched_features
query(matched_features)
target(matched_features)
whichQuery(matched_features)
whichTarget(matched_features)
colnames(matched_features)
matchedData(matched_features)
ms1_subset <- ms1_features[1:100, ]
set.seed(123)
target_df$rtime <- sample(ms1_subset$rtime,
                          nrow(target_df), replace = TRUE) + 2
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+Na]+"),
                       tolerance = 0.005, ppm = 0,
                       toleranceRt = 10)
matched_features <- matchValues(ms1_subset, target_df, param = parm,
                                rtColname = "rtime")
matched_features
matchedData(matched_features)[whichQuery(matched_features), ]
se <- SummarizedExperiment(
  assays = matrix(rnorm(nrow(ms1_features) * 4), ncol = 4,
                  dimnames = list(NULL, c("A", "B", "C", "D"))),
  rowData = ms1_features)
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+Na]+"),
                     tolerance = 0.005, ppm = 0)
matched_features <- matchValues(se, target_df, param = parm)
matched_features
colnames(matched_features)
matchedData(matched_features)
matched_sub <- matched_features[whichQuery(matched_features)]
MetaboAnnotation::query(matched_sub)
library(QFeatures)
qf <- QFeatures(list(features = se))
qf
matched_qf <- matchValues(qf, target_df, param = parm, queryAssay = "features")
matched_qf
colnames(matched_qf)
matchedData(matched_qf)

fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
fl
pest_ms2 <- filterMsLevel(Spectra(fl), 2L)
## subset to selected spectra.
pest_ms2 <- pest_ms2[c(808, 809, 945:955)]
## assign arbitrary *feature IDs* to each spectrum.
pest_ms2$feature_id <- c("FT001", "FT001", "FT002", "FT003", "FT003", "FT003",
                         "FT004", "FT004", "FT004", "FT005", "FT005", "FT006",
                         "FT006")
## assign also *spectra IDs* to each
pest_ms2$spectrum_id <- paste0("sp_", seq_along(pest_ms2))
pest_ms2
load(system.file("extdata", "minimb.RData", package = "MetaboAnnotation"))
minimb
csp <- CompareSpectraParam(requirePrecursor = TRUE, ppm = 10)
mtches <- matchSpectra(pest_ms2, minimb, param = csp)
mtches
mtches[1]
mtches[2]
spectraVariables(mtches)
mtches[2]$target_compound_name
mtches$spectrum_id
mtches_df <- spectraData(mtches, columns = c("spectrum_id", "feature_id",
                                             "score", "target_spectrum_id",
                                             "target_compound_name"))
as.data.frame(mtches_df)
plotSpectraMirror(mtches[2])
scale_int <- function(x, ...) {
  x[, "intensity"] <- x[, "intensity"] / max(x[, "intensity"], na.rm = TRUE)
  x
}
mtches <- addProcessing(mtches, scale_int)
plotSpectraMirror(mtches[2])
mp <- MatchForwardReverseParam(requirePrecursor = TRUE, ppm = 10)
mtches <- matchSpectra(pest_ms2, minimb, param = mp)
mtches
as.data.frame(
  spectraData(mtches, c("spectrum_id", "target_spectrum_id",
                        "target_compound_name", "score", "reverse_score",
                        "presence_ratio")))
select_top_match <- function(x) {
  which.max(x)
}
csp2 <- CompareSpectraParam(ppm = 10, requirePrecursor = FALSE,
                            THRESHFUN = select_top_match)
mtches <- matchSpectra(pest_ms2, minimb, param = csp2)
res <- spectraData(mtches, columns = c("spectrum_id", "target_spectrum_id",
                                       "target_compound_name", "score"))
as.data.frame(res)
res
mtches
plotSpectraMirror(mtches[2],ylab="Normalized Intensity")
title("Azaconazole MS/MS Top Match")
mbank <- MassBankSource("2022.06")
mbank
res <- matchSpectra(
  pest_ms2, mbank,
  param = CompareSpectraParam(requirePrecursor = TRUE, ppm = 10))

res
target(res)
matchedData(res)$target_name
res <- addProcessing(res, scale_int)
plotSpectraMirror(res[2], ylab = "Normalized Intensity")   
colnames(res[2])
?CompareSpectraParam

