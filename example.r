cat("Setting data...\n")

# Set paths
genePath <- "~/MGR/ADNI/genetic/links"
brainPath <- "~/MGR/ADNI/imaging/links"
ref.imgPath <- "/usr/share/fsl/data/standard/MNI152_T1_2mm.nii.gz"
maskPath <- "~/MGR/ADNI/imaging/T1_biascorr_brain_mask.nii.gz"
infoPath <- "~/MGR/ADNI/ADNIMERGE.csv"
outPath <- "~/MGR/OUT_2mm_linear_MNI152_2mm"
# Set path to pre-saved flat ROIs and preliminary analysis results
mockPath <- "~/MGR/OUT_2mm_linear_MNI152_2mm"
mockPath.flatROIs <- paste(mockPath, "/flatROIs-", 2, ".R" , sep = "")
mockPath.pre <- paste(mockPath, "/pre-", 2, ".R" , sep = "")
matPath <- "~/MGR/mdt2mni.mat"

# Stein et al. 2010, Michta 2016
force.snps <- c("rs2132683", "rs713155", "rs476463", "rs2429582", "rs9990343", "rs17314229")

# List imaging files
niiFiles <- list.files(brainPath, full.names = TRUE)[]
# List subjectd IDs of available images
niiIDs <- substr(list.files(brainPath), 6, 15)

# Load subject personal data
subject.info <- read.csv(infoPath)

# Select subjects from ADNI1, with their respective images available, and select ID, gender, age, ethnicity and race
subject.info <- subject.info[which(subject.info$COLPROT == "ADNI1" & subject.info$PTID %in% niiIDs), c("PTID", "PTGENDER", "AGE", "PTETHCAT", "PTRACCAT", "DX_bl")]
# Select unique subjects
subject.info <- unique(subject.info)
# Set row names to subject IDs
rownames(subject.info) <- subject.info$PTID

# Select non-hispanic/latino whites
subject.info <- subject.info[which(subject.info$PTETHCAT=="Not Hisp/Latino" & subject.info$PTRACCAT=="White"),]
# Set gender and age as covariates
covar <- subject.info[,c("PTGENDER", "AGE")]
# Recode gender levels
levels(covar$PTGENDER)[levels(covar$PTGENDER)=="Male"] <- 1
levels(covar$PTGENDER)[levels(covar$PTGENDER)=="Female"] <- 2
# Transpose
covar <- t(covar)
# Sort covariates by subject ID
covar <- covar[, order(colnames(covar))]

time <- system.time(res <- performVGWAS(genePath = genePath, niiFiles = niiFiles, niiIDs = niiIDs, covar = covar, ref.imgPath = ref.imgPath,
                         maskPath = maskPath, subFactor = 1, out.subFactor = 0, matPath = matPath,
                         force.snps = force.snps, outPath = outPath, useModel = MatrixEQTL::modelLINEAR, mockPath.flatROIs = mockPath.flatROIs, mockPath.pre = mockPath.pre))
cat(sprintf("Total time: %.2fs\n", time[3]))
