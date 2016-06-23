cat("Setting data...\n")

# Set paths
genePath <- "~/MGR/ADNI/genetic/links"
brainPath <- "~/MGR/ADNI/imaging/links"
mdtPath <- "~/MGR/ADNI/imaging/ADNI_ICBM9P_MDT.nii"
maskPath <- "~/MGR/ADNI/imaging/T1_biascorr_brain_mask.nii.gz"
infoPath <- "~/MGR/ADNI/ADNIMERGE.csv"
outPath <- "~/MGR/OUT_mock_2mm"
mockPath <- "~/MGR/OUT"

# Stein et al. 2010
force.snps <- c("rs2132683", "rs713155", "rs476463", "rs2429582", "rs9990343", "rs17314229")

system.time(performVGWAS(genePath = genePath, brainPath = brainPath, mdtPath = mdtPath, maskPath = maskPath, infoPath = infoPath, subFactor = 2, force.snps = force.snps, outPath = outPath, sampleSize = 20, mockPath = mockPath))
