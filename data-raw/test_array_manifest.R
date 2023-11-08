## code to prepare `test_array_manifest`
temp <- tempfile()
download.file("https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip",temp, mode="wb")
unzip(temp)
fData_EPIC <- read_csv("MethylationEPIC_v-1-0_B4.csv",
                       skip = 7)
array_manifest = fData_EPIC %>%
  dplyr::mutate(STRAND = rep(BiocGenerics::strand("+"), nrow(fData_EPIC))) %>%
  dplyr::select(MAPINFO, CHR, IlmnID, STRAND)

#Get the first 3k probes of the 21 chromosome
test_array_manifest = array_manifest %>%
  filter(CHR == "21") %>%
  arrange(as.numeric(MAPINFO)) %>%
  slice_head(n = 3000) %>%
  select(-IlmnID) #Remove this column because it takes a lot of space when saving the object, and it is already present in the rownames

usethis::use_data(test_array_manifest, overwrite = TRUE)
