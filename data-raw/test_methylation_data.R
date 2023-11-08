## code to prepare `test_methylation_data` dataset goes here
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

## Create DNAme dataset
sample_size = 30

# Simulate DNAme data
set.seed(123)
distribution_betas = c(rbeta(n = sample_size*nrow(test_array_manifest)/3*2, 5,1), #Methylated distribution - 2 thirds of the distribution
                       rbeta(n = sample_size*nrow(test_array_manifest)/3, 2, 10)) #Unmethylated distribution - 1 third of the dist
m_values = log2(distribution_betas/(1-distribution_betas))

#Make it a data frame
test_methylation_data = matrix(m_values,
                               nrow = nrow(test_array_manifest), ncol = sample_size) %>%
  as.data.frame()

colnames(test_methylation_data) = paste("ID", as.character(1:sample_size), sep = "")
rownames(test_methylation_data) = rownames(test_array_manifest)

usethis::use_data(test_methylation_data, overwrite = TRUE)
