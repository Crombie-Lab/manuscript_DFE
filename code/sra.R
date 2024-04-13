library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# pull fastq file list
files <- data.table::fread("sra/files.txt", header = F) %>%
  dplyr::mutate(read = case_when(grepl(V1, pattern = "_R1_") ~ 1,
                                 grepl(V1, pattern = "_R2_") ~ 2,
                                 TRUE ~ NA_integer_))
r1 <- files %>%
  dplyr::filter(read == 1)

r2 <- files %>%
  dplyr::filter(read == 2)

files2 <- tibble::tibble(filename = r1$V1, filename2 = r2$V1) %>%
  dplyr::mutate(samp = stringr::str_extract(filename, pattern = regex("[^_]*_[^_]*")))

# get the indicies
ind <- readxl::read_excel("sra/20180111_Index_List.xlsx") %>%
  dplyr::mutate(strain = stringr::str_replace(strain, pattern = "\\.", replacement = "_"))

# join the indicies
files3 <- left_join(ind, files2, by = c("strain" = "samp"))

# replace A10_ with RIAIL_ and drop controls
files4 <- files3 %>%
  dplyr::mutate_all(~(stringr::str_replace(., "A10_", "RIAIL_"))) %>%
  dplyr::select(sample_name = strain,
                library = sample_id,
                filename,
                filename2) %>%
  dplyr::filter(sample_name != "CONTROL") %>%
  dplyr::mutate(library_ID = paste(sample_name, library, "paired", "WGS", "IlluminaHiSeq4000", sep = ":"))

# sample + library + strategy + layout + instrument model
# save these
rio::export(files4, file = "sra/proc_metadata.csv")
