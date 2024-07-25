# testing ld
load("data/processed/ld.data.for.heatmaps.rda")

# locus 1
# row 10, col 1
# I_113752 asnd I_7660120
D <- ld.data.for.heatmaps$ld.all.D[10,1]
r <- ld.data.for.heatmaps$ld.all.r[10,1]
r2 <- ld.data.for.heatmaps$ld.all.r2[10,1]

# get genotype matrix
gm <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>%
  dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
  dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
  tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n.na = sum(is.na(genotype))) %>%
  dplyr::filter(n.na<146) %>% # is this the right na filter?
  dplyr::ungroup() %>%
  dplyr::select(-n.na) %>%
  dplyr::mutate(genotype = case_when(genotype == MA563 ~ "L/L",
                                      genotype == MA530 ~ "H/H",
                                      TRUE ~ NA_character_)) %>%
  dplyr::select(-MA530, -MA563) %>%
  tidyr::pivot_wider(names_from = "strain", values_from = "genotype")

go_list <- NULL
for (i in 1:nrow(gm)) {
  go_list[[i]] <- genetics::genotype(as.character(gm[i, 2:ncol(gm)]))
}

# make a dataframe from list
go_df <- data.frame(go_list)

# add informative column names to dataframe
colnames(go_df) <- (gm$snp_id)

# I_113752 asnd I_7660120

test <- go_df %>%
  dplyr::select(I_113752, I_7660120) %>%
  dplyr::filter(!is.na(I_113752), !is.na(I_7660120)) %>%
  dplyr::mutate(pA = sum(I_113752 == "H/H", na.rm = T)/n(),
                pa = sum(I_113752 == "L/L", na.rm = T)/n(),
                pB = sum(I_7660120 == "H/H", na.rm = T)/n(),
                pb = sum(I_7660120 == "L/L", na.rm = T)/n())

pA <- unique(test$pA)
pa <- unique(test$pa)
pB <- unique(test$pB)
pb <- unique(test$pb)
  

-D / sqrt(pA * pa * pB * pb)

my.r <- (-D / sqrt(pA * pa * pB * pb))
my.r
r
r2
my.r2 <- my.r^2


#======================= RILS only
gm_RILs <- gm %>%
   tidyr::pivot_longer(col = 2:ncol(gm), names_to = "strain") %>%
   dplyr::filter(stringr::str_detect(strain, pattern = "RIL")) %>%
  tidyr::pivot_wider(names_from = strain)

# 1) RILs
# make a list of genetics genotype objects from genotype matrix
go_list_RILs <- list()
for (i in 1:nrow(gm_RILs)) {
  go_list_RILs[[i]] <- genetics::genotype(as.character(gm[i, 2:ncol(gm_RILs)]))
}

# make a dataframe from list
go_df_RILs <- data.frame(go_list_RILs)

# add informative column names to dataframe
colnames(go_df_RILs) <- (gm_RILs$snp_id)


test.ril <- go_df_RILs %>%
  dplyr::select(I_113752, I_7660120) %>%
  dplyr::filter(!is.na(I_113752), !is.na(I_7660120)) %>%
  dplyr::mutate(pA = sum(I_113752 == "H/H", na.rm = T)/n(),
                pa = sum(I_113752 == "L/L", na.rm = T)/n(),
                pB = sum(I_7660120 == "H/H", na.rm = T)/n(),
                pb = sum(I_7660120 == "L/L", na.rm = T)/n())

pA.ril <- unique(test.ril$pA)
pa.ril <- unique(test.ril$pa)
pB.ril <- unique(test.ril$pB)
pb.ril <- unique(test.ril$pb)

D.ril <- ld.data.for.heatmaps$ld.RILs.D[10,1]
r.ril <- ld.data.for.heatmaps$ld.RILs.r[10,1]
r2.ril <- ld.data.for.heatmaps$ld.RILs.r2[10,1]

-D.ril / sqrt(pA * pa * pB * pb)

my.r.ril <- (-D.ril / sqrt(pA * pa * pB * pb))
my.r.ril
r.ril

r2.ril
my.r2.ril <- my.r.ril^2
my.r2.ril

cI <- go_df %>%
  dplyr::select(starts_with("I_"))
duplicates <- cI[duplicated(cI), ]

#==============================================================================#
# Test calculation of D | D=p(530A,530B)-p(530A)*p(530B) for loci
# I_113752 asnd I_7660120
#==============================================================================#
geno <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") 

testD <- go_df %>%
  dplyr::select(I_113752, I_7660120) %>%
  dplyr::filter(!is.na(I_113752), !is.na(I_7660120)) %>%
  dplyr::mutate(pA = sum(I_113752 == "H/H", na.rm = T)/n(),
                pB = sum(I_7660120 == "H/H", na.rm = T)/n(),
                pAB = sum(I_7660120 == "H/H" & I_113752 == "H/H", na.rm = T)/n(),
                D = pAB-(pA*pB)) %>%
  distinct(D) %>%
  pull(D)

#==============================================================================#
# Calculate LD with a haplotype object not a genotype object?
#==============================================================================#
# get genotype matrix
gm2 <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>%
  dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
  dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
  tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n.na = sum(is.na(genotype))) %>%
  dplyr::filter(n.na<146) %>% # is this the right na filter?
  dplyr::ungroup() %>%
  dplyr::select(-n.na) %>%
  dplyr::mutate(genotype = case_when(genotype == MA563 ~ "L/L",
                                     genotype == MA530 ~ "H/H",
                                     TRUE ~ NA_character_)) %>%
  dplyr::select(-MA530, -MA563) %>%
  tidyr::pivot_wider(names_from = "strain", values_from = "genotype")

go_list2 <- NULL
for (i in 1:nrow(gm2)) {
  go_list2[[i]] <- genetics::haplotype(as.character(gm2[i, 2:ncol(gm2)]))
}

# make a dataframe from list
go_df2 <- data.frame(go_list2)

# add informative column names to dataframe
colnames(go_df2) <- (gm2$snp_id)


str(go_df2)
glimpse(go_df)
glimpse(go_df2)

# get it in D
ld.all.D.2 <- t(genetics::LD(go_df2)[[2]])
# ugh haplotype is not supported, try diseq?

#==============================================================================#
# haplotype example
#==============================================================================#
# show how genotype and haplotype differ
data1 <- c("C/C", "C/T", "T/C")
data2 <- c("C/C", "T/C", "T/C")
test1 <- genotype( data1 )
test2 <- genotype( data2 )
test3 <- haplotype( data1 )
test4 <- haplotype( data2 )

#==============================================================================#
# Try tsting if the fifth element of the LD output is R2 or not
#==============================================================================#
# get genotype matrix
gm3 <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>%
  dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
  dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
  tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n.na = sum(is.na(genotype))) %>%
  dplyr::filter(n.na<146) %>% # is this the right na filter?
  dplyr::ungroup() %>%
  dplyr::select(-n.na) %>%
  dplyr::mutate(genotype = case_when(genotype == MA563 ~ "L/L",
                                     genotype == MA530 ~ "H/H",
                                     TRUE ~ NA_character_)) %>%
  dplyr::select(-MA530, -MA563) %>%
  tidyr::pivot_wider(names_from = "strain", values_from = "genotype")

go_list3 <- NULL
for (i in 1:nrow(gm3)) {
  go_list3[[i]] <- genetics::genotype(as.character(gm3[i, 2:ncol(gm3)]))
}

# make a dataframe from list
go_df3 <- data.frame(go_list3)

# add informative column names to dataframe
colnames(go_df3) <- (gm3$snp_id)

# get it in r
ld.all<- t(genetics::LD(go_df3))
ld.all.d <- ld.all[[2]]
ld.all.r <- ld.all[[4]]
ld.all.r2 <- ld.all[[5]]

my.ld.all.r2 <- ld.all.r*ld.all.r
test.my.ld <- my.ld.all.r2 - ld.all.r2

# ugh haplotype is not supported, try diseq?

#==============================================================================#
# Test LD D, r and r2 for largest difference r2 value between ALL and RIAILs
#==============================================================================#
# find the pairs of loci with the most missing data between them and calculate LD
# manually
# get genotype matrix
num_lines <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>%
  dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
  dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
  tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
  dplyr::mutate(type = ifelse(grepl(x = strain, pattern = "RIAIL"), "RIAIL", "RIL")) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n.na = sum(is.na(genotype))) %>%
  dplyr::filter(n.na<146) %>% # is this the right na filter?
  dplyr::group_by(snp_id, type) %>%
  dplyr::mutate(n.lines.with.geno.type = sum(!is.na(genotype))) %>%
  dplyr::group_by(snp_id) %>%
  dplyr::mutate(n.lines.with.geno = sum(!is.na(genotype))) %>%
  dplyr::distinct(snp_id, type, .keep_all = T)

