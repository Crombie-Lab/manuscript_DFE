library(tidyverse)

#------------------------------------------------------------------------------#
# Make PLINK files from strain set
#------------------------------------------------------------------------------#
# strain = "all", "ril", "riail"
# geno_path = repo relative path to data/processed/02_DFE_unimputed_genotypes.csv
# out_prefix = prefix for .map .ped files
make_plink_in <- function(geno_path = "data/processed/02_DFE_unimputed_genotypes.csv", strain, out_prefix) {
  if(strain == "all") {
  # get the strains without all missing genotypes
  strains <- data.table::fread(geno_path) %>%
    dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
    dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
    tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
    dplyr::group_by(strain) %>%
    dplyr::mutate(n.na = sum(is.na(genotype))) %>%
    dplyr::filter(n.na<169) %>% # can't have all missing
    dplyr::ungroup() %>%
    dplyr::distinct(strain) %>%
    dplyr::pull(strain)
  }
  if(strain == "ril") {
    # get the strains without all missing genotypes
    strains <- data.table::fread(geno_path) %>%
      dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
      dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
      tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
      dplyr::group_by(strain) %>%
      dplyr::mutate(n.na = sum(is.na(genotype))) %>%
      dplyr::filter(n.na<169) %>% # can't have all missing
      dplyr::ungroup() %>%
      dplyr::filter(grepl(strain, pattern = "RIL_")) %>%
      dplyr::distinct(strain) %>%
      dplyr::pull(strain)
  }
  if(strain == "riail") {
    # get the strains without all missing genotypes
    strains <- data.table::fread(geno_path) %>%
      dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
      dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
      tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
      dplyr::group_by(strain) %>%
      dplyr::mutate(n.na = sum(is.na(genotype))) %>%
      dplyr::filter(n.na<169) %>% # can't have all missing
      dplyr::ungroup() %>%
      dplyr::filter(grepl(strain, pattern = "RIAIL_")) %>%
      dplyr::distinct(strain) %>%
      dplyr::pull(strain)
  }
  
  # Make PLINK .map file
  map <- data.table::fread(geno_path) %>%
    dplyr::select(chrom, full_var_id, pos) %>%
    dplyr::mutate(chrom = dplyr::case_when(chrom == "I" ~ 1,
                                           chrom == "II" ~ 2,
                                           chrom == "III" ~ 3,
                                           chrom == "IV" ~ 4,
                                           chrom == "V" ~ 5,
                                           chrom == "X" ~ 6,
                                           TRUE ~ NA_real_),
                  gen_dist = 0, .before = pos)
  write_tsv(map, file = glue::glue("data/processed/{out_prefix}.map"), col_names = F)
  
  # make the snp-list
  snp_list <- map$full_var_id
  write_file(paste(snp_list, collapse = " "), file = "data/processed/snps.txt")
  
  # Make the .ped file
  #Family ID
  #Sample ID
  #Paternal ID
  #Maternal ID
  #Sex (1=male; 2=female; other=unknown)
  #Affection (0=unknown; 1=unaffected; 2=affected)
  #Genotypes (space or tab separated, 2 for each marker. 0=missing)
  ped <- data.table::fread("data/processed/04_DFE_joined_unimputed_geno_pheno.csv") %>%
    dplyr::filter(grepl(full_id, pattern = "RIAIL|RIL")) %>%
    dplyr::filter(full_id %in% strains) %>% # these are the strains with non-missing genotypes
    dplyr::distinct(full_id, .keep_all = T) %>%
    dplyr::mutate(fam_ID = 100,
                  pat_id = 100,
                  mat_id = 100,
                  sex = 100,
                  aff = 100) %>%
    dplyr::select(fam_ID, samp_id = full_id, pat_id, mat_id, sex, aff, everything(), -block:-ln_ci)
  
  ped[ped==0] <- "1\t1"
  ped[ped==1] <- "2\t2"
  ped[ped==100] <- "0"
  ped[is.na(ped)] <- "0\t0"
  write_tsv(ped, file = glue::glue("data/processed/{out_prefix}.ped"), col_names = F)
}

#------------------------------------------------------------------------------#
# Process PLINK Table output to LD matrix
#------------------------------------------------------------------------------#
#test_plink_proc <- plink_proc(ld = "R2", ld_path = "data/processed/mut_all_allele_ld_r2.ld")
plink_proc <- function(ld, ld_path) {
  # get the ld measures from plink
  ld_df <- data.table::fread(ld_file) %>%
    dplyr::select(SNP_A, SNP_B, ld) %>%
    tidyr::pivot_wider(names_from = SNP_B, values_from = ld) %>% # need to unquote with {{ld}}?
    tibble::column_to_rownames("SNP_A")
  ld_df[upper.tri(ld_df, diag=TRUE)] <- NA
  ldm <- as.matrix(ld_df)
}

#------------------------------------------------------------------------------#
# Intrachromosomal LD calc
#------------------------------------------------------------------------------#
# ld = type of LD output from PLINK c("R" "DP", "D", "R2")
# ld_file = repo relative path to LD file
# plot = T
# filename = NULL, if specified plot will be saved to path
#test_intra <- intra_ld(ld = "DP", ld_file = "data/processed/mut_allele_ld_r.ld", plot = T)
intra_ld <- function(ld, ld_file, plot = T) {
  # get the ld measures from plink
  ld_df <- data.table::fread(ld_file) %>%
  dplyr::select(SNP_A, SNP_B, ld) %>%
  tidyr::pivot_wider(names_from = SNP_B, values_from = ld) %>% # need to unquote with {{ld}}?
  tibble::column_to_rownames("SNP_A")
  ld_df[upper.tri(ld_df, diag=TRUE)] <- NA
  ldm <- as.matrix(ld_df)
  
  # loop through chroms to calculate intraLD
  chroms_list <- NULL
  i = "I"
  for(i in c("I", "II", "III", "IV", "V", "X" )){
    # get I
    chrom_df <- as.data.frame(ldm) %>%
      dplyr::select(starts_with(glue::glue("{i}_"))) %>% # filter columns to chrom 
      tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
      dplyr::filter(stringr::str_detect(POS, pattern = glue::glue("^{i}_"))) %>% #filter to chrom
      dplyr::select(-POS)
    
    chroms_list[[i]] <- tibble::tibble(chrom = i,
                                       ld = ld,
                                       mean_ld = mean(as.matrix(chrom_df), na.rm = T),
                                       sd_ld = sd(as.matrix(chrom_df), na.rm = T),
                                       sem_ld = sd(as.vector(t(chrom_df[!is.na(chrom_df)]))) / sqrt(length(as.vector(t(chrom_df[!is.na(chrom_df)])))))
  }
  
  # combine it
  intra_chrom_ld <- data.table::rbindlist(chroms_list)
  
  if(plot == T) {
  # LD label
  ld_label <- ifelse(ld == "DP", "D'",
                     ifelse(ld == "D", "D",
                            ifelse(ld == "R", "r",
                                   ifelse(ld == "R2", expression(r^2), "wtf?"))))
                                   
  intra_chrom_ld_plot <- ggplot(intra_chrom_ld) +
    aes(x = factor(chrom, levels = c("I", "II", "III", "IV", "V", "X")), y = mean_ld) +
    geom_bar(stat = "identity", width = 0.5, position = "dodge") +
    #scale_fill_manual(labels = c("all strains", "RIAILs", "RILs"), values = c("grey40", "#00B9F1", "#00A875")) +
    geom_errorbar(aes(ymin=mean_ld - sem_ld, ymax=mean_ld + sem_ld), width=.25, position=position_dodge(.5), size = 0.25) +
    theme_bw() +
    labs(x = "Chromosome", y = glue::glue('Mean intra-chromosomal LD ({ld_label})'), fill = "") +
    theme(legend.position = "none")
  intra_chrom_ld_plot
  }
  
  # output
  if(plot == T) {
  out <- list(intra_chrom_ld = intra_chrom_ld, intra_chrom_ld_plot = intra_chrom_ld_plot)
  }
  if(plot == F) {
    out <- intra_chrom_ld
  }
  return(out)
}

#------------------------------------------------------------------------------#
# Interchromosomal LD calc
#------------------------------------------------------------------------------#
# mean interchrom ld function. Pairs is set of patterns to match chroms, see make pairs below function.
#data = "data/processed/mut_allele_ld_r.ld"
#test_inter <- inter_ld(ld = "DP", ld_file = "data/processed/mut_allele_ld_r.ld")
inter_ld <- function(ld, ld_file) {
  # read in data
  # get the ld measures from plink
  # get the ld measures from plink
  ld_df <- data.table::fread(ld_file) %>%
    dplyr::select(SNP_A, SNP_B, ld) %>%
    tidyr::pivot_wider(names_from = SNP_B, values_from = ld) %>% # need to unquote with {{ld}}?
    tibble::column_to_rownames("SNP_A")
  ld_df[upper.tri(ld_df, diag=TRUE)] <- NA
  ldm <- as.matrix(ld_df)
  
  # Make pairs
  pairs <- as.data.frame(t(combn(c("^I_", "^II_", "^III_", "^IV_", "^V_", "^X_"),2)))
  
  # make a list
  lds <- NULL
  # loop through the pairs
  for(i in 1:nrow(pairs)){
    l <- pairs[[i,1]]
    r <- pairs[[i,2]]
    
    # shape data 
    d <- as.data.frame(ldm) %>%
      dplyr::select_if(grepl(l, names(.))) %>% # filter columns to chrom 
      tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
      dplyr::filter(stringr::str_detect(POS, pattern = r)) %>% #filter to chrom
      dplyr::select(-POS)
    
    # calculate mean
    mean_d <- tibble::tibble(chr1 = stringr::str_replace_all(l, pattern = "\\^|_", replacement = ""),
                             chr2 = stringr::str_replace_all(r, pattern = "\\^|_", replacement = ""),
                             mean_ic_ld = mean(as.matrix(d), na.rm = T),
                             sd_ic_ld = sd(as.matrix(d), na.rm = T),
                             sem_ic_ld = sd_ic_ld / sqrt(length(as.vector(t(d[!is.na(d)])))))
    
    # Add to list
    lds[[i]] <- mean_d
  }
  # return the list as dataframe
  lds_df <- data.table::rbindlist(lds) %>%
    dplyr::mutate(chrom_pair = paste0(chr1, "-", chr2))
  
  # return
  return(lds_df)
}

#============================================#
# Make LD heatmap with D
#============================================#
#test <- LDchroms2(ld = "R2", ld_file = "data/processed/mut_all_allele_ld_r2.ld")
#test2 <- LDchroms2(ld = "DP", ld_file = "data/processed/mut_all_allele_ld_r.ld")
LDchroms2 <- function(ld, ld_file){
  # get the parents for plotting using hard coded parent_genos.csv file
  parent_genos <- data.table::fread("data/raw/parent_genos.csv")
  parents <- parent_genos %>%
    dplyr::mutate(loci = stringr::str_replace(loci, pattern = regex("_[^_]*$"), replacement = ""),
                  label = ifelse(MA530 == 1, "MA530", "MA563"))
  
  # get the ld measures from plink
  ld_df <- data.table::fread(ld_file) %>%
    dplyr::select(SNP_A, SNP_B, ld) %>%
    tidyr::pivot_wider(names_from = SNP_B, values_from = ld) %>% # need to unquote with {{ld}}?
    tibble::column_to_rownames("SNP_A")
  ld_df[upper.tri(ld_df, diag=TRUE)] <- NA # diag = TRUE original
  data <- as.matrix(ld_df)
  
  # shape the data for each chromosome
  d <- reshape2::melt(data, na.rm = TRUE)
  
  # get each chrom
  chr.df <- d %>%
    dplyr::mutate(c1 = stringr::str_extract(Var1, pattern = regex("[^_]+")),
                  c2 = stringr::str_extract(Var2, pattern = regex("[^_]+"))) %>%
    dplyr::filter(c1 == c2) %>%
    dplyr::left_join(parents, by = c("Var2" = "loci")) %>%
    dplyr::mutate(Var2 = factor(Var2, levels = levels(d$Var2)))
  
  # get unique chrs
  chrs <- unique(chr.df$c1)
  
  # set the color range for LD type in data and ld label
  if(ld == "D"){
    ld.limits <- c(-0.25,0.25)
  }
  if(ld %in% c("R", "DP")){
    ld.limits <- c(-1,1)
  }
  if(ld == "R2"){
    ld.limits <- c(0,1)
  }
  if(!(ld %in% c("R2", "DP", "R", "D"))){
    stop("type not recognized, Please specify any of R2, DP, R, or D")
  }
  
  # get a plot list
  plot.list <- NULL
  
  # loop through the chrs and make plots
  for(i in unique(chrs)){
    # get loci specific data
    l.i <- parents %>%
      dplyr::filter(grepl(loci, pattern = paste0("^",i, "_"))) %>%
      tidyr::separate(loci, into = c("chrom", "pos"), remove = F) %>%
      dplyr::mutate(pos = as.numeric(pos),
                    min.pos = min(pos),
                    max.pos = max(pos),
                    n.pos = length(unique(pos)),
                    pos.num = 1:n(),
                    std.dist = (max.pos - min.pos)/(n.pos-1),
                    std.pos = case_when(pos == min.pos ~ pos,
                                        (pos != min.pos & pos != max.pos) ~ (min.pos + (pos.num-1)*std.dist),
                                        pos == max.pos ~ max.pos))
    
    # get chrom specific data
    c.i <- chr.df %>%
      dplyr::filter(c1 == i)
    
    # make heatmap plot 
    c.p <- ggplot(c.i, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(na.rm = TRUE) +
      viridis::scale_fill_viridis(option = "viridis", limits = ld.limits) +
      #scales::rescale(x, to = c(-1, 1), from = range(x)) +
      theme_void() +
      coord_equal() +
      theme(legend.position = "none",
            plot.margin = margin(0, 0, 0, 0, "pt")) #upper, right, bottom, left
    
    
    # flip it around
    c.p2_v2 <- ggplotify::as.ggplot(c.p, scale = 1, hjust = 0, vjust = 0.4625, angle = 315) 
    
    # make a line plot
    l.p <- ggplot(l.i) +
      #geom_point(aes(y = 0, x = std.pos, shape = label)) +
      geom_segment(aes(x = std.pos, xend = pos, y = 0, yend = 1, color = label, linetype = label)) +
      scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
      ggnewscale::new_scale_color() +
      geom_segment(aes(x = unique(min.pos), xend = unique(max.pos), y = 1, yend = 1)) +
      theme_void() +
      labs(title = glue::glue("{i}"), color = "Parent") +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      theme(plot.margin = margin(0, 0, 0, 0, "pt")) #upper, right, bottom, left
    
    # put them together
    p.full <- cowplot::plot_grid(l.p, c.p2_v2, ncol = 1, rel_heights = c(1,10), axis = "lrtb", align = "h")
    
    # add to list
    plot.list[[i]] <- p.full
  }
  
  # get nice legends
  l.legend <- cowplot::get_legend(l.p + theme(legend.position = "left") + labs(linetype = "Parent"))
  if(ld == "D"){
    c.legend <- cowplot::get_legend(c.p + theme(legend.position = "left") + labs(fill = expression(r^2)))  
  }
  if(ld == "DP"){
    c.legend <- cowplot::get_legend(c.p + theme(legend.position = "left") + labs(fill = "D'"))
  }
  if(ld == "R"){
    c.legend <- cowplot::get_legend(c.p + theme(legend.position = "left") + labs(fill = "r"))
  }
  if(ld == "R2"){
    c.legend <- cowplot::get_legend(c.p + theme(legend.position = "left") + labs(fill = expression(r^2)))
  }
  
  legend.full <- cowplot::plot_grid(l.legend, c.legend, ncol = 1)
  
  # make a full plot
  p.full.ld <- cowplot::plot_grid(plotlist = plot.list, ncol = 3, nrow = 2)
  p.final <- cowplot::plot_grid(p.full.ld, legend.full, ncol = 2, rel_widths = c(8,1))
  
  # return it
  return(p.final)
}
