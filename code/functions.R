library(tidyverse)

# set colors
colors <- c("grey40", "#00B9F1", "#00A875")
names(colors) <- c("all strains", "RIAILs", "RILs")
colors

#============================================#
# checkImput - test performance of imputation
#              by dropping called variants and
#              imputing without concern for 
#              distnace.
#============================================#
# checkImput <- function(data, strains){
#   # establish list
#   list <- list()
#   # make data a tibble
#   data <- tibble::as_tibble(data)
#   # make a progress bar
#   pb <- progress::progress_bar$new(total = length(unique(strains)))
#   # loop through strains
#   for(i in strains){
#     # filter to strain
#     s <- data %>%
#       dplyr::select(any_of(x = c("chrom", "pos", "var_type", "full_var_id", i)))
#     # loop through chroms
#     for(j in unique(s$chrom)){
#       # filter to chrom and remove missing and setup dataframe
#       sc <- s %>%
#         dplyr::mutate(line = names(.[,5]),
#                       correct = NA,
#                       dist = NA_integer_,
#                       imp_geno = NA_integer_) %>%
#         dplyr::select(chrom:full_var_id, geno = 5, imp_geno, dist, correct, line) %>% 
#         dplyr::filter(chrom == j & !is.na(geno))
#       # handle empty chrom or chrom with a single called variant
#       if(nrow(sc) == 0 | nrow(sc) == 1){
#         next()
#       }
#       # handle non-empty chrom
#       else {
#         for(k in unique(sc$full_var_id)){
#           # handle first var on chrom always impute as the next variant
#           if(k == sc %>%
#              dplyr::slice_head(n = 1) %>%
#              dplyr::pull(full_var_id)) {
#             # get frist var
#             scv <- sc %>%
#               dplyr::slice_head(n = 2) %>%
#               dplyr::mutate(imp_geno = lead(geno),
#                             dist = lead(pos) - pos,
#                             correct = ifelse(imp_geno == geno, T, F)) %>%
#               dplyr::slice_head(n = 1)
#           }
#           # handle last var in chrom
#           if(k == sc %>%
#              dplyr::slice_tail(n = 1) %>%
#              dplyr::pull(full_var_id)) {
#             # get last var
#             scv <- sc %>%
#               dplyr::slice_tail(n = 2) %>%
#               dplyr::mutate(imp_geno = lag(geno),
#                             dist =  pos - lag(pos),
#                             correct = ifelse(imp_geno == geno, T, F)) %>%
#               dplyr::slice_tail(n = 1)
#           }
#           # handle other vars in chrom
#           if(!(k %in% c(sc %>% dplyr::slice_tail(n = 1) %>% dplyr::pull(full_var_id),
#                         sc %>% dplyr::slice_head(n = 1) %>% dplyr::pull(full_var_id)))) {
#             # get other var
#             # get row number of var
#             row <- which(grepl(k, sc$full_var_id))
#             scv <- sc %>%
#               dplyr::slice((row-1):(row+1)) %>%
#               dplyr::mutate(imp_geno_lag = lag(geno),
#                             dist_lag =  pos - lag(pos),
#                             imp_geno_lead = lead(geno),
#                             dist_lag = pos - lag(pos),
#                             dist_lead = lead(pos) - pos,
#                             imp_geno = ifelse(dist_lag >= dist_lead, imp_geno_lead, imp_geno_lag),
#                             dist = ifelse(dist_lag >= dist_lead, dist_lead, dist_lag),
#                             correct = ifelse(imp_geno == geno, T, F)) %>%
#               dplyr::select(1:9) %>%
#               dplyr::slice(2)
#           }
#           # add to output
#           list[[length(list) + 1]] <- scv
#           
#         }
#       }
#     }
#     # add a tick
#     pb$tick()
#   }
#   # tidy up the output
#   out <- data.table::rbindlist(list)
#   # message
#   message("DONE")
#   # return
#   return(out)
# }

# UGh, redo imputation to reflect parental genotypes not just REF or ALT
checkImput2 <- function(data, strains){
  # establish list
  list <- list()
  # make data a tibble
  data <- tibble::as_tibble(data)
  # make a progress bar
  pb <- progress::progress_bar$new(total = length(unique(strains)))
  # loop through strains
  for(i in strains){
    # filter to strain
    s <- data %>%
      dplyr::select(any_of(x = c("chrom", "pos", "var_type", "full_var_id", "MA530", "MA563", i)))
    # loop through chroms
    for(j in unique(s$chrom)){
      # filter to chrom and setup dataframe
      sc <- s %>%
        dplyr::mutate(line = names(.[,7]),
                      correct = NA,
                      dist = NA_integer_,
                      imp_geno = NA_integer_) %>%
        dplyr::select(chrom:MA563, geno = 7, imp_geno, dist, correct, line) %>% 
        dplyr::filter(chrom == j & !is.na(geno))
      # handle empty chrom or chrom with a single called variant
      if(nrow(sc) == 0 | nrow(sc) == 1){
        next()
      }
      # handle non-empty chrom
      else {
        for(k in unique(sc$full_var_id)){
          # handle first var on chrom always impute as the next variant
          if(k == sc %>%
             dplyr::slice_head(n = 1) %>%
             dplyr::pull(full_var_id)) {
            # get frist var
            scv <- sc %>%
              dplyr::slice_head(n = 2) %>%
              dplyr::mutate(imp_geno = case_when(lead(geno) == MA530 ~ lead(MA530),
                                                  lead(geno) == MA563 ~ lead(MA563),
                                                  TRUE ~ NA_integer_),
                            dist = lead(pos) - pos,
                            correct = ifelse(imp_geno == geno, T, F)) %>%
              dplyr::slice_head(n = 1)
          }
          # handle last var in chrom
          if(k == sc %>%
             dplyr::slice_tail(n = 1) %>%
             dplyr::pull(full_var_id)) {
            # get last var
            scv <- sc %>%
              dplyr::slice_tail(n = 2) %>%
              dplyr::mutate(imp_geno = case_when(lag(geno) == MA530 ~ lag(MA530),
                                                 lag(geno) == MA563 ~ lag(MA563),
                                                 TRUE ~ NA_integer_),
                            dist =  pos - lag(pos),
                            correct = ifelse(imp_geno == geno, T, F)) %>%
              dplyr::slice_tail(n = 1)
          }
          # handle other vars in chrom
          if(!(k %in% c(sc %>% dplyr::slice_tail(n = 1) %>% dplyr::pull(full_var_id),
                        sc %>% dplyr::slice_head(n = 1) %>% dplyr::pull(full_var_id)))) {
            # get other var
            # get row number of var
            row <- which(grepl(k, sc$full_var_id))
            scv <- sc %>%
              dplyr::slice((row-1):(row+1)) %>%
              dplyr::mutate(#imp_geno_lag = lag(geno),
                            imp_geno_lag = case_when(lag(geno) == MA530 ~ lag(MA530),
                                                 lag(geno) == MA563 ~ lag(MA563),
                                                 TRUE ~ NA_integer_),
                            dist_lag =  pos - lag(pos),
                            #imp_geno_lead = lead(geno),
                            imp_geno_lead = case_when(lead(geno) == MA530 ~ lead(MA530),
                                                 lead(geno) == MA563 ~ lead(MA563),
                                                 TRUE ~ NA_integer_),
                            dist_lag = pos - lag(pos),
                            dist_lead = lead(pos) - pos,
                            imp_geno = ifelse(dist_lag >= dist_lead, imp_geno_lead, imp_geno_lag),
                            dist = ifelse(dist_lag >= dist_lead, dist_lead, dist_lag),
                            correct = ifelse(imp_geno == geno, T, F)) %>%
              dplyr::select(1:11) %>%
              dplyr::slice(2)
          }
          # add to output
          list[[length(list) + 1]] <- scv
          
        }
      }
    }
    # add a tick
    pb$tick()
  }
  # tidy up the output
  out <- data.table::rbindlist(list)
  # message
  message("DONE")
  # return
  return(out)
}
#============================================#
# Bootstrap to find mutational effects
#============================================#
mutEffBoots <- function(data, nboots=1000, seed = 99) {
  # bootstap lines with replacement to find confidence intervals for the mutational effect of each variant
  mut.effs <- NULL
  pb <- progress::progress_bar$new(total = length(unique(data$variant))*nboots)
  pb$tick(0)
  
  # set the seed
  set.seed(seed)
  # do the loop
  for(i in unique(data$variant)) {
    # filter to variant and remove nas for genotype
    var <- data %>%
      dplyr::filter(variant == i & !is.na(genotype))
    
    # get mutant means
    mut_means <- var %>%
      dplyr::filter(genotype == 1) %>%
      dplyr::pull(reg.line_mean_fitness)
    # get number of samples
    #mut_n <- as.numeric(length(mut_means))
    
    # get ref means
    ref_means <- var %>%
      dplyr::filter(genotype == 0) %>%
      dplyr::pull(reg.line_mean_fitness)
    # get number of samples
    #ref_n <- as.numeric(length(ref_means))
    
    mut_effect_boots <- NULL
    for (j in 1:nboots) {
      # resample
      mut_resample <- sample(mut_means, replace = TRUE)
      ref_resample <- sample(ref_means, replace = TRUE)
      
      # calculate means from resampled data
      mut_mean <- mean(mut_resample)
      ref_mean <- mean(ref_resample)
      
      # find mutational effect
      mut_effect <- mut_mean - ref_mean
      mut_effect_boots[j] <- mut_effect
      
      # add a tick
      pb$tick()
    }
    # get mean and CIs from mut_effect_boots
    boots_dat <- tibble::tibble(variant = i,
                                mut_effect_point_estimate = unique(var$reg.mean_mutant_fitness_effect_mean_scaled),
                                mean_mut_effect_from_boots = mean(mut_effect_boots),
                                lower.05 = unname(quantile(mut_effect_boots, probs = .05)),
                                upper.95 = unname(quantile(mut_effect_boots, probs = .95)),
                                n_mut_lines = as.numeric(length(mut_means)),
                                n_ref_lines = as.numeric(length(ref_means)))
    
    # put these in var output
    mut.effs[[i]] <- boots_dat
  }
  # return the data as a nice dataframe
  out <- data.table::rbindlist(mut.effs)
  return(out)
}  

#=================================================#
# Impute with sepecified distance between markers
#=================================================#
# impute <- function(data, strains, impDist = 10000){
#   # establish list
#   list <- list()
#   # make data a tibble
#   data <- tibble::as_tibble(data)
#   # Create the progress bar
#   pb <- progress::progress_bar$new(total = length(unique(strains))*nrow(data))
#   # loop through strains
#   for(i in strains){
#     # filter to strain
#     s <- data %>%
#       dplyr::select(any_of(x = c("chrom", "pos", "var_type", "full_var_id", i)))
#     # loop through chroms
#     for(j in unique(s$chrom)){
#       # filter to chrom and remove missing and setup dataframe
#       sc <- s %>%
#         dplyr::mutate(line = names(.[,5])) %>%
#         dplyr::select(chrom:full_var_id, geno = 5, line) %>% 
#         dplyr::filter(chrom == j)
#       for(k in unique(sc$full_var_id)){
#         # filter to the variant
#         scv <- sc %>%
#           dplyr::filter(full_var_id == k) %>%
#           dplyr::mutate(imp_geno = NA_integer_,
#                         imp_dist = NA_integer_)
#         # find the closest genotyped variant
#         closest <- sc %>%
#           dplyr::filter(full_var_id != k & !is.na(geno)) %>% 
#           dplyr::mutate(imp_dist = as.integer(abs(scv$pos - pos))) %>%
#           dplyr::arrange(imp_dist) %>%
#           dplyr::slice(1)
#         # handle empty cases
#         if(nrow(closest) == 0){
#           scv_imp <- scv 
#         }
#         else {
#           scv_imp <- scv %>%
#             dplyr::mutate(imp_geno = case_when(is.na(geno) ~ closest$geno,
#                                                TRUE ~ geno),
#                           imp_dist = case_when(is.na(geno) ~ closest$imp_dist,
#                                                TRUE ~ imp_dist))
#         }
#         # add to output
#         list[[length(list) + 1]] <- scv_imp
#         pb$tick()
#       }
#     }
#   }
#   # tidy up the output
#   out1 <- data.table::rbindlist(list)
#   out2 <- out1 %>%
#     dplyr::mutate(imputed = case_when(is.na(geno) & imp_dist <= impDist ~ T,
#                                       TRUE ~ F),
#                   geno = case_when(is.na(geno) & imp_dist <= impDist ~ imp_geno,
#                                    TRUE ~ geno)) %>%
#     dplyr::select(chrom:line, imputed)
#   
#   return(out2)
# }

# # Need to use the parental genotypes!!! Doh
# impute2 <- function(data, strains, impDist = 10000){
#   # establish list
#   list <- list()
#   # make data a tibble
#   data <- tibble::as_tibble(data)
#   # Create the progress bar
#   pb <- progress::progress_bar$new(total = length(unique(strains))*nrow(data))
#   # loop through strains
#   for(i in strains){
#     # filter to strain
#     s <- data %>%
#       dplyr::select(any_of(x = c("chrom", "pos", "var_type", "full_var_id", "MA530", "MA563", i)))
#     # loop through chroms
#     for(j in unique(s$chrom)){
#       # filter to chrom and remove missing and setup dataframe
#       sc <- s %>%
#         dplyr::mutate(line = names(.[,7])) %>%
#         dplyr::select(chrom:MA563, geno = 7, line) %>% 
#         dplyr::filter(chrom == j)
#       for(k in unique(sc$full_var_id)){
#         # filter to the variant
#         scv <- sc %>%
#           dplyr::filter(full_var_id == k) %>%
#           dplyr::mutate(imp_geno = NA_integer_,
#                         imp_dist = NA_integer_)
#         # find the closest genotyped variant
#         closest <- sc %>%
#           dplyr::filter(full_var_id != k & !is.na(geno)) %>% 
#           dplyr::mutate(imp_dist = as.integer(abs(scv$pos - pos))) %>%
#           dplyr::arrange(imp_dist) %>%
#           dplyr::slice(1) %>%
#           dplyr::mutate(parent = case_when(geno == MA530 ~ "MA530",
#                                            geno == MA563 ~ "MA563"))
#         # handle empty cases
#         if(nrow(closest) == 0){
#           scv_imp <- scv 
#         }
#         else {
#           scv_imp <- scv %>%
#             dplyr::mutate(imp_geno = case_when(is.na(geno) & closest$parent == "MA563" ~ as.integer(MA563),
#                                                is.na(geno) & closest$parent == "MA530" ~ as.integer(MA530),
#                                                TRUE ~ geno),
#                           imp_dist = case_when(is.na(geno) ~ closest$imp_dist,
#                                                TRUE ~ imp_dist))
#         }
#         # add to output
#         list[[length(list) + 1]] <- scv_imp
#         pb$tick()
#       }
#     }
#   }
#   # tidy up the output
#   out1 <- data.table::rbindlist(list)
#   out2 <- out1 %>%
#     dplyr::mutate(imputed = case_when(is.na(geno) & imp_dist <= impDist ~ T,
#                                       TRUE ~ F),
#                   geno = case_when(is.na(geno) & imp_dist <= impDist ~ imp_geno,
#                                    TRUE ~ geno)) %>%
#     dplyr::select(chrom:line, imputed)
#   
#   return(list(out1, out2))
# }

impute3 <- function(data, strains){
  # establish chrom list
  c.list <- list()
  # make data a tibble
  data <- tibble::as_tibble(data)
  # Create the progress bar
  pb <- progress::progress_bar$new(format = "  imputing [:bar] :percent eta: :eta",
                                   clear = FALSE,
                                   total = length(unique(strains))*length(unique(data$chrom)))
  # loop through chroms
  for(i in unique(data$chrom)){
    # get chrom dataset
    c <- data %>%
      dplyr::filter(chrom == i) %>%
      dplyr::select(-N2_G0_ANC) # remove ancestor
    
    # make a list missing and present vars per strain
    strains2 <- match(strains,names(c))
    sc.list <- NULL
    for(j in min(strains2):max(strains2)){
      # get missing variants for a strain
      ms <- c[which(is.na(c[j])),]
      
      # get the present variants for a strain
      ps <- c[which(!is.na(c[j])),]
      
      # move on if everything is missing
      if(nrow(ps) == 0 | nrow(ms) == 0) {
        # add a tick for strain and chrom
        pb$tick()
        next()
      }
      
      # get the distance from each missing variant to each present variant
      distm <- abs(outer(ms$pos,ps$pos, `-`))
      
      # get a list of the imputed loci
      i.list <- NULL
      for(k in 1:nrow(distm)){
        # bind the closest present locus as a pair with missing locus. The present
        # locus is always frist (A), this means locus (B) is always missing in a two locus model (AB). 
        l.pair <- dplyr::bind_rows(ps[which.min(distm[k,]),], ms[k,]) 
        
        # calculate haplotype frequencies pAB, pAb, paB, pab: AB=530 ab=563
        AB <- t(l.pair[5][, colSums(is.na(l.pair[5])) == 0])
        ab <- t(l.pair[6][, colSums(is.na(l.pair[6])) == 0])
        Ab <- paste0(AB[1], ab[2])
        aB <- paste0(ab[1], AB[2])
        # cleanup 
        AB <- paste0(AB, collapse = "")
        ab <- paste0(ab, collapse = "")
        # get the frequencies
        freq <- as.data.frame(t(l.pair[7:ncol(l.pair)][, colSums(is.na(l.pair[7:ncol(l.pair)])) == 0])) %>%
          dplyr::mutate(haplo = case_when(paste0(V1, V2) == AB ~ "AB",
                                           paste0(V1, V2) == ab ~ "ab",
                                           paste0(V1, V2) == Ab ~ "Ab",
                                           paste0(V1, V2) == aB ~ "aB",
                                           TRUE ~ "wtf?")) %>%
          dplyr::group_by(haplo) %>%
          dplyr::summarise(n = n(), .groups = "drop") %>%
          dplyr::mutate(p = n/sum(n)) %>%
          dplyr::arrange(haplo) # ensure haplotypes are in the right order
        
        # fix missing if necessary
        if(nrow(freq) != 4){
          ap <- freq$haplo # get haplos present
          ae <- c("ab", "aB", "Ab", "AB") # get haplotypes expected
          am <- which(!(ae %in% ap)) # get index missing
          # make the missing data
          freq.add <- tibble::tibble(allele = ae[am],
                                     n = rep(0, times = length(am)),
                                     p = n)
          # add missing in
          freq <- freq %>%
            dplyr::bind_rows(freq.add) %>%
            dplyr::arrange(haplo) # ensure haplotypes are in the right order
        }
        
        # get the probabilities
        gApB <- freq[[3]][[4]]/(freq[[3]][[4]] + freq[[3]][[3]]) 
        gApb <- freq[[3]][[3]]/(freq[[3]][[3]] + freq[[3]][[4]])
        gapB <- freq[[3]][[2]]/(freq[[3]][[2]] + freq[[3]][[1]])
        gapb <- freq[[3]][[1]]/(freq[[3]][[1]] + freq[[3]][[2]])
        
        # calc haplotype probabilities for  missing locus (B)
        imp <- l.pair[c(2,4:6, j)] %>%
          dplyr::mutate(pMA530 = case_when(.[[5]][[1]] == .[[3]][[1]] ~ gApB, 
                                           .[[5]][[1]] == .[[4]][[1]] ~ gapB,
                                           TRUE ~ NA_real_),
                        pMA563 = case_when(.[[5]][[1]] == .[[3]][[1]] ~ gApb, 
                                           .[[5]][[1]] == .[[4]][[1]] ~ gapb,
                                           TRUE ~ NA_real_),
                        imp2 = case_when(MA530 == 1 ~ rbinom(n = 1, size = 1, prob = pMA530),
                                         MA530 == 0 ~ rbinom(n = 1, size = 1, prob = pMA563),
                                         TRUE ~ NA_integer_), # rbinom just set B to
                        imp.d = distm[k, which.min(distm[k,])],
                        closest.locus = .[[2]][[1]],
                        strain = names(.[5]),
                        chrom = i) %>%
          dplyr::slice_tail(n = 1)
        
        # impute locus B
        imp[[5]] <- imp[[8]]
        #rename geontype column
        colnames(imp)[5] <- "geno"
        # shape output
        ivs <- imp %>%
          dplyr::select(chrom, pos, full_var_id, geno, strain,
                        MA530, MA563, pMA530, pMA563, imp.d, closest.locus)
        
        # add to list
        i.list[[k]] <- ivs
      }
      
      # bind all the imputed loci for strain in this chrom
      isc <- data.table::rbindlist(i.list)
      
      # add a tick for strain and chrom
      pb$tick()
      
      # store data
      s <- names(c[j])
      sc.list[[s]] <- isc
    }
    # wrap up chrom
    sc.df <- data.table::rbindlist(sc.list)
    
    # add to c.list
    c.list[[i]] <- sc.df
  }
  # wrap up all chroms
  full.df <- data.table::rbindlist(c.list)
  
  # return
  return(full.df)
  
}

#============================================#
# Make LD heatmap with D
#============================================#
#data <- ld.RILs.r2
#parents <- pg.2
#type <- "R2"
LDchroms <- function(data, parents, type){
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
  
  # set the color range for LD type in data
  if(type == "D"){
    ld.limits <- c(-0.25,0.25)
  }
  if(type == "r"){
    ld.limits <- c(-1,1)
  }
  if(type %in% c("R2", "D'")){
    ld.limits <- c(0,1)
  }
  if(!(type %in% c("R2", "D'", "r", "D"))){
    stop("type not recognized, Please specify any of R2, D', r, or D")
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
      viridis::scale_fill_viridis(option = "D", limits = ld.limits) +
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
      geom_segment(aes(x = std.pos, xend = pos, y = 0, yend = 1, linetype = label)) +
      geom_segment(aes(x = unique(l.i$min.pos), xend = unique(l.i$max.pos), y = 1, yend = 1)) +
      theme_void() +
      labs(title = paste0("Chromosome ", i)) +
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
  c.legend <- cowplot::get_legend(c.p + theme(legend.position = "left") + labs(fill = glue::glue("  {type}")))
  legend.full <- cowplot::plot_grid(l.legend, c.legend, ncol = 1)
  
  # make a full plot
  p.full.ld <- cowplot::plot_grid(plotlist = plot.list, ncol = 3, nrow = 2)
  p.final <- cowplot::plot_grid(p.full.ld, legend.full, ncol = 2, rel_widths = c(8,1))
  
  # return it
  return(p.final)
}