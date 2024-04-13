library(tidyverse)

# set colors
colors <- c("grey40", "#00B9F1", "#00A875")
names(colors) <- c("all strains", "RIAILs", "RILs")
colors

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