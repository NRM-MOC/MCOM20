plot_ts <- function(mocis_object){
  library(colorspace)
  index <- mocis_contents(mocis_object) %>% 
    filter(!((loc == "ORFJ") & (gen == "CLUP"))) # No Herring in Örefjärden!!!
  n_rows <- nrow(index)
  for (row in 1:n_rows){
    file_name <- paste0(index[row, ], collapse = "_")
    pdf(paste0("figs/", file_name, ".pdf"), width = 5, height = 5)
    plot_mocis(mocis_object, 
               var = index[["var"]][row],
               loc = index[["loc"]][row],
               genus = index[["gen"]][row], 
               cex.mtext = 1.2, 
               newlimit = ifelse((index[["gen"]][row] == "GADU") & (!(index[["var"]][row] %in% c("HG", "PB", "CD"))), 0, NA)) # No limits for organic Cod
    dev.off()
    png(paste0("figs/", file_name, ".png"))
    plot_mocis(mocis_object, 
               var = index[["var"]][row],
               loc = index[["loc"]][row],
               genus = index[["gen"]][row], 
               cex.mtext = 1.2, 
               newlimit = ifelse((index[["gen"]][row] == "GADU") & (!(index[["var"]][row] %in% c("HG", "PB", "CD"))), 0, NA)) # No limits for organic Cod
    dev.off()
  }
}

save_tables <- function(mocis_object, report_year = 2019){
  all_vars <- mocis_tibble(mocis_object)  %>% 
    filter(!((loc == "ORFJ") & (gen == "CLUP"))) %>% # No Herring in Örefjärden!!!
    group_by(var) %>% 
    nest()
  nrows <- nrow(all_vars)
  for (row in 1:nrows){
    data <- all_vars$data[[row]] %>% 
      mutate(var = all_vars$var[row],
             lim = map_dbl(limit, ~.x$tv[1]),
             EKOO = map_dbl(aggdata, ~.x[["EKOO"]][1]),
             NKOO = map_dbl(aggdata, ~.x[["NKOO"]][1]),
             lat = ne2latlong(NKOO, EKOO)[, 1],
             long = ne2latlong(NKOO, EKOO)[, 2],
             station_code = map_chr(aggdata, ~.x[["LOC"]][1]),
             station = mocis_name_station(station_code),
             station = fct_reorder(station, -(-1)^(long > 13.6)*lat),
             slope = map_dbl(linmod10, ~ifelse(is.null(.x[["slope"]]), NA, .x[["slope"]]))
      ) %>% arrange(desc(station))
    herring_data <-  data %>% 
      filter(gen == "CLUP")
    fish_data <- data %>% 
      filter(gen %in% c("PERC", "GADU", "ZOAR"))
    mussel_data <- data %>% 
      filter(gen == "MYTI")
    egg_data <- data %>% 
      filter(gen %in% c("SIGR", "STER", "HAEM"))
    if (nrow(herring_data) > 0){
      mocis_table(herring_data) %>%
        select(-var, -gen, -loc) %>% 
        write_csv(path = paste0("tables/", all_vars$var[row], "_herring_table.csv"))
    }
    if (nrow(fish_data) > 0){
      mocis_table(fish_data) %>%
        select(-var, -gen, -loc) %>% 
        write_csv(path = paste0("tables/", all_vars$var[row], "_fish_table.csv"))
    }
    if (nrow(mussel_data) > 0){
      mocis_table(mussel_data) %>%
        select(-var, -gen, -loc) %>% 
        write_csv(path = paste0("tables/", all_vars$var[row], "_mussel_table.csv"))
    }
    if (nrow(egg_data) > 0){
      mocis_table(egg_data) %>%
        select(-var, -gen, -loc) %>% 
        write_csv(path = paste0("tables/", all_vars$var[row], "_egg_table.csv"))
    }
  }
}


plot_mean_maps <- function(mocis_object, report_year = 2019){
  data <- mocis_tibble(mocis_object) %>% filter(gen == "CLUP", loc != "ORFJ")
  if (nrow(data) == 0) return()
  world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
    st_crop(c(xmin=8, xmax=25, ymin=54, ymax=70))
  #lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf") %>% 
  #  st_crop(c(xmin=8, xmax=25, ymin=54, ymax=70))
  plot_data <- data %>% mutate(
    EKOO = map_dbl(aggdata, ~.x[["EKOO"]][1]),
    NKOO = map_dbl(aggdata, ~.x[["NKOO"]][1]),
    lat = ne2latlong(NKOO, EKOO)[, 1],
    long = ne2latlong(NKOO, EKOO)[, 2],
    long = if_else(loc %in% c("ANGV", "UTLV"), long + 0.5, long),
    lim = map_dbl(limit, ~.x$tv[1]),
    slope = map_dbl(linmod10,  ~ifelse(is.null(.x[["slope"]]), NA, .x[["slope"]])),
    #mean_Conc = map_dbl(tidyaggdata, ~filter(.x, YEAR %in% (report_year-3):report_year) %>% 
    mean_Conc = map_dbl(tidyaggdata, ~filter(.x, YEAR %in% (report_year-2):report_year) %>% 
                          summarise(mean_Conc = exp(mean(log(value)))) %>% 
                          pull(mean_Conc)),
    station = mocis_name_station(loc),
    station = fct_reorder(station, -(-1)^(long > 13.6)*lat)) %>% 
    select(gen, var, station, slope, lim, mean_Conc, lat, long) %>% 
    filter(!is.na(mean_Conc)) %>% 
    arrange(desc(station)) %>% 
    group_by(var) %>% 
    nest()
  vars <- plot_data[["var"]]
  for (idx in (1:length(vars))){
    sub_data <- plot_data %>% pull(data) %>% .[[idx]]
    p <- plot_grid(ggplot(world_map) + 
                     geom_sf(color = "white", fill = "grey") + 
                     #geom_sf(data = lakes, fill = "white", color = "grey", size = 0) +
                     annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
                     annotation_scale() +
                     coord_sf(xlim = c(8, 25), ylim = c(54, 70), expand = FALSE) +
                     geom_point(data = sub_data, 
                                aes(x = long, y = lat, 
                                    fill = mean_Conc,
                                    color = ifelse((mean_Conc < lim)|is.na(lim), "black", "red")), 
                                alpha = .8, size = 5, shape = 21, stroke = 1)+ theme_minimal() + 
                     theme(legend.position = "none") + 
                     ylab("Latitude")+xlab("Longitude") +
                     scale_color_manual(values = c("red" = "red", "black" = "black")) +
                     scale_fill_continuous(high = "#132B43", low = "#56B1F7"),
                   
                   sub_data %>%
                     ggplot(aes(x = station, y = mean_Conc)) + 
                     geom_col(aes(fill = mean_Conc), width = .6, color = "black", size = .5) +  theme_minimal() +
                     coord_flip() + xlab("") +ylab(mocis_get_unit(var = vars[idx], gen = "CLUP")) + 
                     theme(legend.position = "none", axis.text.x = element_text(angle = 45)) +
                     scale_fill_continuous(high = "#132B43", low = "#56B1F7")
    )
    ggsave(p, filename = paste0("figs/", vars[idx], "_mean_map.pdf"), width = 8, height = 5)
    ggsave(p, filename = paste0("figs/", vars[idx], "_mean_map.png"), width = 8, height = 5)
  }
  
  
}

plot_rate_maps <- function(mocis_object, report_year = 2019){
  data <- mocis_tibble(mocis_object) %>% filter(gen == "CLUP", loc != "ORFJ") %>% #2020 addition from Martin
    mutate(prc_all_lod = map_dbl(aggdata, ~mean(.x$all.lod))) %>%
    filter(prc_all_lod < .5)
  ##data <- mocis_tibble(mocis_object) %>% filter(gen == "CLUP", loc != "ORFJ")
  if (nrow(data) == 0) return()
  world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
    st_crop(c(xmin=8, xmax=25, ymin=54, ymax=70))
  #lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf") %>% 
  #  st_crop(c(xmin=8, xmax=25, ymin=54, ymax=70))
  plot_data <- data %>% mutate(
    EKOO = map_dbl(aggdata, ~.x[["EKOO"]][1]),
    NKOO = map_dbl(aggdata, ~.x[["NKOO"]][1]),
    lat = ne2latlong(NKOO, EKOO)[, 1],
    long = ne2latlong(NKOO, EKOO)[, 2],
    long = if_else(loc %in% c("ANGV", "UTLV"), long + 0.5, long),
    lim = map_dbl(limit, ~.x$tv[1]),
    slope = map_dbl(linmod10,  ~ifelse(is.null(.x[["slope"]]), NA, .x[["slope"]])),
    mean_Conc = map_dbl(tidyaggdata, ~filter(.x, YEAR %in% (report_year-3):report_year) %>% 
                          summarise(mean_Conc = exp(mean(log(value)))) %>% 
                          pull(mean_Conc)),
    station = mocis_name_station(loc),
    station = fct_reorder(station, -(-1)^(long > 13.6)*lat)) %>% 
    select(gen, var, station, slope, lim, mean_Conc, lat, long) %>% 
    filter(!is.na(slope)) %>% 
    arrange(desc(station)) %>% 
    group_by(var) %>% 
    nest()
  vars <- plot_data[["var"]]
  for (idx in (1:length(vars))){
    sub_data <- plot_data %>% pull(data) %>% .[[idx]]
    p <- plot_grid(ggplot(world_map) + 
                     geom_sf(color = "white", fill = "grey") + 
                     #geom_sf(data = lakes, fill = "white", color = "grey", size = 0) +
                     annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
                     annotation_scale() +
                     coord_sf(xlim = c(8, 25), ylim = c(54, 70), expand = FALSE) +
                     geom_point(data = sub_data, 
                                aes(x = long, y = lat, 
                                    fill = slope), 
                                alpha = .8, size = 5, shape = 21, stroke = 1)+ theme_minimal() + 
                     scale_fill_gradient2(high = muted("red"), mid = "white", low = muted("blue")) +
                     theme(legend.position = "none") + 
                     ylab("Latitude")+xlab("Longitude"),
                   
                   sub_data %>%
                     ggplot(aes(x = station, y = slope)) + 
                     geom_col(aes(fill = slope), width = .6, color = "black", size = .5) +  theme_minimal() +
                     coord_flip()  +xlab("") +ylab("Yearly relative change, %") + 
                     theme(legend.position = "none", axis.text.x = element_text(angle = 45)) +
                     scale_fill_gradient2(high = muted("red"), mid = "white", low = muted("blue")))
    
    
    ggsave(p, filename = paste0("figs/", vars[idx], "_rate_map.pdf"), width = 8, height = 5)
    ggsave(p, filename = paste0("figs/", vars[idx], "_rate_map.png"), width = 8, height = 5)
  }
}




mocis_map <- function(var, data){
  plot_grid(ggplot(world_map) + 
              geom_sf(color = "white", fill = "grey") + 
              geom_sf(data = lakes, fill = "white", color = "grey", size = 0) +
              annotation_north_arrow(location = "tl", height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
              annotation_scale() +
              coord_sf(xlim = c(8, 25), ylim = c(54, 70), expand = FALSE) +
              geom_point(data = data, 
                         aes(x = long, y = lat, 
                             fill = log(mean_Conc),
                             color = ifelse((mean_Conc < lim)|is.na(lim), "black", "red")), 
                         alpha = .8, size = 5, shape = 21, stroke = 1) + 
              theme(legend.position = "none") + 
              ylab("Latitude")+xlab("Longitude") +
              scale_color_manual(values = c("red" = "red", "black" = "black")),
            data %>%
              ggplot(aes(x = station, y = mean_Conc)) + 
              geom_col(aes(fill = log(mean_Conc)), color = "black") +  
              coord_flip() + xlab("") +ylab(mocis_get_unit(var = "AHCH", gen = "CLUP")) + 
              theme(legend.position = "none", axis.text.x = element_text(angle = 45)))
}


appendix_section <- function(section_contents){
  var <- section_contents$var[1]
  var_name <- section_contents$var_name[1]
  maps <- section_contents$maps[1]
  cat(paste("\n\n##", var_name, "\n\n"))
  if (sum(section_contents$gen == "CLUP") > 0){
    knitr::knit_expand(file = "rmd/CLUPsection.Rmd") %>% cat(sep = "\n")
  }
  if (sum(section_contents$gen %in% c("GADU", "ZOAR", "PERC")) > 0){
    knitr::knit_expand(file = "rmd/FISHsection.Rmd") %>% cat(sep = "\n")
  }
  if (sum(section_contents$gen == "MYTI") > 0){
    knitr::knit_expand(file = "rmd/MUSSELsection.Rmd") %>% cat(sep = "\n")
  }
  if (sum(section_contents$gen %in% c("SIGR", "HAEM", "STER")) > 0){
    knitr::knit_expand(file = "rmd/EGGsection.Rmd") %>% cat(sep = "\n")
  }
}

appendix_chapter <- function(chapter_contents){
  contents <- split(chapter_contents, chapter_contents$var_order)
  cat(paste("\n\n#", chapter_contents$group[1], "\n\n"))
  map(contents, appendix_section)
}

appendix <- function(contents, file_name = "appendix.Rmd", exclude_maps = NULL){
  contents <- mutate(contents, maps = !(var %in% exclude_maps)) %>% 
    split(contents$group_order)
  con <- file(file_name, open = "wt", encoding = "UTF-8")
  sink(con)
  knitr::knit_expand(file = "rmd/header.Rmd") %>% cat(sep = "\n")
  map(contents, appendix_chapter)
  cat("\n\n# References\n\n")
  sink()
  close(con)
}


variables <- function(){
  tibble(files = c("bio", "bio", "bio", "bio", "bio", "met", 
                   "met", "met", "met", "met", "met", "met", "met", "met", "met", 
                   "met", "clc", "clc", "clc", "clc", "clc", "clc", "clc", "clc", 
                   "clc", "clc", "clc", "clc", "clc", "clc", "clc", "dx", "dx", 
                   "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", 
                   "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", "dx", 
                   "dx", "dx", "dx", "dx", "dx", "dx", "dx", "bfr", "bfr", "bfr", 
                   "bfr", "bfr", "bfr", "bfr", "pfc", "pfc", "pfc", "pfc", "pfc", 
                   "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", 
                   "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", "pfc", "pah", "pah", "pah", 
                   "pah", "pah", "pah", "pah", "pah", "pah", "pah", "pah", "pah", 
                   "pah", "pah", "pah", "pah", "to", "to", "to", "to", "to", "to", 
                   "to", "to"), 
         var = c("ALDR", "TOTV", "TOTL", "FPRC", "KOND", 
                 "AG", "AL", "AS", "CD", "CR", "CU", "HG", "NI", "PB", "SE", "ZN", 
                 "DDT", "DDD", "DDE", "AHCH", "BHCH", "LINDA", "HCB", "PCBSUM", 
                 "CB28", "CB52", "CB101", "CB118", "CB118", "CB153", "CB180", 
                 "CB77", "CB81", "CB126", "CB169", "CB123", "CB114", "CB105", 
                 "CB167", "CB156", "CB157", "CB189", "CBEQV", "TCDF", "PECDF1", 
                 "PECDF2", "HXCDF1", "HXCDF2", "HXCDF3", "HXCDF4", "HPCDF1", "HPCDF2", 
                 "OCDF", "TCDD", "PECDD", "HXCDD1", "HXCDD2", "HXCDD3", "HPCDD", 
                 "OCDD", "TCDDEQV", "TCDDEQVW", "BDE47", "BDE99", "BDE100", "BDE153", 
                 "BDE154", "HBCD", "BDE28", "PFHXA", "PFHPA", "PFOA", "PFNA", 
                 "PFDA", "PFUNDA", "PFDODA", "PFTRDA", "PFTEDA", "PFPEDA", "LPFBS", 
                 "PFHXS", "LPFOS", "BPFOS", "PFOS", "LPFDS", "BPFDS", "PFDS", 
                 "LFOSA", "BFOSA", "FOSA", "SUMPAH", "NAP", "ACNE", "FLE", "PA", 
                 "ANT", "FLU", "PYR", "BAA", "CHR", "BBF", "BKF", "BAP", "BGHIP", 
                 "DBAHA", "ICDP", "MBT", "DIBT", "TBT", "MPT", "DIPT", "TPT", 
                 "MOT", "DIOT"), 
         group = c("Biological data", "Biological data", 
                   "Biological data", "Biological data", "Biological data", "Metals", 
                   "Metals", "Metals", "Metals", "Metals", "Metals", "Metals", "Metals", 
                   "Metals", "Metals", "Metals", "Chlorinated Pesticides", "Chlorinated Pesticides", 
                   "Chlorinated Pesticides", "Chlorinated Pesticides", "Chlorinated Pesticides", 
                   "Chlorinated Pesticides", "Chlorinated Pesticides", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", "Polychlorinated Biphenyls", 
                   "Polychlorinated Biphenyls", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Polychlorinated Dibenzodioxins /-furans", "Polychlorinated Dibenzodioxins /-furans", 
                   "Brominated Flame Retardants", "Brominated Flame Retardants", 
                   "Brominated Flame Retardants", "Brominated Flame Retardants", 
                   "Brominated Flame Retardants", "Brominated Flame Retardants", 
                   "Brominated Flame Retardants", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", "Perfluoroalkyl Substances", 
                   "Perfluoroalkyl Substances", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Polycyclic Aromatic Hydrocarbons", 
                   "Polycyclic Aromatic Hydrocarbons", "Organotin compounds", "Organotin compounds", 
                   "Organotin compounds", "Organotin compounds", "Organotin compounds", 
                   "Organotin compounds", "Organotin compounds", "Organotin compounds"
         ), 
         group_order = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                         2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
                         4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
                         5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 
                         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 
                         8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9), 
         var_order = c(1, 
                       2, 3, 4, 6, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 5, 
                       6, 7, 1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
                       17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                       16, 17, 18, 19, 2, 3, 4, 5, 6, 7, 1, 1, 2, 3, 4, 5, 6, 7, 8, 
                       9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 1, 2, 3, 4, 
                       5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 
                       7, 8), 
         var_name = c("Age, years", "Total weight, g", "Total length, cm", 
                      "Fat, %", "Fultons condition factor", "Ag", "Al", "As", "Cd", 
                      "Cr", "Cu", "Hg", "Ni", "Pb", "Se", "Zn", "DDT", "DDD", "DDE", 
                      "alpha-HCH", "beta-HCH", "Lindane", "HCB", "Sum PCBs", "PCB-28", 
                      "PCB-52", "PCB-101", "PCB-118", "PCB-118", "PCB-153", "PCB-180", 
                      "PCB-77", "PCB-81", "PCB-126", "PCB-169", "PCB-123", "PCB-114", 
                      "PCB-105", "PCB-167", "PCB-156", "PCB-157", "PCB-189", "TCDD-equivalents (Sum WHO-PCB-TEQ)", 
                      "2,3,7,8-TCDF", "1,2,3,7,8-PeCDF", "2,3,4,7,8-PeCDF", "1,2,3,4,7,8-HxCDF", 
                      "1,2,3,6,7,8-HxCDF", "2,3,4,6,7,8-HxCDF", "1,2,3,7,8,9-HxCDF", 
                      "1,2,3,4,6,7,8-HpCDF", "1,2,3,4,7,8,9-HpCDF", "1,2,3,4,6,7,8,9-OCDF", 
                      "2,3,7,8-TCDD", "1,2,3,7,8-PeCDD", "1,2,3,4,7,8-HxCDD", "1,2,3,6,7,8-HxCDD", 
                      "1,2,3,7,8,9-HxCDD", "1,2,3,4,6,7,8-HpCDD", "1,2,3,4,6,7,8,9-OCDD", 
                      "TCDD-equivalents (Sum WHO-PCDD/F-TEQ, lw)", "TCDD-equivalents (Sum WHO-PCDD/F-TEQ, ww)", "BDE-47", "BDE-99", 
                      "BDE-100", "BDE-153", "BDE-154", "HBCDD", "BDE-28", "PFHxA", 
                      "PFHpA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFDoDA", "PFTrDA", 
                      "PFTeTA", "PFPeDA", "lin-PFBS", "PFHxS", "lin-PFOS", "br-PFOS", "PFOS", 
                      "lin-PFDS", "br-PFDS", "PFDS", "lin-FOSA", "br-FOSA", "FOSA", 
                      "Sum PAHs", "Naphtalene", "Acenapthene", "Flourene", "Phenantrene", 
                      "Anthracene", "Flouranthene", "Pyrene", "Benso (a) anthracene", 
                      "Chrysene", "Benso (b) fluoranthene", "Benso (k) fluoranthene", 
                      "Benso (a) pyrene", "Benso (g, h, i) perylene", "Dibenso (a , h) anthracene", 
                      "Indeno (1, 2, 3) pyrene", "Monobytultin", "Dibutyltin", "Tributyltin (TBT)", 
                      "Monophenyltin", "Diphenyltin", "Triphenyltin", "Monooctyltin", 
                      "Dioctyltin"))
}
