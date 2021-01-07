library(tidyverse)
library(openxlsx)
library(readxl)
library(MoCiS)
library(MoCiS.tools)
library(ggforce)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(scales)

setwd("Y:/MFO_Privat/MG/Miljögiftsövervakning/MCOM/MCOM2020/R/MCOM20_report/")

report_year <- 2019
# Backgrounds for station map
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  st_crop(c(xmin=1, xmax=36, ymin=53, ymax=70))
#lakes <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf") %>% 
#  st_crop(c(xmin=1, xmax=35, ymin=54, ymax=70))
# Utility functions
source("functions/report_functions.R")
# Aggregate data from MoCiS for figures and tables
load("../MCOM20_MoCiS/mocis_all.Rdata")
mocis_data <- mocis_all %>% 
  map_df(mocis_tibble) %>% 
  unique() %>% 
  select(loc, var, gen, tidyaggdata, linmod10, limit)
# Station data for map and basin
station_table <- read_excel("stations.xlsx")

##
##  Station map
##

map <- ggplot(world_map) + 
  geom_sf(color = "white", aes(fill = (name == "Sweden"))) + 
  geom_sf(color = "grey90", size = 0) +
  #geom_sf(data = lakes, fill = "white", color = "grey90", size = 0) +
  theme_minimal() +
  coord_sf(xlim = c(1, 36), ylim = c(53, 70), expand = FALSE) +
  geom_text_repel(data = filter(station_table, basin != "West Coast"), aes(x = long, y = lat, label = paste0(station, " (", collected, ")")), 
                  size = 2.5, nudge_x = 25- filter(station_table, basin != "West Coast")$long, direction = "y", hjust = 0,
                  segment.size  = 0.5,
                  segment.color = "grey50")  +
  geom_text_repel(data = filter(station_table, basin == "West Coast"), aes(x = long, y = lat, label = paste0(station, " (", collected, ")")), 
                  size = 2.5, nudge_x = 2 - filter(station_table, basin == "West Coast")$long, direction = "y", hjust = 1,
                  segment.size  = 0.5,
                  segment.color = "grey50")  +
  geom_point(data = station_table, aes(x = long, y = lat, color = basin), size = 2) +
  theme(legend.position = c(0.18, .87), legend.title = element_blank()) + xlab("") + ylab("") +
  scale_fill_manual(values = c("grey90", "grey"), guide = FALSE)

ggsave(map, filename = "report_figs/map.png", dpi = 600, width = 5, height = 5)

##
## Basin plots
##

# Contaminants to be included in basin-level analysis
basin_contaminants <- c("BDE47", "CB118", "CB153", "CD", "CU", "DDE", "FOSA", "HBCD", "HCB", "HG", "PB", "PFOS", "PFUNDA", "TCDD", "TCDDEQV", "TCDF")

# Data wrangling
basin_data <- filter(mocis_data, var %in% basin_contaminants,
                     gen == "CLUP") %>% 
  select(tidyaggdata) %>% 
  unnest(cols = tidyaggdata) %>% 
  filter(YEAR > 2008) %>% 
  left_join(station_table, by = c("LOC" = "loc")) %>% 
  select(YEAR, var, value, LOC, station, basin) %>% 
  filter(!(LOC %in% c("ANGV", "UTLV"))) # Exclude spring localities

# Model fitting
basin_fits <- basin_data %>% 
  group_by(var, basin) %>% 
  nest() %>% 
  mutate(lm_fit = map(data, ~lm(log(value) ~ LOC + YEAR, data = .x)),
         lm_fit_full = map(data, ~lm(log(value) ~ LOC * YEAR, data = .x)),
         p_vs_full = map2_dbl(lm_fit, lm_fit_full, ~anova(.x, .y)$"Pr(>F)"[2]),
         lm_summary = map(lm_fit, summary),
         lm_intercept = map_dbl(lm_fit, common_intercept),
         lm_slope = map_dbl(lm_fit, ~coef(.x)["YEAR"]),
         lm_slope_se = map_dbl(lm_fit, ~vcov(.x)["YEAR", "YEAR"] %>% sqrt()),
         lm_slope_CI = map(lm_fit, ~confint(.x)["YEAR", ]),
         lm_slope_upper = map_dbl(lm_slope_CI, ~.x[2]),
         lm_slope_lower = map_dbl(lm_slope_CI, ~.x[1]),
         lm_predict = pmap(list(lm_intercept, lm_slope, lm_slope_lower, lm_slope_upper), 
                           ~tibble(YEAR = seq(2010, 2019, by = 0.5), CONC = exp(..1 + ..2 * YEAR),
                                   UPPER = exp(..1 + ..2 * YEAR + (..4 - ..2) * abs(YEAR - mean(YEAR))),
                                   LOWER = exp(..1 + ..2 * YEAR + (..3 - ..2) * abs(YEAR - mean(YEAR))))),
         table_ci = pmap_chr(list(slope2rate(lm_slope), slope2rate(lm_slope_lower), slope2rate(lm_slope_upper), p_vs_full), 
                             ~if_else(..4 > .05, print_ci(..1, ..2, ..3), as.character(round(..1, 2))))
  ) %>% 
  ungroup() %>% 
  left_join(read_excel("variables.xlsx"), by = "var") %>% 
  mutate(var_name = fct_reorder(var_name, group_order * 1000 + var_order)) %>% 
  select(var_name, var, group, basin, lm_intercept, lm_slope, lm_slope_CI, lm_slope_lower, lm_slope_upper, p_vs_full, lm_predict, table_ci, data)

# Group-wise plots of concentrations
basin_plot_data <- unnest(basin_fits, lm_predict)
basin_met <- basin_plot(basin_plot_data, grp = "Metals", unit = bquote("Concentration: Hg - other metals,"~n*g/g - ~mu*g/g))
basin_pd <- basin_plot(basin_plot_data, grp = "Polychlorinated Dibenzodioxins /-furans", unit = bquote("Concentration,"~p*g/g))
basin_pb <- basin_plot(basin_plot_data, grp = "Polychlorinated Biphenyls", unit = bquote("Concentration,"~mu*g/g))
basin_cp <- basin_plot(basin_plot_data, grp = "Chlorinated Pesticides", unit = bquote("Concentration,"~mu*g/g))
basin_bfr <- basin_plot(basin_plot_data, grp = "Brominated Flame Retardants", unit = bquote("Concentration,"~n*g/g))
basin_ps <- basin_plot(basin_plot_data, grp = "Perfluoroalkyl Substances", unit = bquote("Concentration,"~n*g/g))

h <-  4
w <-  7
ggsave(basin_met, filename  = "report_figs/basin_met.png", dpi = 600, width = w, height = h * 1.8)
ggsave(basin_pd, filename  = "report_figs/basin_pd.png", dpi = 600, width = w, height = h * 1.8)
ggsave(basin_pb, filename  = "report_figs/basin_pb.png", dpi = 600, width = w, height = h)
ggsave(basin_cp, filename  = "report_figs/basin_cp.png", dpi = 600, width = w, height = h)
ggsave(basin_bfr, filename  = "report_figs/basin_bfr.png", dpi = 600, width = w, height = h)
ggsave(basin_ps, filename  = "report_figs/basin_ps.png", dpi = 600, width = w, height = h * 1.8)

ggsave(basin_plot(basin_plot_data, "Perfluoroalkyl Substances", cols = 3, unit = bquote("Concentration,"~n*g/g)), 
       filename = "report_figs/basin_ps3x1.png", width = w * 1.3, height = h)
ggsave(basin_plot(basin_plot_data, "Polychlorinated Dibenzodioxins /-furans", cols = 3, unit = bquote("Concentration,"~p*g/g)), 
       filename = "report_figs/basin_pd3x1.png", width = w * 1.3, height = h)


##
## Basin table
##

basin_fits %>% select(var_name, basin,  table_ci) %>% 
  distinct() %>% 
  pivot_wider(names_from = basin, values_from = table_ci) %>% 
  rename(Contaminant = var_name) %>% 
  write.xlsx(file = "report_figs/basin_table.xlsx")

##
## Heatmaps
##

plot_data <- mocis_data %>% 
  filter(!(var %in% c("CR", "AHCH", "LINDA", "PFHXA", "PFHPA", "LPFBS", "LPFDS", "PFDS"))) %>% # Many LOQs...
  filter(!((loc == "HOLM") & (gen == "ZOAR"))) %>% 
  left_join(read_excel("variables.xlsx"), by = "var") %>% 
  mutate(null.list = map_lgl(linmod10, is.null)) %>% 
  filter(null.list == FALSE) %>% 
  mutate(mean_Conc = map_dbl(tidyaggdata, ~filter(.x, YEAR %in% (report_year-3):report_year) %>% 
                               summarise(mean_Conc = exp(mean(log(value)))) %>% 
                               pull(mean_Conc)),
         lim = map_dbl(limit, ~.x$tv[1]),
         EKOO = map_dbl(tidyaggdata, ~.x[["EKOO"]][1]),
         NKOO = map_dbl(tidyaggdata, ~.x[["NKOO"]][1]),
         lat = ne2latlong(NKOO, EKOO)[, 1],
         long = ne2latlong(NKOO, EKOO)[, 2],
         station_short = map_chr(tidyaggdata, ~.x[["LOC"]][1]),
         station = mocis_name_station(station_short),
         station = fct_reorder(station, -(-1)^(long > 13.6)*lat),
         slope = map_dbl(linmod10, ~ifelse(is.null(.x[["slope"]]), NA, .x[["slope"]])),
         upper = map_dbl(linmod10, ~ifelse(is.null(.x[["upper"]]), NA, .x[["upper"]])),
         lower = map_dbl(linmod10, ~ifelse(is.null(.x[["lower"]]), NA, .x[["lower"]])),
         #p_value = map_dbl(linmod10, ~summary(.x$fit) %>% .[["coefficients"]] %>% .["YEAR", "Pr(>|t|)"]),
         p_value = map_dbl(linmod10, "p"),
         stars = case_when(p_value > 0.05 ~ "", (p_value <= 0.05) & (p_value > 0.01) ~ "*", (p_value <= 0.01) & (p_value > 0.001) ~ "**", p_value <= 0.001 ~ "***"),
         gen_name = mocis_name_species(gen),
         gen_name = fct_relevel(gen_name, "Herring", "Cod", "Perch", "Eelpout", "Blue mussel", "Guillemot", "Common tern", "Eurasian Oystercatcher")
  ) %>% 
  arrange(desc(station)) %>% 
  select(var, var_name, gen, gen_name, station, group, mean_Conc, slope, upper, lower, stars, lim, lat, long, var_order, group_order) %>% 
  group_by(gen, station, group) %>% 
  mutate(all.na = all(is.na(slope))) %>% 
  ungroup() %>% 
  filter(!all.na) %>% 
  mutate(var_name = fct_reorder(var_name, group_order*100 + var_order))

heat_bio <- mocis_heatmap(plot_data %>% filter(gen != "MYTI"), grp = "Biological data")
heat_met <- mocis_heatmap(plot_data, grp = "Metals")
heat_pd <- mocis_heatmap(plot_data, grp = "Polychlorinated Dibenzodioxins /-furans")
heat_pb <- mocis_heatmap(plot_data, grp = "Polychlorinated Biphenyls")
heat_cp <- mocis_heatmap(plot_data, grp = "Chlorinated Pesticides")
heat_bfr <- mocis_heatmap(plot_data, grp = "Brominated Flame Retardants")
heat_ps <- mocis_heatmap(plot_data, grp = "Perfluoroalkyl Substances")
heat_pah <- mocis_heatmap(plot_data, grp = "Polycyclic Aromatic Hydrocarbons")
heat_to <- mocis_heatmap(plot_data, grp = "Organotin compounds")


ggsave(heat_bio, filename = "report_figs/heat_bio.png", width = 8, height = 9, dpi = 600)
ggsave(heat_met, filename = "report_figs/heat_met.png", width = 8, height = 9, dpi = 600)
ggsave(heat_pd, filename = "report_figs/heat_pd.png", width = 8, height = 9, dpi = 600)
ggsave(heat_pb, filename = "report_figs/heat_pb.png", width = 8, height = 11, dpi = 600)
ggsave(heat_cp, filename = "report_figs/heat_cp.png", width = 8, height = 9, dpi = 600)
ggsave(heat_bfr, filename = "report_figs/heat_bfr.png", width = 8, height = 9, dpi = 600)
ggsave(heat_ps, filename = "report_figs/heat_ps.png", width = 8, height = 7, dpi = 600)
ggsave(heat_pah, filename = "report_figs/heat_pah.png", width = 8, height = 4, dpi = 600)
ggsave(heat_to, filename = "report_figs/heat_to.png", width = 8, height = 4, dpi = 600)

#
# Isotopes
#


iso_data <- read_csv("../MCOM20_data/clean_data/full_data_2019.csv", guess_max = 100000) %>% 
  select(YWD, ACCNR, GENUS, LOC, D13CUCD, CUCD, D15NUCD,NUCD) %>% 
  mutate(d13C = D13CUCD - 3.32  + 0.99 * CUCD / NUCD,
         d15N = D15NUCD,
         YEAR = str_sub(ACCNR, 2, 3),
         #YWD = 2 %>%
         YEAR = case_when(
         YWD == 18 ~ "18",
         YWD == 19 ~ "19",
         YWD > 30 ~ YEAR)) %>%
  na.omit() %>% 
  left_join(station_table, by = c("LOC"="loc")) %>% 
  mutate(Species = mocis_name_species(GENUS))


sia_plot <- iso_data  %>% 
  group_by(LOC, YEAR, basin, Species) %>% 
  summarise(d13C = mean(d13C), d15N = mean(d15N)) %>% 
ggplot(aes(x = d13C, y = d15N, color  = Species)) + 
  geom_point(alpha = .8) +
  stat_ellipse(show.legend = FALSE) +
  facet_wrap(~basin) +
  theme(legend.position = c(.85, .25), 
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        panel.border = element_rect(colour = "lightgrey", fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        aspect.ratio = 1) +
  ylab(bquote(delta^{15}*N*","~"\211"))+
  xlab(bquote(delta^{13}*C*","~"\211")) +
  ggtitle("Stable isotopes")



sia_plot_herring <- iso_data %>% filter(GENUS == "CLUP") %>% 
  group_by(LOC, YEAR, basin) %>% 
  summarise(d13C = mean(d13C), d15N = mean(d15N)) %>% 
  ggplot(aes(x = d13C, y = d15N, color  = basin)) +
  geom_point(alpha = .8) +
  stat_ellipse(show.legend = FALSE) +  
  theme(legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        panel.border = element_rect(colour = "lightgrey", fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        aspect.ratio = 1)+
  ylab(bquote(delta^{15}*N*","~"\211"))+
  xlab(bquote(delta^{13}*C*","~"\211")) +
  ggtitle("Stable isotopes, Herring")


ggsave(sia_plot, filename = "report_figs/sia_plot.png", width = 7, height = 6, dpi = 600)
ggsave(sia_plot_herring, filename = "report_figs/sia_plot_herring.png", width = 7, height = 6, dpi = 600)

