read_lab_files_loq <- function(lab_path){
  files <- dir(lab_path)
  files <- subset(files, str_ends(files, "xlsm") & (!str_detect(files, "~")))
  file_paths <- paste0(lab_path, files) 
  # Read files
  # Convert to long format and bind rows
  # Convert back to wide format
  lab_data <- map(file_paths, read_lab_file_loq) %>% 
    map(~pivot_longer(.x, -c("ACCNR", "LABID", "GENUS", "STATION", "TISSUE"), names_to = "VAR", values_to = "VAL")) %>% 
    map_df(bind_rows) %>% 
    select(ACCNR, STATION, VAR, VAL) %>% 
    distinct() %>% 
    pivot_wider(id_cols = c("ACCNR", "STATION"), names_from = VAR, values_from = VAL)
}
read_lab_file_loq <- function(path){
  # Reads an excel lab-protocol, function loqify decides what to do with values
  # reported as negative
  results <- read_excel(path, sheet = "results", skip = 1, na = c("-99.99", "N/A")) %>% 
    rename(ACCNR = ...1, LABID = ...2, GENUS = ...3, STATION = ...4, TISSUE = ...5) %>%
    filter(!is.na(ACCNR)) %>% 
    mutate(LABID = as.character(LABID)) %>% 
    mutate_if(is.character, ~str_replace(.x, "<", "-")) %>% # Check if "<" rather than "-" is used for LOQ
    mutate_at(-(1:5), as.numeric)
  loq <- read_excel(path, sheet = "LOQ", skip = 1, na = c("-99.99", "N/A")) %>% 
    rename(ACCNR = ...1)
  has_LOQ <- setdiff(intersect(names(results), names(loq)), "ACCNR")
  df <- left_join(results, loq, by = "ACCNR", suffix = c("", ".LOQ"))
  for (v in has_LOQ){
    results[[v]] <- loqify(df[[v]], df[[paste0(v, ".LOQ")]])
  }
  results
}
loqify <- function(value, loq){
  loq <- as.numeric(loq)
  transform <- (value < 0)
  transform[is.na(transform)] <- FALSE
  value[transform] <- pmin(value[transform], -loq[transform], na.rm = TRUE)
  value
}

# # Functions replaced by loq-versions
# read_lab_file <- function(path, sheet = "results"){
#   read_excel(path, sheet = sheet, skip = 1, na = c("-99.99", "N/A")) %>% 
#     rename(ACCNR = ...1, LABID = ...2, GENUS = ...3, STATION = ...4, TISSUE = ...5) %>%
#     filter(!is.na(ACCNR)) %>% 
#     mutate(LABID = as.character(LABID)) %>% 
#     mutate_if(is.character, ~str_replace(.x, "<", "-")) %>% # Check if "<" rather than "-" is used for LOQ
#     mutate_at(-(1:5), as.numeric)
# }
# 
# read_lab_files <- function(lab_path){
#   file_paths <- paste0(lab_path, dir(lab_path)[dir(lab_path) %>% str_ends("xlsm")]) # Files ending with xlsm in lab_path
#   # Read files
#   # Convert to long format and join by binding rows
#   # Convert back to wide format
#   lab_data <- map(file_paths, read_lab_file) %>% 
#     map(~pivot_longer(.x, -c("ACCNR", "LABID", "GENUS", "STATION", "TISSUE"), names_to = "VAR", values_to = "VAL")) %>% 
#     map_df(bind_rows) %>% 
#     select(ACCNR, STATION, VAR, VAL) %>% 
#     distinct() %>% 
#     pivot_wider(id_cols = c("ACCNR", "STATION"), names_from = VAR, values_from = VAL)
# }

unpool <- function(ACCNR){
  # Given ACCNR, splits into a table of individual ACCNR and corresponding pooled ACCNR
  if (str_detect(ACCNR, "-")){ #homogenat
    first <- str_sub(ACCNR, 7, 11) %>% as.integer()
    last <- str_sub(ACCNR, 13, 17) %>% as.integer()
    all <- first:last %>% as.character()
    ACCNR_IND <- paste0(str_sub(ACCNR, 1, 6), strrep("0", 5- str_length(all)), all)
    tibble(ACCNR_POOL = ACCNR, ACCNR = ACCNR_IND)
  } 
  else
  {
    tibble(ACCNR_POOL = ACCNR, ACCNR = ACCNR)
  }
}

trans_accnr <- function(ACCNR_ESBASE){
  # Translates from Esbase ACCNR-format
  if (str_sub(ACCNR_ESBASE, 1, 1) %in% c("A", "C", "D", "G", "P")){return(ACCNR_ESBASE)}
  series_digit <- str_sub(ACCNR_ESBASE, 1, 1)
  series_letter <- case_when(
    series_digit == "1"~"A",
    series_digit == "3"~"C",
    series_digit == "4"~"D",
    series_digit == "7"~"G",
    series_digit == "9"~"P"
  )
  year <- ifelse(as.numeric(str_sub(ACCNR_ESBASE, 2, 2)) == "9", 
                 paste0("1", str_sub(ACCNR_ESBASE, 2, 4)),
                 paste0("2", str_sub(ACCNR_ESBASE, 2, 4)))
  number <- str_sub(ACCNR_ESBASE, 5, 9)
  paste0(series_letter, year, "/", number)
}

mean_or_na <- function(x){
  # returns NA rather than NaN if all elements of x are NA
  ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))
}

read_historical_file <- function(file, na_list = c("-9", "-9.9", "")){
  n <- read_excel(file, n_max = 2) %>% ncol()
  read_excel(file, col_types = c(rep("text", 4), rep("numeric", n - 4)), na = na_list) %>% select_if(~!all(is.na(.x)))
}

append_cols <- function(x, y){
  # Adds the columns in y that are not already in x
  dy <- select(y, setdiff(names(y), names(x)))
  bind_cols(x, dy)
}

filter_species <- function(data, x){
  if (str_detect(x, "uria"))
  { filter(data, GENUS %in% c("SIGR", "STER", "HAEM"))}
  else
  { filter(data, GENUS %in% c("CLUP", "GADU", "MYTI", "PERC", "ZOAR"))}
}

stationsreg <- function(){
  # This data.frame was obtained by
  # stationsreg <- read_xlsx("Y:/MFO_Privat/MG/Miljögiftsövervakning/SGU/NRM och SGU koder_version 20191209.xlsx", 
  #                          sheet = "Stationsregister_marina", skip = 1) %>% 
  #   mutate(`ESBase-id` = ifelse(LOC == "49G9", paste0(`ESBase-id`, ", 1293"), `ESBase-id`)) %>% 
  #   separate_rows(`ESBase-id`, sep = ",") %>% 
  #   mutate(`ESBase-id` = as.numeric(`ESBase-id`)) %>% 
  #   select(LOC, NKOO, EKOO, `ESBase-id`, GENUS) %>% 
  #   na.omit() %>% 
  #   distinct() %>% 
  #   mutate(Season = ifelse(LOC %in% c("ANGV", "UTLV"), "V", "H"), # Spring localities
  #          Season = ifelse(GENUS %in% c("HAEM", "STER", "SIGR"), "V", Season) # Eggs always collected in spring
  #   )
  tibble(LOC = c("40G7", "40G7", "46H0", "46H0", "49G9", 
                     "49G9", "49G9", "49G9", "49G9", "49G9", "51G9", "51G9", "51G9", 
                     "ABBE", "ABBE", "ABBE", "ANGV", "BYXE", "BYXE", "FLAD", "FLAD", 
                     "FLAD", "GAFJ", "GAFJ", "HABU", "HABU", "HABU", "HAFJ", "HAFJ", 
                     "HOLM", "KIFJ", "KULL", "LAFJ", "LAFJ", "LAGN", "LAND", "RAFJ", 
                     "UTLA", "UTLA", "VADO", "ANGK", "ORFJ", "ORFJ", "HOLM", "KVFJ", 
                     "UTLV", "UTLV", "KVFJ", "FJBA", "FJBA", "FJBA", "SEGO", "SEGO", 
                     "SEGO", "SEGO", "SEGO", "SEGO", "SEGO", "SEGO", "SEGO", "SEGO", 
                     "SEGO", "FLAD", "FLAD", "FLAD", "NIDI", "FJBA", "FJBA", "FJBA", 
                     "KVFJ", "TJRN", "TJRN", "STKA", "STKA"), 
             NKOO = c(6181395, 6181395, 
                      6523708, 6523708, 6686566, 6686566, 6686566, 6686566, 6686566, 
                      6686566, 6798326, 6798326, 6798326, 6134000, 6134000, 6134000, 
                      6715100, 6365800, 6365800, 6348600, 6348600, 6348600, 7005100, 
                      7005100, 6181700, 6181700, 6181700, 7294000, 7294000, 7073600, 
                      7204900, 6249400, 6852200, 6852200, 6593400, 6510000, 7310900, 
                      6208830, 6208830, 6502000, 6715100, 7039900, 7039900, 7073600, 
                      6434800, 6208830, 6208830, 6434800, 6510100, 6510100, 6510100, 
                      6294700, 6294700, 6294700, 6294700, 6294700, 6294700, 6294700, 
                      6294700, 6294700, 6294700, 6294700, 6348600, 6348600, 6348600, 
                      6368600, 6510100, 6510100, 6510100, 6434800, 6539496, 6539496, 
                      6352800, 6352800), 
             EKOO = c(1606416, 1606416, 1771651, 1771651, 
                      1696248, 1696248, 1696248, 1696248, 1696248, 1696248, 1698277, 
                      1698277, 1698277, 1360700, 1360700, 1360700, 1629400, 1571500, 
                      1571500, 1258800, 1258800, 1258800, 1642800, 1642800, 1404600, 
                      1404600, 1404600, 1825900, 1825900, 1750800, 1759200, 1288200, 
                      1587100, 1587100, 1660100, 1627500, 1802700, 1501600, 1501600, 
                      1218300, 1629400, 1679300, 1679300, 1750800, 1556700, 1501600, 
                      1501600, 1556700, 1236200, 1236200, 1236200, 1664600, 1664600, 
                      1664600, 1664600, 1664600, 1664600, 1664600, 1664600, 1664600, 
                      1664600, 1664600, 1258800, 1258800, 1258800, 1265100, 1236200, 
                      1236200, 1236200, 1556700, 1230900, 1230900, 1631500, 1631500
             ), 
             `ESBase-id` = c(806, 1317, 790, 1332, 783, 784, 2218, 1290, 
                             291, 1293, 1276, 778, 290, 349, 1487, 613, 350, 399, 351, 353, 
                             1090, 517, 355, 2110, 397, 1061, 312, 356, 901, 348, 357, 758, 
                             360, 2107, 358, 337, 1166, 368, 1027, 254, 350, 362, 2108, 348, 
                             255, 368, 1027, 255, 336, 352, 254, 265, 580, 578, 831, 159, 
                             898, 1552, 367, 1385, 1897, 854, 353, 1090, 517, 376, 336, 352, 
                             254, 255, 632, 632, 403, 252), 
             GENUS = c("CLUP", "CLUP", "CLUP", 
                       "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", 
                       "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", 
                       "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", 
                       "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", 
                       "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "CLUP", "PERC", "PERC", 
                       "PERC", "PERC", "CLUP", "CLUP", "ZOAR", "ZOAR", "ZOAR", "ZOAR", 
                       "GADU", "GADU", "GADU", "GADU", "GADU", "GADU", "GADU", "GADU", 
                       "GADU", "GADU", "GADU", "GADU", "GADU", "GADU", "MYTI", "MYTI", 
                       "MYTI", "MYTI", "MYTI", "STER", "HAEM", "SIGR", "SIGR"), 
             Season = c("H", 
                        "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
                        "H", "H", "V", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
                        "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
                        "H", "H", "H", "H", "H", "V", "V", "H", "H", "H", "H", "H", "H", 
                        "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
                        "H", "H", "H", "H", "V", "V", "V", "V"))
}