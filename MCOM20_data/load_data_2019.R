library(tidyverse)
library(readxl)

#setwd("Y:/MFO_Privat/MG/Milj?gifts?vervakning/MCOM/MCOM2020/R/MCOM20_data/")
setwd("Y:/MFO_Privat/MG/MOCIS/MCOM/R/MCOM20_data/")

source("functions/functions_load.R")
#source("Y:/MFO_Privat/MG/Milj?gifts?vervakning/MCOM/MCOM2020/R/MCOM20_data/functions/functions_load.R")

# Read historical data untill 2018. We use the combined csv-file with data generated for the 2018 data analysis (see MCOM2019 repo on github).
historical <- read.csv("data/clean_data_historical/full_data_2018.csv", check.names=FALSE, stringsAsFactors = FALSE) %>%
  mutate(OCDF = ifelse(ACCNR %in% "C88/07008-07017", NA, OCDF)) # removed outlier

# Read data from last years lab files and change unit of ClCs
ClC_names <- c("HCB", "AHCH", "BHCH", "LINDAN", "DDEK", "DDDK", "DDTK", "CB-28", "CB-52", "CB-101", "CB-118", "CB-153", "CB-138", "CB-180")
lab_data <- read_lab_files_loq(lab_path = "data/lab_raw/") %>% 
  mutate_at(ClC_names, ~.x / 1000) %>%


# In order to add biological variables to data, we need to load values for individuals 
# and generate averages for pooled samples. First we construct a table of ACCNRs and their
# corresponding pool ACCNR

pool_table <- lab_data$ACCNR %>% map_df(unpool)

# Biological variables are read and massaged, this should ideally be read directly from Esbase
# Note that species and location short names needs to be derived from latin name and locid, the 
# latter through a table "stationsreg". Two stations, Utlängan and Ängkärsklubb, have spring
# and autumn collection, these are separated by a Season variable corresponding to jan-jul or aug-dec
#
# Note that the file used here does not contain biological data on eggs
bio_data <- read_excel("data/esbase/marina19biodata.xlsx", 
                       col_types = c("text", "text", "date", "numeric", "text", "numeric", "numeric", "numeric", "numeric")) %>% 
  mutate(ACCNR = trans_accnr(ACCNR_ESBASE),
         # Get species short form
         GENUS = str_sub(SPECIES, 1, 4) %>% toupper(), # First 4 letters of latin name capitalized (?)
         GENUS = ifelse(GENUS == "URIA", "SIGR", GENUS), # Exception to above rule
         # MoCiS distinguish between eggs and fish/mussels
         TYPE = ifelse(GENUS %in% c("SIGR", "HAEM", "STER"), "uria", "hav"),
         # Get collection date in YWD-format to comply with MoCiS. Note that collection year and accession year may
         # differ, and it is not clear which one to use
         YWD = as.integer(str_sub(ACCNR, 4, 5)), # Extracts year from accession number (MoCiS does not need week/day)
         Season = ifelse((as.POSIXct(DATE) %>% format("%m")) < "08", "V", "H"),
         SEX = if_else(SEX == "male", "1", if_else(SEX == "female", "2","")),
         SEX = as.numeric(SEX),
         # The following fills in a number of values missing in the Esbase dump
         ##GENUS = ifelse(ACCNR %in% paste0("C2018/0", 7173:7255), "MYTI", GENUS), 
         Season = ifelse(ACCNR %in% paste0("C2019/0", 3689:3712), "H", Season), # FLADEN
         Season = ifelse(ACCNR %in% paste0("C2019/0", 5939:5962), "H", Season), # KINNB?CKFJ?RDEN
         ##YWD = ifelse(ACCNR %in% paste0("C2019/0", 5939:5962), "19", YWD) 
         ##Season = ifelse(ACCNR %in% paste0("C2018/0", 4029:4148), "H", Season) 
  ) %>% 
  # Join with stationsreg to get station short name, note that there seem to be no
  # standard way of giving location id to an accession in Esbase (and e.g. loc_id 254 may correspond
  # to either of Fjällbacka or Väderöarna)
  left_join(stationsreg(), by = c("LOCID_ESBASE" = "ESBase-id", "Season", "GENUS")) %>% 
  
    # Join with table of pool ACCNR and summarise
  right_join(pool_table, by = "ACCNR") %>% 
  group_by(ACCNR_POOL) %>% 
  summarise(LOC = unique(LOC),
            EKOO = unique(EKOO),
            NKOO = unique(NKOO),
            TOTL = mean_or_na(TOTL), 
            TOTV = mean_or_na(TOTV), 
            LVKT = mean_or_na(LVKT), 
            ALDR = mean_or_na(ALDR),
            GENUS = unique(GENUS),
            YWD = unique(YWD),
            NHOM = n()) %>% 
  rename(ACCNR = ACCNR_POOL)

# Bind the new data and change some names to confirm with old conventions
new_data <- left_join(bio_data, lab_data, by = "ACCNR") %>% 
  rename(LINDA = LINDAN) %>%
  #rename(LINDA = LINDAN, 'LTPRC*' = LTPRCI) %>% 
  mutate(TPRC = ifelse(is.na(TPRCI), TPRCIVL, TPRCI),
         MTPRC = MTPRCI, 'FPRC*' = FPRC) %>% 
  select(-TPRCIVL, -TPRCI) %>% 
  mutate(GRU = 1) #%>% # For MoCiS
  ##mutate(LOC = ifelse(ACCNR %in% c("C2018/04298-04309", "C2018/04310-04321"), "46H0", LOC))


# Historical data are here extracted from excel-files constructed by:
# 1) From the source .src file, extract a number of .dat-files (fixed width formatted files)
# 2) Manually, in Excel, import .dat-files and save as .xlsx
#
# The first step uses (unsupported) software written by Anders Bignert, it is unclear why a single .dat-file can not be extracted.
# The second step is highly error-prone and should be done automatically (needs start and end position of fields) if needed again.

##path_hist <- "data/historical/"

# First fish and mussels ("hav")
#
##files_hav <- dir(path_hist) %>% subset(str_detect(., "hav") & !str_detect(., "~"))
# Read historical data files to a list
##data_list_hav <- map(paste0(path_hist, files_hav), read_historical_file) 
# Join historical data files by column binding, note that this assumes that each 
# ACCNR occupies the same row in each of the Excel files (this is why "hav" and "uria" needs to be read separately)
##historical_hav <- reduce(data_list_hav, append_cols) %>% 
##  rename_all(~ifelse(.x %in% c("FPRC*", "TPRC*", "MTPRC*", "LTPRC*"), .x, str_remove(.x, "\\*"))) %>% # Remove all * but those needed by MoCiS
##  rename(LPFDS = 'LPFDS...33') %>% # Variable LPFDS occurs twice, second copy empty
##  rename(LPFBS = PFBS)  %>% # PFBS is, and always was, linear...
##  mutate(TYPE = "hav") %>% 
  # Some YWD = 798 does not translate, changing to year
##  mutate(YWD = ifelse(ACCNR %in% c("C07/07689", "C07/07692", "C07/07693", "C07/07694", 
##                                   "C07/07695", "C07/07696", "C07/07697", "C07/07698", "C07/07699"), "07", YWD)) %>% 
  # Errors in weight/length and last minute fixes
##  mutate(TOTV = ifelse(ACCNR == "C07/06942-06953", 33.6, TOTV),
##         TOTL = ifelse(ACCNR == "C07/06942-06953", 17.6, TOTL),
##         TOTL = ifelse(ACCNR == "C17/07860-07871", 17.63, TOTL),
##         TOTL = ifelse(ACCNR == "C17/10560-10569", 21.8, TOTL),
##         TOTL = ifelse(ACCNR == "C09/01992-02003", 17.06, TOTL),
##         TOTL = ifelse(ACCNR == "C09/02004-02015", 17.46, TOTL),
##         DDTK = ifelse(ACCNR %in% paste0("P98/0", 2042:2053), NA, DDTK),
##         DDTK = ifelse(ACCNR %in% paste0("P02/0", 3939:3948), NA, DDTK),
##         DDTK = ifelse(ACCNR %in% paste0("P03/0", 3435:3444), NA, DDTK),
##         HCB = ifelse(ACCNR == "P01/04580", NA, HCB),
##         BHCH = ifelse(ACCNR %in% paste0("P90/0", 5446:5463), NA, BHCH))

# Then eggs ("uria")
#
##files_uria <- dir(path_hist) %>% subset(str_detect(., "uria") & !str_starts(., "~"))
##data_list_uria <- map(paste0(path_hist, files_uria), read_historical_file) # Read old data files to a list
# Join historical data files by column binding, note that this assumes that each 
# ACCNR occupies the same row in each of the Excel files (this is why "hav" and "uria" needs to be read separately)
##historical_uria <- reduce(data_list_uria, append_cols) %>%
##  rename_all(~ifelse(.x %in% c("FPRC*", "TPRC*", "MTPRC*", "LTPRC*"), .x, str_remove(.x, "\\*"))) %>%
##  rename(LPFDS = 'LPFDS...34') %>% # Variable LPFDS occurs twice, second copy empty
##  rename(LPFBS = PFBS)  %>% # PFBS is, and always was, linear...
##  mutate(TYPE = "uria", GRU = if_else(GENUS %in% c("STER", "HAEM"), 1, GRU))

# Fixing unit of PCBs: REMOVE WHEN FIX IS DONE IN SOURCE
##if(mean(filter(historical_uria, str_sub(ACCNR, 2, 2) == "9") %>% pull(`CB-77`), na.rm = TRUE) < 5){ # Check if fixed
##  historical_uria <- historical_uria %>% mutate(YEAR = as.numeric(str_sub(ACCNR, 2,3)),
##                                                YEAR = ifelse(YEAR < 50, YEAR + 2000, YEAR + 1900)) %>% 
##    mutate(`CB-77` = ifelse(YEAR < 2006, `CB-77` * 1000, `CB-77`),
##           `CB-81` = ifelse(YEAR < 2006, `CB-81` * 1000, `CB-81`),
##           `CB-126` = ifelse(YEAR < 2006, `CB-126` * 1000, `CB-126`),
##           `CB-169` = ifelse((YEAR < 2006) & (YEAR != 2001), `CB-169` * 1000, `CB-169`),
##           `CB-123` = ifelse(YEAR < 2006, `CB-123` * 1000, `CB-123`),
##           `CB-114` = ifelse(YEAR < 2006, `CB-114` * 1000, `CB-114`),
##           `CB-105` = ifelse(YEAR < 2006, `CB-105` * 1000, `CB-105`),
##           `CB-118U` = ifelse(YEAR < 2006, `CB-118U` * 1000 * 1000, `CB-118U`),
##           `CB-167` = ifelse(YEAR < 2006, `CB-167` * 1000, `CB-167`),
##           `CB-157` = ifelse(YEAR < 2006, `CB-157` * 1000, `CB-157`),
##           `CB-189` = ifelse(YEAR < 2006, `CB-189` * 1000, `CB-189`),
##           `CB-156` = ifelse((YEAR < 2006) & (YEAR != 2001), `CB-156` * 1000, `CB-156`),
##           CBEQV = ifelse(YEAR %in% c(2004, 2005), CBEQV * 1000, CBEQV)
##    ) %>% 
##    select(-YEAR)
##}



# Finally bind data together and save
#
full_data <- bind_rows(historical, new_data) %>% ### NB NB NB ALS
  mutate('FPRC*' = ifelse(is.na(.[['FPRC*']]), FPRC, .[['FPRC*']]),
         LTPRC = ifelse(is.na(LTPRC), .[['LTPRC*']], .[['LTPRC*']]))
write_csv(full_data, "clean_data/full_data_2019.csv")

