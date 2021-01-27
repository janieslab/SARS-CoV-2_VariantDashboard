###################################
## Script to prepare SARS-CoV-2  ##
##    Variant FASTA files        ##
##        for mapping            ##
###################################
##    By: Colby T. Ford, Ph.D.   ##
###################################


library(seqinr)
library(dplyr)
library(readr)
library(stringr)


### Parse Variant FASTA File

## 20G 677H Variants
variant_files_677H <- paste0("../table1Anew/",
                             list.files(path = "../table1Anew/",
                                        pattern = ".txt"))
## 20G/501Y Variants
variant_files_501Y <- paste0("../table1Bnew/",
                             list.files(path = "../table1Bnew/",
                                        pattern = ".txt"))


variant_files <- c(variant_files_677H,
                   variant_files_501Y)


variants_df <- data.frame(file = character(),
                          accession = character(),
                          position = integer(),
                          reference = character(),
                          variant = character())

for (file_iter in variant_files){
  
  seqs_iter <- seqinr::read.fasta(file = file_iter)
  
  variants_iter <- data.frame(accession = names(seqs_iter),
                              variant = unlist(getSequence(seqs_iter,
                                                           as.string=T))) %>% 
    mutate(file = file_iter,
           position = str_remove(file, "../table[A-Za-z0-9]+/") %>% 
             str_remove(".txt") %>% 
             str_match("[0-9]+") %>%
             as.integer(),
           # reference = str_remove(file, "../table[A-Za-z0-9]+/") %>% 
           #   str_remove("to[A-Z].txt") %>% 
           #   str_remove("[0-9]+"),
           variant = str_to_upper(variant)) %>% 
    relocate(file,
             accession,
             position,
             # reference,
             variant)
  
  variants_df <- rbind(variants_df, variants_iter)
}

## Add in Reference Nucleotides
refs <- read_csv("explanation_tables.csv")

variants_df <- variants_df %>% 
  left_join(refs %>%
              select(position, reference, mutation),
            by = "position")


## Parse Output Location Files

abbreviations <- read_csv("location_abbreviations.csv")

variants_locations_df <- variants_df %>% 
  mutate(
    abbreviation = case_when(
      str_detect(accession, "_AL_") ~ "AL",
      str_detect(accession, "_AK_") ~ "AK",
      str_detect(accession, "_AZ_") ~ "AZ",
      str_detect(accession, "_AR_") ~ "AR",
      str_detect(accession, "_CA_") ~ "CA",
      str_detect(accession, "_CO_") ~ "CO",
      str_detect(accession, "_CT_") ~ "CT",
      str_detect(accession, "_DE_") ~ "DE",
      str_detect(accession, "_FL_") ~ "FL",
      str_detect(accession, "_GA_") ~ "GA",
      str_detect(accession, "_HI_") ~ "HI",
      str_detect(accession, "_ID_") ~ "ID",
      str_detect(accession, "_IL_") ~ "IL",
      str_detect(accession, "_IN_") ~ "IN",
      str_detect(accession, "_IA_") ~ "IA",
      str_detect(accession, "_KS_") ~ "KS",
      str_detect(accession, "_KY_") ~ "KY",
      str_detect(accession, "_LA_") ~ "LA",
      str_detect(accession, "_ME_") ~ "ME",
      str_detect(accession, "_MD_") ~ "MD",
      str_detect(accession, "_MA_") ~ "MA",
      str_detect(accession, "_MI_") ~ "MI",
      str_detect(accession, "_MN_") ~ "MN",
      str_detect(accession, "_MS_") ~ "MS",
      str_detect(accession, "_MO_") ~ "MO",
      str_detect(accession, "_MT_") ~ "MT",
      str_detect(accession, "_NE_") ~ "NE",
      str_detect(accession, "_NV_") ~ "NV",
      str_detect(accession, "_NH_") ~ "NH",
      str_detect(accession, "_NJ_") ~ "NJ",
      str_detect(accession, "_NM_") ~ "NM",
      str_detect(accession, "_NY_") ~ "NY",
      str_detect(accession, "_NC_") ~ "NC",
      str_detect(accession, "_ND_") ~ "ND",
      str_detect(accession, "_OH_") ~ "OH",
      str_detect(accession, "_OK_") ~ "OK",
      str_detect(accession, "_OR_") ~ "OR",
      str_detect(accession, "_PA_") ~ "PA",
      str_detect(accession, "_RI_") ~ "RI",
      str_detect(accession, "_SC_") ~ "SC",
      str_detect(accession, "_SD_") ~ "SD",
      str_detect(accession, "_TN_") ~ "TN",
      str_detect(accession, "_TX_") ~ "TX",
      str_detect(accession, "_UT_") ~ "UT",
      str_detect(accession, "_VT_") ~ "VT",
      str_detect(accession, "_VA_") ~ "VA",
      str_detect(accession, "_WA_") ~ "WA",
      str_detect(accession, "_WV_") ~ "WV",
      str_detect(accession, "_WI_") ~ "WI",
      str_detect(accession, "_WY_") ~ "WY",
      str_detect(accession, "_Canada_QC_") ~ "Canada_QC",
      str_detect(accession, "_Canada_ON_") ~ "Canada_ON",
      str_detect(accession, "_Canada_BC_") ~ "Canada_BC",
      str_detect(accession, "_Jamaica_") ~ "Jamaica",
      str_detect(accession, "_Mexico_TAM_") ~ "Mexico_TAM",
      str_detect(accession, "_Wuhan_Hu_") ~ "Wuhan_Hu"
    )
  ) %>% left_join(abbreviations, by = "abbreviation")



# variants_locations_df <- variants_df %>%
#   left_join(locations_df, by = "accession")

write_csv(variants_locations_df,
          paste0("variants_locations_df_",Sys.Date(),".csv"))

