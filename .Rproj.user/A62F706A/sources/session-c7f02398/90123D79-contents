library(tidyverse)
library(readxl)
library(openxlsx)
library(ape)

# date when analysis source files were generated
date_var <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")
# date_var = "2023-10-26"
# date_var = "2023-10-28"
# date_var <- date_var -1
# date_var = "2024-07-26"
clsters=NULL
# main_path<-file.path("E:/projects/Baby-Germs/KLEPP/KLEPP")
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/DORA-NGINZA/transmission/DORA_ward.csv"
# wards<-list.files(path = main_path,pattern = "_ward.csv",full.names = T)
# 
# read_csv(wards, col_names = F) %>% 
#   left_join(amr_predicted,by=c("X1"="input_file_name")) %>%
#   replace(is.na(.), "no") %>%
#   # mutate_at(genes, ~replace_na("no"))
#   print(n=40)
# cluster_file = file.path(work_dir,paste("Clusters",date_var,"xlsx",sep = "."))



# klepp_clusters <- file.path("C:/Users/Stanfordk/Documents/Terra-Info/ACIBA") # ACIBA
# klepp_clusters <- file.path("C:/Users/Stanfordk/Documents/Terra-Info/KLEPP") # KLEPP
# ast_path <- file.path(klepp_clusters,"amr_aciba.xlsx") # ACIBA

# "D:/Terra-Informatix/baby-germs/KLEPP/TEMBISA/clusters3009"
# path to work directory
wkdir <- file.path("D:/Terra-Informatix/baby-germs/KLEPP")
# wkdir <- file.path("D:/Terra-Informatix/baby-germs/ACIBA")
# wkdir <- file.path("D:/Terra-Informatix/baby-germs/STAAU")
# wkdir <- file.path("D:/Terra-Informatix/baby-germs/ECOLI")
# wkdir <- file.path("D:/Terra-Informatix/baby-germs/FAECI")
# wkdir <- file.path("D:/Terra-Informatix/baby-germs/FAECA")

# Path to metadata
klepp_clusters <- file.path("E:/projects/Baby-Germs/KLEPP") # KLEPP
# klepp_clusters <- file.path("E:/projects/Baby-Germs/ACIBA") # ACIBA
# klepp_clusters <- file.path("E:/projects/Baby-Germs/STAAU") # STAAU
# klepp_clusters <- file.path("E:/projects/Baby-Germs/ECOLI") # ECOLI
# klepp_clusters <- file.path("E:/projects/Baby-Germs/FAECI") # FAECI
# klepp_clusters <- file.path("E:/projects/Baby-Germs/FAECA") # FAECA


# ASTs and predicted AMR files
# raw_ast_path <- file.path("E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/amr_aciba.xlsx")
ast_path <- file.path(klepp_clusters,"MDRO_profiles.xlsx") # KLEPP
# ast_path <- file.path(klepp_clusters,"amr_aciba.xlsx") # ACIBA
# ast_path <- file.path(klepp_clusters,"amr_profiles.xlsx") # STAAU
# ast_path <- list.files(klepp_clusters,pattern="amr_ecoli.xlsx", recursive = T, full.names = T) # ECOLI
# ast_path <- list.files(klepp_clusters,pattern="^amr_efaecium.xlsx", recursive = T, full.names = T) # FAECI
# ast_path <- list.files(klepp_clusters,pattern="^amr_efaecalis.xlsx", recursive = T, full.names = T) # FAECA

dirs_vec <- list.dirs(path = klepp_clusters)
amr_dir <- dirs_vec[str_detect(dirs_vec,"-hmz/hamronize$")]
# amr_path <- file.path(klepp_clusters,"amr_predicted_genes.xlsx")
amr_path <-list.files(path = amr_dir,pattern = "combined_report_hrz_sorted.tsv",full.names = T)





# # mlst profile
# list.files(klepp_clusters, pattern = "05.mlst", recursive = T)
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/TEMBISA/transmission/05.mlst.xlsx"
# Get profile column ------------------------------------------------------
# MDR profiles
# KLEPP
# mdro <- read_excel(path = ast_path) %>%
#   drop_na(contig_name) %>%
#   dplyr::select(1,12) %>%
#   mutate(ESBL=if_else(! str_detect(profile,"ESBL") | is.na(profile),"no","yes"),
#          MDR=if_else(str_detect(profile,"and MDR") | str_detect(profile,"^MDR"),"yes","no"),
#          Aminoglycosides=if_else(str_detect(profile,"Aminoglycosides"),"yes","no"),
#          CRE=if_else(str_detect(profile,"CR"),"yes","no")) %>%
#   replace(is.na(.), "no") %>%
#   dplyr::select(-2) %>%
#   # dplyr::rename("ESBL_autocolor"=ESBL,
#   #               "MDR_autocolor"=MDR,
#   #               "Aminoglycosides_autocolor"=Aminoglycosides,
#   #               "CR_autocolor"=CR) %>%
#   print(n=100)

# NOW USING THIS OPTION
if (str_detect(ast_path,"STAAU")){
  incl_cols <- c("newid","hospitalname","wardtype","contig","contig_name","profile","Methicillin_resistance","Vancomycin_resistance")
  mdro <- read_excel(ast_path) %>%
    dplyr::select(any_of(incl_cols)) %>%
    dplyr::filter(contig == "yes") %>%
    dplyr::mutate(contig_name =str_remove_all(contig_name,"\"")) %>%
    dplyr::mutate(contig_name = str_remove(contig_name,"_.*|.fasta")) %>%
    dplyr::select(-c(contig,newid)) %>%
    drop_na(contig_name) %>%
    # mutate(ESBL=if_else(! str_detect(profile,"ESBL") | is.na(profile),"no","yes"),
    #        MDR=if_else(str_detect(profile,"and MDR") | str_detect(profile,"^MDR"),"yes","no"),
    #        Aminoglycosides=if_else(str_detect(profile,"Aminoglycosides"),"yes","no"),
    #        CRE=if_else(str_detect(profile,"CR"),"yes","no")) %>%
    replace(is.na(.), "no") %>%
    dplyr::mutate(across(.cols=-c(1,2,3),
                         .fns = ~if_else(. == "no","no","yes"))) #%>% print(n=100)
}else if(str_detect(amr_path,"FAECA|FAECI")){
  incl_cols <- c("newid","hospitalname","wardtype","contig","contig_name","profile")
  mdro <- read_excel(ast_path) %>%
    dplyr::select(any_of(incl_cols)) %>%
    # dplyr::filter(contig == "yes") %>%
    dplyr::mutate(contig_name =str_remove_all(contig_name,"\"")) %>%
    dplyr::mutate(contig_name = str_remove(contig_name,"_.*")) %>%
    dplyr::mutate(contig_name = str_remove(contig_name,"_.*|.fasta")) %>%
    # dplyr::select(-c(contig,newid)) %>%
    drop_na(contig_name) %>%
    mutate(Ampicillin=if_else(str_detect(profile,"ampicillin"),"no","yes"),
           Teicplanin=if_else(str_detect(profile,"teicplanin"),"yes","no"),
           Linezolid=if_else(str_detect(profile,"linezolid"),"yes","no"),
           Vancomycin=if_else(str_detect(profile,"vancomycin"),"yes","no")) %>%
    replace(is.na(.), "no") %>%
    dplyr::select(-c(profile)) #%>% print(n=100)
  
  
  
}else{
  incl_cols <- c("newid","hospitalname","wardtype","contig","contig_name","profile")
  mdro <- read_excel(ast_path) %>%
    dplyr::select(any_of(incl_cols)) %>%
    # dplyr::filter(contig == "yes") %>%
    dplyr::mutate(contig_name =str_remove_all(contig_name,"\"")) %>%
    dplyr::mutate(contig_name = str_remove(contig_name,"_.*")) %>%
    dplyr::mutate(contig_name = str_remove(contig_name,"_.*|.fasta")) %>%
    # dplyr::select(-c(contig,newid)) %>%
    drop_na(contig_name) %>%
    mutate(ESBL=if_else(! str_detect(profile,"ESBL") | is.na(profile),"no","yes"),
           MDR=if_else(str_detect(profile,"MDR") | str_detect(profile,"^MDR"),"yes","no"),
           Aminoglycosides=if_else(str_detect(profile,"Aminoglycosides"),"yes","no"),
           CRE=if_else(str_detect(profile,"CR"),"yes","no")) %>%
    replace(is.na(.), "no") %>%
    dplyr::select(-c(profile)) #%>% print(n=100)
}

# Begin analysis ----------------------------------------------------------

# E:\projects\Baby-Germs\KLEPP\KLEPP\TEMBISA\cluster-analysisZ01
# Get paths for cluster directories
dirs_vec = vector()
dirs_vec <- list.files(path = wkdir, #klepp_clusters, 
           include.dirs = T, 
           full.names = T, 
           pattern = "clusters3009", #"clusters-Z02", #"cluster-analysisZ01", 
           recursive = T)



# tvect <- amr_predicted %>% pull(input_file_name)
# filesList[[1]] %>% filter(label %in% tvect)
# input_dir <- dirs_vec[[6]]

for (i in seq_along(dirs_vec)){
  print(dirs_vec[i])
  
  # get paths for cluster files
  input_dir <- dirs_vec[[i]]
  
  
  # get hospital name
  hosp <- basename(dirname(dirs_vec[[i]]))
  hosp_prefix <- str_split(hosp, pattern = " |-|_")[[1]][1]
  
  # output directory
  output_dir <- file.path(input_dir,"phylogenetic-tree-files")
  if(! dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  search_dir <- dirname(input_dir)
  
  # cluster_file
  cluster_paths <- list.files(path = input_dir,pattern = paste0("Core-SNP-Clusters.",date_var),
                          recursive = T,
                          full.names = T)

  # cluster_files <- map(cluster_paths, readxl::read_excel)
  
  fn_getname <- function(x) { str_remove(basename(dirname(x)),"clusters_") }
  
  
  cluster_dfs <- cluster_paths |>
    purrr::set_names(fn_getname) |>
    map(readxl::read_excel) #|>
    #list_rbind(names_to = "cutoff")
                     
  # str_remove(basename(dirname(cluster_paths[1])),"clusters_")
  
  # View(cluster_files)
  
  # SNP-EPI cluster files
  epi_clusters_paths=list.files(path = input_dir,
                                pattern = paste("Clusters-SNPs-Epi-Definition",date_var,"xlsx",sep = "."),
                                recursive = T,full.names = T)
  
  epi_clusters_dfs <- epi_clusters_paths |>
    purrr::set_names(fn_getname) |>
    map(readxl::read_excel)

  # MLST files
  mlst1=NULL
  mlst_path=list.files(path = input_dir,
                                pattern = paste("MLST_profile",date_var,"xlsx",sep = "."),
                                recursive = T,full.names = T)[1]
  
  
  
  if(! is.na(mlst_path)){
    mlst1 <- read_excel(mlst_path) |>
      dplyr::select(c(1,5,3,2,ST,any_of(Var_02))) 
  }

  
  # mlst1 <- read_excel(mlst_path) |>
  #   dplyr::select(-Clusters)
  
  # mlst_df <- mlst_path |>
  #   purrr::set_names(fn_getname) |>
  #   map(readxl::read_excel) |>
  #   dplyr::select(-Clusters)
  
  # # mlst file
  # mlst1_path<-dirname(dirs_vec[i])
  # mlst1<-list.files(file.path(mlst1_path), pattern = "05.mlst", recursive = T, full.names = T)[1]
  
  # AMR predicted genes
  # ESBL: blaCTX-M-15
  # Carbapenemase: blaOXA-181 and blaOXA-48
  # regex("no|something", ignore_case = TRUE)
  if (str_detect(amr_path,"STAAU")){
    
    amr_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,regex("mecA|mecC", ignore_case=T)) | 
                      str_detect(X2,regex("mec-A|mec-C",ignore_case=T)) | 
                      str_detect(X3,regex("methicillin|Linezolid|Vancomycin|SCCmec|MRSA",ignore_case=T))) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      # dplyr::select(-c(`blaOXA-48_1`,`blaOXA-181_1`)) %>% 
      dplyr::select(!matches("_.*")) %>%
      dplyr::rename("input_file_name" = X1) %>% 
      replace(is.na(.), "no") %>%
      dplyr::mutate(across(.cols=-1,
                           .fns = ~if_else(. == "no","no","yes"))) #%>%
      # print(n=40)
    
  }else if(str_detect(amr_path,"FAECA|FAECI")){
    
    
    # Resistance genes for Enterococci
    # Vancomycin - VanA, VanB, VanC-1, Van-C2/3
    # Aminoglycosides genes - aac(6')-aph(2'')
    
    
    vancomycin_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"VanA") |  str_detect(X3,"VanB") | str_detect(X3,"vancomycin") | str_detect(X3,"VanC")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    aminoglyco_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"aac") |  str_detect(X3,"aminoglycoside") | str_detect(X3,"aad")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    
    A1 <- vancomycin_predicted %>%
      dplyr::select(1)
    
    
    A2 <- aminoglyco_predicted %>%
      dplyr::select(1)
    
    
    A3 <- dplyr::union(A1,A2)
    
    
    amr_predicted <- A3 %>%
      left_join(vancomycin_predicted, by=c("input_file_name")) %>%
      left_join(aminoglyco_predicted, by=c("input_file_name")) %>%
      replace(is.na(.), "no") %>%
      # dplyr::select(!matches("_1")) %>%
      # dplyr::select(!matches("_2")) %>%
      # dplyr::select(!matches("_5")) %>%
      dplyr::mutate(across(.cols=-1,
                           .fns = ~if_else(. == "no","no","yes"))) #%>% print(n=40)

    
  }else if(str_detect(amr_path,"ECOLI")){
  
    # Resistance genes for E. coli
    # ESBL genes – CTX-M-14, CTX-M-15, TEM-1B, SHV-2, ampC
    # Carbapenemases genes – OXA genes
    # Aminoglycosides genes – aada1, aada5
    # Plasmid-mediated quinolone resistance (PMQR) genes – qnrS,
    # Sulphonamides – sul1, sul2
    # Colistin - mcr genes
  
    cr_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"OXA-") | str_detect(X2,"blaOXA-48") | str_detect(X3,"carbapenem")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      # dplyr::select(-c(`blaOXA-48_1`,`blaOXA-181_1`)) %>% 
      dplyr::rename("input_file_name" = X1) #%>% 
      # print(n=40)
    # mutate(`blaOXA-48` = if_else(`blaOXA-48`!="blaOXA-48" | is.na(`blaOXA-48`),"no","yes"),
    #        `blaOXA-181` = if_else(`blaOXA-181` !="blaOXA-181" | is.na(`blaOXA-181`),"no","yes")) %>%
    # print(n=40)
    
    
    
    esbl_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"blaCTX") |  str_detect(X3,"extended") | str_detect(X3,"-1B")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    aminoglyco_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"aac") |  str_detect(X3,"aminoglycoside") | str_detect(X3,"aad")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    quinolone_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"qnrS") |  str_detect(X3,"qnrs") | str_detect(X3,"quinolone")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    sulphonamides_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"sul") |  str_detect(X3,"sulfonamide")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    colistin_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"mcr") |  str_detect(X3,"colistin")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    # amr_predicted <- esbl_predicted %>%
    #   left_join(cr_predicted, by=c("input_file_name")) %>%
    #   mutate(`blaOXA-48` = if_else(`blaOXA-48`!="blaOXA-48" | is.na(`blaOXA-48`),"no","yes"),
    #          `blaOXA-181` = if_else(`blaOXA-181` !="blaOXA-181" | is.na(`blaOXA-181`),"no","yes")) %>%
    #   mutate(`blaCTX-M-15` = if_else(`blaCTX-M-15`!="blaCTX-M-15" | is.na(`blaCTX-M-15`),"no","yes"),
    #          `blaSHV-187` = if_else(`blaSHV-187` !="blaSHV-187" | is.na(`blaSHV-187`),"no","yes"),
    #          `blaSHV-106` = if_else(`blaSHV-106` !="blaSHV-106" | is.na(`blaSHV-106`),"no","yes"),
    #          `blaSHV-38` = if_else(`blaSHV-38` !="blaSHV-38" | is.na(`blaSHV-38`),"no","yes")) %>%
    #   print(n=40)
    
    
    
    
    A1 <- cr_predicted %>% 
      dplyr::select(1)
    
    A2 <- esbl_predicted %>% 
      dplyr::select(1)
    
    A3 <- dplyr::union(A1,A2)
    
    A4 <- aminoglyco_predicted %>% 
      dplyr::select(1) %>%
      dplyr::union(A3)
    
    A5 <- quinolone_predicted %>% 
      dplyr::select(1) %>%
      dplyr::union(A4)
    
    A6 <- sulphonamides_predicted %>% 
      dplyr::select(1) %>%
      dplyr::union(A5)
    
    A7 <- colistin_predicted %>% 
      dplyr::select(1) %>%
      dplyr::union(A6)
    
    
    
    amr_predicted <- A7 %>%
      left_join(esbl_predicted, by=c("input_file_name")) %>%
      left_join(cr_predicted, by=c("input_file_name")) %>%
      left_join(aminoglyco_predicted, by=c("input_file_name")) %>%
      left_join(quinolone_predicted, by=c("input_file_name")) %>%
      left_join(sulphonamides_predicted, by=c("input_file_name")) %>%
      left_join(colistin_predicted, by=c("input_file_name")) %>%
      replace(is.na(.), "no") %>%
      dplyr::select(!matches("_1")) %>%
      dplyr::select(!matches("_2")) %>%
      dplyr::select(!matches("_5")) %>%
      dplyr::mutate(across(.cols=-1,
                           .fns = ~if_else(. == "no","no","yes"))) #%>% print(n=40)
    
  }else {
    cr_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"OXA-") | str_detect(X2,"blaOXA-48") | str_detect(X3,"carbapenem")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) 
    
    
    esbl_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"blaCTX") |  str_detect(X3,"extended") | str_detect(X3,"-1B")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    aminoglyco_predicted <- read_tsv(amr_path, col_names = F) %>% #read_excel(amr_path) %>%
      dplyr::filter(str_detect(X2,"aac") |  str_detect(X3,"aminoglycoside") | str_detect(X3,"aad")) %>%
      dplyr::select(1,2) %>% distinct() %>% 
      pivot_wider(names_from = "X2", values_from = "X2") %>%
      dplyr::rename("input_file_name" = X1) #%>% print(n=40)
    
    
    
    A1 <- cr_predicted %>% 
      dplyr::select(1)
    
    A2 <- esbl_predicted %>% 
      dplyr::select(1)
    
    A3 <- dplyr::union(A1,A2)
    
    A4 <- aminoglyco_predicted %>% 
      dplyr::select(1) %>%
      dplyr::union(A3) 
    
    amr_predicted <- A4 %>%
      left_join(esbl_predicted, by=c("input_file_name")) %>%
      left_join(cr_predicted, by=c("input_file_name")) %>%
      left_join(aminoglyco_predicted, by=c("input_file_name")) %>%
      replace(is.na(.), "no") %>%
      dplyr::select(!matches("_1")) %>%
      dplyr::select(!matches("_2")) %>%
      dplyr::select(!matches("_5")) %>%
      dplyr::mutate(across(.cols=-1,
                           .fns = ~if_else(. == "no","no","yes"))) #%>% print(n=40)
    
  }
    
  
    
  # add clusters file
  # if(is.na(cluster_file)){
  if(all(is.na(cluster_paths))){
    
    # copy and rename consensus tree files
    tree_path <- list.files(dirname(input_dir), pattern = "contree", recursive = T, full.names = T)
    if(length(tree_path) != 0){
      tree_name<- str_split(basename(tree_path),"\\.")[[1]][1]
      file.copy(tree_path,file.path(output_dir,paste0(tree_name,".nwk")))
      nwk_labels <- ape::read.tree(tree_path)$tip.label
      
      # mlst
      # mlst_path1 <- list.files(path = klepp_clusters,pattern = "05.mlst.xlsx", recursive = T, full.names = T)[1]
      # mlst_path1 <- list.files(path = search_dir,pattern = "05.mlst.xlsx", recursive = T, full.names = T)[1]
      
      # mlst_path1 <- list.files(path = search_dir,pattern = "MLST_profile", recursive = T, full.names = T)[1]
      
      mlst_path_vec <- list.files(path = klepp_clusters,pattern = "05.mlst.xlsx", recursive = T, full.names = T)
      mlst_path1 <- mlst_path_vec[str_detect(mlst_path_vec,hosp_prefix)]
      
      mlst1 <- read_excel(mlst_path1) |>
        dplyr::select(c(1,5,3,2)) |>
        dplyr::filter(sampleID %in% nwk_labels) #|>
        # dplyr::rename("sampleID" = FILE)
      write.xlsx(mlst1,file.path(output_dir,"MLST.xlsx"),overwrite = T)
      
      
      if(nrow(amr_predicted)!=0){
        tree_annot <- amr_predicted %>% 
          dplyr::filter(input_file_name %in% nwk_labels) 
        write.xlsx(tree_annot,file.path(output_dir,"AMR-predicted.xlsx"),overwrite = T)
      }
      if(nrow(mdro)!=0){
        tree_ast <- mdro %>% 
          dplyr::filter(contig_name %in% nwk_labels) #%>%
          #left_join(mlst1,by=c("contig_name"="sampleID"))
        write.xlsx(tree_ast,file.path(output_dir,"ASTs.xlsx"),overwrite = T)
      }
    }else{
      # mlst
      # mlst_path1 <- list.files(path = klepp_clusters,pattern = "05.mlst.xlsx", recursive = T, full.names = T)[1]
      # mlst_path1 <- list.files(path = search_dir,pattern = "MLST_profile", recursive = T, full.names = T)[1]
      
      mlst_path_vec <- list.files(path = klepp_clusters,pattern = "05.mlst.xlsx", recursive = T, full.names = T)
      mlst_path1 <- mlst_path_vec[str_detect(mlst_path_vec,hosp_prefix)]
      
      mlst1 <- read_excel(mlst_path1) |>
        dplyr::select(c(1,5,3,2)) 
      write.xlsx(mlst1,file.path(output_dir,"MLST.xlsx"),overwrite = T)
      
      if(nrow(amr_predicted)!=0){
        tree_annot <- amr_predicted 
        write.xlsx(tree_annot,file.path(output_dir,"AMR-predicted.xlsx"),overwrite = T)
      }
      if(nrow(mdro)!=0){
        tree_ast <- mdro 
        write.xlsx(tree_ast,file.path(output_dir,"ASTs.xlsx"),overwrite = T)
      }
    }
    
    # next
  }#else{
   #clsters <- read_excel(cluster_file) %>% 
   #   dplyr::filter(sampleID != "reference") %>%
   #  dplyr::select(1,2) %>%
   #   # dplyr::rename("cluster"=coreSNPcluster)
   #    print(n=40)
   #  }
  if(all(is.na(cluster_paths))){
    next
  }else{
    finalList <- list()
    for(i in seq_along(cluster_dfs)) { 
      var1 <- str_replace(names(cluster_dfs[i]),"cluster-analysis_","") %>%
        str_replace(.,"SNPcutoff","SNP.cluster")
      
      if(i==1){
        y=1
        finalList[[y]] <- mlst1
      }
      
      # add SNPclusters
      v=max(seq_along(finalList))
      if(any(colnames(finalList[[v]]) == var1)){
        y=i-1
        next
      }
      
      if(length(cluster_dfs[[i]]) == 0){
        next
      }
      finalList[[i]] <- finalList[[v]] %>%
        left_join(cluster_dfs[[i]], by=c("sampleID")) %>%
        dplyr::rename("{var1}" := coreSNPcluster)
      
      if(i==max(seq_along(cluster_dfs))){
        break
      }
    }
    
    w=max(seq_along(finalList))
    main_df1 <- finalList[[w]]
  }
 
  
  

  
  if(all(is.na(epi_clusters_paths))){
    next
  }else{
    finalList <- list()
    for(i in seq_along(epi_clusters_dfs)) { 
      var1 <- names(epi_clusters_dfs[i]) %>%
        str_replace(.,".*SNPcutoff","SNP.EPI.cluster") %>%
        str_replace(.,"_Days",".")
      
      if(i==1){
        y=1
        finalList[[y]] <- main_df1
      }
      
      # add SNPclusters
      v=max(seq_along(finalList))
      if(any(colnames(finalList[[v]]) == var1)){
        y=i-1
        next
      }
      
      if(length(epi_clusters_dfs[[i]]) == 0){
        next
      }
      finalList[[i]] <- finalList[[v]] %>%
        left_join(epi_clusters_dfs[[i]], by="sampleID") %>%
        dplyr::rename("{var1}" := Clusters) %>%
        dplyr::select(-c(SNPs,Days))
      
      if(i==max(seq_along(epi_clusters_dfs))){
        break
      }
    }
    
    
    w1=max(seq_along(finalList))
    main_df2 <- finalList[[w1]]
  }
  

  if(exists('main_df2') && is.data.frame(get('main_df2'))){
    df_final <- main_df2 %>% 
      as.data.frame() %>%
      replace(is.na(.), "bg") %>%
      # dplyr::arrange(TakenDate) %>%
      left_join(mdro,by=c("sampleID"="contig_name")) %>%
      left_join(amr_predicted,by=c("sampleID"="input_file_name")) %>%
      dplyr::rename("sampleid"=sampleID) %>%
      replace(is.na(.), "no")
    
    
    write_csv(df_final, paste0(output_dir,"/","metadata", ".csv"))
    
  }
  
  # copy and rename consensus tree files
  tree_path <- list.files(dirname(input_dir), pattern = "contree", recursive = T, full.names = T)
  if(length(tree_path) != 0){
    tree_name<- str_split(basename(tree_path),"\\.")[[1]][1]
    file.copy(tree_path,file.path(output_dir,paste0(tree_name,".nwk")))
  }
 
  main_df2 <- NULL
}
  
  
  # files <- list.files(path = input_dir, pattern = '_SNPcutoff',full.names = T, recursive = T)
  # files <- files[str_detect(files,"xlsx")]
  # # str_split(basename(files),"\\.")[[1]][1]
  # filenames <- str_remove(str_remove(basename(files),"transmission_ntwk_nodes_"),"\\.xlsx")
  # names(files) <- filenames
  # filesList <- purrr::map(files, read_excel)
  # 
  # 
  # genes <- names(amr_predicted)[-1]
  # 
  # # combine different metdata files
  # finalList <- list()
  # for(i in seq_along(filesList)){
  #   clst_name <- names(filesList[i])
  #   df_in <- filesList[[i]] %>% 
  #     as.data.frame() %>%
  #     dplyr::mutate(X1 = dplyr::coalesce(X1,"UNKNOWN"),
  #                   X2 = dplyr::coalesce(X2,"UNKNOWN"),
  #                   X3 = dplyr::coalesce(X3,"UNKNOWN"),
  #                   group=na_if(group,"bg")) %>%
  #     dplyr::mutate(Collectiondate = as.Date(Collectiondate),
  #                   Epiyear=lubridate::epiyear(Collectiondate),
  #                   Epiweek=lubridate::epiweek(Collectiondate),
  #                   Epiweek=paste(Epiyear,Epiweek,sep = ".")) %>%
  #     dplyr::select(-c(Epiyear,shape)) %>%
  #     # dplyr::group_by(Collectiondate) %>%
  #     dplyr::arrange(Collectiondate) %>%
  #     # mutate(denserank = data.table::rleid(Collectiondate)) %>%
  #     # dplyr::mutate(id=dplyr::row_number()) %>%
  #     mutate(Epiweek = factor(Epiweek, levels=unique(Epiweek)[id], ordered = T)) %>%
  #     left_join(mdro,by=c("label"="contig_name")) %>%
  #     left_join(clsters, by=c("label"="sampleID")) %>%
  #     left_join(amr_predicted,by=c("label"="input_file_name")) %>%
  #     dplyr::rename("SNP.cluster"=coreSNPcluster,
  #                   "SNP.EPI.cluster"=group, #
  #                   "sampleid"=label) %>%
  #     dplyr::select(-id) %>%
  #     dplyr::select(sampleid,X1,X2,X3,Collectiondate,Epiweek,SNP.EPI.cluster,SNP.cluster,everything()) %>%
  #     dplyr::rename("Hospital"=X1,
  #                   "Ward-Name"=X2, #
  #                   "Ward-Type"=X3) %>%
  #     replace(is.na(.), "no")
  #     #mutate_at(genes, ~replace_na("no"))
  #   
  #   clst_name <- str_remove_all(clst_name,"transmission_ntwk_snpEpiclusters_")
  #   names(df_in)[names(df_in)=="SNP.EPI.cluster"] <- clst_name
  #   
  #   finalList[[i]] <- df_in
  #   names(finalList)[i] <- names(filesList[i])
  #   
  #   
  # }
  
#   finalList %>%
#     names(.) %>% 
#     purrr::walk(~ write_csv(finalList[[.]], paste0(output_dir,"/", ., ".csv")))
#   
#   # copy and rename consensus tree files
#   tree_path <- list.files(dirname(input_dir), pattern = "contree", recursive = T, full.names = T)
#   tree_name<- str_split(basename(tree_path),"\\.")[[1]][1]
#   file.copy(tree_path,file.path(output_dir,paste0(tree_name,".nwk")))
# 
# }

# End ---------------------------------------------------------------------

# write.xlsx(amr_predicted,"KLEPP_predicted_AMR.xlsx")
