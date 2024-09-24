library(tidyverse)
library(readxl)
library(openxlsx)
library(readr)
path <- "E:/projects/Baby-Germs/KLEPP/KLEPP"
# path <- "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA"
# path <- "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU"
# path <- "E:/projects/Baby-Germs/ECOLI/downstream-analysis/ECOLI"
# path <- "E:/projects/Baby-Germs/FAECA/downstream-analysis/FAECA"
# path <- "E:/projects/Baby-Germs/FAECI/downstream-analysis/FAECI"
# path <- "E:/projects/Baby-Germs/KLEPP/KLEPP"
# path <- "E:/projects/Baby-Germs/KLEPP/KLEPP"

# path <- "D:/Terra-Informatix/baby-germs/KLEPP"

paths <- list.files(path = path, recursive = T, pattern = "metadata.csv",
           full.names = T)

# KLEPP
# paths <- paths[str_detect(paths,"cluster-analysisZ01")]
paths <- paths[str_detect(paths,"clusters-Z03")]
# ACIBA
# paths <- paths[str_detect(paths,"cluster-analysis9")]
# STAAU
# paths <- paths[str_detect(paths,"cluster-analysis6")]
# ECOLI
# paths <- paths[str_detect(paths,"cluster-analysis2")]
#FAECA
# paths <- paths[str_detect(paths,"cluster-analysis2")]
# FAECI
# paths <- paths[str_detect(paths,"cluster-analysis2")]


dat_list <- purrr::map(paths, read_csv) #%>%
  # mutate(across(.fns = as.character))

# surveys <- type_convert(dat_list)

cmd_df <-map_dfr(dat_list, ~.x %>% 
                   mutate(across(everything(), as.character)))

# # bind_rows(dat_list)
# cmd_df <- as.data.frame(do.call(rbind, dat_list)) %>%
# mutate(across(c(1,2,3,4,6,7,11,12,13,14,15),.fns = as.character))


write_csv(cmd_df, file.path(path,"klepp_metadata.csv") )
# write_csv(cmd_df, file.path(path,"aciba_metadata.csv") )
# write_csv(cmd_df, file.path(path,"staau_metadata.csv") )
# write_csv(cmd_df, file.path(path,"ecoli_metadata.csv") )
# write_csv(cmd_df, file.path(path,"faeca_metadata.csv") )
# write_csv(cmd_df, file.path(path,"faeci_metadata.csv") )

# get_scatter <- function(df,hospital){

names(cmd_df)[str_detect(names(cmd_df),"Var_00")] <- "Hospital"
hosp_vec <- cmd_df %>% pull(Hospital) %>% unique()
# hosp_vec <- cmd_df %>% pull(Var_00) %>% unique()

year_id <- lubridate::year(cmd_df$TakenDate)
week_id <- lubridate::week(cmd_df$TakenDate)

cmd_df$Epiweek <- paste(year_id,week_id,sep = ".")

# Ward	WardType	
epi_vars <- c("sampleid",	"Hospital",	"TakenDate",	"Epiweek",	"ST")
snp_epi <- names(cmd_df)[str_detect(names(cmd_df),"SNP.EPI.cluster")]

sel_vec <- c(epi_vars,snp_epi)

for(i in seq_along(hosp_vec)){
  hospital <- hosp_vec[i]
  print(hospital)
  
  plotDF <- cmd_df %>% 
    ungroup() %>%
    dplyr::select(any_of(sel_vec)) %>%
    # dplyr::select(c(1,2,3,4,5,6,7,11,12,13,14,15)) %>%
    group_by(Hospital,Epiweek) %>%
    distinct(sampleid, .keep_all = TRUE) %>%
    dplyr::mutate(n_cases=n()) %>%
    ungroup() %>%
    dplyr::filter(str_detect(Hospital,hospital)) %>%
    pivot_longer(cols = all_of(snp_epi),names_to = "Groups",values_to = "Clusters") %>%
    mutate(Clusters = na_if(Clusters, "bg")) %>%
    dplyr::arrange(TakenDate) %>% 
    mutate(denserank = data.table::rleid(TakenDate)) %>%
    dplyr::mutate(id=dplyr::row_number()) %>%
    dplyr::rename("Cases_count"=n_cases)
  
  
  plotDF$Epiweek = factor(plotDF$Epiweek, levels=unique(plotDF$Epiweek)[plotDF$id], ordered = T)
  plotDF$Cases_count = dplyr::coalesce(plotDF$Cases_count,1)
  plotDF$Cases_count = as.factor(plotDF$Cases_count)
  plotDF$Clusters = as.factor(plotDF$Clusters)
  plotDF$ST = as.factor(plotDF$ST)
  
  
  
  # Plot bar plot and scatterplot - epicurves -------------------------------
  px1 <- plotDF %>% group_by(WardType,Epiweek) %>% dplyr::mutate(count=n())
  p1 <- ggplot(px1, 
               aes(x = Epiweek,
                   y = ST, #WardType, #Clusters,
                   # y = count, #,
                   color=Clusters,
                   shape = WardType
               )) +
    geom_point(alpha = .9, size=3,position=position_jitter(h=0.07,w=0.15)) + 
    facet_grid(Groups ~ .) +
    # scale_fill_brewer(palette = "Dark2") +
    labs(title = hospital, #paste0(""),
         y = "Sequence Types",
         x = "Epiweek") +
    theme_bw() +
    # scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
    # theme(panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank()) +
    theme(axis.title = element_text(size = 10)) +
    theme(axis.text = element_text(size = 7.5)) +
    theme(axis.text.x = element_text(angle=90,size=7.5,vjust = 0.5, hjust = 1)) 
  # +
  #   guides(color = guide_colorbar(title = "Horsepower"),
  #          shape = guide_legend(title = "Weight"), size = guide_legend(title = "Gear")
  #   )
  
  hosp_id <- str_replace_all(hospital," ","-")
  ggsave(file.path(path,paste0(hosp_id,"Scatterplot",".png")),p1,width = 10, height = 8)
}  
  
  
  
  
  







  

  
#    %>% 
#   left_join(clusterSet3, by="sampleID") %>% 
#   dplyr::rename("Epiweek"=epiyearweek) %>% 
#   ungroup() %>%
#   group_by(Epiweek) %>%
#   distinct(sampleID, .keep_all = TRUE) %>%
#   dplyr::mutate(n_cases=n()) %>%
#   ungroup() %>%
#   dplyr::arrange(TakenDate) %>%
#   mutate(denserank = data.table::rleid(TakenDate)) %>%
#   dplyr::mutate(id=dplyr::row_number()) %>%
#   left_join(mlst,by=c("sampleID"="FILE")) %>%
#   dplyr::rename("Cases_count"=n_cases) #%>% print(n=40)
# 
# plotDF$Epiweek = factor(plotDF$Epiweek, levels=unique(plotDF$Epiweek)[plotDF$id], ordered = T)
# plotDF$Cases_count = dplyr::coalesce(plotDF$Cases_count,1)
# plotDF$Cases_count = as.factor(plotDF$Cases_count)
# plotDF$Clusters = as.factor(plotDF$Clusters)
# plotDF$ST = as.factor(plotDF$ST)



# Plot bar plot and scatterplot - epicurves -------------------------------

library(visdat)
# library(ggnewscale)
# vis_dat(airquality)
vis_dat(plotDF)

px1 <- plotDF %>% group_by(WardType,Epiweek) %>% dplyr::mutate(count=n())
p1 <- ggplot(px1, 
             aes(x = Epiweek,
                 y = ST, #WardType, #Clusters,
                 # y = count, #,
                 color=Clusters,
                 shape = WardType
             )) +
  geom_point(alpha = .9, size=4,position=position_jitter(h=0.07,w=0.15)) + 
  # scale_fill_brewer(palette = "Dark2") +
  labs(title = paste0("Clusters based on Epi definition: SNPs <= ",snpco," + ",daysco,"-day window"),
       y = "Sequence Types",
       x = "Epiweek") +
  theme_bw() +
  # scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5, hjust = 1)) 
# +
#   guides(color = guide_colorbar(title = "Horsepower"),
#          shape = guide_legend(title = "Weight"), size = guide_legend(title = "Gear")
#   )


ggsave(file.path(work_dir,paste0("Scatterplot.",date_var,".png")),p1,width = 10, height = 8)

