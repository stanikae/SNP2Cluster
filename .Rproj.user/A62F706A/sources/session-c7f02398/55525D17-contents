library(tidyverse)
# library(viridis)

dir <- file.path("D:/Terra-Informatix/baby-germs/KLEPP")



files_vec <- list.files(path = dir,
           pattern = "metadata.csv",
           recursive = T,
           full.names = T,
           include.dirs = F)


names_vec <- str_replace(str_replace(files_vec, paste0(dir,"/"),""),"/.*","")

fixd_vars_vec <- c("sampleid","ST","TakenDate","Var_00")
core_vec <- c("SNP.cluster11_Days14","SNP.cluster11_Days60","SNP.cluster20_Days14","SNP.cluster20_Days60","SNP.cluster25_Days45")
snp_epi_vec <- c("SNP.EPI.cluster11.14","SNP.EPI.cluster11.60","SNP.EPI.cluster20.14","SNP.EPI.cluster20.60","SNP.EPI.cluster25.45")
 


meta_df <- purrr::map(files_vec,read_csv,col_types = cols(.default = "c"),show_col_types = FALSE) %>%
  dplyr::bind_rows() %>% #names()
  dplyr::select(any_of(fixd_vars_vec),
                # any_of(core_vec),
                any_of(snp_epi_vec)) %>%
  pivot_longer(!all_of(fixd_vars_vec), names_to="definition", values_to="clusters")
# names(metad_list) <- names_vec



meta_df$Epiweek = factor(plotDF$Epiweek, levels=unique(plotDF$Epiweek)[plotDF$id], ordered = T)
meta_df$Cases_count = dplyr::coalesce(plotDF$Cases_count,1)
meta_df$Cases_count = as.factor(plotDF$Cases_count)
meta_df$clusters = as.factor(plotDF$Clusters)
meta_df$ST = as.factor(meta_df$ST)



px1 <- meta_df %>% 
  group_by(Var_00,clusters,definition) %>% 
  dplyr::filter(!is.na(clusters)) %>%
  dplyr::mutate(count=n())

p1 <- ggplot(px1, 
             aes(x = definition,
                 # y = clusters, #WardType, #Clusters,
                 y = count, #,
                 color=clusters,
                 shape = Var_00
             )) +
  geom_point(alpha = .9, size=2,position=position_jitter(h=0.07,w=0.15)) #+ 
  # # scale_fill_brewer(palette = "Dark2") +
  # labs(title = paste0("Clusters based on Epi definition: SNPs <= ",snpco," + ",daysco,"-day window"),
  #      y = "Sequence Types",
  #      x = "Epiweek") +
  # theme_bw() +
  # # scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()) +
  # theme(axis.title = element_text(size = 15)) +
  # theme(axis.text = element_text(size = 12)) +
  # theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5, hjust = 1))


# p1 + stat_ellipse()
p1 + facet_grid(rows = vars(Var_00))


ggplot(meta_df, aes(x=definition,y=clusters, colour = "Var_00", shape ="ST")) +
  geom_point()





