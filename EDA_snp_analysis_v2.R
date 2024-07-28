
# Set-up the environment --------------------------------------------------
src_path <- file.path(getwd())
# src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
source(file.path(src_path,"scripts","pacman_install_packages.R"))


date_var <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")
# date_var <- str_remove_all(date_var,"-")


# Functions ---------------------------------------------------------------

source(file.path(src_path,"functions","load_custom_functions.R"))

# Data load ---------------------------------------------------------------
data_paths <- file.path(src_path,"conf")
# source(file.path(data_paths,"Richael_s_pneumo.R"))
source(file.path(data_paths,"BabyGERMS_kpn_temb.R"))

# START ANALYSIS ----------------------------------------------------------

# Set SNP cut-off and EPI-days
# Run all steps using loop
# for(i in 1:nrow(comparisons)){
  # snpco=comparisons[i,1]
  # daysco=comparisons[i,2]

  # }

  # snpco=20
  # daysco=21
  
  
  # get metadata first -- 20240713 
  # The first 3 columns should have the following - in the order specified below:
  #   1. sample_id
  #   2. collection_date
  #   3. facility_name/hospital_name/community codes/regions/metros (Var_00)
  #   4. Ward_name (or ward type) if 3 is hospital (Var_01)
  
  # get BabyGERMS metadata
  # metdata now included in the dates file: 2023-05-01
  mx <- read_csv(dates_path, col_names = F) %>% ncol()
  datesDF <- read_csv(dates_path, col_names = F) %>% dplyr::select(all_of(mx),4,1,3,2) #,1,2,3,4)
  names(datesDF) <- c("sampleID","TakenDate","Hospital","WardType","Ward")
  datesDF$Ward <- dplyr::coalesce(datesDF$Ward,"UNKNOWN")
  datesDF$WardType <- dplyr::coalesce(datesDF$WardType,"UNKNOWN")
  
  
  # # Richael's data
  # vec_incl_filt <- c("FILE","Sampling_date","Facility code","NAMEFACSCH","SEX","type_vaccine","Serotype","GPSC")
  # mx <- read_csv(dates_path, col_names = F) %>% ncol()
  # datesDF <- read_csv(dates_path, col_names = T) %>%
  #   # dplyr::select(all_of(mx),1,2,3,4)
  #   dplyr::select(any_of(vec_incl_filt))
  
  
  # names(datesDF) <- c("sampleID","TakenDate","Var_00","Facility_name","SeX","Vaccine_type","Serotype","GPSC")
  names(datesDF) <- c("sampleID","TakenDate","Var_00","Var_01","Ward")
  
  # datesDF$Ward <- dplyr::coalesce(datesDF$Ward,"UNKNOWN")
  # datesDF$WardType <- dplyr::coalesce(datesDF$WardType,"UNKNOWN")
  
  # Fix dates - optional 2024-05-29
  # 
  if(lubri_fmt == "mdy"){
    datesDF$TakenDate <- lubridate::mdy(datesDF$TakenDate)
  }else if(lubri_fmt == "dmy"){
    datesDF$TakenDate <- lubridate::dmy(datesDF$TakenDate)
  }else if((lubri_fmt == "ydm")){
    datesDF$TakenDate <- lubridate::ydm(datesDF$TakenDate)
  }else{
    datesDF$TakenDate <- lubridate::ymd(datesDF$TakenDate)
  }
  #  to add more options
  # as.integer(difftime(max(datesDF$TakenDate),min(datesDF$TakenDate),  unit="days"))
  
  
  epiwkDF <- datesDF %>% 
    # dplyr::select(1,4,5) %>% 
    mutate(epiyear = lubridate::epiyear(datesDF$TakenDate)) %>%
    mutate(epiweek = lubridate::epiweek(datesDF$TakenDate)) %>%
    mutate(epiyearweek = paste(epiyear,epiweek,sep=".")) #%>% print()
  
  
  # datesJoin <- epiwkDF %>%  # To make more generic 20240529
  #   group_by(WardType)
  
  
  # Get snp pairwise data
  if(! file.exists(filepath)){
    # next
    stop("Provide path to the SNP distance matrix in the config file")
  }
  
  snpdist_df <- read_csv(file = filepath)
  names(snpdist_df)[1] <- "names"
  snpdist_df$names <- as.character(snpdist_df$names)
  
  
  # mat <- as.matrix(column_to_rownames(snpdist_df,var="names"))
  # str(mat)
  # 
  # # scale data
  # mat <- scale(mat)


  # Get MLST profiles
  if(str_detect(mlst_profile,"\\.csv")){
    mlst <- read_csv(mlst_profile, col_names = T) %>%
      dplyr::select(FILE,ST) %>%
      dplyr::mutate(FILE=as.character(FILE))
  }else{
    mlst <- read_excel(mlst_profile) %>%
      dplyr::select(FILE,ST)%>%
      dplyr::mutate(FILE=as.character(FILE))
  }
    
  
  # # Molten file
  # mdf <- read_csv(moltenpath, col_names = F) #%>%
  #   #dplyr::filter(X3<= snpco & X1 != X2) %>% print(n=40)
  
  
  # Add ST info to the metadata data frame
  incl_vec <- c("sampleID","Hospital","Var_00","Var_01","WardType","ST")
  # incl_vec <- c(1,2,3,4)
  
  annotDF01 <- epiwkDF %>%
    left_join(mlst,by=c("sampleID"="FILE")) %>%
    mutate(ST=as.numeric(as.character(ST))) %>%
    dplyr::select(any_of(incl_vec)) %>%
    tibble::column_to_rownames(var = "sampleID")
  
  # mutate(Clusters=as.numeric(as.character(Clusters))) %>%
  # dplyr::rename("SNP.Epi.Clusters" = Clusters)

  annotDF_subset <- annotDF01
  
  ##### RUN CORE SNP ANALYSIS HERE - 2023-09-27 ###############################
  # Run core analysis per facility/community/hospital
  
  if(any(colnames(epiwkDF)=="Var_01")){
    datesJoin <- epiwkDF %>% group_by(Var_01)
    fc_colnames <- c("sampleID","TakenDate","Var_00","Var_01"
                     ,"ST","epiyear","epiweek","epiyearweek")
  }else{
    datesJoin <- epiwkDF
    fc_colnames <- c("sampleID","TakenDate","Var_00","ST"
                     ,"epiyear","epiweek","epiyearweek")
  }
  
  
  
  fc_df <- datesJoin %>%
    left_join(mlst,by=c("sampleID"="FILE")) %>%
    mutate(ST=as.numeric(as.character(ST))) %>%
    dplyr::select(any_of(fc_colnames))

  main_var <- "Var_00" #names(fc_df)[3]
  facility_vec <- fc_df[[main_var]] %>% unique()
  
  # fc_df %>% group_by(Var_00) %>% summarise(count=n())
  
  for(i in seq_along(facility_vec)){
    fc_val <- facility_vec[[i]]
    
    print(paste0("Starting analysis of: ",fc_val))
    
    fc_df_01 <- fc_df %>%
      dplyr::filter(get({{main_var}}) %in% fc_val)
    
    fc_size <- fc_df_01 %>% group_by(Var_00)%>% summarise(cnt=n()) %>% pull(cnt)
    
    if(fc_size <= 2){
      next
    }
    # set work directory for SNP-EPI clusters
    fc_dir <- str_replace_all(fc_val," ","_")
    work_dir <- file.path(out_dir,fc_dir,paste0("cluster-analysis","_SNPcutoff",snpco,"_Days",daysco))
    if (! dir.exists(work_dir)){
      dir.create(work_dir,recursive = T)
    }
    setwd(work_dir)
    
    
    # Convert SNP-dist df to matrix
    fc_samples_vec <- fc_df_01$sampleID %>% unique()
    fc_df_02 <- snpdist_df %>% 
      dplyr::select("names",any_of(fc_samples_vec)) %>%
      dplyr::filter(names %in% fc_samples_vec)
    
    mat <- as.matrix(column_to_rownames(fc_df_02,var="names"))
    # str(mat)
    
    # scale data
    mat <- scale(mat)
    
    annotDF <- annotDF01[order(match(rownames(annotDF01), colnames(mat))), ,drop = FALSE]
    
    # molten SNP dist
    mdf <- tidyr::pivot_longer(fc_df_02,-names,names_to = "X2", values_to = "X3")
    names(mdf) <- c("X1","X2","X3")
    
    # Enhanced k-means clustering ---------------------------------------------
    # Step 01: Get K
    # 2023-06-03: Now using silhouette optimal cluster calculation based on kmeans
    kmx <- nrow(mat) - 1
    
    if(nrow(mat)>10){
      p4 <- fviz_nbclust(mat, kmeans , method= 'silhouette',nboot = 500)
    }else{
      p4 <- fviz_nbclust(mat, kmeans, k.max = kmx, method= 'silhouette',nboot = 500)
    }
    
    # get optimal number of clusters based on silhoute score
    sigN <- p4$data %>% filter(y==max(p4$data$y)) %>% pull(clusters)
    sigN <- as.numeric(as.character(sigN))
    
    # Use cliuster from hclust
    if(exists("res.hc")){
      sigNN <- length(res.hc$silinfo$clus.avg.widths[(res.hc$silinfo$clus.avg.widths>0.5)])
      sigNN<-sigNN+2
    }else{
      if(sigN>=2){
        sigNN <- sigN+2
      }else{
        sigNN <- 2
      }
      
    }
    
    
    # save graph to file
    p4 <- p4 +
      labs(title = paste(fc_val,"Optimal number of clusters",sigN, sep = ": "))
    ggsave(file.path(work_dir,paste0(fc_val,".","Gap_statistic.",date_var,".png")),p4,width = 10, height = 8)
    
    
    # get core SNP clusters ---------------------------------------------------
    # datesJoin <- fc_df_01
    if(any(colnames(fc_df_01)=="Var_01")){
      datesJoin <- fc_df_01 %>% group_by(Var_01)
    }else{
      datesJoin <- fc_df_01
    }
    
    core_snp_results_list <- run_core_snp_cluster_analysis()
    
    if(length(core_snp_results_list) ==1 && is.na(core_snp_results_list)){ 
      snpClust <- NULL
      snpClustID <- NULL
      df2mat <- NULL
      vec_keep_names <- NULL
      mat_optimal_centres <- NULL
      next
    }
    
    snpClust <- core_snp_results_list[[1]]
    snpClustID <- core_snp_results_list[[2]] 
    df2mat <- core_snp_results_list[[3]]
    vec_keep_names <- core_snp_results_list[[4]] 
    mat_optimal_centres <- core_snp_results_list[[5]]
    
    
    # calculate SNP-Epi clusters - new approach as suggested by Lili-CHARM -------
    
    excl_vec <- c("TakenDate","Date2","epicumsum","CG","num","name","cluster","km_cluster",
                  "Hospital","Ward","Var_01")
    
    
    clusterSet3 <- calculate_SNP_Epi_clusters(snpClust=snpClust,epiwkDF=fc_df_01,mdf=mdf,excl_vec=excl_vec)
    
    
    # datesClusterDF <- datesJoin %>%
    #   dplyr::select(sampleID,TakenDate) %>%
    #   arrange(TakenDate) %>%
    #   mutate(ID2 = dplyr::lead(sampleID,1)) %>%
    #   mutate(Date2 = dplyr::lead(TakenDate,1)) %>%
    #   dplyr::select(sampleID,ID2,TakenDate,Date2,everything()) %>%
    #   mutate(Days = as.numeric(difftime(Date2,TakenDate,units = "days"))) %>%
    #   mutate(epicumsum = if_else(Days <= daysco,1,0)) %>%
    #   mutate(epicumsum = if_else(is.na(epicumsum),0,epicumsum)) %>%
    #   ungroup() %>%
    #   mutate(CG = data.table::rleid(epicumsum)) %>%
    #   dplyr::filter(epicumsum != 0)
    
    
    # Create scatter plots ----------------------------------------------------
    
    if(!is.null(clusterSet3)){
      #stop("Provide path to the SNP distance matrix in the config file")
    
    
    scatter_plots_cmb_list <- list()
    scatter_plots_cmb_list <- create_scatter_plots(datesJoin,clusterSet3,mlst,transmission_type="facility")
    
    plotDF <- scatter_plots_cmb_list[[1]]
    p1 <- scatter_plots_cmb_list[[2]]
    # save scatter plots to file
    ggsave(file.path(work_dir,paste0(fc_val,".scatterplot.",date_var,".png")),
           p1,
           width = 10, 
           height = 8)
    
    
    
    incl_vec2 <- c("sampleID", "Var_00" ,"Hospital","Ward","Var_01","TakenDate","Epiweek","Days","SNPs", "Clusters","ST")

    metadf <- plotDF %>%
      # metadf <- px1 %>%
      ungroup() %>%
      dplyr::select(-epiweek) %>%
      # add_row(sampleID="reference",Var_00=NA,Ward=NA,Var_01=NA,TakenDate=as.Date(refDate),Epiweek=NA,ST=as.factor(refST)) %>%
      dplyr::select(any_of(incl_vec2)) %>% distinct(sampleID,.keep_all = TRUE) %>%
      column_to_rownames(var = "sampleID")
    
    }
    # metadf <- plotDF %>%
    #   # metadf <- px1 %>%
    #   ungroup() %>%
    #   column_to_rownames(var = "sampleID")
    
    # Heatmap graph -----------------------------------------------------------
    
    
    # if(nrow(df2mat) == 0) { next }
    
    if(nrow(df2mat) != 0) {
      mat_snp <-as.matrix(df2mat)
      mat_core <- mat[rownames(mat) %in% vec_keep_names,colnames(mat)%in% vec_keep_names]
      
      # order annotDF_subset based on mat colnames
      annotDF_subset <- annotDF %>% dplyr::filter(get({{main_var}}) == fc_val)
      # annotDF_subset <- annotDF_subset %>% tibble::column_to_rownames(var="sampleID")
      annotDF_subset3 <- annotDF_subset[rownames(annotDF_subset) %in% vec_keep_names,]
      annotDF_subset3 <- annotDF_subset3[order(match(rownames(annotDF_subset3), colnames(mat_snp))), ,drop = FALSE]
    
    }else{
      next
    }
      annotDF_subset2 <- annotDF_subset3 %>%
        mutate(name=rownames(.)) %>%
        mutate(ST=as.numeric(as.character(ST))) %>%
        inner_join(snpClust, by="name") %>%
        dplyr::rename("Core.SNP.clusters"=cluster)
      
      
      
      maxP <-   c(
        RColorBrewer::brewer.pal(12,'Paired'),
        RColorBrewer::brewer.pal(12,'Set3')
      )
      
      
      maxX <-   c(
        RColorBrewer::brewer.pal(9,'Set1'),
        RColorBrewer::brewer.pal(8,'Dark2'),
        RColorBrewer::brewer.pal(2,"Accent")
      )
      
      # vec_incl_filt <- c("Hospital","WardType","ST","Core.SNP.clusters","SNP.Epi.Clusters")
      # vec_excl_filt <- c("km_cluster","Days","SNPs", "Cluster_Cases_count")
      # vec_excl_filt <- c("TakenDate","epiyear","epiweek", "epiyearweek","km_cluster",
                         # "Serotype", "GPSC","Days","SNPs", "Cluster_Cases_count")
      
      vec_incl_filt_hm <- c("Var_00","Hospital","Var_01","ST",
                            "name","Core.SNP.clusters","SNP.Epi.Clusters")
      
      # if(nrow(clusterSet3) != 0){
      if(!is.null(clusterSet3)){
        episnpclust <- clusterSet3 %>%
          mutate(Clusters=as.numeric(as.character(Clusters))) %>%
          # mutate(ST=as.numeric(as.character(ST))) #%>%
          dplyr::rename("SNP.Epi.Clusters" = Clusters)
        
        annotDF_subset2 <- annotDF_subset2 %>%
          left_join(episnpclust,by=c("name"="sampleID")) %>%
          dplyr::select(any_of(vec_incl_filt_hm)) #%>%
          # tibble::column_to_rownames(var = "name")
        
        
        n5 <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% dplyr::n_distinct()
        col5 <- maxX[1:n5]#brewer.pal(n5,"Paired")
        names(col5) <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% unique()
        
        
      }
      
    #   annotDF_subset2 <- annotDF_subset2 %>%
    #     column_to_rownames(var = "name") %>%
    #     dplyr::select(! any_of(vec_excl_filt)) 
    #   
    # }else{
    #   next
    # }
      
      trans_type <- str_to_sentence(transmission_type)
      annotDF_subset2 <- annotDF_subset2 %>%
        column_to_rownames(var = "name") %>%
        dplyr::rename(!!trans_type := Var_00)
      
        # dplyr::select(! any_of(vec_excl_filt))
    
    
    # # display.brewer.all()
    # # order annotDf based on mat colnames
    # maxP <-   c(
    #   RColorBrewer::brewer.pal(12,'Paired'),
    #   RColorBrewer::brewer.pal(12,'Set3')
    # )
    # 
    # # maxP <-   c(
    # #   RColorBrewer::brewer.pal(11,'RdYlBu'),
    # #   RColorBrewer::brewer.pal(11,'BrBG')
    # # )
    # 
    # 
    # maxX <-   c(
    #   RColorBrewer::brewer.pal(9,'Set1'),
    #   RColorBrewer::brewer.pal(8,'Dark2')
    # )
    
      # if(nrow(clusterSet3) != 0){
    n1 <- annotDF_subset2 %>% pull(trans_type) %>% dplyr::n_distinct()
    col1 <- brewer.pal(n1,"Set3")
    names(col1) <- annotDF_subset2 %>% pull(trans_type) %>% unique()
    
    
    if(any(names(annotDF_subset2) == "Var_01")){
      n2 <- annotDF_subset2 %>% pull(Var_01) %>% dplyr::n_distinct()
      col2 <- maxX[1:n2]  #brewer.pal(n2,"YlGnBu")
      names(col2) <- annotDF_subset2 %>% pull(Var_01) %>% unique()
    }
   
    
    n3 <- annotDF_subset2 %>% pull(ST) %>% dplyr::n_distinct()
    col3 <- maxP[1:n3] #brewer.pal(n3,"Accent")
    names(col3) <- annotDF_subset2 %>% pull(ST) %>% unique()
    
    
    n4 <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% dplyr::n_distinct()
    col4 <- maxP[1:n4] #brewer.pal(n4,"Accent")
    names(col4) <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% unique()
    
    
    # n5 <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% dplyr::n_distinct()
    # col5 <- brewer.pal(n5,"Paired")
    # names(col5) <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% unique()
    
    
    
    if (any(names(annotDF_subset2) == "SNP.Epi.Clusters")){
      if(any(names(annotDF_subset2) == "Var_01")){
        annotDF_subset2 <- annotDF_subset2 %>% dplyr::rename("WardType" = Var_01)
        colAnnot2 <- ComplexHeatmap::HeatmapAnnotation(
          df = annotDF_subset2, annotation_height = 7,
          annotation_name_gp = gpar(fontsize = 7),
          col = list(trans_type = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
                     WardType = col2[! is.na(names(col2))],
                     ST = col3[! is.na(names(col3))],
                     Core.SNP.clusters = col4[! is.na(names(col4))],
                     SNP.Epi.Clusters = col5[! is.na(names(col5))]
                     # Ward = annotDF$Ward
          )
          ,border = TRUE
          ,simple_anno_size = unit(0.4, "cm")
        )
      }else{
        colAnnot2 <- ComplexHeatmap::HeatmapAnnotation(
          df = annotDF_subset2, annotation_height = 7,
          annotation_name_gp = gpar(fontsize = 7),
          col = list(trans_type = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
                     # Vaccine_type = col2[! is.na(names(col2))],
                     ST = col3[! is.na(names(col3))],
                     Core.SNP.clusters = col4[! is.na(names(col4))],
                     SNP.Epi.Clusters = col5[! is.na(names(col5))]
                     # Ward = annotDF$Ward
          )
          ,border = TRUE
          ,simple_anno_size = unit(0.4, "cm")
        )
      }
    
    }else{
      colAnnot2 <- ComplexHeatmap::HeatmapAnnotation(
        df = annotDF_subset2, annotation_height = 7,
        annotation_name_gp = gpar(fontsize = 7),
        col = list(Hospital = col1[! is.na(names(col1))], 
                   WardType = col2[! is.na(names(col2))],
                   ST = col3[! is.na(names(col3))],
                   Core.SNP.clusters = col4[! is.na(names(col4))]
                   # SNP.Epi.Clusters = col5[! is.na(names(col5))]
                   # Ward = annotDF$Ward
        )
        ,border = TRUE
        ,simple_anno_size = unit(0.4, "cm")
      )
    }
    
    
    
    if(ncol(mat_optimal_centres) > 2){
      check_n <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% n_distinct() 
      # check_n <- as.numeric(check_n)
      if(check_n > 1){
        nSplit <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>%  dplyr::n_distinct()
        fa = annotDF_subset2 %>% pull(Core.SNP.clusters) #%>% unique()
        dend2 = cluster_within_group(mat_optimal_centres, fa)
      }else{
        dend2=F
        nSplit=NULL 
      }
      
    }else{
      nSplit=NULL 
      dend2=F 
    }
    
    
    p61 <- ComplexHeatmap::Heatmap(as.matrix(mat_snp), column_title = NULL,
                                   # na_col = "black",
                                   # border_gp = gpar(col = "grey", lty = 2),
                                   row_title_gp = gpar(fontsize = 7),
                                   column_title_gp = gpar(fontsize = 7),
                                   top_annotation = colAnnot2, name = "scaled(Means)",
                                   clustering_method_rows = "ward.D2",
                                   clustering_method_columns = "ward.D2",
                                   cluster_columns = dend2, column_split = nSplit,
                                   show_row_names = F,
                                   show_column_dend = F,
                                   show_row_dend = T,
                                   row_dend_side = "right",row_dend_width = unit(0.5, "cm"),
                                   # column_dend_side = "bottom",
                                   column_dend_height = unit(0.5, "cm"),
                                   show_column_names = T,
                                   row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7)
                                   
    )
    
    tidyHeatmap::save_pdf(p61,file.path(work_dir,paste0(fc_val,".coreSNP-clusters_heatmap.",date_var,".pdf")),width = 10, height = 5, units = "in")
    # centers (i.e., the average of each variable for each cluster):
    
    annotDF_subset2 <- NULL
    
    # Write data to file ------------------------------------------------------
    # write EPI+SNP cluster data to file
    # if(nrow(clusterSet3) != 0){
    if(!is.null(clusterSet3)){
      snpepiDF <- clusterSet3 %>% dplyr::select(4,3,2,1)
      openxlsx::write.xlsx(snpepiDF, 
                           file.path(work_dir,
                                     paste(fc_val,"Clusters-SNPs-Epi-Definition",date_var,"xlsx",sep = ".")),
                           overwrite = T)
    }
    
    
    
    # write MLST profile to file
    incl_vec2 <- c("sampleID","Var_00","Hospital","Ward","WardType","TakenDate","Epiweek","ST","Clusters")
    outDF <- plotDF %>% 
      dplyr::select(-ST)%>%
      inner_join(mlst,by=c("sampleID"="FILE")) %>% dplyr::select(any_of(incl_vec2))
    # names(mlst) <- c("sampleid","ST")
    openxlsx::write.xlsx(outDF, 
                         file.path(work_dir,
                                   paste(fc_val,"MLST_profile",date_var,"xlsx",sep = ".")),
                         overwrite = T)
    
    
    
    
    # write Core SNP data to file
    if(nrow(snpClust) != 0){
      snpClusterDF <- snpClust
      vec_cluster_new <- snpClusterDF %>% pull(cluster) %>% unique()
      clust_list <- list()
      names_list <- list()
      for (n in seq_along(vec_cluster_new)){
        var1 <- vec_cluster_new[[n]]
        # clust1 <- names(snpClusterDF$cluster[snpClusterDF$cluster == n])
        clust1 <- snpClusterDF %>% 
          dplyr::filter(cluster == var1) %>%
          pull(name)
        
        names_list[[n]] <- clust1
        
        
        clust_list[[n]] <- snpClusterDF %>% 
          dplyr::filter(cluster == var1) %>%
          dplyr::select(1,2) %>%
          dplyr::rename("sampleID"=name,
                        "coreSNPcluster"=cluster)
      }
      
      
      # Create cluster file to  use for epicurve
      clusterDF <- bind_rows(clust_list)
      # write combinded clusters to file
      openxlsx::write.xlsx(clusterDF, 
                           file.path(work_dir,
                                     paste(fc_val,"Core-SNP-Clusters",date_var,"xlsx",sep = ".")),
                           overwrite = T)
      
      # write core snp clusters to individual sheets
      # silDF <- data.frame(Silhouette_score = res_list[[2]]$silinfo$avg.width) #res.km$silinfo$avg.width)
      # write_csv(silDF,file.path(work_dir,paste0("Silhouette_score",".csv")),col_names = T)
      
      clstrOut <- openxlsx::createWorkbook()
      style <- createStyle(textDecoration = "bold") #, bgFill="Green",fontColour="Black")
      
      for (i in seq_along(clust_list)){
        sheetNm <- paste0("cluster_",i)
        openxlsx::addWorksheet(clstrOut, sheetName = sheetNm)
        abun_data <- clust_list[[i]]
        
        
        if(length(abun_data)==0){
          next
        }else{
          freezePane(clstrOut, sheet = sheetNm, firstRow = T)
          openxlsx::writeDataTable(clstrOut, colNames = T, sheet = sheetNm, headerStyle = style,
                                   x=abun_data)
        }
        
      }
      
      # write to file
      openxlsx::saveWorkbook(clstrOut, 
                             file.path(work_dir,
                                       paste(fc_val,"Samples-Clusters",date_var,"xlsx",sep = ".")),
                             overwrite = T)
      
      
      
      for (i in seq_along(names_list)){
        clstrSamples <- data.frame(names=names_list[[i]])
        
        write_csv(clstrSamples,file.path(work_dir,paste0("cluster_",i,".csv")),col_names = F)
        
      }
    }
    

    # Transmission network analysis -------------------------------------------

    # create outdir
    path_dir <- file.path(dirname(work_dir),paste0("transmission-analysis","_SNPcutoff",snpco,"_Days",daysco))
      #paste0(path,"_SNPcutoff",snpco,"_Days",daysco)
    if (! dir.exists(path_dir)){
      dir.create(path_dir,recursive = T)
    }
    
    setwd(path_dir)
    
    # clusterDF # Core-SNP clusters
    # datesJoin # Dates and other metadata
    
    date_file <- datesJoin
    names(date_file)[1] <- "ID"
    names(date_file)[2] <- "Collectiondate"
    
    date_file <- date_file %>% arrange(Collectiondate)
    
    # Get Core SNP clusters
    if(nrow(clusterDF) != 0){
      coreSNPs <- clusterDF
    }
    
    # Get SNP-EPI clusters 
    # if(nrow(clusterSet3) != 0){
    if(!is.null(clusterSet3)){
      clsters <- clusterSet3 %>%
        dplyr::select(4,1) %>%
        dplyr::rename("names"=sampleID) %>%
        dplyr::rename("cluster"=Clusters)
    }
    
    
    # Get SNP matrix
    snpDistMat <- snpdist_df 
    snpDistMat$names <- as.character(snpDistMat$names)
    # snpDistMat <- snpDistMat %>% filter(names != "reference") %>% dplyr::select(-reference)
    snpDistMat <- as.matrix(column_to_rownames(snpDistMat,var="names"))
    
    
    # subset(snpDistMat,snpDistMat[,1]<100)
    
    names(date_file) <- str_replace_all(names(date_file)," ","_")
    date_file$Collectiondate <- ymd(date_file$Collectiondate)
    dates <- as.Date(date_file$Collectiondate)
    range(dates)
    days <- as.integer(difftime(dates, min(dates), unit="days"))
    
    
    
    # Transmission tree reconstruction using SeqTrack -------------------------
    
    # Perform transmission network analysis per ST
    vec_sts <- date_file %>% dplyr::filter(! is.na(ST)) %>% pull(ST) %>% unique()
    countm <- 1
    for(st in seq_along(vec_sts)){
      m <- vec_sts[st]
      m_ids <- date_file %>% filter(ST==m) %>% pull(ID) %>% unique()
      
      if(length(m_ids) < 3){
        next
      }
      aln_names <- date_file %>% 
        dplyr::filter(ID %in% m_ids) %>%
        arrange(Collectiondate) %>% pull(ID)
      
      
      
      row_idx <- match(aln_names,rownames(snpDistMat))
      col_idx <- match(aln_names,colnames(snpDistMat))
      m_mat <- snpDistMat[row_idx,col_idx]
      
      if(! all(aln_names %in% colnames(m_mat))){
        idx_vec <- which(!aln_names %in% colnames(m_mat))
        aln_names <- aln_names[! seq_along(aln_names) %in% idx_vec]
      }
      
      dates_vec <- date_file %>% 
        dplyr::filter(ID %in% aln_names) %>%
        arrange(Collectiondate) %>% 
        mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
        pull(Collectiondate)
      
      # run seqTrack
      m_mat <- m_mat[,! is.na(colnames(m_mat))]
      m_mat <- m_mat[! is.na(rownames(m_mat)),]
      res <- seqTrack(m_mat, aln_names, dates_vec) #, best="min", annot=T)
      
      res$ances[is.na(res$ances)] <- res$id[which(is.na(res$ances))]
      # rownames_to_column(res,var = "name1")
      # rownames(res)[order(match(res$id,res$ances))]
      res$name1 <- rownames(res)
      res$name2 <- rownames(res)[res$ances]
      
      # add SNPs (pairwise snps)
      # res <- res %>% inner_join(snpDist, by=(c("name1" = "X1","name2" ="X2") ))
      # res<-column_to_rownames(res,"name1")
      # names(res)[names(res)=="X3"] <- "SNPs"
      names(res)[names(res)=="weight"] <- "SNPs"
      res$color <- ifelse(res$SNPs>snpco,"grey","green")
      # # Only include SNPs that meet the cut off in snpco
      # res <- filter(res,SNPs <= 100)
      
      #create edges df
      res1 <-  res %>% dplyr::filter(! rownames(res) == "reference") # Exclude reference 
      edges <- res1 %>% 
        dplyr::select(ances, id, SNPs, everything()) %>%
        dplyr::rename(from=ances) %>%
        dplyr::rename(to=id)
      
      edges <- mutate(edges, width = SNPs/5 + 1)
      edges <- edges %>% select(-width)
      # label <- edges$SNPs
      # edges$label <- label
      edges$label <- edges$SNPs
      edges$label <- as.character(edges$label)
      # edges$from[is.na(edges$from)] <- 29
      edges$from <- as.integer(edges$from)
      
      
      
      
      # date_file <- date_file %>% #group_by(X3) %>% dplyr::count()
      #   mutate(shape=case_when(X3 == "neonatal" ~ "circle",
      #                          X3 == "other" ~ "box",
      #                          TRUE ~ "database")) #%>% print(n=40)
      
      
      # create nodes
      nodes <- res %>%
        #rowid_to_column("test") %>%
        rownames_to_column("label") %>%
        dplyr::select(id, label)
      
      # # add MLST profile to NODES
      # nodes <- nodes %>% 
      #   left_join(mlst, by=c("label" = "FILE" )) #%>%
      # # dplyr::select(-SCHEME) 
      
      
      # names(nodes)[which(names(nodes) == "ST")] <- "group"
      
      len <- nrow(res1)
      
      # nodes <- select(nodes,-shape)
      # nodes <- data.frame(nodes, shape=c(rep('circle',len),'ellipse'))
      
      # add MLST profile to NODES
      nodes <- left_join(nodes,date_file,by=c("label"="ID")) %>%
        dplyr::distinct(id, .keep_all = TRUE)
      
      # Set IDs not in a cluster to background color ----------------------------
      
      # if(file.exists(cluster_file)){
      if(nrow(coreSNPs) != 0){
        ntwrk_out_list <- list()
        c_nodes <- nodes %>%
          left_join(coreSNPs, by=c("label" = "sampleID" ))
        
        # Core-SNP clusters
        names(c_nodes)[which(names(c_nodes) == "coreSNPcluster")] <- "vari"
        
        # run function that generates network graph
        if(all(is.na(c_nodes$vari))){
          next
        }
        ntwrk_out_list <- drw_network(nodes=c_nodes,edges,b=countm,m,type="Core")
        ntwrk <- ntwrk_out_list[[1]]
        # save network graph to file
        ntwrk %>% visSave(file = paste(path_dir,
                                       paste0(fc_val,".CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".html"),
                                       sep = "/"))
        
        
        fnodes <- ntwrk_out_list[[2]]
        write.xlsx(fnodes, paste(path_dir,paste0("CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)
        
        # set color number
        countl <- unique(fnodes$group)[order(unique(fnodes$group))] != "bg"
        countn <- sum(as.numeric(countl))
        countm <- countm+countn
        ntwrk <- NULL
        fnodes <- NULL
      }
      
      
      
      
      # SNP-Epi Clusters
      # if(file.exists(cluster_file_epi_snps)){
      if(nrow(clsters) != 0){
        ntwrk_out_list <- list()
        e_nodes <- nodes %>%
          left_join(clsters, by=c("label" = "names" ))
        
        names(e_nodes)[which(names(e_nodes) == "cluster")] <- "vari"
        
        # run function that generates network graph
        if(all(is.na(e_nodes$vari))){
          next
        }
        ntwrk_out_list <- drw_network(nodes=e_nodes,edges,b=countm,m,type="Epi")
        ntwrk <- ntwrk_out_list[[1]]
        # save network graph to file
        ntwrk %>% visSave(file = paste(path_dir,
                                       paste0("SNP-Epi-clusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".html"),
                                       sep = "/"))
        
        # # write data to file
        # nodes$X1 <- dplyr::coalesce(nodes$X1,"UNKNOWN")
        # nodes$X2 <- dplyr::coalesce(nodes$X2,"UNKNOWN")
        # nodes$X3 <- dplyr::coalesce(nodes$X3,"UNKNOWN")
        fnodes <- ntwrk_out_list[[2]]
        write.xlsx(fnodes, paste(path_dir,paste0("SNP-Epi-clusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)
        
        # set color number
        countl <- unique(fnodes$group)[order(unique(fnodes$group))] != "bg"
        countn <- sum(as.numeric(countl))
        countm <- countm+countn
        ntwrk <- NULL
        fnodes <- NULL
      }
      
    }
    
    
  
}
  


    # # Minimum spanning tree ---------------------------------------------------
    # 
    # library(ape)
    # library(visNetwork)
    # library(networkD3)
    # library(igraph)
    # 
    # mst_dir <- file.path(dirname(out_dir),paste0("minimum-spanning-tree","_SNPcutoff",snpco,"_Days",daysco))
    # 
    # if(! dir.exists(mst_dir)){
    #   dir.create(path = mst_dir, recursive = T)
    # }
    # 
    # # get paths for cluster files
    # cluster_file = file.path(work_dir,paste("Core-SNP-Clusters",date_var,"xlsx",sep = "."))
    # date_file <- date_file
    # core_mat <- snpDistMat
    # 
    # 
    # if(exists('date_file')){
    #   aln_names <- date_file %>% 
    #     arrange(Collectiondate) %>% pull(ID)
    # }else{
    #   stopifnot("MST Mandatory df not available")
    # }
    # 
    # 
    # dates_vec <- date_file %>% 
    #   arrange(Collectiondate) %>% 
    #   mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
    #   pull(Collectiondate)
    # 
    # 
    # 
    # 
    # # Get Core SNP clusters
    # if(file.exists(cluster_file)){
    #   clsters <- read_excel(cluster_file) %>%
    #     dplyr::rename("names"=sampleID,
    #                   "cluster"=coreSNPcluster) %>%
    #     dplyr::filter(names != "reference") %>%
    #     # dplyr::filter(sil_width >= 0.5) %>%
    #     dplyr::select(1,2)
    # }
    # # graph.adjacency depreciated now to use graph_from_adjacency_matrix()
    # mst_out <- ape::mst(core_mat)
    # g <- graph.adjacency(mst_out, mode="undirected", weighted=TRUE)
    # mst_edges <- as_data_frame(g, what="edges")  %>% dplyr::rename("name1"=from,
    #                                                                "name2"=to) 
    # 
    # mst_nodes <- as_data_frame(g, what="vertices") 
    # mst_nodes <- mst_nodes %>% mutate(id=1:nrow(mst_nodes)) %>%
    #   dplyr::select(2,1) %>% 
    #   mutate(id=as.integer(id),
    #          label=name) %>%
    #   inner_join(mlst,by=c("name"="FILE"))%>%
    #   dplyr::rename("group"=ST)
    # 
    # if(exists('clsters') && is.data.frame(get('clsters'))){
    #   mst_nodes <- mst_nodes %>%
    #     dplyr::left_join(clsters,by=c("name"="names")) %>% # 2023-09-30
    #     dplyr::rename("ST"=group,"group"=cluster)
    # }
    # 
    # if(exists('epiwkDF') && is.data.frame(get('epiwkDF'))){
    #   mst_nodes <- mst_nodes %>% 
    #     inner_join(epiwkDF,by=c("name"="sampleID")) %>%
    #     mutate(shape=case_when(WardType == "neonatal" ~ "circle",
    #                            WardType == "other" ~ "box",
    #                            TRUE ~ "diamond"))
    # }else{
    #   mst_nodes <- mst_nodes %>% 
    #     inner_join(datesDF,by=c("name"="sampleID")) 
    # }
    # 
    # 
    # mst_edges <- mst_edges %>%
    #   mutate(from=mst_nodes$id[match(mst_edges$name1,mst_nodes$name)],
    #          to=mst_nodes$id[match(mst_edges$name2,mst_nodes$name)]) %>%
    #   dplyr::select(4,5,3,1,2) %>%
    #   mutate(from=as.integer(from),
    #          to=as.integer(to)) %>%
    #   inner_join(mdf,by=c("name1"="X1","name2"="X2")) %>%
    #   dplyr::rename("label"=X3) %>%
    #   mutate(label=as.character(label))
    # 
    # 
    # 
    # nodes <- mst_nodes
    # edges <- mst_edges
    # 
    # vz <-visNetwork(nodes,edges)
    # visEdges(vz,arrows = NULL,font = list(align="top",size=24))
    # 
    # 
    # # vz %>%
    # #   visEdges(arrows = NULL,font = list(align="top",size=24)) %>%
    # #   visInteraction(navigation = "zoom") %>%
    # #   visInteraction(navigation = "drag") %>%
    # #     visOptions(highlightNearest = TRUE) %>%
    # #   visOptions(collapse = list(enabled = TRUE, keepCoord = TRUE, clusterOptions = list(shape = "circle")))
    # 
    # 
    # 
    # # unique(nodes$group)[order(unique(nodes$group))]
    # groupname=unique(nodes$group)[order(unique(nodes$group))]
    # groupname[is.na(groupname)] <- "bg"
    # nodes$group[is.na(nodes$group)] <- "bg"
    # 
    # clr = c("#ABC2E8","#FFF338","#FFA0A0","#82CD47","#525FE1","#98EECC","#F29727","#D3D04F","#10A19D","#B04759","#D09CFA","#B9E9FC",
    #         "#FFE7A0","#FF6969","#00FFCA","#ECF2FF","#E86A33","#569DAA","#FFE5CA","#FA9884","#A6BB8D","#C8B6A6","#8B1874","#FF78C4")
    # 
    # 
    # 
    # if(any(groupname == "bg")){
    #   groupname<-groupname[groupname != "bg"]
    #   len <- length(groupname) # if reference is included
    # }else{
    #   len <- length(groupname)
    # }
    # 
    # 
    # assign_colors <- data.frame(groupname=groupname, 
    #                             color=clr[1:len])
    # 
    # 
    # # add background color
    # library(glue)
    # bg_row <- c(groupname="bg",color="#F8F4EA")
    # assign_colors <- rbind(assign_colors,bg_row)
    # 
    # gl_code <- list()
    # 
    # for (i in 1:nrow(assign_colors)){
    #   # print(i)
    #   gn <- as.character(assign_colors[i,1])
    #   cl <- as.character(assign_colors[i,2])
    #   gl_code[[i]] <- glue('visGroups(groupname = "{gn}",color = "{cl}", shadow = list(enabled = TRUE))')
    # }
    # 
    # 
    # code_app <- glue(paste(gl_code,collapse = ' %>% '))
    # 
    # vz_groups <- assign_colors %>% pull(1)
    # 
    # ntwk_code <- glue('visNetwork(nodes, edges, height = "1000px", width = "100%",
    #               main = paste0("Minimum spanning tree showing core SNP clusters"),
    #               footer = "*Numbers on edges represent SNP differences between connecting nodes (isolates)") %>%
    #               visNodes(font="12px arial black") %>%
    #               visEdges(arrows = NULL,font = list(align="inside",size=20), label=F) %>%
    #               {code_app} %>%
    #               visLegend(main = "Core.SNP.Clusters",width = 0.1, position = "right") %>%
    #               visClusteringByGroup(groups = vz_groups, label = "cluster : ") %>%
    #               visHierarchicalLayout() %>%
    #               visLayout(randomSeed = 12)')
    # 
    # #font = list(align="inside",size=20)
    # 
    # network <- eval(parse(text=ntwk_code))
    # # Export
    # network %>% visSave(file = paste(mst_dir,
    #                                  paste0("minimum_spanning_tree","_SNPcutoff",snpco,"_Days",daysco,".html"),
    #                                  sep = "/"))
    # 
    # 
    # # Save data files
    # nodesDF <- nodes #%>%
    # #dplyr::rename("ST"=group)
    # 
    # write.xlsx(nodesDF,file = file.path(mst_dir,
    #                                     paste0("minimum_spanning_tree_data","_SNPcutoff",snpco,"_Days",daysco,".xlsx")), overwrite = T)
    # 
    # 
    # 
    # 
    
    
    
    

# END ---------------------------------------------------------------------

  
  
 
 



  
 # # Plot bar plot and scatterplot - epicurves -------------------------------
 #  
 #  library(visdat)
 #  # library(ggnewscale)
 #  # vis_dat(airquality)
 #  vis_dat(plotDF)
 #  
 #  px1 <- plotDF %>% group_by(WardType,Epiweek) %>% dplyr::mutate(count=n())
 #  px1 <- plotDF %>% group_by(Facility_code.x,Epiweek) %>% dplyr::mutate(count=n())
 #  
 #  p1 <- ggplot(px1, 
 #         aes(x = Epiweek,
 #             y = ST, #WardType, #Clusters,
 #             # y = count, #,
 #             color=Clusters,
 #             shape = WardType
 #             )) +
 #    geom_point(alpha = .9, size=4,position=position_jitter(h=0.07,w=0.15)) + 
 #    # scale_fill_brewer(palette = "Dark2") +
 #    labs(title = paste0("Clusters based on Epi definition: SNPs <= ",snpco," + ",daysco,"-day window"),
 #         y = "Sequence Types",
 #         x = "Epiweek") +
 #    theme_bw() +
 #    # scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
 #    theme(panel.grid.major = element_blank(),
 #          panel.grid.minor = element_blank()) +
 #    theme(axis.title = element_text(size = 15)) +
 #    theme(axis.text = element_text(size = 12)) +
 #    theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5, hjust = 1)) 
 #  # +
 #  #   guides(color = guide_colorbar(title = "Horsepower"),
 #  #          shape = guide_legend(title = "Weight"), size = guide_legend(title = "Gear")
 #  #   )
 #   
 #  
 #  ggsave(file.path(work_dir,paste0("Scatterplot.",date_var,".png")),p1,width = 10, height = 8)
 #  
 #  
 #  
 #  # # px <- plotDF %>% 
 #  # #   naniar::bind_shadow()
 #  # px1 <- px1 %>% mutate(ST = if_else(!is.na(Clusters),ST,NA))
 #  # st_shape <- length(levels(px1$ST))
 #  # # jitter <- position_jitter(width = 0.1, height = 0.2)
 #  # p2 <- ggplot(px1,aes(x = Epiweek,
 #  #                     fill = Clusters)) +
 #  #   geom_bar() +
 #  #   # ggnewscale::new_scale_color() +
 #  #   # scale_fill_brewer(palette = "Dark2")
 #  #   geom_point(aes(x=Epiweek,y=count,shape=ST),color="black",size = 4) + #, position=jitter #position=position_jitter(h=0.0,w=0.5)
 #  #   # scale_color_manual(values = c("blue", "green","khaki")) +
 #  #   scale_shape_manual(values=c(seq(2,length.out=st_shape)))+
 #  #   # scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
 #  #   # scale_size_manual(values=c(2,3,4))
 #  #   labs(title = paste0("Clusters based on Epi definition: SNPs <= ",snpco," + ",daysco,"-day window"),
 #  #        y = "Counts",
 #  #        x = "Epiweek") +
 #  #  
 #  #   theme_bw() +
 #  #   scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
 #  #   theme(panel.grid.major = element_blank(),
 #  #         panel.grid.minor = element_blank()) +
 #  #   theme(axis.title = element_text(size = 15)) +
 #  #   theme(axis.text = element_text(size = 12)) +
 #  #   theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5, hjust = 1)) 
 #  # 
 #  # 
 #  # ggsave(file.path(work_dir,paste0("Barplot_epicurve.",date_var,".png")),p2,width = 10, height = 8)
 #  # 
  
  # incl_vec2 <- c("sampleID", "Hospital","Ward","WardType","TakenDate","Epiweek","Days","SNPs", "Clusters","ST")
  # 
  # metadf <- plotDF %>%
  # # metadf <- px1 %>%
  #   ungroup() %>%
  #   dplyr::select(-epiweek) %>%
  #   add_row(sampleID="reference",Hospital=NA,Ward=NA,WardType=NA,TakenDate=as.Date(refDate),Epiweek=NA,ST=as.factor(refST)) %>%
  #   dplyr::select(all_of(incl_vec2)) %>% distinct(sampleID,.keep_all = TRUE) %>%
  #   column_to_rownames(var = "sampleID") 
  # 
  # 
  # 
  
  

   
  
  
  

  
  
  # # Transmission network analysis -------------------------------------------
  # 
  # # # get paths for cluster files
  # cluster_file = file.path(work_dir,paste("Core-SNP-Clusters",date_var,"xlsx",sep = "."))
  # cluster_file_epi_snps = file.path(work_dir,paste("Clusters-SNPs-Epi-Definition",date_var,"xlsx",sep = "."))
  # 
  # # create outdir
  # path_dir <- paste0(path,"_SNPcutoff",snpco,"_Days",daysco)
  # if (! dir.exists(path_dir)){
  #   dir.create(path_dir,recursive = T)
  # }
  # 
  # setwd(path_dir)
  # # create date_file
  # date_file <- read_csv(dates_path, col_names = F) %>% dplyr::select(all_of(mx),1,2,3,4)
  # date_file <- datesDF
  # names(date_file)[1] <- "ID"
  # names(date_file)[6] <- "Collectiondate"
  # 
  # # date_file <- date_file %>% 
  # #   add_row(ID="reference",X2=NA,X3=NA,Collectiondate=as_date("2011-01-01"))
  # date_file <- date_file %>% arrange(Collectiondate)
  # 
  # 
  # # get mlst profile
  # # mlst <- read_excel(mlst_profile) %>% 
  # #   dplyr::select(1,2,3) %>% 
  # #   filter(FILE %in% date_file$ID) #%>% print(n=35)
  # 
  # 
  # mlst <- mlst %>% dplyr::filter(FILE %in% date_file$ID)
  # 
  # # mlst <- mlst %>% add_row(FILE = "reference", SCHEME = "kpneumoniae", ST=as.character(refST))
  # 
  # 
  # # add SNP-EPI clusters file
  # if(file.exists(cluster_file_epi_snps)){
  #   clsters <- read_excel(cluster_file_epi_snps) %>% 
  #     dplyr::filter(sampleID != "reference") %>%
  #     # dplyr::select(1,2) %>%
  #     dplyr::select(4,1) %>%
  #     dplyr::rename("names"=sampleID) %>%
  #     dplyr::rename("cluster"=Clusters)
  # }
  # 
  # 
  # 
  # # Get Core SNP clusters
  # if(file.exists(cluster_file)){
  #   coreSNPs <- read_excel(cluster_file)
  #   # coreSNP_list <- get_core_snp_clusters(m=mat, max=sigNN+2, dates=datesJoin, snpco=snpco, orig=T)
  #   # coreSNPs <- coreSNP_list[[1]] %>% 
  #   #   dplyr::filter(cluster != "reference") %>%
  #   #   dplyr::select(1,2) %>%
  #   #   dplyr::rename("sampleID"=name) %>%
  #   #   dplyr::rename("coreSNPcluster"=cluster)
  # }
  # 
  # 
  # 
  # # get SNP distance matrix (molten format)
  # snpDist <- read_csv(file.path(snp_dist),col_names = F)
  # # snpDist <- filter(snpDist, X3<=100) %>% #X1!=X2)
  # #     dplyr::filter(X1 != "reference") %>%
  # #     dplyr::filter(X2 != "reference")
  # 
  # # snpDist %>% pivot_wider()
  # # 
  # # mdf_vec1 <- snpDist %>% pull(X1) %>% unique()
  # # mdf_vec2 <- snpDist %>% pull(X2) %>% unique()
  # # mdf_vec <- c(mdf_vec1,mdf_vec2) %>% sort() %>% unique()
  # # # # # # filter
  # # mat1 <- snpDistMat[,colnames(snpDistMat) %in% mdf_vec]
  # # mat <- mat1[rownames(mat1) %in% mdf_vec,]
  # 
  # 
  # 
  # # # Get SNP alignment file
  # # dna <- adegenet::fasta2DNAbin(file = file.path(aln_path))
  # 
  # 
  # # get SNP matrix
  # snpDistMat <- read_csv(core_snps)
  # names(snpDistMat)[1] <- "names"
  # snpDistMat$names <- as.character(snpDistMat$names)
  # # snpDistMat <- snpDistMat %>% filter(names != "reference") %>% dplyr::select(-reference)
  # snpDistMat <- as.matrix(column_to_rownames(snpDistMat,var="names"))
  # 
  # 
  # # subset(snpDistMat,snpDistMat[,1]<100)
  # 
  # names(date_file) <- str_replace_all(names(date_file)," ","_")
  # date_file$Collectiondate <- ymd(date_file$Collectiondate)
  # dates <- as.Date(date_file$Collectiondate)
  # range(dates)
  # days <- as.integer(difftime(dates, min(dates), unit="days"))
  # 
  # 
  # # Transmission tree reconstruction using SeqTrack -------------------------
  # 
  # # Perform transmission network analysis per ST
  # vec_sts <- mlst %>% dplyr::filter(ST != "NOVEL") %>% distinct(ST) %>% pull(ST) 
  # countm <- 1
  # for(st in seq_along(vec_sts)){
  #   m <- vec_sts[st]
  #   m_ids <- mlst %>% filter(ST==m) %>% pull(FILE)
  #   
  #   if(length(m_ids) < 3){
  #     next
  #   }
  #   aln_names <- date_file %>% 
  #     dplyr::filter(ID %in% m_ids) %>%
  #     arrange(Collectiondate) %>% pull(ID)
  #   
  #  
  #   
  #   row_idx <- match(aln_names,rownames(snpDistMat))
  #   col_idx <- match(aln_names,colnames(snpDistMat))
  #   m_mat <- snpDistMat[row_idx,col_idx]
  #   
  #   if(! all(aln_names %in% colnames(m_mat))){
  #     idx_vec <- which(!aln_names %in% colnames(m_mat))
  #     aln_names <- aln_names[! seq_along(aln_names) %in% idx_vec]
  #   }
  #     
  #   dates_vec <- date_file %>% 
  #     dplyr::filter(ID %in% aln_names) %>%
  #     arrange(Collectiondate) %>% 
  #     mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
  #     pull(Collectiondate)
  #   
  #   # run seqTrack
  #   m_mat <- m_mat[,! is.na(colnames(m_mat))]
  #   m_mat <- m_mat[! is.na(rownames(m_mat)),]
  #   res <- seqTrack(m_mat, aln_names, dates_vec) #, best="min", annot=T)
  #   
  #   res$ances[is.na(res$ances)] <- res$id[which(is.na(res$ances))]
  #   # rownames_to_column(res,var = "name1")
  #   # rownames(res)[order(match(res$id,res$ances))]
  #   res$name1 <- rownames(res)
  #   res$name2 <- rownames(res)[res$ances]
  #   
  #   # add SNPs (pairwise snps)
  #   # res <- res %>% inner_join(snpDist, by=(c("name1" = "X1","name2" ="X2") ))
  #   # res<-column_to_rownames(res,"name1")
  #   # names(res)[names(res)=="X3"] <- "SNPs"
  #   names(res)[names(res)=="weight"] <- "SNPs"
  #   res$color <- ifelse(res$SNPs>snpco,"grey","green")
  #   # # Only include SNPs that meet the cut off in snpco
  #   # res <- filter(res,SNPs <= 100)
  #   
  #   #create edges df
  #   res1 <-  res %>% dplyr::filter(! rownames(res) == "reference") # Exclude reference 
  #   edges <- res1 %>% 
  #     dplyr::select(ances, id, SNPs, everything()) %>%
  #     dplyr::rename(from=ances) %>%
  #     dplyr::rename(to=id)
  #   
  #   edges <- mutate(edges, width = SNPs/5 + 1)
  #   edges <- edges %>% select(-width)
  #   # label <- edges$SNPs
  #   # edges$label <- label
  #   edges$label <- edges$SNPs
  #   edges$label <- as.character(edges$label)
  #   # edges$from[is.na(edges$from)] <- 29
  #   edges$from <- as.integer(edges$from)
  #   
  #   
  #   
  #   
  #   # date_file <- date_file %>% #group_by(X3) %>% dplyr::count()
  #   #   mutate(shape=case_when(X3 == "neonatal" ~ "circle",
  #   #                          X3 == "other" ~ "box",
  #   #                          TRUE ~ "database")) #%>% print(n=40)
  #   
  #   
  #   # create nodes
  #   nodes <- res %>%
  #     #rowid_to_column("test") %>%
  #     rownames_to_column("label") %>%
  #     dplyr::select(id, label)
  #   
  #   # add MLST profile to NODES
  #   nodes <- nodes %>% 
  #     left_join(mlst, by=c("label" = "FILE" )) #%>%
  #     # dplyr::select(-SCHEME) 
  #   
  #   
  #   # names(nodes)[which(names(nodes) == "ST")] <- "group"
  #   
  #   len <- nrow(res1)
  #   
  #   # nodes <- select(nodes,-shape)
  #   # nodes <- data.frame(nodes, shape=c(rep('circle',len),'ellipse'))
  #   nodes <- left_join(nodes,date_file,by=c("label"="ID")) %>%
  #     dplyr::distinct(id, .keep_all = TRUE)
  # 
  #   # Set IDs not in a cluster to background color ----------------------------
  #   
  #   if(file.exists(cluster_file)){
  #     ntwrk_out_list <- list()
  #     c_nodes <- nodes %>%
  #       left_join(coreSNPs, by=c("label" = "sampleID" ))
  #     
  #     # Core-SNP clusters
  #     names(c_nodes)[which(names(c_nodes) == "coreSNPcluster")] <- "vari"
  #     
  #     # run function that generates network graph
  #     if(all(is.na(c_nodes$vari))){
  #       next
  #     }
  #     ntwrk_out_list <- drw_network(nodes=c_nodes,edges,b=countm,m,type="Core")
  #     ntwrk <- ntwrk_out_list[[1]]
  #     # save network graph to file
  #     ntwrk %>% visSave(file = paste(path_dir,
  #                                    paste0("CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".html"),
  #                                    sep = "/"))
  #     
  # 
  #     fnodes <- ntwrk_out_list[[2]]
  #     write.xlsx(fnodes, paste(path_dir,paste0("CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)
  #     
  #     # set color number
  #     countl <- unique(fnodes$group)[order(unique(fnodes$group))] != "bg"
  #     countn <- sum(as.numeric(countl))
  #     countm <- countm+countn
  #     ntwrk <- NULL
  #     fnodes <- NULL
  #   }
  #   
  #   
  #   
  #   
  #   # SNP-Epi Clusters
  #   if(file.exists(cluster_file_epi_snps)){
  #     ntwrk_out_list <- list()
  #     e_nodes <- nodes %>%
  #       left_join(clsters, by=c("label" = "names" ))
  #     
  #     names(e_nodes)[which(names(e_nodes) == "cluster")] <- "vari"
  #     
  #     # run function that generates network graph
  #     if(all(is.na(e_nodes$vari))){
  #       next
  #     }
  #     ntwrk_out_list <- drw_network(nodes=e_nodes,edges,b=countm,m,type="Epi")
  #     ntwrk <- ntwrk_out_list[[1]]
  #     # save network graph to file
  #     ntwrk %>% visSave(file = paste(path_dir,
  #                                    paste0("SNP-Epi-clusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".html"),
  #                                    sep = "/"))
  #     
  #     # # write data to file
  #     # nodes$X1 <- dplyr::coalesce(nodes$X1,"UNKNOWN")
  #     # nodes$X2 <- dplyr::coalesce(nodes$X2,"UNKNOWN")
  #     # nodes$X3 <- dplyr::coalesce(nodes$X3,"UNKNOWN")
  #     fnodes <- ntwrk_out_list[[2]]
  #     write.xlsx(fnodes, paste(path_dir,paste0("SNP-Epi-clusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)
  #     
  #     # set color number
  #     countl <- unique(fnodes$group)[order(unique(fnodes$group))] != "bg"
  #     countn <- sum(as.numeric(countl))
  #     countm <- countm+countn
  #     ntwrk <- NULL
  #     fnodes <- NULL
  #   }
  #   
  #   
  #   # if(any(names(nodes) == "vari")){
  #   #   names(nodes)[which(names(nodes) == "vari")] <- "group"
  #   # }
  #   # 
  #   # 
  #   # if(all(is.na(nodes$group))){
  #   #   next
  #   # }
  #   # 
  #   # nodes$group[is.na(nodes$group)] <- "bg"
  #   
  #   # # run function that generates network graph
  #   # ntwrk <- drw_network(nodes,edges,b=countm,m,type="Core")
  #   # # save network graph to file
  #   # ntwrk %>% visSave(file = paste(path_dir,
  #   #                                paste0("CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".html"),
  #   #                                sep = "/"))
  #   
  #   # # write data to file
  #   # nodes$X1 <- dplyr::coalesce(nodes$X1,"UNKNOWN")
  #   # nodes$X2 <- dplyr::coalesce(nodes$X2,"UNKNOWN")
  #   # nodes$X3 <- dplyr::coalesce(nodes$X3,"UNKNOWN")
  #   # 
  #   # write.xlsx(nodes, paste(path_dir,paste0("CoreSNPclusters","_ST",m,"_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)
  #   # 
  #   # # set color number
  #   # countl <- unique(nodes$group)[order(unique(nodes$group))] != "bg"
  #   # countn <- sum(as.numeric(countl))
  #   # countm <- countm+countn
  #   # ntwrk <- NULL
  # }
  # 
  # # Minimum spanning tree ---------------------------------------------------
  # 
  # library(ape)
  # library(visNetwork)
  # library(networkD3)
  # library(igraph)
  # 
  # mst_dir <- file.path(dirname(out_dir),paste0("minimum-spanning-tree","_SNPcutoff",snpco,"_Days",daysco))
  # 
  # if(! dir.exists(mst_dir)){
  #   dir.create(path = mst_dir, recursive = T)
  # }
  # 
  # # get paths for cluster files
  # cluster_file = file.path(work_dir,paste("Core-SNP-Clusters",date_var,"xlsx",sep = "."))
  # 
  # # # prepare data matrix
  # # date_file <- read_csv(dates_path, col_names = F) %>% dplyr::select(all_of(mx),1,2,3,4)
  # # names(date_file)[1] <- "ID"
  # # names(date_file)[5] <- "Collectiondate"
  # # date_file <- date_file %>% arrange(Collectiondate)
  # date_file <- date_file
  # 
  # # aln_names <- date_file %>% 
  # #   dplyr::filter(ID %in% m_ids) %>%
  # #   arrange(Collectiondate) %>% pull(ID)
  # 
  # # get SNP matrix
  # # snpDistMat <- read_csv(core_snps)
  # # names(snpDistMat)[1] <- "names"
  # # snpDistMat$names <- as.character(snpDistMat$names)
  # # snpDistMat <- snpDistMat %>% filter(names != "reference") %>% dplyr::select(-reference)
  # # snpDistMat <- as.matrix(column_to_rownames(snpDistMat,var="names"))
  # # # row_idx <- match(aln_names,rownames(snpDistMat))
  # # # col_idx <- match(aln_names,colnames(snpDistMat))
  # # # core_mat <- snpDistMat[row_idx,col_idx]
  # core_mat <- snpDistMat
  # 
  # 
  # if(exists('date_file')){
  #   aln_names <- date_file %>% 
  #     arrange(Collectiondate) %>% pull(ID)
  # }else{
  #   stopifnot("MST Mandatory df not available")
  # }
  # 
  # 
  # dates_vec <- date_file %>% 
  #   arrange(Collectiondate) %>% 
  #   mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
  #   pull(Collectiondate)
  # 
  # 
  # 
  # 
  # # Get Core SNP clusters
  # if(file.exists(cluster_file)){
  #   clsters <- read_excel(cluster_file) %>%
  #     dplyr::rename("names"=sampleID,
  #                   "cluster"=coreSNPcluster) %>%
  #     dplyr::filter(names != "reference") %>%
  #     # dplyr::filter(sil_width >= 0.5) %>%
  #     dplyr::select(1,2)
  # }
  # # graph.adjacency depreciated now to use graph_from_adjacency_matrix()
  # mst_out <- ape::mst(core_mat)
  # g <- graph.adjacency(mst_out, mode="undirected", weighted=TRUE)
  # mst_edges <- as_data_frame(g, what="edges")  %>% dplyr::rename("name1"=from,
  #                                                                "name2"=to) 
  # 
  # mst_nodes <- as_data_frame(g, what="vertices") 
  # mst_nodes <- mst_nodes %>% mutate(id=1:nrow(mst_nodes)) %>%
  #   dplyr::select(2,1) %>% 
  #   mutate(id=as.integer(id),
  #          label=name) %>%
  #   inner_join(mlst,by=c("name"="FILE"))%>%
  #   dplyr::rename("group"=ST)
  # 
  # if(exists('clsters') && is.data.frame(get('clsters'))){
  #   mst_nodes <- mst_nodes %>%
  #     dplyr::left_join(clsters,by=c("name"="names")) %>% # 2023-09-30
  #     dplyr::rename("ST"=group,"group"=cluster)
  # }
  # 
  # if(exists('epiwkDF') && is.data.frame(get('epiwkDF'))){
  #   mst_nodes <- mst_nodes %>% 
  #     inner_join(epiwkDF,by=c("name"="sampleID")) %>%
  #     mutate(shape=case_when(WardType == "neonatal" ~ "circle",
  #                            WardType == "other" ~ "box",
  #                            TRUE ~ "diamond"))
  # }else{
  #   mst_nodes <- mst_nodes %>% 
  #     inner_join(datesDF,by=c("name"="sampleID")) 
  # }
  # 
  # 
  # mst_edges <- mst_edges %>%
  #   mutate(from=mst_nodes$id[match(mst_edges$name1,mst_nodes$name)],
  #          to=mst_nodes$id[match(mst_edges$name2,mst_nodes$name)]) %>%
  #   dplyr::select(4,5,3,1,2) %>%
  #   mutate(from=as.integer(from),
  #          to=as.integer(to)) %>%
  #   inner_join(mdf,by=c("name1"="X1","name2"="X2")) %>%
  #   dplyr::rename("label"=X3) %>%
  #   mutate(label=as.character(label))
  # 
  # 
  # 
  # nodes <- mst_nodes
  # edges <- mst_edges
  # 
  # vz <-visNetwork(nodes,edges)
  # visEdges(vz,arrows = NULL,font = list(align="top",size=24))
  # 
  # 
  # # vz %>%
  # #   visEdges(arrows = NULL,font = list(align="top",size=24)) %>%
  # #   visInteraction(navigation = "zoom") %>%
  # #   visInteraction(navigation = "drag") %>%
  # #     visOptions(highlightNearest = TRUE) %>%
  # #   visOptions(collapse = list(enabled = TRUE, keepCoord = TRUE, clusterOptions = list(shape = "circle")))
  # 
  # 
  # 
  # # unique(nodes$group)[order(unique(nodes$group))]
  # groupname=unique(nodes$group)[order(unique(nodes$group))]
  # groupname[is.na(groupname)] <- "bg"
  # nodes$group[is.na(nodes$group)] <- "bg"
  # 
  # clr = c("#ABC2E8","#FFF338","#FFA0A0","#82CD47","#525FE1","#98EECC","#F29727","#D3D04F","#10A19D","#B04759","#D09CFA","#B9E9FC",
  #         "#FFE7A0","#FF6969","#00FFCA","#ECF2FF","#E86A33","#569DAA","#FFE5CA","#FA9884","#A6BB8D","#C8B6A6","#8B1874","#FF78C4")
  # 
  # 
  # 
  # if(any(groupname == "bg")){
  #   groupname<-groupname[groupname != "bg"]
  #   len <- length(groupname) # if reference is included
  # }else{
  #   len <- length(groupname)
  # }
  # 
  # 
  # assign_colors <- data.frame(groupname=groupname, 
  #                             color=clr[1:len])
  # 
  # 
  # # add background color
  # library(glue)
  # bg_row <- c(groupname="bg",color="#F8F4EA")
  # assign_colors <- rbind(assign_colors,bg_row)
  # 
  # gl_code <- list()
  # 
  # for (i in 1:nrow(assign_colors)){
  #   # print(i)
  #   gn <- as.character(assign_colors[i,1])
  #   cl <- as.character(assign_colors[i,2])
  #   gl_code[[i]] <- glue('visGroups(groupname = "{gn}",color = "{cl}", shadow = list(enabled = TRUE))')
  # }
  # 
  # 
  # code_app <- glue(paste(gl_code,collapse = ' %>% '))
  # 
  # vz_groups <- assign_colors %>% pull(1)
  # 
  # ntwk_code <- glue('visNetwork(nodes, edges, height = "1000px", width = "100%",
  #                 main = paste0("Minimum spanning tree showing core SNP clusters"),
  #                 footer = "*Numbers on edges represent SNP differences between connecting nodes (isolates)") %>%
  #                 visNodes(font="12px arial black") %>%
  #                 visEdges(arrows = NULL,font = list(align="inside",size=20), label=F) %>%
  #                 {code_app} %>%
  #                 visLegend(main = "Core.SNP.Clusters",width = 0.1, position = "right") %>%
  #                 visClusteringByGroup(groups = vz_groups, label = "cluster : ") %>%
  #                 visHierarchicalLayout() %>%
  #                 visLayout(randomSeed = 12)')
  # 
  # #font = list(align="inside",size=20)
  # 
  # network <- eval(parse(text=ntwk_code))
  # # Export
  # network %>% visSave(file = paste(mst_dir,
  #                                  paste0("minimum_spanning_tree","_SNPcutoff",snpco,"_Days",daysco,".html"),
  #                                  sep = "/"))
  # 
  # 
  # # Save data files
  # nodesDF <- nodes #%>%
  #   #dplyr::rename("ST"=group)
  # 
  # write.xlsx(nodesDF,file = file.path(mst_dir,
  #                                     paste0("minimum_spanning_tree_data","_SNPcutoff",snpco,"_Days",daysco,".xlsx")), overwrite = T)
  # 
  
  
# }








# End ---------------------------------------------------------------------


# # main = paste0("Minimum spanning tree(s) showing Core-SNP Clusters:"," SNPcutoff: ","{snpco}")) %>%
# ntwk_code <- glue('visNetwork(nodes, edges, height = "1000px", width = "100%",
#            main = paste0("Minimum spanning tree(s) showing all sequence types in the present study"),footer = "*Numbers on edges represent SNP differences between isolates") %>%
#   visNodes(font="12px arial black") %>%
#   visEdges(arrows = NULL,font = list(align="inside",size=20), label=F) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%
#   visHierarchicalLayout() %>%
#   {code_app} %>%
#   visLegend(main = "Sequence.Type",width = 0.1, position = "right") %>%
#   visLayout(randomSeed = 12) ')
# 
# 
# vz <-visNetwork(nodes,edges)
# visEdges(vz,arrows = NULL,font = list(align="top",size=24))
# 
# 
# vz %>%
#   visEdges(arrows = NULL,font = list(align="top",size=24)) %>%
#   visInteraction(navigation = "zoom") %>%
#   visInteraction(navigation = "drag") %>%
#   visOptions(highlightNearest = TRUE) %>%
#   visOptions(collapse = list(enabled = TRUE, keepCoord = TRUE, clusterOptions = list(shape = "circle")))
# 
# vz %>%
# visGroups(groupname = "-",color = "#ABC2E8", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "1",color = "#FFF338", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "2",color = "#FFA0A0", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "267",color = "#82CD47", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "374",color = "#525FE1", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "424",color = "#98EECC", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "79",color = "#F29727", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "bg",color = "#F8F4EA", shadow = list(enabled = TRUE)) %>%
#   visClusteringByGroup(groups = groupname, label = "ST : ")
# 
# 
# 
# vz %>%
#   visGroups(groupname = "1",color = "#ABC2E8", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "2",color = "#FFF338", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "3",color = "#FFA0A0", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "4",color = "#82CD47", shadow = list(enabled = TRUE)) %>% 
#   visGroups(groupname = "bg",color = "#F8F4EA", shadow = list(enabled = TRUE)) %>%
#   visClusteringByGroup(groups = groupname, label = "ST : ")
# 
# 
# 
# 
# 
# ntwk_code <- glue('visNetwork(nodes, edges, height = "1000px", width = "100%",
#                   main = paste0("Minimum spanning tree(s) showing all sequence types in the present study"),
#                   footer = "*Numbers on edges represent SNP differences between isolates") %>%
#                   visNodes(font="12px arial black") %>%
#                   visEdges(arrows = NULL,font = list(align="inside",size=20), label=F) %>%
#                   {code_app} %>%
#                   visLegend(main = "Core.SNP.Clusters",width = 0.1, position = "right") %>%
#                   visClusteringByGroup(groups = groupname, label = "cluster : ") %>%
#                   visHierarchicalLayout()')
# 
# network <- eval(parse(text=ntwk_code))
# visHierarchicalLayout() %>%
# visLayout(randomSeed = 12) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%

  
  
  
  
  
  
  
  
  
  











# res.km <- eclust(mat, FUNcluster = "kmeans",k.max = sigNN+2,nboot = 500,nstart = 25) #nstart = 25,
# 
# 
# # Get optimal clusters only
# vec_widths <- res.km$silinfo$clus.avg.widths
# vec_optimal <- abs(res.km$silinfo$clus.avg.widths) >= 0.5
# vec_cluster <- res.km$silinfo$widths %>% pull(cluster) %>% unique() %>% as.character()
# names(vec_widths) <- vec_cluster
# vec_keep <- names(vec_widths[vec_optimal])
# vec_excl <- names(vec_widths[!vec_optimal])
# vec_km_cluster <- res.km$cluster
# 
# 
# df_dates <- datesJoin %>% 
#   ungroup() %>%
#   dplyr::select(sampleID,TakenDate)
# 
# 
# 
# # get core SNP clusters -- with total SNPS <= to cut-off
# maxV=0
# list_clusters <-list()
# # keep_df3 <- data.frame()
# for(i in seq_along(vec_keep)){
#   if(is.null(i)){next}
#   
#   grp <- vec_keep[[i]]
#   names_vec <- res.km$silinfo$widths %>% 
#     dplyr::filter(cluster %in% grp) %>%
#     rownames_to_column(var="ID") %>%
#     pull(ID)
#   
#   
#   clster_ids_df <- res.km$silinfo$widths %>% 
#     dplyr::filter(cluster %in% grp) %>%
#     rownames_to_column(var="ID") %>%
#     dplyr::select(1,2) %>%
#     dplyr::rename("km_cluster"=cluster)
# 
#   # ## This is the original logic - which randomly compares isolates in same km cluster 
#   # and creates SNP clusters from least pairwise SNP differences up to cut-off
#   # 2023-09-27
#   
#   # 
#   # comps <- combn(seq_along(names_vec),
#   #                m=2, simplify = F,
#   #                FUN = function(x)names_vec[x])
#   # 
#   # compare_df <- as.data.frame(do.call(rbind, comps))
#   # names(compare_df) <- c("X1","X2")
#   # 
#   # 
#   # check_df1 <- compare_df %>%
#   #   inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
#   #   dplyr::filter(X1 != X2 & X3 <=snpco )
#   
#   
#   # if(nrow(check_df1) > 0){
#   #   keep_df1 <- check_df1 %>%
#   #     arrange(X3) %>%
#   #     mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
#   #     dplyr::filter(snpcumsum != 0) %>%
#   #     mutate(num = as.numeric(paste0(X3,snpcumsum)))
#   # }else{
#   #   next
#   # }
#   
#   check_df1 <- clster_ids_df %>% 
#     dplyr::inner_join(df_dates,by=c("ID"="sampleID")) %>% 
#     dplyr::arrange(TakenDate) %>% 
#     mutate(ID2 = dplyr::lead(ID,1)) %>% 
#     # dplyr::select(-c(2,3)) %>%
#     dplyr::select(ID,ID2,everything()) %>%
#     inner_join(mdf,by=c("ID"="X1","ID2"="X2")) %>%
#     dplyr::filter(ID != ID2 & X3 <=snpco ) %>%
#     dplyr::rename("X1"=ID,"X2"=ID2)
#     
#   
#   if(nrow(check_df1) > 0){
#     keep_df1 <- check_df1 %>%
#       dplyr::arrange(TakenDate) %>%
#       mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
#       dplyr::filter(snpcumsum != 0) %>%
#       mutate(num = as.numeric(paste0(X3,snpcumsum))) %>%
#       dplyr::select(-TakenDate)
#   }else{
#     next
#   }
#   
#   
#   if(nrow(keep_df1)>0){
#     
#     keep_df2 <- keep_df1 %>%
#       ungroup() %>%
#       pivot_longer(cols = c(X1,X2),values_to = "sampleID") %>%
#       inner_join(df_dates,by="sampleID") %>%
#       dplyr::rename("SNPs"=X3) %>%
#       distinct(sampleID, .keep_all = TRUE) %>%
#       # mutate(snpcumsum2 = cumsum_group(SNPs, snpco)) %>%
#       # dplyr::filter(snpcumsum2 != 0) %>%
#       # dplyr::rename("clst"=snpcumsum2) 
#       dplyr::rename("clst"=snpcumsum)
#     
#     if(length(keep_df2$SNPs[keep_df2$SNPs == 0])>=3){
#       
#       # Use option below when using TakenDates - 2023-09-27
#       # df_filtx <- keep_df2 %>% 
#       #   arrange(SNPs,TakenDate) %>%
#       #   mutate(ID2=lead(sampleID)) %>%
#       #   mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
#       #   left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>% 
#       #   mutate(rn=row_number())
#       
#       
#       df_filtx <- keep_df2 %>% 
#         arrange(TakenDate) %>%
#         mutate(ID2=lead(sampleID)) %>%
#         mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
#         left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>% 
#         mutate(rn=row_number())
#       
#       # filter 1
#       list_checks <- list()
#       list_checks <-connect_samples(df_filtx)
#       
#       # filter 2
#       df_filty <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
#       
#       
#       df_filty <- df_filty %>%
#         mutate(snpcumsum3 = cumsum_group(as.numeric(V3), snpco)) 
#       
#       
#       grps_vec <- df_filty %>% pull(snpcumsum3) %>% unique()
#       
#       list_dfs <- list()
#       for(z in seq_along(grps_vec)){
#         
#         if(z==1){
#           k=z
#         }else{
#           k=c(1:z)
#         }
#         
#         
#         if(z==1){
#           grp_df <- df_filty %>% 
#             dplyr::filter(snpcumsum3 %in% k)
#           list_dfs[[k]] <- grp_df
#           nm_vec <- unique(c(grp_df$V1,grp_df$V2))
#         }
#         
#         
#         j=z+1
#         
#         keep <- df_filty %>% 
#           dplyr::filter(snpcumsum3 == j) %>%
#           dplyr::filter(! V1 %in% nm_vec) %>%
#           dplyr::filter(! V2 %in% nm_vec)
#         
#         if(nrow(keep) != 0){
#           list_dfs[[j]] <- keep
#         }
#         
#         
#         bnd_df <- bind_rows(list_dfs)
#         nm_vec <- unique(c(bnd_df$V1,bnd_df$V2))
#         
#         if(j==length(grps_vec)){
#           break
#         }
#       }
#       
#       
#       df_filtz <- bind_rows(list_dfs) %>%
#         mutate(Clusters = data.table::rleid(snpcumsum3)+maxV) %>%
#         dplyr::rename("sampleID"=V1,
#                       "ID2"=V2,
#                       "SNPs"=V3,
#                       "clst"=snpcumsum3) %>%
#         mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
#         mutate(across(all_of(c("SNPs","clst")), as.numeric))
#       
#       
#       # get SNP pairs excluded in by above filter 1
#       a_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V1) %>% pull(sampleID)
#       b_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V2) %>% pull(sampleID)
#       c_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V1) %>% pull(ID2)
#       d_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V2) %>% pull(ID2)
#       
#       m_vec <- c(a_vec,b_vec,c_vec,d_vec) %>% unique()
#       len<-abs(length(m_vec))
#       if(len%%2==0 & len !=0){
#         list_comp <- combn(seq_along(m_vec),
#                            m=2, simplify = F,
#                            FUN = function(x)m_vec[x])
#         
#         df_comps <- as.data.frame(do.call(rbind, list_comp))
#         names(df_comps) <- c("X1","X2")
#         
#         maxV=max(df_filtz$Clusters)   
#         
#         df_snps <- df_comps %>%
#           inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
#           dplyr::filter(X1 != X2 & X3 <snpco ) %>%
#           arrange(X3) %>%
#           mutate(snpcumsum4 = cumsum_group(as.numeric(X3), snpco)) %>%
#           dplyr::filter(snpcumsum4 != 0) %>%
#           mutate(Clusters = data.table::rleid(snpcumsum4)+maxV) %>%
#           dplyr::rename("sampleID"=X1,
#                         "ID2"=X2,
#                         "SNPs"=X3,
#                         "clst"=snpcumsum4) %>%
#           mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
#           mutate(across(all_of(c("SNPs","clst")), as.numeric))
#       }else{
#         df_snps<-NULL
#       }
#       
#       
#       
#       maxV=0
#       
#       keep_df2 <- dplyr::bind_rows(list(df_filtz,df_snps)) %>%
#         pivot_longer(cols = c(sampleID,ID2),values_to = "sampleID") %>%
#         dplyr::group_by(Clusters) %>%
#         dplyr::distinct(sampleID, .keep_all = TRUE) %>% 
#         ungroup()
#       
#       # maxV=max(keep_df3$Clusters) 
#       
#     }
#     
#     if(any(names(keep_df2) == "Clusters")){
#       # if(! is.null(keep_df3)){
#       keep_df <- keep_df2
#     }else{
#       keep_df <- keep_df2 %>%
#         mutate(Clusters = data.table::rleid(clst)+maxV) %>%
#         group_by(Clusters) %>%
#         mutate(cnt=n()) %>%
#         dplyr::filter(cnt > 1) %>%
#         dplyr::select(SNPs,clst,Clusters,name,sampleID) 
#     }
#     
#     
#   }else{
#     next
#   }
#   
#   maxV=max(keep_df$Clusters)
#   list_clusters[[i]] <- keep_df %>% 
#     dplyr::select(sampleID,Clusters) %>%
#     dplyr::rename("name"=sampleID,
#                   "cluster"=Clusters) %>%
#     inner_join(clster_ids_df,by=c("name"="ID"))
#   
#   # keep_df3=NULL
#   keep_df2=NULL
#   keep_df=NULL
#   
# }


# snpClust<-NULL
# snpClustID<-NULL
# 
# if(length(vec_keep) == 0){
#   next
# }
# 
# 
# snpClust <- bind_rows(list_clusters) %>%
#   mutate(km_cluster=as.character(km_cluster),
#          cluster=as.character(cluster)) 



# 
# 
# 
# 
# daysDF <-  epiwkDF %>% 
#   dplyr::select(sampleID,TakenDate) %>%
#   arrange(TakenDate) %>%
#   mutate(ID2 = dplyr::lead(sampleID,1)) %>% 
#   mutate(Date2 = dplyr::lead(TakenDate,1)) %>% 
#   dplyr::select(sampleID,ID2,TakenDate,Date2,everything()) %>%
#   mutate(Days = as.numeric(difftime(Date2,TakenDate,units = "days")))
# 
# 
# 
# 
# # datesClusterDF <- epiwkDF %>% 
# #   dplyr::select(sampleID,TakenDate) %>%
# #   arrange(TakenDate) %>%
# #   mutate(ID2 = dplyr::lead(sampleID,1)) %>% 
# #   mutate(Date2 = dplyr::lead(TakenDate,1)) %>% 
# #   dplyr::select(sampleID,ID2,TakenDate,Date2,everything()) %>%
# #   mutate(Days = as.numeric(difftime(Date2,TakenDate,units = "days"))) %>%
# #   mutate(epicumsum = cumsum_group(Days, daysco)) %>%
# #   ungroup() %>%
# #   mutate(CG = data.table::rleid(epicumsum)) #%>% print(n=40)
# 
# 
# if(nrow(datesClusterDF) == 0){
#   next
# }
# # Get clusters meeting SNP and Epi definition
# excl_vec <- c("TakenDate","Date2","epicumsum","CG","snpcumsum","num","name")
# 
# 
# clusterSet3 <- datesClusterDF %>%
#   inner_join(mdf,by=c("sampleID"="X1","ID2"="X2")) 
# 
# if(nrow(clusterSet3) == 0){ next }
# 
# clusterSet3 <- clusterSet3 %>% 
#   mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
#   dplyr::filter(snpcumsum != 0) %>%
#   mutate(num = as.numeric(paste0(CG,snpcumsum))) %>%
#   mutate(Clusters = data.table::rleid(num)) %>%
#   ungroup() %>%
#   pivot_longer(cols = c(sampleID,ID2),values_to = "sampleID") %>%
#   select(! all_of(excl_vec)) %>%
#   group_by(Clusters) %>%
#   dplyr::rename("SNPs"=X3) %>%
#   # dplyr::select(3,4) %>%
#   distinct(sampleID, .keep_all = TRUE) %>%
#   dplyr::mutate(n_clusters=n()) %>%
#   mutate(Clusters=as.factor(as.character(Clusters))) %>%
#   dplyr::rename("Cluster_Cases_count"=n_clusters) #%>% print(n=40)
# 
# 
# if(nrow(clusterSet3) == 0){
#   next
# }
# 
# 
# datesJoin <- epiwkDF %>%
#   group_by(WardType) #%>% print(n=40)


  
  # # Core-SNP-Analysis -------------------------------------------------------
  # 
  # # Enhanced hierarchical clustering
  # # Color labels using k-means clusters
  # 
  # # res.hc <- eclust(mat, "hclust", scale="none", nboot = 500) # compute hclust
  # # 
  # # # dendrogam
  # # (p3 <- fviz_dend(res.hc,
  # #                  k_colors=NULL,
  # #                  rect = T,
  # #                  # k=sigN,
  # #                  # label_cols = km.clust[res.hc$order],
  # #                  cex = 0.6))
  # # 
  # # # to relook: https://stackoverflow.com/questions/38862303/customize-ggplot2-axis-labels-with-different-colors
  # # # p3 +  theme(legend.position=c(.9, .6))
  # # 
  # # p3 <- p3 +
  # #   labs(title = "Cluster dendrogram",
  # #        subtitle = "Core-SNP distance clustering",
  # #        #caption ="Label colors showing kmeans cluster groups",
  # #        y = "Height",
  # #        x = "") +
  # #   theme(axis.title = element_text(size = 15)) +
  # #   theme(axis.text = element_text(size = 12))
  # # 
  # # ggsave(file.path(work_dir,paste0("Cluster_dendrogram.",date_var,".png")),p3,width = 10, height = 8)
  # 
  # 
  # # display.brewer.all()
  # 
  # # scale SNP differences
  # # mat <- scale(mat)
  # 
  # incl_vec <- c("Hospital","WardType","Clusters","ST")
  # annotDF <- metadf %>%
  #   dplyr::select(all_of(incl_vec)) %>%
  #   mutate(Clusters=as.numeric(as.character(Clusters))) %>%
  #   mutate(ST=as.numeric(as.character(ST))) %>%
  #   dplyr::rename("SNP.Epi.Clusters" = Clusters)
  #   # dplyr::rename("EpiClusters" = Clusters)
  #   # mutate(Clusters = if_else(is.na(as.numeric(Clusters)),0,as.numeric(as.character(Clusters)))) %>% #na_if(as.double(Clusters),0)) %>%
  # 
  # annotDF_subset <- annotDF 
  # 
  # 
  # # order annotDf based on mat colnames
  # maxP <-   c(
  #   RColorBrewer::brewer.pal(12,'Paired'),
  #   RColorBrewer::brewer.pal(12,'Set3')
  # )
  # 
  # 
  # maxX <-   c(
  #   RColorBrewer::brewer.pal(9,'Set1'),
  #   RColorBrewer::brewer.pal(8,'Dark2')
  # )
  # 
  # 
  # annotDF <- annotDF[order(match(rownames(annotDF), colnames(mat))), ,drop = FALSE]
  # # 
  # 
  # # n1 <- annotDF %>% pull(Hospital) %>% dplyr::n_distinct()
  # # col1 <- maxP[1:n1] #brewer.pal(n1,"Set3")
  # # names(col1) <- annotDF %>% pull(Hospital) %>% unique()
  # # 
  # # 
  # # n2 <- annotDF %>% pull(WardType) %>% dplyr::n_distinct()
  # # col2 <- maxP[1:n2] #brewer.pal(n2,"Paired")
  # # names(col2) <- annotDF %>% pull(WardType) %>% unique()
  # # 
  # # 
  # # n3 <- annotDF %>% pull(SNP.Epi.Clusters) %>% dplyr::n_distinct()
  # # col3 <- maxX[1:n3] #brewer.pal(n3,"Paired")
  # # names(col3) <- annotDF %>% pull(SNP.Epi.Clusters) %>% unique()
  # # 
  # # 
  # # n4 <- annotDF %>% pull(ST) %>% dplyr::n_distinct()
  # # col4 <- maxP[1:n4] #brewer.pal(n4,"Accent")
  # # names(col4) <- annotDF %>% pull(ST) %>% unique()
  # # 
  # # 
  # # 
  # # colAnnot <- ComplexHeatmap::HeatmapAnnotation(
  # #   df = annotDF, annotation_height = 8,
  # #   annotation_name_gp = gpar(fontsize = 8),
  # #   # Hospital = annotDF$Hospital,
  # #   # Ward = annotDF$Ward,
  # #   # EpiClusters = annotDF$clusterGrp,
  # #   col = list(Hospital = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
  # #              WardType = col2[! is.na(names(col2))],
  # #              SNP.Epi.Clusters = col3[! is.na(names(col3))],
  # #              ST = col4[! is.na(names(col4))]
  # #              # Ward = annotDF$Ward
  # #   )
  # #   ,border = TRUE
  # #   ,simple_anno_size = unit(0.5, "cm")
  # # )
  # # 
  # 
  # 
  # # Enhanced k-means clustering
  # # 2023-06-03: Now using silhouette optimal cluster calculation based on kmeans
  # kmx <- nrow(mat) - 1
  # 
  # if(nrow(mat)>10){
  #   p4 <- fviz_nbclust(mat, kmeans , method= 'silhouette',nboot = 500)
  # }else{
  #   p4 <- fviz_nbclust(mat, kmeans, k.max = kmx, method= 'silhouette',nboot = 500)
  # }
  # 
  # # get optimal number of clusters based on silhoute score
  # sigN <- p4$data %>% filter(y==max(p4$data$y)) %>% pull(clusters)
  # sigN <- as.numeric(as.character(sigN))
  # 
  # # Use cliuster from hclust
  # if(exists("res.hc")){
  #   sigNN <- length(res.hc$silinfo$clus.avg.widths[(res.hc$silinfo$clus.avg.widths>0.5)])
  # }else{
  #   sigNN <- sigN
  # }
  # 
  # 
  # # save graph to file
  # p4 <- p4 +
  #   labs(title = paste("Optimal number of clusters",sigN, sep = ": "))
  # ggsave(file.path(work_dir,paste0("Gap_statistic.",date_var,".png")),p4,width = 10, height = 8)
  # 
  # 
  # res.km <- eclust(mat, FUNcluster = "kmeans",k.max = sigNN+2,nboot = 500,nstart = 25) #nstart = 25,
  # # p6 <- ComplexHeatmap::Heatmap(as.matrix(res.km$centers), 
  # #                         top_annotation = colAnnot, name = "ClusterMeans")
  # # 
  # # tidyHeatmap::save_pdf(p6,file.path(work_dir,paste0("cluster_heatmap.",date_var,".pdf")),width = 10, height = 5, units = "in")
  # 
  # 
  # # Get optimal clusters only
  # vec_widths <- res.km$silinfo$clus.avg.widths
  # vec_optimal <- abs(res.km$silinfo$clus.avg.widths) >= 0.5
  # vec_cluster <- res.km$silinfo$widths %>% pull(cluster) %>% unique() %>% as.character()
  # names(vec_widths) <- vec_cluster
  # vec_keep <- names(vec_widths[vec_optimal])
  # vec_excl <- names(vec_widths[!vec_optimal])
  # vec_km_cluster <- res.km$cluster
  # 
  # 
  # df_dates <- datesJoin %>% 
  #   ungroup() %>%
  #   dplyr::select(sampleID,TakenDate)
  # 
  # 
  # 
  # maxV=0
  # list_clusters <-list()
  # # keep_df3 <- data.frame()
  # for(i in seq_along(vec_keep)){
  #   if(is.null(i)){next}
  #   
  #   grp <- vec_keep[[i]]
  #   names_vec <- res.km$silinfo$widths %>% 
  #     dplyr::filter(cluster %in% grp) %>%
  #     rownames_to_column(var="ID") %>%
  #     pull(ID)
  #   
  #   
  #   clster_ids_df <- res.km$silinfo$widths %>% 
  #     dplyr::filter(cluster %in% grp) %>%
  #     rownames_to_column(var="ID") %>%
  #     dplyr::select(1,2) %>%
  #     dplyr::rename("km_cluster"=cluster)
  # 
  #   comps <- combn(seq_along(names_vec),
  #                  m=2, simplify = F,
  #                  FUN = function(x)names_vec[x])
  #   
  #   compare_df <- as.data.frame(do.call(rbind, comps))
  #   names(compare_df) <- c("X1","X2")
  #   
  #   
  #   check_df1 <- compare_df %>%
  #     inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
  #     dplyr::filter(X1 != X2 & X3 <=snpco ) 
  #   
  #   if(nrow(check_df1) > 0){
  #     keep_df1 <- check_df1 %>%
  #       arrange(X3) %>%
  #       mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
  #       dplyr::filter(snpcumsum != 0) %>%
  #       mutate(num = as.numeric(paste0(X3,snpcumsum)))
  #   }else{
  #     next
  #   }
  #   
  #   
  #   if(nrow(keep_df1)>0){
  # 
  #     keep_df2 <- keep_df1 %>%
  #       ungroup() %>%
  #       pivot_longer(cols = c(X1,X2),values_to = "sampleID") %>%
  #       inner_join(df_dates,by="sampleID") %>%
  #       dplyr::rename("SNPs"=X3) %>%
  #       distinct(sampleID, .keep_all = TRUE) %>%
  #       # mutate(snpcumsum2 = cumsum_group(SNPs, snpco)) %>%
  #       # dplyr::filter(snpcumsum2 != 0) %>%
  #       # dplyr::rename("clst"=snpcumsum2) 
  #       dplyr::rename("clst"=snpcumsum)
  #       
  #     if(length(keep_df2$SNPs[keep_df2$SNPs == 0])>=3){
  #     
  #       df_filtx <- keep_df2 %>% 
  #         arrange(SNPs,TakenDate) %>%
  #         mutate(ID2=lead(sampleID)) %>%
  #         mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
  #         left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>% 
  #         mutate(rn=row_number())
  #       
  #       # filter 1
  #       list_checks <- list()
  #       list_checks <-connect_samples(df_filtx)
  #       
  #       # filter 2
  #       df_filty <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
  #       
  #       
  #       df_filty <- df_filty %>%
  #         mutate(snpcumsum3 = cumsum_group(as.numeric(V3), snpco)) 
  #       
  #       
  #       grps_vec <- df_filty %>% pull(snpcumsum3) %>% unique()
  #       
  #       list_dfs <- list()
  #       for(z in seq_along(grps_vec)){
  #         
  #         if(z==1){
  #           k=z
  #         }else{
  #           k=c(1:z)
  #         }
  #         
  #         
  #         if(z==1){
  #           grp_df <- df_filty %>% 
  #             dplyr::filter(snpcumsum3 %in% k)
  #           list_dfs[[k]] <- grp_df
  #           nm_vec <- unique(c(grp_df$V1,grp_df$V2))
  #         }
  #         
  #           
  #         j=z+1
  #      
  #         keep <- df_filty %>% 
  #           dplyr::filter(snpcumsum3 == j) %>%
  #           dplyr::filter(! V1 %in% nm_vec) %>%
  #           dplyr::filter(! V2 %in% nm_vec)
  #         
  #         if(nrow(keep) != 0){
  #           list_dfs[[j]] <- keep
  #         }
  #         
  #         
  #         bnd_df <- bind_rows(list_dfs)
  #         nm_vec <- unique(c(bnd_df$V1,bnd_df$V2))
  #         
  #         if(j==length(grps_vec)){
  #           break
  #         }
  #       }
  #       
  #       
  #       df_filtz <- bind_rows(list_dfs) %>%
  #         mutate(Clusters = data.table::rleid(snpcumsum3)+maxV) %>%
  #         dplyr::rename("sampleID"=V1,
  #                       "ID2"=V2,
  #                       "SNPs"=V3,
  #                       "clst"=snpcumsum3) %>%
  #         mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
  #         mutate(across(all_of(c("SNPs","clst")), as.numeric))
  #         
  #       
  #       # get SNP pairs excluded in by above filter 1
  #       a_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V1) %>% pull(sampleID)
  #       b_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V2) %>% pull(sampleID)
  #       c_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V1) %>% pull(ID2)
  #       d_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V2) %>% pull(ID2)
  #       
  #       m_vec <- c(a_vec,b_vec,c_vec,d_vec) %>% unique()
  #       len<-abs(length(m_vec))
  #       if(len%%2==0 & len !=0){
  #         list_comp <- combn(seq_along(m_vec),
  #                            m=2, simplify = F,
  #                            FUN = function(x)m_vec[x])
  #         
  #         df_comps <- as.data.frame(do.call(rbind, list_comp))
  #         names(df_comps) <- c("X1","X2")
  #         
  #         maxV=max(df_filtz$Clusters)   
  #         
  #         df_snps <- df_comps %>%
  #           inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
  #           dplyr::filter(X1 != X2 & X3 <snpco ) %>%
  #           arrange(X3) %>%
  #           mutate(snpcumsum4 = cumsum_group(as.numeric(X3), snpco)) %>%
  #           dplyr::filter(snpcumsum4 != 0) %>%
  #           mutate(Clusters = data.table::rleid(snpcumsum4)+maxV) %>%
  #           dplyr::rename("sampleID"=X1,
  #                         "ID2"=X2,
  #                         "SNPs"=X3,
  #                         "clst"=snpcumsum4) %>%
  #           mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
  #           mutate(across(all_of(c("SNPs","clst")), as.numeric))
  #       }else{
  #         df_snps<-NULL
  #       }
  #       
  #     
  #        
  #       maxV=0
  #       
  #       keep_df2 <- dplyr::bind_rows(list(df_filtz,df_snps)) %>%
  #         pivot_longer(cols = c(sampleID,ID2),values_to = "sampleID") %>%
  #         dplyr::group_by(Clusters) %>%
  #         dplyr::distinct(sampleID, .keep_all = TRUE) %>% 
  #         ungroup()
  #       
  #       # maxV=max(keep_df3$Clusters) 
  #       
  #     }
  #     
  #     if(any(names(keep_df2) == "Clusters")){
  #     # if(! is.null(keep_df3)){
  #       keep_df <- keep_df2
  #     }else{
  #       keep_df <- keep_df2 %>%
  #         mutate(Clusters = data.table::rleid(clst)+maxV) %>%
  #         group_by(Clusters) %>%
  #         mutate(cnt=n()) %>%
  #         dplyr::filter(cnt > 1) %>%
  #         dplyr::select(SNPs,clst,Clusters,name,sampleID) 
  #     }
  #    
  #     
  #   }else{
  #     next
  #   }
  #   
  #   maxV=max(keep_df$Clusters)
  #   list_clusters[[i]] <- keep_df %>% 
  #     dplyr::select(sampleID,Clusters) %>%
  #     dplyr::rename("name"=sampleID,
  #                   "cluster"=Clusters) %>%
  #     inner_join(clster_ids_df,by=c("name"="ID"))
  #   
  #   # keep_df3=NULL
  #   keep_df2=NULL
  #   keep_df=NULL
  #   
  # }
  # 
  # 
  # # if(exists('epiwkDF') && is.data.frame(get('epiwkDF'))){
  # 
  # snpClust<-NULL
  # snpClustID<-NULL
  # 
  # if(length(vec_keep) == 0){
  #   next
  # }
  # 
  # 
  # snpClust <- bind_rows(list_clusters) %>%
  #   mutate(km_cluster=as.character(km_cluster),
  #          cluster=as.character(cluster)) 
  
  
  # snpClustID <- bind_rows(list_clusters) %>%
  #   mutate(km_cluster=as.character(km_cluster),
  #          cluster=as.character(cluster)) %>%
  #   dplyr::select(2,3) %>%
  #   distinct(.,.keep_all = T) 
  
  
  # snpClustID <- snpClust %>%
  #   mutate(km_cluster=as.character(km_cluster),
  #          cluster=as.character(cluster)) %>%
  #   dplyr::select(2,3) %>%
  #   distinct(.,.keep_all = T) 
  # 
  # vec_keep_names <- snpClust %>% pull(name) %>% unique() 
  # # mat_optimal_centres <- res.km$centers[,colnames(res.km$centers) %in% vec_keep_names]
  # centers_mat <- res_list[[2]]$centers
  # mat_optimal_centres <- centers_mat[,colnames(centers_mat) %in% vec_keep_names]
  # # mat_optimal_centres <- mat_optimal_centres[! rownames(mat_optimal_centres) %in% vec_excl,]
  # 
  # 
  # 
  # df2mat <- as.data.frame(mat_optimal_centres) %>%
  #   rownames_to_column(var = "rowID") %>%
  #   inner_join(snpClustID, by=c("rowID"="km_cluster")) %>%
  #   column_to_rownames(var = "cluster") %>%
  #   dplyr::select(-rowID)
  # 
  # if(nrow(df2mat) == 0) { next }
  # 
  # mat_snp <-as.matrix(df2mat)
  # 
  # 
  # 
  # mat_core <- mat[rownames(mat) %in% vec_keep_names,colnames(mat)%in% vec_keep_names]
  # 
  # 
  # 
  # 
  # 
  # # order annotDF_subset based on mat colnames
  # annotDF_subset <- annotDF
  # annotDF_subset <- annotDF_subset %>% tibble::column_to_rownames(var="sampleID")
  # annotDF_subset2 <- annotDF_subset[rownames(annotDF_subset) %in% vec_keep_names,]
  # annotDF_subset2 <- annotDF_subset2[order(match(rownames(annotDF_subset2), colnames(mat_snp))), ,drop = FALSE]
  # 
  # 
  # episnpclust <- clusterSet3 %>%
  #   mutate(Clusters=as.numeric(as.character(Clusters))) %>%
  #   # mutate(ST=as.numeric(as.character(ST))) #%>%
  #   dplyr::rename("SNP.Epi.Clusters" = Clusters)
  # 
  # vec_incl_filt <- c("Hospital","WardType","ST","Core.SNP.clusters","SNP.Epi.Clusters") 
  # 
  # annotDF_subset2 <- annotDF_subset2 %>%
  #   mutate(name=rownames(.)) %>%
  #   mutate(ST=as.numeric(as.character(ST))) %>%
  #   inner_join(snpClust, by="name") %>%
  #   dplyr::rename("Core.SNP.clusters"=cluster) %>%
  #   left_join(episnpclust,by=c("name"="sampleID")) %>%
  #   column_to_rownames(var = "name") %>%
  #   dplyr::select(all_of(vec_incl_filt)) 
  # 
  # 
  # 
  # 
  # 
  # n1 <- annotDF_subset2 %>% pull(Hospital) %>% dplyr::n_distinct()
  # col1 <- brewer.pal(n1,"Set3")
  # names(col1) <- annotDF_subset2 %>% pull(Hospital) %>% unique()
  # 
  # 
  # n2 <- annotDF_subset2 %>% pull(WardType) %>% dplyr::n_distinct()
  # col2 <- brewer.pal(n2,"Paired")
  # names(col2) <- annotDF_subset2 %>% pull(WardType) %>% unique()
  # 
  # 
  # n3 <- annotDF_subset2 %>% pull(ST) %>% dplyr::n_distinct()
  # col3 <- brewer.pal(n3,"Accent")
  # names(col3) <- annotDF_subset2 %>% pull(ST) %>% unique()
  # 
  # 
  # n4 <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% dplyr::n_distinct()
  # col4 <- maxP[1:n4] #brewer.pal(n4,"Accent")
  # names(col4) <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% unique()
  # 
  # 
  # n5 <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% dplyr::n_distinct()
  # col5 <- brewer.pal(n5,"Paired")
  # names(col5) <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% unique()
  # 
  # 
  # 
  # colAnnot2 <- ComplexHeatmap::HeatmapAnnotation(
  #   df = annotDF_subset2, annotation_height = 7,
  #   annotation_name_gp = gpar(fontsize = 7),
  #   # Hospital = annotDF$Hospital,
  #   # Ward = annotDF$Ward,
  #   # EpiClusters = annotDF$clusterGrp,
  #   col = list(Hospital = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
  #              WardType = col2[! is.na(names(col2))],
  #              ST = col3[! is.na(names(col3))],
  #              core.SNP.clusters = col4[! is.na(names(col4))],
  #              SNP.Epi.Clusters = col5[! is.na(names(col5))]
  #              # Ward = annotDF$Ward
  #   )
  #   ,border = TRUE
  #   ,simple_anno_size = unit(0.4, "cm")
  # )
  # 
  # 
  # 
  # if(ncol(mat_optimal_centres) > 2){
  #   check_n <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% n_distinct() 
  #   # check_n <- as.numeric(check_n)
  #   if(check_n > 1){
  #     nSplit <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>%  dplyr::n_distinct()
  #     fa = annotDF_subset2 %>% pull(Core.SNP.clusters) #%>% unique()
  #     dend2 = cluster_within_group(mat_optimal_centres, fa)
  #   }else{
  #     dend2=F
  #     nSplit=NULL 
  #   }
  #   
  # }else{
  #   nSplit=NULL 
  #   dend2=F 
  # }
  # 
  # 
  # p61 <- ComplexHeatmap::Heatmap(as.matrix(mat_snp), column_title = NULL,
  #                                # na_col = "black",
  #                                # border_gp = gpar(col = "grey", lty = 2),
  #                                row_title_gp = gpar(fontsize = 7),
  #                                column_title_gp = gpar(fontsize = 7),
  #                                top_annotation = colAnnot2, name = "scaled(Means)",
  #                         clustering_method_rows = "ward.D2",
  #                         clustering_method_columns = "ward.D2",
  #                         cluster_columns = dend2, column_split = nSplit,
  #                         show_row_names = F,
  #                         show_column_dend = T,
  #                         show_row_dend = T,
  #                         row_dend_side = "right",row_dend_width = unit(0.5, "cm"),
  #                         column_dend_side = "bottom",column_dend_height = unit(0.5, "cm"),
  #                         show_column_names = T,
  #                         row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7)
  #                         
  #                         )
  # 
  # tidyHeatmap::save_pdf(p61,file.path(work_dir,paste0("coreSNP-clusters_heatmap.",date_var,".pdf")),width = 10, height = 5, units = "in")
  # # centers (i.e., the average of each variable for each cluster):
  

  # Get samples names from each cluster -------------------------------------
  # this will be used when we run the annotation step
  # metadf <- pcaMeta
  # snpClusterDF <- res.km$silinfo$widths[rownames(res.km$silinfo$width) %in% vec_keep_names ,] #res.km$silinfo$widths
  



# END ---------------------------------------------------------------------



# aln_names <- date_file %>% arrange(Collectiondate) %>% pull(ID)
# dates_vec <- date_file %>% 
#   arrange(Collectiondate) %>% 
#   mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
#   pull(Collectiondate)
# 
# 
# row_idx <- match(aln_names,rownames(snpDistMat))
# col_idx <- match(aln_names,colnames(snpDistMat))
# snpDistMat <- snpDistMat[row_idx,col_idx]
# 
# 
# # run seqTrack
# res <- seqTrack(snpDistMat, aln_names, dates_vec) #, best="min", annot=T)


# res$ances[is.na(res$ances)] <- res$id[which(is.na(res$ances))]
# # rownames_to_column(res,var = "name1")
# # rownames(res)[order(match(res$id,res$ances))]
# res$name1 <- rownames(res)
# res$name2 <- rownames(res)[res$ances]
# 
# # add SNPs (pairwise snps)
# # res <- res %>% inner_join(snpDist, by=(c("name1" = "X1","name2" ="X2") ))
# # res<-column_to_rownames(res,"name1")
# # names(res)[names(res)=="X3"] <- "SNPs"
# names(res)[names(res)=="weight"] <- "SNPs"
# res$color <- ifelse(res$SNPs>snpco,"grey","green")
# # # Only include SNPs that meet the cut off in snpco
# # res <- filter(res,SNPs <= 100)
# 
# #create edges df
# res1 <-  res %>% dplyr::filter(! rownames(res) == "reference") # Exclude reference 
# edges <- res1 %>% 
#   dplyr::select(ances, id, SNPs, everything()) %>%
#   dplyr::rename(from=ances) %>%
#   dplyr::rename(to=id)
# 
# edges <- mutate(edges, width = SNPs/5 + 1)
# edges <- edges %>% select(-width)
# # label <- edges$SNPs
# # edges$label <- label
# edges$label <- edges$SNPs
# edges$label <- as.character(edges$label)
# # edges$from[is.na(edges$from)] <- 29
# edges$from <- as.integer(edges$from)
# 
# 
# 
# 
# date_file <- date_file %>% #group_by(X3) %>% dplyr::count()
#   mutate(shape=case_when(X3 == "neonatal" ~ "circle",
#                          X3 == "other" ~ "box",
#                          TRUE ~ "database")) #%>% print(n=40)
# 
# # node list ---------------------------------------------------------------
# 
# nodes <- res %>%
#   #rowid_to_column("test") %>%
#   rownames_to_column("label") %>%
#   dplyr::select(id, label)
# 
# # add MLST profile to NODES
# nodes <- nodes %>% 
#   left_join(mlst, by=c("label" = "FILE" )) %>%
#   dplyr::select(-SCHEME) %>%
#   # left_join(epi_snp_clusters, by=c("label" = "sampleID"))
#   left_join(clsters, by=c("label" = "names" ))
# 
# # names(nodes)[which(names(nodes) == "ST")] <- "group"
# names(nodes)[which(names(nodes) == "cluster")] <- "group"
# len <- nrow(res1)
# # nodes <- select(nodes,-shape)
# # nodes <- data.frame(nodes, shape=c(rep('circle',len),'ellipse'))
# nodes <- left_join(nodes,date_file,by=c("label"="ID")) %>%
#   dplyr::distinct(id, .keep_all = TRUE)
# 
# 
# # Set IDs not in a cluster to background color ----------------------------
# 
# nodes$group[is.na(nodes$group)] <- "bg"


# # interactive network graphs
# library(visNetwork)
# library(networkD3)
# # basic plot 1
# vz <-visNetwork(nodes, edges)
# visEdges(vz,arrows = "to",font = list(align="top",size=24))

# unique(nodes$group)[order(unique(nodes$group))]


# network %>% visSave(file = paste(path_dir,
#                                  paste0("transmission_ntwk_snpEpiclusters","_SNPcutoff",snpco,"_Days",daysco,".html"),
#                                  sep = "/"))
# 
# 
# # write data to file
# nodes$X1 <- dplyr::coalesce(nodes$X1,"UNKNOWN")
# nodes$X2 <- dplyr::coalesce(nodes$X2,"UNKNOWN")
# nodes$X3 <- dplyr::coalesce(nodes$X3,"UNKNOWN")
# 
# write.xlsx(nodes, paste(path_dir,paste0("transmission_ntwk_snpEpiclusters","_SNPcutoff",snpco,"_Days",daysco,".xlsx"),sep = "/"), overwrite = T, asTable = T)



# # Transmission network analysis (STs) -------------------------------------
# 
# nodes <- res %>%
#   #rowid_to_column("test") %>%
#   rownames_to_column("label") %>%
#   dplyr::select(id, label)
# 
# # add MLST profile to NODES
# nodes <- nodes %>%
#   left_join(mlst, by=c("label" = "FILE" )) %>%
#   dplyr::select(-SCHEME) #%>%
#   #left_join(clsters, by=c("label" = "names" ))
# 
# names(nodes)[which(names(nodes) == "ST")] <- "group"
# # names(nodes)[which(names(nodes) == "cluster")] <- "group"
# len <- nrow(res1)
# nodes <- left_join(nodes,date_file,by=c("label"="ID"))
# 
# 
# 
# # Set IDs not in a cluster to background color ----------------------------
# 
# nodes$group[is.na(nodes$group)] <- "bg"
# 
# write.xlsx(nodes, paste(path_dir,"transmission_network_metadata_STs.xlsx",sep = "/"), overwrite = T, asTable = T)
# 
# 
# # interactive network graphs
# library(visNetwork)
# library(networkD3)
# # basic plot 1
# vz <-visNetwork(nodes, edges)
# visEdges(vz,arrows = "to",font = list(align="top",size=24))
# 
# 
# 
# 
# unique(nodes$group)[order(unique(nodes$group))]
# 
# groupname=unique(nodes$group)[order(unique(nodes$group))]
# 
# # clr = c("#ABC2E8","#FFF338","#FFA0A0","#C3ACD0","#C68B59","#A6BB8D","#10A19D","brown","#D09CFA","#B9E9FC","#FFE7A0","#FF6969","#00FFCA","#82CD47","#ECF2FF")
# clr = c("#ABC2E8","#FFF338","#FFA0A0","#FF78C4","#525FE1","#98EECC","#F29727","#D3D04F","#10A19D","#B04759","#D09CFA","#B9E9FC",
#         "#FFE7A0","#FF6969","#00FFCA","#ECF2FF","#E86A33","#569DAA","#FFE5CA","#FA9884","#A6BB8D","#C8B6A6","#8B1874","#82CD47")
# 
# groupname=unique(nodes$group)[order(unique(nodes$group))]
# 
# if(any(groupname == "bg")){
#   groupname<-groupname[groupname != "bg"]
#   len <- length(groupname) # if reference is included
# }else{
#   len <- length(groupname)
# }
# 
# 
# assign_colors <- data.frame(groupname=groupname, 
#                             color=clr[1:len]
# )
# # add background color
# library(glue)
# bg_row <- c(groupname="bg",color="#F8F4EA")
# assign_colors <- rbind(assign_colors,bg_row)
# 
# gl_code <- list()
# 
# for (i in 1:nrow(assign_colors)){
#   # print(i)
#   gn <- as.character(assign_colors[i,1])
#   cl <- as.character(assign_colors[i,2])
#   gl_code[[i]] <- glue('visGroups(groupname = "{gn}",color = "{cl}", shadow = list(enabled = TRUE))')
# }
# 
# 
# code_app <- glue(paste(gl_code,collapse = ' %>% '))
# 
# 
# ntwk_code <- glue('visNetwork(nodes, edges, height = "500px", width = "100%",
#   main = paste0("Sequence Types:"," SNPcutoff: ",{snpco})) %>%
#   visEdges(arrows = "to",font = list(align="top",size=14, color=ifelse(edges$weight>11,"grey","green"))) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%
#   visHierarchicalLayout() %>%
#   {code_app} %>%
#   visLegend(main = "Sequence.Types",width = 0.1, position = "right")
#   ')
# 
# 
# network <- eval(parse(text=ntwk_code))
# 
# # Export
# # Use visSave() for save network in .html file, and visExport() to save as .png with shiny :
# network %>% visSave(file = paste(path_dir,
#                                  paste0("transmission_ntwk_STs","_SNPcutoff",snpco,"_Days",daysco,".html"),
#                                  sep = "/")
# )
# 





# # Transmission network analysis (coreSNPs) --------------------------------
# snpob <- 100
# 
# # IDs to incl
# incl_ids <- date_file %>% 
#   arrange(Collectiondate) %>% pull(ID)
# 
# # # Molten file
# # moltdf <- read_csv(moltenpath, col_names = F) %>%
# #   dplyr::filter(X1 %in% incl_ids | X2 %in% incl_ids)
# #   dplyr::filter(X3<= snpob ) %>% 
# #   dplyr::arrange(X1) %>% 
# #   pivot_wider(names_from = "X2", values_from = "X3")
# #   print(n=40)
# 
# # get SNP matrix
# snpDistMat <- read_csv(core_snps)
# names(snpDistMat)[1] <- "names"
# snpDistMat$names <- as.character(snpDistMat$names)
# snpDistMat <- snpDistMat %>% filter(names != "reference") %>% dplyr::select(-reference)
# snpDistMat <- as.matrix(column_to_rownames(snpDistMat,var="names"))
# 
# 
# # subset(snpDistMat,snpDistMat[,1]<100)
# 
# names(date_file) <- str_replace_all(names(date_file)," ","_")
# date_file$Collectiondate <- ymd(date_file$Collectiondate)
# dates <- as.Date(date_file$Collectiondate)
# range(dates)
# days <- as.integer(difftime(dates, min(dates), unit="days"))
# 
# 
# aln_names <- date_file %>% 
#   arrange(Collectiondate) %>% pull(ID)
# 
# dates_vec <- date_file %>% 
#   arrange(Collectiondate) %>% 
#   mutate(Collectiondate=as.POSIXct(Collectiondate)) %>%
#   pull(Collectiondate)
# 
# row_idx <- match(aln_names,rownames(snpDistMat))
# col_idx <- match(aln_names,colnames(snpDistMat))
# core_mat <- snpDistMat[row_idx,col_idx]
# 
# # run seqTrack
# res <- seqTrack(core_mat, aln_names, dates_vec) #, best="min", annot=T)
# 
# 
# 
# res$ances[is.na(res$ances)] <- res$id[which(is.na(res$ances))]
# res$name1 <- rownames(res)
# res$name2 <- rownames(res)[res$ances]
# names(res)[names(res)=="weight"] <- "SNPs"
# res$color <- ifelse(res$SNPs>snpco,"grey","green")
# # res <- res %>% dplyr::filter(SNPs <= 100)
# 
# 
# # newIDs <- res %>% 
# #   dplyr::filter(SNPs <= 100) %>%
# #   rownames_to_column(var = "IDs") %>% 
# #   pull(IDs)
# # 
# # core_mat1 <- core_mat[, match(newIDs, colnames(core_mat))] 
# # core_mat1[match(newIDs, rownames(core_mat1)),]
# 
# 
# # # add clusters file
# clsters <- read_excel(cluster_file) %>%
#   dplyr::rename("names"=sampleID,
#                 "cluster"=coreSNPcluster) %>%
#   dplyr::filter(names != "reference") %>%
#   # dplyr::filter(sil_width >= 0.5) %>%
#   dplyr::select(1,2)
# 
# #create edges df
# res1 <-  res %>% dplyr::filter(! rownames(res) == "reference") # Exclude reference 
# edges <- res1 %>% 
#   dplyr::select(ances, id, SNPs, everything()) %>%
#   dplyr::rename(from=ances) %>%
#   dplyr::rename(to=id)
# 
# edges <- mutate(edges, width = SNPs/5 + 1)
# edges <- edges %>% select(-width)
# # label <- edges$SNPs
# # edges$label <- label
# edges$label <- edges$SNPs
# edges$label <- as.character(edges$label)
# # edges$from[is.na(edges$from)] <- 29
# edges$from <- as.integer(edges$from)
# 
# # create nodes DF
# nodes <- res %>%
#   #rowid_to_column("test") %>%
#   rownames_to_column("label") %>%
#   dplyr::select(id, label)
# 
# # add MLST profile to NODES
# nodes <- nodes %>% 
#   left_join(mlst, by=c("label" = "FILE" )) %>%
#   dplyr::select(-SCHEME) %>%
#   # left_join(epi_snp_clusters, by=c("label" = "sampleID"))
#   left_join(clsters, by=c("label" = "names" ))
# 
# # names(nodes)[which(names(nodes) == "ST")] <- "group"
# names(nodes)[which(names(nodes) == "cluster")] <- "group"
# len <- nrow(res1)
# # nodes <- select(nodes,-shape)
# # nodes <- data.frame(nodes, shape=c(rep('circle',len),'ellipse'))
# nodes <- left_join(nodes,date_file,by=c("label"="ID")) 
# 
# 
# 
# # Set IDs not in a cluster to background color ----------------------------
# 
# nodes$group[is.na(nodes$group)] <- "bg"
# 
# 
# #%>% print(n=40)
# 
# 
# # MLST_PROFILE <-read_excel(mlst_profile) %>%
# #   inner_join(nodes, by=c("FILE" = "label" )) %>%
# #   dplyr::select(1:10)
# write.xlsx(nodes, paste(path_dir,"transmission_network_coreSNPs.xlsx",sep = "/"), overwrite = T, asTable = T)
# 
# 
# # interactive network graphs
# library(visNetwork)
# library(networkD3)
# # basic plot 1
# vz <-visNetwork(nodes, edges)
# visEdges(vz,arrows = "to",font = list(align="top",size=24))
# 
# 
# 
# 
# unique(nodes$group)[order(unique(nodes$group))]
# 
# 
# # clr = c("#ABC2E8","#FFF338","#FFA0A0","#C3ACD0","#C68B59","#A6BB8D","#10A19D","brown","#D09CFA","#B9E9FC","#FFE7A0","#FF6969","#00FFCA","#82CD47","#ECF2FF")
# clr = c("#ABC2E8","#FFF338","#FFA0A0","#FF78C4","#525FE1","#98EECC","#F29727","#D3D04F","#10A19D","#B04759","#D09CFA","#B9E9FC",
#         "#FFE7A0","#FF6969","#00FFCA","#ECF2FF","#E86A33","#569DAA","#FFE5CA","#FA9884","#A6BB8D","#C8B6A6","#8B1874","#82CD47")
# 
# groupname=unique(nodes$group)[order(unique(nodes$group))]
# 
# if(any(groupname == "bg")){
#   groupname<-groupname[groupname != "bg"]
#   len <- length(groupname) # if reference is included
# }else{
#   len <- length(groupname)
# }
# 
# 
# assign_colors <- data.frame(groupname=groupname, 
#                             color=clr[1:len]
# )
# 
# # add background color
# library(glue)
# bg_row <- c(groupname="bg",color="#F8F4EA")
# assign_colors <- rbind(assign_colors,bg_row)
# 
# gl_code <- list()
# 
# for (i in 1:nrow(assign_colors)){
#   # print(i)
#   gn <- as.character(assign_colors[i,1])
#   cl <- as.character(assign_colors[i,2])
#   gl_code[[i]] <- glue('visGroups(groupname = "{gn}",color = "{cl}", shadow = list(enabled = TRUE))')
# }
# 
# 
# code_app <- glue(paste(gl_code,collapse = ' %>% '))
# 
# 
# ntwk_code <- glue('visNetwork(nodes, edges, height = "500px", width = "100%",
#   main = paste0("core SNP clusters:"," SNPcutoff: ",{snpco}),footer = "*bg=background, i.e. cases not included in clusters") %>% 
#   visEdges(arrows = "to",font = list(align="top",size=14, color=ifelse(edges$weight>11,"grey","green"))) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%
#   visHierarchicalLayout() %>%
#   {code_app} %>%
#   visLegend(main = "core.SNP.Clusters",width = 0.1, position = "right")
#   ')
# 
# # visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>% 
# # visHierarchicalLayout(sortMethod = "directed", direction = "UD", levelSeparation = 100) %>% 
# network <- eval(parse(text=ntwk_code))
# 
# 
# # Export
# # Use visSave() for save network in .html file, and visExport() to save as .png with shiny :
# # network <- visNetwork(nodes, edges, width = "100%")
# 
# network %>% visSave(file = paste(path_dir,
#                                  paste0("transmission_ntwk_coreSNPs","_SNPcutoff",snpco,"_Days",daysco,".html"),
#                                  sep = "/")
# )

# }


# visConfigure(graph=network,enabled = TRUE)
