
# Set-up the environment --------------------------------------------------
src_path <- file.path(getwd())
# src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
source(file.path(src_path,"scripts","pacman_install_packages.R"))


date_var <- as.Date(date(), format = "%a %b %d %H:%M:%S %Y")
# date_var <- str_remove_all(date_var,"-")


# Functions ---------------------------------------------------------------

source(file.path(src_path,"functions","load_custom_functions.R"))

# Data load ---------------------------------------------------------------
# data_paths <- file.path(src_path,"conf")
# source(file.path(data_paths,"Richael_s_pneumo.R"))
# source(file.path(data_paths,"BabyGERMS_kpn_temb.R"))
# source(file.path(data_paths,"Thabo_p_aerogenosa.R"))
# START ANALYSIS ----------------------------------------------------------

# Set SNP cut-off and EPI-days
# Run all steps using loop
# for(i in 1:nrow(comparisons)){
  # snpco=comparisons[i,1]
  # daysco=comparisons[i,2]

  # }

  # snpco=20
  # daysco=21
  
  
  # # get metadata first -- 20240713 
  # # The first 3 columns should have the following - in the order specified below:
  # #   1. sample_id
  # #   2. collection_date
  # #   3. facility_name/hospital_name/community codes/regions/metros (Main_var)
  # #   4. Ward_name (or ward type) if 3 is hospital (Var_01)
  # 
  # # get BabyGERMS metadata
  # # metdata now included in the dates file: 2023-05-01
  # mx <- read_csv(dates_path, col_names = F) %>% ncol()
  # datesDF <- read_csv(dates_path, col_names = F) %>% dplyr::select(all_of(mx),4,1,3,2) #,1,2,3,4)
  # names(datesDF) <- c("sampleID","TakenDate","Hospital","WardType","Ward")
  # datesDF$Ward <- dplyr::coalesce(datesDF$Ward,"UNKNOWN")
  # datesDF$WardType <- dplyr::coalesce(datesDF$WardType,"UNKNOWN")
  # 
  # 
  # # # Richael's data
  # # vec_incl_filt <- c("FILE","Sampling_date","Facility code","NAMEFACSCH","SEX","type_vaccine","Serotype","GPSC")
  # # mx <- read_csv(dates_path, col_names = F) %>% ncol()
  # # datesDF <- read_csv(dates_path, col_names = T) %>%
  # #   # dplyr::select(all_of(mx),1,2,3,4)
  # #   dplyr::select(any_of(vec_incl_filt))
  # 
  # 
  # # # Thabo's metadata
  # # # metdata now included in the dates file: 2023-05-01
  # # var_sel <- c("sampleid","biosample","collection_year","collection_date","Country")
  # # # var_sel <- c("biosample","collection_year","Country")
  # # mx <- read_excel(dates_path) %>% dplyr::filter(Host=="Human") %>% ncol()
  # # datesDF <- read_excel(dates_path) %>% 
  # #   dplyr::filter(Host=="Human") %>% 
  # #   dplyr::select(any_of(var_sel)) #%>%
  # #   # dplyr::filter(collection_date != "Missing")
  # # 
  # # datesDF <- dplyr::select(datesDF,sampleid,collection_date,Country,collection_year,biosample)
  # # 
  # # names(datesDF) <- c("sampleID","TakenDate",Var_01,"Year","biosample") #,"Var_01","Ward")
  # # 
  # # 
  # # # names(datesDF) <- c("sampleID","TakenDate","Var_00","Facility_name","SeX","Vaccine_type","Serotype","GPSC")
  # # names(datesDF) <- c("sampleID","TakenDate",Main_var,Var_01,"Ward")
  # # 
  # # # datesDF$Ward <- dplyr::coalesce(datesDF$Ward,"UNKNOWN")
  # # # datesDF$WardType <- dplyr::coalesce(datesDF$WardType,"UNKNOWN")
  # 
  
  
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
  st_remove <- c("ST")
  incl_vec <- c("sampleID",Main_var,Var_01,Var_02,"biosample","ST")
  
  if(clust_type == "Core"){
    epiwkDF <- datesDF
    
    annotDF01 <- epiwkDF %>%
      dplyr::select(!any_of(st_remove)) %>%
      left_join(mlst,by=c("sampleID"="FILE")) %>%
      mutate(ST=as.numeric(as.character(ST))) %>%
      dplyr::select(any_of(incl_vec)) %>%
      tibble::column_to_rownames(var = "sampleID")
    
  }else{
    # Fix dates - optional 2024-05-29
    # 
    if(lubri_fmt == "mdy"){
      datesDF <- datesDF %>%
        dplyr::mutate({{Var_02}}:=lubridate::mdy(.data[[Var_02]]))
    }else if(lubri_fmt == "dmy"){
      datesDF <- datesDF %>%
        dplyr::mutate({{Var_02}}:=lubridate::dmy(.data[[Var_02]]))
    }else if((lubri_fmt == "ydm")){
       datesDF <- datesDF %>%
        dplyr::mutate({{Var_02}}:=lubridate::ydm(.data[[Var_02]]))
    }else{
      datesDF <- datesDF %>%
        dplyr::mutate({{Var_02}}:=lubridate::ymd(.data[[Var_02]]))
    }
    #  to add more options
    # as.integer(difftime(max(datesDF$TakenDate),min(datesDF$TakenDate),  unit="days"))
    
    
    epiwkDF <- datesDF %>% 
      # dplyr::select(1,4,5) %>% 
      mutate(epiyear = lubridate::epiyear(.data[[Var_02]])) %>%
      mutate(epiweek = lubridate::epiweek(.data[[Var_02]])) %>%
      mutate(epiyearweek = paste(epiyear,epiweek,sep=".")) #%>% print()
    
    
    # datesJoin <- epiwkDF %>%  # To make more generic 20240529
    #   group_by(WardType)
    
    # # Add ST info to the metadata data frame
    # incl_vec <- c("sampleID","Hospital",Main_var,"Var_01","WardType","ST")
    # incl_vec <- c(1,2,3,4)
    
    annotDF01 <- epiwkDF %>%
      left_join(mlst,by=c("sampleID"="FILE")) %>%
      mutate(ST=as.numeric(as.character(ST))) %>%
      dplyr::select(any_of(incl_vec)) %>%
      tibble::column_to_rownames(var = "sampleID")
    
    # mutate(Clusters=as.numeric(as.character(Clusters))) %>%
    # dplyr::rename("SNP.Epi.Clusters" = Clusters)
    
  }
  


  annotDF_subset <- annotDF01
  
  ##### RUN CORE SNP ANALYSIS HERE - 2023-09-27 ###############################
  # Run core analysis per facility/community/hospital
  
  if(! exists("Main_var")){
    # next
    stop("Define main variable in the config file")
  }
  
  
  if(any(colnames(epiwkDF)==Var_01)){
    datesJoin <- epiwkDF %>% group_by(.data[[Var_01]])
    fc_colnames <- c("sampleID",Var_02,Main_var,Var_01,
                     "ST","epiyear","epiweek","epiyearweek")
  }else{
    datesJoin <- epiwkDF
    fc_colnames <- c("sampleID",Var_02,Main_var,"ST",
                     "epiyear","epiweek","epiyearweek","biosample")
  }
  
  
  
  fc_df <- datesJoin %>%
    left_join(mlst,by=c("sampleID"="FILE")) %>%
    mutate(ST=as.numeric(as.character(ST))) %>%
    dplyr::select(any_of(fc_colnames))

  main_var <- Main_var #names(fc_df)[3]
  facility_vec <- fc_df[[main_var]] %>% unique()
  
  # fc_df %>% group_by(Var_00) %>% summarise(count=n())
  
  for(i in seq_along(facility_vec)){
    fc_val <- facility_vec[[i]]
    
    print(paste0("Starting analysis of: ",fc_val))
    
    fc_df_01 <- fc_df %>%
      dplyr::filter(get({{main_var}}) %in% fc_val)
    
    fc_size <- fc_df_01 %>% group_by(.data[[Main_var]])%>% summarise(cnt=n()) %>% pull(cnt)
    
    if(fc_size <= 2){
      next
    }
    # set work directory for SNP-EPI clusters
    fc_dir <- str_replace_all(fc_val," ","_")
    if(any(toupper(c("Transmission","Community")) %in% toupper(clust_type))){
      work_dir <- file.path(out_dir,fc_dir,
                            paste0("cluster-analysis","_SNPcutoff",snpco,"_Days",daysco))
     
    }else{
      work_dir <- file.path(out_dir,fc_dir,
                            paste0("cluster-analysis","_SNPcutoff",snpco))
    }
  
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
    kmx <- nrow(mat) - 2
    
    if(kmx < 2){
      next
    }
    
    if(nrow(mat)>15){
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
    if(any(colnames(fc_df_01)==Var_01)){
      datesJoin <- fc_df_01 %>% group_by(.data[[Var_01]])
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
    if(toupper(clust_type) == toupper("Transmission")){
      
      excl_vec <- c("Date2","epicumsum","CG","num","name","cluster","km_cluster",
                    "Hospital","Ward")
      
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
    }
   
    
    
    # Create scatter plots ----------------------------------------------------
    if(exists("clusterSet3")){
      if(!is.null(clusterSet3)){
        #stop("Provide path to the SNP distance matrix in the config file")
        
        
        scatter_plots_cmb_list <- list()
        scatter_plots_cmb_list <- create_scatter_plots(datesJoin,clusterSet3,mlst,transmission_lvl=trans_lvl)
        
        plotDF <- scatter_plots_cmb_list[[1]]
        p1 <- scatter_plots_cmb_list[[2]]
        # save scatter plots to file
        ggsave(file.path(work_dir,paste0(fc_val,".scatterplot.",date_var,".png")),
               p1,
               width = 10, 
               height = 8)
        
        
        
        incl_vec2 <- c("sampleID", Main_var ,"Hospital","Ward",Var_01,Var_02,"Epiweek","Days","SNPs", "Clusters","ST")
        
        metadf <- plotDF %>%
          # metadf <- px1 %>%
          ungroup() %>%
          dplyr::select(-epiweek) %>%
          # add_row(sampleID="reference",Var_00=NA,Ward=NA,Var_01=NA,TakenDate=as.Date(refDate),Epiweek=NA,ST=as.factor(refST)) %>%
          dplyr::select(any_of(incl_vec2)) %>% distinct(sampleID,.keep_all = TRUE) %>%
          column_to_rownames(var = "sampleID")
        
      }
    }
    
    
    # metadf <- plotDF %>%
    #   # metadf <- px1 %>%
    #   ungroup() %>%
    #   column_to_rownames(var = "sampleID")
    
    # Heatmap graph -----------------------------------------------------------
    maxP <-   c(
      RColorBrewer::brewer.pal(12,'Paired'),
      RColorBrewer::brewer.pal(12,'Set3')
    )
    
    
    maxX <-   c(
      RColorBrewer::brewer.pal(9,'Set1'),
      RColorBrewer::brewer.pal(8,'Dark2'),
      RColorBrewer::brewer.pal(2,"Accent")
    )
    
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
      # mutate(ST=as.numeric(as.character({{Main_var}}))) %>%
      inner_join(snpClust, by="name") %>%
      dplyr::rename("Core.SNP.clusters"=cluster)
    
    if(any(names(annotDF_subset3) == 'ST')){
      annotDF_subset2 <- annotDF_subset2 %>%
        mutate(ST=as.numeric(as.character(ST)))
      
      # Set colors for ST
      n3 <- annotDF_subset2 %>% pull(ST) %>% dplyr::n_distinct()
      col3 <- maxX[1:n3]#brewer.pal(n5,"Paired")
      names(col3) <- annotDF_subset2 %>% pull(ST) %>% unique()
    }
      
      
      
   
    
    # vec_incl_filt <- c("Hospital","WardType","ST","Core.SNP.clusters","SNP.Epi.Clusters")
    # vec_excl_filt <- c("km_cluster","Days","SNPs", "Cluster_Cases_count")
    # vec_excl_filt <- c("TakenDate","epiyear","epiweek", "epiyearweek","km_cluster",
    # "Serotype", "GPSC","Days","SNPs", "Cluster_Cases_count")
    
    vec_incl_filt_hm <- c(Main_var,"Hospital",Var_01,Var_02,"ST",
                          "name","Core.SNP.clusters","SNP.Epi.Clusters")
    
    # if(nrow(clusterSet3) != 0){
    if(exists("clusterSet3")){
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
    }
  
      
    #   annotDF_subset2 <- annotDF_subset2 %>%
    #     column_to_rownames(var = "name") %>%
    #     dplyr::select(! any_of(vec_excl_filt)) 
    #   
    # }else{
    #   next
    # }
      
    m_variable <- str_to_sentence(Main_var)
    annotDF_subset2 <- annotDF_subset2 %>%
      column_to_rownames(var = "name") %>%
      dplyr::rename(!!m_variable := Main_var)
    
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
    n1 <- annotDF_subset2 %>% pull(m_variable) %>% dplyr::n_distinct()
    col1 <- brewer.pal(n1,"Set3")
    names(col1) <- annotDF_subset2 %>% pull(m_variable) %>% unique()
    
    
    if(any(names(annotDF_subset2) == Var_01)){
      n2 <- annotDF_subset2 %>% pull(Var_01) %>% dplyr::n_distinct()
      col2 <- maxX[1:n2]  #brewer.pal(n2,"YlGnBu")
      names(col2) <- annotDF_subset2 %>% pull(Var_01) %>% unique()
    }
    
    
    # if(any(names(annotDF_subset2) == Var_02)){
    #   mutate(epiyear = lubridate::epiyear(.data[[Var_04]]))
    #   annotDF_subset2 %>%
    #     dplyr::mutate(epiyear=lubridate::year(.data[[Var_02]]))
    #   
    #   n6 <- annotDF_subset2 %>% pull(Var_02) %>% dplyr::n_distinct()
    #   col2 <- maxX[1:n2]  #brewer.pal(n2,"YlGnBu")
    #   names(col2) <- annotDF_subset2 %>% pull(Var_01) %>% unique()
    # }
    
    
    
    # n3 <- annotDF_subset2 %>% pull(Main_var) %>% dplyr::n_distinct()
    # col3 <- maxP[1:n3] #brewer.pal(n3,"Accent")
    # names(col3) <- annotDF_subset2 %>% pull(Main_var) %>% unique()
    
    
    n4 <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% dplyr::n_distinct()
    col4 <- maxP[1:n4] #brewer.pal(n4,"Accent")
    names(col4) <- annotDF_subset2 %>% pull(Core.SNP.clusters) %>% unique()
    
    
    # n5 <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% dplyr::n_distinct()
    # col5 <- brewer.pal(n5,"Paired")
    # names(col5) <- annotDF_subset2 %>% pull(SNP.Epi.Clusters) %>% unique()
    
    
    
    cols_hm <- c(m_variable,Var_01,"ST","Core.SNP.clusters","SNP.Epi.Clusters")
    annotDF_subset2 <- annotDF_subset2 %>%
      dplyr::select(any_of(cols_hm))
    
    if (any(names(annotDF_subset2) == "SNP.Epi.Clusters")){
      if(any(names(annotDF_subset2) == Var_01)){
        annotDF_subset2 <- annotDF_subset2 #%>% dplyr::rename("WardType" = Var_01)
        colAnnot2 <- ComplexHeatmap::HeatmapAnnotation(
          df = annotDF_subset2, annotation_height = 7,
          annotation_name_gp = gpar(fontsize = 7),
          col = list(m_variable = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
                     Var_01 = col2[! is.na(names(col2))],
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
          col = list(m_variable = col1[! is.na(names(col1))], #c("DORA NGINZA HOSPITAL" = "#8DD3C7"),
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
        col = list(m_variable = col1[! is.na(names(col1))], 
                   # WardType = col2[! is.na(names(col2))],
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
    if(exists("clusterSet3")){
      if(!is.null(clusterSet3)){
        snpepiDF <- clusterSet3 %>% dplyr::select(4,3,2,1)
        openxlsx::write.xlsx(snpepiDF, 
                             file.path(work_dir,
                                       paste(fc_val,"Clusters-SNPs-Epi-Definition",date_var,"xlsx",sep = ".")),
                             overwrite = T)
      }
    }
  
    
    
    
    # write MLST profile to file
    st_remove <- c("ST")
    incl_vec2 <- c("sampleID",Main_var,"Hospital","Ward","WardType",Var_02,"Epiweek","ST","Clusters")
    
    if(clust_type == "Transmission"){
      if(!is.null(clusterSet3)){
          outDF <- plotDF %>% 
          dplyr::select(!any_of(st_remove)) %>%
          inner_join(mlst,by=c("sampleID"="FILE")) %>% dplyr::select(any_of(incl_vec2))
      # names(mlst) <- c("sampleid","ST")
      }
     
    }else{
      outDF <- fc_df_01 %>%
        dplyr::select(!any_of(st_remove)) %>%
        inner_join(mlst,by=c("sampleID"="FILE")) %>% 
        left_join(snpClust,by=c("sampleID"="name"))
        # dplyr::select(any_of(incl_vec2))
      
    }
   
    if(exists("outDF")){
       openxlsx::write.xlsx(outDF, 
                            file.path(work_dir,
                                      paste(fc_val,"MLST_profile",date_var,"xlsx",sep = ".")),
                            overwrite = T)
     }
   
    
    
    
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
    if(clust_type == "Core"){
      next
    }
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
    
    
    if(clust_type == "Core"){
      next
    }
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
      if(exists("clsters")){
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

