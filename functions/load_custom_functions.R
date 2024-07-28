# Functions ---------------------------------------------------------------

# https://gist.github.com/jgilfillan/23336d0f5bcfffe6a71d0bdd634d023e
# group rows based on cumsum with reset
# modified 03/05/2023


# Function 00 -------------------------------------------------------------

cumsum_group <- function(x, threshold) {
  #x <- na.omit(x)
  x[is.na(x)] <- -555
  grp <- 0
  cumsum <- 0
  group <- 1
  result <- numeric()

  for (i in 1:length(x)) {
    # cumsum <- cumsum + abs(x[i]) 
    # 2024-07-21: Breaking clusters based on cumsum can potentially 
    # break cluster chains. Now applying SNP cut-off pair-wisely.
    # This is based on feedback received from the Klebs CoP session - 20240716
    
    cumsum <- abs(x[i])

    if (abs(cumsum) > threshold) {
      if(sum(abs(x[i]),abs(x[i-1])) <= threshold){
        group <- group + 1
        cumsum <- x[i]
      } else if(abs(x[i]) <= threshold){
        group <- group + 1
        cumsum <- x[i]

      }else{
        grp <- group
        group <- 0 #group + 1
        cumsum <- 0 #x[i]
      }

    }

    result = c(result, group)
    if(group == 0 & i != 1){
      group <- grp + 1
    }else if(group == 0 & i == 1){
      group = 1
    }


  }

  return (result)
}


# Function 01 -------------------------------------------------------------

# Default SNP cut-off threshold == 20

connect_samples <- function(x,threshold=20){
  rn_idx <- x$rn
  
  vec_id1 <- x %>% pull(sampleID)
  vec_id2 <- x %>% pull(ID2)
  vec_x3 <- x %>% pull(X3)
  vec_rn <- x %>% pull(rn)
  
  
  vec_id3 <- vector()
  pos_idx = 0
  list_checks <- list()
  grpCNT = 0
  # df_check <- data.frame()
  vec_excl <- vector()
  # vec_checks <- vector()
  added_combn_list <- list()
  
  for(i in rn_idx){
    
    if(pos_idx !=0){
      i = pos_idx
    }
    
    if(as.integer(vec_x3[[i]]) <= threshold){
      grpCNT <- grpCNT + 1
      id1 <- vec_id1[[i]]
      val2 <- vec_id2[[i]]
      snp_val <- vec_x3[[i]]
      
      vec_in_pair <- c(id1,val2)
      
      if(i!=1){
        for(j in seq_along(list_checks)){
          # if(j>=1){
          #   j <- j+1
          # }
          
          id_pair_vec <- c(list_checks[[j]][1],list_checks[[j]][2])
          if(all(vec_in_pair %in% id_pair_vec)){
            # do nothing
            added_combn_list[[j]] <- vec_in_pair
            id1<- NULL
            val2<-NULL
            snp_val <- NULL
            grpCNT<-NULL
            id_pair_vec <-NULL
            vec_in_pair <- NULL
            next
          }#else{
          #}
        }
        
        # if(is.null(vec_in_pair)){
        #   next
        # }
        vec_col1_vals <- unlist(map(list_checks, 1)) %>% unique()
        vec_col2_vals <- unlist(map(list_checks, 2)) %>% unique()
        # if(any(vec_col1_vals == id1 & vec_col2_vals == val2)){
        if(any(vec_col1_vals == id1)){ #& vec_col2_vals == val2)){
          pos1 <- which(vec_in_pair == id1)
          id1_incl <- vec_in_pair[[pos1]]
          idA <- id1
          idB <- val2
          # pos_idx = 0
        }else if(any(vec_col2_vals == val2)){
          pos1 <- which(vec_in_pair == val2)
          id1_incl <- vec_in_pair[[pos1]]
          idA <- val2
          idB <- id1
          # pos_idx = 0
        }else{
          # pos_idx = 0
          idA <- id1
          idB <- val2
        }
        
        vec_id3=as.vector(c(idA,idB,snp_val,grpCNT))
        pos_idx = 0
        
      }else{
        vec_id3=as.vector(c(vec_in_pair,snp_val,grpCNT))
        pos_idx = 0
      }
      
      
      
    }else{
      grpCNT= grpCNT + 2
    }
    
    
    list_checks[[i]] <- vec_id3
    id1<- NULL
    val2<-NULL
    snp_val <- NULL
    vec_id3 <- NULL
    grpCNT <- 0
  }
  
  
  return(list_checks)
  
}




# Function 02 -------------------------------------------------------------

# Only use actual value of sigN for K -- now adding 2 directly to sigNN in main script -- 2024-07-13

# get_core_snp_clusters <- function(m=mat, k=sigN, max=sigNN, dates=NULL, snpco=snpco, orig=TRUE){
get_core_snp_clusters <- function(m, k, max, dates=NULL, snpco, orig=TRUE){
  library(factoextra)
  
  if(nrow(m)<= max){
    max=nrow(m)-1
  }
  
  # calculate kmeans and group closely related samples
  res.km <- eclust(m, 
                   FUNcluster = "kmeans",
                   k=k,
                   k.max = max,
                   nboot = 500,
                   nstart = 25) #nstart = 25,
  
  
  # Get optimal clusters only
  vec_widths <- res.km$silinfo$clus.avg.widths
  if(is.null(vec_widths)){
    return(NA)
  }
  vec_optimal <- abs(res.km$silinfo$clus.avg.widths) >= 0.5
  vec_cluster <- res.km$silinfo$widths %>% pull(cluster) %>% unique() %>% as.character()
  names(vec_widths) <- vec_cluster
  
  vec_keep <- names(vec_widths[vec_optimal])
  vec_excl <- names(vec_widths[!vec_optimal])
  vec_km_cluster <- res.km$cluster
  
  
  # df_dates <- dates %>% 
  #   ungroup() %>%
  #   dplyr::select(sampleID,TakenDate)
  
  
  snp_data_list <- list()
  # get core SNP clusters -- with total SNPS <= to cut-off
  maxV=0
  list_clusters <-list()
  # keep_df3 <- data.frame()
  for(i in seq_along(vec_keep)){
    if(is.null(i)){next}
    
    grp <- vec_keep[[i]]
    names_vec <- res.km$silinfo$widths %>% 
      dplyr::filter(cluster %in% grp) %>%
      rownames_to_column(var="ID") %>%
      pull(ID)
    
    if(length(names_vec)<=2){next}
    
    clster_ids_df <- res.km$silinfo$widths %>% 
      dplyr::filter(cluster %in% grp) %>%
      rownames_to_column(var="ID") %>%
      dplyr::select(1,2) %>%
      dplyr::rename("km_cluster"=cluster)
    
    
    if(orig){
      # ## This is the original logic - which randomly compares isolates in same km cluster 
      # and creates SNP clusters from least pairwise SNP differences up to cut-off
      # 2023-09-27
      
      comps <- combn(seq_along(names_vec),
                     m=2, simplify = F,
                     FUN = function(x)names_vec[x])
      
      compare_df <- as.data.frame(do.call(rbind, comps))
      names(compare_df) <- c("X1","X2")
      
      # filter 1: apply snp cut-off
      check_df1 <- compare_df %>%
        inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
        dplyr::filter(X1 != X2 & X3 <=snpco )
      
      # filter 2: group isolates within K-clusters 
      # into SNP clusters based on provided SNP threshold
      
      if(nrow(check_df1) > 0){
        keep_df1 <- check_df1 %>%
          arrange(X3) %>%
          mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
          dplyr::filter(snpcumsum != 0) %>%
          mutate(num = as.numeric(paste0(X3,snpcumsum)))
      }else{
        next
      }
      
      
      if(nrow(keep_df1)>0){
        keep_df2 <- keep_df1 %>%
          ungroup() %>%
          dplyr::rename("clst"=snpcumsum) %>%
          dplyr::rename("SNPs"=X3) 
        
        
        
        df_filtx <- keep_df2 %>%
          arrange(SNPs) %>%
          mutate(rn=row_number())
        
        names(df_filtx) <- c("sampleID","ID2","X3","clust","num","rn")
        
        # filter 3: (main filter)
        list_checks <- list()
        list_checks <-connect_samples(df_filtx,threshold=snpco)
        
        
        # filter 4
        len01 <- length(list_checks)
        if(len01 == 1){
          len02 <- len01
        }else{
          len02 <- len01 - 1
        }
        
        if(is.null(list_checks[[len02]])){
          df_001 <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
          df_000 <- tail(df_001,1) %>% dplyr::mutate(V4="1")
          df_0000 <- df_001[1:nrow(df_001)-1,]
          df_002 <- bind_rows(df_000,df_0000)
        }else{
          df_002 <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
        }
        
        # df_filty <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
        df_filty <- df_002 
        
        df_filty <- df_filty %>%
          # group_by(V4) %>%
          mutate(snpcumsum8 = as.character(cumsum_group(as.numeric(V3), snpco))) %>%
          dplyr::mutate(clstGRP=str_c(V4,snpcumsum8)) %>%
          mutate(snpcumsum3 = data.table::rleid(clstGRP))     #20240529
        
        
        grps_vec <- df_filty %>% pull(snpcumsum3) %>% unique()
        
        list_dfs <- list()
        for(z in seq_along(grps_vec)){
          
          if(z==1){
            k=z
          }else{
            k=c(1:z)
          }
          
          
          if(z==1){
            grp_df <- df_filty %>% 
              dplyr::filter(snpcumsum3 %in% k)
            list_dfs[[k]] <- grp_df
            nm_vec <- unique(c(grp_df$V1,grp_df$V2))
          }
          
          
          j=z+1
          
          keep <- df_filty %>% 
            dplyr::filter(snpcumsum3 == j) %>%
            dplyr::filter(! V1 %in% nm_vec) %>%
            dplyr::filter(! V2 %in% nm_vec)
          
          if(nrow(keep) != 0){
            list_dfs[[j]] <- keep
          }
          
          
          bnd_df <- bind_rows(list_dfs)
          nm_vec <- unique(c(bnd_df$V1,bnd_df$V2))
          
          if(j==length(grps_vec)){
            break
          }
        }
        
        
        
        df_filtz <- bind_rows(list_dfs) %>%
          mutate(Clusters = data.table::rleid(snpcumsum3)+maxV) %>%
          dplyr::rename("sampleID"=V1,
                        "ID2"=V2,
                        "SNPs"=V3,
                        "clst"=snpcumsum3) %>%
          mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
          mutate(across(all_of(c("SNPs","clst")), as.numeric))
        
        
        x_vec1 <- df_filtx %>% pull(sampleID) %>% unique()
        x_vec2 <- df_filtx %>% pull(ID2) %>% unique()
        x_vec3 <- c(x_vec1,x_vec2) %>% unique()
        
        z_vec1 <- df_filtz %>% pull(sampleID) %>% unique()
        z_vec2 <- df_filtz %>% pull(ID2) %>% unique()
        z_vec3 <- c(z_vec1,z_vec2) %>% unique()
        
        # m_vec <- c(a_vec,b_vec,c_vec,d_vec) %>% unique()
        
        m_vec <- x_vec3[!x_vec3 %in% z_vec3]
        
        len<-abs(length(m_vec))
        if(len%%2==0 & len !=0){
          list_comp <- combn(seq_along(m_vec),
                             m=2, simplify = F,
                             FUN = function(x)m_vec[x])
          
          df_comps <- as.data.frame(do.call(rbind, list_comp))
          names(df_comps) <- c("X1","X2")
          
          maxV=max(df_filtz$Clusters)   
          
          df_snps <- df_comps %>%
            inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
            # dplyr::filter(X1 != X2 & X3 <snpco ) %>%
            arrange(X3) %>%
            mutate(snpcumsum4 = cumsum_group(as.numeric(X3), snpco)) %>%
            dplyr::filter(snpcumsum4 != 0) %>%
            mutate(Clusters = data.table::rleid(snpcumsum4)+maxV) %>%
            dplyr::rename("sampleID"=X1,
                          "ID2"=X2,
                          "SNPs"=X3,
                          "clst"=snpcumsum4) %>%
            mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
            mutate(across(all_of(c("SNPs","clst")), as.numeric))
          
          if(nrow(df_snps)==0){
            df_snps<-NULL
          }
          
        }else{
          df_snps<-NULL
        }
        
        
        maxV=0
        
        keep_df3 <- dplyr::bind_rows(list(df_filtz,df_snps)) %>%
          pivot_longer(cols = c(sampleID,ID2),values_to = "sampleID") %>%
          dplyr::group_by(Clusters) %>%
          dplyr::distinct(sampleID, .keep_all = TRUE) %>% 
          ungroup()
        
        
        if(any(names(keep_df3) == "Clusters")){
          # if(! is.null(keep_df3)){
          keep_df <- keep_df3
        }else{
          keep_df <- keep_df3 %>%
            mutate(Clusters = data.table::rleid(clst)+maxV) %>%
            group_by(Clusters) %>%
            mutate(cnt=n()) %>%
            dplyr::filter(cnt > 1) %>%
            dplyr::select(SNPs,clst,Clusters,name,sampleID) 
        }
        
        
      }else{
        next
      }
      
      
      
    
  
    }else{
      
      if(is.null(dates)){
        next
      }
      
      df_dates <- dates %>% 
        ungroup() %>%
        dplyr::select(sampleID,TakenDate)
      
      
      check_df1 <- clster_ids_df %>% 
        dplyr::inner_join(df_dates,by=c("ID"="sampleID")) %>% 
        dplyr::arrange(TakenDate) %>% 
        mutate(ID2 = dplyr::lead(ID,1)) %>% 
        # dplyr::select(-c(2,3)) %>%
        dplyr::select(ID,ID2,everything()) %>%
        inner_join(mdf,by=c("ID"="X1","ID2"="X2")) %>%
        dplyr::filter(ID != ID2 & X3 <=snpco ) %>%
        dplyr::rename("X1"=ID,"X2"=ID2)
      
      if(nrow(check_df1) > 0){
        keep_df1 <- check_df1 %>%
          dplyr::arrange(TakenDate) %>%
          mutate(snpcumsum = cumsum_group(X3, snpco)) %>%
          dplyr::filter(snpcumsum != 0) %>%
          mutate(num = as.numeric(paste0(X3,snpcumsum))) %>%
          dplyr::select(-TakenDate)
      }else{
        next
      }
    }
    
    
    maxV=max(keep_df$Clusters)
        list_clusters[[i]] <- keep_df %>%
          dplyr::select(sampleID,Clusters) %>%
          dplyr::rename("name"=sampleID,
                        "cluster"=Clusters) %>%
          inner_join(clster_ids_df,by=c("name"="ID"))

        # keep_df3=NULL
        keep_df2=NULL
        keep_df=NULL
    
    }
    if(length(list_clusters) != 0){
      snpClust <- bind_rows(list_clusters) %>%
        mutate(km_cluster=as.character(km_cluster),
               cluster=as.character(cluster))
    }else{
      return(NA)
    }
    # return(snpClust) # output is a data frame
    if(nrow(snpClust)!=0){
      snp_data_list[[1]] <- snpClust
      snp_data_list[[2]] <- res.km
    }
    if(length(snp_data_list) != 0){
      return(snp_data_list)
    }else{
      return(NA)
    }
    
    
  }
      ##########################################################
      
      
    
    
    
#     
#     if(nrow(keep_df1)>0){
#       
#       keep_df2 <- keep_df1 %>%
#         ungroup() %>%
#         pivot_longer(cols = c(X1,X2),values_to = "sampleID") %>%
#         # inner_join(df_dates,by="sampleID") %>%
#         dplyr::rename("SNPs"=X3) %>%
#         distinct(sampleID, .keep_all = TRUE) %>%
#         # mutate(snpcumsum2 = cumsum_group(SNPs, snpco)) %>%
#         # dplyr::filter(snpcumsum2 != 0) %>%
#         # dplyr::rename("clst"=snpcumsum2) 
#         dplyr::rename("clst"=snpcumsum)
#       
#       if(length(keep_df2$SNPs[keep_df2$SNPs == 0])>=3){ 
#         
#         if(orig){
#           df_filtx <- keep_df2 %>%
#             # arrange(SNPs,TakenDate) %>% # 2023-10-06
#             arrange(SNPs) %>%
#             mutate(ID2=lead(sampleID)) %>%
#             mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
#             left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>%
#             mutate(rn=row_number())
#         }else{
#           # Use option below when using TakenDates - 2023-09-27
#           df_filtx <- keep_df2 %>% 
#             arrange(TakenDate) %>%
#             mutate(ID2=lead(sampleID)) %>%
#             mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
#             left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>% 
#             mutate(rn=row_number())
#         }
#         
#         # filter 1
#         list_checks <- list()
#         list_checks <-connect_samples(df_filtx)
#         
#         # df_filtx %>% 
#         # dplyr::mutate(snpt=cumsum_group(as.numeric(X3),snpco)) %>% 
#         # print(n=40)
#         
#         # filter 2
#         len01 <- length(list_checks)
#         len02 <- len01 - 1
#         if(is.null(list_checks[[len02]])){
#           df_001 <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
#           df_000 <- tail(df_001,1) %>% dplyr::mutate(V4="1")
#           df_0000 <- df_001[1:nrow(df_001)-1,]
#           df_002 <- bind_rows(df_000,df_0000)
#         }else{
#           df_002 <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
#         }
#         
#         # df_filty <- as.data.frame(do.call(rbind, list_checks)) %>% distinct()
#         df_filty <- df_002
#         
#         
#         df_filty <- df_filty %>%
#           # group_by(V4) %>%
#           mutate(snpcumsum8 = as.character(cumsum_group(as.numeric(V3), snpco))) %>%
#           dplyr::mutate(clstGRP=str_c(V4,snpcumsum8)) %>%
#           mutate(snpcumsum3 = data.table::rleid(clstGRP))     #20240529
#         
#         
#         
#         grps_vec <- df_filty %>% pull(snpcumsum3) %>% unique()
#         
#         list_dfs <- list()
#         for(z in seq_along(grps_vec)){
#           
#           if(z==1){
#             k=z
#           }else{
#             k=c(1:z)
#           }
#           
#           
#           if(z==1){
#             grp_df <- df_filty %>% 
#               dplyr::filter(snpcumsum3 %in% k)
#             list_dfs[[k]] <- grp_df
#             nm_vec <- unique(c(grp_df$V1,grp_df$V2))
#           }
#           
#           
#           j=z+1
#           
#           keep <- df_filty %>% 
#             dplyr::filter(snpcumsum3 == j) %>%
#             dplyr::filter(! V1 %in% nm_vec) %>%
#             dplyr::filter(! V2 %in% nm_vec)
#           
#           if(nrow(keep) != 0){
#             list_dfs[[j]] <- keep
#           }
#           
#           
#           bnd_df <- bind_rows(list_dfs)
#           nm_vec <- unique(c(bnd_df$V1,bnd_df$V2))
#           
#           if(j==length(grps_vec)){
#             break
#           }
#         }
#         
#         
#         df_filtz <- bind_rows(list_dfs) %>%
#           mutate(Clusters = data.table::rleid(snpcumsum3)+maxV) %>%
#           dplyr::rename("sampleID"=V1,
#                         "ID2"=V2,
#                         "SNPs"=V3,
#                         "clst"=snpcumsum3) %>%
#           mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
#           mutate(across(all_of(c("SNPs","clst")), as.numeric))
#         
#         
#         # get SNP pairs excluded in by above filter 1
#         # a_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V1) %>% pull(sampleID)
#         # b_vec <- df_filtx %>% filter(! sampleID %in% df_filty$V2) %>% pull(sampleID)
#         # c_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V1) %>% pull(ID2)
#         # d_vec <- df_filtx %>% filter(! ID2 %in% df_filty$V2) %>% pull(ID2)
#         
#         
#         x_vec1 <- df_filtx %>% pull(sampleID) %>% unique()
#         x_vec2 <- df_filtx %>% pull(ID2) %>% unique()
#         x_vec3 <- c(x_vec1,x_vec2) %>% unique()
#         
#         z_vec1 <- df_filtz %>% pull(sampleID) %>% unique()
#         z_vec2 <- df_filtz %>% pull(ID2) %>% unique()
#         z_vec3 <- c(z_vec1,z_vec2) %>% unique()
#         
#         # m_vec <- c(a_vec,b_vec,c_vec,d_vec) %>% unique()
#         
#         m_vec <- x_vec3[!x_vec3 %in% z_vec3]
#         
#         len<-abs(length(m_vec))
#         if(len%%2==0 & len !=0){
#           list_comp <- combn(seq_along(m_vec),
#                              m=2, simplify = F,
#                              FUN = function(x)m_vec[x])
#           
#           df_comps <- as.data.frame(do.call(rbind, list_comp))
#           names(df_comps) <- c("X1","X2")
#           
#           maxV=max(df_filtz$Clusters)   
#           
#           df_snps <- df_comps %>%
#             inner_join(mdf,by=c("X1"="X1","X2"="X2")) %>%
#             # dplyr::filter(X1 != X2 & X3 <snpco ) %>%
#             arrange(X3) %>%
#             mutate(snpcumsum4 = cumsum_group(as.numeric(X3), snpco)) %>%
#             dplyr::filter(snpcumsum4 != 0) %>%
#             mutate(Clusters = data.table::rleid(snpcumsum4)+maxV) %>%
#             dplyr::rename("sampleID"=X1,
#                           "ID2"=X2,
#                           "SNPs"=X3,
#                           "clst"=snpcumsum4) %>%
#             mutate(across(all_of(c("sampleID","ID2")), as.character)) %>%
#             mutate(across(all_of(c("SNPs","clst")), as.numeric))
#             
#             if(nrow(df_snps)==0){
#               df_snps<-NULL
#             }
#           
#         }else{
#           df_snps<-NULL
#         }
#         
#         
#         maxV=0
#         
#         keep_df2 <- dplyr::bind_rows(list(df_filtz,df_snps)) %>%
#           pivot_longer(cols = c(sampleID,ID2),values_to = "sampleID") %>%
#           dplyr::group_by(Clusters) %>%
#           dplyr::distinct(sampleID, .keep_all = TRUE) %>% 
#           ungroup()
#         
#         # maxV=max(keep_df3$Clusters) 
#         
#       }#else{ 
#        # next
#       #}
#       
#       if(any(names(keep_df2) == "Clusters")){
#         # if(! is.null(keep_df3)){
#         keep_df <- keep_df2
#       }else{
#         keep_df <- keep_df2 %>%
#           mutate(Clusters = data.table::rleid(clst)+maxV) %>%
#           group_by(Clusters) %>%
#           mutate(cnt=n()) %>%
#           dplyr::filter(cnt > 1) %>%
#           dplyr::select(SNPs,clst,Clusters,name,sampleID) 
#       }
#       
#       
#     }else{
#       next
#     }
#     
#     maxV=max(keep_df$Clusters)
#     list_clusters[[i]] <- keep_df %>% 
#       dplyr::select(sampleID,Clusters) %>%
#       dplyr::rename("name"=sampleID,
#                     "cluster"=Clusters) %>%
#       inner_join(clster_ids_df,by=c("name"="ID"))
#     
#     # keep_df3=NULL
#     keep_df2=NULL
#     keep_df=NULL
#     
#   }
#   # return(list_clusters)
#   if(length(list_clusters) != 0){
#     snpClust <- bind_rows(list_clusters) %>%
#       mutate(km_cluster=as.character(km_cluster),
#              cluster=as.character(cluster))
#   }else{
#     return(NA)
#   }
#   # return(snpClust) # output is a data frame
#   if(nrow(snpClust)!=0){
#     snp_data_list[[1]] <- snpClust
#     snp_data_list[[2]] <- res.km
#   }
#   if(length(snp_data_list) != 0){
#     return(snp_data_list)
#   }else{
#     return(NA)
#   }
#   
#   
# }
# 


# Function 04 -------------------------------------------------------------


drw_network <- function(nodes,edges,b=1,st,type="Core"){
  ntwrk_list <- list()
  
  if(any(names(nodes) == "vari")){
    names(nodes)[which(names(nodes) == "vari")] <- "group"
  }
  
  
  # if(all(is.na(nodes$group))){
  #   next
  # }
  nodes$group <- as.character(nodes$group)
  nodes$group[is.na(nodes$group)] <- "bg"
  
  st=as.name(st)
  clr = c("#ABC2E8","#FFF338","#82CD47","#FF78C4","#525FE1","#98EECC","#F29727","#D3D04F","#10A19D","#B04759","#D09CFA","#B9E9FC",
          "#FFE7A0","#FF6969","#00FFCA","#ECF2FF","#E86A33","#569DAA","#FFE5CA","#FA9884","#A6BB8D","#C8B6A6","#8B1874","#FFA0A0")
  # "#C68B59",,"#C3ACD0"
  
  y=length(clr)
  groupname=unique(nodes$group)[order(unique(nodes$group))]
  
  if(any(groupname == "bg", na.rm = T)){
    groupname<-groupname[groupname != "bg"]
    len <- length(groupname) # if reference is included
  }else{
    len <- length(groupname)
  }
  
  c=(len+b)-1
  if(len+b > y){
    b=1
    c=len
  }
  
  assign_colors <- data.frame(groupname=groupname, 
                              color=clr[b:c]
  )
  
  # add background color
  library(glue)
  bg_row <- c(groupname="bg",color="#F8F4EA")
  assign_colors <- rbind(assign_colors,bg_row)
  
  gl_code <- list()
  
  for (i in 1:nrow(assign_colors)){
    # print(i)
    gn <- as.character(assign_colors[i,1])
    cl <- as.character(assign_colors[i,2])
    gl_code[[i]] <- glue('visGroups(groupname = "{gn}",color = "{cl}", shadow = list(enabled = TRUE))')
  }
  
  
  code_app <- glue(paste(gl_code,collapse = ' %>% '))
  
  
  # visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%
  # visOptions(highlightNearest = TRUE, nodesIdSelection = F, selectedBy = "group") %>%
  
  vz_groups <- assign_colors %>% pull(1)
  
  if(type=="Core"){
    ntwk_code <- glue('visNetwork(nodes, edges, height = "500px", width = "100%",
    main = paste0("ST",{st},", ","Core-SNP Clusters:"," SNPcutoff: ",{snpco}),footer = "*bg=background, i.e. cases not included in clusters") %>%
    visEdges(arrows = "to",font = list(align="top",size=14, color=ifelse(edges$weight>11,"grey","green"))) %>%
    
    {code_app} %>%
    visLegend(main = "Core.SNP.Clusters",width = 0.1, position = "right") %>%
    visClusteringByGroup(groups = vz_groups, label = "cluster : ") %>%
    visHierarchicalLayout() %>%
    visLayout(randomSeed = 12)
    ')
  }else{
    ntwk_code <- glue('visNetwork(nodes, edges, height = "500px", width = "100%",
    main = paste0("ST",{st},", ","Core-SNP + EpiData Clusters:"," SNPcutoff: ",{snpco},", Days: ",{daysco}),footer = "*bg=background, i.e. cases not included in clusters") %>%
    visEdges(arrows = "to",font = list(align="top",size=14, color=ifelse(edges$weight>11,"grey","green"))) %>%
    
    {code_app} %>%
    visLegend(main = "SNP+Epi Clusters",width = 0.1, position = "right") %>%
    visClusteringByGroup(groups = vz_groups, label = "cluster : ") %>%
    visHierarchicalLayout() %>%
    visLayout(randomSeed = 12)
    ')
  }
  
  
  # write data to file
  nodes$X1 <- dplyr::coalesce(nodes$X1,"UNKNOWN")
  nodes$X2 <- dplyr::coalesce(nodes$X2,"UNKNOWN")
  nodes$X3 <- dplyr::coalesce(nodes$X3,"UNKNOWN")
  
  network <- eval(parse(text=ntwk_code))
  ntwrk_list[[1]] <- network
  ntwrk_list[[2]] <- nodes
  return(ntwrk_list)
}




# Function 05 -------------------------------------------------------------

run_core_snp_cluster_analysis <- function(){
  core_snp_results_list <- list()
  # get core SNP clusters
  res_list <- get_core_snp_clusters(m=mat, k=sigN, max=sigNN, dates=datesJoin, snpco=snpco, orig=T)
  
  snpClust <- res_list[[1]]
  
  if(!exists('snpClust')){
    return(NA)
  }
  
  # if(exists('snpClust') && is.na(snpClust)){return(NA)}
  if(exists('snpClust')){
    if(length(snpClust) == 1 && is.na(snpClust)){
      # print("Hello")
      return(NA)
    }
  }
  
  # if(is.na(snpClust)){
  #   return(NA)
  # }
  
  if(nrow(snpClust) == 0){
    return(NA)
  }#else if(is.na(snpClust)){next}
  
  
  # create matrix for heatmap
  snpClustID <- snpClust %>%
    mutate(km_cluster=as.character(km_cluster),
           cluster=as.character(cluster)) %>%
    dplyr::select(2,3) %>%
    distinct(.,.keep_all = T) 
  
  vec_keep_names <- snpClust %>% pull(name) %>% unique() 
  # mat_optimal_centres <- res.km$centers[,colnames(res.km$centers) %in% vec_keep_names]
  centers_mat <- res_list[[2]]$centers
  
  # if(is.na(centers_mat)){next}
  
  mat_optimal_centres <- centers_mat[,colnames(centers_mat) %in% vec_keep_names]
  # mat_optimal_centres <- mat_optimal_centres[! rownames(mat_optimal_centres) %in% vec_excl,]
  
  df2mat <- as.data.frame(mat_optimal_centres) %>%
    rownames_to_column(var = "rowID") %>%
    inner_join(snpClustID, by=c("rowID"="km_cluster")) %>%
    column_to_rownames(var = "cluster") %>%
    dplyr::select(-rowID)
  
  
  core_snp_results_list[[1]] <- snpClust
  core_snp_results_list[[2]] <- snpClustID
  core_snp_results_list[[3]] <- df2mat
  core_snp_results_list[[4]]<- vec_keep_names
  core_snp_results_list[[5]] <- mat_optimal_centres
  
  return(core_snp_results_list)
}


# Function 06 -------------------------------------------------------------

calculate_SNP_Epi_clusters <- function(snpClust,epiwkDF,mdf,excl_vec){
  
  vec_clust <- snpClust %>% pull(cluster) %>% unique()
  
  episnp_list = list()
  for(k in vec_clust){
    episnp_list[[k]] <- snpClust %>% dplyr::filter(cluster %in% k) %>%
      inner_join(epiwkDF, by=c("name"="sampleID")) %>%
      dplyr::arrange(TakenDate) %>%
      mutate(ID2 = dplyr::lead(name,1)) %>%
      mutate(Date2 = dplyr::lead(TakenDate,1)) %>%
      dplyr::select(name,ID2,TakenDate,Date2,everything()) %>%
      mutate(Days = as.numeric(difftime(Date2,TakenDate,units = "days"))) %>%
      mutate(epicumsum = if_else(Days <= daysco,1,0)) %>%
      mutate(epicumsum = if_else(is.na(epicumsum),0,epicumsum)) %>%
      ungroup() %>%
      inner_join(mdf,by=c("name"="X1","ID2"="X2")) %>%
      mutate(CG = data.table::rleid(epicumsum)) %>%
      dplyr::filter(epicumsum != 0) #%>%
    # print(n=40)
    
  }
  
  
  
  # excl_vec <- c("TakenDate","Date2","epicumsum","CG","num","name","cluster","km_cluster",
  #               "Hospital","Ward","WardType")
  
  clusterSet3 <- bind_rows(episnp_list) %>%
    mutate(num = as.numeric(paste0(cluster,CG))) %>%
    mutate(Clusters = data.table::rleid(num)) %>%
    ungroup() %>%
    pivot_longer(cols = c(name,ID2),values_to = "sampleID") %>%
    # select(! all_of(excl_vec)) %>%  #20240529
    select(! any_of(excl_vec)) %>%
    group_by(Clusters) %>%
    dplyr::rename("SNPs"=X3) %>%
    distinct(sampleID, .keep_all = TRUE) %>%
    dplyr::mutate(n_clusters=n()) %>%
    mutate(Clusters=as.factor(as.character(Clusters))) %>%
    dplyr::rename("Cluster_Cases_count"=n_clusters)
  
  
  
  clusterSet3 <- clusterSet3 %>%
    dplyr::select("sampleID","Days","SNPs","Clusters","Cluster_Cases_count")
  
  if(nrow(clusterSet3) == 0){ 
    return(NULL) 
  }else{
    return(clusterSet3)
  }
  
}



# Function 07 -------------------------------------------------------------

create_scatter_plots <- function(datesJoin,clusterSet3,mlst,transmission_type="facility"){
  
  scatter_plots_list <- list()
  
  plotDF <- datesJoin %>% 
    left_join(clusterSet3, by="sampleID") %>% 
    dplyr::rename("Epiweek"=epiyearweek) %>%
    ungroup() %>%
    group_by(Epiweek) %>%
    distinct(sampleID, .keep_all = TRUE) %>%
    dplyr::mutate(n_cases=n()) %>%
    ungroup() %>%
    dplyr::arrange(TakenDate) %>%
    mutate(denserank = data.table::rleid(TakenDate)) %>%
    dplyr::mutate(id=dplyr::row_number()) %>%
    # left_join(mlst,by=c("sampleID"="FILE")) %>%
    dplyr::rename("Cases_count"=n_cases)
  
  plotDF$Epiweek = factor(plotDF$Epiweek, levels=unique(plotDF$Epiweek)[plotDF$id], ordered = T)
  plotDF$Cases_count = dplyr::coalesce(plotDF$Cases_count,1)
  plotDF$Cases_count = as.factor(plotDF$Cases_count)
  plotDF$Clusters = as.factor(plotDF$Clusters)
  plotDF$ST = as.factor(plotDF$ST)
  
  
  
  # create scatter plots
  st_list <- list()
  st_final_list <- list()
  scatter_plots <- list()
  scatter_plots_cmb_list <- list()
  if(transmission_type == "community"){
    
    st_vec <- levels(plotDF$ST)
    st_vec <- st_vec[! str_detect(st_vec,"NOVEL")]
    for(i in seq_along(st_vec)){
      st_val <- st_vec[[i]]
      
      
      st_df <- plotDF %>%
        dplyr::filter(ST %in% st_val)
      
      st_to_exclude_vec <- st_df %>% 
        group_by(ST,Clusters) %>% 
        dplyr::summarise(grp_cnt=n()) %>%
        dplyr::filter(grp_cnt <2) %>%
        ungroup() %>%
        # mutate(Clusters=as.character(Clusters)) %>%
        pull(Clusters)
      
      st_df <- st_df %>% dplyr::filter(! Clusters %in% st_to_exclude_vec)
      
      if(nrow(st_df)>=2){
        st_list[[i]] <- st_df
        
        if(sum(st_df$Cluster_Cases_count, na.rm = T) > 0){
          cluster_n <-st_df %>% dplyr::filter(!is.na(Clusters)) %>% dplyr::summarise(n())
          if(cluster_n >= 2){
            st_final_list[[i]] <- st_df
            px1 <- st_df %>% 
              group_by(Facility_code,Epiweek) %>% 
              dplyr::mutate(count=n())
            
            p1 <- ggplot(px1, 
                         aes(x = Epiweek,
                             y = Facility_code , #WardType, #Clusters,
                             # y = count, #,
                             color=Clusters,
                             shape = ST #WardType
                         )) +
              geom_point(alpha = .9, size=4,position=position_jitter(h=0.07,w=0.15)) + 
              # scale_fill_brewer(palette = "Dark2") +
              labs(title = paste0("Clusters based on Epi definition: SNPs <= ",snpco," + ",daysco,"-day window"),
                   y = "Facility code",
                   x = "Epiweek") +
              theme_bw() +
              # scale_y_continuous(breaks = ~round(unique(pretty(.)))) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) +
              theme(axis.title = element_text(size = 15)) +
              theme(axis.text = element_text(size = 12)) +
              theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5, hjust = 1)) 
            
            
            ggsave(file.path(work_dir,paste0("Scatterplot.ST.",st_val,date_var,".png")),p1,width = 10, height = 8)
            
            scatter_plots[[i]] <- p1
          }
          
        }
        
      }else{
        next
      }
    }
  }else{
    if(all(c("Var_01","Epiweek") %in% colnames(plotDF))){
      px1 <- plotDF %>% group_by(Var_01,Epiweek) %>% dplyr::mutate(count=n())
      
      p1 <- ggplot(px1, 
                   aes(x = Epiweek,
                       y = ST, #WardType, #Clusters,
                       # y = count, #,
                       color=Clusters,
                       shape = Var_01
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
      
    }else{
      px1 <- plotDF %>% group_by(Epiweek) %>% dplyr::mutate(count=n())
      
      p1 <- ggplot(px1, 
                   aes(x = Epiweek,
                       y = ST, 
                       color=Clusters,
                       # shape = WardType
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
    }
    
    
   
    
    scatter_plots[[1]] <- p1
  }
  
  
  scatter_plots_list[[1]] <- plotDF
  scatter_plots_cmb_list <- c(scatter_plots_list, scatter_plots)
  # scatter_plots_list[[2]] <- 
  
  return(scatter_plots_cmb_list)
  
}
