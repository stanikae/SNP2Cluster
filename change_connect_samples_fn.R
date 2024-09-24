
keep_df1 <- check_df1 %>%
  arrange(X3) %>%
  mutate(rn=row_number())



keep_df2 <- keep_df1 %>%
  ungroup() %>%
  pivot_longer(cols = c(X1,X2),values_to = "sampleID") %>%
  # inner_join(df_dates,by="sampleID") %>%
  dplyr::rename("SNPs"=X3) %>%
  distinct(sampleID, .keep_all = TRUE)


names(keep_df1) <- c("sampleID","ID2","X3","rn")


df_filtx <- keep_df2 %>%
  # arrange(SNPs,TakenDate) %>% # 2023-10-06
  arrange(SNPs) %>%
  mutate(ID2=lead(sampleID)) %>%
  mutate(ID2=replace_na(ID2,dplyr::first(sampleID))) %>%
  left_join(mdf,by=c("sampleID"="X1","ID2"="X2")) %>%
  mutate(rn=row_number())


connect_samples(keep_df1)



l_dfs <- list(data.frame(ID1=vec_id1,ID2=vec_id2),
data.frame(ID1=vec_id2,ID2=vec_id1))
bind_rows(l_dfs) %>% arrange(ID1,ID2) %>% dist

connect_samples <- function(x,threshold=snpco){
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

