library(RANN)
library(dplyr)
library(purrr)
library(stringr)
library(future.apply)

## tdens is the target phenotype we want to analyze.
## Examples are:
## CD8p, meaning cells showing ONLY CD8p marker
## CD8p_total, meaning cells that shows, AT LEAST, CD8p marker
## CD8p_CD68p_PD1neg_total, meaning all cells that shows CD8p AND CD68p, NOT PD1, plus other markers

line_selector_on_pheno=function(tdens, pheno_count){
  if(grepl("total", tdens)){
    all_marker=strsplit(tdens, "_")[[1]]
    ###we need to excluded eventual negative marker
    
    excl_marker=NULL
    if(grepl("neg", tdens)){
      excl_marker=all_marker[grepl("neg", all_marker)]
      excl_marker=sapply(excl_marker, function(x){
        return(str_replace (x, "neg", ""))
      })
      
      for(j in excl_marker){
        all_marker=all_marker[!grepl(j, all_marker)]
      }
    }
    
    all_marker=all_marker[!grepl("total", all_marker)]
    true_list=rep(TRUE, length(pheno_count))
    for(sing_marker in all_marker){
      if(length(sing_marker)!=0){
        true_list=true_list & grepl(sing_marker, pheno_count, fixed=TRUE)        
      }
    }
    
    for(sing_marker in excl_marker){
      true_list=true_list & !grepl(sing_marker, pheno_count, fixed=TRUE)
    }
  }
  
  else{
    true_list= tdens==pheno_count
  }
  return(true_list)
}



calc_triplets_mutual=function(samp_data, threshold=40, df_name_points){
  
  ## duplicate the matrix
  ## when using triples the combinations become more complicated
  df_name_points=rbind(df_name_points,
                       data.frame(Var1=df_name_points$Var2, Var2=df_name_points$Var1, Var3=df_name_points$Var3),
                       data.frame(Var1=df_name_points$Var3, Var2=df_name_points$Var1, Var3=df_name_points$Var2)
                       )
  
  
  ## we need to transform these columns in character
  df_name_points$Var1=as.character(df_name_points$Var1)
  df_name_points$Var2=as.character(df_name_points$Var2)
  df_name_points$Var3=as.character(df_name_points$Var3)
  
  
  d_data=samp_data
  d_data$ID=samp_data$sample
  d_data$type=samp_data$phenotype
  d_data$xpos=samp_data$nucleus.x
  d_data$ypos=samp_data$nucleus.y
  d_data=d_data[, c("ID", "type", "xpos", "ypos", "tissue.type") ]
  ## the data frame here MUST HAVE the following columns: "ID, type, xpos, ypos
  ## type must be as FACTORS
  ###d_data$type=as.factor(d_data$type)
  
  plan(multisession)
  cell_in_range=future_lapply(1:nrow(df_name_points), function(i, df_name_points, d_data){
    sp=df_name_points[i, 1]
    tp=df_name_points[i, 2]
    tp2=df_name_points[i, 3]
    
    ### checking the neighbors from sp to tp
    t_data=d_data[line_selector_on_pheno(sp, d_data$type) |
                    line_selector_on_pheno(tp, d_data$type), ]
    
    t_data$type[line_selector_on_pheno(sp, t_data$type)]=sp
    t_data$type[line_selector_on_pheno(tp, t_data$type)]=tp
    
    t_data$type=factor(t_data$type)
    
    tdf=count_in_cutoff(t_data, threshold, start = sp, target = tp,
                        a_max_neigh=1)
    
    if(nrow(tdf)==0){
      tdf=t_data[, c("ID", "type", "xpos", "ypos", "tissue.type")]
      tdf[, tp]=0
    }
    
    ### checking the neighbors from sp to tp2
    t_data=d_data[line_selector_on_pheno(sp, d_data$type) |
                    line_selector_on_pheno(tp2, d_data$type), ]
    
    t_data$type[line_selector_on_pheno(sp, t_data$type)]=sp
    t_data$type[line_selector_on_pheno(tp2, t_data$type)]=tp2
    
    t_data$type=factor(t_data$type)
    
    tdf2=count_in_cutoff(t_data, threshold, start = sp, target = tp2,
                         a_max_neigh=1)
    
    if(nrow(tdf2)==0){
      tdf2=t_data[1, c("ID", "type", "xpos", "ypos", "tissue.type")]
      tdf2[, tp2]=0
    }
    ### we need now to merge the two results
    tdf_m=merge(x=tdf, y=tdf2, by.x=c("ID", "type", "xpos", "ypos", "tissue.type"), by.y=c("ID", "type", "xpos", "ypos", "tissue.type"))
    ## remove the ".x" ".y"
    colnames(tdf_m)=str_remove_all(colnames(tdf_m), "\\.x|\\.y")
    
    tdf_m[, tp2]= as.numeric(tdf_m[, tp] & tdf_m[,tp2])
    tdf_m=tdf_m[, c("ID", "type", "xpos", "ypos", "tissue.type", tp2)]
    
    ## we need to generate a column name
    a_cname=paste0(tp, "_to_", tp2)
    
    colnames(tdf_m)=c("sample", "phenotype", "xpos", "ypos", "tissue.type", a_cname)
    ## we remove, from the results, some of the columns
    tdf_m=tdf_m[, c( "tissue.type", "phenotype", "sample", a_cname )]
    
    
    return(tdf_m)
    
  }, df_name_points=df_name_points, d_data=d_data)
  
  names(cell_in_range)=sapply( 1:nrow(df_name_points), function(i, df_name_points){
    sp=df_name_points[i, 1]
    tp=df_name_points[i, 2]
    tp2=df_name_points[i, 3]
    
    return(paste0(sp, "_to_", tp, "_to_", tp2))
  }, df_name_points=df_name_points)
  
  
  
  #From the cell to cell statistics we generate the  mutual statistics
  mutual_freq=list()
  
  for(i in 1:(nrow(df_name_points)/3) ){
    
    sp=df_name_points[i, 1]
    tp=df_name_points[i, 2]
    tp2=df_name_points[i, 3]
    
    
    tdf1=cell_in_range[[paste0(sp, "_to_", tp, "_to_", tp2)]]
    
    tdf2=cell_in_range[[paste0(tp, "_to_", sp, "_to_", tp2)]]
    
    tdf3=cell_in_range[[paste0(tp2, "_to_", sp, "_to_", tp)]]
    
    
    a_j_str="___"
    
    a_stat1=tdf1 %>% group_by(sample, tissue.type) %>% summarise_at(paste0(tp, "_to_", tp2), list(sum=sum, length=length))
    a_stat1$samp_tissue=paste0(a_stat1$sample, a_j_str, a_stat1$tissue.type)
    
    a_stat2=tdf2 %>% group_by(sample, tissue.type) %>% summarise_at(paste0(sp, "_to_", tp2), list(sum=sum, length=length))
    a_stat2$samp_tissue=paste0(a_stat2$sample, a_j_str, a_stat2$tissue.type)
    
    a_stat3=tdf3 %>% group_by(sample, tissue.type) %>% summarise_at(paste0(sp, "_to_", tp), list(sum=sum, length=length))
    a_stat3$samp_tissue=paste0(a_stat3$sample, a_j_str, a_stat3$tissue.type)
    
    
    a_stat1_nt=merge(a_stat1[, c("samp_tissue", "sum", "length")], a_stat2[, c("samp_tissue", "sum", "length")], by="samp_tissue", all=TRUE)
    a_stat1_nt=merge(a_stat1_nt, a_stat3[, c("samp_tissue", "sum", "length")], by="samp_tissue", all=TRUE)
    
    ## in the merge there are some missing values, we need to fix this for the column sample.x and tissue.type.x
    a_stat1_nt$sample=sapply(a_stat1_nt$samp_tissue, function(x){
      return(str_split(x, a_j_str)[[1]][1])
    })
    
    a_stat1_nt$tissue.type=sapply(a_stat1_nt$samp_tissue, function(x){
      return(str_split(x, a_j_str)[[1]][2])
    })
    
    
    a_stat1_nt[is.na(a_stat1_nt)]=0
    a_stat1_nt$freq=(a_stat1_nt$sum.x + a_stat1_nt$sum.y + a_stat1_nt$sum) /(a_stat1_nt$length.x + a_stat1_nt$length.y + a_stat1_nt$length)  
    a_stat1_nt=a_stat1_nt[, c("sample", "tissue.type", "freq")]
    mutual_freq[[paste0(sp, "_mutual_", tp, "_mutual_", tp2)]]=a_stat1_nt
    
    
  }
  
  return(mutual_freq)
}

