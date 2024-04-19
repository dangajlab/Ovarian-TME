library(dplyr)
library(vroom)
library(here)
library(stringr)
library(tidyr)
source("triplets_code.r")


### creating a "temp" directory if it does not exists already
tempdir="temp"
if (!dir.exists(here(tempdir))) {dir.create(here(tempdir))}

### we set the temporary vroom directory to our generated temp directory
Sys.setenv(VROOM_TEMP_PATH = tempdir)

### we select the columns we want to keep, for the immunoclassification we will need only the Annotation.ID, tag, Tissue.Category, but
### for sake of other analyses we keep also others
pick_columns_df=c("Annotation.ID", "tag", "Tissue.Category", "Cell.X.Position", "Cell.Y.Position", "Phenotype.CD8",
                  "Phenotype.CD11c", "Phenotype.CD68", "Phenotype.CK", "Phenotype.PD1", "Phenotype.PDL1")


## We read all the data files of UPPEN.
datadir="data/" ###"../uppen_imcol_analysis/uppen_raw_distances/rawdata/"
all_df=lapply(list.files(datadir), function(i, pick_columns){
  cons_data=vroom(paste0(datadir, i), delim ="\t", show_col_types = FALSE)
  colnames(cons_data)=str_replace_all(colnames(cons_data), " ", ".")
  ## fill with "not-detected" the missing columns
  cons_data[, pick_columns[!pick_columns_df %in% colnames(cons_data)]]="not-detected"
  return(cons_data[, pick_columns_df])
}, pick_columns=pick_columns_df)

all_df=do.call(rbind, all_df)



### we remove the "Other tissue type
all_df=all_df[all_df$Tissue.Category!="Other", ]

### we count the number of CD8+ for each Annotatation.ID and for each Tissue.Category
counts_cd8=all_df %>%  group_by(tag, Annotation.ID, Tissue.Category) %>%  summarize(cd8_count = sum(grepl("CD8[+]", Phenotype.CD8))) %>% pivot_wider(values_from = cd8_count, names_from = Tissue.Category)
counts_cd8=counts_cd8 %>% rename(cd8p_stroma=Stroma, cd8p_tumor=Tumor)
counts_cd8[is.na(counts_cd8)]=0

### we compute the density --------->>>>>>>>>>>> PER ROI <<<<<<<<<<<<--------------
### note also that the size of the ROI matters, much smaller or much larger ROI will completely skew the immunoclassification

### a ROI is 0.64 mm2
roi_surface=0.64
counts_cd8 <- counts_cd8 %>%
  mutate(cd8p_stroma = cd8p_stroma / roi_surface,
         cd8p_tumor = cd8p_tumor / roi_surface)

### we compute the fraction of ROI that are inflammed in tumor and stroma using the threshold of 21 cells/mm2
perc_roi_infl= counts_cd8 %>% group_by(tag) %>% summarize( inf_tumor=sum(cd8p_tumor>21)/n(), inf_stroma=sum(cd8p_stroma>21)/n(), total=n()  )


### Assigning the immuno type according to the inflammation level
perc_roi_infl$immuno_cat= sapply(1:nrow(perc_roi_infl), function(n){
  
  x=perc_roi_infl$inf_tumor[n]
  y=perc_roi_infl$inf_stroma[n]
  if(x>=0.7){
    return("purely_inflammed")
  }
  else if(x<=0.699999 & x>=0.5){
    return("mixed_inflammed")
  }
  else if(x<=0.499999 & y>=0.1){
    return("excluded")
  }
  else if(x<=0.1 & y<0.1){
    return("desert")
  }
  else{return("mixed_inflammed")}
  
})


############### Triplets code
#### For the code to work we need a single column for the phenotypes
replace_entries <- function(x) {
  ifelse(grepl("-", x), "", x)
}

replace_na <- function(x) {
  ifelse(is.na(x), "", x)
}


### The columns "Phenotype.XX" (imagine: Phenotype.CD8) containts the following entries:
### "XX+"  "XX-", we convert every "XX-".  With the following command we convert every entry that has a "-" into an empty string "".
df <- mutate_at(all_df, vars(starts_with("Phenotype")), ~replace_entries(.))

#### note that there are several entries on the phenotype columns that are NA, we replace them with an empty string
df <- mutate_at(df, vars(starts_with("Phenotype")), ~replace_na(.))


### we paste together all the phenotypes
df$phenotype=paste0(df$Phenotype.CD8, df$Phenotype.CD11c, df$Phenotype.CD68, df$Phenotype.CK, df$Phenotype.PD1, df$Phenotype.PDL1)

### we convert the [+] with a "p_",  this is done because it helps another routine to select a specific phenotype.
### note that this add an underscore at the end, we will take care of this later
df$phenotype=str_replace_all(df$phenotype, "[+]", "p_")

### empty phenotypes are DAPIp,  we have to add an underscore for the moment
df$phenotype[df$phenotype==""]="DAPIp_"

###  in all the entries we have an additional underscore at the end, we remove it
df$phenotype=substr(df$phenotype, 1, nchar(df$phenotype)-1)


### We perform the analysis on triple niches.
### In this case we start from CD8p_PD1neg_total, there MUST be a CD8p_PD1p (pure), then we need also
### any APC cell.   An APC is either a CD11cp (pure) or a CD68p (pure).

### for this purpose we need to rename CD11cp and CD68p PURE as APCp
### we perform this on a copy of the datraframe "df"
df2=df
df2$phenotype[df2$phenotype=="CD11cp" | df2$phenotype=="CD68p"]="APCp"

### the triplets routine assume a dataframe with the following columns: sample, nucleus.x, nucleus.y, phenotype
df2=df2 %>% select(tag, Annotation.ID, Tissue.Category, Cell.X.Position, Cell.Y.Position, phenotype)
df2= df2 %>% rename(sample=tag, nucleus.x=Cell.X.Position, nucleus.y=Cell.Y.Position, tissue.type=Tissue.Category)

### define the triplets we want to investigate
df_name_points=data.frame(Var1="CD8p_PD1neg_total", Var2="CD8p_PD1p", Var3="APCp")

### call the routine that compute the mutual interaction of the triplets
### NOTE: if there are no interactions at all the sample will NOT be reported
res_mutual=calc_triplets_mutual(df2, threshold=100, df_name_points)




