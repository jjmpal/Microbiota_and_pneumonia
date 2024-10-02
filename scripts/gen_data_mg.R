
#--------------------#
# Generate the data  #
#--------------------#

#Generates the data for the 'HELIUS microbiota and infections' replication study
#This reads in metagenome data. Requires finnriskR_env and R4.3


#---------------------------------------------------------------------------
#Set up

#Libraries
library(data.table)
library(microbiome)
suppressPackageStartupMessages(library(tidyverse))


#I was not able to install yingtools2 dependence due to outdated environment.
#However some of their functions were usable. 
source("scripts/functions_yingtools2.R")

#Paths from external file
path.list <- fread("paths.txt") %>% 
  group_by(path_name) %>%
  transmute(named_vec = list(path)) %>%
  deframe()

#Output file names
data <- "data_mg"
phy_spe_out        <- str_glue("{data}/pseq_species.Rds")
phy_all_out        <- str_glue("{data}/pseq_all.Rds")
phy_gen_out        <- str_glue("{data}/pseq_genera.Rds")

beta_dist_out <- str_glue("{data}/beta_dist.Rds")
beta_ord_out  <- str_glue("{data}/beta_ord.Rds")
df_final_out  <- str_glue("{data}/df_final.Rds")


# Outcomes
outcomes_nicenames <- c(BL1ST_J10_PNEUMONIA = "Pneumonia", BL1ST_INPAT_INFECT_HELIUS_LOWER_RESPIRATORY_TRACT = "Pneumonia_helius", BL1ST_PNEUMONIA_HELIUS = "Pneumonia_helius_new")
outcomes <- names(outcomes_nicenames)
outcomes_age <- paste0(outcomes, "_AGE")
#outcomes_year <- paste0(outcomes, "_YEAR")

#covariates
covs_nicenames <- 
  c(BL_AGE            = "Age",                #Age at the time of sample collection
    MEN               = "Men",                #Is person a men
    BMI               = "BMI",
    SMOKING           = "Smoking",            #Three levels: Current smoker: 2, ex-smoker: 1, never: 0 ; self defined
    ALKI2_FR02        = "Alchol",             #Alcohol usage     ##Confirm that is OK to use this!!
    Q57X              = "Physical activity",  #Physical activity
    BL_USE_RX_J01     = "Antibiotics",        #Prior antibiotics
    PREVAL_DIAB       = "Diabetes",           #Prevalent diabetes
    PREVAL_CVD        = "CVD",                #Prevalent CVD
    PREVAL_CR_ANYCANC = "Cancer",             #Prevalent cancer
    HYPERT_AHA        = "Hypertension",       #Hypertension at baseline
    PREVAL_PULMONARY  = "Pulmonary",          #PREVAL_COPD == 1 | PREVAL_ASTHMA == 1; self defined
    PREVAL_GASTRO     = "Gastrointestinal",   #PREVAL_LIVERDIS == 1 | PREVAL_IBD == 1 | PREVAL_K11_COELIAC ==1; self defined
    GRAVID            = "Gravid")             #2=pregnant, 1=no (only women)
covs_names <- names(covs_nicenames)

#Buryrate.producers
butyrate.producers <- c("Butyricimonas", "Odoribacter", "Anaerostipes", "Anaerobutyricum", "Agathobacter", "Butyrivibrio", "Coprococcus", "Roseburia","Shuttleworthia","Butyricicoccus", "Faecalibacterium","Flavonifractor","Pseudoflavonifractor","Oscillibacter", "Subdoligranulum")
butyrate.producers.species <- c("Subdoligranulum variabile", "Eubacterium_G ventriosum")

#-----------------------------------------------------------------------
#Load

# Load data
df0 <- fread(path.list$pheno_file)
df0_old <- fread(path.list$pheno_file_old) %>% 
  select(Barcode, DEATH_INFECT_HELIUS, any_of(outcomes),any_of(outcomes_age))
covs0 <-  fread(path.list$covs_file)
P.all <- readRDS(path.list$phymg_file) %>% mia::makePhyloseqFromTreeSummarizedExperiment()

#Number of individuals in files
print("Number of individuals in new and old phenotype file, covariate file and phyloseq file, respectively.")
dim(df0)[1]
dim(df0_old)[1]     
dim(covs0)[1]
P.all
print("\n")

#Test with subset of data
#randsamples <- sample(colnames(otu_table(P.all)), 50)
#P.all <- prune_samples(randsamples, P.all)
#P.all

#--------------------------------------------------------------------------
#Preprocess phenotype data

# Covariates
covs <- covs0 %>% 
  mutate(
    #Covariates, by condition
    SMOKING = case_when(CURR_SMOKE==1 ~ 2L, 
                        EX_SMOKE ==1 ~ 1L, 
                        EX_SMOKE ==0 & CURR_SMOKE ==0 ~ 0L ),
    #alcohol = if_else(ALKI2_FR02 > 0, 1L, 0L),  #Isn't the continuous variable better?
    PREVAL_PULMONARY = if_else(PREVAL_COPD == 1 | PREVAL_ASTHMA == 1L, 1L, 0L),
    PREVAL_GASTRO = if_else(PREVAL_LIVERDIS == 1 | PREVAL_IBD == 1 | PREVAL_K11_COELIAC ==1, 1L, 0L )
  ) %>% 
  select("Barcode", all_of(covs_names)) %>%
  mutate_if(is.integer, as.factor)

# Other variables

df <- df0 %>% 
  select(Barcode, contains("DEATH"), any_of(outcomes), any_of(outcomes_age)) %>%
  filter(Barcode!="") 
  
print("Number of individuals at phenotype file, FR02 and FR07")
dim(df0)[1]
print("Number of individuals with Barcode (Microbiome data), only in FR02")
dim(df)[1]

df <- df %>%
  left_join(df0_old, by="Barcode") %>%
  left_join(covs, by="Barcode") %>%
  #Calculate time differences
  mutate(DEATH_AGEDIFF = DEATH_AGE - BL_AGE) %>%
  mutate_at(outcomes_age, list(DIFF =~ . - BL_AGE)) %>%
  rename_at(vars(ends_with("AGE_DIFF")), ~str_replace(., "AGE_DIFF","AGEDIFF")) %>%
  #Define variable for competing risk analysis
  mutate_at(outcomes, 
            list(EVENTDEATH =~ case_when(. == 1 ~ "event",
                                         DEATH == 1 ~ "death",
                                         TRUE ~ "no event"))) %>%
  mutate_at(vars(contains("EVENTDEATH")), ~fct_relevel(.,"no event", "event", "death")) %>%
  mutate_at(outcomes, as.factor)

print('Number of individuals with all phenotype data combined')
dim(df)[1]
print('Remove all individiduals with missing covariate values')
df <- df %>% drop_na(all_of(covs_names[covs_names!="GRAVID"]))
dim(df)[1] 
print('Remove all pregnant individuals')
df <- df %>% filter(GRAVID==1 | is.na(GRAVID))
dim(df)[1]

print('Summary: processed phentype data')
summary(df)

#--------------------------------------------------------------------------
#Preprocess microbiome data

# Final phyloseq object let us use species level features
print("Agglomerating to species level...")
set.seed(42)
P <- phyloseq::tax_glom(P.all, "Species") 
print("P, after agglomeration")
P

# Set row names
df <- df %>% column_to_rownames("Barcode") # also converts to data.frame; below did not work for tibble
#df <- as.data.frame(df)
#rownames(df) <- df$Barcode


# Match statistics
print("Data frame entries matching phyloseq")
print(
  c(mean=mean(rownames(df) %in% colnames(otu_table(P))),
    N=sum(rownames(df) %in% colnames(otu_table(P))))
)

print("phyloseq matching data frame entries")
print(
  c(mean=mean(colnames(otu_table(P)) %in% rownames(df)),
    N=sum(colnames(otu_table(P)) %in% rownames(df)))
)

# Subset to samples that we have in phenodata and in P
coms <- intersect(rownames(df), colnames(otu_table(P)))

# Also require more than N reads (there are some very low read count samples; suitable cutoff is around 1000-2000 reads, 2-6% of samples)
coms <- intersect(coms, names(which(colSums(otu_table(P)) > 1000)))

P <- prune_samples(coms, P)
P.all <- prune_samples(coms, P.all)
sample_data(P) <- df[rownames(sample_data(P)), ]
sample_data(P.all) <- df[colnames(otu_table(P.all)), ]
print("P, after pruning")
P

df <- sample_data(P) %>% as_tibble()
df$Barcode <- rownames(sample_data(P))

print('df dimensions after including only individuals in physeq data')
print(df %>% dim())

#-------------------------------------------------------------------------------------------------
# Agglomerate to Genus level

print("Agglomerating to Genus level...")
set.seed(42)
P.genus <- phyloseq::tax_glom(P.all, "Genus")
print("P.genus, after agglomeration:")
P.genus

#------------------------------------------------------------------------------------------------
# Butyrate producers

#We quantified the relative abundance of butyrate-producing bacteria by
#calculating the relative abundance of 17 bacteria that are known to be
#the most abundant drivers of butyrate production.  This might take a
#while to run.

print("Abundance of butyrate producers...")
  
#Clean the tax_table so that we can recognise butyrate producers and create new 'taxon level'
tax.table.buty <- as.data.frame(tax_table(P)) %>% 
  mutate_all(~str_remove(., ".__")) %>% 
  mutate_at("Species", ~str_replace(., "([^ ^_]+)_.* ([^_]*).*", "\\1 \\2")) %>%
  mutate_at("Genus",   ~str_replace(., "^([^_]+)_.*$", "\\1")) %>%
  mutate(butyrateproducers = if_else(Genus %in% butyrate.producers, Genus, 
                                     if_else(Species %in% butyrate.producers.species, Species, "Other"))) 
P.buty <- P
tax_table(P.buty) <- as.matrix(tax.table.buty)

print("Butyrate producers:")
as.data.table(tax.table.buty) %>%
  filter(butyrateproducers != "Other") 

#AK: get.otu.melt is function from problematic 'yingtools2' package. 
#Calculate abundance
butyrate <- get.otu.melt(P.buty, filter.zero = F) %>%
  filter(butyrateproducers != "Other") %>%
  group_by(sample)%>%                       
    summarize(sum = sum(pctseqs)) %>%
    dplyr::rename(Barcode=sample)               

#Add to data frame and create tertiles
df <- left_join(df, butyrate) %>%
  mutate(butyrate = sum*100) %>% # convert fractions to percents
  mutate(tertiles_butyrate = ntile(`sum`, 3)) %>%
  mutate(tertiles_butyrate = as.factor(if_else(tertiles_butyrate == 1, 'Low butyrate', if_else(tertiles_butyrate == 2, 'Intermediate butyrate', 'High butyrate')))) %>%
  mutate(tertiles_butyrate = fct_relevel(tertiles_butyrate, "Low butyrate", "Intermediate butyrate", "High butyrate")) %>%
  dplyr::select(-sum)

#CLR-transformation for butyrate producers
#We performed a CLR-transformation to correct for the compositional nature of microbiome data
P.buty.clr <- microbiome::aggregate_taxa(P.buty, level = "butyrateproducers")
P.buty.clr <- microbiome::transform(P.buty.clr, "clr")

butyrate.clr <- get.otu.melt(P.buty.clr, filter.zero = F) %>%
  filter(butyrateproducers != "Other") %>%
  group_by(sample)%>%
    summarize(butyrate.clr = sum(numseqs)) %>% # calculate the CLR-transformed abundance of the butyrate-producing bacteria per sample
    mutate(tertiles_butyrate_clr = ntile(`butyrate.clr`, 3)) %>%
    mutate(tertiles_butyrate_clr = as.factor(if_else(tertiles_butyrate_clr == 1, 'Low butyrate', if_else(tertiles_butyrate_clr == 2, 'Intermediate butyrate', 'High butyrate')))) %>%
    mutate(tertiles_butyrate_clr = fct_relevel(tertiles_butyrate_clr, "Low butyrate", "Intermediate butyrate", "High butyrate")) %>%
    dplyr::rename(Barcode=sample)

df <- left_join(df, butyrate.clr)

print("Dimensions and column names")
df %>% dim()
df %>% names()

#-------------------------------------------------------------------------------------------------
#Alpha diversity

#Calculate alpha diversity
print("Calculating alpha diversity...")

set.seed(88)
P.alpha <- rarefy_even_depth(prune_samples(sample_names(P), P.all), rngseed = 88)      
P.alpha <- microbiome::aggregate_taxa(P.alpha, level = "Species") # Genus or ta7?
alpha <- estimate_richness(P.alpha, measures="Shannon")
alpha$Barcode <- row.names(alpha)                               #AK sample -> Barcode in all occurances
alpha$Barcode <- gsub("X","",as.character(alpha$Barcode))
alpha$Barcode <- str_replace(alpha$Barcode, "\\.","-")        #AK 
alpha <- alpha %>% dplyr::select(Barcode, Shannon)

df <- left_join(df, alpha) %>%
  filter(!is.na(Shannon)) %>% # Added by AK - some individuals did not have microbiome data, even they had 'Barcode'
  mutate(tertiles = ntile(`Shannon`, 3)) %>%
  mutate(tertiles = as.factor(if_else(tertiles == 1, 'Low diversity', if_else(tertiles == 2, 'Intermediate diversity', 'High diversity')))) %>%
  mutate(tertiles = fct_relevel(tertiles, "Low diversity", "Intermediate diversity", "High diversity"))

print("Dimensions and column names")
df %>% dim()
df %>% names()


#--------------------------------------------------------------------------------------------------
# Beta diversity

print("Calculating Beta diversity...")

#Calculate beta diversity
P.comp <- microbiome::transform(P, "compositional")
set.seed(88)
dist <- phyloseq::distance(P.comp, method = "bray") 
ord <- ordinate(P.comp, method = "PCoA", distance = "bray")

print("Summary, ordination")
ord %>% summary()

#-------------------------------------------------------------------------------------------------
#Write out

print("The final df: Summary")
df %>% summary()

print("Writing outputfiles...")

#Write out microbiome data
saveRDS(P, phy_spe_out)
saveRDS(P.all, phy_all_out)
saveRDS(P.genus, phy_gen_out)
#Beta diversity
saveRDS(dist, beta_dist_out) 
saveRDS(ord, beta_ord_out) 
#Write out the processed data frame
saveRDS(df, df_final_out)



