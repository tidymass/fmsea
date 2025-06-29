no_function()

###20220307

# object mass_dataset class
# metabolite_set is the pathway database

library(massdataset)
library(metid)
library(metpath)

setwd(masstools::get_project_wd())
setwd("demo_data/denmark_project/metabolome/peaks/")

load("expression_data")
load("sample_info")
load("variable_info")
load("phenotype_info")

idx <- c("subject_id", setdiff(colnames(phenotype_info), colnames(sample_info)))

phenotype_info <-
  phenotype_info[,idx]

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_info, by = "subject_id")

sample_info <- 
  sample_info %>% 
  dplyr::mutate(class = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    stringr::str_detect(sample_id, "blk") ~ "Blank",
    TRUE ~ "Subject"
  ))

object <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info %>% dplyr::select(-na_freq)
  )

sum(is.na(object))

qc_id <- 
  object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(class == "QC") %>% 
  pull(sample_id)

subject_id <- 
  object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(class == "Subject") %>% 
  pull(sample_id)

object <-
  object %>% 
  mutate_variable_na_freq(according_to_samples = qc_id) %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  rename(na_freq_qc = na_freq) %>% 
  mutate_variable_na_freq(according_to_samples = subject_id) %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  rename(na_freq_subject = na_freq)

###remove NA > 20 in QC
object <- 
  object %>% 
  filter(na_freq_qc < 0.2)

###remove samples which are after delivery
object <- 
  object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  filter(!postdelivery)


###calculate the correlation between ga and metabolic features
plot(density(unlist(object[1,,drop = TRUE])))
plot(density(unlist(object[4,,drop = TRUE])))

object <- 
  object %>% 
  `+`(1) %>% 
  log(2)

plot(density(unlist(object[4,,drop = TRUE])))

sum(is.na(object))

###remove QC samples
object <-
  object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(!is.na(g_stage))

# ###linear regression
ga <- object@sample_info$g_stage

pheno <-
  data.frame(sample_id = colnames(object),
             stringsAsFactors = FALSE) %>%
  dplyr::left_join(object@sample_info, by = "sample_id")

age <-
  pheno$age

sum(is.na(age))

bmi <-
  pheno$pre_bmi

sum(is.na(bmi))

parity <-
  pheno$earlier_preg_with_children

sum(is.na(parity))

library(furrr)
library(future)

plan(future::multiprocess())

p_cor <-
  furrr::future_map(
    as.data.frame(t(object@expression_data)),
    .f = function(x) {
      temp_data <-
        data.frame(ga,
                   x,
                   age,
                   bmi,
                   parity,
                   stringsAsFactors = FALSE)

      adjusted_x = lm(data = temp_data, formula = x ~ age + bmi + parity)
      adjusted_x = adjusted_x$residuals
      temp_data =
        data.frame(temp_data, adjusted_x, stringsAsFactors = FALSE)

      temp = cor.test(temp_data$adjusted_x, temp_data$ga, method = "spearman")
      return(c(p = temp$p.value, cor = unname(temp$estimate)))
    },
    .progress = TRUE
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

p_cor <-
  p_cor %>%
  tibble::rownames_to_column(var = "variable_id")

p_cor %>%
  dplyr::filter(cor > 0.8)

plot(ga, as.numeric(object["M286T556_POS",,drop = TRUE]))


p_cor %>%
  dplyr::filter(cor < -0.7)

plot(ga, as.numeric(object["M465T448_POS",,drop = TRUE]))


object <- 
  object %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::mutate(cor = p_cor$cor,
                cor_p = p_cor$p)

save(object, file = "object")

# library(metID)
# 
# annotation_pos =
#   metID::identify_metabolites(
#     ms1.data = "ms1_data_pos.csv",
#     ms1.match.ppm = 25,
#     database = "hmdbMS1Database0.0.1",
#     polarity = "positive",
#     threads = 5,
#     candidate.num = 1000
#   )
# save(annotation_pos, file = "annotation_pos")
load("annotation_pos")

# annotation_neg =
#   metID::identify_metabolites(
#     ms1.data = "ms1_data_neg.csv",
#     ms1.match.ppm = 25,
#     database = "hmdbMS1Database0.0.1",
#     polarity = "negative",
#     threads = 5,
#     candidate.num = 1000
#   )
# save(annotation_neg, file = "annotation_neg")
load("annotation_neg")













