library(tidyverse)
library(vegan)

vtable = read_csv("../generated_data/virus_table_filtered.csv") %>%
    select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) # drop empty columns

vmeta = read_csv("../raw_data/virus_meta.csv")
lib_meta = read_csv("../raw_data/lib_metadata.csv")

virus_name_arr = colnames(select(vtable, -lib_id))
lib_id_arr = vtable$lib_id

vmat = t(select(vtable, all_of(virus_name_arr)))
colnames(vmat) = vtable$lib_id

vmat_df = as.data.frame(vmat) %>%
    mutate(virus_name = rownames(.), .before=1) %>%
    left_join(vmeta, by="virus_name")

# viral load
total_RPM = colSums(vmat)
total_RPM_df = data.frame(lib_id = names(total_RPM), viral_load=log10(total_RPM+1))

# count the numbers of virus of concern species per individual
voc_df = vmat_df %>%
    group_by(virus_of_concern) %>% 
    summarise(across(all_of(lib_id_arr), function(x){sum(x>0)})) %>%
    select(-virus_of_concern) %>%
    t() %>%
    as.data.frame()
colnames(voc_df) = c("n_non_voc", "n_voc")
voc_df$lib_id = rownames(voc_df)

# Shannon index
shannon = diversity(t(vmat))
shannon_df = data.frame(lib_id = names(shannon), shannon=shannon)

# count virus family incidence
vfi_df = vmat_df %>%
    group_by(virus_family) %>% 
    summarise(across(all_of(lib_id_arr), function(x){ifelse(sum(x>0)>0, 1, 0)}))
vf_arr = vfi_df$virus_family

vfi_df = vfi_df %>% select(-virus_family) %>%
    t() %>%
    as.data.frame()
colnames(vfi_df) = vf_arr
vfi_df$lib_id = rownames(vfi_df)

lib_meta_new = lib_meta %>%
    left_join(total_RPM_df, by="lib_id") %>%
    left_join(voc_df, by="lib_id") %>%
    mutate(n_virus = n_non_voc + n_voc) %>%
    left_join(shannon_df, by="lib_id") %>%
    left_join(vfi_df, by="lib_id")

write_csv(lib_meta_new, "../generated_data/lib_meta_w_div.csv")

lib_meta_new %>% 
    select(genus, n_virus, shannon) %>% 
    write_csv("../output/alpha_diversity.csv")

glm_fit = glm(n_virus ~ site + genus, family="poisson", data=lib_meta_new)
drop1(glm_fit, test="Chi")

glm_optimal = glm(n_virus ~ genus, family="poisson", data=lib_meta_new)

glm_pred = predict(glm_optimal, terms="genus", se.fit=T, type="terms" )

pred_df = tibble(
    lib_id = lib_meta_new$lib_id,
    genus = lib_meta_new$genus,
    fitted = glm_pred$fit[,1],
    se = glm_pred$se.fit[,1]) %>%
    group_by(genus) %>%
    summarise(fitted = mean(fitted), se = 1.98*mean(se))

write_csv(pred_df, "../output/nvirus_glm_effect_size.csv")

# evaluate prevalence of each virus family
lib_meta_new %>%
    group_by(genus) %>% 
    summarise(n_libs=n(), across(all_of(unique(vmat_df$virus_family)), sum)) %>%
    write_csv("../output/virus_family_by_genus.csv")

lib_meta_new %>%
    group_by(species) %>% 
    summarise(n_libs=n(), across(all_of(unique(vmat_df$virus_family)), sum)) %>%
    write_csv("../output/virus_family_by_species.csv")

lib_meta_new %>%
    count(species, site) %>%
    arrange(site) %>%
    write_csv("../output/sample_stat.csv")
