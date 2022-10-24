# This script do the RPM calculation and index-hopping correction
# The assumed work dir is <project_root>/R

library(tidyverse)

reads_mapping = read_csv("../raw_data/raw_reads_mapping.csv")
lane_id_table = read_csv("../raw_data/lane_id_table.csv")

# Transpose the table
virus_name_arr = reads_mapping$virus_name
lib_id_arr = colnames(select(reads_mapping, -virus_name))

reads_mapping_t = reads_mapping %>%
    select(-virus_name) %>%
    t() %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(lib_id = lib_id_arr, .before=1)

colnames(reads_mapping_t) = c("lib_id", virus_name_arr)

# Merge lane id info into the transposed table
reads_mapping_t_lane = inner_join(lane_id_table, reads_mapping_t, by="lib_id")
# check if the two tables are correctly joined
flag = setequal(reads_mapping_t_lane$lib_id, lib_id_arr) 
if (!flag) { warning("join failed.") }

reads_sum_by_lane = reads_mapping_t_lane %>%
    group_by(lane_id) %>%
    summarise(across(all_of(virus_name_arr), sum))

# set threshold to reduce false positive due to index hopping
threshold_matrix = reads_mapping_t_lane %>%
    select(lib_id, lane_id) %>%
    left_join(reads_sum_by_lane, by="lane_id") %>%
    select(all_of(virus_name_arr)) %>%
    mutate_all(function(x) {x*0.001}) # threshold value set at 0.1%

raw_reads_matrix = select(reads_mapping_t_lane, all_of(virus_name_arr))
mask_matrix = raw_reads_matrix > threshold_matrix

filtered_count_lane = sum((raw_reads_matrix <= threshold_matrix) & (raw_reads_matrix != 0))
warning(str_glue("{filtered_count_lane} possible false positives deleted (index hopping filter)"))

filtered_reads_matrix = mask_matrix * raw_reads_matrix

reads_mapping_t_lane_filtered = bind_cols(
    select(reads_mapping_t_lane, lib_id, lane_id), 
    filtered_reads_matrix)

# load norRNA reads count
norRNA_reads = read_tsv("../raw_data/norRNA.readcount")
lane_filtered_with_reads = inner_join(norRNA_reads, reads_mapping_t_lane_filtered, by="lib_id")
# check if the two tables are correctly joined
flag = assertthat::are_equal(nrow(lane_filtered_with_reads), nrow(reads_mapping_t_lane_filtered)) 
if (!flag) { warning("join failed.") }

# calc RPM and apply RPM > 1 filter
virus_table_filtered = lane_filtered_with_reads %>%
    mutate(across(all_of(virus_name_arr), function(x) {x/norRNA_reads * 1e6})) %>%
    mutate(across(all_of(virus_name_arr), function(x) {ifelse(x>1, x, 0)})) %>%
    select(-norRNA_reads, -lane_id)

RPM_mat = lane_filtered_with_reads %>% 
    mutate(across(all_of(virus_name_arr), function(x) {x/norRNA_reads * 1e6})) %>% 
    select(all_of(virus_name_arr))
filtered_count_RPM = sum(RPM_mat >0 & RPM_mat <= 1)
warning(str_glue("{filtered_count_RPM} possible false positives deleted (RPM filter)"))

# output
write_csv(virus_table_filtered, "../generated_data/virus_table_filtered.csv")