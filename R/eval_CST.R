library(tidyverse)

vtable = read_csv("../generated_data/virus_table_filtered.csv") %>%
    select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) # drop empty columns

vmeta = read_csv("../raw_data/virus_meta.csv")
lib_meta = read_csv("../raw_data/lib_metadata.csv") %>%
    mutate(species = str_replace(species, " ", "_"))

virus_name_arr = colnames(select(vtable, -lib_id))
lib_id_arr = vtable$lib_id

num_host_species = vtable %>%
    left_join(select(lib_meta, lib_id, species), by="lib_id") %>%
    select(-lib_id) %>%
    group_by(species) %>%
    summarise_all(function(x) {ifelse(sum(x)>0, 1, 0)}) %>%
    select(-species) %>%
    colSums()

num_host_genus = vtable %>%
    left_join(select(lib_meta, lib_id, genus), by="lib_id") %>%
    select(-lib_id) %>%
    group_by(genus) %>%
    summarise_all(function(x) {ifelse(sum(x)>0, 1, 0)}) %>%
    select(-genus) %>%
    colSums()

virus_stat = tibble(
    virus_name = names(num_host_species),
    num_host_species = num_host_species,
    num_host_genus = num_host_genus,
    is_multi_host_s = num_host_species > 1,
    is_multi_host_g = num_host_genus > 1
) %>%
    left_join(vmeta, by="virus_name")

write_csv(virus_stat, "../output/virus_stat.csv")

#### Network analysis ####
require(igraph)
vh_g = vtable %>%
    left_join(lib_meta, by="lib_id") %>%
    select(species, all_of(virus_name_arr)) %>%
    pivot_longer(cols=all_of(virus_name_arr), names_to="virus_name", values_to="RPM") %>%
    filter(RPM > 0) %>%
    select(-RPM) %>%
    table() %>%
    as.matrix() %>%
    graph.incidence(weighted = T) 
bp_g = bipartite.projection(vh_g)
v_g = bp_g$proj2
h_g = bp_g$proj1

vm = tibble(virus_name=names(V(vh_g))) %>% 
    left_join(vmeta) %>% 
    mutate(virus_name_abbr = ifelse(is.na(virus_name_abbr), virus_name, virus_name_abbr)) %>%
    replace_na(list(virus_of_concern="Host"))
cm = c("No"="green", "Yes"="red", "Host"="blue")
pdf("../output/network.pdf", width=10, height=10)
plot(vh_g, vertex.color=cm[vm$virus_of_concern], vertex.label=vm$virus_name_abbr, vertex.size=8)
dev.off()
vb = cluster_edge_betweenness(v_g)

#### GLM model of virus sharing between hosts ####
host_tree = ape::read.tree("../raw_data/bat_tree_timetree.nwk")
host_dist = cophenetic(host_tree)/2
host_dist_df = host_dist %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("host_sp1", "host_sp2"), value.name = "host_dist")

host_tree_COI = ape::read.tree("../raw_data/bat_tree_COI.nwk")
host_dist_COI = cophenetic(host_tree_COI)
host_dist_df_COI = host_dist_COI %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("lib1", "lib2"), value.name = "host_dist_COI")


site_loc = read_csv("../raw_data/geo_locations.csv") %>%
    select(-site_name) %>%
    as.data.frame()
rownames(site_loc) = site_loc$site
site_loc$site = NULL
site_xy = SoDA::geoXY(site_loc$lat, site_loc$lon, unit=1000) # unit=km
rownames(site_xy) = rownames(site_loc)
site_dist_df = dist(site_xy) %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("site1", "site2"), value.name = "geo_dist")

time_dist = as.matrix(dist(lib_meta$year))
rownames(time_dist) = lib_meta$lib_id
colnames(time_dist) = lib_meta$lib_id
time_dist_df = time_dist %>%
    reshape2::melt(varnames = c("lib1", "lib2"), value.name = "time_dist")

lmeta = select(lib_meta, lib_id, species, site)

tbl_host_sp = vtable %>%
    left_join(select(lib_meta, lib_id, species), by="lib_id") %>%
    pivot_longer(cols=all_of(virus_name_arr), names_to="virus_name", values_to="RPM") %>%
    filter(RPM > 0) %>%
    distinct(species, virus_name) %>%
    table() 
vs_host_sp_df = proxy::dist(tbl_host_sp, function(x, y) {sum((x>0)&(y>0))}) %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("host_sp1", "host_sp2"), value.name = "shared_virus") %>%
    left_join(host_dist_df, c("host_sp1", "host_sp2")) %>%
    filter(host_sp1 != host_sp2)

ggplot(vs_host_sp_df, aes(x=host_dist, y=shared_virus)) +
    geom_jitter(width=0.5, height=0.1, alpha=0.7) +
    geom_smooth(method="glm", formula=y~x, method.args=list(family="poisson")) +
    xlab("Host divergence time (Mya)") +
    ylab("Shared number of virus") +
    theme_bw()
ggsave("../output/shared_virus_by_host_sp.pdf", width=10.1, height=8.9, units="cm")

tbl = vtable %>%
    pivot_longer(cols=all_of(virus_name_arr), names_to="virus_name", values_to="RPM") %>%
    filter(RPM > 0) %>%
    select(-RPM) %>%
    table() 

vs_df = proxy::dist(tbl, function(x, y) {sum((x>0)&(y>0))}) %>%
    as.matrix() %>%
    reshape2::melt(varnames = c("lib1", "lib2"), value.name = "shared_virus") %>%
    left_join(select(lmeta, lib_id, site1=site, host_sp1=species), by=c("lib1"="lib_id")) %>%
    left_join(select(lmeta, lib_id, site2=site, host_sp2=species), by=c("lib2"="lib_id"))

tmp = vs_df %>%
    left_join(host_dist_df, c("host_sp1", "host_sp2")) %>%
    left_join(host_dist_df_COI, c("lib1", "lib2")) %>%
    left_join(site_dist_df, c("site1", "site2")) %>%
    left_join(time_dist_df, c("lib1", "lib2")) %>%
    filter(lib1 != lib2)

model_terms = list(
    host_c = c(NA, "host_dist", "host_dist_COI"),
    site_c = c(NA, "geo_dist"),
    time_c = c(NA, "time_dist")
)
models_comb = expand.grid(model_terms)

f = function(row) {
    c = paste(na.omit(row), collapse = ' + ')
    paste0(ifelse(c=="", "1", c))
}
models = apply(models_comb, 1, f) %>% 
    data.frame(Formula = ., stringsAsFactors = F) %>% 
    mutate(Formula = paste('shared_virus ~', Formula))

fitted = models %>%
    rowwise() %>%
    mutate(Fit = list(glm(as.formula(Formula), family="poisson", data=tmp)))

sel = MuMIn::model.sel(fitted$Fit)
supported_models = filter(sel, delta<2)


M0 = glm(shared_virus ~ host_dist_COI + geo_dist, family="poisson", data=tmp)
M1 = glm(shared_virus ~ host_dist_COI + geo_dist + time_dist, family="poisson", data=tmp)

np = 601
host_x = (0:(np-1))/(np-1)
geo_x = (0:(np-1))/(np-1) * 600
host_pred = predict(M0, type="response", se.fit = T, newdata = tibble(host_dist_COI=host_x, geo_dist=rep(0,np)))
geo_pred =  predict(M0, type="response", se.fit=T, newdata = tibble(host_dist_COI=rep(0,np), geo_dist=geo_x))

x_h = c(host_x,rev(host_x))
y_h = c(host_pred$fit + 2*host_pred$se.fit, rev(host_pred$fit-2*host_pred$se.fit))
g1 = ggplot() + 
    
    geom_jitter(aes(x=tmp$host_dist_COI, y=tmp$shared_virus), height=0.1, alpha=0.9, color="gray", shape=1) +
    geom_polygon(aes(x_h,y_h), fill="blue", alpha=0.5) +
    geom_line(aes(x=host_x, y=host_pred$fit)) +
    xlab("Phylogenetic distance") +
    ylab("Shared num. of virus") +
    theme_bw()

x_g = c(geo_x,rev(geo_x))
y_g = c(geo_pred$fit + 2*geo_pred$se.fit, rev(geo_pred$fit-2*geo_pred$se.fit))
g2 = ggplot() + 
    
    geom_jitter(aes(x=tmp$geo_dist, y=tmp$shared_virus), height=0.1, alpha=0.9, color="gray", shape=1) +
    geom_polygon(aes(x_g,y_g), fill="blue", alpha=0.5) +
    geom_line(aes(x=geo_x, y=geo_pred$fit)) +
    xlab("Geographic distance (km)") +
    ylab("Shared num. of virus") +
    theme_bw()
    
cowplot::plot_grid(g1, g2, nrow=2)
ggsave("../output/glm_est_host_and_geo.pdf", width=3, height=4.6)