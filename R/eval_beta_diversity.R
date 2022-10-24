library(tidyverse)
library(ecodist)

vtable = read_csv("../generated_data/virus_table_filtered.csv") %>%
    select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) # drop empty columns

vmeta = read_csv("../raw_data/virus_meta.csv")
lib_meta = read_csv("../raw_data/lib_metadata.csv")

virus_name_arr = colnames(select(vtable, -lib_id))
lib_id_arr = vtable$lib_id

vmat = t(select(vtable, all_of(virus_name_arr)))
colnames(vmat) = vtable$lib_id

nonempty_libs = names(which(colSums(vmat)>0))
vmat_ne = vmat[,nonempty_libs]
vmat_vegan = t(vmat_ne)

pcoa = cmdscale(vegan::vegdist(vmat_vegan, method="bray"), eig=T, k=3)

host_genus = filter(lib_meta, lib_id %in% nonempty_libs)$genus
sample_site = filter(lib_meta, lib_id %in% nonempty_libs)$site
sample_year = filter(lib_meta, lib_id %in% nonempty_libs)$year
# draw PCoA
ggplot() + 
    geom_point(aes(x=pcoa$points[,1] , y=pcoa$points[,2], color=host_genus, shape=sample_site), size=4) +
    xlab(paste0("PcoA1 (", round(pcoa$eig[1]/sum(pcoa$eig) * 1000) /10, "%)")) + 
    ylab(paste0("PCoA2 (", round(pcoa$eig[2]/sum(pcoa$eig) * 1000) /10, "%)")) + 
    theme_bw()

ggsave("../output/virus_pcoa.pdf", width=6.24, height=4.58)

# perform PERMANOVA
dbrda_fit = vegan::dbrda(vmat_vegan ~ sample_site + host_genus + sample_year, distance = "bray") 
anova(dbrda_fit, by="margin") %>%
    capture.output(file="../output/dbRDA_bray_marginal_effect.txt")
anova(dbrda_fit, by="term") %>%
    capture.output(file="../output/dbRDA_bray_type_I_SS.txt")

# perform Mantel test
site_loc = read_csv("../raw_data/geo_locations.csv") %>%
    select(-site_name) %>%
    as.data.frame()
rownames(site_loc) = site_loc$site
site_loc$site = NULL

site_xy = SoDA::geoXY(site_loc$lat, site_loc$lon, unit=1000) # unit=km
rownames(site_xy) = rownames(site_loc)
site_dist = as.matrix(dist(site_xy))

host_tree = ape::read.tree("../raw_data/bat_tree_timetree.nwk")
host_dist = cophenetic(host_tree)/2
host_tree_COI = ape::read.tree("../raw_data/bat_tree_COI.nwk")
host_dist_COI = cophenetic(host_tree_COI)

time_dist = as.matrix(dist(lib_meta$year))
rownames(time_dist) = lib_meta$lib_id
colnames(time_dist) = lib_meta$lib_id

bc_virus = bcdist(vmat_vegan)
vmat_sp = tibble(lib_id=rownames(vmat_vegan)) %>% 
    left_join(select(lib_meta, lib_id, species, site), by="lib_id")
sp = str_replace(vmat_sp$species, " ", "_")
st = vmat_sp$site
lid = vmat_sp$lib_id


hd = as.dist(host_dist[sp,sp])
hd_COI = as.dist(host_dist_COI[lid,lid])
s_d = as.dist(site_dist[st,st])
td = as.dist(time_dist[lid,lid])

tmp_data = data.frame(
    VirusDist = as.vector(bc_virus), 
    SiteDist = as.vector(s_d), 
    HostDist = as.vector(hd), 
    HostDistCOI = as.vector(hd_COI), 
    TimeDist = as.vector(td)
)

site_m = mantel(VirusDist ~ SiteDist + HostDist + TimeDist, data=tmp_data, mrank=T)
host_m = mantel(VirusDist ~ HostDist + SiteDist + TimeDist, data=tmp_data, mrank=T)
time_m = mantel(VirusDist ~ TimeDist + SiteDist + HostDist, data=tmp_data, mrank=T)
mantel_df = as.data.frame(rbind(host_m, site_m, time_m))
rownames(mantel_df) = c("host", "site", "time")
write_csv(mantel_df, "../output/partial_mantel.csv")

site_m_COI = mantel(VirusDist ~ SiteDist + HostDistCOI + TimeDist, data=tmp_data, mrank=T)
host_m_COI = mantel(VirusDist ~ HostDistCOI + SiteDist + TimeDist, data=tmp_data, mrank=T)
time_m_COI = mantel(VirusDist ~ TimeDist + SiteDist + HostDistCOI, data=tmp_data, mrank=T)
mantel_df_COI = as.data.frame(rbind(host_m_COI, site_m_COI, time_m_COI))
rownames(mantel_df_COI) = c("host", "site", "time")
write_csv(mantel_df_COI, "../output/partial_mantel_COI.csv")

mrm_res = MRM(VirusDist ~ SiteDist + HostDist + TimeDist, data=tmp_data, nperm=10000, mrank=T)
capture.output(print(mrm_res), file="../output/MRM_result.txt")

mrm_res_COI = MRM(VirusDist ~ SiteDist + HostDistCOI + TimeDist, data=tmp_data, nperm=10000, mrank=T)
capture.output(print(mrm_res), file="../output/MRM_result_COI.txt")

host_pm = pmgram(bc_virus, hd, s_d, nperm = 10000, stepsize=20)
site_pm = pmgram(bc_virus, s_d, hd, nperm = 10000, stepsize=50)
write_csv(host_pm$mgram, "../output/host_pmg.csv")
write_csv(site_pm$mgram, "../output/site_pmg.csv")

host_pm_COI = pmgram(bc_virus, hd_COI, s_d, nperm = 10000)
site_pm_COI = pmgram(bc_virus, s_d, hd_COI, nperm = 10000, stepsize=50)
write_csv(host_pm$mgram, "../output/host_pmg_COI.csv")
write_csv(site_pm$mgram, "../output/site_pmg_COI.csv")