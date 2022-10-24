library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

vtable = read_csv("../generated_data/virus_table_filtered.csv") %>%
    select_if(function(x) { ifelse(is.character(x), T, sum(x)>0)}) # drop empty columns

virus_name_arr = colnames(select(vtable, -lib_id))
lib_id_arr = vtable$lib_id

vmeta = read_csv("../raw_data/virus_meta.csv") %>%
    filter(virus_name %in% virus_name_arr)
lib_meta = read_csv("../raw_data/lib_metadata.csv")

vmat = t(select(vtable, all_of(virus_name_arr)))
colnames(vmat) = vtable$lib_id

# reordering
vmeta_ord = arrange(vmeta, order)
meta_ord = arrange(lib_meta, order)
vmat_ord = log10(vmat[vmeta_ord$virus_name, meta_ord$lib_id]+1)

# set color mapping
genus = sort(c("Rhinolophus", "Aselliscus", "Hipposideros", "Rousettus", "Cynopterus", "Eonycteris"))
genus_color = brewer.pal(6, "Set1")
site = c("LS","BS","CX","ZK","WD","ML")
site_color = brewer.pal(6, "Set2")
family = unique(vmeta_ord$virus_family)
family_color = brewer.pal(length(family), "Set3")
species = sort(unique(meta_ord$species))
species_color = c(brewer.pal(10, "Paired"), brewer.pal(5, "Dark2"))
RPM_col_fun = colorRamp2(c(0, ceiling(max(vmat_ord))), c("#F4F2ED", "#DD1332"))

# heatmap annotations
col_ano = columnAnnotation(
    host_genus = anno_simple(meta_ord$genus, col=setNames(genus_color, genus), height=unit(0.3, "cm")),
    host_species = anno_simple(meta_ord$species, col=setNames(species_color, species), height=unit(0.3, "cm")),
    sample_site = anno_simple(meta_ord$site, col=setNames(site_color, site), height=unit(0.3, "cm")),
    gap = unit(1, "mm")
)
row_ano = rowAnnotation(
    virus_family = anno_simple(vmeta_ord$virus_family, col=setNames(family_color, family), width=unit(0.3, "cm"))
)

# draw heatmap
pdf(file="../output/virus_heatmap.pdf", width=12, height=8)
ht = Heatmap(vmat_ord, 
             cluster_columns = F, 
             cluster_rows = F, 
             row_order = 1:nrow(vmat_ord),
             column_order = 1:ncol(vmat_ord),
             column_split = factor(meta_ord$genus, levels=c("Aselliscus", "Hipposideros", "Rhinolophus", "Rousettus", "Eonycteris", "Cynopterus"), ordered=T),
             row_split = factor(vmeta_ord$virus_family, levels=c("Coronaviridae","Paramyxoviridae","Astroviridae","Caliciviridae","Picornaviridae","Reoviridae","Cicoviridae","Parvoviridae","Anelloviridae","Adenoviridae","Polyomaviridae"), ordered = T),
             
             col = RPM_col_fun, 
             top_annotation = col_ano,
             left_annotation = row_ano,
             column_labels = rep("", ncol(vmat_ord)))

l1 = Legend(labels = genus, title="Host genus", legend_gp=gpar(fill=genus_color))
l2 = Legend(labels = species, title="Host species", legend_gp=gpar(fill=species_color))
l3 = Legend(labels = site, title="Sample site", legend_gp=gpar(fill=site_color))
l4 = Legend(labels = family, title="Virus family", legend_gp=gpar(fill=family_color))
draw(ht, annotation_legend_list = list(l1,l2,l3,l4))
dev.off()