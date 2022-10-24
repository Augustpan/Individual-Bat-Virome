library(tidyverse)
library(cowplot)

be_files = c(
    "SARS-CoV-2.1-ANALYSIS/foldx_be.csv",
    "SARS-CoV-2.2-ANALYSIS/foldx_be.csv",
    "SARS-CoV-2.3-ANALYSIS/foldx_be.csv",
    "Mod1-ANALYSIS/foldx_be.csv",
    "Mod2-ANALYSIS/foldx_be.csv",
    "Mod3-ANALYSIS/foldx_be.csv"
)

rmsd_files = c(
    "SARS-CoV-2.1-ANALYSIS/RMSD.csv",
    "SARS-CoV-2.2-ANALYSIS/RMSD.csv",
    "SARS-CoV-2.3-ANALYSIS/RMSD.csv",
    "Mod1-ANALYSIS/RMSD.csv",
    "Mod2-ANALYSIS/RMSD.csv",
    "Mod3-ANALYSIS/RMSD.csv"
)

be_all = NULL
for (filename in be_files) {
    mod = str_split(filename, "-ANALYSIS/", simplify = T)[1]
    be = read_tsv(filename) %>%
        mutate(model=mod, frame=0:(n()-1), time=frame/10, .before=1)
    be_all = rbind(be_all, be)
}

rmsd_all = NULL
for (filename in rmsd_files) {
    mod = str_split(filename, "-ANALYSIS/", simplify = T)[1]
    rmsd = read_csv(filename) %>%
        mutate(time = time/10) %>%
        mutate(model=mod, frame=time*10, .before=1)
    rmsd_all = rbind(rmsd_all, rmsd)
}

linewidth = 0.3
g1 = ggplot(be_all, aes(x=time, y=BE, color=model)) +
    geom_line(size=linewidth) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    xlab("Time [ns]") +
    ylab("Binding energy\n[kcal/mol]") +
    theme_bw() +
    theme(legend.position="none")

g2 = ggplot(be_all, aes(x=model, y=BE, fill=model)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), size=linewidth) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size) + 
    xlab("") +
    ylab("") + 
    theme_bw() + theme(legend.position="none")

p1 = plot_grid(g1, g2, rel_widths=c(2, 1))

g1 = ggplot(rmsd_all, aes(x=time, y=RMSD_INTER, color=model)) +
    geom_line(size=linewidth) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    xlab("Time [ns]") +
    ylab("RMSD [Å]") +
    theme_bw() + theme(legend.position="none")
g2 = ggplot(rmsd_all, aes(x=model, y=RMSD_INTER, fill=model)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), size=linewidth) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size, shape=1, alpha=0.1) + 
    xlab("") +
    ylab("") +
    theme_bw() + theme(legend.position="none")

g3 = ggplot(rmsd_all, aes(x=time, y=RMSD_RBD, color=model)) +
    geom_line(size=linewidth) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    xlab("Time [ns]") +
    ylab("RMSD [Å]") +
    theme_bw() + theme(legend.position="none")
g4 = ggplot(rmsd_all, aes(x=model, y=RMSD_RBD, fill=model)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), size=linewidth) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size, shape=1, alpha=0.1) + 
    xlab("") +
    ylab("") +
    theme_bw() + theme(legend.position="none")

p2=plot_grid(g1,g2,g3,g4, nrow=2, ncol=2, rel_widths = c(2, 1))

ggsave("be.pdf", plot=p1, width=7, height=2)
ggsave("rmsd.pdf", plot=p2, width=7, height=4)