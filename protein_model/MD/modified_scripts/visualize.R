library(tidyverse)
library(cowplot)

pt_size = 0.4

# FoldX binding energy
foldx_be = read_tsv("FOLDX/FoldX.dat", col_names=c("BE"))

g1 = ggplot(foldx_be, aes(x=1:length(BE)/10, y=BE)) +
    #geom_point(size=pt_size) +
    geom_line() +
    geom_hline(yintercept=mean(foldx_be$BE), linetype = "dashed", color="gray") +
    xlab("Time [ns]") +
    ylab("Binding energy\n[kcal/mol]")

g2 = ggplot(foldx_be, aes(x="", y=BE)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size) + 
    xlab("") +
    ylab("")

p = plot_grid(g1, g2, rel_widths=c(3.5, 1))
ggsave("foldx_be.pdf", plot=p, width=1.718045 * 3, height=1.8)
write_csv(foldx_be, "foldx_be.csv")

# Salt bridge
SB = read_fwf("PLUMED/SALT_BRIDGE.BAR")
colnames(SB) = c("SB_name", "freq")
SB$SB_name = factor(SB$SB_name, levels=c("D30-K417K", "E35-Q493K", "D38-Q493K", "K31-E35", "D38-K353"), ordered=T)
p = ggplot(SB, aes(SB_name, freq)) +
    geom_col() + 
    xlab("Salt bridges") +
    ylab("Frequency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("sb_freq.pdf", plot=p, width=2, height=3)

# RMSD
RMSD_raw = read_lines("PLUMED/COLVAR_RMSD", skip=1)
RMSD_mat = str_split(RMSD_raw, " ", simplify=T)[,2:5]
colnames(RMSD_mat) = c("time", "RMSD_ACE", "RMSD_RBD", "RMSD_INTER")
RMSD = RMSD_mat %>%
    as_tibble() %>%
    mutate_all(as.double)

g1 = ggplot(RMSD, aes(x=time/10, y=RMSD_ACE)) +
    #geom_point(size=pt_size) +
    geom_line() +
    geom_hline(yintercept=mean(RMSD$RMSD_ACE), linetype = "dashed", color="gray") +
    xlab("Time [ns]") +
    ylab("RMSD [Å]") +
    ggtitle("hACE2")

g2 = ggplot(RMSD, aes(x="", y=RMSD_ACE)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size, shape=1, alpha=0.1) + 
    xlab("") +
    ylab("") +
    ggtitle("") 

g3 = ggplot(RMSD, aes(x=time/10, y=RMSD_RBD)) +
    #geom_point(size=pt_size) +
    geom_line() +
    geom_hline(yintercept=mean(RMSD$RMSD_RBD), linetype = "dashed", color="gray") +
    xlab("Time [ns]") +
    ylab("RMSD [Å]") +
    ggtitle("RBD")

g4 = ggplot(RMSD, aes(x="", y=RMSD_RBD)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size, shape=1, alpha=0.1) + 
    xlab("") +
    ylab("") +
    ggtitle("")

g5 = ggplot(RMSD, aes(x=time/10, y=RMSD_INTER)) +
    #geom_point(size=pt_size) +
    geom_line() +
    geom_hline(yintercept=mean(RMSD$RMSD_INTER), linetype = "dashed", color="gray") +
    xlab("Time [ns]") +
    ylab("RMSD [Å]") +
    ggtitle("RBD-hACE2 interface")
g6 = ggplot(RMSD, aes(x="", y=RMSD_INTER)) +
    geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) +
    #geom_jitter(height = 0, width = 0.2, size=pt_size, shape=1, alpha=0.1) + 
    xlab("") +
    ylab("") +
    ggtitle("")

p=plot_grid(g1,g2,g3,g4,g5,g6, nrow=3, ncol=2, rel_widths = c(3.5, 1))
ggsave("backbone_rmsd.pdf", plot=p, width=1.718045 * 3, height=6)
write_csv(RMSD, "RMSD.csv")