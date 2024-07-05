# carregar libs ----
{
library("ggplot2")
library("ggpubr")
library("dplyr")
}


#colorir pontos ----
cores <- c( "P1" = "#3900A4FF",
            "P2" = "#0100BFFF",
            "P3" = "#008AEFFF",
            "P4" = "#01C43BFF",
            "P5" = "#316C00FF",
            "P6" = "#A4A200FF")



# ler tabela ou criar tabela de entrada ----
div_elev_tbl <- readRDS(file = "/home/heron/prjcts/paranaiba/tabela_alphaDiv_elevation.RDS")


# extraindo tabela por ponto amostral par ao manuscrito
elevation_plot <- div_elev_tbl %>%
  ggplot(aes(x = Elevation,
             y = Values
             )) +
  geom_smooth(method = "lm",
              se = T) +
  stat_cor(method = "pearson",
           label.x = 600,
           label.y = 30
  ) +
  geom_point(aes(
    group = agrupador,
    col = agrupador),
    size = 4) +
  scale_colour_manual(values = cores,
                      breaks = names(cores),
                      name = "Sampling site") +
  scale_x_continuous(breaks = c(seq(500, 800, 50)),transform = "reverse") +
  scale_y_continuous(breaks = c(seq(0,35,5)),
                     name = "alpha Diversity (Num. species)") +
  theme_bw(base_size = 14)



# salvar gráfico ----
ggsave(file = "/home/heron/prjcts/paranaiba/results/figs/corr_elevation_div.pdf",
       plot = elevation_plot,
       device = "pdf",
       width = 16,
       height = 10,
       units = "cm",
       limitsize = FALSE,
       dpi = 300)
