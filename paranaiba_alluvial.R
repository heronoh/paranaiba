# carregar libs ----
{
  library("ggplot2")
  library("ggpubr")
  library("dplyr")
  library("ggalluvial")
  library("ggrepel")
}

#
# #colorir pontos ----
# cores <- c( "P1" = "#3900A4FF",
#             "P2" = "#0100BFFF",
#             "P3" = "#008AEFFF",
#             "P4" = "#01C43BFF",
#             "P5" = "#316C00FF",
#             "P6" = "#A4A200FF")
#

# cores <- viridis::turbo(n = 8)[2:6]
cores <- viridis::viridis(n = 7)[2:6]
# cores <- viridis::viridis(n = 5)
#read table ---

# entre no google sheets, na planilha https://docs.google.com/spreadsheets/d/152W3J1of0GEKej2yLQ4qAhcfAHG7zR1wQPgHWUmD50E/edit?usp=sharing
# "Arquivo" -> "Download" -> "Comma separated values (.csv)"

#no Rstudio, vá no canto inferior direito em "Upload"
# salve este arquivo no seu diretório de trabalho e altere o caminho a seguir até apontar pra ele, usando o tab

                                      #AUTO COMPLETE AQUI
comp_analysis_res_tbl <- readxl::read_xlsx(path = "/home/heron/prjcts/paranaiba/results/Paranaíba-Complete_analysis_results-2024-09-05.xlsx",
                                           col_names = TRUE)




# Convert the data to a long format using make_long()

sample_names <- comp_analysis_res_tbl$`Metadata 1` %>% unique()

df_allu <-
  comp_analysis_res_tbl %>%
  dplyr::filter(Type %in% c("Sample")) %>%
  dplyr::mutate("Sample Name" = `Metadata 1`) %>%
  dplyr::mutate("ID clean abd on sample" = 0) %>%
  dplyr::group_by(`Sample Name`,`Curated ID`) %>%
  dplyr::summarize(`ID clean abd on sample` = sum(`ASV absolute abundance`)/mean(`Total clean sample abd.`)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(`Curated ID` %in% c("Pimelodus maculatus",
                                    "Prochilodus lineatus",
                                    "Leporinus octofasciatus",
                                    "Pimelodella sp.3",
                                    "Astyanax lacustris"
                                    )) %>%
  tidyr::pivot_wider(id_cols = "Curated ID",
                     names_from = "Sample Name",
                     values_from = "ID clean abd on sample",
                     values_fill = 0) %>%
  tidyr::pivot_longer(cols = any_of(sample_names),
                      names_to = "Sample Name",
                      values_to = "ID clean abd on sample") %>%
  mutate(`Curated ID` = factor(`Curated ID`,levels = rev(c("Pimelodus maculatus",
                                                       "Prochilodus lineatus",
                                                       "Leporinus octofasciatus",
                                                       "Pimelodella sp.3",
                                                       "Astyanax lacustris"))))


options(scipen = 999)


allu_plot <-
  df_allu %>%
  ggplot(
    aes(y = `ID clean abd on sample`,
        x = `Sample Name`,
        alluvium  = `Curated ID`,
        stratum  = `Curated ID`,
        label  = `Curated ID`)) +
  geom_alluvium(aes(fill = `Curated ID`,
                    col = "#888888"), width = 1/2) +
  geom_stratum(aes(fill = `Curated ID`,
                   col = "#888888"),
               width = 1/6
               # , color = "grey"
  ) +
  # geom_label() +
  # geom_label(stat = "stratum", aes(label = after_stat(stratum)),label.size=0.25, size = 2) +
  # geom_label(stat = "stratum", aes(label = after_stat(stratum)),label.size=0.25, size = 2) +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum),
                                        fill = (`Curated ID`)),
                   col = "#000000",
                   alpha = 0.75,
                   label.size = 0.1,label.padding = 0.1,
                   size = 2) +
  scale_fill_manual(values = cores) +
  scale_colour_manual(values = cores) +
  ggtitle("Top five species on all sites") +
  xlab("Sample Sites") +
  ylab("eDNA RRA on sample") +
  theme_bw() +
  guides(fill = "none",
         col = "none")

allu_plot

ggsave(plot = allu_plot,
       # file = "/home/danielc/projetos/paranaiba/paranaíba--top5_sps--alluvial.pdf",
       file = "/home/heron/prjcts/paranaiba/results/figs/paranaíba--top5_sps--alluvial.pdf",
       device = "pdf",
       units = "cm",
       width = 20,
       height = 12,
       dpi = 300,
       limitsize = FALSE)
