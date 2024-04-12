# Load R libs ----
{
  library(DECIPHER)
  library(dada2)
  library(dplyr)
  # library(phyloseq)
  library(ggplot2)
  library(gridExtra)
  library(ggtree)
  library(Biostrings)
  library(phangorn)
  library(msa)
}

# Load fasta seqs ----

ASVs_seqs <- readDNAStringSet("/home/guiberger/Projetos/Paranaíba/5AlinhamentoPB/PB5.fasta",
                               format = "FASTA") %>%
  unique()

# edit seq names if needed
names(ASVs_seqs) <- names(ASVs_seqs) %>% stringr::str_remove(pattern = "_[:alpha:].*.$")



# orient nucleotides do avoid revcomps ----
ASVs_seqs_ori <- DECIPHER::OrientNucleotides(myXStringSet = ASVs_seqs) %>%
  sort(decreasing = T)

# ASVs_seqs_align <- AlignSeqs(myXStringSet = ASVs_seqs_ori,
#                            iterations = 700,
#                            refinements = 700,
#                            verbose = TRUE)
#
# BrowseSeqs(ASVs_seqs_align,colorPatterns = F,highlight = 0)
# BrowseSeqs(ASVs_seqs_align,colorPatterns = F)

# writeXStringSet(x = ASVs_seqs_align, filepath = "/home/guiberger/Projetos/Paranaíba/5AlinhamentoPB/align_pb5.fasta",
                # format = "FASTA")



# ASVs_seqs_aligniple sequence alignment ----
ASVs_seqs_align <- msa::msa(inputSeqs = ASVs_seqs_ori,
                       method="ClustalW",
                       substitutionMatrix = "clustalw",
                       type="dna",
                       order ="input")


# calculate distance matrix ----
#convert object to phy
ASVs_seqs_align_phy <- as.phyDat(ASVs_seqs_align, type="DNA", names=names(ASVs_seqs_ori))
#calculate dist
ASVs_seqs_dist <- dist.ml(x = ASVs_seqs_align_phy,
                          model ="JC69")

# built tree ----
ASVs_tree <- phangorn::NJ(x = ASVs_seqs_dist) # Note, tip order != sequence order
ASVs_tree


#sequence, names and tips ----
tips_labels <- c(as.character(ASVs_seqs_align))

add.tips(tree = ASVs_tree, tips = tips_labels,where = 10, edge.length = NULL)

fit = pml(ASVs_tree, data=ASVs_seqs_align_phy)

## negative edges length changed to 0!

# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    # rearrangement = "stochastic", control = pml.control(trace = 0))


# write/read tree ----
ape::write.tree(phy = ASVs_tree,
                file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_tree.nwk")

ASVs_tree <- read.tree(file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_tree.nwk")




#add labels ----


#read metadata table ----
metadata_tbl <- read.csv(file =  "/home/heron/prjcts/paranaiba/results/paranaiba-todas_infos_ASVs-2024-04-12.csv",
         header = T,sep = ",",check.names = F)


  tips_metadata <- metadata_tbl %>%
    mutate("tip_label" = stringr::str_remove(pattern = ">",`ASV header`)) %>%
    mutate("BLASTn pseudo-score" = round(`BLASTn pseudo-score`,digits = 2)) %>%
    relocate("tip_label" )

# rownames(tips_metadata) <- tips_metadata$node






# plot tree ----

options(ignore.negative.edge=TRUE)

  tree_plot <- ggtree(tr = ASVs_tree,
                    branch.length = 3,
                    ladderize = T)  %<+%
    tips_metadata +                           #adicionar metadados
  geom_tiplab(align=T,
              linesize=0.5)  +
  geom_treescale(width = 0.4)  +
  theme_tree2() +

    # BLASTn pseudo-score
    geom_tiplab(
    aes(label = `BLASTn pseudo-score`,
        fill = `BLASTn pseudo-score`),
    # offset = 0.04,
    geom = "label",
    # size = 3,
    linetype = "blank" ,
    align = TRUE) +
    scale_fill_gradientn(name = "BLASTn pseudo-score",
                         colours = c("white","red","yellow","green","dark green"),
                         values = c(0.6,1),
                         na.value ="white") +

    # `Curated ID`
    geom_tiplab(
      aes(label = `Curated ID`, col = `Curated ID`),
      # offset = 0.05,
      # size = 3,
      linetype = "blank" ,
      geom = "text",
      align = TRUE) +
    scale_color_manual(values = viridis::turbo(n=10)) +


    # Order
    geom_tiplab(
      aes(label = `Order (BLASTn)`, col = `Order (BLASTn)`),
      # offset = 0.09,
      # size = 3,
      linetype = "blank" ,
      geom = "text",
      align = TRUE) +
    scale_color_manual(values = viridis::turbo(n=10)) +

    # Family
    geom_tiplab(
      aes(label = `Family (BLASTn)`, col = `Family (BLASTn)`),
      # offset = 0.12,
      # size = 3,
      linetype = "blank" ,
      geom = "text",
      align = TRUE) +
    scale_color_manual(values = viridis::turbo(n=105)) +
    guides(col="none")



tree_plot

  # ggsave(file = paste0(figs_path,"/",project_name,"-ASVs_tree.pdf",collapse = ""),
  ggsave(file = paste0("/home/heron/prjcts/paranaiba/results/figs","/",project_name,"-ASVs_tree.pdf",collapse = ""),
         # ggsave(file = paste0(results_path,"/",
         #                      unique(smp_abd_ID_Final$Project),"/",
         #                      unique(smp_abd_ID_Final$Project),"-ASV_length_by_sample-ALL-ASVs.pdf",collapse = ""),
         plot = tree_plot,
         device = "pdf",
         width = 140,
         height = 120,
         limitsize=FALSE,
         units = "cm",
         dpi = 300)





  #referencias ----
  #https://yulab-smu.top/treedata-book/chapter7.html
  #https://www.youtube.com/watch?v=3swFCSt2_x4
  #https://f1000research.com/articles/5-1492/v1
  #https://www.molecularecologist.com/2017/02/08/phylogenetic-trees-in-r-using-ggtree/
  #https://epirhandbook.com/en/phylogenetic-trees-1.html
  #https://boopsboops.blogspot.com/2010/10/negative-branch-lengths-in-neighbour.html
