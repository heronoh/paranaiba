# Load R libs ----
{
  library(DECIPHER)
  library(dada2)
  library(dplyr)
  library(phytools)
  library(ggplot2)
  library(gridExtra)
  library(ggtree)
  library(Biostrings)
  library(phangorn)
  library(msa)
  library(ape)
  library(mergeTrees)
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
                       method = "ClustalW",
                       substitutionMatrix = "clustalw",
                       type = "dna",
                       order = "input")



# MaximumLlikelyhood tree ----

# https://klausvigo.github.io/phangorn/articles/MLbyHand.html

library(ape)
library(phangorn)
# fdir <- system.file("extdata/trees", package = "phangorn")
# primates <- read.phyDat(file.path(fdir, "primates.dna"),
#                         format = "interleaved")

#convert object to phy
# ASVs_seqs_align_phy
primates <- as.phyDat(ASVs_seqs_align,
                                 type = "DNA",
                                 names = names(ASVs_seqs_align))



dm <- dist.ml(primates)
treeNJ  <- NJ(dm)
fit <- pml(treeNJ, data=primates)
fit
methods(class="pml")
fitJC  <- optim.pml(fit, rearrangement="NNI")
logLik(fitJC)
fitF81 <- update(fitJC, k=4, inv=0.2, bf=baseFreq(primates))
fitF81
fitGTR <- optim.pml(fitF81, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

anova(fitJC, fitGTR)

SH.test(fitGTR, fitJC)
AIC(fitJC)
AIC(fitGTR)

bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE,
                    control = pml.control(trace = 0))

plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

cnet <- consensusNet(bs, p=0.2)
plot(cnet, show.edge.label=TRUE)


ASVs_tree <-  fitJC$tree

fitJC$


                  # calculate distance matrix ----

                  #calculate dist
                  ASVs_seqs_dist <- dist.ml(x = ASVs_seqs_align_phy,
                                            model = "JC69")

                  # built tree ----
                  ASVs_tree <- phangorn::NJ(x = ASVs_seqs_dist) # Note, tip order != sequence order
                  ASVs_tree


                  #sequence, names and tips ----
                  tips_labels <- c(as.character(ASVs_seqs_align))

                  add.tips(tree = ASVs_tree, tips = tips_labels,where = 10, edge.length = NULL)

                  fit = pml(ASVs_tree, data = ASVs_seqs_align_phy)

                  ## negative edges length changed to 0!

                  # fitGTR <- update(fit, k=4, inv=0.2)
                  # fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                      # rearrangement = "stochastic", control = pml.control(trace = 0))


# write/read tree ----
# ape::write.tree(phy = bs,
                # file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_bootstraps.nwk")

                  # write/read tree ----
ape::write.tree(phy = ASVs_tree,
                file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_tree.nwk")

# ASVs_tree_boot <- read.tree(file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_bootstraps.nwk")
ASVs_tree <- read.tree(file = "/home/heron/prjcts/paranaiba/results/paranaiba_ASVs_tree.nwk")


data(woodmouse)


ASVs_seqs_align

ASVs_seqs_align_bin <- ape::as.DNAbin(ASVs_seqs_align)

f <- function(x) nj(dist.dna(x))


#bootstra ----
bs<-boot.phylo(phy = ASVs_tree, ASVs_seqs_align_bin, f, quiet = TRUE)


ASVs_tree$node.label <- bs

ASVs_tree

bs_tibble <- tibble(
  # hard-code node ID: internal nodes start after tip nodes,
  # and phy$node.label is in the same order as internal nodes
  "node" = 1:Nnode(ASVs_tree) + Ntip(ASVs_tree),
  # Don't print BS < 50%
  # (or to show all BS, just do `bootstrap = phy$node.label`)
  "bootstrap" = ifelse(ASVs_tree < 50, "", ASVs_tree$node.label))

#read metadata table ----
metadata_tbl <- read.csv(file =  "/home/heron/prjcts/paranaiba/results/paranaiba-todas_infos_ASVs-2024-04-12.csv",
         header = T,sep = ",",check.names = F)


  tips_metadata <- metadata_tbl %>%
    mutate("tip_label" = stringr::str_remove(pattern = ">",`ASV header`)) %>%
    mutate("BLASTn pseudo-score" = round(`BLASTn pseudo-score`,digits = 2)) %>%
    relocate("tip_label" )

# rownames(tips_metadata) <- tips_metadata$node


  ASVs_tree$edge
  ASVs_tree$edge.length
  ASVs_tree$tip.label
  ASVs_tree$Nnode


# rotate nodes to match taxa? ----
  # http://blog.phytools.org/2015/04/finding-closest-set-of-node-rotations.html

  #nao funcionou
  # tree<-pbtree(n=26,tip.label=LETTERS)
  # plotTree(tree)


  ## random set of 100 rotations
  ## tree all scrambled up
  ## the objective function going to zero indicated fully
## unscrambled

  # plotTree(tree)

  # ## the objective function going to zero indicated fully
  # ## unscrambled
  # x<-setNames(1:Ntip(tree),LETTERS)
  # unscrambled<-minRotate(tree,x)
  #
  #
  #
  #


  # ordered_tip_names <- tips_metadata %>%
  #   arrange(`Order (BLASTn)`) %>%
  #   pull(tip_label)
  #
  #
  # ASVs_tree$tip.label
  #
  #
  # tip_order <- setNames(match(ASVs_tree$tip.label,ordered_tip_names),ordered_tip_names)
  #
  #
  # ASVs_tree <- phytools::minRotate(tree = ASVs_tree,
  #                                  x = tip_order )
  #
  #
  #
  # ASVs_tree <- phytools::untangle(ASVs_tree,method="reorder")
  #
  #
# plot tree ----

options(ignore.negative.edge = TRUE)

#set anotations position relative to tip
tip_alignment <- TRUE



  tree_plot <- ggtree(tr = ASVs_tree,
                    branch.length = 3,
                    ladderize = T)  %<+%
    tips_metadata +                           #adicionar metadados
  geom_tiplab(align = tip_alignment,
              linesize = 0.5)  +
  geom_treescale(width = 0.4)  +
  theme_tree2() +

    # `Curated ID`              #####################
  geom_tiplab(
    aes(label = `Curated ID`, col = `Curated ID`),
    offset = 0.027,
    linetype = "blank" ,
    geom = "text",
    align = tip_alignment) +
    scale_color_manual(values = viridis::turbo(n = 10)) +

    # BLASTn pseudo-score      #####################
  geom_tiplab(
    aes(label = `BLASTn pseudo-score`,
        fill = `BLASTn pseudo-score`),
    offset = 0.053,
    geom = "label",
    size = 2,
    linetype = "blank" ,
    align = tip_alignment) +
    scale_fill_gradientn(name = "BLASTn pseudo-score",
                         colours = c("white","red","yellow","green","dark green"),
    values = c(0.6,1),
    na.value = "white") +

    # Order      #####################
  geom_tiplab(
    aes(label = `Order (BLASTn)`, col = `Order (BLASTn)`),
    offset = 0.12,
    linetype = "blank" ,
    geom = "text",
    align = tip_alignment) +
    scale_color_manual(values = viridis::turbo(n = 10)) +
    # Family      #####################
  geom_tiplab(
    aes(label = `Family (BLASTn)`, col = `Family (BLASTn)`),
    offset = 0.09,
    linetype = "blank" ,
    geom = "text",
    align = tip_alignment) +
    scale_color_manual(values = viridis::turbo(n = 105)) +
    guides(col = "none") +
    geom_nodelab(aes(label = nodesupport))
    %<+% bs_tibble +
    # now we can show BS values using geom_text()
    geom_text(aes(label=bootstrap), hjust=-.25, size = 3)



tree_plot


# sava tree plot ----

  # ggsave(file = paste0("/home/heron/prjcts/paranaiba/results/figs","/",project_name,"-ASVs_tree.pdf",collapse = ""),
  ggsave(file = paste0("/home/heron/prjcts/paranaiba/results/figs/Paranaiba-ASVs_tree.pdf",
                       collapse = ""),
         plot = tree_plot,
         device = "pdf",
         width = 150,
         height = 120,
         limitsize = FALSE,
         units = "cm",
         dpi = 300)



# referencias ----
  #https://yulab-smu.top/treedata-book/chapter7.html
  #https://www.youtube.com/watch?v=3swFCSt2_x4
  #https://f1000research.com/articles/5-1492/v1
  #https://www.molecularecologist.com/2017/02/08/phylogenetic-trees-in-r-using-ggtree/
  #https://epirhandbook.com/en/phylogenetic-trees-1.html
  #https://boopsboops.blogspot.com/2010/10/negative-branch-lengths-in-neighbour.html
#https://cran.r-project.org/web/packages/phangorn/vignettes/IntertwiningTreesAndNetworks.html
