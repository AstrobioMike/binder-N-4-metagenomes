library(KEGGREST)
library(tidyverse)

N_KO_norm_cov_tab <- read.table("data/N_KO_norm_cov.tsv", sep="\t", header=TRUE)

N_KO_info_tab <- read.table("data/All_N_KO_info.tsv", sep = "\t", header = TRUE)

N_gene_level_tab <- read.table("data/N-gene-normalized-master-tab.tsv", sep = "\t", header = TRUE, na.strings = "")

get_KO_info <- function( target_KO, KO_KEGG_tab = N_KO_norm_cov_tab ) {

    # bulding link
    curr_KO_link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", target_KO)

    # getting info on KEGG term
    try(curr_info <- keggGet(target_KO), silent = TRUE)

    # checking if it was found at KEGG
    if ( ! exists("curr_info") ) {

        cat(paste0("\n  It seems '", target_KO, "' wasn't found at KEGG. It may have been removed, or may have never existed.\n\n  This should be the link if it were there if you wanna take a look:\n\n      ", curr_KO_link, "\n\n"))
        # couldn't figure out a better way to just report this and not give an error while exiting
        stop_no_error <- function() {
            opt <- options(show.error.messages = FALSE)
            on.exit(options(opt))
            stop()
        }

        stop_no_error()
    }

    # parsing some info
    if ( length(curr_info[[1]]$ENTRY) > 0 ) { curr_KO_ID <- curr_info[[1]]$ENTRY %>% as.vector } else { curr_KO_ID <- NA }
    if ( length(curr_info[[1]]$NAME) > 0 ) { curr_KO_name <- curr_info[[1]]$NAME %>% as.vector } else { curr_KO_name <- NA }
    if ( length(curr_info[[1]]$DEFINITION) > 0 ) { curr_KO_def <- curr_info[[1]]$DEFINITION %>% as.vector } else { curr_KO_def <- NA }

    # seeing if it's in our table
    if ( target_KO %in% (KO_KEGG_tab %>% pull(KO_ID)) ) {

      curr_in_data <- "Yes"

    } else {

      curr_in_data <- "No"

    }

    # reporting info
    cat("\n  KO ID          :  ", curr_KO_ID, "\n")
    cat("  KO name        :  ", curr_KO_name, "\n")
    cat("  KO definition  :  ", curr_KO_def, "\n")
    cat("  KO link        :  ", curr_KO_link, "\n")
    cat("  In our N data? :  ", curr_in_data, "\n\n")

}

plot_KO <- function( target, tab = N_KO_norm_cov_tab ) {

    # subsetting to target KO
    sub_tab <- tab %>% filter(label == target) %>% select(c(1:5, 10, 11, 12))

    # making long formatted version
    sub_long <- sub_tab %>% select(1:5) %>% pivot_longer(-KO_ID, names_to = "Depth", values_to = "norm_cov") %>% data.frame(check.names = FALSE)

    # making labels
    num_uniq_genes <- sub_tab$num_uniq_genes
    ko_id <- sub_tab$KO_ID
    ko_name <- sub_tab$KO_name
    ko_def <- sub_tab$KO_def
    main_lab <- paste0(ko_id, " (", num_uniq_genes, " unique gene copies)")
    sub_lab <- paste0(ko_name, " - ", ko_def)

    # making plot
    plot <- ggplot() + geom_bar(data = sub_long, aes(y = norm_cov, x = Depth, fill = Depth), color = "black", stat = "identity", width = 0.75) +
        theme_bw() + labs(title = main_lab, subtitle = sub_lab, y = "Coverage per Million (CPM)") +
        theme(axis.title.y = element_text(face = "bold", size = 15), axis.text.y = element_text(size=14)) +
        theme(axis.title.x = element_text(face = "bold", size = 15), axis.text.x = element_text(size=15)) +
        theme(legend.position = "none") + theme(plot.title = element_text(size = 16, face = "bold")) +
        theme(plot.subtitle = element_text(size = 14))

    return(plot)
}

plot_all_N_KOs <- function( tab = N_KO_norm_cov_tab ) {

    # adding unique number of genes to KO ID
    ids <- paste0(tab$KO_ID, " (", tab$num_uniq_genes, ")")
    tab$KO_ID <- ids

    # converting to long format
    long_tab <- tab %>% select(1:5) %>% pivot_longer(-KO_ID, names_to = "Depth", values_to = "norm_cov")

    # plotting

    plot <- ggplot() +
      geom_bar(data=long_tab, aes(y=norm_cov, x=Depth, fill=Depth), stat="identity", width=0.75) +
      facet_wrap(. ~ KO_ID, scales="free_y", ncol = 10) + theme_bw() + labs(title = "All Nitrogen-related KO coverages", y = "Coverage per Million (CPM)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank()) +
      theme(legend.position="bottom", legend.title=element_text(face="bold")) + labs(fill = "Depth") +
      theme(strip.text = element_text(size = 10, face="bold"))

    return(plot)
}

get_seq <- function( target_gene_ID, type = "AA", seqs_tab = N_gene_level_tab ) {

    header <- paste0(">", target_gene_ID)

    if ( type == "AA" ) {
        seq <- seqs_tab %>% filter(gene_ID == target_gene_ID) %>% pull(aa_sequence)
    } else {
        seq <- seqs_tab %>% filter(gene_ID == target_gene_ID) %>% pull(dna_sequence)
    }

    paste(header, seq, sep = "\n")
}

plot_KO_tax <- function(target_KO_label, tab = N_gene_level_tab, info_tab = N_KO_norm_cov_tab) {

    # getting the KO ID
    target_KO <- info_tab %>% filter(label == target_KO_label) %>% pull(KO_ID)

    # subset target KO
    sub_tab <- tab %>% filter(KO_ID == target_KO) %>% select(-c("dna_sequence", "aa_sequence"))

    # splitting bacteria, archaea, and completely unclassified
    bac_sub <- sub_tab %>% filter(domain == "Bacteria")
    arc_sub <- sub_tab %>% filter(domain == "Archaea")
    completely_unclassified <- sub_tab %>% filter(domain == "NA")

    # grouping and summing
    bac_grouped <- bac_sub %>% group_by(phylum) %>% summarize(S1=sum(S1), S2=sum(S2), S3=sum(S3), S4=sum(S4)) %>% data.frame()
    arc_grouped <- arc_sub %>% group_by(phylum) %>% summarize(S1=sum(S1), S2=sum(S2), S3=sum(S3), S4=sum(S4)) %>% data.frame()
    completely_unclassified_grouped <- completely_unclassified %>% group_by(phylum) %>% summarize(S1=sum(S1), S2=sum(S2), S3=sum(S3), S4=sum(S4)) %>% data.frame()


    # converting from factors to vectors
    bac_grouped$phylum <- as.vector(bac_grouped$phylum)
    arc_grouped$phylum <- as.vector(arc_grouped$phylum)
    completely_unclassified_grouped$phylum <- as.vector(completely_unclassified_grouped$phylum)

    # setting NAs to unclassifieds
    bac_grouped$phylum[bac_grouped$phylum == "NA"] <- "Unclassified Bacteria"
    arc_grouped$phylum[arc_grouped$phylum == "NA"] <- "Unclassified Archaea"
    completely_unclassified_grouped$phylum[completely_unclassified_grouped$phylum == "NA"] <- "Unclassified"

    # adding "(A)" in front of archaea phyla just to make them easier to spot
    if ( nrow(arc_grouped) > 0 ) {
        arc_grouped$phylum <- paste0("(A) ", arc_grouped$phylum)
    }

    # combining all together
    combined_tab <- rbind(bac_grouped, arc_grouped, completely_unclassified_grouped)

    # filtering to only those contributing at least 1% of any depth's minimum total (and adding an "other" to track the remaining)
    num_tab <- combined_tab %>% column_to_rownames(var = "phylum")
    min_vals <- colSums(num_tab) * 0.01

    keepers <- vector()

    for ( col in 1:4 ) {

        for ( row in 1:nrow(num_tab) ) {

            curr_tax <- row.names(num_tab)[row]

            if ( num_tab[row, col] > min_vals[col] ) {

                if (! curr_tax %in% keepers ) {

                    keepers <- c(keepers, curr_tax)

                }

            }

        }

    }

    # parsing those we want and dont want
    keepers_tab <- combined_tab %>% filter(phylum %in% keepers)
    others_tab <- combined_tab %>% filter(!phylum %in% keepers)

    # combining all the "others" into one row
    others_row <- data.frame("phylum" = "Other*", "S1" = sum(others_tab$S1), "S2" = sum(others_tab$S2), "S3" = sum(others_tab$S3), "S4" = sum(others_tab$S4))

    # adding others to make final table
    final_tab <- rbind(keepers_tab, others_row)

    # getting factor order for plotting, want bacteria, unclassified bacteria, archaea, unclassified archaea, other
    all_IDs <- final_tab$phylum

    # getting all bacterial phyla
    classified_bac_IDs <- all_IDs[!grepl("Unclassified", all_IDs) & !grepl("(A)", all_IDs) & all_IDs != "Other*"]
    # getting all archaeal phyla
    classified_arc_IDs <- all_IDs[grepl("(A)", all_IDs)]

    # generating factor order (should handle things fine when there is no "Unclassified Archaea" for example)
    factor_order <- c(classified_bac_IDs, "Unclassified Bacteria", classified_arc_IDs, "Unclassified Archaea", "Unclassified", "Other*")

    # setting factor order for plotting
    final_tab$phylum <- factor(final_tab$phylum, levels = factor_order)

    # converting to long format
    sub_long <- final_tab %>% pivot_longer(-phylum, names_to = "Depth", values_to = "CPM")

    # making titles and such
    num_unique_genes <- info_tab %>% filter(KO_ID == target_KO) %>% pull(num_uniq_genes)
    ko_name <- info_tab %>% filter(KO_ID == target_KO) %>% pull(KO_name)
    ko_def <- info_tab %>% filter(KO_ID == target_KO) %>% pull(KO_def)

    main_lab <- paste0("Phyla distribution for ", target_KO, " (", num_unique_genes, " unique genes)")
    sub_lab <- paste0(ko_name, " - ", ko_def)

    ggplot() + geom_bar(data = sub_long, aes(y = CPM, x = Depth, fill = phylum), color = "black", position = "dodge", stat = "identity", width = 0.75) +
        theme_bw() + labs(title = main_lab, subtitle = sub_lab, y = "Coverage per Million (CPM)") +
        theme(axis.title.y = element_text(face = "bold", size = 15), axis.text.y = element_text(size=14)) +
        theme(axis.title.x = element_text(face = "bold", size = 15), axis.text.x = element_text(size=15)) +
        theme(plot.title = element_text(size = 16, face = "bold")) +
        theme(plot.subtitle = element_text(size = 14)) + labs(fill = "Phylum") +
        theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))
}

subset_tax_tab <- function(target_KO_label, tab = N_gene_level_tab, info_tab = N_KO_norm_cov_tab) {

    # getting the KO ID
    target_KO <- info_tab %>% filter(label == target_KO_label) %>% pull(KO_ID)

    # subsetting table
    tab %>% filter(KO_ID == target_KO) %>% select(1:16)
}

# example sending to shinyapps.io
# rsconnect::deployApp('N-shiny/')
