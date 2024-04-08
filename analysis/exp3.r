
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)

########################################################
# Plot showing the time comparison
########################################################

time_df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/time.csv")

titles <- c("aquatic" = "Aquatic Reads (n=200k)", 
            "human_gut" = "Human Gut Reads (n=200k)", 
            "soil" = "Soil Reads (n=200k)")
plot <- ggplot(time_df, aes(x=region, y=time, group=interaction(method, digestion), color=method)) +
        geom_line(aes(linetype=digestion)) +
        geom_point() +
        facet_wrap(~dataset, labeller=as_labeller(titles)) +
        theme_bw() +
        theme(strip.text = element_text(face = "bold"),
              legend.position="right") +
        labs(x="Region", y="Time (s)") +
        scale_y_continuous(breaks=seq(0,90,10)) +
        scale_color_discrete(name="Method", labels=c("Cliffy", "Kraken2")) +
        scale_linetype_manual(name="Digestion Approach",
                              values = c("no_digestion"="longdash", 
                                         "dna_minimizers" = "dashed",
                                         "minimizers"="dotted",
                                         "na"="solid"),
                              labels=c("DNA minimizers", "Minimizers",  "No digestion"),
                              breaks=c("dna_minimizers", "minimizers", "no_digestion")) 
                
plot
ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/time.pdf", 
       plot, 
       width=8, 
       height=3,
       dpi=500)

########################################################
# Plot showing the accuracy comparison
########################################################

acc_df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/accuracy.csv")

# Define the plotting function
plot_acc_func <- function(curr_dataset) {
        # subset the dataframe to relavant dataset
        df_temp <- acc_df[acc_df$dataset == curr_dataset, ]
        
        # subset to only relevant kraken approaches
        method_names <- c("cliffy", "kraken2_with_vp_m31", "kraken2_without_vp_m31", "kraken2_with_vp_m27", "kraken2_without_vp_m27", "kraken2_with_vp_m23", "kraken2_without_vp_m23")
        df_subset <- df_temp %>% filter((method %in% method_names))
        df_subset$method <- factor(df_subset$method, levels = c("cliffy", "kraken2_with_vp_m31", "kraken2_with_vp_m27", "kraken2_with_vp_m23", "kraken2_without_vp_m31", "kraken2_without_vp_m27", "kraken2_without_vp_m23"))
        method_new_names <- c("Cliffy (TP/n)", "Kraken2, k=31 (TP+VP/n)", "Kraken2, k=27 (TP+VP/n)", "Kraken2, k=23 (TP+VP/n)", "Kraken2, k=31 (TP/n)", "Kraken2, k=27 (TP/n)", "Kraken2, k=23 (TP/n)")
        
        # set the order of the clades since that will be x-axis
        levels_order <- c("genus", "family", "order", "class", "phylum", "domain")
        df_subset$clade <- factor(df_subset$clade, levels=levels_order, labels=c("Genus", "Family", "Order", "Class", "Phylum", "Domain"))
        
        # title of facets
        titles <- c("V1_V2" = "V1-V2 Reads (n=200k)",
                    "V3_V4" = "V3-V4 Reads (n=200k)",
                    "V4_V4" = "V4 Reads (n=200k)",
                    "V4_V5" = "V4-V5 Reads (n=200k)")
        
        # make the plot
        plot <- ggplot(df_subset, aes(x=clade, y=sensitivity, group=interaction(method, digestion, query_approach))) +
                geom_line(aes(color=method, linetype=query_approach)) +
                geom_point(aes(shape=digestion, color=method, fill=method)) +
                theme_bw() +
                theme(axis.text.x=element_text(angle=45, vjust = 1., hjust=1,size=12,color="black"),
                      axis.title.x=element_text(size=14),
                      axis.title.y=element_text(size=14),
                      axis.text.y=element_text(size=12, color="black"),
                      strip.text=element_text(face="bold", size=12),
                      legend.title=element_text(size=12),
                      legend.text=element_text(size=10)) +
                scale_y_continuous(breaks=seq(0.70, 1.0, 0.05)) +
                labs(x="Major Clade", y="Metric") +
                facet_wrap(~region, labeller=as_labeller(titles)) +
                scale_color_discrete(name="Method (Metric)", labels=method_new_names) +
                scale_fill_discrete(name="Method (Metric)", labels=method_new_names) +
                scale_linetype_discrete(name="Query Approach",
                                      labels=c("LCA", "Approximate Doc. Listing")) +
                scale_shape_manual(name="Digestion Approach",
                                      values = c("no_digestion"=24, 
                                                 "dna_minimizers" = 22,
                                                 "minimizers"=23,
                                                 "na"=21),
                                      labels=c("DNA minimizers", "Minimizers",  "No digestion"),
                                      breaks=c("dna_minimizers", "minimizers", "no_digestion")) 
        
        plot    
        
        # save plot with dataset name
        ggsave(paste0("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/accuracy_", curr_dataset, ".pdf"), 
               plot, 
               width=8, 
               height=6,
               dpi=500)
}

# make a plot for each dataset
lapply(unique(acc_df$dataset), plot_acc_func)


##########################################################
# Plot showing the # of increases in each direction
##########################################################

inc_df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/minimizers_increases.csv", 
               header=FALSE)
colnames(inc_df) <- c("pos", "incfromleft", "incfromright")
inc_df['total'] <- inc_df['incfromleft'] + inc_df['incfromright']

df_filt <- inc_df %>% filter(pos<25)
df_long <- df_filt %>% pivot_longer(cols = c(incfromleft, incfromright), names_to = "variable", values_to = "value")

ggplot(df_long, aes(x = factor(pos), y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = "Position", y = "Value") +
        theme_minimal()

plot1 <- ggplot(df_long, aes(x=factor(pos),y=value, fill=variable)) +
         theme_bw() +
         geom_bar(stat="identity", position=position_dodge(), color="black") +
         scale_x_discrete(breaks=seq(0,25,1)) +
         scale_y_continuous(breaks=seq(0,5000000,500000), labels = scales::scientific) +
         labs(x="Number of Monotonic Increases", y="Number of Profiles") +
         theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
              axis.title.x=element_text(size=12),
              axis.title.y=element_text(size=12),
              axis.text=element_text(size=12, color="black")) +
         scale_fill_discrete(name="Number of Increases",
                             labels=c("Left to Right Increases", "Right to Left Increases")) 
plot1

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/minimizers_inc_plot.pdf", 
       plot=plot1, 
       dpi=800, 
       device="pdf", 
       width=8, 
       height=4)

########################################################
# Plot showing the abundances for each dataset
########################################################
abund_df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/abundance.csv")

# Define the plotting function
plot_abund_func <- function(curr_dataset) {
        # subset the dataframe to relevant dataset
        df_subset <- abund_df[abund_df$dataset == "human_gut", ]
        
        # set the order of the methods since that will be x-axis
        levels_order <- c("truth", "docprofiles", "kraken2")
        df_subset$method <- factor(df_subset$method, levels=levels_order, labels=c("Truth", "Cliffy", "Kraken2"))
        df_subset$genus <- factor(df_subset$genus, levels = c(setdiff(unique(df_subset$genus), "Others"), "Others"))

        # title of facets
        titles <- c("V1_V2" = "V1-V2 Reads (n=200k)",
                    "V3_V4" = "V3-V4 Reads (n=200k)",
                    "V4_V4" = "V4 Reads (n=200k)",
                    "V4_V5" = "V4-V5 Reads (n=200k)")
        
        # define color mapping
        color_mapping <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(df_subset$genus)))
        names(color_mapping) <- sort(unique(df_subset$genus))
        color_mapping["Others"] <- "gray"
        
        # make the plot
        plot <- ggplot(df_subset, aes(x=method, y=abundance, fill=genus)) +
                geom_bar(stat="identity", colour="black", size=0.1) +
                scale_fill_manual(values = color_mapping) +
                theme_bw() +
                facet_wrap(~region, labeller=as_labeller(titles)) +
                theme(axis.text.x=element_text(angle=0,size=12,color="black"),
                      axis.title.x=element_text(size=14),
                      axis.title.y=element_text(size=14),
                      axis.text.y=element_text(size=12, color="black"),
                      strip.text=element_text(face="bold", size=12),
                      legend.title=element_text(size=12),
                      legend.text=element_text(size=10),
                      legend.position="none") +
                labs(x="Method", y="Abundance") 
        plot    
        
        # save plot with dataset name
        ggsave(paste0("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/abundance_", curr_dataset, ".pdf"), 
               plot, 
               width=8, 
               height=6,
               dpi=500)
}

# make a plot for each dataset
lapply(unique(abund_df$dataset), plot_abund_func)

########################################################
# Plot showing the bray-curtis distances
########################################################

bc_df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/abundance_bc_distances.csv")

# title of facets
titles <- c("human_gut" = "Human Gut Reads (n=200k)",
            "aquatic" = "Aquatic Reads (n=200k)",
            "soil" = "Soil Reads (n=200k)")

plot <- ggplot(bc_df, aes(x=region, y=bcdist, group=method)) +
        geom_line(aes(color=method), size=0.75) +
        geom_point(aes(color=method), size=1.5) +
        theme_bw() +
        facet_wrap(~dataset, labeller=as_labeller(titles)) +
        theme(axis.text.x=element_text(angle=0,size=12,color="black"),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14),
              axis.text.y=element_text(size=12, color="black"),
              strip.text=element_text(face="bold", size=12),
              legend.title=element_text(size=14),
              legend.text=element_text(size=12)) +
        labs(x="Region", y="Bray-Curtis Dissimilarity") +
        scale_color_discrete(name="Method", labels=c("Cliffy", "Kraken2")) +
        scale_x_discrete(labels=c("V1-V2", "V3-V4", "V4", "V4-V5")) +
        scale_y_continuous(breaks=seq(0,1,0.02), limits=c(0,0.20))
plot
ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/bc_distances.pdf", 
       plot, 
       width=10, 
       height=4,
       dpi=500)