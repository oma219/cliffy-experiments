

library(ggplot2)
library(dplyr)

df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/accuracy_debugging.csv", header=TRUE)

########################################################
# Plot #1: Density of minimizers with no hits
########################################################
df <- df %>% group_by(dataset, region, type) %>% mutate(total_no_hits = sum(no_hits))
df$percent_no_hits <- df$no_hits / df$total_no_hits

titles <- c("aquatic.V1_V2" = "Aquatic V1-V2 Reads (n=200k)", 
            "aquatic.V4_V4" = "Aquatic V4 Reads (n=200k)", 
            "soil.V1_V2" = "Soil V1-V2 Reads (n=200k)",
            "soil.V4_V4" = "Soil V4 Reads (n=200k)")

plot1 <- ggplot(df, aes(x=x, ymin=0, ymax=percent_no_hits, fill=type)) +
       geom_ribbon(alpha = 0.4) +
       geom_line(aes(y=percent_no_hits, color=type)) +
       labs(x = "# of Minimizers with No Hits", y = "Density", fill = "Method") +
       theme_bw() +
       facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
       theme(axis.text=element_text(size=10, color="black"),
             strip.text=element_text(size=12, face="bold"),
             axis.title.x=element_text(size=12),
             axis.title.y=element_text(size=12),
             legend.text=element_text(size=10),
             legend.title=element_text(size=12)) +
       scale_fill_discrete(name="Method") +
       scale_color_discrete(name="Method") 

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_1.pdf", 
       plot1, 
       width=8, 
       height=6,
       dpi=800)

########################################################
# Plot #2: Density of minimizers to correct genus
########################################################    
df <- df %>% group_by(dataset, region, type) %>% mutate(total_hits_to_correct_genus = sum(hits_to_correct_genus))
df$percent_correct_genus_hits <- df$hits_to_correct_genus / df$total_hits_to_correct_genus

plot2 <-ggplot(df, aes(x=x, ymin=0, ymax=percent_correct_genus_hits, fill=type)) +
      geom_ribbon(alpha = 0.4) +
      geom_line(aes(y=percent_correct_genus_hits, color=type)) +
      labs(x = "# of Minimizers with Correct Genus Hits", y = "Density", fill = "Method") +
      theme_bw() +
      facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
      theme(axis.text=element_text(size=10, color="black"),
            strip.text=element_text(size=12, face="bold"),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            legend.text=element_text(size=10),
            legend.title=element_text(size=12)) +
      scale_fill_discrete(name="Method") +
      scale_y_continuous(limits=c(0,0.05)) +
      scale_color_discrete(name="Method")

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_2.pdf", 
       plot2, 
       width=8, 
       height=6,
       dpi=800)

#############################################################
# Plot #3: Density of lengths that correspond to correct doc
############################################################# 

df <- df %>% group_by(dataset, region, type) %>% mutate(total_lengths_contain_doc = sum(lengths_contain_doc))
df$percent_lengths_contain_doc <- df$lengths_contain_doc / df$total_lengths_contain_doc

plot3 <-ggplot(df, aes(x=x, ymin=0, ymax=percent_lengths_contain_doc, fill=type)) +
        geom_ribbon(alpha = 0.4) +
        geom_line(aes(y=percent_lengths_contain_doc, color=type)) +
        labs(x = "Length of Exact Matching Containing Correct Doc", y = "Density", fill = "Method") +
        theme_bw() +
        facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
        theme(axis.text=element_text(size=10, color="black"),
              strip.text=element_text(size=12, face="bold"),
              axis.title.x=element_text(size=12),
              axis.title.y=element_text(size=12),
              legend.text=element_text(size=10),
              legend.title=element_text(size=12)) +
        scale_fill_discrete(name="Method") +
        xlim(0, 150) +
        scale_color_discrete(name="Method")

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_3.pdf", 
       plot3, 
       width=8, 
       height=6,
       dpi=800)

#########################################################################
# Plot #4: Density of lengths that correspond to correct doc exclusively
#########################################################################

df <- df %>% group_by(dataset, region, type) %>% mutate(total_lengths_exclusive_doc = sum(lengths_exclusive_doc))
df$percent_lengths_exclusive_doc <- df$lengths_exclusive_doc / df$total_lengths_exclusive_doc

plot4 <-ggplot(df, aes(x=x, ymin=0, ymax=percent_lengths_exclusive_doc, fill=type)) +
  geom_ribbon(alpha = 0.4) +
  geom_line(aes(y=percent_lengths_exclusive_doc, color=type)) +
  labs(x = "Length of Exact Matching Exclusive to Correct Doc", y = "Density", fill = "Method") +
  theme_bw() +
  facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
  theme(axis.text=element_text(size=10, color="black"),
        strip.text=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  scale_fill_discrete(name="Method") +
  #xlim(0, 150) +
  scale_color_discrete(name="Method")

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_4.pdf", 
       plot4, 
       width=8, 
       height=6,
       dpi=800)

#########################################################################
# Plot #5: Density of minimizer hits to each level of taxonomy
#########################################################################
labels <- c("Genus", "Family", "Order", "Class", "Phylum", "Domain", "Root", "No Hit")

titles <- c("aquatic.V1_V2" = "Aquatic V1-V2 Reads (n=200k)", 
            "aquatic.V4_V4" = "Aquatic V4 Reads (n=200k)", 
            "soil.V1_V2" = "Soil V1-V2 Reads (n=200k)",
            "soil.V4_V4" = "Soil V4 Reads (n=200k)")

df_filt <- df %>% filter(database == "silva_m31")

plot5 <-ggplot(df_filt, aes(x=x, ymin=0, ymax=percent_to_clade, fill=type)) +
        geom_ribbon(alpha = 0.4) +
        geom_line(aes(y=percent_to_clade, color=type)) +
        labs(x = "Major Clade", y = "Density", fill = "Method") +
        theme_bw() +
        facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
        theme(axis.text.y=element_text(size=10, color="black"),
              axis.text.x=element_text(size=10, color="black", angle=40, hjust=1),
              strip.text=element_text(size=12, face="bold"),
              axis.title.x=element_text(size=12),
              axis.title.y=element_text(size=12),
              legend.text=element_text(size=10),
              legend.title=element_text(size=12)) +
        scale_fill_discrete(name="Method") +
        scale_x_continuous(breaks = 0:7, labels = labels, limits=c(0,7)) +
        scale_color_discrete(name="Method")
plot5
ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_5.pdf", 
       plot5, 
       width=8, 
       height=6,
       dpi=800)

#########################################################################
# Plot #6: Avg size of each taxa that minimizer hits to
#########################################################################
labels <- c("Genus", "Family", "Order", "Class", "Phylum", "Domain")

plot6 <-ggplot(df, aes(x=x, ymin=0, ymax=avg_clade_size, fill=type)) +
  geom_ribbon(alpha = 0.4) +
  geom_line(aes(y=avg_clade_size, color=type)) +
  labs(x = "Major Clade", y = "Density", fill = "Method") +
  theme_bw() +
  facet_wrap(~interaction(dataset,region), labeller=as_labeller(titles)) +
  theme(axis.text.y=element_text(size=10, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=40, hjust=1),
        strip.text=element_text(size=12, face="bold"),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  scale_fill_discrete(name="Method") +
  scale_x_continuous(breaks = 0:5, labels = labels, limits=c(0,5)) +
  scale_color_discrete(name="Method")

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3_analysis/plot_6.pdf", 
       plot6, 
       width=8, 
       height=6,
       dpi=800)