library(ggplot2)
library(tidyr)
library(dplyr)

########################################################
# Plot monotonic increases for each digestion type
########################################################

df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_1/combined_increases.csv", 
               header=FALSE,
               col.names=c("digestion_type", "x", "inc_from_left", "inc_from_right")) 

# pivot longer and filter
df_long <- df %>% pivot_longer(cols = c(inc_from_left, inc_from_right), names_to = "variable", values_to = "value")
df_filt <- df_long %>% filter(x<25)

titles <- c("dna_minimizers" = "DNA Minimizers", 
            "minimizers" = "Minimizers",
            "full_text" = "Full Text")

plot1 <- ggplot(df_filt, aes(x=factor(x), y=value, fill=variable)) +
         theme_bw() +
         facet_wrap(~digestion_type, labeller=as_labeller(titles), ncol=1) +
         geom_bar(stat="identity", position=position_dodge(), color="black") +
         scale_x_discrete(breaks=seq(0,25,2)) +
         #scale_y_continuous(breaks=seq(0,5000000,500000), labels = scales::scientific) +
         labs(x="Number of Monotonic Increases", y="Number of Profiles") +
         theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
               axis.title.x=element_text(size=12),
               axis.title.y=element_text(size=12),
               axis.text=element_text(size=10, color="black"),
               strip.text=element_text(size=12, face="bold")) +
         scale_fill_discrete(name="Number of Increases", labels=c("Left to Right Increases", "Right to Left Increases")) 
plot1

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_1/minimizers_inc_plot.pdf", 
       plot=plot1, 
       dpi=800, 
       device="pdf", 
       width=6, 
       height=8)

########################################################
# Plot CDF functions for monotonic increases
########################################################

df <- read.csv("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_1/combined_increases.csv", 
               header=FALSE,
               col.names=c("digestion_type", "x", "inc_from_left", "inc_from_right")) 

df$total <- df$inc_from_left + df$inc_from_right
df <- df %>%
      arrange(digestion_type, x) %>%
      group_by(digestion_type) %>%
      mutate(cdf = cumsum(total) / sum(total))
df <- df %>% filter(x<25)

titles <- c("dna_minimizers" = "DNA Minimizers", 
            "minimizers" = "Minimizers",
            "full_text" = "Full Text")

plot2 <- ggplot(df, aes(x=factor(x), y=cdf, group=digestion_type)) +
          theme_bw() +
          facet_wrap(~digestion_type, labeller=as_labeller(titles)) +
          geom_line() +
          geom_point() +
          geom_area(alpha=0.4) +
          scale_x_discrete(breaks=seq(0,25,2)) +
          #scale_y_continuous(breaks=seq(0,5000000,500000), labels = scales::scientific) +
          labs(x="Number of Monotonic Increases", y="Percent") +
          theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                axis.text=element_text(size=10, color="black"),
                strip.text=element_text(size=12, face="bold"))  
plot2

ggsave("/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_1/minimizers_inc_plot_cdfs.pdf", 
       plot=plot2, 
       dpi=800, 
       device="pdf", 
       width=10, 
       height=4)