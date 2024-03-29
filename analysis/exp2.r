
library(ggplot2)
library(viridis)

# load data
df <- read.delim("/Users/omarahmed/Downloads/work_dir/chromap_paper/exp_2/humangut_abundance.tsv", 
                 header=FALSE, 
                 sep = "\t",
                 col.names=c("lineage", "read.count", "fraction"))
df$genus <- sapply(df$lineage, function(x) tail(strsplit(x, ";")[[1]], n = 1))

# compute cumulative percentages (top of each rectangle)
df$ymax <- cumsum(df$fraction)

# compute the bottom of each rectangle
df$ymin <- c(0, head(df$ymax, n=-1))

# Reorder genus by size of fraction (descrending order)
df$genus <- with(df, reorder(genus, -fraction))

# Sample from viridis color pallete for each genus
set.seed(123)
palette <- viridis(length(unique(df$genus)))
random_colors <- sample(palette, length(unique(df$genus)))

# Set random colors for each genus so it is easier to distinguish them
colors <- setNames(random_colors, levels(df$genus))

# Use this vector of colors in the plot
plot <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2, fill=genus)) +
        geom_rect(stat="identity") +
        scale_fill_manual(values=colors, 
                          guide="legend",
                          labels=paste(levels(df$genus), ": ", df$fraction, sep="")) +
        coord_polar(theta="y") +
        xlim(c(0, 4)) +
        theme_void() +
        theme(legend.position = "bottom",
              legend.text = element_text(size=8)) +
        labs(fill="Genus")
plot
ggsave(plot=plot,
       filename="/Users/omarahmed/Downloads/work_dir/chromap_paper/exp_2/humangut_abundance.pdf", 
       device="pdf",
       width=10, 
       height=15, 
       dpi=300)


# Get a list of all *_abundance.tsv files
files <- list.files(path = "/Users/omarahmed/Downloads/work_dir/chromap_paper/exp_2",
                    pattern = "*_abundance.tsv", 
                    full.names = TRUE)

# Define a function to read and process a file
process_file <- function(file) {
  # Load data
  df <- read.delim(file, 
                   header = FALSE, 
                   sep = "\t",
                   col.names = c("lineage", "read.count", "fraction"))
  df$genus <- sapply(df$lineage, function(x) tail(strsplit(x, ";")[[1]], n = 1))

  # Compute cumulative percentages (top of each rectangle)
  df$ymax <- cumsum(df$fraction)

  # Compute the bottom of each rectangle
  df$ymin <- c(0, head(df$ymax, n=-1))

  # Reorder genus by size of fraction (descrending order)
  df$genus <- with(df, reorder(genus, -fraction))

  # Sample from viridis color pallete for each genus
  set.seed(123)
  palette <- viridis(length(unique(df$genus)))
  random_colors <- sample(palette, length(unique(df$genus)))

  # Set random colors for each genus so it is easier to distinguish them
  colors <- setNames(random_colors, levels(df$genus))

  # Use this vector of colors in the plot
  plot <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2, fill=genus)) +
          geom_rect(stat="identity") +
          scale_fill_manual(values=colors, 
                            guide="legend",
                            labels=paste(levels(df$genus), ": ", df$fraction, sep="")) +
          coord_polar(theta="y") +
          xlim(c(0, 4)) +
          theme_void() +
          theme(legend.position = "bottom",
                legend.text = element_text(size=8)) +
          labs(fill="Genus")

  # Save the plot
  ggsave(plot=plot,
         filename=paste0(gsub(".tsv", "", file), ".pdf"), 
         device="pdf",
         width=11, 
         height=15, 
         dpi=300)

  # Print a success message
  print(paste("Processed file:", file))
}

# Use lapply to read and process each file
res <- lapply(files, process_file)
