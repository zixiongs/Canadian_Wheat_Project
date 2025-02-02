library(tidyverse)
library(patchwork)
library(broom) # for the tidy() function for summarising model output
source("functions/fn-plotting.R") # has my functions
#library(cowplot) 

biomass_2022 <- readRDS("data/data-input/biomass_2022.rds")

current_genotype <- biomass_2022 |> filter(genotype == "AAC_Goodwin")
aov1 <- aov(above_ground_biomass ~ temperature_level*co2_level, data = current_genotype)
#summary(aov1)

p_vals <- tidy(aov1)$p.value[1:3]

sig_labs <- paste(
  "T:", sig_level(p_vals[1]),
  "\nCO2:", sig_level(p_vals[2]),
  "\nT*CO2:", sig_level(p_vals[3])
)


biomass_2022 |> 
  filter(genotype == "AAC_Goodwin") |> 
  ggplot(aes(temperature_level, above_ground_biomass, fill = co2_level)) +
  geom_boxplot(width = 0.24, outliers = FALSE) +
  geom_point(position = position_dodge(width = 0.24)) +
  theme_classic() +
  #facet_wrap(~genotype) +
  scale_fill_manual(values = c("white", "gray")) +
  labs(y="Total above ground biomass (g)", size = 5) +
  scale_y_continuous(breaks=seq(5,40,5), limits = c(0,45),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size=12), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black")),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Add ticks at the bottom x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.7, y=42, label="(a)",size = 5) +
  annotate("text", x= 1.2, y=37, label=sig_labs, size = 4.2, hjust = 0) +
  theme(legend.position = c(0.8, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Goodwin") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")


# Writing functions -------------------------------------------------------

create_biomass_plot <- function(data, genotype_name, plot_label = "(a)") {
  genotype_data <- data |> 
    filter(genotype == genotype_name)
  
  # Then perform ANOVA on the filtered data
  aov_result <- aov(above_ground_biomass ~ temperature_level*co2_level, data = genotype_data)
  
  p_vals <- tidy(aov_result)$p.value[1:3]
  
  sig_labs <- paste(
    "T:", sig_level(p_vals[1]),
    "\nCO2:", sig_level(p_vals[2]),
    "\nT*CO2:", sig_level(p_vals[3])
  )
  
  genotype_data |> 
    ggplot(aes(temperature_level, above_ground_biomass, fill = co2_level)) +
    geom_boxplot(width = 0.24, outliers = FALSE) +
    geom_point(position = position_dodge(width = 0.24)) +
    theme_classic() +
    scale_fill_manual(values = c("white", "gray")) +
    labs(y = "Total above ground biomass (g)", size = 5) +
    scale_y_continuous(
      breaks = seq(5, 40, 5), 
      limits = c(0, 45),
      expand = c(0, 0), 
      sec.axis = dup_axis(name = NULL, labels = NULL)
    ) +
    theme(
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(
        size = 12, 
        margin = margin(3, 3, 3, 3, "pt"), 
        colour = rep("black", 4)
      ),
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 12, 
        margin = margin(3, 3, 3, 3, "pt"), 
        colour = rep("black", 6)
      )
    ) +
    theme(
      plot.background = element_rect(fill = "white"),
      plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"),
      axis.line = element_line(colour = "black", size = 0.45),
      axis.ticks.length = unit(0, "cm"),
      axis.ticks.length.x.top = unit(0, "cm"),
      axis.ticks.length.y.right = unit(0, "cm"),
      axis.ticks.length.x.bottom = unit(0.25, "cm"),
      axis.ticks.length.y.left = unit(0.25, "cm")
    ) +
    annotate("text", x= 0.7, y=42, label="(a)",size = 5) +
    annotate("text", x= 1.2, y=37, label=sig_labs, size = 4.2, hjust = 0) +
    theme(
      legend.position = c(0.8, 0.85),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.key = element_blank()
    ) +
    ggtitle(genotype_name) +
    theme(plot.title = element_text(
      vjust = 0.5, 
      hjust = 0.5, 
      size = 12, 
      colour = "red"
    )) +
    geom_segment(
      aes(x = Inf, xend = -Inf, y = Inf, yend = Inf),
      color = "black", 
      size = 0.5, 
      lineend = "round"
    )
}

create_biomass_plot(biomass_2022, "A1614-KE07", "(a)")

distinct(biomass_2022,genotype)


# Writing for-loop --------------------------------------------------------

genotypes <- distinct(biomass_2022, genotype)
labels <- letters[1:nrow(genotypes)]
plot_list <- list()

# Loop through genotypes and create plots
for(i in 1:nrow(genotypes)) {
  current_genotype <- genotypes$genotype[i]
  current_label <- paste0("(", labels[i], ")")
  
  # Create plot and store in list
  plot_list[[i]] <- create_biomass_plot(
    data = biomass_2022,
    genotype_name = current_genotype,
    plot_label = current_label
  )
}

# Don't show y-axis title except for the left most plot
show_y_labels <- rep(c(TRUE, FALSE, FALSE, FALSE), 5)
plot_list_modified <- map2(plot_list, show_y_labels, function(p, show) {
  if (!show) {
    p + theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
  } else {
    p
  }
})


p1 <- patchwork::wrap_plots(plot_list_modified, nrow = 5, ncol = 4)
pdf("outputs/plots/p1.pdf", width = 16, height = 16)
print(p1)
dev.off()
