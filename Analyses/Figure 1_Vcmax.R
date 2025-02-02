library(dplyr)
library(ggplot2)
library(GGally)
library(MASS)
library(plotrix)
library(gplots)
library(plyr)
library(RColorBrewer)
library(cowplot)
library(ggthemes)
library(gridExtra)
library(grid)
library(lattice)
library(lsmeans)
library(emmeans)
library(nlme)
library(lme4)
library(multcomp)
library(forcats)
library(tidyverse)
library(ggpattern)
library(lmerTest)
library(car)
#------------------------------------------------------------------------------------
Database <- read.csv("Biotron2023.csv")
str(Database)
Vcmax <- dplyr::select(Database, T_Treat, CO2, Genotype, Vcmax)
Vcmax <- na.omit(Vcmax)

VcmaxMean <- ddply(Vcmax, c("T_Treat", "CO2", "Genotype"), summarise,
                          N    = length(Vcmax),
                          mean = mean(Vcmax),
                          sd   = sd(Vcmax),
                          se   = sd / sqrt(N))
names(VcmaxMean) <- c("T_Treat", "CO2", "Genotype", "N", "Vcmax", "sd.Vcmax", "se.Vcmax")

#ToptASpMeanSuccMean <- ddply(ToptASpMean, c("Elevation", "Succ"), summarise,
                     #N    = length(ToptA),
                     #mean = mean(ToptA),
                    # sd   = sd(ToptA),
                     #se   = sd / sqrt(N))
#names(ToptASpMeanSuccMean) <- c("Elevation", "Succ", "N", "ToptA", "sd.ToptA", "se.ToptA")

ToptASpMean$Elevation <- factor(ToptASpMean$Elevation, levels = c("HE", "ME", "LE"))

a1 <- ggplot(ToptASpMean, aes(x=Elevation, y=ToptA, fill = Succ)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y = expression(italic('T')[optA] ~ '(' * degree*C * ')')) +
  scale_y_continuous(breaks=seq(15,35,5), limits = c(15,35),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size=12), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black"))) + 
  theme(axis.title.x = element_text(size=12, margin=margin(3,3,3,3,"pt")), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill="white"), plot.margin = unit(c(0.5, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=34, label="(a)",size = 5) +
  theme(legend.position = c(0.85, 0.9),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank())
a1
#-------------------------------------------------------------

m3a <- aov(ToptA ~ Elevation*Succ, data = ToptASpMean)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                     ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#------------------------
str(Database)
ToptA_obs <- dplyr::select(Database, Site, Elevation, Plot, Line, Species, Succ, sunit, ToptA_obs)
ToptA_obs <- na.omit(ToptA_obs)

model <- lmer(ToptA_obs ~ Elevation * Succ + (1 | Species) + (1 | sunit), data = ToptA_obs)
Anova(model, test="F")

marginalm3a = lsmeans(model,
                      ~ Elevation*Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------------
Database <- read.csv("Project_Database_no_FS.csv")
str(Database)

Aopt <- dplyr::select(Database, Site, Elevation, Plot, Line, Species, Succ, Aopt_obs)
Aopt <- na.omit(Aopt)

AoptSpMean <- ddply(Aopt, c("Site", "Elevation", "Succ", "Species"), summarise,
                    N    = length(Aopt_obs),
                    mean = mean(Aopt_obs),
                    sd   = sd(Aopt_obs),
                    se   = sd / sqrt(N))
names(AoptSpMean) <- c("Site", "Elevation", "Succ", "Species", "N", "Aopt", "sd.Aopt", "se.Aopt")

#AoptSpMeanSuccMean <- ddply(AoptSpMean, c("Elevation", "Succ"), summarise,
                            #N    = length(Aopt),
                            #mean = mean(Aopt),
                            #sd   = sd(Aopt),
                            #se   = sd / sqrt(N))
#names(AoptSpMeanSuccMean) <- c("Elevation", "Succ", "N", "Aopt", "sd.Aopt", "se.Aopt")

AoptSpMean$Elevation <- factor(AoptSpMean$Elevation, levels = c("HE", "ME", "LE"))

a2 <- ggplot(AoptSpMean, aes(x=Elevation, y=Aopt, fill = Succ)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('A')[opt]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(0,25,5), limits = c(0,25),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size=12), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black"))) + 
  theme(axis.title.x = element_text(size=12, margin=margin(3,3,3,3,"pt")), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill="white"), plot.margin = unit(c(0.5, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=24, label="(b)",size = 5) +
  theme(legend.position = c(0.65, 0.9),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank())
a2
#-------------------------------------------------------------
#-------------------------------------------------------------
m3a <- aov(Aopt ~ Elevation*Succ, data = AoptSpMean)
summary(m3a)

#plot(m3a)
#bartlett.test(ToptV ~ interaction(Elevation, Succ), data = ToptVSpMean)
#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#-----------------------------------------------------------
str(Database)

Aopt <- dplyr::select(Database, Site, Elevation, Plot, Line, Species, Succ, Aopt_obs)
Aopt <- na.omit(Aopt)

Aopt <- dplyr::mutate(Aopt, LogAopt = log10(Aopt_obs))

model1 <- lmer(LogAopt ~ Elevation * Succ + (1 | Species), data = Aopt)
Anova(model1, test="F")

#--------------------------------------------------------------
p4 = plot_grid(a1, a2, nrow = 1, ncol = 2)
p4

# Open a PDF device
pdf("p4.pdf", width = 8, height = 4)

# Plot and save the plot
print(p4)

# Close the PDF device
dev.off()

