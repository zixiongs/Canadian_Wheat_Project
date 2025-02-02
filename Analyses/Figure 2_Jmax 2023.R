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
Jmax <- dplyr::select(Database, T_Treat, CO2, Genotype, Jmax)
Jmax <- na.omit(Jmax)

JmaxA1614_KE07 <- dplyr::filter(Jmax, Genotype == "A1614_KE07")

a1 <- ggplot(JmaxA1614_KE07, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
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
  annotate("text", x= 0.65, y=385, label="(a)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("A1614_KE07") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a1
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxA1614_KE07)
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
JmaxAAC_Brandon <- dplyr::filter(Jmax, Genotype == "AAC_Brandon")

a2 <- ggplot(JmaxAAC_Brandon, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(b)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Brandon") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a2
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Brandon)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------------
#--------------------------------------------------------------
JmaxAAC_Concord <- dplyr::filter(Jmax, Genotype == "AAC_Concord")

a3 <- ggplot(JmaxAAC_Concord, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(c)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Concord") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a3
#-------------------------------------------------------------
m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Concord)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------------
JmaxAAC_Goodwin <- dplyr::filter(Jmax, Genotype == "AAC_Goodwin")

a4 <- ggplot(JmaxAAC_Goodwin, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(d)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Goodwin") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a4
#-------------------------------------------------------------
m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Goodwin)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------------
JmaxAAC_Grainland <- dplyr::filter(Jmax, Genotype == "AAC_Grainland")

a5 <- ggplot(JmaxAAC_Grainland, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
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
  annotate("text", x= 0.65, y=385, label="(e)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Grainland") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a5
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Grainland)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#------------------------------------------------------------
JmaxAAC_Hodge <- dplyr::filter(Jmax, Genotype == "AAC_Hodge")

a6 <- ggplot(JmaxAAC_Hodge, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(f)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Hodge") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a6
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Hodge)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#-----------------------------------------------------------
JmaxAAC_Leroy <- dplyr::filter(Jmax, Genotype == "AAC_Leroy")

a7 <- ggplot(JmaxAAC_Leroy, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(g)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Leroy") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a7
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Leroy)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2

#----------------------------------------------------------
JmaxAAC_Magnet <- dplyr::filter(Jmax, Genotype == "AAC_Magnet")

a8 <- ggplot(JmaxAAC_Magnet, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(h)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Magnet") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a8
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Magnet)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2

#---------------------------------------------------------
JmaxAAC_Prevail <- dplyr::filter(Jmax, Genotype == "AAC_Prevail")

a9 <- ggplot(JmaxAAC_Prevail, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size = 12), 
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
  annotate("text", x= 0.65, y=385, label="(i)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Prevail") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a9
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Prevail)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#---------------------------------------------------------
JmaxAAC_Warman <- dplyr::filter(Jmax, Genotype == "AAC_Warman")

a10 <- ggplot(JmaxAAC_Warman, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(j)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Warman") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a10
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Warman)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#---------------------------------------------------------
JmaxAAC_Wheatland <- dplyr::filter(Jmax, Genotype == "AAC_Wheatland")

a11 <- ggplot(JmaxAAC_Wheatland, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(k)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("AAC_Wheatland") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a11
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxAAC_Wheatland)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------
JmaxBW1094 <- dplyr::filter(Jmax, Genotype == "BW1094")

a12 <- ggplot(JmaxBW1094, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(l)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("BW1094") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a12
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxBW1094)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------
JmaxBW5089 <- dplyr::filter(Jmax, Genotype == "BW5089")

a13 <- ggplot(JmaxBW5089, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size = 12), 
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
  annotate("text", x= 0.65, y=385, label="(m)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("BW5089") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a13
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxBW5089)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------
JmaxDT2016 <- dplyr::filter(Jmax, Genotype == "DT2016")

a14 <- ggplot(JmaxDT2016, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(n)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("DT2016") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a14
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxDT2016)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#--------------------------------------------------------
JmaxDT2030 <- dplyr::filter(Jmax, Genotype == "DT2030")

a15 <- ggplot(JmaxDT2030, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(o)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("DT2030") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a15
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxDT2030)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#-------------------------------------------------------
JmaxHY2136 <- dplyr::filter(Jmax, Genotype == "HY2136")

a16 <- ggplot(JmaxHY2136, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
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
  annotate("text", x= 0.65, y=385, label="(p)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("HY2136") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a16
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxHY2136)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#-------------------------------------------------------
JmaxPT4002 <- dplyr::filter(Jmax, Genotype == "PT4002")

a17 <- ggplot(JmaxPT4002, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  labs(x=expression('Treatment temperature'~'('*~degree*C*')')) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black")),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Add ticks at the bottom x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=385, label="(q)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("PT4002") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a17
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxPT4002)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#-------------------------------------------------------
JmaxStettler <- dplyr::filter(Jmax, Genotype == "Stettler")

a18 <- ggplot(JmaxStettler, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  labs(x=expression('Treatment temperature'~'('*~degree*C*')')) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black")),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Add ticks at the bottom x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=385, label="(r)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("Stettler") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a18
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxStettler)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#------------------------------------------------------
JmaxStrongfield <- dplyr::filter(Jmax, Genotype == "Strongfield")

a19 <- ggplot(JmaxStrongfield, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  labs(x=expression('Treatment temperature'~'('*~degree*C*')')) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black")),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Add ticks at the bottom x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=385, label="(s)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("Strongfield") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a19
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxStrongfield)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#------------------------------------------------------
JmaxTranscend <- dplyr::filter(Jmax, Genotype == "Transcend")

a20 <- ggplot(JmaxTranscend, aes(x=T_Treat, y=Jmax, fill = CO2)) +
  geom_boxplot(width=0.24) +
  geom_point(position=position_dodge(width=0.24)) +
  scale_fill_manual(values = c("white", "gray")) +
  theme_classic() +
  labs(y=expression(italic('J')[max25]~(mu*'mol'~'m'^-2~'s'^-1)), size = 5) +
  labs(x=expression('Treatment temperature'~'('*~degree*C*')')) +
  scale_y_continuous(breaks=seq(100,400,50), limits = c(100,400),expand=c(0,0), sec.axis = dup_axis(name = NULL, labels = NULL)) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size =12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black")),
        axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 12, margin=margin(3,3,3,3,"pt"), colour=c("black", "black", "black", "black", "black", "black"))) +
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"), 
        axis.line = element_line(colour = "black", size = 0.45),
        axis.ticks.length = unit(0, "cm"),  # Remove ticks on all axes
        axis.ticks.length.x.top = unit(0, "cm"),  # Remove ticks on the top x-axis
        axis.ticks.length.y.right = unit(0, "cm"),
        axis.ticks.length.x.bottom = unit(0.25, "cm"),  # Add ticks at the bottom x-axis
        axis.ticks.length.y.left = unit(0.25, "cm")) +
  annotate("text", x= 0.65, y=385, label="(t)",size = 5) +
  theme(legend.position = c(0.65, 0.85),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  ggtitle("Transcend") +
  theme(plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 12, colour = "red")) +
  geom_segment(aes(x = Inf, xend = -Inf, y = Inf, yend = Inf), 
               color = "black", size = 0.5, lineend = "round")
a20
#-------------------------------------------------------------

m3a <- aov(Jmax ~ T_Treat*CO2, data = JmaxTranscend)
summary(m3a)

#plot(m3a)

marginalm3a = lsmeans(m3a,
                      ~ Elevation+Succ)
CLDm2 = cld(marginalm3a,
            alpha=0.05,
            Letters=letters,
            adjust="Bonferroni")
CLDm2
#------------------------------------------------------
p4 = plot_grid(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, nrow = 5, ncol = 4)
p4

# Open a PDF device
pdf("p4.pdf", width = 16, height = 16)

# Plot and save the plot
print(p4)

# Close the PDF device
dev.off()

