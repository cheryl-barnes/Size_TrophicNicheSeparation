# The code below estimates linear relationships among fork length (cm), gape height (mm), and gape width (mm) for Pacific Halibut (PH) and Arrowtooth Flounder (ATF) and tests for species-specific differences in allometric growth. Specimens were collected using hook-and-line near Juneau, AK in 2015. See Barnes et al. (2021) for details.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(ggplot2)
###########################################################  
# Import gape height, gape width, and fork length measurements (obtained from freshly dead fishes):
Data = read.csv("Data/PH_ATF_Morphometrics_SEAK.csv")
  Data$Species = as.factor(Data$Species)
  levels(Data$Species) = c(" ATF", " PH")

# Test for species-specific differences in allometric growth:
# gape height-fork length:
ANCOVA_gh = aov(Final.Fork.Length..cm. ~ -1 + Gape.Height..mm. + Species, data=Data)
  summary(ANCOVA_gh)
# gape width-fork length:
ANCOVA_gw = aov(Final.Fork.Length..cm. ~ -1 + Gape.Width..mm. + Species, data=Data)
  summary(ANCOVA_gw)
# gape width-gape height:
ANCOVA_both = aov(Gape.Height..mm. ~ -1 + Gape.Width..mm. + Species, data=Data)
  summary(ANCOVA_both)

########################################################### 
# Quantify species-specific linear relationships among size metrics (Table 1):

# Pacific Halibut:
PHdata = subset(Data, Species==" PH")
PH.fl.gh = lm(Gape.Height..mm. ~ -1 + Final.Fork.Length..cm., PHdata)
  summary(PH.fl.gh)
PH.fl.gw = lm(Gape.Width..mm. ~ -1 + Final.Fork.Length..cm., PHdata)
  summary(PH.fl.gw)
PH.gh.gw = lm(Gape.Width..mm. ~ -1 + Gape.Height..mm., PHdata)
  summary(PH.gh.gw)

PHdata$predGH = predict(PH.fl.gh) # predicted gape height (GH) based on fork length (FL)
PHdata$predGW = predict(PH.fl.gw) # predicted gape width (GW) based on fork length (FL)

# Arrowtooth Flounder:
ATFdata = subset(Data, Species==" ATF")
ATF.fl.gh = lm(Gape.Height..mm. ~ -1 + Final.Fork.Length..cm., ATFdata)
  summary(ATF.fl.gh)
ATF.fl.gw = lm(Gape.Width..mm. ~ -1 + Final.Fork.Length..cm., ATFdata)
  summary(ATF.fl.gw)
ATF.gh.gw = lm(Gape.Width..mm. ~ -1 + Gape.Height..mm., ATFdata)
  summary(ATF.gh.gw)

ATFdata$predGH = predict(ATF.fl.gh) # predicted gape height (GH) based on fork length (FL)
ATFdata$predGW = predict(ATF.fl.gw) # predicted gape width (GW) based on fork length (FL)

# Recombine the two data frames and plot:
PH_ATFplot = rbind(PHdata, ATFdata)

plot.settings = function() {
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18),
        panel.background = element_rect(colour="black", fill="white", size=1, line="solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(vjust=-0.13, size=14),
        axis.title.y = element_text(vjust=2.0, size=14),
        axis.text = element_text(family="Arial", size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.key.width = unit(1.9, "mm"),
        legend.key.height = unit(4.5, "mm"),
        strip.background = element_rect(colour="white"),
        strip.text = element_text(size=12)) }

# Fig. 2 (left):
LM.FL.GH = ggplot(Data, aes(x=Final.Fork.Length..cm., y=Gape.Height..mm., colour = Species, group = Species)) +
  plot.settings() +
  geom_point() + 
  geom_smooth(aes(color = Species, fill = Species), method = "lm", se = T) +
  scale_color_manual(breaks = c("PH", "ATF"), values = c("royalblue", "indianred2", "indianred2", "royalblue")) +
  scale_fill_manual(values = c("royalblue", "indianred2", "indianred2", "royalblue")) +
  xlab("Fork Length (cm)") + ylab("Gape Height (mm)") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(50, 150), breaks=c(56,76,96,116,136))
ggsave(filename="Plots/LM_FL_GH.png", plot=LM.FL.GH, width=3, height=4, units="in")

# Fig. 2 (right):
LM.FL.GW = ggplot(Data, aes(x=Final.Fork.Length..cm., y=Gape.Width..mm., color = Species)) +
  geom_point() + 
  geom_smooth(aes(color = Species, fill = Species), method = "lm", se = T) +
  plot.settings() +
  scale_color_manual(values = c("royalblue", "indianred2", "indianred2")) +
  scale_fill_manual(values = c("royalblue", "indianred2")) +
  xlab("Fork Length (cm)") + ylab("Gape Width (mm)") +
  theme(legend.position = c(0.85, 0.08)) +
  scale_y_continuous(limits = c(50, 150), breaks=c(56,76,96,116,136))
ggsave(filename="Plots/LM_FL_GW.png", plot=LM.FL.GW, width=3, height=4, units="in")