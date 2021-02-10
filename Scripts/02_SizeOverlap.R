# The code below estimates linear relationships among fork length (cm), gape height (mm), and gape width (mm) for Pacific Halibut (PH) and Arrowtooth Flounder (ATF) and tests for species-specific differences in allometric growth. Specimens were collected using hook-and-line near Juneau, AK in 2015 and 2016. See Barnes et al. (2021) for details.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(overlapping)
require(ggplot2)
####################################################################
# Import and prep field data:
PredData = read.csv("Data/PH_ATF_FieldData_SEAK.csv")

# Remove unknown species and those not measured:
PredData = subset(PredData, Species != "UNK" & Length..cm. != "NA")

# Remove spring and fall months with very few samples (focus = summer):
PredData = subset(PredData, Month != 5 & Month != 9)

# Consolidate individual capture locations into sampling sites:
PredData$Location = with(PredData, 
  ifelse(Grid..  < 1, "Chatham.Strait",
  ifelse(Grid.. <= 3, "Lynn.Canal",
  ifelse(Grid.. <= 7, "Favorite-Saginaw.Channel",
  ifelse(Grid.. <= 8, "Pt.Howard",
  ifelse(Grid.. <= 9, "Funter.Bay",
  ifelse(Grid.. <= 10, "Other", 
  ifelse(Grid.. <= 11, "Pt.Couverden", 
  ifelse(Grid.. <= 13, "Icy.Strait", "Other")))))))))

# Remove locations outside of spatial extent, order others north to south, and rename:
LFD = subset(PredData, Location != "Other" & Location != "Chatham.Strait")
LFD$Location = ordered(LFD$Location, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard", "Funter.Bay", "Pt.Couverden", "Icy.Strait"))
levels(LFD$Location) = c("Lynn Canal", "Favorite-Saginaw Channels", "Point Howard", "Funter Bay", "Point Couverden", "Icy Strait")

# Select columns of interest:
LFD = LFD[ , c("Location","Species","PredatorID","Length..cm.", "Year")]
  LFD = na.omit(unique(LFD, incomparables = FALSE))
LFD$FL = round(LFD$Length..cm.) # round fork length to nearest whole number
####################################################################
# Predict gape from length (see '01_AllometricGrowth' for linear relationships):
PH.data = subset(LFD, Species == "PH")
  PH.data$GH = round((PH.data$Length..cm. * 1.11497), 0)
  PH.data$GW = round((PH.data$Length..cm. * 1.19663), 0)
ATF.data = subset(LFD, Species == "ATF")  
  ATF.data$GH = round((ATF.data$Length..cm. * 2.14768), 0)
  ATF.data$GW = round((ATF.data$Length..cm. * 2.06384), 0)
All.data = rbind(PH.data, ATF.data) # recombine
  All.data$Species = as.factor(All.data$Species)
####################################################################
# Estimate overlap in size by species, year, site, and metric (values added to Fig. S2):
PH_ATF_Syr.fl = data.frame()
PH_ATF_Syr.gh = data.frame()
PH_ATF_Syr.gw = data.frame()

for (i in unique(All.data$Year)) {
  LFD_yr = subset(All.data, Year == i)
for (j in unique(LFD_yr$Location)) {
  LFD_yr_loc = subset(LFD_yr, Location == j)
  
FL.yr_vect = list(ATF = LFD_yr_loc$FL[LFD_yr_loc$Species == "ATF"], PH = LFD_yr_loc$FL[LFD_yr_loc$Species == "PH"]) 
  FL.yr_over = overlap(FL.yr_vect, nbins = 1000)
GH.yr_vect = list(ATF = LFD_yr_loc$GH[LFD_yr_loc$Species == "ATF"], PH = LFD_yr_loc$GH[LFD_yr_loc$Species == "PH"]) 
  GH.yr_over = overlap(GH.yr_vect, nbins = 1000)
GW.yr_vect = list(ATF = LFD_yr_loc$GW[LFD_yr_loc$Species == "ATF"], PH = LFD_yr_loc$GW[LFD_yr_loc$Species == "PH"]) 
  GW.yr_over = overlap(GW.yr_vect, nbins = 1000)

PH_ATF_Syr.fl = rbind(PH_ATF_Syr.fl, data.frame(unique(LFD_yr_loc$Year), unique(LFD_yr_loc$Location), FL.yr_over$OV))
PH_ATF_Syr.gh = rbind(PH_ATF_Syr.gh, data.frame(unique(LFD_yr_loc$Year), unique(LFD_yr_loc$Location), GH.yr_over$OV))
PH_ATF_Syr.gw = rbind(PH_ATF_Syr.gw, data.frame(unique(LFD_yr_loc$Year), unique(LFD_yr_loc$Location), GW.yr_over$OV))
  }
}

colnames(PH_ATF_Syr.fl) = c("Year", "Location", "S")
colnames(PH_ATF_Syr.gh) = c("Year", "Location", "S")
colnames(PH_ATF_Syr.gw) = c("Year", "Location", "S")

# Test for differences in overlap (S) by year, location, and size metric:
PH_ATF_Syr.fl_test = aov(S ~ Year, data = PH_ATF_Syr.fl)
  summary(PH_ATF_Syr.fl_test) # no difference between years --> pool
PH_ATF_Syr.gh_test = aov(S ~ Year, data = PH_ATF_Syr.gh)
  summary(PH_ATF_Syr.gh_test) # no difference between years --> pool
PH_ATF_Syr.gw_test = aov(S ~ Year, data = PH_ATF_Syr.gw)
  summary(PH_ATF_Syr.gw_test) # no difference between years --> pool
####################################################################
# Estimate overlap in size by species, site, and metric (years pooled; values added to Fig. 3):
PH_ATF_S.fl = data.frame()
PH_ATF_S.gh = data.frame()
PH_ATF_S.gw = data.frame()
  
for (j in unique(All.data$Location)) {
  LFD_loc = subset(All.data, Location == j)
    
FL_vect = list(ATF = LFD_loc$FL[LFD_loc$Species == "ATF"], PH = LFD_loc$FL[LFD_loc$Species == "PH"])      
  FL_over = overlap(FL_vect, nbins = 1000)
GH_vect = list(ATF = LFD_loc$GH[LFD_loc$Species == "ATF"], PH = LFD_loc$GH[LFD_loc$Species == "PH"])      
  GH_over = overlap(GH_vect, nbins = 1000)
GW_vect = list(ATF = LFD_loc$GW[LFD_loc$Species == "ATF"], PH = LFD_loc$GW[LFD_loc$Species == "PH"]) 
  GW_over = overlap(GW_vect, nbins = 1000)
      
PH_ATF_S.fl = rbind(PH_ATF_S.fl, data.frame(unique(LFD_loc$Location), FL_over$OV))
PH_ATF_S.gh = rbind(PH_ATF_S.gh, data.frame(unique(LFD_loc$Location), GH_over$OV))
PH_ATF_S.gw = rbind(PH_ATF_S.gw, data.frame(unique(LFD_loc$Location), GW_over$OV))
    }
  
colnames(PH_ATF_S.fl) = c("Location", "S")
colnames(PH_ATF_S.gh) = c("Location", "S")
colnames(PH_ATF_S.gw) = c("Location", "S")
  
# Plot size frequency distributions:
levels(All.data$Species) = c(" ATF", " PH")

plot.settings = function() {
  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18),
          panel.background = element_rect(colour="black", size=0.5, line="solid"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_text(vjust = -0.13, size = 12, family="Arial"),
          axis.title.y = element_text(vjust = 2.5, size = 12.5, family="Arial"),
          axis.text = element_text(family="Arial", size=10),
          legend.title = element_blank(),
          legend.text = element_text(family = "Arial", size = 12), legend.text.align = 0,
          legend.background = element_rect(fill = "transparent"),
          legend.key.width = unit(4, "mm"),
          legend.key.height = unit(4, "mm"),
          panel.spacing.y = unit(0.5, "lines"),
          strip.background = element_rect(colour="white", fill = "white"),
          strip.text.x = element_text(family = "Arial", size = 11.5, hjust=0),
          strip.text.y = element_text(family = "Arial", size = 8, hjust=0)) }

# Fig. 3 (left):  
FL.dist = ggplot(All.data, aes(x = FL, fill = Species)) +   
  geom_density(adjust = 0.5, alpha = 0.75) +
  facet_wrap(~ Location, ncol = 1) +
  labs(x = "Fork Length (cm)", y = " Relative Frequency") +
  plot.settings() + 
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.12), breaks = seq(0, 0.1, 0.05)) +
  scale_x_continuous(expand = c(0, 0), limits = c(30,165), breaks = seq(30, 160, 20)) 
ggsave(filename = "Plots/FL_dist.png", plot = FL.dist, width = 3.3, height = 5.7, units = "in")

# Fig. S2 (left):
FL.dist_yr = ggplot(All.data, aes(x = FL, fill = Species)) +   
  geom_density(adjust = 0.5, alpha = 0.75) +
  plot.settings() + 
  facet_grid(Location ~ Year) +
  labs(x = "Fork Length (cm)", y = " Relative Frequency") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.12), breaks = seq(0, 0.1, 0.05)) +
  scale_x_continuous(expand = c(0, 0), limits = c(30,165), breaks = seq(30, 160, 40)) 
ggsave(filename = "Plots/FL_dist_yr.png", plot = FL.dist_yr, width = 4, height = 6, units = "in")

# Fig. 3 (center):
GH.dist = ggplot(All.data, aes(x = GH, fill = Species)) +  
  geom_density(adjust = 0.5, alpha = 0.75) +
  plot.settings() + 
  facet_wrap(~ Location, ncol = 1) +
  labs(x = "Gape Height (mm)", y = "Relative Frequency") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.068), breaks = seq(0, 0.06, 0.03)) +
  scale_x_continuous(expand = c(0, 0), limits = c(30,202), breaks = seq(30, 190, 20)) 
ggsave(filename = "Plots/GH_dist.png", plot = GH.dist, width = 3.3, height = 5.7, units = "in")

# Fig. S2 (center):
GH.dist_yr = ggplot(All.data, aes(x = GH, fill = Species)) +  
  geom_density(adjust = 0.5, alpha = 0.75) +
  plot.settings() +
  facet_grid(Location ~ Year) +
  labs(x = "Gape Height (mm)", y = "Relative Frequency") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.068), breaks = seq(0, 0.06, 0.03)) +
  scale_x_continuous(expand = c(0, 0), limits = c(30,202), breaks = seq(30, 190, 40)) 
ggsave(filename = "Plots/GH_dist_yr.png", plot = GH.dist_yr, width = 4, height = 6, units = "in")

# Fig. 3 (right):
GW.dist = ggplot(All.data, aes(x = GW, fill = Species)) +   
  geom_density(adjust = 0.5, alpha = 0.75) +
  facet_wrap(~ Location, ncol = 1) +
  plot.settings() +
  labs(x = "Gape Width (mm)", y = "Relative Frequency") +
  theme(legend.position = c(0.11, 0.07), legend.text = element_text(size=9)) +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.068), breaks = seq(0, 0.06, 0.03)) +
  scale_x_continuous(expand = c(0, 0), limits = c(30,202), breaks = seq(30, 190, 20)) 
ggsave(filename = "Plots/GW_dist.png", plot = GW.dist, width = 3.3, height = 5.7, units = "in")

# Fig. S2 (right):
GW.dist_yr = ggplot(All.data, aes(x = GW, fill = Species)) +   
  geom_density(adjust = 0.5, alpha = 0.75) +
  plot.settings() +
  facet_grid(Location ~ Year) +
  labs(x = "Gape Width (mm)", y = "Relative Frequency") +
  theme(legend.position = c(0.89, 0.97), legend.text = element_text(size=9)) +
  scale_fill_manual(breaks = c(" ATF", " PH"), values = c("blue", "red")) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.068), breaks = seq(0, 0.07, 0.03)) +
  scale_x_continuous(expand = c(0, 0), limits = c(29,202), breaks = seq(30, 190, 40)) 
ggsave(filename = "Plots/GW_dist_yr.png", plot = GW.dist_yr, width = 4, height = 6, units = "in")