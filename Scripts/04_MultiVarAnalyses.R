# The code below tests for multivariate differences in diet compositions of Pacific Halibut (PH) and Arrowtooth Flounder (ATF) caught in Southeast Alaska (2015 and 2016). Specimens were collected using hook-and-line near Juneau, AK in 2015. See Barnes et al. (2021) for details. Specific size classes analyzed: fork length (FL) 60-69 cm; gape heights (GH) 96-115 mm and 116-135 mm; gape widths (GW) 96-115 mm and 116-135 mm.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(reshape2)
require(vegan)
require(ggplot2)

##################################################################
# Prepare diet data:
DietData = read.csv("Data/PH_ATF_DietData_SEAK.csv")
  DietData$Year = as.factor(DietData$Year)

# Remove spring and fall months with very few samples (focus = summer):
DietData = subset(DietData, Month != 5 & Month != 9)

# Consolidate individual capture locations into sampling sites:
DietData$Location = with(DietData, 
                          ifelse(Grid.. < 1, "Chatham.Strait",
                          ifelse(Grid.. <= 3, "Lynn.Canal",
                          ifelse(Grid.. <= 7, "Favorite-Saginaw.Channel",
                          ifelse(Grid.. <= 8, "Pt.Howard",
                          ifelse(Grid.. <= 9, "Funter.Bay",
                          ifelse(Grid.. <= 10, "Other", 
                          ifelse(Grid.. <= 11, "Pt.Couverden", 
                          ifelse(Grid.. <= 13, "Icy.Strait", "Other")))))))))

# Remove locations outside of spatial extent, order others north to south, and rename:
DietData = subset(DietData, Location != "Other" & Location != "Chatham.Strait")
DietData$Location = ordered(DietData$Location, levels=c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard", "Funter.Bay", "Pt.Couverden", "Icy.Strait"))

# Remove unmeasured predators:
DietData = subset(DietData, Length..cm. != "NA")

# Assign size bins based on length:
DietData$FL = round(DietData$Length..cm., 0)
DietData$FL_Bin = cut(DietData$FL, breaks = c(0, 29, 39, 49, 59, 69, 79, 89, 99, 109, 119, 129, 139, 149, 159), include.lowest = F)
levels(DietData$FL_Bin) = c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100-109", "110-119", "120-129", "130-139", "140-149", "150-159")

# Remove sites undersampled for ATF:
DietData_site = subset(DietData, Location != "Funter.Bay" & Location != "Pt.Couverden" & Location != "Icy.Strait")

# Predict gape sizes from gape-length relationships:
PH = subset(DietData_site, Species=="PH")
ATF = subset(DietData_site, Species=="ATF")  

# Predict gape height from fork length:
ATF$GH = ATF$Length..cm. * 2.14768
PH$GH = PH$Length..cm. * 1.11497

# Predict gape width from fork length:
ATF$GW = ATF$Length..cm. * 2.06384
PH$GW = PH$Length..cm. * 1.19663

# Recombine PH and ATF data:
DietData_comp = rbind(PH, ATF)

# Round gapes to nearest mm:
DietData_comp$GapeHeight = round(DietData_comp$GH, 0)
DietData_comp$GapeWidth = round(DietData_comp$GW, 0)  

# Assign gape size bins:
DietData_comp$GH_Bin = cut(DietData_comp$GapeHeight, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = FALSE)
levels(DietData_comp$GH_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

DietData_comp$GW_Bin = cut(DietData_comp$GapeWidth, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = FALSE)
levels(DietData_comp$GW_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

####################################################################
# Multivariate analyses according to fork length (FL), 60-69 cm:
DietData_FL60 = subset(DietData_comp, FL_Bin == "60-69")

# Summarize diet data:
propData_FL60 = DietData_FL60 %>%
  group_by(Species, Year, Month, Location, ID_final_Mod) %>%
  summarise(preyWT = sum(Prey.Mass..g.)) %>%
  mutate(propWT = preyWT / sum (preyWT))

# Transform data for normality:
propData_FL60 = as.data.frame(propData_FL60)
  propData_FL60$transWT = propData_FL60$propWT ^ 0.25

# Select specific columns of interest:
propData_FL60 = propData_FL60[,c("Species", "Year", "Month", "Location", "ID_final_Mod", "propWT", "transWT")]

# Set up input variables:
WideData_FL60 = dcast(propData_FL60, Species + Year + Month + Location ~ ID_final_Mod, value.var = "transWT", fun.aggregate = mean, fill=0)
  WideData_FL60 = na.omit(WideData_FL60)
WideData_FL60$Month = as.factor(WideData_FL60$Month)
droplevels(WideData_FL60$Location)

ncol(WideData_FL60); colnames(WideData_FL60)[5]; colnames(WideData_FL60)[24]
  RespVar_FL60 = WideData_FL60[ ,5:24]
ExpVar_FL60 = WideData_FL60[ ,c("Species", "Year", "Month", "Location")]

# Calculate Morisita-Horn dissimilarity index:  
M.H_FL60 = vegdist(RespVar_FL60, method = "horn") # all data

# Test for differences in multivariate dispersion among factor levels:
disp1_FL60 = betadisper(M.H_FL60, WideData_FL60$Species)
disp2_FL60 = betadisper(M.H_FL60, WideData_FL60$Year)
disp3_FL60 = betadisper(M.H_FL60, WideData_FL60$Month)
disp4_FL60 = betadisper(M.H_FL60, WideData_FL60$Location)

permutest(disp1_FL60) # no diff. disp. between species
permutest(disp2_FL60) # no diff. disp. between years
permutest(disp3_FL60) # no diff. disp. among months
permutest(disp4_FL60) # no diff. disp. among locations

# Develop individual models with single explanatory variable:
Spp_input_FL60 = adonis2(M.H_FL60 ~ Species, data = WideData_FL60, method = "horn", permutations = 999)
Spp_input_FL60

Yr_input_FL60 = adonis2(M.H_FL60 ~ Year, data = WideData_FL60, method = "horn", permutations = 999)
Yr_input_FL60

Mo_input_FL60 = adonis2(M.H_FL60 ~ Month, data = WideData_FL60, method = "horn", permutations = 999)
Mo_input_FL60

Loc_input_FL60 = adonis2(M.H_FL60 ~ Location, data = WideData_FL60, method = "horn", permutations = 999)
Loc_input_FL60

# Full model with explanatory variable order according to pseudo F and p:
Perm_full_FL60 = adonis2(M.H_FL60 ~ Year + Species + Month + Location, data=WideData_FL60, method = "horn", by = "margin", permutations = 9999); Perm_full_FL60

# Drop non-significant main effects:
Perm_mod_FL60 = adonis2(M.H_FL60 ~ Year + Species + Month , data=WideData_FL60, method = "horn", by = "margin", permutations = 9999); Perm_mod_FL60

Perm_final_FL60 = adonis2(M.H_FL60 ~ Year + Species, data=WideData_FL60, method = "horn", by = "margin", permutations = 9999); Perm_final_FL60

# Visualize differences in diet compositions:
veganCovEllipse = function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{ theta = (0:npoints) * 2 * pi/npoints
  Circle = cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov))) }

plot.theme = function() {
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        panel.spacing = unit(0.5, "lines"),
        legend.title = element_blank(),
        legend.position = c(0.885, 0.932),
        legend.background = element_blank(),
        legend.text.align = 0, legend.text = element_text(family="Arial", size=10),
        legend.key.size = unit(0, 'lines'),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-1, size=12),
        axis.title.y = element_text(family="Arial", vjust=1, size=12.5)) }

# Stress plot:
NMDS.data_FL60 = metaMDS(RespVar_FL60, distance = "horn", k = 2, trymax = 100, maxint = 1000); NMDS.data_FL60
  # stress < 0.05 = excellent representation in reduced dimensions
  # stress < 0.1 = great, < 0.2 = good/ok, stress < 0.3 = poor

pdf("Plots/StressPlot_PH-ATF_FL60.pdf")
stressplot(NMDS.data_FL60, M.H_FL60)
dev.off()

# PH & ATF diets according to fork length only (interspecific comparison):
NMDS_FL60 = data.frame(NMDS1_FL60 = NMDS.data_FL60$points[,1], NMDS2_FL60 = NMDS.data_FL60$points[,2], group=as.factor(ExpVar_FL60$Species))

df_FL60 = data.frame()
for(g in NMDS_FL60$group) {
  df_FL60 = rbind(df_FL60, cbind(as.data.frame(with(NMDS_FL60[NMDS_FL60$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_FL60,NMDS2_FL60), wt=rep(1/length(NMDS1_FL60), length(NMDS1_FL60)))$cov, center=c(mean(NMDS1_FL60), mean(NMDS2_FL60))))), group=g)) }

Centroids_FL60 = df_FL60 %>%
  group_by(group) %>%
  mutate(mean.NMDS1_FL60 = mean(NMDS1_FL60)) %>%
  mutate(mean.NMDS2_FL60 = mean(NMDS2_FL60))

# Figure S3 a):
NMDS_FL60 = ggplot(data = NMDS_FL60, aes(x=NMDS1_FL60, y=NMDS2_FL60)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_FL60, aes(x=NMDS1_FL60, y=NMDS2_FL60, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
geom_point(data=Centroids_FL60, aes(x=mean.NMDS1_FL60, y=mean.NMDS2_FL60, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("blue", "red"), name = "Species") +
  scale_fill_manual(values = c("blue", "red"), name = "Species") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_FL60.png", plot=NMDS_FL60, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to fork length and year:
NMDS_FL60_yr = data.frame(NMDS1_FL60_yr = NMDS.data_FL60$points[,1], NMDS2_FL60_yr = NMDS.data_FL60$points[,2], group=as.factor(ExpVar_FL60$Year))

df_FL60_yr = data.frame()
for(g in NMDS_FL60_yr$group) {
  df_FL60_yr = rbind(df_FL60_yr, cbind(as.data.frame(with(NMDS_FL60_yr[NMDS_FL60_yr$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_FL60_yr,NMDS2_FL60_yr), wt=rep(1/length(NMDS1_FL60_yr), length(NMDS1_FL60_yr)))$cov, center=c(mean(NMDS1_FL60_yr), mean(NMDS2_FL60_yr))))), group=g)) }

Centroids_FL60_yr = df_FL60_yr %>%
  group_by(group) %>%
  mutate(mean.NMDS1_FL60_yr = mean(NMDS1_FL60_yr)) %>%
  mutate(mean.NMDS2_FL60_yr = mean(NMDS2_FL60_yr))

# Figure S3 a):
NMDS_FL60_yr = ggplot(data = NMDS_FL60_yr, aes(x=NMDS1_FL60_yr, y=NMDS2_FL60_yr)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_FL60_yr, aes(x=NMDS1_FL60_yr, y=NMDS2_FL60_yr, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
geom_point(data=Centroids_FL60_yr, aes(x=mean.NMDS1_FL60_yr, y=mean.NMDS2_FL60_yr, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("gray20", "gray70"), name = "Year") +
  scale_fill_manual(values = c("gray20", "gray70"), name = "Year") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_FL60_yr.png", plot=NMDS_FL60_yr, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to fork length and month:
NMDS_FL60_mo = data.frame(NMDS1_FL60_mo = NMDS.data_FL60$points[,1], NMDS2_FL60_mo = NMDS.data_FL60$points[,2], group=as.factor(ExpVar_FL60$Month))
levels(NMDS_FL60_mo$group) = c("Jun", "Jul", "Aug")  

  df_FL60_mo = data.frame()
for(g in levels(NMDS_FL60_mo$group)) {
  df_FL60_mo = rbind(df_FL60_mo, cbind(as.data.frame(with(NMDS_FL60_mo[NMDS_FL60_mo$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_FL60_mo,NMDS2_FL60_mo), wt=rep(1/length(NMDS1_FL60_mo), length(NMDS1_FL60_mo)))$cov, center=c(mean(NMDS1_FL60_mo), mean(NMDS2_FL60_mo))))), group=g)) }

Centroids_FL60_mo = df_FL60_mo %>%
  group_by(group) %>%
  mutate(mean.NMDS1_FL60_mo = mean(NMDS1_FL60_mo)) %>%
  mutate(mean.NMDS2_FL60_mo = mean(NMDS2_FL60_mo))

# Figure S3 a):
NMDS_FL60_mo = ggplot(data = NMDS_FL60_mo, aes(x=NMDS1_FL60_mo, y=NMDS2_FL60_mo)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_FL60_mo, aes(x=NMDS1_FL60_mo, y=NMDS2_FL60_mo, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
geom_point(data=Centroids_FL60_mo, aes(x=mean.NMDS1_FL60_mo, y=mean.NMDS2_FL60_mo, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  scale_fill_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.895, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_FL60_mo.png", plot=NMDS_FL60_mo, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to fork length and site:
NMDS_FL60_site = data.frame(NMDS1_FL60_site = NMDS.data_FL60$points[,1], NMDS2_FL60_site = NMDS.data_FL60$points[,2], group=as.factor(as.character(ExpVar_FL60$Location))) # site cannot be an ordered factor

df_FL60_site = data.frame()
for(g in levels(NMDS_FL60_site$group)) {
  df_FL60_site = rbind(df_FL60_site, cbind(as.data.frame(with(NMDS_FL60_site[NMDS_FL60_site$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_FL60_site,NMDS2_FL60_site), wt=rep(1/length(NMDS1_FL60_site), length(NMDS1_FL60_site)))$cov, center=c(mean(NMDS1_FL60_site), mean(NMDS2_FL60_site))))), group=g)) }

Centroids_FL60_site = df_FL60_site %>%
  group_by(group) %>%
  mutate(mean.NMDS1_FL60_site = mean(NMDS1_FL60_site)) %>%
  mutate(mean.NMDS2_FL60_site = mean(NMDS2_FL60_site))

NMDS_FL60_site$group = ordered(NMDS_FL60_site$group, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard"))

# Figure S3 a):
NMDS_FL60_site = ggplot(data = NMDS_FL60_site, aes(x=NMDS1_FL60_site, y=NMDS2_FL60_site)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_FL60_site, aes(x=NMDS1_FL60_site, y=NMDS2_FL60_site, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
geom_point(data=Centroids_FL60_site, aes(x=mean.NMDS1_FL60_site, y=mean.NMDS2_FL60_site, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  scale_fill_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.61, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_FL60_site.png", plot=NMDS_FL60_site, dpi=500, width=4, height=3.25, units="in")

####################################################################
# Multivariate analyses according to gape height (GH), 96-115 mm:
DietData_GH96 = subset(DietData_comp, GH_Bin == "96-115")

# Summarize diet data:
propData_GH96 = DietData_GH96 %>%
  group_by(Species, Year, Month, Location, ID_final_Mod) %>%
  summarise(preyWT = sum(Prey.Mass..g.)) %>%
  mutate(propWT = preyWT / sum (preyWT))

# Transform data for normality:
propData_GH96 = as.data.frame(propData_GH96)
propData_GH96$transWT = propData_GH96$propWT ^ 0.25

# Select specific columns of interest:
propData_GH96 = propData_GH96[,c("Species", "Year", "Month", "Location", "ID_final_Mod", "propWT", "transWT")]

# Set up input variables:
WideData_GH96 = dcast(propData_GH96, Species + Year + Month + Location ~ ID_final_Mod, value.var = "transWT", fun.aggregate = mean, fill=0)
WideData_GH96 = na.omit(WideData_GH96)
WideData_GH96$Month = as.factor(WideData_GH96$Month)
droplevels(WideData_GH96$Location)

ncol(WideData_GH96); colnames(WideData_GH96)[5]; colnames(WideData_GH96)[24]
RespVar_GH96 = WideData_GH96[ ,5:24]
ExpVar_GH96 = WideData_GH96[ ,c("Species", "Year", "Month", "Location")]

# Calculate Morisita-Horn dissimilarity index:  
M.H_GH96 = vegdist(RespVar_GH96, method = "horn") # all data

# Test for differences in multivariate dispersion among factor levels:
disp1_GH96 = betadisper(M.H_GH96, WideData_GH96$Species)
disp2_GH96 = betadisper(M.H_GH96, WideData_GH96$Year)
disp3_GH96 = betadisper(M.H_GH96, WideData_GH96$Month)
disp4_GH96 = betadisper(M.H_GH96, WideData_GH96$Location)

permutest(disp1_GH96) # diff. in disp. between species (ATF > PH)
  # do not include in PERMANOVA
permutest(disp2_GH96) # no diff. disp. between years
permutest(disp3_GH96) # no diff. disp. among months
permutest(disp4_GH96) # no diff. disp. among locations

# Develop individual models with single explanatory variable:
Yr_input_GH96 = adonis2(M.H_GH96 ~ Year, data = WideData_GH96, method = "horn", permutations = 999)
Yr_input_GH96

Mo_input_GH96 = adonis2(M.H_GH96 ~ Month, data = WideData_GH96, method = "horn", permutations = 999)
Mo_input_GH96

Loc_input_GH96 = adonis2(M.H_GH96 ~ Location, data = WideData_GH96, method = "horn", permutations = 999)
Loc_input_GH96

# Full model with explanatory variable order according to pseudo F and p:
Perm_full_GH96 = adonis2(M.H_GH96 ~ Location + Month + Year, data=WideData_GH96, method = "horn", by = "margin", permutations = 9999); Perm_full_GH96

# Drop non-significant main effects:
Perm_mod_GH96 = adonis2(M.H_GH96 ~ Location + Month, data=WideData_GH96, method = "horn", by = "margin", permutations = 9999); Perm_mod_GH96

Perm_final_GH96 = adonis2(M.H_GH96 ~ Location, data=WideData_GH96, method = "horn", by = "margin", permutations = 9999); Perm_final_GH96

# Stress plot:
NMDS.data_GH96 = metaMDS(RespVar_GH96, distance = "horn", k = 2, trymax = 100, maxint = 1000); NMDS.data_GH96
  # stress < 0.05 = excellent representation in reduced dimensions
  # stress < 0.1 = great, < 0.2 = good/ok, stress < 0.3 = poor

pdf("Plots/StressPlot_PH-ATF_GH96.pdf")
stressplot(NMDS.data_GH96, M.H_GH96)
dev.off()

# PH & ATF diets according to gape width only (interspecific comparison):
NMDS_GH96 = data.frame(NMDS1_GH96 = NMDS.data_GH96$points[,1], NMDS2_GH96 = NMDS.data_GH96$points[,2], group=as.factor(ExpVar_GH96$Species))

  df_GH96 = data.frame()
for(g in NMDS_GH96$group) {
  df_GH96 = rbind(df_GH96, cbind(as.data.frame(with(NMDS_GH96[NMDS_GH96$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH96,NMDS2_GH96), wt=rep(1/length(NMDS1_GH96), length(NMDS1_GH96)))$cov, center=c(mean(NMDS1_GH96), mean(NMDS2_GH96))))), group=g)) }

Centroids_GH96 = df_GH96 %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH96 = mean(NMDS1_GH96)) %>%
  mutate(mean.NMDS2_GH96 = mean(NMDS2_GH96))

# Figure S3 b):
NMDS_GH96 = ggplot(data = NMDS_GH96, aes(x=NMDS1_GH96, y=NMDS2_GH96)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH96, aes(x=NMDS1_GH96, y=NMDS2_GH96, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH96, aes(x=mean.NMDS1_GH96, y=mean.NMDS2_GH96, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("blue", "red"), name = "Species") +
  scale_fill_manual(values = c("blue", "red"), name = "Species") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH96.png", plot=NMDS_GH96, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and year:
NMDS_GH96_yr = data.frame(NMDS1_GH96_yr = NMDS.data_GH96$points[,1], NMDS2_GH96_yr = NMDS.data_GH96$points[,2], group=as.factor(ExpVar_GH96$Year))

df_GH96_yr = data.frame()
for(g in NMDS_GH96_yr$group) {
  df_GH96_yr = rbind(df_GH96_yr, cbind(as.data.frame(with(NMDS_GH96_yr[NMDS_GH96_yr$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH96_yr,NMDS2_GH96_yr), wt=rep(1/length(NMDS1_GH96_yr), length(NMDS1_GH96_yr)))$cov, center=c(mean(NMDS1_GH96_yr), mean(NMDS2_GH96_yr))))), group=g)) }

Centroids_GH96_yr = df_GH96_yr %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH96_yr = mean(NMDS1_GH96_yr)) %>%
  mutate(mean.NMDS2_GH96_yr = mean(NMDS2_GH96_yr))

# Figure S3 b):
NMDS_GH96_yr = ggplot(data = NMDS_GH96_yr, aes(x=NMDS1_GH96_yr, y=NMDS2_GH96_yr)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH96_yr, aes(x=NMDS1_GH96_yr, y=NMDS2_GH96_yr, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH96_yr, aes(x=mean.NMDS1_GH96_yr, y=mean.NMDS2_GH96_yr, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("gray20", "gray70"), name = "Year") +
  scale_fill_manual(values = c("gray20", "gray70"), name = "Year") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH96_yr.png", plot=NMDS_GH96_yr, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and month:
NMDS_GH96_mo = data.frame(NMDS1_GH96_mo = NMDS.data_GH96$points[,1], NMDS2_GH96_mo = NMDS.data_GH96$points[,2], group=as.factor(ExpVar_GH96$Month))
levels(NMDS_GH96_mo$group) = c("Jun", "Jul", "Aug")  

df_GH96_mo = data.frame()
for(g in levels(NMDS_GH96_mo$group)) {
  df_GH96_mo = rbind(df_GH96_mo, cbind(as.data.frame(with(NMDS_GH96_mo[NMDS_GH96_mo$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH96_mo,NMDS2_GH96_mo), wt=rep(1/length(NMDS1_GH96_mo), length(NMDS1_GH96_mo)))$cov, center=c(mean(NMDS1_GH96_mo), mean(NMDS2_GH96_mo))))), group=g)) }

Centroids_GH96_mo = df_GH96_mo %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH96_mo = mean(NMDS1_GH96_mo)) %>%
  mutate(mean.NMDS2_GH96_mo = mean(NMDS2_GH96_mo))

# Figure S3 b):
NMDS_GH96_mo = ggplot(data = NMDS_GH96_mo, aes(x=NMDS1_GH96_mo, y=NMDS2_GH96_mo)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH96_mo, aes(x=NMDS1_GH96_mo, y=NMDS2_GH96_mo, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH96_mo, aes(x=mean.NMDS1_GH96_mo, y=mean.NMDS2_GH96_mo, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  scale_fill_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.895, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH96_mo.png", plot=NMDS_GH96_mo, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and site:
NMDS_GH96_site = data.frame(NMDS1_GH96_site = NMDS.data_GH96$points[,1], NMDS2_GH96_site = NMDS.data_GH96$points[,2], group=as.factor(as.character(ExpVar_GH96$Location))) # site cannot be an ordered factor

  df_GH96_site = data.frame()
for(g in levels(NMDS_GH96_site$group)) {
  df_GH96_site = rbind(df_GH96_site, cbind(as.data.frame(with(NMDS_GH96_site[NMDS_GH96_site$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH96_site,NMDS2_GH96_site), wt=rep(1/length(NMDS1_GH96_site), length(NMDS1_GH96_site)))$cov, center=c(mean(NMDS1_GH96_site), mean(NMDS2_GH96_site))))), group=g)) }

Centroids_GH96_site = df_GH96_site %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH96_site = mean(NMDS1_GH96_site)) %>%
  mutate(mean.NMDS2_GH96_site = mean(NMDS2_GH96_site))

NMDS_GH96_site$group = ordered(NMDS_GH96_site$group, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard"))

# Figure S3 b):
NMDS_GH96_site = ggplot(data = NMDS_GH96_site, aes(x=NMDS1_GH96_site, y=NMDS2_GH96_site)) + 
  geom_point (aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH96_site, aes(x=NMDS1_GH96_site, y=NMDS2_GH96_site, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH96_site, aes(x=mean.NMDS1_GH96_site, y=mean.NMDS2_GH96_site, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  scale_fill_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.61, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH96_site.png", plot=NMDS_GH96_site, dpi=500, width=4, height=3.25, units="in")

####################################################################
# Multivariate analyses according to gape height(GH), 116-135 mm:
DietData_GH116 = subset(DietData_comp, GH_Bin == "116-135")

# Summarize diet data:
propData_GH116 = DietData_GH116 %>%
  group_by(Species, Year, Month, Location, ID_final_Mod) %>%
  summarise(preyWT = sum(Prey.Mass..g.)) %>%
  mutate(propWT = preyWT / sum (preyWT))

# Transform data for normality:
propData_GH116 = as.data.frame(propData_GH116)
propData_GH116$transWT = propData_GH116$propWT ^ 0.25

# Select specific columns of interest:
propData_GH116 = propData_GH116[,c("Species", "Year", "Month", "Location", "ID_final_Mod", "propWT", "transWT")]

# Set up input variables:
WideData_GH116 = dcast(propData_GH116, Species + Year + Month + Location ~ ID_final_Mod, value.var = "transWT", fun.aggregate = mean, fill=0)
WideData_GH116 = na.omit(WideData_GH116)
WideData_GH116$Month = as.factor(WideData_GH116$Month)
droplevels(WideData_GH116$Location)

ncol(WideData_GH116); colnames(WideData_GH116)[5]; colnames(WideData_GH116)[20]
RespVar_GH116 = WideData_GH116[ ,5:20]
ExpVar_GH116 = WideData_GH116[ ,c("Species", "Year", "Month", "Location")]

# Calculate Morisita-Horn dissimilarity index:  
M.H_GH116 = vegdist(RespVar_GH116, method = "horn") # all data

# Test for differences in multivariate dispersion among factor levels:
disp1_GH116 = betadisper(M.H_GH116, WideData_GH116$Species)
disp2_GH116 = betadisper(M.H_GH116, WideData_GH116$Year)
disp3_GH116 = betadisper(M.H_GH116, WideData_GH116$Month)
disp4_GH116 = betadisper(M.H_GH116, WideData_GH116$Location)

permutest(disp1_GH116) # no diff. disp. between species
permutest(disp2_GH116) # no diff. disp. between years
permutest(disp3_GH116) # no diff. disp. among months
permutest(disp4_GH116) # no diff. disp. among locations

# Develop individual models with single explanatory variable:
Spp_input_GH116 = adonis2(M.H_GH116 ~ Species, data = WideData_GH116, method = "horn", permutations = 999)
Spp_input_GH116

Yr_input_GH116 = adonis2(M.H_GH116 ~ Year, data = WideData_GH116, method = "horn", permutations = 999)
Yr_input_GH116

Mo_input_GH116 = adonis2(M.H_GH116 ~ Month, data = WideData_GH116, method = "horn", permutations = 999)
Mo_input_GH116

Loc_input_GH116 = adonis2(M.H_GH116 ~ Location, data = WideData_GH116, method = "horn", permutations = 999)
Loc_input_GH116

# Full model with explanatory variable order according to pseudo F and p:
Perm_full_GH116 = adonis2(M.H_GH116 ~ Species + Month + Year + Location, data=WideData_GH116, method = "horn", by = "margin", permutations = 9999); Perm_full_GH116

# Drop non-significant main effects:
Perm_mod_GH116 = adonis2(M.H_GH116 ~ Species + Month + Year, data=WideData_GH116, method = "horn", by = "margin", permutations = 9999); Perm_mod_GH116

Perm_final_GH116 = adonis2(M.H_GH116 ~ Species + Month, data=WideData_GH116, method = "horn", by = "margin", permutations = 9999); Perm_final_GH116

# Stress plot:
NMDS.data_GH116 = metaMDS(RespVar_GH116, distance = "horn", k = 2, trymax = 100, maxint = 1000); NMDS.data_GH116
# stress < 0.05 = excellent representation in reduced dimensions
# stress < 0.1 = great, < 0.2 = good/ok, stress < 0.3 = poor

pdf("Plots/StressPlot_PH-ATF_GH116.pdf")
stressplot(NMDS.data_GH116, M.H_GH116)
dev.off()

# PH & ATF diets according to gape width only (interspecific comparison):
NMDS_GH116 = data.frame(NMDS1_GH116 = NMDS.data_GH116$points[,1], NMDS2_GH116 = NMDS.data_GH116$points[,2], group=as.factor(ExpVar_GH116$Species))

df_GH116 = data.frame()
for(g in NMDS_GH116$group) {
  df_GH116 = rbind(df_GH116, cbind(as.data.frame(with(NMDS_GH116[NMDS_GH116$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH116,NMDS2_GH116), wt=rep(1/length(NMDS1_GH116), length(NMDS1_GH116)))$cov, center=c(mean(NMDS1_GH116), mean(NMDS2_GH116))))), group=g)) }

Centroids_GH116 = df_GH116 %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH116 = mean(NMDS1_GH116)) %>%
  mutate(mean.NMDS2_GH116 = mean(NMDS2_GH116))

# Figure S3 c):
NMDS_GH116 = ggplot(data = NMDS_GH116, aes(x=NMDS1_GH116, y=NMDS2_GH116)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH116, aes(x=NMDS1_GH116, y=NMDS2_GH116, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH116, aes(x=mean.NMDS1_GH116, y=mean.NMDS2_GH116, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("blue", "red"), name = "Species") +
  scale_fill_manual(values = c("blue", "red"), name = "Species") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH116116.png", plot=NMDS_GH116, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and year:
NMDS_GH116_yr = data.frame(NMDS1_GH116_yr = NMDS.data_GH116$points[,1], NMDS2_GH116_yr = NMDS.data_GH116$points[,2], group=as.factor(ExpVar_GH116$Year))

df_GH116_yr = data.frame()
for(g in NMDS_GH116_yr$group) {
  df_GH116_yr = rbind(df_GH116_yr, cbind(as.data.frame(with(NMDS_GH116_yr[NMDS_GH116_yr$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH116_yr,NMDS2_GH116_yr), wt=rep(1/length(NMDS1_GH116_yr), length(NMDS1_GH116_yr)))$cov, center=c(mean(NMDS1_GH116_yr), mean(NMDS2_GH116_yr))))), group=g)) }

Centroids_GH116_yr = df_GH116_yr %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH116_yr = mean(NMDS1_GH116_yr)) %>%
  mutate(mean.NMDS2_GH116_yr = mean(NMDS2_GH116_yr))

# Figure S3 c):
NMDS_GH116_yr = ggplot(data = NMDS_GH116_yr, aes(x=NMDS1_GH116_yr, y=NMDS2_GH116_yr)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH116_yr, aes(x=NMDS1_GH116_yr, y=NMDS2_GH116_yr, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH116_yr, aes(x=mean.NMDS1_GH116_yr, y=mean.NMDS2_GH116_yr, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("gray20", "gray70"), name = "Year") +
  scale_fill_manual(values = c("gray20", "gray70"), name = "Year") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH116_yr.png", plot=NMDS_GH116_yr, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and month:
NMDS_GH116_mo = data.frame(NMDS1_GH116_mo = NMDS.data_GH116$points[,1], NMDS2_GH116_mo = NMDS.data_GH116$points[,2], group=as.factor(ExpVar_GH116$Month))
levels(NMDS_GH116_mo$group) = c("Jun", "Jul", "Aug")  

df_GH116_mo = data.frame()
for(g in levels(NMDS_GH116_mo$group)) {
  df_GH116_mo = rbind(df_GH116_mo, cbind(as.data.frame(with(NMDS_GH116_mo[NMDS_GH116_mo$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH116_mo,NMDS2_GH116_mo), wt=rep(1/length(NMDS1_GH116_mo), length(NMDS1_GH116_mo)))$cov, center=c(mean(NMDS1_GH116_mo), mean(NMDS2_GH116_mo))))), group=g)) }

Centroids_GH116_mo = df_GH116_mo %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH116_mo = mean(NMDS1_GH116_mo)) %>%
  mutate(mean.NMDS2_GH116_mo = mean(NMDS2_GH116_mo))

# Figure S3 c):
NMDS_GH116_mo = ggplot(data = NMDS_GH116_mo, aes(x=NMDS1_GH116_mo, y=NMDS2_GH116_mo)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH116_mo, aes(x=NMDS1_GH116_mo, y=NMDS2_GH116_mo, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH116_mo, aes(x=mean.NMDS1_GH116_mo, y=mean.NMDS2_GH116_mo, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  scale_fill_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.895, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH116_mo.png", plot=NMDS_GH116_mo, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and site:
NMDS_GH116_site = data.frame(NMDS1_GH116_site = NMDS.data_GH116$points[,1], NMDS2_GH116_site = NMDS.data_GH116$points[,2], group=as.factor(as.character(ExpVar_GH116$Location))) # site cannot be an ordered factor

df_GH116_site = data.frame()
for(g in levels(NMDS_GH116_site$group)) {
  df_GH116_site = rbind(df_GH116_site, cbind(as.data.frame(with(NMDS_GH116_site[NMDS_GH116_site$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GH116_site,NMDS2_GH116_site), wt=rep(1/length(NMDS1_GH116_site), length(NMDS1_GH116_site)))$cov, center=c(mean(NMDS1_GH116_site), mean(NMDS2_GH116_site))))), group=g)) }

Centroids_GH116_site = df_GH116_site %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GH116_site = mean(NMDS1_GH116_site)) %>%
  mutate(mean.NMDS2_GH116_site = mean(NMDS2_GH116_site))

# Figure S3 c):
NMDS_GH116_site$group = ordered(NMDS_GH116_site$group, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard"))
NMDS_GH116_site = ggplot(data = NMDS_GH116_site, aes(x=NMDS1_GH116_site, y=NMDS2_GH116_site)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GH116_site, aes(x=NMDS1_GH116_site, y=NMDS2_GH116_site, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GH116_site, aes(x=mean.NMDS1_GH116_site, y=mean.NMDS2_GH116_site, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  scale_fill_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.61, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GH116_site.png", plot=NMDS_GH116_site, dpi=500, width=4, height=3.25, units="in")

####################################################################
# Multivariate analyses according to gape width (GW), 96-115 mm:
DietData_GW96 = subset(DietData_comp, GW_Bin == "96-115")

# Summarize diet data:
propData_GW96 = DietData_GW96 %>%
  group_by(Species, Year, Month, Location, ID_final_Mod) %>%
  summarise(preyWT = sum(Prey.Mass..g.)) %>%
  mutate(propWT = preyWT / sum (preyWT))

# Transform data for normality:
propData_GW96 = as.data.frame(propData_GW96)
propData_GW96$transWT = propData_GW96$propWT ^ 0.25

# Select specific columns of interest:
propData_GW96 = propData_GW96[,c("Species", "Year", "Month", "Location", "ID_final_Mod", "propWT", "transWT")]

# Set up input variables:
WideData_GW96 = dcast(propData_GW96, Species + Year + Month + Location ~ ID_final_Mod, value.var = "transWT", fun.aggregate = mean, fill=0)
WideData_GW96 = na.omit(WideData_GW96)
WideData_GW96$Month = as.factor(WideData_GW96$Month)
droplevels(WideData_GW96$Location)

ncol(WideData_GW96); colnames(WideData_GW96)[5]; colnames(WideData_GW96)[24]
RespVar_GW96 = WideData_GW96[ ,5:24]
ExpVar_GW96 = WideData_GW96[ ,c("Species", "Year", "Month", "Location")]

# Calculate Morisita-Horn dissimilarity index:  
M.H_GW96 = vegdist(RespVar_GW96, method = "horn") # all data

# Test for differences in multivariate dispersion among factor levels:
disp1_GW96 = betadisper(M.H_GW96, WideData_GW96$Species)
disp2_GW96 = betadisper(M.H_GW96, WideData_GW96$Year)
disp3_GW96 = betadisper(M.H_GW96, WideData_GW96$Month)
disp4_GW96 = betadisper(M.H_GW96, WideData_GW96$Location)

permutest(disp1_GW96) # diff. in disp. between species (ATF > PH)
# do not include in PERMANOVA
permutest(disp2_GW96) # no diff. disp. between years
permutest(disp3_GW96) # no diff. disp. among months
permutest(disp4_GW96) # no diff. disp. among locations

# Develop individual models with single explanatory variable:
Yr_input_GW96 = adonis2(M.H_GW96 ~ Year, data = WideData_GW96, method = "horn", permutations = 999)
Yr_input_GW96

Mo_input_GW96 = adonis2(M.H_GW96 ~ Month, data = WideData_GW96, method = "horn", permutations = 999)
Mo_input_GW96

Loc_input_GW96 = adonis2(M.H_GW96 ~ Location, data = WideData_GW96, method = "horn", permutations = 999)
Loc_input_GW96

# Full model with explanatory variable order according to pseudo F and p:
Perm_full_GW96 = adonis2(M.H_GW96 ~ Month + Location + Year, data=WideData_GW96, method = "horn", by = "margin", permutations = 9999); Perm_full_GW96

# Drop non-significant main effects:
Perm_mod_GW96 = adonis2(M.H_GW96 ~ Month + Location, data=WideData_GW96, method = "horn", by = "margin", permutations = 9999); Perm_mod_GW96

Perm_final_GW96 = adonis2(M.H_GW96 ~ Month, data=WideData_GW96, method = "horn", by = "margin", permutations = 9999); Perm_final_GW96

# Stress plot:
NMDS.data_GW96 = metaMDS(RespVar_GW96, distance = "horn", k = 2, trymax = 100, maxint = 1000); NMDS.data_GW96
# stress < 0.05 = excellent representation in reduced dimensions
# stress < 0.1 = great, < 0.2 = good/ok, stress < 0.3 = poor

pdf("Plots/StressPlot_PH-ATF_GW96.pdf")
stressplot(NMDS.data_GW96, M.H_GW96)
dev.off()

# PH & ATF diets according to gape width only (interspecific comparison):
NMDS_GW96 = data.frame(NMDS1_GW96 = NMDS.data_GW96$points[,1], NMDS2_GW96 = NMDS.data_GW96$points[,2], group=as.factor(ExpVar_GW96$Species))

df_GW96 = data.frame()
for(g in NMDS_GW96$group) {
  df_GW96 = rbind(df_GW96, cbind(as.data.frame(with(NMDS_GW96[NMDS_GW96$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW96,NMDS2_GW96), wt=rep(1/length(NMDS1_GW96), length(NMDS1_GW96)))$cov, center=c(mean(NMDS1_GW96), mean(NMDS2_GW96))))), group=g)) }

Centroids_GW96 = df_GW96 %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW96 = mean(NMDS1_GW96)) %>%
  mutate(mean.NMDS2_GW96 = mean(NMDS2_GW96))

# Figure S3 d):
NMDS_GW96 = ggplot(data = NMDS_GW96, aes(x=NMDS1_GW96, y=NMDS2_GW96)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW96, aes(x=NMDS1_GW96, y=NMDS2_GW96, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW96, aes(x=mean.NMDS1_GW96, y=mean.NMDS2_GW96, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("blue", "red"), name = "Species") +
  scale_fill_manual(values = c("blue", "red"), name = "Species") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW96.png", plot=NMDS_GW96, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and year:
NMDS_GW96_yr = data.frame(NMDS1_GW96_yr = NMDS.data_GW96$points[,1], NMDS2_GW96_yr = NMDS.data_GW96$points[,2], group=as.factor(ExpVar_GW96$Year))

df_GW96_yr = data.frame()
for(g in NMDS_GW96_yr$group) {
  df_GW96_yr = rbind(df_GW96_yr, cbind(as.data.frame(with(NMDS_GW96_yr[NMDS_GW96_yr$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW96_yr,NMDS2_GW96_yr), wt=rep(1/length(NMDS1_GW96_yr), length(NMDS1_GW96_yr)))$cov, center=c(mean(NMDS1_GW96_yr), mean(NMDS2_GW96_yr))))), group=g)) }

Centroids_GW96_yr = df_GW96_yr %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW96_yr = mean(NMDS1_GW96_yr)) %>%
  mutate(mean.NMDS2_GW96_yr = mean(NMDS2_GW96_yr))

# Figure S3 d):
NMDS_GW96_yr = ggplot(data = NMDS_GW96_yr, aes(x=NMDS1_GW96_yr, y=NMDS2_GW96_yr)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW96_yr, aes(x=NMDS1_GW96_yr, y=NMDS2_GW96_yr, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW96_yr, aes(x=mean.NMDS1_GW96_yr, y=mean.NMDS2_GW96_yr, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("gray20", "gray70"), name = "Year") +
  scale_fill_manual(values = c("gray20", "gray70"), name = "Year") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW96_yr.png", plot=NMDS_GW96_yr, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and month:
NMDS_GW96_mo = data.frame(NMDS1_GW96_mo = NMDS.data_GW96$points[,1], NMDS2_GW96_mo = NMDS.data_GW96$points[,2], group=as.factor(ExpVar_GW96$Month))
levels(NMDS_GW96_mo$group) = c("Jun", "Jul", "Aug")  

df_GW96_mo = data.frame()
for(g in levels(NMDS_GW96_mo$group)) {
  df_GW96_mo = rbind(df_GW96_mo, cbind(as.data.frame(with(NMDS_GW96_mo[NMDS_GW96_mo$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW96_mo,NMDS2_GW96_mo), wt=rep(1/length(NMDS1_GW96_mo), length(NMDS1_GW96_mo)))$cov, center=c(mean(NMDS1_GW96_mo), mean(NMDS2_GW96_mo))))), group=g)) }

Centroids_GW96_mo = df_GW96_mo %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW96_mo = mean(NMDS1_GW96_mo)) %>%
  mutate(mean.NMDS2_GW96_mo = mean(NMDS2_GW96_mo))

# Figure S3 d):
NMDS_GW96_mo = ggplot(data = NMDS_GW96_mo, aes(x=NMDS1_GW96_mo, y=NMDS2_GW96_mo)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW96_mo, aes(x=NMDS1_GW96_mo, y=NMDS2_GW96_mo, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW96_mo, aes(x=mean.NMDS1_GW96_mo, y=mean.NMDS2_GW96_mo, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  scale_fill_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.895, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW96_mo.png", plot=NMDS_GW96_mo, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and site:
NMDS_GW96_site = data.frame(NMDS1_GW96_site = NMDS.data_GW96$points[,1], NMDS2_GW96_site = NMDS.data_GW96$points[,2], group=as.factor(as.character(ExpVar_GW96$Location))) # site cannot be an ordered factor

df_GW96_site = data.frame()
for(g in levels(NMDS_GW96_site$group)) {
  df_GW96_site = rbind(df_GW96_site, cbind(as.data.frame(with(NMDS_GW96_site[NMDS_GW96_site$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW96_site,NMDS2_GW96_site), wt=rep(1/length(NMDS1_GW96_site), length(NMDS1_GW96_site)))$cov, center=c(mean(NMDS1_GW96_site), mean(NMDS2_GW96_site))))), group=g)) }

Centroids_GW96_site = df_GW96_site %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW96_site = mean(NMDS1_GW96_site)) %>%
  mutate(mean.NMDS2_GW96_site = mean(NMDS2_GW96_site))

# Figure S3 d):
NMDS_GW96_site$group = ordered(NMDS_GW96_site$group, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard"))
NMDS_GW96_site = ggplot(data = NMDS_GW96_site, aes(x=NMDS1_GW96_site, y=NMDS2_GW96_site)) + 
  geom_point (aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW96_site, aes(x=NMDS1_GW96_site, y=NMDS2_GW96_site, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW96_site, aes(x=mean.NMDS1_GW96_site, y=mean.NMDS2_GW96_site, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  scale_fill_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.61, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW96_site.png", plot=NMDS_GW96_site, dpi=500, width=4, height=3.25, units="in")

####################################################################
# Multivariate analyses according to gape width (GW), 116-135 mm:
DietData_GW116 = subset(DietData_comp, GW_Bin == "116-135")

# Summarize diet data:
propData_GW116 = DietData_GW116 %>%
  group_by(Species, Year, Month, Location, ID_final_Mod) %>%
  summarise(preyWT = sum(Prey.Mass..g.)) %>%
  mutate(propWT = preyWT / sum (preyWT))

# Transform data for normality:
propData_GW116 = as.data.frame(propData_GW116)
propData_GW116$transWT = propData_GW116$propWT ^ 0.25

# Select specific columns of interest:
propData_GW116 = propData_GW116[,c("Species", "Year", "Month", "Location", "ID_final_Mod", "propWT", "transWT")]

# Set up input variables:
WideData_GW116 = dcast(propData_GW116, Species + Year + Month + Location ~ ID_final_Mod, value.var = "transWT", fun.aggregate = mean, fill=0)
WideData_GW116 = na.omit(WideData_GW116)
WideData_GW116$Month = as.factor(WideData_GW116$Month)
droplevels(WideData_GW116$Location)

ncol(WideData_GW116); colnames(WideData_GW116)[5]; colnames(WideData_GW116)[20]
RespVar_GW116 = WideData_GW116[ ,5:20]
ExpVar_GW116 = WideData_GW116[ ,c("Species", "Year", "Month", "Location")]

# Calculate Morisita-Horn dissimilarity index:  
M.H_GW116 = vegdist(RespVar_GW116, method = "horn") # all data

# Test for differences in multivariate dispersion among factor levels:
disp1_GW116 = betadisper(M.H_GW116, WideData_GW116$Species)
disp2_GW116 = betadisper(M.H_GW116, WideData_GW116$Year)
disp3_GW116 = betadisper(M.H_GW116, WideData_GW116$Month)
disp4_GW116 = betadisper(M.H_GW116, WideData_GW116$Location)

permutest(disp1_GW116) # no diff. disp. between species
permutest(disp2_GW116) # no diff. disp. between years
permutest(disp3_GW116) # no diff. disp. among months
permutest(disp4_GW116) # no diff. disp. among locations

# Develop individual models with single explanatory variable:
Spp_input_GW116 = adonis2(M.H_GW116 ~ Species, data = WideData_GW116, method = "horn", permutations = 999)
Spp_input_GW116

Yr_input_GW116 = adonis2(M.H_GW116 ~ Year, data = WideData_GW116, method = "horn", permutations = 999)
Yr_input_GW116

Mo_input_GW116 = adonis2(M.H_GW116 ~ Month, data = WideData_GW116, method = "horn", permutations = 999)
Mo_input_GW116

Loc_input_GW116 = adonis2(M.H_GW116 ~ Location, data = WideData_GW116, method = "horn", permutations = 999)
Loc_input_GW116

# Full model with explanatory variable order according to pseudo F and p:
Perm_full_GW116 = adonis2(M.H_GW116 ~ Species + Month + Location + Year, data=WideData_GW116, method = "horn", by = "margin", permutations = 9999); Perm_full_GW116

# Drop non-significant main effects:
Perm_mod_GW116 = adonis2(M.H_GW116 ~ Species + Month + Location, data=WideData_GW116, method = "horn", by = "margin", permutations = 9999); Perm_mod_GW116

Perm_final_GW116 = adonis2(M.H_GW116 ~ Species + Month, data=WideData_GW116, method = "horn", by = "margin", permutations = 9999); Perm_final_GW116

# Stress plot:
NMDS.data_GW116 = metaMDS(RespVar_GW116, distance = "horn", k = 2, trymax = 100, maxint = 1000); NMDS.data_GW116
  # stress < 0.05 = excellent representation in reduced dimensions
  # stress < 0.1 = great, < 0.2 = good/ok, stress < 0.3 = poor

pdf("Plots/StressPlot_PH-ATF_GW116.pdf")
stressplot(NMDS.data_GW116, M.H_GW116)
dev.off()

# PH & ATF diets according to gape width only (interspecific comparison):
NMDS_GW116 = data.frame(NMDS1_GW116 = NMDS.data_GW116$points[,1], NMDS2_GW116 = NMDS.data_GW116$points[,2], group=as.factor(ExpVar_GW116$Species))

df_GW116 = data.frame()
for(g in NMDS_GW116$group) {
  df_GW116 = rbind(df_GW116, cbind(as.data.frame(with(NMDS_GW116[NMDS_GW116$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW116,NMDS2_GW116), wt=rep(1/length(NMDS1_GW116), length(NMDS1_GW116)))$cov, center=c(mean(NMDS1_GW116), mean(NMDS2_GW116))))), group=g)) }

Centroids_GW116 = df_GW116 %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW116 = mean(NMDS1_GW116)) %>%
  mutate(mean.NMDS2_GW116 = mean(NMDS2_GW116))

# Figure S3 e):
NMDS_GW116 = ggplot(data = NMDS_GW116, aes(x=NMDS1_GW116, y=NMDS2_GW116)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW116, aes(x=NMDS1_GW116, y=NMDS2_GW116, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW116, aes(x=mean.NMDS1_GW116, y=mean.NMDS2_GW116, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("blue", "red"), name = "Species") +
  scale_fill_manual(values = c("blue", "red"), name = "Species") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW116116.png", plot=NMDS_GW116, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and year:
NMDS_GW116_yr = data.frame(NMDS1_GW116_yr = NMDS.data_GW116$points[,1], NMDS2_GW116_yr = NMDS.data_GW116$points[,2], group=as.factor(ExpVar_GW116$Year))

df_GW116_yr = data.frame()
for(g in NMDS_GW116_yr$group) {
  df_GW116_yr = rbind(df_GW116_yr, cbind(as.data.frame(with(NMDS_GW116_yr[NMDS_GW116_yr$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW116_yr,NMDS2_GW116_yr), wt=rep(1/length(NMDS1_GW116_yr), length(NMDS1_GW116_yr)))$cov, center=c(mean(NMDS1_GW116_yr), mean(NMDS2_GW116_yr))))), group=g)) }

Centroids_GW116_yr = df_GW116_yr %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW116_yr = mean(NMDS1_GW116_yr)) %>%
  mutate(mean.NMDS2_GW116_yr = mean(NMDS2_GW116_yr))

# Figure S3 e):
NMDS_GW116_yr = ggplot(data = NMDS_GW116_yr, aes(x=NMDS1_GW116_yr, y=NMDS2_GW116_yr)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW116_yr, aes(x=NMDS1_GW116_yr, y=NMDS2_GW116_yr, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW116_yr, aes(x=mean.NMDS1_GW116_yr, y=mean.NMDS2_GW116_yr, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("gray20", "gray70"), name = "Year") +
  scale_fill_manual(values = c("gray20", "gray70"), name = "Year") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.885, 0.932)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW116_yr.png", plot=NMDS_GW116_yr, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and month:
NMDS_GW116_mo = data.frame(NMDS1_GW116_mo = NMDS.data_GW116$points[,1], NMDS2_GW116_mo = NMDS.data_GW116$points[,2], group=as.factor(ExpVar_GW116$Month))
levels(NMDS_GW116_mo$group) = c("Jun", "Jul", "Aug")  

df_GW116_mo = data.frame()
for(g in levels(NMDS_GW116_mo$group)) {
  df_GW116_mo = rbind(df_GW116_mo, cbind(as.data.frame(with(NMDS_GW116_mo[NMDS_GW116_mo$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW116_mo,NMDS2_GW116_mo), wt=rep(1/length(NMDS1_GW116_mo), length(NMDS1_GW116_mo)))$cov, center=c(mean(NMDS1_GW116_mo), mean(NMDS2_GW116_mo))))), group=g)) }

Centroids_GW116_mo = df_GW116_mo %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW116_mo = mean(NMDS1_GW116_mo)) %>%
  mutate(mean.NMDS2_GW116_mo = mean(NMDS2_GW116_mo))

# Figure S3 e):
NMDS_GW116_mo = ggplot(data = NMDS_GW116_mo, aes(x=NMDS1_GW116_mo, y=NMDS2_GW116_mo)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW116_mo, aes(x=NMDS1_GW116_mo, y=NMDS2_GW116_mo, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW116_mo, aes(x=mean.NMDS1_GW116_mo, y=mean.NMDS2_GW116_mo, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  scale_fill_manual(values = c("greenyellow", "turquoise1", "lightseagreen"), name = "Month") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.895, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW116_mo.png", plot=NMDS_GW116_mo, dpi=500, width=4, height=3.25, units="in")

# PH & ATF diets according to gape width and site:
NMDS_GW116_site = data.frame(NMDS1_GW116_site = NMDS.data_GW116$points[,1], NMDS2_GW116_site = NMDS.data_GW116$points[,2], group=as.factor(as.character(ExpVar_GW116$Location))) # site cannot be an ordered factor

df_GW116_site = data.frame()
for(g in levels(NMDS_GW116_site$group)) {
  df_GW116_site = rbind(df_GW116_site, cbind(as.data.frame(with(NMDS_GW116_site[NMDS_GW116_site$group==g,], veganCovEllipse(cov.wt(cbind(NMDS1_GW116_site,NMDS2_GW116_site), wt=rep(1/length(NMDS1_GW116_site), length(NMDS1_GW116_site)))$cov, center=c(mean(NMDS1_GW116_site), mean(NMDS2_GW116_site))))), group=g)) }

Centroids_GW116_site = df_GW116_site %>%
  group_by(group) %>%
  mutate(mean.NMDS1_GW116_site = mean(NMDS1_GW116_site)) %>%
  mutate(mean.NMDS2_GW116_site = mean(NMDS2_GW116_site))

# Figure S3 e):
NMDS_GW116_site$group = ordered(NMDS_GW116_site$group, levels = c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard"))
NMDS_GW116_site = ggplot(data = NMDS_GW116_site, aes(x=NMDS1_GW116_site, y=NMDS2_GW116_site)) + 
  geom_point(aes(color = group, fill = group), size=0.5) +
  geom_polygon(data=df_GW116_site, aes(x=NMDS1_GW116_site, y=NMDS2_GW116_site, fill=group, colour=group), alpha=0.05, linetype=1, lwd = 0.25, show.legend = F) +
  geom_point(data=Centroids_GW116_site, aes(x=mean.NMDS1_GW116_site, y=mean.NMDS2_GW116_site, color = group), shape = 3, size = 2, show.legend = F) +
  scale_color_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  scale_fill_manual(values = c("darkorange", "forestgreen", "darkorchid3"), name = "Location") +
  plot.theme() +
  coord_fixed() +
  ggtitle("") +
  labs(x="NMDS 1", y="NMDS 2") +
  guides(fill=guide_legend(nrow=c(4))) +
  theme(legend.position = c(0.61, 0.91)) +
  scale_x_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0","")) +
  scale_y_continuous(limits = c(-2.03, 2.03), breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5), labels=c("","-1.0","","0.0","","1.0",""))
ggsave(filename="Plots/NMDS_PH_ATF_GW116_site.png", plot=NMDS_GW116_site, dpi=500, width=4, height=3.25, units="in")