# The code below uses quantile regression to assess predator-prey size spectra for Pacific Halibut (PH) and Arrowtooth Flounder (ATF) caught in Southeast Alaska (2015 and 2016). Specimens were collected using hook-and-line near Juneau, AK in 2015. See Barnes et al. (2021) for details.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(ggplot2)
require(quantreg)

#############################################################
# Import and prepare diet data:
PreyData = read.csv("Data/PH_ATF_DietData_SEAK.csv")

# Remove spring and fall months with very few samples (focus = summer):
PreyData$Year = as.factor(PreyData$Year)
PreyData = subset(PreyData, Month != 5 & Month != 9)

# Consolidate individual capture locations into sampling sites:
PreyData$Location = with(PreyData, 
                        ifelse(Grid..  < 1, "Chatham.Strait",
                        ifelse(Grid.. <= 3, "Lynn.Canal",
                        ifelse(Grid.. <= 7, "Favorite-Saginaw.Channel",
                        ifelse(Grid.. <= 8, "Pt.Howard",
                        ifelse(Grid.. <= 9, "Funter.Bay",
                        ifelse(Grid.. <= 10, "Other", 
                        ifelse(Grid.. <= 11, "Pt.Couverden", 
                        ifelse(Grid.. <= 13, "Icy.Strait", "Other")))))))))

# Remove locations outside of spatial extent, order others north to south, and rename:
PreyData = subset(PreyData, Location != "Other" & Location != "Chatham.Strait")
PreyData$Location = ordered(PreyData$Location, levels=c("Lynn.Canal", "Favorite-Saginaw.Channel", "Pt.Howard", "Funter.Bay", "Pt.Couverden", "Icy.Strait"))

# Remove unmeasured predators:
PreyData = subset(PreyData, Length..cm. != "NA")

# Assign size bins based on length:
PreyData$FL = round(PreyData$Length..cm., 0)
PreyData$FL_Bin = cut(PreyData$FL, breaks = c(0, 29, 39, 49, 59, 69, 79, 89, 99, 109, 119, 129, 139, 149, 159), include.lowest = F)
levels(PreyData$FL_Bin) = c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100-109", "110-119", "120-129", "130-139", "140-149", "150-159")

# Select only non-empty stomachs for cumulative prey curves:
WithContents = subset(PreyData, Prey.Mass..g. > 0)

# Predict gape from length (see '01_AllometricGrowth' for linear relationships):
PH.data = subset(WithContents, Species == "PH")
  PH.data$GH = round((PH.data$Length..cm. * 1.11497), 0)
  PH.data$GW = round((PH.data$Length..cm. * 1.19663), 0)
ATF.data = subset(WithContents, Species == "ATF")  
  ATF.data$GH = round((ATF.data$Length..cm. * 2.14768), 0)
  ATF.data$GW = round((ATF.data$Length..cm. * 2.06384), 0)
All.data = rbind(PH.data, ATF.data) # recombine

# Assign size bins based on gape:
All.data$GH_Bin = cut(All.data$GH, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = F)
levels(All.data$GH_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

All.data$GW_Bin = cut(All.data$GW, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = F)
levels(All.data$GW_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

# Order prey items based upon phylogeny and exclude "other" category:
All.data$ID_final_Mod = as.factor(All.data$ID_final_Mod)
All.data$ID_final_Mod = ordered(All.data$ID_final_Mod, levels = c("Mollusca", "Teuthida", "Octopodidae", "Crustacea", "Pandalidae", "Lithodidae", "Majoidea", "Chionoecetes bairdi", "Metacarcinus magister", "Paguroidea", "Chondrichthyes", "Teleostei", "Clupea pallasii", "Salmoniformes", "Gadus", "Gadus chalcogrammus", "Gadus macrocephalus", "Scorpaeniformes", "Pleuronectiformes", "Benthic Material", "Other"))
All.data = subset(All.data, ID_final_Mod!="Other" & ID_final_Mod!="Benthic Material")

#############################################################
# Prey sizes: standard lengths for fishes:

# Calculate standard length (mm) for herring without heads:
herring = subset(All.data, ID_final_Mod == "Clupea pallasii")
herring = herring[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "Total.L..mm.", "Std.L..mm.", "Head.L..mm.", "Std.L..mm....No.Head")]

# sample sizes (direct measurements):
nrow(subset(herring, Total.L..mm. > 0))
nrow(subset(herring, Std.L..mm. > 0))

# Quantify relationship between standard length and head length (to predict length from prey without heads):
H.SL_herr = lm(herring$Head.L..mm. ~ herring$Std.L..mm. - 1) 
  summary(H.SL_herr)

# Subtract head length from std length (with head):
herring$Std.L..mm.NoHead = herring$Std.L..mm. - (herring$Std.L..mm. * H.SL_herr$coefficients[1])
  herring$diff = herring$Std.L..mm. - herring$Std.L..mm.NoHead
NH.SL_herr = lm(herring$diff ~ herring$Std.L..mm.NoHead - 1) 
  summary(NH.SL_herr)

# Predict standard length from herring prey without heads (many in this condition):
herring$Std.L..mm....WithHead = herring$Std.L..mm....No.Head + (herring$Std.L..mm....No.Head * NH.SL_herr$coefficients[1])

# Estimate std length from total length:
SL.TL_herr = lm(herring$Std.L..mm. ~ herring$Total.L..mm. - 1)
  summary(SL.TL_herr)

# Group best estimates of herring std. length (mm):
herring$SL = with(herring, 
  ifelse(Total.L..mm. > 0 & is.na(Std.L..mm.), 
      (Total.L..mm. * SL.TL_herr$coefficients[1]), Std.L..mm.))
herring$SL = with(herring, 
  ifelse(is.na(SL), Std.L..mm....WithHead, SL))


# Estimate pollock std length (mm) from otolith length (mm):
pollock = subset(All.data, ID_final_Mod == "Gadus chalcogrammus")
pollock = pollock[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "Total.L..mm.", "Std.L..mm.", "Head.L..mm.", "Std.L..mm....No.Head", "Oto.1.L..mm.", "Oto.2.L..mm.", "Oto.Cond...1", "Oto.Cond...2")]

# Estimate std length from total length:
SL.TL_poll = lm(pollock$Std.L..mm. ~ pollock$Total.L..mm. - 1)
  summary(SL.TL_poll)

# Select the otolith with the greater length (i.e., less digested):
pollock$OtoL = pmax(pollock$Oto.1.L..mm., pollock$Oto.2.L..mm., na.rm = T)

# Quantify relationship between standard length and otolith length:
poll.oto = subset(pollock, Oto.Cond...1 != "fair" | Oto.Cond...2 != "fair")
poll.oto = subset(poll.oto, Oto.Cond...1 != "poor" | Oto.Cond...2 != "poor")
poll.oto = subset(poll.oto, Oto.Cond...1 != "poor")

SL.OL_poll = lm(poll.oto$Std.L..mm. ~ poll.oto$OtoL - 1)
  summary(SL.OL_poll)

pollock$SL = with(pollock,
  ifelse(Total.L..mm. > 0 & is.na(Std.L..mm.), 
        (Total.L..mm. * SL.TL_poll$coefficients[1]), Std.L..mm.))
pollock$SL = with(pollock,
  ifelse(is.na(SL), (OtoL * SL.OL_poll$coefficients[1]), SL))


# Estimate P. cod std length (mm) from otolith length (mm):
cod = subset(All.data, ID_final_Mod == "Gadus macrocephalus")
cod = cod[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "Total.L..mm.", "Std.L..mm.", "Head.L..mm.", "Std.L..mm....No.Head", "Oto.1.L..mm.", "Oto.2.L..mm.", "Oto.Cond...1", "Oto.Cond...2")]

# Estimate std length from total length:
SL_TL_cod = lm(cod$Std.L..mm. ~ cod$Total.L..mm. - 1)
  summary(SL_TL_cod)

# Select otolith with greater length:
cod$OtoL = pmax(cod$Oto.1.L..mm., cod$Oto.2.L..mm., na.rm = T)

# Quantify the relationship between standard length and otolith length:
cod.oto = subset(cod, Oto.Cond...1 != "fair" | Oto.Cond...2 != "fair")
cod.oto = subset(cod.oto, Oto.Cond...1 != "poor" | Oto.Cond...2 != "poor")
cod.oto = subset(cod.oto, Oto.Cond...1 != "poor")

SL_OL_cod = lm(cod.oto$Std.L..mm. ~ cod.oto$OtoL - 1)
  summary(SL_OL_cod)

cod$SL = with(cod,
   ifelse(Total.L..mm. > 0 & is.na(Std.L..mm.), 
         (Total.L..mm. * SL_TL_cod$coefficients[1]), Std.L..mm.))
cod$SL = with(cod,
   ifelse(is.na(SL), (OtoL * SL_OL_cod$coefficients[1]), SL))


# Estimate std length (mm) for all other fishes:
fish = subset(All.data, 
                ID_final_Mod == "Teleostei" | 
                ID_final_Mod == "Chondrichthyes" | 
                ID_final_Mod == "Salmoniformes" | 
                ID_final_Mod == "Gadus" | 
                ID_final_Mod == "Scorpaeniformes" | 
                ID_final_Mod == "Pleuronectiformes")

fish = fish[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "Std.L..mm.")]

fish = subset(fish, Std.L..mm. > 0)
fish$SL = fish$Std.L..mm.

#############################################################
# Prey sizes: body mass (kg) for cephalopods:

# Octopuses 'std length' = body mass (kg):
  # Reference: Buckley et al. (unpub. data) https://www.afsc.noaa.gov/Quarterly/ond2011/divrptsREFM1.htm
octo = subset(All.data, ID_final_Mod == "Octopodidae")
octo = octo[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "LHL")]
  octo$SL = 0.0038 * (octo$LHL ^ 3.0419)


# Squids 'std length' = body mass (kg), estimated from dorsal mantle length (mm), estimated from URL (mm)
  # Reference: Sinclair et al. (2015; assuming B. magister) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4545836/)
squid = subset(All.data, ID_final_Mod == "Teuthida")
squid = squid[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "URL")]
  squid$DML = (45.47 * squid$URL) - 0.72
  squid$SL = (0.00008101 * (squid$DML ^ 2.816)) / 1000

#############################################################
# Prey sizes: carapace width (mm) for crabs:
# Crabs 'std length' = carapace width (mm), measured directly:
crabs = subset(All.data, 
              ID_final_Mod == "Lithodidae" | 
              ID_final_Mod == "Majoidea" | 
              ID_final_Mod == "Chionoecetes bairdi" | 
              ID_final_Mod == "Metacarcinus magister" | 
              ID_final_Mod == "Paguroidea")
crabs = crabs[, c("Species", "PreyData.PredatorID", "Length..mm.", "GH", "GW", "PreyID", "ID_final_Mod", "Carape.W..mm.")]
crabs = subset(crabs, Carape.W..mm. > 0)
crabs$SL = crabs$Carape.W..mm.

#############################################################
# Combine all prey size data:
PreySizes = herring %>% 
  full_join(., pollock) %>% 
  full_join(., cod) %>% 
  full_join(., fish) %>% 
  full_join(., octo) %>% 
  full_join(., squid) %>% 
  full_join(., crabs)

# Combine gadids for analyses:
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Gadus chalcogrammus"] = "Gadus"
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Gadus macrocephalus"] = "Gadus"
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Salmoniformes"] = "Teleostei"

# Combine all true crabs:
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Majoidea"] = "Brachyura"
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Chionoecetes bairdi"] = "Brachyura"
levels(PreySizes$ID_final_Mod)[levels(PreySizes$ID_final_Mod)=="Metacarcinus magister"] = "Brachyura"

# Group all fishes:
PreySizes$TaxaType = with(PreySizes,
    ifelse(ID_final_Mod == "Clupea pallasii", " fishes",
    ifelse(ID_final_Mod == "Teleostei", " fishes",
    ifelse(ID_final_Mod == "Salmoniformes", " fishes",
    ifelse(ID_final_Mod == "Scorpaeniformes", " fishes",
    ifelse(ID_final_Mod == "Gadus", " fishes",
    ifelse(ID_final_Mod == "Octopodidae", " cephalopods",
    ifelse(ID_final_Mod == "Teuthida", " cephalopods",
    ifelse(ID_final_Mod == "Lithodidae", " crabs",
    ifelse(ID_final_Mod == "Brachyura", " crabs",
    ifelse(ID_final_Mod == "Paguroidea", " crabs", NA)))))))))))

# Extract only those with direct or predicted measurements:
PreySizes = subset(PreySizes, SL > 0)

PreySizes$TaxaType = as.factor(PreySizes$TaxaType)
  PreySizes$TaxaType = ordered(PreySizes$TaxaType, levels = c(" fishes", " crabs", " cephalopods"))

#############################################################
# Estimate quantile regression parameters for Pacific Halibut (according to fork length, gape height, and gape width):

# Fish prey (Table 5):
FishOnly = subset(PreySizes, TaxaType == " fishes")
FishOnly = subset(FishOnly, ID_final_Mod != "Teleostei")
  PH_fishes = subset(FishOnly, Species == "PH")
rqPH_fish_FL = rq(Length..mm. ~ SL, data = PH_fishes, tau = c(0.10, 0.90))
  summary(rqPH_fish_FL, se = "ker")
rqPH_fish_GH = rq(GH ~ SL, data = PH_fishes, tau = c(0.10, 0.90))
  summary(rqPH_fish_GH, se = "ker")
rqPH_fish_GW = rq(GW ~ SL, data = PH_fishes, tau = c(0.10, 0.90))
  summary(rqPH_fish_GW, se = "ker")
  
# Fish prey size by predator gape width (PH only; Fig. 7 [top]):
SizeSpectra.fishes = ggplot(PH_fishes, aes(x = GW, y = SL, col = ID_final_Mod)) + 
  geom_point(size=1.0) +
  scale_color_manual(values = c("Clupea pallasii"="cadetblue1", "Gadus"="blue", "Scorpaeniformes"="purple2"), labels = c(expression(paste(italic(" Clupea pallasii"))), expression(paste(italic(" Gadus spp."))), " Scorpaeniformes")) +
  geom_quantile(quantiles = c(0.10, 0.90), lwd = 0.75, color = "blue", show.legend = F) +
  theme_bw() +
  expand_limits(x=0, y=0) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_blank(),
        legend.key.size = unit(0, 'lines'),
        legend.key.height = unit(4, "mm"),
        legend.position = c(0.28, 0.895),
        legend.text = element_text(family="Arial", size=9),
        legend.text.align = 0,
        strip.text = element_blank()) +
  labs(x="Predator Gape Width (mm)", y="Fish Standard Length (mm)") +
  scale_x_continuous(limits = c(0,205)) +
  scale_y_continuous(limits = c(0,650))
ggsave(SizeSpectra.fishes, filename="Plots/SizeSpectra_PH_fishes.png", width=3, height=3, units="in")


# Cephalopod prey (Table 5):
CephOnly = subset(PreySizes, TaxaType == " cephalopods")
  PH_ceph = subset(CephOnly, Species == "PH")
rqPH_ceph_FL = rq(Length..mm. ~ SL, data = PH_ceph, tau = c(0.10, 0.90))
  summary(rqPH_ceph_FL, se = "ker")
rqPH_ceph_GH = rq(GH ~ SL, data = PH_ceph, tau = c(0.10, 0.90))
  summary(rqPH_ceph_GH, se = "ker")
rqPH_ceph_GW = rq(GW ~ SL, data = PH_ceph, tau = c(0.10, 0.90))
  summary(rqPH_ceph_GW, se = "ker")

# Cephalopod prey size by predator gape width (PH only; Fig. 7 [center]):
SizeSpectra.cephalopods = ggplot(CephOnly, aes(x = GW, y = SL, col = ID_final_Mod)) + 
  geom_point(size=1.0) +
  scale_color_manual(values = c("Teuthida"="brown4", "Octopodidae"="firebrick3"), labels = c(" Teuthida", " Octopodidae")) +
# geom_quantile(quantiles = c(0.10, 0.90), lwd = 0.75, show.legend = F, colour = "brown3") +
  theme_bw() +
  expand_limits(x=0, y=0) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_blank(),
        legend.key.size = unit(0, 'lines'),
        legend.key.height = unit(4, "mm"),
        legend.position = c(0.22, 0.93),
        legend.text = element_text(family="Arial", size=9),
        legend.text.align = 0,
        strip.text = element_blank()) +
  labs(x="Predator Gape Width (mm)", y="Cephalopod Body Mass (kg)") +
  scale_x_continuous(limits = c(0,205)) +
  scale_y_continuous(limits = c(0,3.5))
ggsave(SizeSpectra.cephalopods, filename="Plots/SizeSpectra_PH_ceph.png", width=3, height=3, units="in")


# Crab prey (Table 5):
CrabsOnly = subset(PreySizes, TaxaType == " crabs")
  PH_crabs = subset(CrabsOnly, Species == "PH")
rqPH_crabs_FL = rq(Length..mm. ~ SL, data = PH_crabs, tau = c(0.10, 0.90))
  summary(rqPH_crabs_FL, se = "ker")
rqPH_crabs_GH = rq(GH ~ SL, data = PH_crabs, tau = c(0.10, 0.90))
  summary(rqPH_crabs_GH, se = "ker")
rqPH_crabs_GW = rq(GW ~ SL, data = PH_crabs, tau = c(0.10, 0.90))
  summary(rqPH_crabs_GW, se = "ker")

# Crab prey size by predator gape width (PH only; Fig. 7 [center]):
SizeSpectra.crabs = ggplot(CrabsOnly, aes(x = GW, y = SL, col = ID_final_Mod)) + 
  geom_point(size=1.0) +
  scale_color_manual(values = c("Lithodidae"="chocolate", "Brachyura"="orange", "Paguroidea"="yellow2"), labels = c(" Lithodidae", " Brachyura", " Paguroidea")) +
  geom_quantile(quantiles = c(0.10, 0.90), lwd = 0.75, show.legend = F, colour = "orange") +
  theme_bw() +
  expand_limits(x=0, y=0) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_blank(),
        legend.key.size = unit(0, 'lines'),
        legend.key.height = unit(4, "mm"),
        legend.position = c(0.21 , 0.895),
        legend.text = element_text(family="Arial", size=9),
        legend.text.align = 0,
        strip.text = element_blank()) +
  labs(x="Predator Gape Width (mm)", y="Crab Carapace Width (mm)") +
  scale_x_continuous(limits = c(0,205)) +
  scale_y_continuous(limits = c(0,75))
ggsave(SizeSpectra.crabs, filename="Plots/SizeSpectra_PH_crabs.png", width=3, height=3, units="in")