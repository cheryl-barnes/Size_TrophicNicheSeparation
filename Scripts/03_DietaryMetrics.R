# The code below estimates a variety of dietary metrics (e.g., prey richness, diversity, evenness) for Pacific Halibut (PH) and Arrowtooth Flounder (ATF). Specimens were collected using hook-and-line near Juneau, AK in 2015 and 2016. See Barnes et al. (2021) for details.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(ggplot2)
require(vegan)
require(tidyr)
require(reshape2)

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

# Identify range of lengths sampled (for results text):
DietData %>% 
  group_by(Species) %>%
  summarize(range(Length..cm.))

# Calculate sample sizes:
Field.N = unique(DietData[ ,c("Year", "Month", "Location", "Type", "Species", "PredatorData.PredatorID", "FL_Bin")])
  table(Field.N$Species)
  table(Field.N$Year, Field.N$Species)
  table(Field.N$Location, Field.N$Species)
  table(Field.N$Type, Field.N$Species)
  table(Field.N$Month, Field.N$Year, Field.N$Species)
  table(Field.N$FL_Bin, Field.N$Species)

# Order prey items based upon phylogeny:
DietData$ID_final_Mod = ordered(DietData$ID_final_Mod, levels = c("Mollusca", "Teuthida", "Octopodidae", "Crustacea", "Pandalidae", "Lithodidae", "Majoidea", "Chionoecetes bairdi", "Metacarcinus magister", "Paguroidea", "Chondrichthyes", "Teleostei", "Clupea pallasii", "Salmoniformes", "Gadus", "Gadus chalcogrammus", "Gadus macrocephalus", "Scorpaeniformes", "Pleuronectiformes", "Benthic Material", "Other"))

# Predict gape from length (see '01_AllometricGrowth' for linear relationships):
PH.data = subset(DietData, Species == "PH")
  PH.data$GH = round((PH.data$Length..cm. * 1.11497), 0)
  PH.data$GW = round((PH.data$Length..cm. * 1.19663), 0)
ATF.data = subset(DietData, Species == "ATF")  
  ATF.data$GH = round((ATF.data$Length..cm. * 2.14768), 0)
  ATF.data$GW = round((ATF.data$Length..cm. * 2.06384), 0)
All.data = rbind(PH.data, ATF.data) # recombine

# Assign size bins based on gape:
All.data$GH_Bin = cut(All.data$GH, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = F)
levels(All.data$GH_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

All.data$GW_Bin = cut(All.data$GW, breaks = c(0, 35, 55, 75, 95, 115, 135, 155, 175, 195, 215), include.lowest = F)
levels(All.data$GW_Bin) = c("<36", "36-55", "56-75", "76-95", "96-115", "116-135", "136-155", "156-175", "176-195", ">195")

# Calculate sample sizes by size metric and class:
Diet.N = unique(All.data[ , c("Year", "Month", "Location", "Species", "PredatorData.PredatorID", "FL_Bin", "GH_Bin", "GW_Bin")])
table(Diet.N$FL_Bin, Diet.N$Species)
table(Diet.N$GH_Bin, Diet.N$Species)
table(Diet.N$GW_Bin, Diet.N$Species)
table(Diet.N$Location, Diet.N$GW_Bin, Diet.N$Species)

# Remove sites with insufficient samples for ATF:
All.data_sites = subset(All.data, Location != "Funter.Bay" & Location != "Pt.Couverden" & Location != "Icy.Strait")

# Order predators:
All.data_sites$Species = as.factor(All.data_sites$Species)
###################################################################
# Summary dietary metrics:

##########
# Trophic level (TL) - by predator (all sizes, comparable sites, and time periods combined):
# Proportions by weight (to be multiplied to prey trophic levels):
propWT = All.data_sites %>%
  group_by(Species, ID_final) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT = as.data.frame(propWT)

# Assign trophic levels to prey >= 1% total weight (Aydin et al. 2007; to be multiplied to proportions by weight):
propWT$TrophLevel = with(propWT,
ifelse(ID_final == "Acantholithodes hispidus", 3.1, 
ifelse(ID_final == "Amphipoda", 2.5,
ifelse(ID_final == "Argis dentata", 2.5,
ifelse(ID_final == "Atheresthes stomias", 3.9,
ifelse(ID_final == "Berryteuthis magister", 3.7,
ifelse(ID_final == "Bivalvia", 2.5,
ifelse(ID_final == "Bothidae", 3.5,
ifelse(ID_final == "Caridea", 2.9,
ifelse(ID_final == "Cephalopoda", 3.75,
ifelse(ID_final == "Chionoecetes", 3.4,
ifelse(ID_final == "Chionoecetes bairdi", 3.4,   
ifelse(ID_final == "Chondrichthyes", 4.47,
ifelse(ID_final == "Chorilia longipes", 3.1,
ifelse(ID_final == "Clupea pallasii", 3.5, 
ifelse(ID_final == "Cnidaria", 3.5,
ifelse(ID_final == "Cottidae", 3.9,
ifelse(ID_final == "Crangon", 2.5,
ifelse(ID_final == "Craniata", 3.65,
ifelse(ID_final == "Derbesia marina", 1.0,
ifelse(ID_final == "Diogenidae", 3.1,
ifelse(ID_final == "Echinoidea", 2.0,
ifelse(ID_final == "Elassochirus cavimanus", 3.1,
ifelse(ID_final == "Enteroctopus dofleini", 3.8,
ifelse(ID_final == "Gadus", 3.65,
ifelse(ID_final == "Gadus chalcogrammus", 3.6,
ifelse(ID_final == "Gadus macrocephalus", 3.7,
ifelse(ID_final == "Gallus gallus", 1.0,
ifelse(ID_final == "Gastropoda", 2.9,
ifelse(ID_final == "Heptacarpus brevirostris", 3.1,
ifelse(ID_final == "Homo sapiens", 0.0,
ifelse(ID_final == "Hyas lyratus", 3.1,
ifelse(ID_final == "Isopoda", 2.5,
ifelse(ID_final == "Lithodes aequispinus", 3.4,
ifelse(ID_final == "Lithodidae", 3.4,
ifelse(ID_final == "Lopholithodes foraminatus", 3.4,
ifelse(ID_final == "Lycodes", 3.6,
ifelse(ID_final == "Melonchela clathriata", 2.5,
ifelse(ID_final == "Metacarcinus magister", 3.1,  
ifelse(ID_final == "Mollusca", 2.7,
ifelse(ID_final == "Mycale loveni", 2.5,
ifelse(ID_final == "Mytilidae", 2.5,
ifelse(ID_final == "Ochrophyta", 1.0,
ifelse(ID_final == "Octopodidae", 3.8,
ifelse(ID_final == "Odontopyxis trispinosa", 3.8,
ifelse(ID_final == "Oncorhynchus gorbuscha", 3.8,
ifelse(ID_final == "Oncorhynchus keta", 3.8,
ifelse(ID_final == "Ophiuroidea", 2.2, NA
))))))))))))))))))))))))))))))))))))))))))))))))

propWT$TrophLevel = with(propWT,       
ifelse(ID_final == "Paguridae", 3.1,
ifelse(ID_final == "Paguroidea", 3.1, 
ifelse(ID_final == "Pandalidae", 2.9,
ifelse(ID_final == "Pandalopsis dispar", 2.9,
ifelse(ID_final == "Pandalus", 2.9,
ifelse(ID_final == "Pandalus borealis", 2.9,
ifelse(ID_final == "Pandalus goniurus", 2.9,
ifelse(ID_final == "Pandalus platyceros", 2.9,
ifelse(ID_final == "Paralithodes camtschaticus", 3.4,
ifelse(ID_final == "Pasiphaea pacifica", 2.5,
ifelse(ID_final == "Pasiphaeidea", 2.5,
ifelse(ID_final == "Phyllospadix", 1.0,
ifelse(ID_final == "Plantae", 1.0,
ifelse(ID_final == "Pleocyemata", 2.5,
ifelse(ID_final == "Pleuronectiformes", 3.5,
ifelse(ID_final == "Polychaeta", 2.5,
ifelse(ID_final == "Porifera", 2.5,
ifelse(ID_final == "Rajidae", 4.47,
ifelse(ID_final == "Rhinolithodes wosnessenskii", 3.1,
ifelse(ID_final == "Rhodophyta", 1.0,
ifelse(ID_final == "Rock", 0.0,
ifelse(ID_final == "Rock and Shell Hash", 0.0,
ifelse(ID_final == "Rossia pacifica", 3.7,
ifelse(ID_final == "Salmonidae", 3.8,
ifelse(ID_final == "Scorpaeniformes", 3.73,
ifelse(ID_final == "Sebastidae", 3.8,
ifelse(ID_final == "Soranthera ulvoidea", 1.0,
ifelse(ID_final == "Stichaeus punctatus", 3.5,
ifelse(ID_final == "Teleostei", 3.65,
ifelse(ID_final == "Teuthida", 3.7,
ifelse(ID_final == "Unidentified Inorganic Material", 0.0,
ifelse(ID_final == "Unidentified Organic Matter", 0.0, TrophLevel)))))))))))))))))))))))))))))))))

# Predator trophic level:
propWT$TL = with(propWT, TrophLevel * propWT)

propWT %>%
  group_by(Species) %>%
  summarise(sum(TL)+1)

#######################################
# Trophic level (TL) - by predator and size class according to length (all comparable sites, months, and years combined); Table 2:
# Proportions by weight (to be multiplied to prey trophic levels):
propWT_FL = All.data_sites %>%
  group_by(Species, FL_Bin, ID_final) %>%   
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT_FL = as.data.frame(propWT_FL)

# Assign trophic levels to prey >= 1% total weight (Aydin et al. 2007; to be multiplied to proportions by weight):
propWT_FL$TrophLevel = with(propWT_FL,
ifelse(ID_final == "Acantholithodes hispidus", 3.1, 
ifelse(ID_final == "Amphipoda", 2.5,
ifelse(ID_final == "Argis dentata", 2.5,
ifelse(ID_final == "Atheresthes stomias", 3.9,
ifelse(ID_final == "Berryteuthis magister", 3.7,
ifelse(ID_final == "Bivalvia", 2.5,
ifelse(ID_final == "Bothidae", 3.5,
ifelse(ID_final == "Caridea", 2.9,
ifelse(ID_final == "Cephalopoda", 3.75,
ifelse(ID_final == "Chionoecetes", 3.4,
ifelse(ID_final == "Chionoecetes bairdi", 3.4,   
ifelse(ID_final == "Chondrichthyes", 4.47,
ifelse(ID_final == "Chorilia longipes", 3.1,
ifelse(ID_final == "Clupea pallasii", 3.5, 
ifelse(ID_final == "Cnidaria", 3.5,
ifelse(ID_final == "Cottidae", 3.9,
ifelse(ID_final == "Crangon", 2.5,
ifelse(ID_final == "Craniata", 3.65,
ifelse(ID_final == "Derbesia marina", 1.0,
ifelse(ID_final == "Diogenidae", 3.1,
ifelse(ID_final == "Echinoidea", 2.0,
ifelse(ID_final == "Elassochirus cavimanus", 3.1,
ifelse(ID_final == "Enteroctopus dofleini", 3.8,
ifelse(ID_final == "Gadus", 3.65,
ifelse(ID_final == "Gadus chalcogrammus", 3.6,
ifelse(ID_final == "Gadus macrocephalus", 3.7,
ifelse(ID_final == "Gallus gallus", 1.0,
ifelse(ID_final == "Gastropoda", 2.9,
ifelse(ID_final == "Heptacarpus brevirostris", 3.1,
ifelse(ID_final == "Homo sapiens", 0.0,
ifelse(ID_final == "Hyas lyratus", 3.1,
ifelse(ID_final == "Isopoda", 2.5,
ifelse(ID_final == "Lithodes aequispinus", 3.4,
ifelse(ID_final == "Lithodidae", 3.4,
ifelse(ID_final == "Lopholithodes foraminatus", 3.4,
ifelse(ID_final == "Lycodes", 3.6,
ifelse(ID_final == "Melonchela clathriata", 2.5,
ifelse(ID_final == "Metacarcinus magister", 3.1,  
ifelse(ID_final == "Mollusca", 2.7,
ifelse(ID_final == "Mycale loveni", 2.5,
ifelse(ID_final == "Mytilidae", 2.5,
ifelse(ID_final == "Ochrophyta", 1.0,
ifelse(ID_final == "Octopodidae", 3.8,
ifelse(ID_final == "Odontopyxis trispinosa", 3.8,
ifelse(ID_final == "Oncorhynchus gorbuscha", 3.8,
ifelse(ID_final == "Oncorhynchus keta", 3.8,
ifelse(ID_final == "Ophiuroidea", 2.2, NA
))))))))))))))))))))))))))))))))))))))))))))))))
propWT_FL$TrophLevel = with(propWT_FL,       
ifelse(ID_final == "Paguridae", 3.1,
ifelse(ID_final == "Paguroidea", 3.1, 
ifelse(ID_final == "Pandalidae", 2.9,
ifelse(ID_final == "Pandalopsis dispar", 2.9,
ifelse(ID_final == "Pandalus", 2.9,
ifelse(ID_final == "Pandalus borealis", 2.9,
ifelse(ID_final == "Pandalus goniurus", 2.9,
ifelse(ID_final == "Pandalus platyceros", 2.9,
ifelse(ID_final == "Paralithodes camtschaticus", 3.4,
ifelse(ID_final == "Pasiphaea pacifica", 2.5,
ifelse(ID_final == "Pasiphaeidea", 2.5,
ifelse(ID_final == "Phyllospadix", 1.0,
ifelse(ID_final == "Plantae", 1.0,
ifelse(ID_final == "Pleocyemata", 2.5,
ifelse(ID_final == "Pleuronectiformes", 3.5,
ifelse(ID_final == "Polychaeta", 2.5,
ifelse(ID_final == "Porifera", 2.5,
ifelse(ID_final == "Rajidae", 4.47,
ifelse(ID_final == "Rhinolithodes wosnessenskii", 3.1,
ifelse(ID_final == "Rhodophyta", 1.0,
ifelse(ID_final == "Rock", 0.0,
ifelse(ID_final == "Rock and Shell Hash", 0.0,
ifelse(ID_final == "Rossia pacifica", 3.7,
ifelse(ID_final == "Salmonidae", 3.8,
ifelse(ID_final == "Scorpaeniformes", 3.73,
ifelse(ID_final == "Sebastidae", 3.8,
ifelse(ID_final == "Soranthera ulvoidea", 1.0,
ifelse(ID_final == "Stichaeus punctatus", 3.5,
ifelse(ID_final == "Teleostei", 3.65,
ifelse(ID_final == "Teuthida", 3.7,
ifelse(ID_final == "Unidentified Inorganic Material", 0.0,
ifelse(ID_final == "Unidentified Organic Matter", 0.0, TrophLevel)))))))))))))))))))))))))))))))))

propWT_FL$TL = propWT_FL$TrophLevel * propWT_FL$propWT
TL_FL = na.omit(propWT_FL)

# Keep only groups with sufficient sample sizes (n >= 20):
TL_FL = subset(TL_FL, 
              Species == "ATF" & FL_Bin == "50-59" | 
              Species == "ATF" & FL_Bin == "60-69" | 
              Species == "PH" & FL_Bin == "60-69" | 
              Species == "PH" & FL_Bin == "70-79" | 
              Species == "PH" & FL_Bin == "80-89" | 
              Species == "PH" & FL_Bin == "90-99" | 
              Species == "PH" & FL_Bin == "100-109")

# Predator trophic level by fork length size class:  
TL_FL %>%
  group_by(Species, FL_Bin) %>%
  summarise(sum(TL)+1) 
#######################################
# Trophic level (TL) - by predator and size class according to gape height (all comparable sites, months, and years combined); Table 3:
# Proportions by weight (to be multiplied to prey trophic levels):
propWT_GH = All.data_sites %>%
  group_by(Species, GH_Bin, ID_final) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT_GH = as.data.frame(propWT_GH)

# Assign trophic levels to prey >= 1% total weight (Aydin et al. 2007; to be multiplied to proportions by weight):
propWT_GH$TrophLevel = with(propWT_GH,
ifelse(ID_final == "Acantholithodes hispidus", 3.1, 
ifelse(ID_final == "Amphipoda", 2.5,
ifelse(ID_final == "Argis dentata", 2.5,
ifelse(ID_final == "Atheresthes stomias", 3.9,
ifelse(ID_final == "Berryteuthis magister", 3.7,
ifelse(ID_final == "Bivalvia", 2.5,
ifelse(ID_final == "Bothidae", 3.5,
ifelse(ID_final == "Caridea", 2.9,
ifelse(ID_final == "Cephalopoda", 3.75,
ifelse(ID_final == "Chionoecetes", 3.4,
ifelse(ID_final == "Chionoecetes bairdi", 3.4,   
ifelse(ID_final == "Chondrichthyes", 4.47,
ifelse(ID_final == "Chorilia longipes", 3.1,
ifelse(ID_final == "Clupea pallasii", 3.5, 
ifelse(ID_final == "Cnidaria", 3.5,
ifelse(ID_final == "Cottidae", 3.9,
ifelse(ID_final == "Crangon", 2.5,
ifelse(ID_final == "Craniata", 3.65,
ifelse(ID_final == "Derbesia marina", 1.0,
ifelse(ID_final == "Diogenidae", 3.1,
ifelse(ID_final == "Echinoidea", 2.0,
ifelse(ID_final == "Elassochirus cavimanus", 3.1,
ifelse(ID_final == "Enteroctopus dofleini", 3.8,
ifelse(ID_final == "Gadus", 3.65,
ifelse(ID_final == "Gadus chalcogrammus", 3.6,
ifelse(ID_final == "Gadus macrocephalus", 3.7,
ifelse(ID_final == "Gallus gallus", 1.0,
ifelse(ID_final == "Gastropoda", 2.9,
ifelse(ID_final == "Heptacarpus brevirostris", 3.1,
ifelse(ID_final == "Homo sapiens", 0.0,
ifelse(ID_final == "Hyas lyratus", 3.1,
ifelse(ID_final == "Isopoda", 2.5,
ifelse(ID_final == "Lithodes aequispinus", 3.4,
ifelse(ID_final == "Lithodidae", 3.4,
ifelse(ID_final == "Lopholithodes foraminatus", 3.4,
ifelse(ID_final == "Lycodes", 3.6,
ifelse(ID_final == "Melonchela clathriata", 2.5,
ifelse(ID_final == "Metacarcinus magister", 3.1,  
ifelse(ID_final == "Mollusca", 2.7,
ifelse(ID_final == "Mycale loveni", 2.5,
ifelse(ID_final == "Mytilidae", 2.5,
ifelse(ID_final == "Ochrophyta", 1.0,
ifelse(ID_final == "Octopodidae", 3.8,
ifelse(ID_final == "Odontopyxis trispinosa", 3.8,
ifelse(ID_final == "Oncorhynchus gorbuscha", 3.8,
ifelse(ID_final == "Oncorhynchus keta", 3.8,
ifelse(ID_final == "Ophiuroidea", 2.2, NA
))))))))))))))))))))))))))))))))))))))))))))))))
propWT_GH$TrophLevel = with(propWT_GH,       
ifelse(ID_final == "Paguridae", 3.1,
ifelse(ID_final == "Paguroidea", 3.1, 
ifelse(ID_final == "Pandalidae", 2.9,
ifelse(ID_final == "Pandalopsis dispar", 2.9,
ifelse(ID_final == "Pandalus", 2.9,
ifelse(ID_final == "Pandalus borealis", 2.9,
ifelse(ID_final == "Pandalus goniurus", 2.9,
ifelse(ID_final == "Pandalus platyceros", 2.9,
ifelse(ID_final == "Paralithodes camtschaticus", 3.4,
ifelse(ID_final == "Pasiphaea pacifica", 2.5,
ifelse(ID_final == "Pasiphaeidea", 2.5,
ifelse(ID_final == "Phyllospadix", 1.0,
ifelse(ID_final == "Plantae", 1.0,
ifelse(ID_final == "Pleocyemata", 2.5,
ifelse(ID_final == "Pleuronectiformes", 3.5,
ifelse(ID_final == "Polychaeta", 2.5,
ifelse(ID_final == "Porifera", 2.5,
ifelse(ID_final == "Rajidae", 4.47,
ifelse(ID_final == "Rhinolithodes wosnessenskii", 3.1,
ifelse(ID_final == "Rhodophyta", 1.0,
ifelse(ID_final == "Rock", 0.0,
ifelse(ID_final == "Rock and Shell Hash", 0.0,
ifelse(ID_final == "Rossia pacifica", 3.7,
ifelse(ID_final == "Salmonidae", 3.8,
ifelse(ID_final == "Scorpaeniformes", 3.73,
ifelse(ID_final == "Sebastidae", 3.8,
ifelse(ID_final == "Soranthera ulvoidea", 1.0,
ifelse(ID_final == "Stichaeus punctatus", 3.5,
ifelse(ID_final == "Teleostei", 3.65,
ifelse(ID_final == "Teuthida", 3.7,
ifelse(ID_final == "Unidentified Inorganic Material", 0.0,
ifelse(ID_final == "Unidentified Organic Matter", 0.0, TrophLevel)))))))))))))))))))))))))))))))))

propWT_GH$TL = propWT_GH$TrophLevel * propWT_GH$propWT
TL_GH = na.omit(propWT_GH)

# Keep only groups with sufficient sample sizes (n >= 20):
TL_GH = subset(TL_GH, 
                 Species == "ATF" & GH_Bin == "96-115" | 
                 Species == "ATF" & GH_Bin == "116-135" | 
                 Species == "PH" & GH_Bin == "56-75" | 
                 Species == "PH" & GH_Bin == "76-95" | 
                 Species == "PH" & GH_Bin == "96-115" | 
                 Species == "PH" & GH_Bin == "116-135")

# Predator trophic level by gape height size class:
TL_GH %>%
  group_by(Species, GH_Bin) %>%
  summarise(sum(TL)+1) 
#######################################
# Trophic level (TL) - by predator and size class according to gape width (all comparable sites, months, and years combined); Table 3:
# Proportions by weight (to be multiplied to prey trophic levels):
propWT_GW = All.data_sites %>%
  group_by(Species, GW_Bin, ID_final) %>%   
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT_GW = as.data.frame(propWT_GW)

# Assign trophic levels to prey >= 1% total weight (Aydin et al. 2007; to be multiplied to proportions by weight):
propWT_GW$TrophLevel = with(propWT_GW,
ifelse(ID_final == "Acantholithodes hispidus", 3.1, 
ifelse(ID_final == "Amphipoda", 2.5,
ifelse(ID_final == "Argis dentata", 2.5,
ifelse(ID_final == "Atheresthes stomias", 3.9,
ifelse(ID_final == "Berryteuthis magister", 3.7,
ifelse(ID_final == "Bivalvia", 2.5,
ifelse(ID_final == "Bothidae", 3.5,
ifelse(ID_final == "Caridea", 2.9,
ifelse(ID_final == "Cephalopoda", 3.75,
ifelse(ID_final == "Chionoecetes", 3.4,
ifelse(ID_final == "Chionoecetes bairdi", 3.4,   
ifelse(ID_final == "Chondrichthyes", 4.47,
ifelse(ID_final == "Chorilia longipes", 3.1,
ifelse(ID_final == "Clupea pallasii", 3.5, 
ifelse(ID_final == "Cnidaria", 3.5,
ifelse(ID_final == "Cottidae", 3.9,
ifelse(ID_final == "Crangon", 2.5,
ifelse(ID_final == "Craniata", 3.65,
ifelse(ID_final == "Derbesia marina", 1.0,
ifelse(ID_final == "Diogenidae", 3.1,
ifelse(ID_final == "Echinoidea", 2.0,
ifelse(ID_final == "Elassochirus cavimanus", 3.1,
ifelse(ID_final == "Enteroctopus dofleini", 3.8,
ifelse(ID_final == "Gadus", 3.65,
ifelse(ID_final == "Gadus chalcogrammus", 3.6,
ifelse(ID_final == "Gadus macrocephalus", 3.7,
ifelse(ID_final == "Gallus gallus", 1.0,
ifelse(ID_final == "Gastropoda", 2.9,
ifelse(ID_final == "Heptacarpus brevirostris", 3.1,
ifelse(ID_final == "Homo sapiens", 0.0,
ifelse(ID_final == "Hyas lyratus", 3.1,
ifelse(ID_final == "Isopoda", 2.5,
ifelse(ID_final == "Lithodes aequispinus", 3.4,
ifelse(ID_final == "Lithodidae", 3.4,
ifelse(ID_final == "Lopholithodes foraminatus", 3.4,
ifelse(ID_final == "Lycodes", 3.6,
ifelse(ID_final == "Melonchela clathriata", 2.5,
ifelse(ID_final == "Metacarcinus magister", 3.1,  
ifelse(ID_final == "Mollusca", 2.7,
ifelse(ID_final == "Mycale loveni", 2.5,
ifelse(ID_final == "Mytilidae", 2.5,
ifelse(ID_final == "Ochrophyta", 1.0,
ifelse(ID_final == "Octopodidae", 3.8,
ifelse(ID_final == "Odontopyxis trispinosa", 3.8,
ifelse(ID_final == "Oncorhynchus gorbuscha", 3.8,
ifelse(ID_final == "Oncorhynchus keta", 3.8,
ifelse(ID_final == "Ophiuroidea", 2.2, NA
))))))))))))))))))))))))))))))))))))))))))))))))
propWT_GW$TrophLevel = with(propWT_GW,       
ifelse(ID_final == "Paguridae", 3.1,
ifelse(ID_final == "Paguroidea", 3.1, 
ifelse(ID_final == "Pandalidae", 2.9,
ifelse(ID_final == "Pandalopsis dispar", 2.9,
ifelse(ID_final == "Pandalus", 2.9,
ifelse(ID_final == "Pandalus borealis", 2.9,
ifelse(ID_final == "Pandalus goniurus", 2.9,
ifelse(ID_final == "Pandalus platyceros", 2.9,
ifelse(ID_final == "Paralithodes camtschaticus", 3.4,
ifelse(ID_final == "Pasiphaea pacifica", 2.5,
ifelse(ID_final == "Pasiphaeidea", 2.5,
ifelse(ID_final == "Phyllospadix", 1.0,
ifelse(ID_final == "Plantae", 1.0,
ifelse(ID_final == "Pleocyemata", 2.5,
ifelse(ID_final == "Pleuronectiformes", 3.5,
ifelse(ID_final == "Polychaeta", 2.5,
ifelse(ID_final == "Porifera", 2.5,
ifelse(ID_final == "Rajidae", 4.47,
ifelse(ID_final == "Rhinolithodes wosnessenskii", 3.1,
ifelse(ID_final == "Rhodophyta", 1.0,
ifelse(ID_final == "Rock", 0.0,
ifelse(ID_final == "Rock and Shell Hash", 0.0,
ifelse(ID_final == "Rossia pacifica", 3.7,
ifelse(ID_final == "Salmonidae", 3.8,
ifelse(ID_final == "Scorpaeniformes", 3.73,
ifelse(ID_final == "Sebastidae", 3.8,
ifelse(ID_final == "Soranthera ulvoidea", 1.0,
ifelse(ID_final == "Stichaeus punctatus", 3.5,
ifelse(ID_final == "Teleostei", 3.65,
ifelse(ID_final == "Teuthida", 3.7,
ifelse(ID_final == "Unidentified Inorganic Material", 0.0,
ifelse(ID_final == "Unidentified Organic Matter", 0.0, TrophLevel)))))))))))))))))))))))))))))))))

propWT_GW$TL = propWT_GW$TrophLevel * propWT_GW$propWT
TL_GW = na.omit(propWT_GW)

# Keep only groups with sufficient sample sizes (n >= 20):
TL_GW = subset(TL_GW, 
                 Species == "ATF" & GW_Bin == "96-115" | 
                 Species == "ATF" & GW_Bin == "116-135" | 
                 Species == "PH" & GW_Bin == "56-75" | 
                 Species == "PH" & GW_Bin == "76-95" | 
                 Species == "PH" & GW_Bin == "96-115" | 
                 Species == "PH" & GW_Bin == "116-135")

# Predator trophic level by gape width size class:
TL_GW %>%
  group_by(Species, GW_Bin) %>%
  summarise(sum(TL)+1)  

##########
# Prey richness, diversity (H'), and evenness (J'), by predator (comparable sites only):

# Exclude non-prey items:
All.data_div = subset(All.data_sites, ID_final != "Homo sapiens" & ID_final != "Rock" & ID_final != "Rock and Shell Hash" & ID_final != "Unidentified Inorganic Material")

# ALL SIZES:
prey.data = All.data_div %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
prey.all = spread(prey.data, key = ID_final, value = sumWT)
prey.all[is.na(prey.all)] = 0

Richness = ifelse(prey.all[,2:ncol(prey.all)] > 0, 1, 0)
  sum(Richness[1,1:ncol(Richness)]) # ATF prey richness
  sum(Richness[2,1:ncol(Richness)]) # PH prey richness

# ATF #
H_ATF = subset(prey.all, Species=="ATF")
  H_ATF = H_ATF[,2:ncol(H_ATF)]
Hprime_ATF = diversity(H_ATF); Hprime_ATF
J_ATF = Hprime_ATF/log(ncol(H_ATF)); J_ATF

# PH #
H_PH = subset(prey.all, Species=="PH")
  H_PH = as.data.frame(H_PH[,2:ncol(H_PH)])
Hprime_PH = diversity(H_PH); Hprime_PH
J_PH = Hprime_PH/log(ncol(H_PH)); J_PH


# Table 2:
# FL:50 - 59 cm:
FL.50 = subset(All.data_div, FL_Bin == "50-59")
H_50 = FL.50 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.50 = spread(H_50, key = ID_final, value = sumWT)
H.50[is.na(H.50)] = 0

Richness.50 = ifelse(H.50[,2:ncol(H.50)] > 0, 1, 0)
  sum(Richness.50[1,1:ncol(Richness.50)]) # ATF prey richness only

# ATF only #
H_ATF.50 = subset(H.50, Species=="ATF")
  H_ATF.50 = H_ATF.50[,2:ncol(H.50)]
Hprime_ATF.50 = diversity(H_ATF.50); Hprime_ATF.50
J_ATF.50 = Hprime_ATF.50/log(ncol(H_ATF.50)); J_ATF.50


# FL:60 - 69 cm:
FL.60 = subset(All.data_div, FL_Bin == "60-69")
H_60 = FL.60 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.60 = spread(H_60, key = ID_final, value = sumWT)
H.60[is.na(H.60)] = 0

Richness.60 = ifelse(H.60[,2:ncol(H.60)] > 0, 1, 0)
  sum(Richness.60[1,1:ncol(Richness.60)]) # ATF prey richness
  sum(Richness.60[2,1:ncol(Richness.60)]) # PH prey richness

# ATF #
H_ATF.60 = subset(H.60, Species=="ATF")
  H_ATF.60 = H_ATF.60[,2:ncol(H.60)]
Hprime_ATF.60 = diversity(H_ATF.60); Hprime_ATF.60
J_ATF.60 = Hprime_ATF.60/log(ncol(H_ATF.60)); J_ATF.60

# PH #
H_PH.60 = subset(H.60, Species=="PH")
  H_PH.60 = H_PH.60[,2:ncol(H.60)]
Hprime_PH.60 = diversity(H_PH.60); Hprime_PH.60
J_PH.60 = Hprime_PH.60/log(ncol(H_PH.60)); J_PH.60


# FL:70 - 79 cm:
FL.70 = subset(All.data_div, FL_Bin == "70-79")
H_70 = FL.70 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.70 = spread(H_70, key = ID_final, value = sumWT)
H.70[is.na(H.70)] = 0

Richness.70 = ifelse(H.70[,2:ncol(H.70)] > 0, 1, 0)
sum(Richness.70[2,1:ncol(Richness.70)]) # PH prey richness only

# PH only #
H_PH.70 = subset(H.70, Species=="PH")
  H_PH.70 = H_PH.70[,2:ncol(H.70)]
Hprime_PH.70 = diversity(H_PH.70); Hprime_PH.70
J_PH.70 = Hprime_PH.70/log(ncol(H_PH.70)); J_PH.70


# FL:80 - 89 cm:
FL.80 = subset(All.data_div, FL_Bin == "80-89")
H_80 = FL.80 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.80 = spread(H_80, key = ID_final, value = sumWT)
H.80[is.na(H.80)] = 0

Richness.80 = ifelse(H.80[,2:ncol(H.80)] > 0, 1, 0)
sum(Richness.80[2,1:ncol(Richness.80)]) # PH prey richness pnly

# PH only #
H_PH.80 = subset(H.80, Species=="PH")
  H_PH.80 = H_PH.80[,2:ncol(H.80)]
Hprime_PH.80 = diversity(H_PH.80); Hprime_PH.80
J_PH.80 = Hprime_PH.80/log(ncol(H_PH.80)); J_PH.80


# FL:90 - 99 cm:
FL.90 = subset(All.data_div, FL_Bin == "90-99")
H_90 = FL.90 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.90 = spread(H_90, key = ID_final, value = sumWT)
H.90[is.na(H.90)] = 0

Richness.90 = ifelse(H.90[,2:ncol(H.90)] > 0, 1, 0)
sum(Richness.90[1,1:ncol(Richness.90)]) # PH prey richness only

# PH only #
H_PH.90 = subset(H.90, Species=="PH")
  H_PH.90 = H_PH.90[,2:ncol(H.90)]
Hprime_PH.90 = diversity(H_PH.90); Hprime_PH.90
J_PH.90 = Hprime_PH.90/log(ncol(H_PH.90)); J_PH.90


# FL:100 - 109 cm:
FL.100 = subset(All.data_div, FL_Bin == "100-109")
H_100 = FL.100 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.100 = spread(H_100, key = ID_final, value = sumWT)
H.100[is.na(H.100)] = 0

Richness.100 = ifelse(H.100[,2:ncol(H.100)] > 0, 1, 0)
sum(Richness.100[1,1:ncol(Richness.100)]) # PH prey richness only

# PH only #
H_PH.100 = subset(H.100, Species=="PH")
  H_PH.100 = H_PH.100[,2:ncol(H.100)]
Hprime_PH.100 = diversity(H_PH.100); Hprime_PH.100
J_PH.100 = Hprime_PH.100/log(ncol(H_PH.100)); J_PH.100


# All fork length bins, ATF (FL:50 - 69 cm):
FL.all = subset(All.data_div, FL_Bin == "50-59" | FL_Bin == "60-69")
H_all = FL.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[1,1:ncol(Richness.all)]) # ATF prey richness only

# ATF only #
H_ATF.all = subset(H.all, Species=="ATF")
  H_ATF.all = H_ATF.all[,2:ncol(H.all)]
Hprime_ATF.all = diversity(H_ATF.all); Hprime_ATF.all
J_ATF.all = Hprime_ATF.all/log(ncol(H_ATF.all)); J_ATF.all


# All fork length bins, PH (FL:60 - 109 cm):
FL.all = subset(All.data_div, 
                  FL_Bin == "60-69" | 
                  FL_Bin == "70-79" | 
                  FL_Bin == "70-79" | 
                  FL_Bin == "70-79" |
                  FL_Bin == "80-89" | 
                  FL_Bin == "90-99" | 
                  FL_Bin == "100-109")
H_all = FL.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[2,1:ncol(Richness.all)]) # PH prey richness only

# PH only #
H_PH.all = subset(H.all, Species=="PH")
  H_PH.all = H_PH.all[,2:ncol(H.all)]
Hprime_PH.all = diversity(H_PH.all); Hprime_PH.all
J_PH.all = Hprime_PH.all/log(ncol(H_PH.all)); J_PH.all

##########################################################
# Table 3:
# GH: 56 - 75 mm:
GH.56 = subset(All.data_div, GH_Bin == "56-75")
H_56 = GH.56 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.56 = spread(H_56, key = ID_final, value = sumWT)
H.56[is.na(H.56)] = 0

Richness.56 = ifelse(H.56[,2:ncol(H.56)] > 0, 1, 0)
sum(Richness.56[1,1:ncol(Richness.56)]) # PH prey richness only

# PH only #
H_PH.56 = subset(H.56, Species=="PH")
  H_PH.56 = H_PH.56[,2:ncol(H.56)]
Hprime_PH.56 = diversity(H_PH.56); Hprime_PH.56
J_PH.56 = Hprime_PH.56/log(ncol(H_PH.56)); J_PH.56


# GH: 76 - 95 mm:
GH.76 = subset(All.data_div, GH_Bin == "76-95")
H_76 = GH.76 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.76 = spread(H_76, key = ID_final, value = sumWT)
H.76[is.na(H.76)] = 0

Richness.76 = ifelse(H.76[,2:ncol(H.76)] > 0, 1, 0)
sum(Richness.76[2,1:ncol(Richness.76)]) # PH prey richness only

# PH only #
H_PH.76 = subset(H.76, Species=="PH")
  H_PH.76 = H_PH.76[,2:ncol(H.76)]
Hprime_PH.76 = diversity(H_PH.76); Hprime_PH.76
J_PH.76 = Hprime_PH.76/log(ncol(H_PH.76)); J_PH.76


# GH: 96 - 115 mm:
GH.96 = subset(All.data_div, GH_Bin == "96-115")
H_96 = GH.96 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.96 = spread(H_96, key = ID_final, value = sumWT)
H.96[is.na(H.96)] = 0

Richness.96 = ifelse(H.96[,2:ncol(H.96)] > 0, 1, 0)
sum(Richness.96[1,1:ncol(Richness.96)]) # ATF prey richness
sum(Richness.96[2,1:ncol(Richness.96)]) # PH prey richness

# ATF #
H_ATF.96 = subset(H.96, Species=="ATF")
  H_ATF.96 = H_ATF.96[,2:ncol(H.96)]
Hprime_ATF.96 = diversity(H_ATF.96); Hprime_ATF.96
J_ATF.96 = Hprime_ATF.96/log(ncol(H_ATF.96)); J_ATF.96

# PH #
H_PH.96 = subset(H.96, Species=="PH")
  H_PH.96 = H_PH.96[,2:ncol(H.96)]
Hprime_PH.96 = diversity(H_PH.96); Hprime_PH.96
J_PH.96 = Hprime_PH.96/log(ncol(H_PH.96)); J_PH.96


# GH: 116 - 135 mm:
GH.116 = subset(All.data_div, GH_Bin == "116-135")
H_116 = GH.116 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.116 = spread(H_116, key = ID_final, value = sumWT)
H.116[is.na(H.116)] = 0

Richness.116 = ifelse(H.116[,2:ncol(H.116)] > 0, 1, 0)
sum(Richness.116[1,1:ncol(Richness.116)]) # ATF prey richness
sum(Richness.116[2,1:ncol(Richness.116)]) # PH prey richness

# ATF #
H_ATF.116 = subset(H.116, Species=="ATF")
  H_ATF.116 = H_ATF.116[,2:ncol(H.116)]
Hprime_ATF.116 = diversity(H_ATF.116); Hprime_ATF.116
J_ATF.116 = Hprime_ATF.116/log(ncol(H_ATF.116)); J_ATF.116

# PH #
H_PH.116 = subset(H.116, Species=="PH")
  H_PH.116 = H_PH.116[,2:ncol(H.116)]
Hprime_PH.116 = diversity(H_PH.116); Hprime_PH.116
J_PH.116 = Hprime_PH.116/log(ncol(H_PH.116)); J_PH.116


# All gape height bins, ATF (GH: 96 - 135 cm):
GH.all = subset(All.data_div, GH_Bin == "96-115" | GH_Bin == "116-135")
H_all = GH.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[1,1:ncol(Richness.all)]) # ATF prey richness only

# ATF only #
H_ATF.all = subset(H.all, Species=="ATF")
  H_ATF.all = H_ATF.all[,2:ncol(H.all)]
Hprime_ATF.all = diversity(H_ATF.all); Hprime_ATF.all
J_ATF.all = Hprime_ATF.all/log(ncol(H_ATF.all)); J_ATF.all


# All gape height bins, PH (GH:56 - 135 cm):
GH.all = subset(All.data_div, 
                 GH_Bin == "56-75" |
                 GH_Bin == "76-95" |
                 GH_Bin == "96-115" |
                 GH_Bin == "116-125")
H_all = GH.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[2,1:ncol(Richness.all)]) # PH prey richness only

# PH only #
H_PH.all = subset(H.all, Species=="PH")
  H_PH.all = H_PH.all[,2:ncol(H.all)]
Hprime_PH.all = diversity(H_PH.all); Hprime_PH.all
J_PH.all = Hprime_PH.all/log(ncol(H_PH.all)); J_PH.all

##########################################################
# Table 3:
# GW: 56 - 75 mm:
GW.56 = subset(All.data_div, GW_Bin == "56-75")
H_56 = GW.56 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.56 = spread(H_56, key = ID_final, value = sumWT)
H.56[is.na(H.56)] = 0

Richness.56 = ifelse(H.56[,2:ncol(H.56)] > 0, 1, 0)
sum(Richness.56[1,1:ncol(Richness.56)]) # PH prey richness only

# PH only #
H_PH.56 = subset(H.56, Species=="PH")
  H_PH.56 = H_PH.56[,2:ncol(H.56)]
Hprime_PH.56 = diversity(H_PH.56); Hprime_PH.56
J_PH.56 = Hprime_PH.56/log(ncol(H_PH.56)); J_PH.56


# GW: 76 - 95 mm:
GW.76 = subset(All.data_div, GW_Bin == "76-95")
H_76 = GW.76 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.76 = spread(H_76, key = ID_final, value = sumWT)
H.76[is.na(H.76)] = 0

Richness.76 = ifelse(H.76[,2:ncol(H.76)] > 0, 1, 0)
sum(Richness.76[2,1:ncol(Richness.76)]) # PH prey richness only

# PH only #
H_PH.76 = subset(H.76, Species=="PH")
  H_PH.76 = H_PH.76[,2:ncol(H.76)]
Hprime_PH.76 = diversity(H_PH.76); Hprime_PH.76
J_PH.76 = Hprime_PH.76/log(ncol(H_PH.76)); J_PH.76


# GW: 96 - 115 mm:
GW.96 = subset(All.data_div, GW_Bin == "96-115")
H_96 = GW.96 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.96 = spread(H_96, key = ID_final, value = sumWT)
H.96[is.na(H.96)] = 0

Richness.96 = ifelse(H.96[,2:ncol(H.96)] > 0, 1, 0)
sum(Richness.96[1,1:ncol(Richness.96)]) # ATF prey richness
sum(Richness.96[2,1:ncol(Richness.96)]) # PH prey richness

# ATF #
H_ATF.96 = subset(H.96, Species=="ATF")
  H_ATF.96 = H_ATF.96[,2:ncol(H.96)]
Hprime_ATF.96 = diversity(H_ATF.96); Hprime_ATF.96
J_ATF.96 = Hprime_ATF.96/log(ncol(H_ATF.96)); J_ATF.96

# PH #
H_PH.96 = subset(H.96, Species=="PH")
  H_PH.96 = H_PH.96[,2:ncol(H.96)]
Hprime_PH.96 = diversity(H_PH.96); Hprime_PH.96
J_PH.96 = Hprime_PH.96/log(ncol(H_PH.96)); J_PH.96


# GW: 116 - 135 mm:
GW.116 = subset(All.data_div, GW_Bin == "116-135")
H_116 = GW.116 %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.116 = spread(H_116, key = ID_final, value = sumWT)
H.116[is.na(H.116)] = 0

Richness.116 = ifelse(H.116[,2:ncol(H.116)] > 0, 1, 0)
sum(Richness.116[1,1:ncol(Richness.116)]) # ATF prey richness
sum(Richness.116[2,1:ncol(Richness.116)]) # PH prey richness

# ATF #
H_ATF.116 = subset(H.116, Species=="ATF")
  H_ATF.116 = H_ATF.116[,2:ncol(H.116)]
Hprime_ATF.116 = diversity(H_ATF.116); Hprime_ATF.116
J_ATF.116 = Hprime_ATF.116/log(ncol(H_ATF.116)); J_ATF.116

# PH #
H_PH.116 = subset(H.116, Species=="PH")
  H_PH.116 = H_PH.116[,2:ncol(H.116)]
Hprime_PH.116 = diversity(H_PH.116); Hprime_PH.116
J_PH.116 = Hprime_PH.116/log(ncol(H_PH.116)); J_PH.116


# All gape height bins, ATF (GW:96 - 135 cm):
GW.all = subset(All.data_div, GW_Bin == "96-115" | GW_Bin == "116-135")
H_all = GW.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[1,1:ncol(Richness.all)]) # ATF prey richness only

# ATF only #
H_ATF.all = subset(H.all, Species=="ATF")
  H_ATF.all = H_ATF.all[,2:ncol(H_ATF.all)]
Hprime_ATF.all = diversity(H_ATF.all); Hprime_ATF.all
J_ATF.all = Hprime_ATF.all/log(ncol(H_ATF.all)); J_ATF.all


# All gape height bins, PH (GW:56 - 135 cm):
GW.all = subset(All.data_div, 
                GW_Bin == "56-75" |
                  GW_Bin == "76-95" |
                  GW_Bin == "96-115" |
                  GW_Bin == "116-125")
H_all = GW.all %>%
  group_by(Species, ID_final) %>%
  summarise(sumWT = sum(Prey.Mass..g.))
H.all = spread(H_all, key = ID_final, value = sumWT)
H.all[is.na(H.all)] = 0

Richness.all = ifelse(H.all[,2:ncol(H.all)] > 0, 1, 0)
sum(Richness.all[2,1:ncol(Richness.all)]) # PH prey richness only

# PH only #
H_PH.all = subset(H.all, Species=="PH")
  H_PH.all = H_PH.all[,2:ncol(H_PH.all)]
Hprime_PH.all = diversity(H_PH.all); Hprime_PH.all
J_PH.all = Hprime_PH.all/log(ncol(H_PH.all)); J_PH.all

################################################################
# Fig. 4: Diet compositions by size metric and class - (all sites, years, and months combined):
# top panel - fork length:
All.data$Species = as.factor(All.data$Species)
propWT.fl = All.data %>%
  group_by(Species, FL_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.fl = as.data.frame(propWT.fl)

Predator = unique(propWT.fl$Species)
FL = unique(propWT.fl$FL_Bin)
PreyTaxa = unique(propWT.fl$ID_final_Mod)
combinations = expand.grid(Species = Predator, FL_Bin = FL, ID_final_Mod = PreyTaxa)

data.fl = full_join(propWT.fl, combinations, by = c("Species" = "Species", "FL_Bin" = "FL_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, FL_Bin, ID_final_Mod)

# Remove species-size class combinations with insufficient diet samples (n >= 20):
data.fl = subset(data.fl, Species != "ATF" | 
                                  FL_Bin != "30-39" & 
                                  FL_Bin != "40-49" & 
                
                                  FL_Bin != "70-79" & 
                                  FL_Bin != "80-89" & 
                                  FL_Bin != "90-99" & 
                                  FL_Bin != "100-109" & 
                                  FL_Bin != "110-119" & 
                                  FL_Bin != "120-129" & 
                                  FL_Bin != "130-139" & 
                                  FL_Bin != "140-149" & 
                                  FL_Bin != "150-159")

data.fl = subset(data.fl, Species != "PH" | 
                                  FL_Bin != "30-39" & 
                                  FL_Bin != "40-49" & 
                                  FL_Bin != "50-59" & 
                
                                  FL_Bin != "110-119" & 
                                  FL_Bin != "120-129" & 
                                  FL_Bin != "130-139" & 
                                  FL_Bin != "140-149" & 
                                  FL_Bin != "150-159")

levels(data.fl$Species) = c("Arrowtooth Flounder", "Pacific Halibut")
PropPrey.Spp.FL = ggplot(data.fl, aes(x = FL_Bin, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Species, nrow=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill="transparent"),
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, colour="white", size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        strip.text = element_text(vjust=1.0, family="Arial", size=12),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(5,18,5,5)) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=FALSE) +
  labs(x="Fork Length (cm)", y="Proportions of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
  scale_x_discrete(expand = c(0.004, 0.004), labels = c("50", "60", "70", "80", "90", "100")) 
ggsave(PropPrey.Spp.FL, filename="Plots/PropPrey_Spp_FL.png", dpi=500, width=6, height=3.75, units="in")


# middle panel - gape height:
propWT.gh = All.data %>%
  group_by(Species, GH_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.gh = as.data.frame(propWT.gh)

Predator = unique(propWT.gh$Species)
GH = unique(propWT.gh$GH_Bin)
PreyTaxa = unique(propWT.gh$ID_final_Mod)
combinations = expand.grid(Species = Predator, GH_Bin = GH, ID_final_Mod = PreyTaxa)

data.gh = full_join(propWT.gh, combinations, by = c("Species" = "Species", "GH_Bin" = "GH_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, GH_Bin, ID_final_Mod)

# Remove species-size class combinations with insufficient diet samples (n >= 20):
data.gh = subset(data.gh, Species != "ATF" | 
                GH_Bin != "<36" &
                GH_Bin != "36-55" &
                GH_Bin != "56-75" & 
                GH_Bin != "76-95" & 
                
                GH_Bin != "136-155" & 
                GH_Bin != "156-175" & 
                GH_Bin != "176-195" &
                GH_Bin != ">195")

data.gh = subset(data.gh, Species != "PH" | 
                GH_Bin != "<36" &
                GH_Bin != "36-55" &
                
                GH_Bin != "136-155" & 
                GH_Bin != "156-175" &              
                GH_Bin != "176-195" &
                GH_Bin != ">195")

levels(data.gh$Species) = c("Arrowtooth Flounder", "Pacific Halibut")
PropPrey.Spp.GH = ggplot(data.gh, aes(x = GH_Bin, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Species, nrow=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill="transparent"),
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, colour="white", size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        strip.text = element_text(vjust=1.0, family="Arial", size=12),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(5,18,5,5)) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=FALSE) +
  labs(x="Gape Height (mm)", y="Proportions of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
  scale_x_discrete(expand = c(0.004, 0.004), labels = c("56", "76", "96", "116")) 
ggsave(PropPrey.Spp.GH, filename="Plots/PropPrey_Spp_GH.png", dpi=500, width=6, height=3.75, units="in")


# bottom panel - gape width:
propWT.gw = All.data %>%
  group_by(Species, GW_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.gw = as.data.frame(propWT.gw)

Predator = unique(propWT.gw$Species)
GW = unique(propWT.gw$GW_Bin)
PreyTaxa = unique(propWT.gw$ID_final_Mod)
combinations = expand.grid(Species = Predator, GW_Bin = GW, ID_final_Mod = PreyTaxa)

data.gw = full_join(propWT.gw, combinations, by = c("Species" = "Species", "GW_Bin" = "GW_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, GW_Bin, ID_final_Mod)

# Remove species-size class combinations with insufficient diet samples (n >= 20):
data.gw = subset(data.gw, Species != "ATF" | 
                GW_Bin != "<36" &
                GW_Bin != "36-55" &
                GW_Bin != "56-75" & 
                GW_Bin != "76-95" & 
                
                GW_Bin != "136-155" & 
                GW_Bin != "156-175" & 
                GW_Bin != "176-195" &
                GW_Bin != ">195")

data.gw = subset(data.gw, Species != "PH" | 
                GW_Bin != "<36" &
                GW_Bin != "36-55" &
                
                GW_Bin != "136-155" & 
                GW_Bin != "156-175" &              
                GW_Bin != "176-195" &
                GW_Bin != ">195")

levels(data.gw$Species) = c("Arrowtooth Flounder", "Pacific Halibut")
PropPrey.Spp.GW = ggplot(data.gw, aes(x = GW_Bin, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Species, nrow=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.background = element_rect(fill="transparent"),
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(family="Arial", vjust=-0.13, size=12),
        axis.title.y = element_text(family="Arial", vjust = 2.5, colour="white", size=12.5),
        strip.background = element_rect(colour="white",fill="white"),
        strip.text = element_text(vjust=1.0, family="Arial", size=12),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(5,18,5,5)) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=FALSE) +
  labs(x="Gape Width (mm)", y="Proportions of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
  scale_x_discrete(expand = c(0.004, 0.004), labels = c("56", "76", "96", "116")) 
ggsave(PropPrey.Spp.GW, filename="Plots/PropPrey_Spp_GW.png", dpi=500, width=6, height=3.75, units="in")

################################################################
# Calculate Schoener's index of dietary overlap (Schoener 1968):
# Fig. 6: Diet compositions by sampling site - (all sizes, years, and months combined):
propLoc = All.data_sites %>%
  group_by(Species, Location, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 

propLoc_wide = dcast(propLoc, Location + ID_final_Mod ~ Species, value.var = "propWT", fill = 0)
  propLoc_wide$overlap = with(propLoc_wide, abs(PH - ATF))
D_Loc = propLoc_wide %>%
  group_by(Location) %>%
  summarize(n = sum(overlap)) 
D_Loc$D = round((1 - 0.5 * D_Loc$n), digits=3); D_Loc$D 


levels(All.data$Location) = c("Lynn Canal", "Fav.-Sag. Channels", "Point Howard", "Funter Bay", "Point Couverden", "Icy Strait")

propWT = All.data %>%
  group_by(Species, Location, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT = as.data.frame(propWT)

Predator = unique(propWT$Species)
Site = unique(propWT$Location)
PreyTaxa = unique(propWT$ID_final_Mod)
combinations = expand.grid(Species = Predator, Location = Site, ID_final_Mod = PreyTaxa)

data = full_join(propWT, combinations, by = c("Species" = "Species", "Location" = "Location", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, Location, ID_final_Mod)

data = subset(data, Species != "ATF" | Location != "Funter Bay" & Location != "Point Couverden" & Location != "Icy Strait")

PropPrey.Spp.Loc = ggplot(data, aes(x = Species, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_bar(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Location, ncol=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.title.y = element_text(vjust = 2.5, size = 12.5, family="Arial"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white",fill="white")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=F) +
  labs(x="", y="Proportion of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks=c(0.0,0.25,0.50,0.75,1.0), labels=c("0.0","","0.5","","1.0")) +
  scale_x_discrete(expand = c(0.001, 0.001)) +
  theme(legend.background = element_rect(fill="transparent")) 
ggsave(PropPrey.Spp.Loc, filename="Plots/PropPrey_Spp_Loc.png", dpi=500, width=3.75, height=7.5, units="in", bg = "transparent")

##########################################################
# Fig. S4: by size class and metric - (all comparable sites, years, and months combined):
# left panel - 60-69 cm fork length:
All.data_sites$Species = as.factor(All.data_sites$Species)
propWT.fl = All.data_sites %>%
  group_by(Species, FL_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.fl = as.data.frame(propWT.fl)

Predator = unique(propWT.fl$Species)
FL = unique(propWT.fl$FL_Bin)
PreyTaxa = unique(propWT.fl$ID_final_Mod)
combinations = expand.grid(Species = Predator, FL_Bin = FL, ID_final_Mod = PreyTaxa)

data.fl = full_join(propWT.fl, combinations, by = c("Species" = "Species", "FL_Bin" = "FL_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, FL_Bin, ID_final_Mod)

FL_60 = subset(data.fl, FL_Bin == "60-69")
  levels(FL_60$Species) = c("ATF", "PH")

PropPrey.Spp.FL60 = ggplot(FL_60, aes(x = Species, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.title.y = element_text(vjust = 2.5, size = 12.5, family="Arial"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white",fill="white")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=F) +
  labs(x="", y="Proportion of Prey by Weight") +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_x_discrete(expand = c(0.001, 0.001)) +
  theme(legend.background = element_rect(fill="transparent")) 
ggsave(PropPrey.Spp.FL60, filename="Plots/PropPrey_Spp_FL60.png", dpi=500, width=3.75, height=4.5, units="in")

FL_60_wide = dcast(FL_60, ID_final_Mod ~ Species, value.var = "propWT", fill = 0)
FL_60_wide$overlap = with(FL_60_wide, abs(PH - ATF))
1 - 0.5 * sum(FL_60_wide$overlap)


# upper-right panel - 96-135 mm gape height:
propWT.gh = All.data_sites %>%
  group_by(Species, GH_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.gh = as.data.frame(propWT.gh)

Predator = unique(propWT.gh$Species)
GH = unique(propWT.gh$GH_Bin)
PreyTaxa = unique(propWT.gh$ID_final_Mod)
combinations = expand.grid(Species = Predator, GH_Bin = GH, ID_final_Mod = PreyTaxa)

data.gh = full_join(propWT.gh, combinations, by = c("Species" = "Species", "GH_Bin" = "GH_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, GH_Bin, ID_final_Mod)

GH.96_135 = subset(data.gh, 
                    GH_Bin == "96-115" | 
                    GH_Bin == "116-135")
  levels(GH.96_135$Species) = c("ATF", "PH")
  
PropPrey.Spp.GH.96_135 = ggplot(GH.96_135, aes(x = GH_Bin, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Species, nrow=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.title.y = element_text(vjust = 2.5, size = 12.5, family="Arial"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white",fill="white")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=F) +
  labs(x="Gape Height (mm)", y="Proportions of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
  scale_x_discrete(expand = c(0.004, 0.004), labels = c("56", "76", "96", "116")) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(plot.margin = margin(5,18,5,5))
ggsave(PropPrey.Spp.GH.96_135, filename="Plots/PropPrey_Spp_GH96_135.png", dpi=500, width=6, height=3.75, units="in")

GH.96_135_wide = dcast(GH.96_135, GH_Bin + ID_final_Mod ~ Species, value.var = "propWT", fill = 0)
GH.96_135_wide$overlap = with(GH.96_135_wide, abs(PH - ATF))
D_GH = GH.96_135_wide %>%
  group_by(GH_Bin) %>%
  summarize(n = sum(overlap)) 
D_GH$D = round((1 - 0.5 * D_GH$n), digits=3); D_GH$D


# lower-right panel - 96-135 mm gape width:
propWT.gw = All.data_sites %>%
  group_by(Species, GW_Bin, ID_final_Mod) %>% 
  summarise(n = sum(Prey.Mass..g.)) %>%
  mutate(propWT = n / sum(n)) 
propWT.gw = as.data.frame(propWT.gw)

Predator = unique(propWT.gw$Species)
GW = unique(propWT.gw$GW_Bin)
PreyTaxa = unique(propWT.gw$ID_final_Mod)
combinations = expand.grid(Species = Predator, GW_Bin = GW, ID_final_Mod = PreyTaxa)

data.gw = full_join(propWT.gw, combinations, by = c("Species" = "Species", "GW_Bin" = "GW_Bin", "ID_final_Mod" = "ID_final_Mod")) %>%
  mutate(propWT = ifelse(is.na(propWT), 0, propWT)) %>%
  arrange(Species, GW_Bin, ID_final_Mod)

GW.96_135 = subset(data.gw, 
                     GW_Bin == "96-115" | 
                     GW_Bin == "116-135")
levels(GW.96_135$Species) = c("ATF", "PH")

PropPrey.Spp.GW.96_135 = ggplot(GW.96_135, aes(x = GW_Bin, y = propWT, group = ID_final_Mod, fill = ID_final_Mod, order = -as.numeric(ID_final_Mod))) +
  geom_area(position = "stack", stat = "identity", na.rm=T) +
  theme_bw() +
  facet_wrap(~ Species, nrow=1) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(colour="black", size=1, line="solid"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.text = element_text(family="Arial", size=8),
        axis.title.y = element_text(vjust = 2.5, size = 12.5, family="Arial"),
        axis.text = element_text(family="Arial", size=11.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white",fill="white")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_fill_manual(values = c("Mollusca"="lightpink3", "Teuthida"="brown4", "Octopodidae"="firebrick3", "Crustacea"="red", "Pandalidae"="chocolate", "Lithodidae"="chocolate1", "Majoidea"="orange", "Chionoecetes bairdi"="gold", "Metacarcinus magister"="khaki1", "Paguroidea"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupea pallasii"="cadetblue1", "Salmoniformes"="turquoise3", "Gadus"="blue", "Gadus chalcogrammus"="mediumblue", "Gadus macrocephalus"="blue4", "Scorpaeniformes"="purple2", "Pleuronectiformes"="mediumvioletred", "Benthic Material"="gray91", "Other"="azure4"), name="Prey Taxa", labels = c(" Mollusca", " Teuthida", " Octopodidae", " Crustacea", " Pandalidae", " Lithodidae", " Majoidea", expression(paste(italic(" Chionoecetes spp."))), expression(paste(italic(" Metacarcinus magister"))), " Paguroidea", " Chondrichthyes", " Teleostei", expression(paste(italic(" Clupea pallasii"))), " Salmoniformes", expression(paste(italic(" Gadus spp."))), expression(paste(italic("   G. chalcogrammus"))), expression(paste(italic("   G. macrocephalus"))), " Scorpaeniformes", " Pleuronectiformes", " Benthic Material", " Other"), drop=F) +
  labs(x="Gape Width (mm)", y="Proportions of Prey by Weight") +
  scale_y_continuous(expand = c(0.003, 0.003), breaks = c(0.00,0.25,0.50,0.75,1.00)) +
  scale_x_discrete(expand = c(0.004, 0.004), labels = c("56", "76", "96", "116")) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(plot.margin = margin(5,18,5,5))
ggsave(PropPrey.Spp.GW.96_135, filename="Plots/PropPrey_Spp_GW96_135.png", dpi=500, width=6, height=3.75, units="in")

GW.96_135_wide = dcast(GW.96_135, GW_Bin + ID_final_Mod ~ Species, value.var = "propWT", fill = 0)
GW.96_135_wide$overlap = with(GW.96_135_wide, abs(PH - ATF))
D_GW = GW.96_135_wide %>%
  group_by(GW_Bin) %>%
  summarize(n = sum(overlap)) 
D_GW$D = round((1 - 0.5 * D_GW$n), digits=3); D_GW$D