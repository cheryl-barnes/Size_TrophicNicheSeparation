# The code below cumulative prey curves from diets of Pacific Halibut (PH) and Arrowtooth Flounder (ATF). Predators were collected using hook-and-line near Juneau, AK in 2015 and 2016. See Barnes et al. (2021) for details.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# Reference:
# Barnes, C.L., A.H. Beaudreau, and R.N. Yamada (2021). The Role of size in trophic niche separation between two groundfish predators in Alaskan waters. Marine and Coastal Fisheries. doi:10.1002/mcf2.10141

setwd("~/Desktop/SEAK_FieldStudy/")
require(dplyr)
require(reshape2)
require(vegan)
###################################################################### 
# Import and prepare diet data:
PreyData = read.csv("Data/PH_ATF_DietData_SEAK.csv")

# Remove spring and fall months with very few samples (focus = summer):
PreyData$Year = as.factor(PreyData$Year)
PreyData = subset(PreyData, Month != 5 & Month != 9)

# Consolidate individual capture locations into sampling sites:
PreyData$Location = with(PreyData, 
ifelse(Grid.. < 1, "Chatham.Strait",
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

# Format long to wide:
PreyData_wide = dcast(All.data, Month + Day + Year + Location + PredatorData.PredatorID + Species + FL_Bin + GH_Bin + GW_Bin ~ ID_final_Mod, value.var = "Prey.Mass..g.", fun=sum)

###################################################################
# Construct cumulative prey curves to determine whether or not sufficient samples were collected for any one sampling group:
PH_data = subset(PreyData_wide, Species == "PH")
ATF_data = subset(PreyData_wide, Species == "ATF")
###################################################################
# Fig. S1 a) by year:

### 2015:
PH_2015 = subset(PH_data, Year == 2015)
ATF_2015 = subset(ATF_data, Year == 2015)

# Remove 'Species' column from data frames:
PH_2015sub = ifelse(PH_2015[10:ncol(PH_2015)] > 1, 1, 0)
ATF_2015sub = ifelse(ATF_2015[10:ncol(ATF_2015)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_2015 = specaccum(PH_2015sub, method="random", permutations=100)
curve.ATF_2015 = specaccum(ATF_2015sub, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-2015.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_2015, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_2015, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("2015: All Sites, Sizes, and Months", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 2016:
PH_2016 = subset(PH_data, Year == 2016)
ATF_2016 = subset(ATF_data, Year == 2016)

# Remove 'Species' column from data frames:
PH_2016sub = ifelse(PH_2016[10:ncol(PH_2016)] > 1, 1, 0)
ATF_2016sub = ifelse(ATF_2016[10:ncol(ATF_2016)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_2016 = specaccum(PH_2016sub, method="random", permutations=100)
curve.ATF_2016 = specaccum(ATF_2016sub, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-2016.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_2016, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_2016, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("2016: All Sites, Sizes, and Months", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()
###################################################################
# Fig. S1 b) by site:
PH_LC = subset(PH_data, Location == "Lynn.Canal")
PH_FS = subset(PH_data, Location == "Favorite-Saginaw.Channel")
PH_PH = subset(PH_data, Location == "Pt.Howard")
PH_FB = subset(PH_data, Location == "Funter.Bay")
PH_PC = subset(PH_data, Location == "Pt.Couverden")
PH_IS = subset(PH_data, Location == "Icy.Strait")

ATF_LC = subset(ATF_data, Location == "Lynn.Canal")
ATF_FS = subset(ATF_data, Location == "Favorite-Saginaw.Channel")
ATF_PH = subset(ATF_data, Location == "Pt.Howard")
ATF_FB = subset(ATF_data, Location == "Funter.Bay")
ATF_PC = subset(ATF_data, Location == "Pt.Couverden")
ATF_IS = subset(ATF_data, Location == "Icy.Strait")


### Lynn Canal:
# Remove 'Species' column from data frames:
PH.sub_LC = ifelse(PH_LC[10:ncol(PH_LC)] > 1, 1, 0)
ATF.sub_LC = ifelse(ATF_LC[10:ncol(ATF_LC)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_LC = specaccum(PH.sub_LC, method="random", permutations=100)
curve.ATF_LC = specaccum(ATF.sub_LC, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_LC.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_LC, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_LC, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Lynn Canal; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### Favorite-Saginaw Channels:
# Remove 'Species' column from data frames:
PH.sub_FS = ifelse(PH_FS[10:ncol(PH_FS)] > 1, 1, 0)
ATF.sub_FS = ifelse(ATF_FS[10:ncol(ATF_FS)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_FS = specaccum(PH.sub_FS, method="random", permutations=100)
curve.ATF_FS = specaccum(ATF.sub_FS, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_FS.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_FS, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_FS, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Favorite-Saginaw Channels; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### Point Howard:
# Remove 'Species' column from data frames:
PH.sub_PH = ifelse(PH_PH[10:ncol(PH_PH)] > 1, 1, 0)
ATF.sub_PH = ifelse(ATF_PH[10:ncol(ATF_PH)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_PH = specaccum(PH.sub_PH, method="random", permutations=100)
curve.ATF_PH = specaccum(ATF.sub_PH, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_PH.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_PH, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_PH, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Point Howard; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### Funter Bay:
# Remove 'Species' column from data frames:
PH.sub_FB = ifelse(PH_FB[10:ncol(PH_FB)] > 1, 1, 0)
ATF.sub_FB = ifelse(ATF_FB[10:ncol(ATF_FB)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_FB = specaccum(PH.sub_FB, method="random", permutations=100)
curve.ATF_FB = specaccum(ATF.sub_FB, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_FB.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_FB, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_FB, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Funter Bay; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### Point Couverden:
# Remove 'Species' column from data frames:
  PH.sub_PC = ifelse(PH_PC[10:ncol(PH_PC)] > 1, 1, 0)
# ATF.sub_PC = ifelse(ATF_PC[10:ncol(ATF_PC)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_PC = specaccum(PH.sub_PC, method="random", permutations=100)
# curve.ATF_PC = specaccum(ATF.sub_PC, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_PC.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_PC, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_PC, 
    # ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
    # col = "blue",
    # xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Point Couverden; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### Icy Strait:
# Remove 'Species' column from data frames:
PH.sub_IS = ifelse(PH_IS[10:ncol(PH_IS)] > 1, 1, 0)
ATF.sub_IS = ifelse(ATF_IS[10:ncol(ATF_IS)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_IS = specaccum(PH.sub_IS, method="random", permutations=100)
curve.ATF_IS = specaccum(ATF.sub_IS, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-Site_IS.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_IS, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_IS, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("Icy Strait; All Sizes and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()

###################################################################
# Fig. S1 c) by fork length bin:
PH_60 = subset(PH_data, FL_Bin == "60-69")
PH_70 = subset(PH_data, FL_Bin == "70-79")
PH_80 = subset(PH_data, FL_Bin == "80-89")
PH_90 = subset(PH_data, FL_Bin == "90-99")
PH_100 = subset(PH_data, FL_Bin == "100-109")

ATF_50 = subset(ATF_data, FL_Bin == "50-59")
ATF_60 = subset(ATF_data, FL_Bin == "60-69")


### 50-59 cm:
# Remove 'Species' column from data frames:
# PH.sub_50 = ifelse(PH_50[10:ncol(PH_50)] > 1, 1, 0)
  ATF.sub_50 = ifelse(ATF_50[10:ncol(ATF_50)] > 1, 1, 0)

# Species accumulation curves:
# curve.PH_50 = specaccum(PH.sub_50, method="random", permutations=100)
  curve.ATF_50 = specaccum(ATF.sub_50, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL50.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
# plot(curve.PH_50, 
     # main="",
     # ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), # ci.lty = 0, 
     # col = "red",
     # xlab = "", ylab="", 
     # xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_50, 
     main="",
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlab = "", ylab="",
     xlim = c(0,100), ylim = c(0,25))

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 50-59 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 60-69 cm:
# Remove 'Species' column from data frames:
PH.sub_60 = ifelse(PH_60[10:ncol(PH_60)] > 1, 1, 0)
ATF.sub_60 = ifelse(ATF_60[10:ncol(ATF_60)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_60 = specaccum(PH.sub_60, method="random", permutations=100)
curve.ATF_60 = specaccum(ATF.sub_60, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL60.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_60, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_60, 
ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
col = "blue",
xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 60-69 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 70-79 cm:
# Remove 'Species' column from data frames:
  PH.sub_70 = ifelse(PH_70[10:ncol(PH_70)] > 1, 1, 0)
# ATF.sub_70 = ifelse(ATF_70[10:ncol(ATF_70)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_70 = specaccum(PH.sub_70, method="random", permutations=100)
# curve.ATF_70 = specaccum(ATF.sub_70, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL70.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_70, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_70, 
     # ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     # col = "blue",
     # xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 70-79 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 80-89 cm:
# Remove 'Species' column from data frames:
PH.sub_80 = ifelse(PH_80[10:ncol(PH_80)] > 1, 1, 0)
# ATF.sub_80 = ifelse(ATF_80[10:ncol(ATF_80)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_80 = specaccum(PH.sub_80, method="random", permutations=100)
# curve.ATF_80 = specaccum(ATF.sub_80, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL80.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_80, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_80, 
# ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
# col = "blue",
# xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 80-89 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 90-99 cm:
# Remove 'Species' column from data frames:
  PH.sub_90 = ifelse(PH_90[10:ncol(PH_90)] > 1, 1, 0)
# ATF.sub_90 = ifelse(ATF_90[10:ncol(ATF_90)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_90 = specaccum(PH.sub_90, method="random", permutations=100)
# curve.ATF_90 = specaccum(ATF.sub_90, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL90.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_90, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_90, 
# ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
# col = "blue",
# xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 90-99 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 100-109 cm:
# Remove 'Species' column from data frames:
  PH.sub_100 = ifelse(PH_100[10:ncol(PH_100)] > 1, 1, 0)
# ATF.sub_100 = ifelse(ATF_100[10:ncol(ATF_100)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_100 = specaccum(PH.sub_100, method="random", permutations=100)
# curve.ATF_100 = specaccum(ATF.sub_100, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-FL100.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_100, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_100, 
# ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
# col = "blue",
# xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("FL: 100-109 cm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()

###################################################################
# Fig. S1 d) by gape height bin:
PH_56 = subset(PH_data, GH_Bin == "56-75")
PH_76 = subset(PH_data, GH_Bin == "76-95")
PH_96 = subset(PH_data, GH_Bin == "96-115")
PH_116 = subset(PH_data, GH_Bin == "116-135")
PH_136 = subset(PH_data, GH_Bin == "136-155")

ATF_56 = subset(ATF_data, GH_Bin == "56-75")
ATF_76 = subset(ATF_data, GH_Bin == "76-95")
ATF_96 = subset(ATF_data, GH_Bin == "96-115")
ATF_116 = subset(ATF_data, GH_Bin == "116-135")
ATF_136 = subset(ATF_data, GH_Bin == "136-155")


### 56-75 mm:
# Remove 'Species' column from data frames:
  PH.sub_56 = ifelse(PH_56[10:ncol(PH_56)] > 1, 1, 0)
# ATF.sub_56 = ifelse(ATF_56[10:ncol(ATF_56)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_56 = specaccum(PH.sub_56, method="random", permutations=100)
# curve.ATF_56 = specaccum(ATF.sub_56, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GH56.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_56, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_56, 
     # ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     # col = "blue",
     # xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GH: 56-75 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 76-95 mm:
# Remove 'Species' column from data frames:
  PH.sub_76 = ifelse(PH_76[10:ncol(PH_76)] > 1, 1, 0)
# ATF.sub_76 = ifelse(ATF_76[10:ncol(ATF_76)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_76 = specaccum(PH.sub_76, method="random", permutations=100)
# curve.ATF_76 = specaccum(ATF.sub_76, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GH76.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_76, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_76, 
     # ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     # col = "blue",
     # xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GH: 76-95 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 96-115 mm:
# Remove 'Species' column from data frames:
PH.sub_96 = ifelse(PH_96[10:ncol(PH_96)] > 1, 1, 0)
ATF.sub_96 = ifelse(ATF_96[10:ncol(ATF_96)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_96 = specaccum(PH.sub_96, method="random", permutations=100)
curve.ATF_96 = specaccum(ATF.sub_96, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GH96.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_96, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_96, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GH: 96-115 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 116-135 mm:
# Remove 'Species' column from data frames:
PH.sub_116 = ifelse(PH_116[10:ncol(PH_116)] > 1, 1, 0)
ATF.sub_116 = ifelse(ATF_116[10:ncol(ATF_116)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_116 = specaccum(PH.sub_116, method="random", permutations=100)
curve.ATF_116 = specaccum(ATF.sub_116, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GH116.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_116, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_116, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GH: 116-135 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 136-155 mm:
# Remove 'Species' column from data frames:
PH.sub_136 = ifelse(PH_136[10:ncol(PH_136)] > 1, 1, 0)
ATF.sub_136 = ifelse(ATF_136[10:ncol(ATF_136)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_136 = specaccum(PH.sub_136, method="random", permutations=100)
curve.ATF_136 = specaccum(ATF.sub_136, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GH136.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_136, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_136, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GH: 136-155 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()

###################################################################
# Fig. S1 e) by gape width bin:
PH_56 = subset(PH_data, GW_Bin == "56-75")
PH_76 = subset(PH_data, GW_Bin == "76-95")
PH_96 = subset(PH_data, GW_Bin == "96-115")
PH_116 = subset(PH_data, GW_Bin == "116-135")
PH_136 = subset(PH_data, GW_Bin == "136-155")

ATF_56 = subset(ATF_data, GW_Bin == "56-75")
ATF_76 = subset(ATF_data, GW_Bin == "76-95")
ATF_96 = subset(ATF_data, GW_Bin == "96-115")
ATF_116 = subset(ATF_data, GW_Bin == "116-135")
ATF_136 = subset(ATF_data, GW_Bin == "136-155")


### 56-75 mm:
# Remove 'Species' column from data frames:
  PH.sub_56 = ifelse(PH_56[10:ncol(PH_56)] > 1, 1, 0)
# ATF.sub_56 = ifelse(ATF_56[10:ncol(ATF_56)] > 1, 1, 0)

# Species accumulation curves:
  curve.PH_56 = specaccum(PH.sub_56, method="random", permutations=100)
# curve.ATF_56 = specaccum(ATF.sub_56, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GW56.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_56, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
# plot(curve.ATF_56, 
# ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
# col = "blue",
# xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GW: 56-75 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 76-95 mm:
# Remove 'Species' column from data frames:
PH.sub_76 = ifelse(PH_76[10:ncol(PH_76)] > 1, 1, 0)
ATF.sub_76 = ifelse(ATF_76[10:ncol(ATF_76)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_76 = specaccum(PH.sub_76, method="random", permutations=100)
curve.ATF_76 = specaccum(ATF.sub_76, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GW76.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_76, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_76, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GW: 76-95 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 96-115 mm:
# Remove 'Species' column from data frames:
PH.sub_96 = ifelse(PH_96[10:ncol(PH_96)] > 1, 1, 0)
ATF.sub_96 = ifelse(ATF_96[10:ncol(ATF_96)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_96 = specaccum(PH.sub_96, method="random", permutations=100)
curve.ATF_96 = specaccum(ATF.sub_96, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GW96.pdf",width=6.5,height=4.25)
par(mfrow=c(1, 1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_96, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_96, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GW: 96-115 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 116-135 mm:
# Remove 'Species' column from data frames:
PH.sub_116 = ifelse(PH_116[10:ncol(PH_116)] > 1, 1, 0)
ATF.sub_116 = ifelse(ATF_116[10:ncol(ATF_116)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_116 = specaccum(PH.sub_116, method="random", permutations=100)
curve.ATF_116 = specaccum(ATF.sub_116, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GW116.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_116, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_116, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GW: 116-135 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()


### 136-155 mm:
# Remove 'Species' column from data frames:
PH.sub_136 = ifelse(PH_136[10:ncol(PH_136)] > 1, 1, 0)
ATF.sub_136 = ifelse(ATF_136[10:ncol(ATF_136)] > 1, 1, 0)

# Species accumulation curves:
curve.PH_136 = specaccum(PH.sub_136, method="random", permutations=100)
curve.ATF_136 = specaccum(ATF.sub_136, method="random", permutations=100)

pdf("Plots/CumulativePreyCurve-GW136.pdf",width=6.5,height=4.25)
par(mfrow=c(1,1), mar=c(2,2,3,1), omi=c(0.4,0.4,0,0))

# Plot PH curve with CI:
plot(curve.PH_136, 
     main="",
     ci.type = "polygon", ci.col = rgb(1, 158/255, 115/255, 0.5), ci.lty = 0, 
     col = "red",
     xlab = "", ylab="", 
     xlim = c(0,100), ylim = c(0,25))

# Add ATF curve and CIT to the same plot:
plot(curve.ATF_136, 
     ci.type = "polygon", ci.col = rgb(86/255, 180/255, 233/255, 0.5), ci.lty = 0, 
     col = "blue",
     xlim = c(0,100), ylim = c(0,25), add = TRUE)

mtext("Cumulative Prey Taxa", side=2, line=0.5, outer=T)
mtext("No. Stomachs Sampled", side=1, line=-0.0, padj=0.5, outer=T)
mtext("GW: 136-155 mm; All Sites and Time Periods", side=3, padj=-0.5, outer=F)
legend("bottomright", inset=0.005, legend=c("PH", "ATF"), bty="n", lwd=1, col=c("red", "blue"), lty=1, x.intersp=0.5, y.intersp=1.2, cex=0.75, seg.len=1.5)
dev.off()