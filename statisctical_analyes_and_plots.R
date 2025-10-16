install.packages("car")    
install.packages("ggpubr")
install.packages("emmeans")
install.packages("multcomp")
install.packages("multcompView")
install.packages("heplots")
install.packages("lsr")
install.packages("forcats")

library(lsr)
library(car)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)
library(ggpubr)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

all.data      <- read_csv("Github/UU_NP_Soil/data/UU_NP_Soil_data")
sand.data     <- all.data %>% filter(substrate == "Sand")
sandsoil.data <- all.data %>% filter(year == "2")

pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(family = "Helvetica", size = 0, face = "italic"),
        panel.border = element_rect(linewidth = 3, fill = NA),
        #     axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.ticks.length = unit(0.25, "cm"))




# Nutrients in sand, two years
shapiro_test_results <- by(sand.data, list(sand.data$trt, sand.data$year, sand.data$spp), function(subgroup) shapiro.test(subgroup$vcmax))
print(shapiro_test_results)

# traits

sand.data<-sand.data %>%
  mutate(
    CcostN    = (c.root*roots.wt/100)/((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    RootC.wt  = (c.root*roots.wt/100),
    PlantN.wt = ((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    Narea     = (n.leaf/100*plant.lma),
    Comb1 = paste(spp, year, sep = " - Year ") 
    )

all.data <- all.data %>%
  mutate(
    CcostN    = (c.root*roots.wt/100)/((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    RootC.wt  = (c.root*roots.wt/100),
    PlantN.wt = ((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    Narea     = (n.leaf/100*plant.lma)
  )

sandsoil.data <- sandsoil.data %>%
  mutate(
    CcostN    = (c.root*roots.wt/100)/((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    RootC.wt  = (c.root*roots.wt/100),
    PlantN.wt = ((n.root*roots.wt/100)+(n.stem*stem.wt/100)+(n.leaf*leaves.wt/100)),
    Narea     = (n.leaf/100*plant.lma),
    Comb = paste(trt, substrate, sep = " - ")
  )

 ########### statistical analyses ###########

       ## Whole-plant nitrogen biomass 

plantNmodel        <- lm(PlantN.wt ~ trt * spp + year, data = sand.data)
plantNmodel_SI     <- lm(PlantN.wt ~ trt * spp + spp * year + trt * year, data = sand.data)


dev.off()
par(mfrow = c(3, 2))
plot(plantNmodel, main = "Plant N")
qqnorm(residuals(plantNmodel))
qqline(residuals(plantNmodel))
densityPlot(residuals(plantNmodel))
shapiro.test(residuals(plantNmodel)) # all normal
outlierTest(plantNmodel)             # all normal

summary(plantNmodel)
Anova(plantNmodel)

summary(plantNmodel_SI)
Anova(plantNmodel_SI)


emmeans(plantNmodel, pairwise~trt)
emmeans(plantNmodel, pairwise~spp)
emmeans(plantNmodel, pairwise~trt * spp)


       ## Root carbon ####

RootCmodel           <- lm(RootC.wt ~ trt * spp + year, data = sand.data)
RootCmodel_SI        <- lm(RootC.wt ~ trt * spp + spp * year + trt * year, data = sand.data)

dev.off()
par(mfrow = c(3, 2))
plot(RootCmodel, main = "Root C")
qqnorm(residuals(RootCmodel))
qqline(residuals(RootCmodel))
densityPlot(residuals(RootCmodel))
shapiro.test(residuals(RootCmodel)) # all normal
outlierTest(RootCmodel)             # all normal

summary(RootCmodel)
Anova(RootCmodel)

summary(RootCmodel_SI)
Anova(RootCmodel_SI)

emmeans(RootCmodel, pairwise~trt)
emmeans(RootCmodel, pairwise~spp)
emmeans(RootCmodel, pairwise~spp*trt)



                ## Carbon cost for nitrogen acquisition 
CcostNmodel        <- lm(CcostN ~ trt * spp + year, data = sand.data)
CcostNmodel_SI     <- lm(CcostN ~ trt * spp + spp * year + trt * year, data = sand.data)

tukey_interaction <- lsmeans(CcostNmodel, pairwise ~ trt:spp, adjust = "tukey")
print(tukey_interaction)


dev.off()
par(mfrow = c(3, 2))
plot(CcostNmodel, main = "C cost N")
qqnorm(residuals(CcostNmodel))
qqline(residuals(CcostNmodel))
densityPlot(residuals(CcostNmodel))
shapiro.test(residuals(CcostNmodel)) # all normal
outlierTest(CcostNmodel)             # all normal

summary(CcostNmodel)
Anova(CcostNmodel)
r.squaredGLMM(CcostNmodel)

summary(CcostNmodel_SI)
Anova(CcostNmodel_SI)

emmeans(CcostNmodel, pairwise~trt)
emmeans(CcostNmodel, pairwise~spp)
emmeans(CcostNmodel, ~ trt * spp)


 



        ## Vcmax ####

Vcmaxmodel      <- lm(vcmax ~ trt * spp + year, data = sand.data)
Vcmaxmodel_SI   <- lm(vcmax ~ trt * spp + spp * year + trt * year, data = sand.data)


dev.off()
par(mfrow = c(3, 2))
plot(Vcmaxmodel, main = "vcmax")
qqnorm(residuals(Vcmaxmodel))
qqline(residuals(Vcmaxmodel))
densityPlot(residuals(Vcmaxmodel))
shapiro.test(residuals(Vcmaxmodel)) # all normal
outlierTest(Vcmaxmodel)             # all normal

summary(Vcmaxmodel)
Anova(Vcmaxmodel)

summary(Vcmaxmodel_SI)
Anova(Vcmaxmodel_SI)


emmeans(Vcmaxmodel, pairwise~trt, type = "response")
emmeans(Vcmaxmodel, pairwise~spp, type = "response")
emmeans(Vcmaxmodel, ~ trt * spp)


    ## Leaf N area ####
Nareamodel   <- lm(Narea ~ trt * spp + year, data = sand.data)
logNareamodel<- lm(log(Narea) ~ trt * spp + year, data = sand.data)
Nareamodel_SI  <- lm(Narea ~ trt * spp + spp * year + trt * year, data = sand.data)


dev.off()
par(mfrow = c(3, 2))
plot(Nareamodel, main = "Leaf N area")
qqnorm(residuals(Nareamodel))
qqline(residuals(Nareamodel))
densityPlot(residuals(Nareamodel))
shapiro.test(residuals(Nareamodel)) # not normal
outlierTest(Nareamodel)             # all normal

dev.off()
par(mfrow = c(3, 2))
plot(logNareamodel, main = "Leaf N area")
qqnorm(residuals(logNareamodel))
qqline(residuals(logNareamodel))
densityPlot(residuals(logNareamodel))
shapiro.test(residuals(logNareamodel)) # all normal
outlierTest(logNareamodel)             # all normal


summary(logNareamodel)
Anova(logNareamodel)
summary(Nareamodel)
Anova(Nareamodel)

shapiro.test(residuals(Nareamodel_SI))
summary(Nareamodel_SI)
Anova(Nareamodel_SI)


emmeans(Nareamodel, pairwise~trt)
emmeans(Nareamodel, pairwise~spp)
emmeans(Nareamodel, ~ trt * spp)


       ## Chi ####
Chimodel     <- lm(chi ~ trt * spp + year, data = sand.data)
Chimodel_SI  <- lm(chi ~ trt * spp + spp * year + trt * year, data = sand.data)


dev.off()
par(mfrow = c(3, 2))
plot(Chimodel, main = "Chi")
qqnorm(residuals(Chimodel))
qqline(residuals(Chimodel))
densityPlot(residuals(Chimodel))
shapiro.test(residuals(Chimodel)) # all normal
outlierTest(Chimodel)             # all normal

summary(Chimodel)
Anova(Chimodel)
summary(Chimodel_SI)
Anova(Chimodel_SI)


emmeans(Chimodel, pairwise~trt)
emmeans(Chimodel, pairwise~spp)
emmeans(Chimodel, ~ trt * spp)


################## cross table emmeans ######

##  Get emmeans by Treatment by Year by Species
emm_rootC  <- emmeans(RootCmodel , ~ trt * year * spp)
emm_plantN <- emmeans(plantNmodel, ~ trt * year * spp)
emm_CcostN <- emmeans(CcostNmodel, ~ trt * year * spp)

emm_vcmax <- emmeans(Vcmaxmodel, ~ trt * year * spp)
emm_Narea <- emmeans(Nareamodel, ~ trt * year * spp)
emm_chi   <- emmeans(Chimodel  , ~ trt * year * spp)


#  Extract Holcus lanatus only
emm_holcus <- subset(emm_species, spp == "Holcus lanatus")
emm_solanum <- subset(emm_species, spp == "Solanum dulcamara")

#  Pairwise comparisons for Holcus lanatus
pairs_holcus <- pairs(emm_holcus, adjust = "Tukey")
pairs_solanum <- pairs(emm_solanum, adjust = "Tukey")

emm_holcus
emm_solanum
pairs_holcus
pairs_solanum

emm_rootC
emm_plantN
emm_CcostN

emm_vcmax
emm_Narea
emm_chi


########################################

lighter_shades <- c( "#80C0C0",  "#C5A289","#008080","#8B4513")

Sand.CcostN.png <- sand.data %>%
  mutate(trt=fct_relevel(trt, "Low", "Medium", "High"))%>%
  ggplot(aes(x=trt, y=CcostN, fill = Comb1))+
  geom_boxplot(position = position_dodge(0.75), outlier.shape = NA) +  # Hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.75), size = 2) + 
  scale_fill_manual(values = setNames(lighter_shades, c("Solanum dulcamara - Year 1", "Holcus lanatus - Year 1","Solanum dulcamara - Year 2", "Holcus lanatus - Year 2"))) +
  pubtheme +
  guides(fill = guide_legend(title = "Species"), color = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y=expression("C cost to acquire N (gC gN "^"-1"~")" ),x="Nutrient treatment")

Sand.plantN.png <- sand.data %>%
  mutate(trt=fct_relevel(trt, "Low", "Medium", "High"))%>%
  ggplot(aes(x=trt, y=PlantN.wt, fill = Comb1))+
  geom_boxplot(position = position_dodge(0.75), outlier.shape = NA) +  # Hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.75), size = 2) + 
  scale_fill_manual(values = setNames(lighter_shades, c("Solanum dulcamara - Year 1", "Holcus lanatus - Year 1","Solanum dulcamara - Year 2", "Holcus lanatus - Year 2"))) +
  pubtheme +
  guides(fill = guide_legend(title = "Species"), color = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y=expression("Whole-plant nitrogen biomass (gN)" ),x="Nutrient treatment")

Sand.RootC.png <- sand.data %>%
  mutate(trt=fct_relevel(trt, "Low", "Medium", "High"))%>%
  ggplot(aes(x=trt, y=RootC.wt, fill = Comb1))+
  geom_boxplot(position = position_dodge(0.75), outlier.shape = NA) +  # Hide outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.75), size = 2) + 
  scale_fill_manual(values = setNames(lighter_shades, c("Solanum dulcamara - Year 1", "Holcus lanatus - Year 1","Solanum dulcamara - Year 2", "Holcus lanatus - Year 2"))) +
  pubtheme +
  guides(fill = guide_legend(title = "Species"), color = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y=expression("Root carbon biomass (gC)" ),x="Nutrient treatment")

# Create the top row with two smaller plots
top_row <- (Sand.plantN.png + theme(legend.position = "none") | 
              Sand.RootC.png + theme(legend.position = "none")) + 
  plot_layout(widths = c(1, 1))  # Equal width for both plots

# Add the wide bottom plot with the legend
bottom_row <- Sand.CcostN.png + 
  theme(legend.position = "right")  # Place legend within the bottom plot

# Combine the rows
final_plot <- top_row /
  bottom_row +  # Combine rows
  plot_layout(heights = c(1, 1.1)) + # Adjust heights and collect the legend
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")  # Increase tag size and make bold
    )
  )

# Save the plot
png("../fig1_Sandcarboncosts.png",
    height = 13, width = 20, units = "in", res = 600)
final_plot + 
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")  # Increase tag size and make bold
    )
  )
dev.off()
tiff("../fig1_Sandcarboncosts.tiff", height = 13, width = 20, units = "in", res = 600, compression = "lzw")
print(final_plot)
dev.off()

pdf("fig1_Sandcarboncosts.pdf", width = 13, height = 20, bg = "white", colormodel = "cmyk")
print(final_plot)  # Ensure the plot is written to the file
dev.off()


o

# Save the plot
png("../fig2_Sandleaftraits.png",
    height = 13, width = 20, units = "in", res = 600)
final_plot + 
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold")  # Increase tag size and make bold
    )
  )
dev.off()
tiff("../fig2_Sandleaftraits.tiff", height = 13, width = 20, units = "in", res = 600, compression = "lzw")
print(final_plot)
dev.off()

##### Sand vs soils. ####



sandsoil.data <- sandsoil.data %>% mutate(Comb = paste(trt, substrate, sep = " - ") )

comb_colorssoil <- c("Sand" = "#8c8c8c",
                     "Microp"        = "#1F77B4", 
                     "Reijerscamp"    = "#D62728"
)

sandsoilvcmax     <-lm(vcmax ~ n.mg*spp + substrate, data= sandsoil.data)
sandsoilvcmax_SI  <-lm(vcmax ~ p.mg*spp + substrate, data= sandsoil.data)

dev.off()
par(mfrow = c(3, 2))
plot(sandsoilvcmax, main = "Vcmax")
qqnorm(residuals(sandsoilvcmax))
qqline(residuals(sandsoilvcmax))
densityPlot(residuals(sandsoilvcmax))
shapiro.test(residuals(sandsoilvcmax)) # all normal
outlierTest(sandsoilvcmax)             # all normal

summary(sandsoilvcmax)
Anova(sandsoilvcmax)

test(emtrends(sandsoilvcmax, ~ spp, "n.mg"))
vcmax.regline <- data.frame(emmeans(sandsoilvcmax, ~spp, "n.mg",
                                       at = list(n.mg = seq(20, 90, 0.1),
                                                 type = "response")))


sandsoilchi     <-lm(chi ~ n.mg*spp + substrate, data= sandsoil.data)

sandsoilchi_SI  <-lm(chi ~ p.mg*spp + substrate, data= sandsoil.data)

dev.off()
par(mfrow = c(3, 2))
plot(sandsoilchi, main = "Chi")
qqnorm(residuals(sandsoilchi))
qqline(residuals(sandsoilchi))
densityPlot(residuals(sandsoilchi))
shapiro.test(residuals(sandsoilchi)) # all normal
outlierTest(sandsoilchi)             # all normal

summary(sandsoilchi)
Anova(sandsoilchi)

test(emtrends(sandsoilchi, ~spp, "n.mg"))
chi.regline <- data.frame(emmeans(sandsoilchi, ~spp, "n.mg",
                                     at = list(n.mg = seq(20, 90, 0.1),
                                               type = "response")))


sandsoilCcostN    <-lm(CcostN ~ n.mg*spp + substrate, data= sandsoil.data)
sandsoilCcostN_SI <-lm(CcostN ~ p.mg*spp + substrate, data= sandsoil.data)

dev.off()
par(mfrow = c(3, 2))
plot(sandsoilCcostN, main = "CcostN")
qqnorm(residuals(sandsoilCcostN))
qqline(residuals(sandsoilCcostN))
densityPlot(residuals(sandsoilCcostN))
shapiro.test(residuals(sandsoilCcostN)) # all normal
outlierTest(sandsoilCcostN)             # one outlier, no changes made

test(emtrends(sandsoilCcostN, ~spp, "n.mg"))
CcostN.regline <- data.frame(emmeans(sandsoilCcostN, ~spp, "n.mg",
                                     at = list(n.mg = seq(20, 90, 0.1),
                                               type = "response")))




sandsoil.data %>%
  mutate(Comb=fct_relevel(Comb, "Low - Sand", "Medium - Sand", "Natural - Microp", "Natural - Reijerscamp", "High - Sand"))%>%
  ggplot(aes(x = spp, y = vcmax)) +
  geom_boxplot(aes(fill = Comb)) +
  scale_fill_manual(values = comb_colorssoil) +
  theme_bw(base_size = 13) + 
  guides(fill = guide_legend(title = "Combination Levels")) +
  labs(y = expression("V"["cmax"]~"(mumol m"^"-2"~"s"^"-1"~")"),
       x = "Species")

sandsoil.data %>%
  mutate(Comb=fct_relevel(Comb, "Low - Sand", "Medium - Sand", "Natural - Microp", "Natural - Reijerscamp", "High - Sand"))%>%
  ggplot(aes(x = spp, y = chi)) +
  geom_boxplot(aes(fill = Comb)) +
  scale_fill_manual(values = comb_colorssoil) +
  theme_bw(base_size = 13) + 
  guides(fill = guide_legend(title = "Combination Levels")) +
  labs(y = expression("Isotope derived " ~chi~" " ),
       x = "Species")

NutrientNatural %>%
  mutate(Comb=fct_relevel(Comb, "Low - Sand", "Medium - Sand", "Natural - Microp", "Natural - Reijerscamp", "High - Sand"))%>%
  ggplot(aes(x = Species, y = CcostN)) +
  geom_boxplot(aes(fill = Comb)) +
  scale_fill_manual(values = comb_colorsNN) +
  theme_bw(base_size = 13) + 
  guides(fill = guide_legend(title = "Combination Levels")) +
  labs(y = expression("Carbon cost to acquire nitrogen (gC gN "^"-1"~")"  ),
       x = "Species")

# --- Vcmax ---

sandsoil.vcmax.png <- ggplot(
  sandsoil.data,
  aes(x = n.mg, y = vcmax, shape = spp, color = substrate)
) +
  # points
  geom_jitter(size = 3, alpha = 0.8, width = 0.6, height = 0) +
  
  # fitted regression lines (from emmeans or lmer predictions)
  geom_smooth(
    data = vcmax.regline,
    aes(x = n.mg, y = emmean, linetype = spp),
    linewidth = 1,
    se = FALSE,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # manual styling
  scale_color_manual(values = comb_colorssoil) +
  scale_shape_manual(values = c(16, 17)) +        # solid circle + triangle
  scale_linetype_manual(values = c("solid", "dashed")) +
  
  # tidy theme & axis limits
  pubtheme +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  
  # labels
  labs(
    x = "Inorganic N (mg / pot)",
    y = expression(V[cmax]~"(µmol m"^{-2}*" s"^{-1}*")"),
    color = "Substrate",
    shape = "Species",
    linetype = "Species",
    tag = "A"
  )


sandsoil.chi.png <- sandsoil.data %>%
  mutate(Comb = fct_relevel(
    Comb, "Low - Sand", "Medium - Sand", "Natural - Microp",
    "High - Sand", "Natural - Reijerscamp"
  )) %>%
  ggplot(aes(x = n.mg, y = chi, shape = spp, color = substrate)) +
  geom_jitter(size = 3, alpha = 0.8, width = 0.6, height = 0) +
  geom_smooth(
    data = chi.regline,
    aes(x = n.mg, y = emmean, linetype = spp),
    linewidth = 1, se = FALSE, color = "black", inherit.aes = FALSE
  ) +
  scale_color_manual(values = comb_colorssoil) +
  scale_shape_manual(values = c(16, 17)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  pubtheme +
  scale_x_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.05))) +
  labs(
    y = "Isotope derived \u03C7",
    x = "Inorganic N (mg / pot)",
    color = "Substrate",
    shape = "Species",
    linetype = "Species",
    tag = "B"
  )


# --- CcostN ---

sandsoil.CcostN.png <- sandsoil.data %>%
  mutate(Comb = fct_relevel(
    Comb, "Low - Sand", "Medium - Sand", "Natural - Microp", "High - Sand", "Natural - Reijerscamp"
  )) %>%
  ggplot(aes(x = n.mg, y = CcostN, shape = spp, color = substrate)) + 
  geom_jitter(size = 3, alpha = 0.8) +
  geom_smooth(
    data = CcostN.regline,
    aes(x = n.mg, y = emmean, linetype = spp),
    linewidth = 1, se = FALSE, color = "black"
  ) +
  scale_color_manual(values = comb_colorssoil) +
  scale_shape_manual(values = c(16, 17)) +
  pubtheme +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    y = "Carbon cost to acquire nitrogen",
    x = "Inorganic N (mg / pot)",
    color = "Substrate",
    shape = "Species",
    linetype = "Species",
    tag = "C"
  )




combined_plot<-ggarrange(
  sandsoil.vcmax.png,
  sandsoil.chi.png,
  sandsoil.CcostN.png, 
  common.legend = TRUE,
  legend = "right",
  ncol = 1, nrow = 3,  # All plots in one row
  labels = c("A", "B", "C"),  # Add labels
  font.label = list(size = 18, face = "bold")
)

combined_plot <- ggarrange(
  sandsoil.vcmax.png,
  sandsoil.chi.png,
  sandsoil.CcostN.png,
  ncol = 1, 
  nrow = 3,
  common.legend = TRUE,
  legend = "right",
  label.x = 0.03,    # Adjust as needed
  label.y = 0.97,    # Adjust as needed
  align = "hv"
)


tiff("../fig3_whole.tiff", height = 13, width = 20, units = "in", res = 600, compression = "lzw")
print(combined_plot)
dev.off()

