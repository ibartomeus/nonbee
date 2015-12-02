#this file will have 4 analysis:
#1) % visits per guild across crops 
#2) Efficiency of bees/non bees
#3) Fruit set in relation with each guild visitation
#4) General relationship of each guild with landscape. 

#read data and load packages needed----
library(reshape)
library(dplyr)
library(nlme)
library(lme4)
library(multcomp)
#library(plotrix)

vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }

d <- read.table("Data/Site_data2.txt", h = TRUE)
str(d)

sv <- read.table("Data/Efficiency_data2.txt")
str(sv)

#1) % visits per guild across crops----
str(d)
#check how many (%) each group represents
length(which(rowSums(table(d$author_crop, d$ant_percent)[,-1]) != 0)) #12 /37
length(which(rowSums(table(d$author_crop, d$other_hymenoptera_percent)[,-1]) != 0)) #17
length(which(rowSums(table(d$author_crop, d$coleoptera_percent)[,-1]) != 0)) #25
length(which(rowSums(table(d$author_crop, d$other_diptera_percent)[,-1]) != 0)) #26
length(which(rowSums(table(d$author_crop, d$hemiptera_percent)[,-1]) != 0)) #9
length(which(rowSums(table(d$author_crop, d$lepidoptera_percent)[,-1]) != 0)) #19
length(which(rowSums(table(d$author_crop, d$other_percent)[,-1]) != 0)) #9
length(which(rowSums(table(d$author_crop, d$syrphidae_percent)[,-1]) != 0)) #31
#We decide to group at bee/ honey bee / non bee level.

#As suggested by Lucas, we use the lnme framework, to avoid treating all systems 
# as equal in the boxplots, regardless of sample size.

m_v_nb <- lme(nonbee_percent ~ 1, random = ~1|author_crop, data = d)
summary(m_v_nb)
i_mvnb <- intervals(m_v_nb, level=0.95, which="fixed")

m_v_a <- lme(apis_percent ~ 1, random = ~1|author_crop, data = d)
summary(m_v_a)
i_mva <- intervals(m_v_a, level=0.95, which="fixed")

m_v_b <- lme(bee_percent ~ 1, random = ~1|author_crop, data = d)
summary(m_v_b)
i_mvb <- intervals(m_v_b, level=0.95, which="fixed")
str(i_mvb)

#Fig 2 (See Below for its inclusion in a panel)
#pdf("Figures/Final_final/figure2.pdf", width=4.6, height=4.6)
#dotchart(c(i_mvnb$fixed[1,2], i_mvb$fixed[1,2], i_mva$fixed[1,2]),
#         #labels = c("non-bee", "other bees", "honey bee"), xlim = c(15,50),
#         labels = c("", "", ""), xlim = c(15,50),
#         xlab = "% visits")

#lines(x = c(i_mva$fixed[1,1],i_mva$fixed[1,3]), y = c(3,3))
#lines(x = c(i_mvb$fixed[1,1],i_mvb$fixed[1,3]), y = c(2,2))
#lines(x = c(i_mvnb$fixed[1,1],i_mvnb$fixed[1,3]), y = c(1,1))
#mtext("B", side = 3, adj = 0)
#invisible(dev.off())

#Fig 2
#excluded one site with NO visitors only to calculate percentages
d2 <- d[-which(d$author_crop == "smitha_Coffea_canephora" & d$bee_percent == 0),]

dat <- d2[,c("author_crop", "ant_percent", "apis_percent", "bee_percent",
             "other_hymenoptera_percent", "coleoptera_percent", "other_diptera_percent",
             "hemiptera_percent", "lepidoptera_percent", "other_percent", 
             "syrphidae_percent")] %>%
  group_by(author_crop) %>%
  summarise_each(funs(mean))

head(dat)
rowSums(dat[,-1])

dat2 <- t(as.matrix(as.data.frame(dat[,-1])))
colnames(dat2) <- dat$author_crop
dat2
dat2 <- dat2[c("apis_percent", "bee_percent", "syrphidae_percent",
               "other_diptera_percent", "other_hymenoptera_percent",
               "lepidoptera_percent", "coleoptera_percent", "hemiptera_percent",
               "ant_percent", "other_percent"),]
dat2 <- dat2[, order(dat2[1,]+dat2[2,], decreasing = TRUE)]

#change col names for plotting
original <- colnames(dat2)
colnames(dat2) <- c("Lowland coffee", "Grapefruit",              
                    "Watermelon A", "Oilseed rape A",      
                    "Watermelon B", "Field bean",                   
                    "Sunflower A", "Highland coffee",               
                    "Strawberry", "Kiwi",              
                    "Apple A", "Sunflower B",        
                    "Oilseed rape B", "Buckwheat A",    
                    "Almond", "Pear",         
                    "Oilseed rape C", "Oilseed rape D",
                    "Oilseed rape E", "Onion (seed)",                  
                    "Oilseed rape F", "Oilseed rape G",             
                    "Turnip rape J", "Apple B",                        
                    "Oilseed rape K", "Buckwheat B",            
                    "Turnip rape L", "Mango A",               
                    "Oilseed rape M", "Carrot (seed)",                
                    "Cherry", "Mango B",         
                    "Oilseed rape N", "Mango C",
                    "Custard apple A", "Custard apple B", 
                    "Soursop")
coded <- colnames(dat2)
cbind(original, coded)
#Fig 2 in final paper
pdf("Figures/Final_final/figure2.pdf", width=8, height=4.6)
par(mar = c(11,6.1,4.1,6.5))
clrs <- c("darkblue", "blue", "yellow", "cornsilk", "goldenrod",
          "orange", "chocolate", "brown", "darkred", "red")
barplot(dat2, beside = FALSE, ylab = "% visits", cex.names = 0.8
        , las = 2, col = clrs)
abline (h = 0)
par(xpd=TRUE)
leg <- c("honey bee", "other bee", "Syrphidae", "other Diptera", "Hymenoptera",
         "Lepidoptera", "Coleoptera", "Hemiptera", "ants",
         "other")
legend(x = 45, y = 100, legend = rev(leg), cex = 0.5, fill = rev(clrs))
par(mar = c(5,1.1,4.1,2.1))
par(xpd=FALSE)
invisible(dev.off())

#2) Efficiency of bees/non bees----
str(sv)

crop <- levels(sv$author_crop)
names_ <- crop
names_ <- c("Oilseed rape N", "Apple B", "Field bean", "Onion (seed)",
            "Avocado", "Turnip rape L", "Carrot (seed)", "Kiwi", "Radish", "Mango A", "Watermelon A")
#Fig SUp Mat in final paper
pdf("Figures/Final_final/figure5.pdf", width=8, height=8)
par(mfrow = c(3,4))
for(i in 1:length(levels(sv$author_crop))){
  temp <- subset(sv, author_crop == crop[i])
  temp <- droplevels(temp)
  boxplot(temp$Mean ~ temp$guild, las = 2, varwidth = TRUE,
          main = names_[i], ylab = "Effectiveness")
}
par(mfrow = c(1,1))
invisible(dev.off())

#We use also here the nlme framework.
head(sv)
m <- lme(Mean_sc ~ guild, random = ~1|author_crop, 
         data= sv, na.action = na.omit)
plot(m)
summary(m)
anova(m)
#fetch posthoc diferences
con <- glht(m, linfct = mcp(guild = "Tukey"))
summary(con)
summary(con,test = adjusted("hochberg"))

#intervals
i_m <- intervals(m, level=0.95, which="fixed")
i_m

#Fig 4: See below for the figure in pannels
#pdf("Figures/Final_final/figure4.pdf", width=4.6, height=4.6)
#dotchart(c(i_m$fixed[c(2,1,3),2]),
#         #labels = c("non-bee","other bees" , "honey bee"), 
#         labels = c("","" , ""), 
#         xlim = c(-2,2.5),
#         xlab = "Z- per visit effectiveness")

#lines(x = c(i_m$fixed[2,1],i_m$fixed[2,3]), y = c(1,1))
#lines(x = c(i_m$fixed[1,1],i_m$fixed[1,3]), y = c(2,2))
#lines(x = c(i_m$fixed[3,1],i_m$fixed[3,3]), y = c(3,3))
#text(x = 2, y = 1, labels = "a")
#text(x = 2, y = 2, labels = "b")
#text(x = 2, y = 3, labels = "a")
#mtext("C", side = 3, adj = 0)
#invisible(dev.off())

#Initially in Suplement: Visits * eff #end up in main paper
sv5 <- read.table("Data/total_efficiency_z.txt", header = TRUE)

m <- lme(Mean_sc ~ variable, random = ~1|author_crop, 
         data= sv5, na.action = na.omit)
plot(m)
summary(m)
anova(m)
#fetch posthoc diferences
con <- glht(m, linfct = mcp(variable = "Tukey"))
summary(con)
summary(con,test = adjusted("hochberg"))

i_m2 <- intervals(m, level=0.95, which="fixed")
i_m2

#Fig 1----
pdf("Figures/Final_final/figure1f.pdf", width=12, height=6)
#as per joint figure
par(mfrow = c(1,3))
par(mar = c(5.1,4.1,4.1,0.3), cex = 1.4, xpd=FALSE)
dotchart(c(i_m2$fixed[c(3,2,1),2]),
         labels = c("non-bee", "other bees", "honey bee"), 
         xlim = c(-1,1),
         xlab = "Z-total effectiveness")

lines(x = c(i_m2$fixed[3,1],i_m2$fixed[3,3]), y = c(1,1))
lines(x = c(i_m2$fixed[2,1],i_m2$fixed[2,3]), y = c(2,2))
lines(x = c(i_m2$fixed[1,1],i_m2$fixed[1,3]), y = c(3,3))
text(x = 0.9, y = 1, labels = "a")
text(x = 0.9, y = 2, labels = "a")
text(x = 0.9, y = 3, labels = "a")
mtext("A", side = 3, adj = 0)
par(mar = c(5.1,4.1,4.1,2.1), cex = 1.4, xpd=FALSE)
#panelB (from above)
dotchart(c(i_mvnb$fixed[1,2], i_mvb$fixed[1,2], i_mva$fixed[1,2]),
         #labels = c("non-bee", "other bees", "honey bee"), xlim = c(15,50),
         labels = c("", "", ""), xlim = c(10,50),
         xlab = "% visits")

lines(x = c(i_mva$fixed[1,1],i_mva$fixed[1,3]), y = c(3,3))
lines(x = c(i_mvb$fixed[1,1],i_mvb$fixed[1,3]), y = c(2,2))
lines(x = c(i_mvnb$fixed[1,1],i_mvnb$fixed[1,3]), y = c(1,1))
mtext("B", side = 3, adj = 0)
#panelC (from above)
dotchart(c(i_m$fixed[c(2,1,3),2]),
         #labels = c("non-bee","other bees" , "honey bee"), 
         labels = c("","" , ""), 
         xlim = c(-2,2.5),
         xlab = "Z-per visit effectiveness")

lines(x = c(i_m$fixed[2,1],i_m$fixed[2,3]), y = c(1,1))
lines(x = c(i_m$fixed[1,1],i_m$fixed[1,3]), y = c(2,2))
lines(x = c(i_m$fixed[3,1],i_m$fixed[3,3]), y = c(3,3))
text(x = 2, y = 1, labels = "a")
text(x = 2, y = 2, labels = "b")
text(x = 2, y = 3, labels = "a")
mtext("C", side = 3, adj = 0)

invisible(dev.off())

#3) Fruit set in relation with each guild visitation----

#check correlation
plot(d$nonbee_sc ~ d$bee_sc) #not strong...
plot(d$nonbee_sc ~ d$apis_sc) #not strong...
cor(d[,c("bee_sc", "nonbee_sc", "apis_sc")])
cor.test(d$nonbee_sc, d$bee_sc) #so its correlated, but weakly...
cor.test(d$nonbee_sc, d$apis_sc)
cor.test(d$bee_sc, d$apis_sc)

#models
#null model
m0 <- lme(data = d, fruitset_sc ~ ., random = ~1|author_crop, 
         na.action = na.exclude, weights = varIdent(~ 1 | author_crop)) 
plot(m0)
AIC(m0)

m <- lme(data = d, fruitset_sc ~ nonbee_sc + bee_sc + apis_sc, random = ~1|author_crop, 
         na.action = na.exclude, weights = varIdent(~ 1 | author_crop)) 
plot(m)
AIC(m)
summary(m) 
anova(m)
vif.mer(m)
m$coefficients
intervals(m)
m$apVar 

#Fig 3A (not presented because intervals can not bee computed (but they use to be before updating R!))
#see below for equivalent lmer4 model.
#dotchart(c(fixef(m)[2], fixef(m)[3],fixef(m)[4]),
 #        labels = c("non-bees", "other bees", "honey bees"), 
 #        xlim = c(-0.5,0.5),
 #        xlab = "slope (fruit set)")
#ci <- intervals(m)
#lines(x = c(ci[[1]][2,1],ci[[1]][2,3]), y = c(1,1))
#lines(x = c(ci[[1]][3,1],ci[[1]][3,3]), y = c(2,2))
#lines(x = c(ci[[1]][4,1],ci[[1]][4,3]), y = c(3,3))
#abline(v = 0)

#interactions.
m1 <- lme(data = d, fruitset_sc ~ nonbee_sc * bee_sc * apis_sc, random = ~1|author_crop, 
          na.action = na.omit) 
AIC(m1)
summary(m1)
m2 <- lme(data = d, fruitset_sc ~ nonbee_sc * bee_sc + apis_sc, random = ~1|author_crop, 
          na.action = na.omit) 
summary(m2)
AIC(m2)
m3 <- lme(data = d, fruitset_sc ~ nonbee_sc + bee_sc * apis_sc, random = ~1|author_crop, 
          na.action = na.omit) 
AIC(m3)
m4 <- lme(data = d, fruitset_sc ~ nonbee_sc * apis_sc + bee_sc, random = ~1|author_crop, 
          na.action = na.omit) 
AIC(m4)
#best model m

#random slopes---
m5 <- lmer(data = d, fruitset_sc ~ 
             nonbee_sc + apis_sc + bee_sc + (1 + nonbee_sc + apis_sc + bee_sc | author_crop), 
           na.action = na.omit)
AIC(m)
AIC(m5) #worst AIC...
m6 <- lmer(data = d, fruitset_sc ~ 
             nonbee_sc + apis_sc + bee_sc + (1 + bee_sc | author_crop), 
           na.action = na.omit)
AIC(m6)
m7 <- lmer(data = d, fruitset_sc ~ 
             nonbee_sc + apis_sc + bee_sc + (1 + apis_sc | author_crop), 
           na.action = na.omit)
AIC(m7)
m8 <- lmer(data = d, fruitset_sc ~ 
             nonbee_sc + apis_sc + bee_sc + (1 + nonbee_sc | author_crop), 
           na.action = na.omit)
AIC(m8)
#All AIC worst m

#non linear effects for hb
m9 <- lme(data = d, fruitset_sc ~ nonbee_sc + bee_sc + apis_sc + I(apis_sc^2), random = ~1|author_crop, 
          na.action = na.exclude, weights = varIdent(~ 1 | author_crop)) 
plot(m9)
summary(m9)
AIC(m9)


#Check is its the same with lmer 
m <- lmer(data = d, fruitset_sc ~ nonbee_sc + bee_sc + apis_sc + (1|author_crop),
          na.action = na.omit)
ci <- confint(m)
ci
summary(m)$coefficients

#Fig 3A. See below for the pannel)
#pdf("Figures/Final_final/figure6.pdf", width=4.6, height=4.6)
#pdf("Figures/Final_final/figure2.pdf", width=6, height=4.6)
#par(mfrow = c(1,2))
#par(mar = c(5.1,4.1,4.1,0))
#dotchart(summary(m)$coefficients[2:4,1],
#         labels = c("non-bee", "other bees", "honey bee"), 
#         xlim = c(-0.25,0.5),
#         xlab = "slope (fruit set)")

#lines(x = c(ci[4,1],ci[4,2]), y = c(1,1))
#lines(x = c(ci[5,1],ci[5,2]), y = c(2,2))
#lines(x = c(ci[6,1],ci[6,2]), y = c(3,3))
#abline(v = 0)
#mtext("A", side = 3, adj = 0)
#par(mar = c(5.1,4.1,4.1,2.1))
#invisible(dev.off())

#all summaries
summary(m0)
summary(m)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)
summary(m7)
summary(m8)

#4) General relationship of each guild with landscape----

#Final model
m_ <- lme(data = d, bee_sc ~ distance_sc, random = ~1|author_crop, 
         na.action = na.omit,  weights = varFixed(~ distance_sc))
m1_ <- lme(data = d, nonbee_sc ~ distance_sc, random = ~1|author_crop, 
          na.action = na.omit, weights = varFixed(~ distance_sc))
m2_ <- lme(data = d, apis_sc ~ distance_sc, random = ~1|author_crop, 
          na.action = na.omit,  weights = varFixed(~ distance_sc))

plot(m_)
plot(m1_)
plot(m2_) #all plots look nice enough
summary(m_)
summary(m1_)
summary(m2_)

#Fig 3B in final. See below for pannel
#pdf("Figures/Final_final/figure7.pdf", width=4.6, height=4.6)
#dotchart(c(fixef(m1_)[2], fixef(m_)[2],fixef(m2_)[2]),
#         #labels = c("non-bee", "other bees", "honey bee"), 
#         labels = c("", "", ""), 
#         xlim = c(-0.5,0.5),
#         xlab = "slope (isolation from natural areas)")

#ci_ <- intervals(m_)
#ci1_ <- intervals(m1_)
#ci2_ <- intervals(m2_)
#lines(x = c(ci1_[[1]][2,1],ci1_[[1]][2,3]), y = c(1,1))
#lines(x = c(ci_[[1]][2,1],ci_[[1]][2,3]), y = c(2,2))
#lines(x = c(ci2_[[1]][2,1],ci2_[[1]][2,3]), y = c(3,3))
#abline(v = 0)
#mtext("B", side = 3, adj = 0)
#invisible(dev.off())


#Fig 3----
pdf("Figures/Final_final/figure3f.pdf", width=8, height=6)
par(mfrow = c(1,2))
par(mar = c(5.1,4.1,4.1,0.1), cex = 1.4)
dotchart(summary(m)$coefficients[2:4,1],
         labels = c("non-bee", "other bees", "honey bee"), 
         xlim = c(-0.25,0.5),
         xlab = "slope (fruit set)")

lines(x = c(ci[4,1],ci[4,2]), y = c(1,1))
lines(x = c(ci[5,1],ci[5,2]), y = c(2,2))
lines(x = c(ci[6,1],ci[6,2]), y = c(3,3))
abline(v = 0, lty = 2)
mtext("A", side = 3, adj = 0)
par(mar = c(5.1,4.1,4.1,2.1))
#paneB
dotchart(c(fixef(m1_)[2], fixef(m_)[2],fixef(m2_)[2]),
         #labels = c("non-bee", "other bees", "honey bee"), 
         labels = c("", "", ""), 
         xlim = c(-0.6,0.6),
         xlab = "slope (isolation)")

ci_ <- intervals(m_)
ci1_ <- intervals(m1_)
ci2_ <- intervals(m2_)
lines(x = c(ci1_[[1]][2,1],ci1_[[1]][2,3]), y = c(1,1))
lines(x = c(ci_[[1]][2,1],ci_[[1]][2,3]), y = c(2,2))
lines(x = c(ci2_[[1]][2,1],ci2_[[1]][2,3]), y = c(3,3))
abline(v = 0, lty = 2)
mtext("B", side = 3, adj = 0)

invisible(dev.off())