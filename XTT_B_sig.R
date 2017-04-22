setwd("G:/Light/Thomas")
XTTresults <- read.csv("XTT-ASC-100317-7.5min-24_48h.csv")
XTTclean <- subset(XTTresults, Data != 0 & Row != "A" & Column != c(1,6,7,12))
windows(4,4)
p <- ggplot(XTTclean, aes(factor(Harvesting.time..h.), Fold.Change)) + 
  geom_boxplot(outlier.shape = 3, aes(group = Harvesting.time..h.)) +
  coord_cartesian(ylim = c(0.7, 1.5)) +
  theme_bw() + theme(panel.grid = element_blank()) 
p + geom_jitter(width = 0.1, height = 0)



'label.df <- data.frame(Group = c("S1", "S2"),
Value = c(6, 9))

p + geom_text(data = label.df, label = "***")'