#Tidy up workspace just for debugging
#rm(list =ls())

###################### Boxplot with significance levels ########################
# This script generates a boxplot for plate reader data in the 96-well format  #
# to compare the fold change between irradiated and control samples.           #
# A stutent's t-test is done and p-value levels between groups are shown       #
# with stars over the groups. (* p < 0.05, ** p< 0.01, *** p< 0.001)           #
# Important: Only graphs with 3 or 4 groups (eg. Control, 24h, 48h and 72h)    #
# can be processed currently!                                                  #
# To prepare the data for imput please use Excel_Export_R.xlsx                 #
# The constants in the next part, control the basics of the plot (like the     #    
# titel but there are also options to exclude columns and rows of the plates). #
# Please change these values to meet your needs.                               #
################################################################################


###################### Change values here to fit your desire ###################

# y-scaling
ymini = 0.5 
ymaxi = 1.5

# Rows and columns not to use 
# (eg. boarders of plate or boarder between irradiated/control)
skiprows <- c("A", "H")
skipcolumns <- c(1,6,7,12)

# Main-, sub- and axis titel of plot
maintitel <- "XTT assay"
subtitel <- "ADSC irratiated 30min (453nm, 23mW/cm?), 3 plates each time point"
xaxis <- "Time since irradiation [h]"
yaxis <- "Fold change (irradiated/control)"

# Labels for the ticks on the x axis 
xaxisticks <- c("CTRL","24","48", "72")

###################### End of setup section ####################################
################################################################################

# Load ggplot2 for plotting and readxl for excel import
library("ggplot2")
library("readxl")

# Load data
XTTpath <- file.choose()
XTTresults <- read_excel(XTTpath, 2)

# Name Columns corretly
for(i in grep("lab", names(XTTresults), ignore.case=T)) {
  if(any(XTTresults[i] == LETTERS[1:8])) {
    names(XTTresults)[i] <- "Row"
  } else if(any(XTTresults[i] == c(1:12))) {
    names(XTTresults)[i] <- "Column"      
  }
}

names(XTTresults)[grep("data", names(XTTresults), ignore.case=T)] <- "Data"
names(XTTresults)[grep("harvest", names(XTTresults), ignore.case=T)] <- "Harvest"
names(XTTresults)[grep("fold", names(XTTresults), ignore.case=T)] <- "Fold.Change"

# Clean data
XTTclean1 <- subset(XTTresults, Data != 0)
XTTclean2 <- subset(XTTclean1, !is.na(Data))
XTTclean3 <- subset(XTTclean2, Row != skiprows) #uncomment to exclude righ and left  wells
sel <- XTTclean3[,"Column"] %in% skipcolumns
XTTclean <- XTTclean3[!sel,] #uncomment to exclude top and bottom and border wells

# Set types of columns correctly
XTTclean$Harvest <- as.factor(XTTclean$Harvest)
XTTclean$Column <- as.factor(XTTclean$Column)
XTTclean$Row <- as.factor(XTTclean$Row)

if(is.factor(XTTclean$Fold.Change)){
  XTTclean$Fold.Change <- as.numeric(levels(XTTclean$Fold.Change))[XTTclean$Fold.Change]
} else {
  XTTclean$Fold.Change <- as.numeric(XTTclean$Fold.Change)
}

# Test significance
tt <- pairwise.t.test(x = XTTclean$Fold.Change, 
                      g = XTTclean$Harvest)

# Assign stars according to significance level
get.stars <- function(pvalue) {
  strs <- c()
  if(is.na(pvalue)) {
    strs <- ""
  } else {
    if(pvalue >= 0.05) {
      strs <- "" 
    } else if(pvalue < 0.05 & pvalue >= 0.01) {
      strs <- "*"
    } else if(pvalue < 0.01 & pvalue >= 0.001) {
      strs <- "**"
    } else {
      strs <- "***"
    }
  }
  return(strs)
}

stars <- c() #reset vector

for(i in seq_along(tt$p.value)) {
  if(!is.na(tt$p.value[i])) {
    stars <- append(stars, get.stars(tt$p.value[i]))
  }
  
}

# get number of measurments
lvls <- levels(XTTclean$Harvest)

# generate labels for number of measurments
lvlslabel <- c() #reset vector
for(i in lvls) {
  lvlslabel <- append(lvlslabel, paste("n = ", nrow(subset(XTTclean, 
                                                           Harvest == i))))
}

######### Positioning ######### 
# Get maximal FC value and add offset as basis for postion of lines
MaxFC <- max(XTTclean$Fold.Change) + 0.325 
# Check if MaxFC exceeds upper point of y-axis
# set it back if needed
if(MaxFC > ymaxi) MaxFC <- ymaxi

# Postition of stars
if(length(lvls) == 3) {
  slabel.df <- data.frame(Harvest = c(1.5, 2, 2.5),
                          Fold.Change = c(
                            MaxFC - 0.05, 
                            MaxFC, 
                            MaxFC - 0.1 )) 
} else if(length(lvls) == 4) {
  slabel.df <- data.frame(Harvest = c(1.5, #1 postion of stars from top to bottom
                                      2,  #2
                                      2.5, #3
                                      2.5, #4                                  
                                      3, #5
                                      3.5), #6
                          Fold.Change = c(MaxFC - 0.1, #1
                                          MaxFC - 0.05, #2
                                          MaxFC, #3
                                          MaxFC - 0.2, #4
                                          MaxFC - 0.15, #5                                          
                                          MaxFC - 0.25)) #6
}

# Postion of number of measurments
nlabel.df <- data.frame(Harvest = lvls,
                        Fold.Change = ymini)

# Postition of lines
if(length(lvls) == 3) {
  df1 <- data.frame(a = c(1, 1, 3, 3), 
                    b = c(MaxFC - 0.05, 
                          MaxFC - 0.025,  
                          MaxFC - 0.025, 
                          MaxFC - 0.05))
  df2 <- data.frame(a = c(1, 1, 2, 2), 
                    b = c(MaxFC - 0.1, 
                          MaxFC - 0.075, 
                          MaxFC - 0.075, 
                          MaxFC - 0.1))
  df3 <- data.frame(a = c(2, 2, 3, 3), 
                    b = c(MaxFC - 0.15, 
                          MaxFC - 0.125, 
                          MaxFC - 0.125, 
                          MaxFC - 0.15))
} else if(length(lvls) == 4) {
  df1 <- data.frame(a = c(1, 1, 4, 4), 
                    b = c(MaxFC - 0.05, 
                          MaxFC - 0.025, 
                          MaxFC - 0.025, 
                          MaxFC - 0.05))
  df2 <- data.frame(a = c(1, 1, 3, 3), 
                    b = c(MaxFC - 0.1, 
                          MaxFC - 0.075, 
                          MaxFC - 0.075, 
                          MaxFC - 0.1))
  df3 <- data.frame(a = c(1, 1, 2, 2), 
                    b = c(MaxFC - 0.15, 
                          MaxFC - 0.125, 
                          MaxFC - 0.125, 
                          MaxFC - 0.15))
  df4 <- data.frame(a = c(2, 2, 4, 4),
                    b = c(MaxFC - 0.2, 
                          MaxFC - 0.175, 
                          MaxFC - 0.175, 
                          MaxFC - 0.2))
  df5 <- data.frame(a = c(2, 2, 3, 3),
                    b = c(MaxFC - 0.25, 
                          MaxFC - 0.225, 
                          MaxFC - 0.225, 
                          MaxFC - 0.25))
  df6 <- data.frame(a = c(3, 3, 4, 4),
                    b = c(MaxFC - 0.3, 
                          MaxFC - 0.275, 
                          MaxFC - 0.275, 
                          MaxFC - 0.3))
  dfbig <- list(df3, df2, df1, df4, df5, df6)
}

######### Plot ################
# Set plot size
windows(8,8)

# Basic setup of plot
p <- ggplot(XTTclean, aes(Harvest, Fold.Change)) + 
  geom_boxplot(outlier.shape = 3, aes(group = Harvest)) +
  coord_cartesian(ylim = c(ymini, ymaxi)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  ggtitle(maintitel, subtitel) + ylab(yaxis) + xlab(xaxis) +
  geom_jitter(width = 0.1, height = 0) +
  scale_y_continuous(breaks=seq(ymini,ymaxi,0.1)) +
  geom_text(data = nlabel.df, label = lvlslabel)   # uncomment to lable number of observations

# Draw in datapoints as second layer and lines with stars
# Check number of groups (3 or 4)
if(length(lvls) == 3) {
  p + 
    geom_line(data = df1, aes(x = a, y = b)) + 
    geom_line(data = df2, aes(x = a, y = b)) +  
    geom_line(data = df3, aes(x = a, y = b)) +
    scale_x_discrete(labels=xaxisticks[1:3]) + # uncomment to manualy set x-axis ticks 
    geom_text(data = slabel.df, label = stars[1:3]) #Stars 1st: 1vs2, 2nd: 1vs3 3rd: 2vs3
} else if(length(lvls) == 4) {
  # basic setup
  p <- p + 
    scale_x_discrete(labels= xaxisticks) + # uncomment to manualy set x-axis ticks 
    geom_text(data = slabel.df, label = stars[1:6])
  # draw only lines for significant differences 
  for(i in seq_along(stars)) {
    if(stars[i] != "") {
      p <- p + geom_line(data = dfbig[[i]], aes(x = a, y = b))
    }
  }
  p
} 
