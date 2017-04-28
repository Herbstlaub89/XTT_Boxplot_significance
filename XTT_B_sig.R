###################### Boxplot with significance levels ########################
# This script generates a boxplot for plate reader data in the 96-well format  #
# to compare the fold change between irradiated and control samples.           #
# A stutent's t-test is done and p-value levels between groups are shown       #
# with stars over the groups (* p < 0.05, ** p< 0.01, *** p< 0.001).           #
# To prepare the data for imput please use "Excel_Export_R.xlsx".              #
# The constants in the first part control the basics of the plot (like the     #    
# titel but there are also options to exclude columns and rows of the plates). #
# Please change these values to meet your needs.                               #
###################### Change values here to fit your desire ###################

# Scaling of y-axis
ymini = 0.5 
ymaxi = 1.5

# Rows and columns not to use 
# (eg. boarders of plate or boarder between irradiated/control)
skiprows <- c("A", "H")
skipcolumns <- c(1,6,7,12)

# Main-, sub- and axis titel of plot
maintitel <- "XTT assay"
subtitel <- "ADSC irratiated 30min (453nm, 23mW/cm), 3 plates each time point"
xaxis <- "Time since irradiation [h]"
yaxis <- "Fold change (irradiated/control)"

# Labels for the ticks on the x axis 
xaxisticks <- c("CTRL","24","48", "72")

# Outliers-test, set to 0 to disable
maxoutliers <- 1 # maximal outliers to remove in each direction
# (highest AND lowest value are checked!)

senivity <- 1 # smaller values increase senitivity of outlier detection
# value of 1 means if highest value > second highst value + iqr*factor: outlier
# the same applies for the smallest value

###################### End of setup section ####################################

####  Load ggplot2 for plotting and readxl for excel import #### 
library("ggplot2")
library("readxl")

#### Load data ####
XTTpath <- file.choose()
XTTresults <- read_excel(XTTpath, 2)

#### Name Columns corretly #### 
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

#### Clean data #### 
XTTclean1 <- subset(XTTresults, Data != 0)
XTTclean2 <- subset(XTTclean1, !is.na(Data))
XTTclean3 <- subset(XTTclean2, Row != skiprows) #uncomment to exclude righ and left  wells
sel <- XTTclean3[,"Column"] %in% skipcolumns
XTTclean <- XTTclean3[!sel,] #uncomment to exclude top and bottom and border wells

#### Set types of columns correctly #### 
XTTclean$Harvest <- as.factor(XTTclean$Harvest)
XTTclean$Column <- as.factor(XTTclean$Column)
XTTclean$Row <- as.factor(XTTclean$Row)

if(is.factor(XTTclean$Fold.Change)){
  XTTclean$Fold.Change <- as.numeric(levels(XTTclean$Fold.Change))[XTTclean$Fold.Change]
} else {
  XTTclean$Fold.Change <- as.numeric(XTTclean$Fold.Change)
}

#### add row index for referencing ###

XTTclean$id <- c(1:nrow(XTTclean))

#### get number of measurments #### 
lvls <- levels(XTTclean$Harvest)

#### outliers test for each level ####
reflist <- XTTclean[,c("id", "Fold.Change", "Harvest")]

for(lvl in lvls) {
  outliers <- 0
  while(outliers < maxoutliers ) {
    idlist <- c()
    outliers <- outliers +1
    subref <- reflist[reflist$Harvest == lvl,]
    subref <- subref[order(subref$Fold.Change),]
    iqr <- IQR(subref$Fold.Change)
    if(subref[1,2] < (subref[2,2] - iqr*senivity)) {
      idlist <- append(idlist, subref$id[1])
    }
    if(subref[nrow(subref),2] > (subref[nrow(subref)-1,2] + iqr*senivity)) {
      idlist <- append(idlist, subref$id[nrow(subref)])
    }
    XTTclean <- XTTclean[!(XTTclean$id %in% idlist),]
  }
}


#### t-test to get significant differences #### 
tt <- pairwise.t.test(x = XTTclean$Fold.Change, 
                      g = XTTclean$Harvest)

#### function that assigns stars according to significance level #### 
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

#### generate dataframe from p.values of t-test #### 
stars <- data.frame("col" = integer(), 
                    "row" = integer(), 
                    "value" = integer(), 
                    "poscol" = integer(), 
                    "posrow" = integer()) #reset vector

for(i in colnames(tt$p.value)) {
  for(j in rownames(tt$p.value)) {
    
    stars[nrow(stars)+1,] <- c(i,j,get.stars(tt$p.value[j,i]),which(i==lvls),which(j==lvls)) 
  }   
}
#### drop entries that are not significant #### 
stars <- stars[stars$value != "",]

#### generate labels for number of measurments #### 
lvlslabel <- c() #reset vector
for(i in lvls) {
  lvlslabel <- append(lvlslabel, 
                      paste("n = ", 
                            nrow(subset(XTTclean, Harvest == i))))
}

####  Get maximal FC value and add offset as basis for postion of lines #### 
MaxFC <- max(XTTclean$Fold.Change) + 0.325 

####  Check if MaxFC exceeds upper point of y-axis, set it back if needed #### 
if(MaxFC > ymaxi) MaxFC <- ymaxi

####  Postion of number of measurments #### 
nlabel.df <- data.frame(Harvest = lvls,
                        Fold.Change = ymini)

######### Plot ################
# Set plot size
windows(8,8)

#####  Basic setup of plot #### 
p <- ggplot(XTTclean, aes(Harvest, Fold.Change)) + 
  geom_boxplot(outlier.shape = 3, aes(group = Harvest)) +
  coord_cartesian(ylim = c(ymini, ymaxi)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  ggtitle(maintitel, subtitel) + ylab(yaxis) + xlab(xaxis) +
  # geom_jitter(width = 0.1, height = 0) + #uncomment to see individual datapoints
  scale_y_continuous(breaks=seq(ymini,ymaxi,0.1)) + 
  scale_x_discrete(labels= xaxisticks)  # uncomment to manualy set x-axis ticks
  # geom_text(data = nlabel.df, label = lvlslabel)   # uncomment to lable number of observations
  
  
  ######### postion of  lines and stars ######### 
# create list of dataframes for coordinates of lines
lines <- vector("list", nrow(stars))
starlabel <- data.frame("label" = integer(), "x" = integer(), "y" = integer())
# loop over stars dataframe and generate coordinates for lines
for(i in seq_along(stars[,1])) {
  lines[[i]] <- data.frame(xli =  as.numeric(c(rep(stars$poscol[i],2),
                                               rep(stars$posrow[i],2))), 
                           # i*step size, gap between lines, offset (horizontal length of line)
                           yli = c(MaxFC - 2*i*0.025, 
                                   MaxFC - 2*i*0.025+0.025,  
                                   MaxFC - 2*i*0.025+0.025, 
                                   MaxFC - 2*i*0.025))
  # loop over stars dataframe and generate coordinates for stars, x is between compared groups(mean), and y depends on MaxFC and i 
  starlabel[i,] <- c(stars$value[i], mean(c(as.numeric(stars$poscol[i]), as.numeric(stars$posrow[i]))), (MaxFC - i*0.025-i*0.025+0.035))
  
}

##### draw lines #### 
for(i in seq_along(lines)) {
  p <- p + geom_line(data = lines[[i]], aes(x = xli, y = yli))
} 

####  draw stars #### 
p + geom_text(data = starlabel, aes(label = label, x = as.numeric(x), y = as.numeric(y)))

#### Histogram as diagnosis plot ####
windows(8,8)
ggplot(XTTclean, aes(x = Fold.Change, color = Harvest, fill = Harvest)) + geom_histogram(alpha = 0.8, binwidth =  0.05) + facet_grid(Harvest~.)
