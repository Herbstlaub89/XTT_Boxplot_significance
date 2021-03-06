###################### Boxplot with significance levels ########################
# This script generates a boxplot for plate reader data in the 96-well format  #
# to compare the fold change between irradiated and control samples.           #
# A stutent's t-test is done and p-value levels between groups are shown       #
# with stars over the groups (* p < 0.05, ** p < 0.01, *** p < 0.001).         #
# To prepare the data for imput please use "Excel_Export_R.xlsx".              #
# The constants in the first part control the basics of the plot (like the     #    
# titel but there are also options to exclude columns and rows of the plates). #
# Please change these values to meet your needs.                               #
###################### Change values here to fit your desire ###################

# Scaling of y-axis
ymini = 0.5
ymaxi = 1.5

# Rows and columns not to use 
# (eg. borders of plate or border between irradiated/control)
skiprows <- c("A", "H") 
skipcolumns <- c(1,6,7,12) 

# Main-, sub- and axis titel of plot
maintitel <- "XTT assay of ADSC´s"
subtitel <- "transfered medium \nirradiated 30 min (453 nm, 23 mW/cm)"
xaxis <- "Time since medium transfer [h]"
yaxis <- "Fold change (irradiated/control)"

# Labels for the ticks on the x axis 
xaxisticks <- c("CTRL","24","48", "72")

### Balance group sizes? 
# Otherwise a controlgroup with more observations can lead 
# to overestimation of significance
balanceGroups = TRUE

### Outliers-test, set to 0 to disable
# maximal outliers to remove in each direction and group
# (highest AND lowest value are checked)
maxoutliers <- 1

# smaller values increase sensitivity of outlier detection
# if highest value > (75% quantile + iqr * threshold) -> outlier
# if lowest value < (25% quantile - iqr * threshold) -> outlier
# a threshold of 1.5 is a good starting point 
threshold <- 1

### Show histogram as diagnosis plot?
# useful to manualy detect outliers and afterwards set the maxoutliers variable accordingly
showHist = TRUE

### Show jittered points
# shows individual data points in plot
showJitter = FALSE

###################### End of setup section ####################################

# Set Random seed (mostly affects the group balancing)
set.seed(42)

CTRLcols = c(7,8,9,10,11,12)

# Load ggplot2 for plotting, readxl for excel import and dplyr for cleaning of data  
library(ggplot2)
library(readxl)
library(dplyr) 

#### Plot theme ####

xtheme <- function() {
  theme(
    # legend
    legend.position = "right", legend.title = element_text(face = "bold", colour = "#000000", size = 10),
    legend.background = element_rect(fill = "#FFFFFF"),
    legend.key = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
    legend.text = element_text(colour = "#000000", size = 10),
    #legend.title=element_blank(),
    # plot background
    plot.title = element_text(colour = "black", face = "bold", size = 18, vjust = 1),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white"),
    # Axis 
    axis.text = element_text(colour = "#000000", size=14),
    axis.text.x = element_text(colour = "#000000", size=12),
    axis.title = element_text(colour = "black", face = "bold", size = 12),
    axis.title.y = element_text(colour = "black", face = "bold", size = 12, vjust=1),
    axis.ticks = element_line(colour = "black"),
    # panel
    panel.grid.major.x = element_line(colour = "#FFFFFF"),
    panel.grid.minor.x = element_line(colour = "#FFFFFF"), #element_blank(),
    panel.grid.major.y = element_line(colour = "#FFFFFF"), #element_blank(),
    panel.grid.minor.y = element_line(colour = "#FFFFFF"), #element_blank(),
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(fill = "#333333")
  )
}

#### Load data ####
XTTpath <- file.choose()
XTTraw <- read_excel(XTTpath, 2)

#### Name Columns correctly #### 
for(i in grep("lab", names(XTTraw), ignore.case=T)) {
  if(any(XTTraw[i] == LETTERS[1:8])) {
    names(XTTraw)[i] <- "Row"
  } else if(any(XTTraw[i] == c(1:12))) {
    names(XTTraw)[i] <- "Column"      
  }
}

XTTraw <- XTTraw %>%
  mutate(Harvest = ifelse(Column %in% CTRLcols,
                          yes = 0,
                          no = `Harvesting time (h)`)
  )

names(XTTraw)[grep("data", names(XTTraw), ignore.case=T)] <- "Data"


#### Clean data #### 
XTTclean <- XTTraw %>%
  filter(Data != 0 &
           !is.na(Data) &
           !(Row %in% skiprows) &
           !(Column %in% skipcolumns))

#### Set types of columns correctly #### 
XTTclean$Harvest <- as.factor(XTTclean$Harvest)
XTTclean$Column <- as.factor(XTTclean$Column)
XTTclean$Row <- as.factor(XTTclean$Row)

#### add row index for referencing ###
XTTclean$id <- c(1:nrow(XTTclean))

#### get number of measurments #### 
lvls <- levels(XTTclean$Harvest)

#### outliers test for each level ####
reflist <- XTTclean[,c("id", "Data", "Harvest")] # list containing only needed data for outlier test
outlCount <- 0
idlist <- c()
for(lvl in lvls) {
  outliers <- 0
  while(outliers < maxoutliers ) {
    outliers <- outliers +1
    subref <- reflist[reflist$Harvest == lvl & !(reflist$id %in% idlist),]
    subref <- subref[order(subref$Data),]
    iqr <- IQR(subref$Data)
    if(subref[1,2] < (quantile(subref$Data, 0.25) - iqr*threshold)) {
      idlist <- append(idlist, subref$id[1])
      outlCount <- outlCount + 1
    }
    if(subref[nrow(subref),2] > (quantile(subref$Data, 0.75) + iqr*threshold)) {
      idlist <- append(idlist, subref$id[nrow(subref)])
      outlCount <- outlCount + 1
    }
  }
}
XTTout <- XTTclean[!(XTTclean$id %in% idlist),]

#### calculate averages for control groups ####
avg <- XTTout %>% 
  group_by(Harvest, PlateID) %>% 
  summarise(Avg = mean(Data))

#### calculate FoldChanges ####
suppressWarnings( # suppress warning because of uninitialised column: 'FoldChange' 
  for(i in 1:nrow(XTTout)) {
    XTTout$FoldChange[i] <- XTTout$Data[i]/as.numeric(avg$Avg[
      avg$Harvest == 0 &
        avg$PlateID == XTTout$PlateID[i]])
  }
)

#### generate labels for number of measurments #### 
lvlslabel <- c() #reset vector
for(i in lvls) {
  lvlslabel <- append(
    lvlslabel, 
    nrow(subset(XTTout, Harvest == i))
    )
}

#### Balance group-sizes #####
if(balanceGroups) {
  # get size of smallest group
  gmin <- min(lvlslabel) 
  # random sample from each group
  XTTplot <- XTTout %>% 
    group_by(Harvest) %>%
    sample_n(gmin)
} else {
  XTTplot <- XTTout
}

#### t-test to get significant differences #### 
tt <- pairwise.t.test(x = XTTplot$FoldChange, 
                      g = XTTplot$Harvest)

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

####  Get maximal FC value and add offset as basis for postion of lines #### 
MaxFC <- max(XTTplot$FoldChange) + 0.325 

####  Check if MaxFC exceeds upper point of y-axis, set it back if needed #### 
if(MaxFC > (ymaxi + 0.05)) MaxFC <- (ymaxi + 0.05)

####  Postion of number of measurments #### 
nlabel.df <- data.frame(Harvest = lvls,
                        FoldChange = ymini)

######### Plot ################
# Set plot size
windows(8,8)

#####  Basic setup of plot #### 
p <- ggplot(XTTplot, aes(Harvest, FoldChange)) + 
  geom_boxplot(outlier.shape = 3, aes(group = Harvest)) +
  coord_cartesian(ylim = c(ymini, ymaxi)) +
  xtheme()+
  ggtitle(maintitel, subtitel) + ylab(yaxis) + xlab(xaxis) +
  #geom_text(data = nlabel.df, label = paste("n = ",lvlslabel)) + #number of observations/group
  #does not work if balanceGroups is TRUE because original number is shown
  scale_y_continuous(breaks=seq(ymini,ymaxi,0.1)) + 
  scale_x_discrete(labels= xaxisticks)  # uncomment to manualy set x-axis ticks

######### postion of  lines and stars ######### 
# create list of dataframes for coordinates of lines
lines <- vector("list", nrow(stars))
starlabel <- data.frame("label" = integer(), "x" = integer(), "y" = integer())

# loop over stars dataframe and generate coordinates for lines
for(i in seq_along(stars[,1])) {
  lines[[i]] <- data.frame(xli =  as.numeric(c(rep(stars$poscol[i],2),
                                               rep(stars$posrow[i],2))), 
                           # i*step size, 
                           # gap between lines, 
                           # offset (horizontal length of line)
                           yli = c(MaxFC - 2*i*0.025, 
                                   MaxFC - 2*i*0.025+0.025,  
                                   MaxFC - 2*i*0.025+0.025, 
                                   MaxFC - 2*i*0.025))
  
  # loop over stars dataframe and generate coordinates for stars, 
  # x is between compared groups(mean), 
  # and y depends on MaxFC and i 
  starlabel[i,] <- c(stars$value[i], 
                     mean(c(as.numeric(stars$poscol[i]), 
                            as.numeric(stars$posrow[i]))), 
                     (MaxFC - i*0.025-i*0.025+0.035))
}

##### draw lines #### 
for(i in seq_along(lines)) {
  p <- p + geom_line(data = lines[[i]], aes(x = xli, y = yli))
} 

#### draw stars, show plot #### 
p <- p + geom_text(data = starlabel, 
                   aes(label = label, 
                       x = as.numeric(x), 
                       y = as.numeric(y)))


#### Jittered points ####
if(showJitter) {
 p <- p + geom_point(shape = 1, position = position_jitter(w = 0.1, h = 0))
}

#### Histogram as diagnosis plot ####
if(showHist) {
  hist <- ggplot(XTTplot, 
         aes(x = FoldChange)) +
    geom_histogram(alpha = 0.8, binwidth =  0.05, color = "black") +
    facet_grid(Harvest~.) +
    xtheme() +
    scale_x_continuous(breaks=seq(ymini,ymaxi,0.1)) +
    geom_vline(xintercept = 1)
  gridExtra::grid.arrange(p, hist, ncol = 2)
  g <- gridExtra::arrangeGrob(p, hist, ncol = 2)
} else {
  p
}

XTTplot %>% 
  group_by(Harvest) %>% 
  summarise(meanFC = mean(FoldChange), 
            sdFC = sd(FoldChange),
            meanRaw = mean(Data),
            sdRaw = sd(Data),
            n = n())

if(outlCount > 0) {
  cat(paste(sep = "","\n", outlCount, " outliers were found and removed.\n",
                       "The following list contains observations considered outliers:\n\n"))
  data.frame(XTTclean[idlist,])
}
ggsave(paste(tools::file_path_sans_ext(XTTpath),".png", sep =""),plot = p,height=12,width=12,dpi=600,units="cm")
