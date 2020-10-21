# PSA-ReD Plot Generator v1.0.2. Use this script to make your own PSA-ReD plots.  
# Copyright (C) 2020, Joost Geenen
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version.

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

### ----------------------------------------------------------------------- ###
### --------------  Installing and / or loading libraries. ---------------- ### 


### -------------- Do not change the part in between / below -------------- ###
### ----------------------------------------------------------------------- ###
if (!require('ggplot2')) {
  install.packages("ggplot2")
}
library('ggplot2')

if (!require('MASS')) {
  install.packages("MASS")
}
library('MASS')

if (!require('directlabels')) {
  install.packages("directlabels")
}
library('directlabels')

if (!require('grid')) {
  install.packages("grid")
}
library('grid')

if (!require('gridExtra')) {
  install.packages("gridExtra")
}
library('gridExtra')

if (!require('reshape2')) {
  install.packages("reshape2")
}
library('reshape2')

ProcessContourData <- function(kde.data) {
  # Processess kde data for plotting the contours.  
  # 
  # Args:
  #  kde.data: a list containing the kernel density data,
  #  which is generated using the 'data' which the user loaded. 
  #
  # Returns: 
  #  A dataframe containing the cumulative density per x and y.
  
  kde.dx <- diff(kde.data$x[1:2]) 
  kde.dy <- diff(kde.data$y[1:2]) 
  kde.sz <- sort(kde.data$z)       
  kde_c1 <- cumsum(kde.sz) * kde.dx * kde.dy 
  
  dimnames(kde.data$z) <- list(kde.data$x, kde.data$y)
  kde_dc <- melt(kde.data$z)   
  kde_dc$contour.levels <- approx(kde.sz, 1 - kde_c1, kde_dc$value)$y 
  
  plot.data.contour <- cbind.data.frame(kde_dc[, 1], kde_dc[, 2], kde_dc[, 4])
  colnames(plot.data.contour)[1] <- "x"
  colnames(plot.data.contour)[2] <- "y"
  colnames(plot.data.contour)[3] <- "contour.levels"
  return(plot.data.contour)
}

GenerateNormalisedDensity <- function(kde.data) {
  # Normalises the density to 1 and generates a matrix for plotting the density
  # 
  # Args:
  #  kde.data: a list containing the kernel density data,
  #  which is generated using the 'data' which the user loaded. 
  # 
  # Returns: 
  #  A dataframe containing the relative density per x and y.
  x <- kde.data$x
  y <- kde.data$y
  
  df.density   <- expand.grid(X = x, Y = y)
  df.density$Z <- as.numeric(unlist(kde.data$z))
  
  df.density.normalised   <- df.density
  df.density.normalised$Z <- df.density.normalised$Z * 
    (1 / max(df.density.normalised$Z))
  
  return(df.density.normalised)
}


GeneratePlot <- function(df.density.normalised, 
                         legend.title,
                         plot.data.contour,
                         contour.levels,
                         font.size,
                         font.face,
                         font.family,
                         plot.type,
                         grayscale,
                         x.axis.title,
                         y.axis.title,
                         x.range,
                         y.range,
                         clipping,
                         extend.panel, 
                         WTP.thresholds,
                         basecase,
                         average.PSA) {
  # Generates the ggplot object that can be plotted. 
  #
  # Args:
  #   df.density.normalised: A dataframe consisting of the normalised relative 
  #                          density per x and y.    
  #   legend.title:          The title of the legend. 
  #   plot.data.contour:     A dataframe containing the cumulative density 
  #                          per x and y. 
  #   contour.levels:        A vector containing the specified contour levels. 
  #   font.size:             The font size of characters in the plot
  #   font.face:             The font face (ie, bold, italic) of the characters
  #                          in the plot
  #   font.family:           The font family (ie, sans) of the characters
  #                          in the plot
  #   plot.type              A string containing the plot type, which can be
  #                          "combined", "density" or "contour"
  #   grayscale               A boolean specifying the use of a graysacele  
  #   x.axis.title:          The x-axis title. 
  #   y.axis.title:          The y-axis title.
  #   x.range:               The range of values on the x-axis, 
  #                          as a vector (min, max)
  #   y.range:               The range of values on the y-axis,
  #                          as a vector (min, max)
  #   clipping:              A Boolean specifying whether contour labels
  #                          may be clipped from the plot area.
  #   extend.panel:          A Boolean specifying whether the grid panel 
  #                          may be extended. 
  #   WTP.thresholds:        A vector containing WTP thresholds to draw. 
  #   basecase:              A vector as: (incremental QALYs, incremental costs)                       
  #   average.PSA:           A vector as: (incremental QALYs, incremental costs)
  #
  # Returns:
  #   A gtable object with the plot data.  
  
  if (grayscale == TRUE){
    color.palette <- "Greys"
  } else {
    color.palette <- "Spectral"
  }
  
  if (missing(x.range)) {
    x.range <- c(min(df.density.normalised$X), max(df.density.normalised$X))
  }
  
  if (missing(y.range)) {
    y.range <- c(min(df.density.normalised$Y), max(df.density.normalised$Y))
  }
  
  
  if (!is.null(average.PSA) | !is.null(basecase) | !is.null(WTP.thresholds)) {
    density.barheight <- NULL
  } else {
    density.barheight <- 15
  }
  
  element.list <- list()
  label.list1  <- list("far.from.others.borders", 
                       "calc.boxes", 
                       "enlarge.box",
                       rot = 0,
                       hjust = 0, 
                       vjust = 0, 
                       box.color = NA,
                       fill = "transparent", 
                       "draw.rects")
  
  if (!is.null(WTP.thresholds)) {
    #init segment dataframe
    segment.df <- as.data.frame(matrix(0, ncol = 6, 
                                       nrow = length(WTP.thresholds)))
    min.x <- min(df.density.normalised$X)
    min.y <- min(df.density.normalised$Y)
    max.x <- max(df.density.normalised$X)
    max.y <- max(df.density.normalised$Y)
    colnames(segment.df) <- c("segment,start.x", "segment.end.x",
                              "min.y", "max.y", "WTP Threshold", "i")
    # calculate the coordinates of each WTP segment 
    for(i in 1:length(WTP.thresholds)) {
      segment.start.x <- min.y / WTP.thresholds[i]
      segment.end.x <- max.y / WTP.thresholds[i]
      
      if (segment.end.x > max.x){
        segment.end.x <- max.x 
        max.y <- WTP.thresholds[i] * max.x
      }
      if(segment.start.x < min.x) {
        segment.start.x <- min.x
        min.y <- WTP.thresholds[i] * min.x
      }
      
      # if a segment is entirely out of the plot window:
      if( segment.start.x > max.x){
        segment.start.x <- min.x
        segment.end.x <- min.x
        min.y <- min(df.density.normalised$Y)
        max.y <- min(df.density.normalised$Y)
      }
      segment.df[i,] <- c(segment.start.x, segment.end.x, min.y, max.y, 
                          WTP.thresholds[i], i)
    }  
    WTP <- factor(segment.df[, 6], labels = as.character(segment.df[, 5]))
  } 
  
  # Generate plot using ggplot() call
  
  if (plot.type == "contour") {
    #plot only contour
    plot.contour <- ggplot(data = df.density.normalised, 
                           aes(x = X, y = Y, z = Z)) + 
      theme_minimal()  +
      geom_contour(aes(z  = plot.data.contour$contour.levels),
                   breaks = rev(contour.levels),
                   size   = 0.5,
                   colour = "black") +
      theme(panel.grid.major = element_line(colour = "gray30", size = 0.25),
            panel.grid.minor = element_line(colour = "gray30", size = 0.25),
            panel.ontop      = TRUE,
            text             = element_text(size = font.size, 
                                            family = font.family, 
                                            face   = font.face),
            legend.spacing.y = unit(0.15, "cm")) +
      labs(x = x.axis.title, 
           y = y.axis.title) +
      guides(fill = guide_colourbar(barheight = density.barheight))
    
  } else if (plot.type == "density") {  
    #plot only density
    plot.contour <- ggplot(data = df.density.normalised, 
                           aes(x = X, y = Y, z = Z)) + 
      geom_tile(aes(fill = Z), alpha = 1) + 
      scale_fill_distiller(name      = legend.title,
                           palette   = color.palette, 
                           direction = -1,
                           guide     = "colourbar") + 
      theme_minimal()  +
      theme(panel.grid.major = element_line(colour = "gray30", size = 0.25),
            panel.grid.minor = element_line(colour = "gray30", size = 0.25),
            panel.ontop      = TRUE,
            text             = element_text(size = font.size, 
                                            family = font.family, 
                                            face   = font.face),
            legend.spacing.y = unit(0.15, "cm")) +
      labs(x = x.axis.title, 
           y = y.axis.title) +
      guides(fill = guide_colourbar(barheight = density.barheight))
    
  } else {
    #default combined plot, displaying both density and contours  
    plot.contour <- ggplot(data = df.density.normalised, 
                           aes(x = X, y = Y, z = Z)) + 
      geom_tile(aes(fill = Z), alpha = 1) + 
      scale_fill_distiller(name      = legend.title,
                           palette   = color.palette, 
                           direction = -1,
                           guide     = "colourbar") + 
      theme_minimal()  +
      geom_contour(aes(z  = plot.data.contour$contour.levels),
                   breaks = rev(contour.levels),
                   size   = 0.5,
                   colour = "black") +
      theme(panel.grid.major = element_line(colour = "gray30", size = 0.25),
            panel.grid.minor = element_line(colour = "gray30", size = 0.25),
            panel.ontop      = TRUE,
            text             = element_text(size = font.size, 
                                            family = font.family, 
                                            face   = font.face),
            legend.spacing.y = unit(0.15, "cm")) +
      labs(x = x.axis.title, 
           y = y.axis.title) +
      guides(fill = guide_colourbar(barheight = density.barheight))
  }
  
  if (!is.null(WTP.thresholds)) {
    element.list <- append(element.list, geom_segment(data = segment.df, 
                                                      aes(x    = segment.df[, 1],
                                                          xend = segment.df[, 2],
                                                          y    = segment.df[, 3],
                                                          yend = segment.df[, 4], 
                                                          linetype = WTP),
                                                      color = "black",
                                                      size  = 0.5,
                                                      inherit.aes = F)) 
  } 
  
  if (!is.null(average.PSA) & !is.null(basecase)) {
    point.df <- rbind.data.frame(average.PSA, basecase)
    point.df <- cbind(point.df, c("Average", "Base case"))
    colnames(point.df) <- c("x", "y", "Type")
    points.data <- factor(c(1, 2), labels = as.character(point.df$Type))
    
    element.list <- append(element.list, geom_point(inherit.aes = F,
                                                    data = point.df,
                                                    aes(x = x,
                                                        y = y,
                                                        group = "Type",
                                                        shape = points.data),
                                                    color = "black",
                                                    fill  = "red",
                                                    size  = 3)) 
    element.list <- append(element.list,
                           scale_shape_manual(name = "Points", 
                                              values = c(23, 21), 
                                              labels = c("Average",
                                                         "Base case")))
  } else if (!is.null(average.PSA) & is.null(basecase)) {
    point.df <- rbind.data.frame(average.PSA)
    point.df <- cbind(point.df, c("Average"))
    colnames(point.df) <- c("x", "y", "Type")
    points.data <- factor(c(1), labels = as.character(point.df$Type))
    
    element.list <- append(element.list, geom_point(inherit.aes = F,
                                                    data = point.df,
                                                    aes(x = x,
                                                        y = y,
                                                        group = "Type",
                                                        shape = points.data),
                                                    color = "black",
                                                    fill  = "red",
                                                    size  = 3)) 
    element.list <- append(element.list,
                           scale_shape_manual(name = "Points", 
                                              values = c(23), 
                                              labels = c("Average")))
  } else if (is.null(average.PSA) & !is.null(basecase)) {
    point.df <- rbind.data.frame(basecase)
    point.df <- cbind(point.df, c("Base case"))
    colnames(point.df) <- c("x", "y", "Type")
    points.data <- factor(c(1), labels = as.character(point.df$Type))
    
    element.list <- append(element.list, geom_point(inherit.aes = F,
                                                    data = point.df,
                                                    aes(x = x,
                                                        y = y,
                                                        group = "Type",
                                                        shape = points.data),
                                                    color = "black",
                                                    fill  = "red",
                                                    size  = 3)) 
    element.list <- append(element.list,
                           scale_shape_manual(name = "Points", 
                                              values = c(21), 
                                              labels = c("Base case")))
  }
  
  element.list <- append(element.list, geom_dl(aes(label = ..level.., 
                                                   x = plot.data.contour$x,
                                                   y = plot.data.contour$y, 
                                                   z = plot.data.contour$contour.levels, 
                                                   fontface = "bold"),
                                               inherit.aes = F,
                                               color = "gray15",
                                               cex = 0.75,
                                               method = label.list1,
                                               stat = "contour", 
                                               breaks = rev(contour.levels)))
  
  if (clipping == FALSE) {
    # limits rendering to given coordinates, 
    # does not exclude (clip) data during plot generation
    element.list <- append(element.list, coord_cartesian(xlim = x.range,
                                                         ylim = y.range,
                                                         expand = FALSE))
  } else {
    #clip datapoints to fall within a range
    element.list <- append(element.list, xlim(x.range))
    element.list <- append(element.list, ylim(y.range))                           
  }
  
  #add elements stored in list
  plot.contour <- plot.contour + element.list
  
  if (extend.panel == TRUE) {
    # extend the plot panel outside of the density area
    # so that contour labels are not partially cropped. 
    plot.contour <- ggplot_gtable(ggplot_build(plot.contour))
    plot.contour$layout$clip[plot.contour$layout$name == "panel"] <- "off"
  }
  return(plot.contour)
}

message("Copyright (C) 2020, Joost Geenen\n",
        "This program comes with ABSOLUTELY NO WARRANTY.\n",
        "This is free software, and you are welcome to redistribute it",
        " under certain conditions.\n","You should have received a copy of the GNU",
        " General Public License along with this program.\n",
        "If not, see <https://www.gnu.org/licenses/>.")
### ----------------------------------------------------------------------- ###
### -------------- Do not change the part in between / above -------------- ###



### ----------------------------------------------------------------------- ###
### ------------------ Setup your data and plot settings ------------------ ###

# Your raw data should have the following characteristics:

# - Saved as a .csv file
# - The first column should be incremental effects (QALYs / LYs / etc)
# - The second column should be incremental costs

# Then, prepare to load your data:

# - set the filepath of your raw datafile, do not forget the .csv
#   when your datafile is in the same folder as this script, 
#   only the filename is needed (e.g. "your_data.scv")
# - set whether your data has a header (ie, column names)
# - set the type of decimal seperator (, or . within " ")
# - set the column separator used for your .csv file (within " "). 
# You can see the seperator when opening your csv file with, for example, Excel. 
filepath   <- "eHTA case study data, used for figure 2.csv"
has.header <- FALSE
decimal.separator <- ","
column.separator  <- ";"

# We recommend 10.000 rows, although somewhere between 
# 1.000 and 100.000 will produce proper figures. 
# More than 10.000 will slow various computations whilst 
# it does not add information to your plots. 
# We therefore recommend a maximum number of rows of 10.000 
# If you would like to limit your rows to this number, 
# set the 'limit_rows' variable to TRUE and specify 
# the number of rows in number_rows: 
# If your data has less than 10.000 rows, leave this setting at FALSE. 
limit.rows <- FALSE
number.rows <- 10000

# The number of bins determines the granularity of the plot. 
# Warning: Larger bin numbers require more RAM. 
# 1000 bins typically produces images without pixelation. 
# More bins does not provide better images whilst it increases computation time.
# Less bins (eg, 100, 500), provide decent figure with limited RAM usage and 
# reduced computation time. 
# We therefore recommend using 1000 bins. 

# set bin.number here
bin.number <- 1000

# You can specify the following plot characteristics: 
# - specify desired contour levels. 
# - set as, for example: contour.levels <- c(0.95, 0.5, 0.1)
contour.levels <- c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)

# - Specify WTP thresholds to add to the plot. 
#   Example: WTP.thresholds <- c(30000, 80000)
#   If you don't want to add these, run "WTP.thresholds <- c()"
WTP.thresholds <- c(20000, 50000, 80000)

# - If you want to add the base-case results to your plot, 
#   set as: "basecase <- c(<basecase incr. QALYs>, <basecase incr. costs>)
#   Example: basecase <- c(0.5, 2000)
#   If you don't want to add this, set: "basecase <- c()"
basecase <- c(0.0978, 4946)

# - If you want to add a marker with the average of the PSA iterations, 
#   set add.average.PSA <- TRUE
add.average.PSA <- TRUE

# - Set Font Title for the Plot
x.axis.title <- "Incremental QALYs"
y.axis.title <- "Incremental Costs"
font.size    <- 14

# - Set font characteristics for the plot
font.face    <- "bold"         ## choose bold or plain or italix
font.family  <- "sans"         ## 
legend.title <- "Density"      ## Set your desired legend title

# - Set plot type. "combined" yields contours as well as density, 
# "density" provides only the density plot and "contour" provides
# only the contours. Default = "combined"
plot.type <- "combined"

# - Set color-gradient or a grayscale gradient (e.g. for printing)
#  It is not recommended to use grayscale, default = FALSE
grayscale <- FALSE

# Now, proceed running the following parts, 
# parts within a 'do not change the part in between / below (or above)
# should just be run but do not require input from the user. P

### -------------------- End of data and plot settings -------------------- ###
### ----------------------------------------------------------------------- ###


### -------------- Do not change the part in between / below -------------- ###
### ----------------------------------------------------------------------- ###

# Load data
data <- read.csv(file   = filepath, 
                 header = has.header, 
                 dec    = decimal.separator, 
                 sep    = column.separator)

print(paste("your data has", as.character(nrow(data)), "rows"))
if (limit.rows == TRUE & nrow(data) > number.rows){
  data <- data[1:number.rows, ]
}

# Perform kernel density estimation
kde.data <- kde2d(data[, 1], data[, 2], n = bin.number)
plot.data.contour     <- ProcessContourData(kde.data)
df.density.normalised <- GenerateNormalisedDensity(kde.data)

#calculate average PSA results
if (add.average.PSA == TRUE) {
  average.PSA <- c(mean(data[, 1]), mean(data[, 2]))
} else {
  average.PSA <- NULL
}

# Generate the ggplot plot object
contour.plot <- GeneratePlot(df.density.normalised, 
                             legend.title,
                             plot.data.contour,
                             contour.levels,
                             font.size,
                             font.face,
                             font.family,
                             plot.type,
                             grayscale,
                             x.axis.title,
                             y.axis.title,
                             clipping = FALSE,
                             extend.panel = TRUE,
                             WTP.thresholds = WTP.thresholds,
                             basecase = basecase,
                             average.PSA = average.PSA)

# display the plot
grid.newpage()
grid.draw(contour.plot)


### ----------------------------------------------------------------------- ###
### -------------- Do not change the part in between / above -------------- ###



### -------------------- Options for Saving your plot --------------------- ###
### ----------------------------------------------------------------------- ###

# If you wish to save this plot, enter the required filename and settings here

# set your filepath (ie, the folder)
filepath <- "C://your_path/saved_PSAReD_plots/"

# set the file name of the new plot. 
# WARNING: It will overwrite files / plots with the same name!
filename <- "Figurename.PNG"

# set plot dpi (resolution)
dpi <- 600

# Use size of the "Plots" Panel in Rstudio?
# "TRUE" will guarantee that your figure is saved exactly like you see it now. 
# Another user, with a differently sized "Plots" panel, 
# will then however get a different plot using the same data
# Selecting "FALSE" allows you to specifiy your own, fixed figure size.
# This figure will be different from the one you see in your "Plots" panel
# But you will then reproduce this exact figure using the size values values. 
use.my.panel.size <- FALSE

# If you have set use.my.panel.size to FALSE, set the plot size in inches:
plot.width  <- 8.47
plot.height <- 4.25

# Saving is now set-up, run the part below. 

### -------------- Do not change the part in between / below -------------- ###
### ----------------------------------------------------------------------- ###

if (use.my.panel.size == TRUE) {
  ggsave(contour.plot, path = filepath, filename = filename, dpi = dpi)
} else {
  ggsave(contour.plot, path = filepath, filename = filename, dpi = dpi,
         width = plot.width, height = plot.height, units = "in")
}
### ----------------------------------------------------------------------- ###
### -------------- Do not change the part in between / above -------------- ###



### -------------- Options for zooming in on a specific area -------------- ###
### ----------------------------------------------------------------------- ###

# If you wish, you can zoom in on a particular area of the plot. 
# set x- and y range as c(min, max)
x.range <- c(0, 0.5)
y.range <- c(1500, 8000)

# To prevent the contour labels from clipping of the sides, 
# the drawing panel is extended. 
# This however, has some consequences when zooming. 

# Select ONE the following and then go to "Now, run the part below":

# 1 If your contour labels all lie within the plot, set: 
clip <- FALSE

# If your contour labels (partly) fall outside the plot, you have 2 choices:
# 1: Clip them off, which may yield a pretier plot. set:
clip <- FALSE

# 2: Extend panel and clip data to fit the panel, this may yield a pretier plot. 
# set:
clip <- TRUE

# Now, run the part below. 

### -------------- Do not change the part in between / below -------------- ###
### ----------------------------------------------------------------------- ###
if (clip == FALSE) {
  extend.panel <- FALSE
  clipping     <- FALSE
} else{
  extend.panel <- TRUE
  clipping     <- TRUE
}

zoomed.plot <- GeneratePlot(df.density.normalised, 
                            legend.title,
                            plot.data.contour,
                            contour.levels,
                            font.size,
                            font.face,
                            font.family,
                            plot.type,
                            grayscale,
                            x.axis.title,
                            y.axis.title,
                            x.range,
                            y.range,
                            clipping,
                            extend.panel,
                            WTP.thresholds,
                            basecase,
                            average.PSA)
grid.newpage()
grid.draw(zoomed.plot)
### ----------------------------------------------------------------------- ###
### ---------- Do not change the part in between / above ------------------ ###



### -------------- Options for zooming in on a specific area -------------- ###
### ----------------------------------------------------------------------- ###

# If you wish to save this plot, enter the required filename and settings here

# set your filepath (ie, the folder)
filepath  <- "C://your_path/saved_PSAReD_plots/"

# set the file name of the new plot. 
# WARNING: It will overwrite files / plots with the same name!
filename <- "Figurename.PNG"

# set plot dpi (resolution)
dpi <- 600

# Use size of the "Plots" Panel in Rstudio?
# "TRUE" will guarantee that your figure is saved exactly like you see it now. 
# Another user, with a differently sized "Plots" panel, will then however 
# get a different plot using the same data
# Selecting "FALSE" allows you to specifiy your own, fixed figure size.
# This figure will be different from the one you see in your "Plots" panel
# But you will then reproduce this exact figure using the size values values.
use.my.panel.size <- FALSE

# If you have set use.my.panel.size to FALSE, set the plot size in inches:
plot.width  <- 8.47
plot.height <- 4.25

# Saving is now set-up, run the part below. 

### -------------- Do not change the part in between / below -------------- ###
### ----------------------------------------------------------------------- ###

if (use.my.panel.size == TRUE) {
  ggsave(zoomed.plot, path = filepath, filename = filename, dpi = dpi)
} else {
  ggsave(zoomed.plot, path = filepath, filename = filename, dpi = dpi,
         width = plot.width, height = plot.height, units = "in")
}
### ----------------------------------------------------------------------- ###
### -------------- Do not change the part in between / above -------------- ###
