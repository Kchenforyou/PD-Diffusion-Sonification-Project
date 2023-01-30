#NOTE: Loop from CSV file 91-end!!!!

library("ggplot2")
library("ggforce")
library("gifski")

# it's hardcoded for now... I know, not a good practice...
setwd("C:/Users/kchen/Documents/PD_diffusion_project/CSV_files6/")

# need to figure out how to transfer coordinates from SS to here.
#plot(-400:400,-400:400,type="n",xlab="",ylab="",main="Test draw.circle")
# a <- -500:500
# b <- -500:500
# p <- data.frame(cbind(a,b))

# object to access csv files of SkewerProcess array.
files <- Sys.glob(file.path(getwd(), 'SkewerProcess*csv'))

# change directory to save PNG files
setwd("C:/Users/kchen/Documents/PD_diffusion_project/CSV_files6/Frames")

# save current path to create gifs below
# png_path = getwd()

# create ggplot object
gg <-  ggplot() + 
  geom_blank() + theme_bw() +
  theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
  theme(axis.line = element_line(color = 'black')) +
  theme(legend.position = "none")

# this is to initialize the for-loop
num_frame = 500 # initially 500
# this is to help create labels like "001"
digits = floor(log10(num_frame)) + 1

for (i in 1:num_frame){
  # store csv file
  sp_test <- read.csv(files[i])
  
  # radius and coordinates for circle
  r <- sqrt(sp_test$x1)
  xc <- sp_test$x5
  yc <- sp_test$x6
  
  # for colors
  red = sp_test$x2
  green = sp_test$x3
  blue = sp_test$x4
  # opacity = sp_test$x5  no longer needed
  
  # create circles using ggplot
  gg +
  geom_rect(aes(xmin = -400, xmax = 400,
                ymin = -400, ymax = 400),
            fill = "white") +
  geom_circle(aes(x0 = x5, y0 = x6, 
                   r = sqrt(x1), 
                  fill = I(rgb(red = x2, green = x3, # I() solved color-wobble
                             blue = x4))),
              alpha = 0.4,
              data = sp_test) +
  coord_fixed()
  
  # write a short function that converts a number
  # into a character string with 0's. Hint, use logarithm?
  
  # create labels by attaching 0's
  zeros = strrep("0", digits - floor(log10(i)))
  
  # save as png
  pic <- paste0("SkewerProcess", zeros, i, ".png")
  ggsave(pic, width = 7, height = 7, dpi = 700)
}

# If you want high quality GIFs, use this code block.
# A warning: This will take a long time to run.
# Might want to experiment with 500 by 500?
png_files <- Sys.glob(file.path(getwd(), 'Skewer*.png'))
gifski(png_files[1:10], gif_file = "PD_5.gif",
       width = 600, height = 600,
       delay = 0.0625)
#3:02:23

# Noah's sample code:

# in r commands like lines, polygon, circle, 
# use 'col = rgb('red', 'green', 'blue', 'opacity')'
# note each arg above is just a number between 0-1 (continuous).
# sample code:
# polygon(X[i,1] + c(-S*a, rev(S)*a), 
# y = X[i,2] + (c(1,J,rev(J),1)-1)*yRatio, 
# col=rgb(Clrs[i,1],Clrs[i,2],Clrs[i,3], opacity), border=NA)
# S is seq. of herd populations for spindle.