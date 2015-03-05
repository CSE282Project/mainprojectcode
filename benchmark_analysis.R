#Set working directory
setwd("/Users/sjroth/Documents/CSE282/Final/mainprojectcode/")

#Load data
data <- read.csv('benchmark.csv',header = F)
names(data) <- c('Algorithm','Length','Motif','Recursive.Count',
                 'Running.Time','Weight','Elapsed.Time')

#Label the algorithms
data$Algorithm[data$Algorithm == 0] <- 'Greedy'
data$Algorithm[data$Algorithm == 1] <- 'Dynamic Programming'
data$Algorithm[data$Algorithm == 2] <- 'Greedy 2'

#Begin plotting
#Pick the data for plotting.
data.plot <- data[data$Length == 600,c(1,3:7)]

library('ggplot2')
p1 <-
  ggplot(data, aes(x=Motif, y=Recursive.Count, colour=Algorithm)) +
  geom_line() +
  ggtitle("Motif Length vs. Recursive Count") +
  ylab("Recursive Count")

p2 <-
  ggplot(data, aes(x=Motif, y=Running.Time, colour=Algorithm)) +
  geom_line() +
  ggtitle("Motif Length vs. Running Time") +
  ylab("Running Time")

p3 <-
  ggplot(data, aes(x=Motif, y=Weight, colour=Algorithm)) +
  geom_line() +
  ggtitle("Motif Length vs. Weight") 

p4 <-
  ggplot(data, aes(x=Motif, y=Elapsed.Time, colour=Algorithm)) +
  geom_line() +
  ggtitle("Motif Length vs. Elapsed Time") +
  ylab("Elapsed Time")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p1,p2,p3,p4,cols=2)

