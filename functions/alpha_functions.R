require(alphahull)

# alpha shapes
ashape_to_SPLDF <- function(x, proj4string=NA)
{
  if(class(x) != 'ashape')
    stop('this function only works with `ashape` class objects')
  
  # convert ashape edges to DF
  x.as.df <- as.data.frame(x$edges)
  
  # convert each edge to a line segment
  l.list <- list()
  for(i in 1:nrow(x.as.df))
  {
    # extract line start and end points as 1x2 matrices
    p1 <- cbind(x.as.df$x1[i], x.as.df$y1[i])
    p2 <- cbind(x.as.df$x2[i], x.as.df$y2[i])
    # row-bind into 2x3 matrix
    l.list[[i]] <- Line(rbind(p1, p2))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- Lines(l.list, ID=1)
  
  # copy over CRS data from original point data
  l.spl <- SpatialLines(list(l), proj4string=CRS(as.character(NA)))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <- SpatialLinesDataFrame(l.spl, data=data.frame(id=1), match.ID=FALSE)
  
  return(l.spldf)
}



# alpha hulls
ahull_to_SPLDF <- function(x, proj4string=NA)
{
  if(class(x) != 'ahull')
    stop('this function only works with `ahull` class objects')
  
  # convert ashape edges to DF
  x.ah.df <- as.data.frame(x$arcs)
  
  # convert each arc to a line segment
  l.list <- list()
  for(i in 1:nrow(x.ah.df))
  {
    # extract row i
    row_i <- x.ah.df[i,]
    
    # extract elements for arc()
    v <- c(row_i$v.x, row_i$v.y)
    theta <- row_i$theta
    r <- row_i$r
    cc <- c(row_i$c1, row_i$c2)
    # from arc()
    angles <- anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- cc[1] + r * cos(seqang)
    y <- cc[2] + r * sin(seqang)
    
    # convert to line segment
    l.list[[i]] <- Line(cbind(x,y))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- Lines(l.list, ID=1)
  
  # copy over CRS data from original point data
  l.spl <- SpatialLines(list(l), proj4string=CRS(as.character(NA)))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <- SpatialLinesDataFrame(l.spl, data=data.frame(id=1), match.ID=FALSE)
  
  return(l.spldf)
}
