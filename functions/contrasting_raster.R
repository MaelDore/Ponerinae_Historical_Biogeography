### Function to contrast a raster between a min and max value ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com


# Input = Raster Stack/Layer

# Output = Raster Stack/Layer with values framed between zmin and zmax



### Function to contrast a raster between a min and max value
contrasting_raster <- function(x, zmin, zmax)
{
  # Create raster with null values for terrestrial areas
  continental_mask <- (calc(x, fun = sum) >= 0) - 1
  # plot(continental_mask)
  
  # Create raster with NA values outside of the species range (assuming a sum of zero is outside the range)
  community_mask_PA <- (calc(x, fun = sum) > 0) - 1
  community_mask <- community_mask_PA
  community_mask[community_mask[] < 0] <- NA
  # plot(community_mask)
  
  x[x[] <= zmin] <- zmin  # Fix low values
  x[x[] >= zmax] <- zmax  # Fix high values
  
  x <- mask(x, mask = community_mask)  # Cut out values that are outside species range
  
  y <- continental_mask + zmin  # Create final new raster from continental mask with baseline = zmin
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}
