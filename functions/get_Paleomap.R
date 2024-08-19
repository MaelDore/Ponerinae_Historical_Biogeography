
### paleoMap::getmap by Sonja Rothkugel <rothkugelsonja@aol.de>

###########################get_Paleomap#############################
#' get_Paleomap
#' 
#' Downloads a map of a specific age or a list of maps of specific ages from a model specified by the user.
#' Available models and ages can be found at 
#' https://github.com/GPlates/gplates_web_service_doc/wiki/Reconstruction-Models.
#' 
#' @usage get_Paleomap(ma, model = "SETON2012", show.plates = FALSE, 
#'                   save.as = NULL, colland = "#66666660", 
#'                   colsea = "#00509010", 
#'                   do.plot = TRUE, ...)
#' 
#' @param ma numeric. Age in ma(million years ago). Can also be a vector of ages (vector of numeric age values).
#' @param model character. The model the map should be created with. "SETON2012" (default), 
#' "MULLER2016", "GOLONKA", "PALEOMAP" or "MATTHEWS2016".
#' @param get.plate.polygons boolean. If the user wants to get the continental plate polygons or not. By default get.plate.polygons = FALSE.
#' @param get.plate.boundaries boolean. If the user wants to get the continental plate boundaries as lines, or not. By default get.plate.boundaries = FALSE.
#' @param show.plates boolean. If the user wants to display the continental plate boundaries (lines) on the plot, or not. By default show.plates = FALSE.
#' @param save.as character. The format the plots should be saved. "tiff", "pdf", "jpeg" or "png". 
#' By default save.as = NULL, plots are only shown and are not automatically saved as a file.
#' @param colland character. The color of the land masses. By default colland = "#66666660".
#' @param colsea character. The color of the sea. By default colsea = "#00509010".
#' @param do.plot logical. If a plot of the map is created or not. By default do.plot = TRUE. 
#' @param verbose logical. To display progress in loading maps
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to plot, such as main = "my own title" or main.col = "red".
#' @export
#' @examples
#' \dontrun{
#' 
#' library(mapast)
#' 
#' #with continental plates
#' map <- get_Paleomap(ma = 100, model = "SETON2012", show.plates = T)
#' coastlines <- map[[1]][[1]]
#' plates <- map[[2]][[1]]
#' 
#' #without continental plates
#' coastline <- get_Paleomap(ma = 100, model = "SETON2012")
#' 
#' #save map directly as pdf
#' get_Paleomap(ma = 100, model = "SETON2012", save.as="pdf")
#' 
#' #save multiple maps in one pdf
#' pdf("get_Paleomap_multi.pdf")
#' par(mfrow = c(2, 1))
#' get_Paleomap(ma = c(1, 2))
#' dev.off()
#' par(mfrow = c(1, 1))
#' 
#'}

get_Paleomap <- function(ma,                             # Set of ages
                         model = "SETON2012",            # Plate model to use
                         get.plate.polygons = FALSE,     # Obtain plate polygons
                         get.plate.boundaries = FALSE,   # Obtain plate boundaries as lines
                         show.plates = FALSE,            # Show plate boundaries (lines) on plot
                         save.as = NULL,                 # The format the plots should be saved. "tiff", "pdf", "jpeg" or "png".
                         colland = "#66666660",          # Color for lands on plot
                         colsea = "#00509010",           # Color for seas on plot
                         do.plot = TRUE,
                         verbose = TRUE, ...)            # Print progress every map
{
  #create empty list of shapes and plates
  shapes <- base::list()
  plates <- base::list()
  plate_polygons_all_ages <- base::list()
  
  #go through the ages given by the user
  for(ages in 1:base::length(ma))
  {
    #url for api request
    url <- base::paste0("https://gws.gplates.org/reconstruct/coastlines/?&time=", ma[ages], "&model=", model)
    #try to get the requested map of the model
    #print error message if map is not available
    #err: boolean if there was an error getting the map
    err <- FALSE
    shape <- tryCatch(
      {
        # Old version with rgdal
        # rgdal::readOGR(url, verbose = FALSE)
        
        # New version with sf objects and geojsonsf
        geojsonsf::geojson_sf(geojson = url)

      }, error = function(e) {
        err <- TRUE
        message(base::paste0("There is no map for ", ma[ages], " mya in ", model, " model. Please check the spelling, the age and the model you chose."))
        stop()
      }
    )
    
    # Convert to sf to SpatialPolygonsDataFrame
    shape <- as(object = shape, Class = "Spatial")
    
    #add metadata to the shape file
    #age and model
    shape@data$age <- ma[ages]
    shape@data$model <- model
    shape@data$name <- paste0("Coastlines_", ma[ages], "My_", model)
    
    #save shape in list of shapes
    shapes[[ages]] <- shape
    
    #errplate: boolean if there was an error getting the plates
    errplate <- FALSE
    err_plate_polygons <- FALSE
    
    # If user wants to get the plate boundaries (lines)
    if(get.plate.boundaries)
    {
      #check if model is golonka or paleomap: they dont have plates. Throw warning.
      if(model == "GOLONKA" || model == "PALEOMAP"){
        base::warning(base::paste0("No plate boundaries available for model ", model, "."))
        errplate <- TRUE
      }else{
        #create url for api request and try to get plates
        plateurl <- base::paste0("http://gws.gplates.org/topology/plate_boundaries/?time=", ma[ages], "&model=", model)
        platebounds <- tryCatch(
          {
            # Old version with rgdal
            # rgdal::readOGR(plateurl, verbose = FALSE)
            
            # New version with sf objects and geojsonsf
            geojsonsf::geojson_sf(geojson = plateurl)
            
          }, error = function(e){
            errplate <- TRUE
            #if cannot be loaded throw a warning that there was a problem loading the map
            base::warning(base::paste0("No Plate Boundaries available for ", ma[ages], " mya in ", model, " model. Please check the spelling, the age and the model you chose."))
            stop()
          }
        )
        
        # Convert to sf to SpatialPolygonsDataFrame
        platebounds <- as(object = platebounds, Class = "Spatial")
        
        #add metadata to plates, the age and the model
        platebounds@data$age <- ma[ages]
        platebounds@data$model <- model
        platebounds@data$name <- paste0("Plate_boundaries_", ma[ages], "My_", model)

        #save plate in list of plates
        plates[[ages]] <- platebounds
      }
    }
    
    # If user wants to get the plate polygons
    if(get.plate.polygons)
    {
      # Check if model is golonka or paleomap: they dont have plates. Throw warning.
      if (model == "GOLONKA" || model == "PALEOMAP")
      {
        base::warning(base::paste0("No plate boundaries available for model ", model, "."))
        err_plate_polygons <- TRUE
      } else {
        
        # Create url for api request and try to get plates
        plate_polygons_url <- base::paste0("https://gws.gplates.org/topology/plate_polygons/?&time=", ma[ages], "&model=", model)
        
        plate_polygons <- tryCatch(
          {
            # Old version with rgdal
            # rgdal::readOGR(plate_polygons_url, verbose = FALSE)
            
            # New version with sf objects and geojsonsf
            geojsonsf::geojson_sf(geojson = plate_polygons_url)
            
          }, error = function(e){
            err_plate_polygons <- TRUE
            #if cannot be loaded throw a warning that there was a problem loading the map
            base::warning(base::paste0("No Plate Polygons available for ", ma[ages], " mya in ", model, " model. Please check the spelling, the age and the model you chose."))
            stop()
          }
        )
        
        # Convert to sf to SpatialPolygonsDataFrame
        plate_polygons <- as(object = plate_polygons, Class = "Spatial")
        
        # Quick plot to check results
        # plot(plate_polygons)
        
        # Add metadata to plates, the age and the model
        plate_polygons@data$age <- ma[ages]
        plate_polygons@data$model <- model
        plate_polygons@data$name <- paste0("Plate_polygons_", ma[ages], "My_", model)
        
        # Save plate in list of plate polygons
        plate_polygons_all_ages[[ages]] <- plate_polygons
      }
    }
    
    #if there was no error getting the map, create parameter for plotting and plot the map
    if(!err & do.plot)
    {
      #getting final parameter list for plot
      #default parameter list for plotting
      graphparams.def <- base::list(x = shape, col = "white", border = FALSE
                                    , xlim = c(-180, 180), ylim = c(-90, 90)
                                    , xaxs = "i", yaxs = "i")
      #list of user defined graphical parameter
      graphparams.user <- base::list(...)
      names_graphparams.user <- base::as.vector(base::names(graphparams.user))
      names_graphparams.def <- base::as.vector(base::names(graphparams.def))
      #remove default parameter if user specifies a different value
      for( param in names_graphparams.user){
        if(param %in% names_graphparams.def) graphparams.def <- graphparams.def[ - base::which(base::names(graphparams.def) == param)] 
      }
      #complete new list of plotting parameters, including default and user specified ones
      graphparams <- c(graphparams.def, graphparams.user)
      
      # if user does not set plot=FALSE plot the shape file
      if (do.plot)
      {
        #if user specifies save.as, save the plot
        if(!base::is.null(save.as)){
          if(save.as == "tiff"){
            grDevices::tiff(base::paste0("get_Paleomap-", ma[ages], "mya_", model, ".tiff"), 
                            height = 10.5, width = 17, units = "cm", res = 300)
          }
          if(save.as == "jpeg"){
            grDevices::jpeg(base::paste0("get_Paleomap-", ma[ages], "mya_", model, ".jpeg"), 
                            height = 10.5, width = 17, units = "cm", res = 300)
          }
          if(save.as == "png"){
            grDevices::png(base::paste0("get_Paleomap-", ma[ages], "mya_", model, ".png"), 
                           height = 10.5, width = 17, units = "cm", res = 300)
          }
        }
        #define the size of the margin of the plot and save the former definition
        def.mar <- graphics::par("mar")
        #define the margin size
        graphics::par(mar = c(2, 2, 2, 2))
        #do a first plot with the graphical parameters set by the user
        base::do.call(sp::plot, graphparams)
        #if the user set show.plates true and there was no problem loading the plates add the plates to the plot
        if (!errplate && show.plates && get.plate.boundaries)
        {
          sp::plot(platebounds, add = T, col = "#66666680")
        }
        #draw a rectangle showing the sea
        graphics::rect(xleft = -180, xright = 180, ybottom = -90, 
                       ytop = 90, col = colsea, 
                       border = FALSE)
        #add x-axis and x-axis label
        graphics::axis(side = 1, pos = -84, lwd = 0, xaxp = c(180, -180, 4), col.ticks = "darkgrey", 
                       col.axis = "darkgrey", cex.axis = 0.6)
        graphics::axis(side = 1, pos = -89, lwd = 0, at = 0 , labels = "Longitude", col.ticks = "darkgrey", 
                       col.axis = "darkgrey", cex.axis = 0.6)
        #add y-axis and y-axis label
        graphics::axis(side = 2, pos = -175, lwd = 0, yaxp = c(90, -90, 4), col.ticks = "darkgrey", 
                       col.axis = "darkgrey", cex.axis = 0.6, las = 1)
        graphics::axis(side = 2, pos = -178, lwd = 0, at = 0 , labels = "Latitude", col.ticks = "darkgrey", 
                       col.axis = "darkgrey", cex.axis = 0.6)
        #add model and age info top right of the plot
        graphics::axis(side = 3, pos = 89, lwd = 0, at = 135 , labels = model, col.ticks = "darkgrey", 
                       col.axis = "darkgrey", cex.axis = 1)
        graphics::axis(side = 3, pos = 81, lwd = 0, at = 135 , labels = base::paste(ma[ages], " mya", sep = ""), 
                       col.ticks = "darkgrey", col.axis = "darkgrey", cex.axis = 0.7)
        #add the landmasses to the plot
        sp::plot(shape, col = colland, border = FALSE, add = TRUE)
        #restore the former graphical mar parameters
        graphics::par(mar = def.mar)
        #if the user wants to have the maps as pdf
        if(!base::is.null(save.as) && save.as == "pdf"){
          filename <- base::paste0("get_Paleomap-", ma[ages], "mya_", model, ".pdf")
          grDevices::pdf(filename, width= 6.885417, height = 4.291667)
          #define the size of the margin of the plot and save the former definition
          graphics::par(mar = c(2, 2, 2, 2))
          #do a first plot with the graphical parameters set by the user
          plotmap <- base::do.call(sp::plot, graphparams)
          if (!errplate && show.plates && get.plate.boundaries)
          {
            plotmap <- sp::plot(platebounds, add = T, col = "#66666680")
          }
          #draw a rectangle showing the sea
          plotmap <- graphics::rect(xleft = -180, xright = 180, ybottom = -90, 
                                    ytop = 90, col = colsea, 
                                    border = FALSE)
          #add x-axis and x-axis label
          plotmap <- graphics::axis(side = 1, pos = -84, lwd = 0, xaxp = c(180, -180, 4), col.ticks = "darkgrey", 
                                    col.axis = "darkgrey", cex.axis = 0.6)
          plotmap <- graphics::axis(side = 1, pos = -89, lwd = 0, at = 0 , labels = "Longitude", col.ticks = "darkgrey", 
                                    col.axis = "darkgrey", cex.axis = 0.6)
          #add y-axis and y-axis label
          plotmap <- graphics::axis(side = 2, pos = -175, lwd = 0, yaxp = c(90, -90, 4), col.ticks = "darkgrey", 
                                    col.axis = "darkgrey", cex.axis = 0.6, las = 1)
          plotmap <- graphics::axis(side = 2, pos = -178, lwd = 0, at = 0 , labels = "Latitude", col.ticks = "darkgrey", 
                                    col.axis = "darkgrey", cex.axis = 0.6)
          #add model and age info top right of the plot
          plotmap <- graphics::axis(side = 3, pos = 89, lwd = 0, at = 135 , labels = model, col.ticks = "darkgrey", 
                                    col.axis = "darkgrey", cex.axis = 1)
          plotmap <- graphics::axis(side = 3, pos = 81, lwd = 0, at = 135 , labels = paste(ma[ages], " mya", sep = ""), 
                                    col.ticks = "darkgrey", col.axis = "darkgrey", cex.axis = 0.7)
          #add the landmasses to the plot
          plotmap <- sp::plot(shape, col = colland, border = FALSE, add = TRUE)
          base::print(plotmap)
          #restore the former graphical mar parameters
          graphics::par(mar = def.mar)
        }
        #if the user wanted to save the plot close the device after plotting
        if(!base::is.null(save.as)){
          grDevices::dev.off()
        }
      }
    }
    
    ## Print progress
    if ((ages %% 1 == 0) & verbose)
    {
      cat(paste0(Sys.time(), " - Maps loaded for model ",model," at age ", ma[ages]," My - n°", ages, "/", length(ma),"\n"))
    }
  }
  
  
  ## Make a list of coastlines

  # Extract names to identify elements of the shapes list
  coastlines_names <- unlist(lapply(X = shapes, FUN = function (x) { x@data$name[1] } ))
  names(shapes) <- coastlines_names
  
  output_list <- base::list(shapes)
  names(output_list) <- c("Coastlines")
    
  ## Add plate polygons if needed
  if(!err_plate_polygons && get.plate.polygons)
  {
    
    # Extract names to identify elements of the plates list
    plates_polygons_names <- unlist(lapply(X = plate_polygons_all_ages, FUN = function (x) { x@data$name[1] } ))
    names(plate_polygons_all_ages) <- plates_polygons_names
    
    # Append plate boundaries to the output list
    old_names <- names(output_list)
    output_list <- append(x = output_list, values = list(plate_polygons_all_ages))
    new_names <- c(old_names, "Plate_polygons")
    names(output_list) <- new_names
  }
  
  ## Add plate boundaries (lines) if needed
  if(!errplate && get.plate.boundaries)
  {
    # return(base::list(shapes, plates))
    
    # Extract names to identify elements of the plates list
    plates_names <- unlist(lapply(X = plates, FUN = function (x) { x@data$name[1] } ))
    names(plates) <- plates_names
    
    # Append plate boundaries to the output list
    old_names <- names(output_list)
    output_list <- append(x = output_list, values = list(plates))
    new_names <- c(old_names, "Plate_boundaries")
    names(output_list) <- new_names
  }
  
  ## Return the shape file/ list of shapes
  return(output_list)
}