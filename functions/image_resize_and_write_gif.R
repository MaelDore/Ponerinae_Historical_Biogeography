
### Original functions from magick by Jeroen Ooms (jeroen@berkeley.edu)
# magick::image_write_gif
# magick::write_png_files

# Changes: 
  # Allows to resize the output GIF to any size
  # Print progress

# Workflow:
#  1/ Convert to PNG
#  2/ Resize all images
#  3/ Export temporary PNG
#  4/ Read PNG
#  5/ Export as GIF

image_resize_and_write_gif <- function (image, path, width = NULL, height = NULL, delay = 1/10, ...) 
{
  # Convert and export temporary PNG
  png_paths <- resize_and_write_png_files(image = image, width = width, height = height)
  # Import first image to get dimensions
  png_tmp <- magick::image_read(path = png_paths[1])
  
  # Extract info of dimensions
  info <- magick::image_info(png_tmp)
  width <- info$width[1]
  height <- info$height[1]
  
  # Remove temporary directory
  on.exit(unlink(png_paths))
  
  message("Exporting final GIF..")
  
  # Export in GIF
  gifski::gifski(png_files = png_paths, gif_file = path, width = width, 
                 height = height, delay = delay, ...)
}

resize_and_write_png_files <- function(image, width, height)
{
  # Create temporary directory
  imgdir <- tempfile('tmppng')
  dir.create(imgdir)
  
  # Convert to PNG
  image <- magick::image_convert(image, format = 'png')
  # Extract info
  info <- magick::image_info(image)
  
  # Case with only height or width provided => Call an error
  if (is.null(width) != is.null(height))
  {
    stop("Width and Height should both be provided or ignored altogether")
  }
  
  # Case with height and width not provided => Use the dimension of the first frame
  if (is.null(width) & is.null(height))
  {
    width <- info$width[1]
    height <- info$height[1]
    
    if(length(unique(info$width)) > 1 || length(unique(info$height)) > 1)
    {
      message(sprintf("Images are not the same size. Resizing to dimension of the first frame: %dx%d...", width, height))
      image <- image_resize(image, sprintf('%dx%d!', width, height))
    } else {
      message(sprintf("Using intial dimensions of images. GIF is %dx%d.", width, height))
    }
  }
  
  # Case with height and width provided => Resize all frames
  if (!is.null(width) & !is.null(height))
  {
    initial_widths <- info$width
    initial_heights <- info$height
    
    any(initial_widths != width)
    any(initial_heights != height)
    
    if(any(initial_widths != width) || any(initial_heights != height))
    {
      message(sprintf("Images are resized to %dx%d...", width, height))
      image <- image_resize(image, sprintf('%dx%d!', width, height))
    } else {
      message(sprintf("Dimensions provided are the intial dimensions of images. GIF is %dx%d.", width, height))
    }
  }

  message("Exporting temporary PNG..")
  
  # Export temporary PNG
  vapply(X = seq_along(image), 
         FUN = function(i)
         {
            tmppng <- file.path(imgdir, sprintf("tmpimg_%05d.png", i))
            magick::image_write(image[i], path = tmppng, format = 'png')
         },
         FUN.VALUE =  character(1))
}



magick::image_write_gif 

function (image, path = NULL, delay = 1/10, ...) 
{
  png_files <- write_png_files(image)
  on.exit(unlink(png_files))
  info <- image_info(image)
  width <- info$width[1]
  height <- info$height[1]
  gifski::gifski(png_files, gif_file = path, width = width, 
                 height = height, delay = delay, ...)
}


write_png_files <- function(image){
  imgdir <- tempfile('tmppng')
  dir.create(imgdir)
  image <- image_convert(image, format = 'png')
  info <- image_info(image)
  width <- info$width[1]
  height <- info$height[1]
  if(length(unique(info$width)) > 1 || length(unique(info$height)) > 1){
    message(sprintf("Images are not the same size. Resizing to %dx%d...", width, height))
    image <- image_resize(image, sprintf('%dx%d!', width, height))
  }
  vapply(seq_along(image), function(i){
    tmppng <- file.path(imgdir, sprintf("tmpimg_%05d.png", i))
    image_write(image[i], path = tmppng, format = 'png')
  }, character(1))
}