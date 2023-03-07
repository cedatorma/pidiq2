##pidiq_calc2() - Alternate version of PIDIQ (with image segmentation for
##extracting plants and eliminating background noise which can be applied to
##multiple-well/plate screen images).

#Requires:
#The mcgrittr "%>%" pipe...
#magick (https://github.com/ropensci/magick)
#EBImage (https://bioconductor.org/packages/release/bioc/html/EBImage.html) packages.



##
#Function Parameters:
##

  #img_file - path to plant/plate spray image.
  #output_dir - path to output directory for storing processed images and saving
            #PIDIQ stats results tables (default: './test_pidiq2').

##Pixel H,S,V cutoffs and image segmentation parameters:

#H,S,V cutoffs for healthy/diseased pixels,
#provided in the pidiq_presets.csv file:
  #filter - name of filtering preset corresponding to the 'name' column of the presets file.
  #assumed to be found in the current working directory where the script is called...


##OLD:
#Note: For development purposes, spectral filtering defaults are currently
#hard-coded based on A. thaliana HTP plate images
#WIP: incorporating user-supplied H,S,V values...

  #sat_min - minimum pixel saturation cutoff for spectral filtering
            #(H, and V will be added in later versions) (default: 0.02).
######


  #gsize_cut - minimum segmentation group size threshold for distinguishing plant from
            #spurious background pixels which pass filtering, if not supplied, then
            #a segmentation group size distribution plot will be generated to help you
            #to interactively select the value (default: NULL).
  #test_plt - plot preliminary segmentation group images (bw and watershed) (default: FALSE).


##For multiple-well seedling spray images (single plant image assumed by default),
##User-Supplied plate-coordinates for seedling wells:
  #col_lab - indices for well columns (default: NULL).
  #row_lab - indices for well rows (default: NULL).


#Print progress-messages to console:
  #msg = FALSE


##
#Description:
##
#PIDIQv2 aims to overcome some limitations of the previous PIDIQ macro and
#accounts for variability in well/seedling/ecotype growth on a given plate/spray,
#increasing the reliability and accuracy of resulting PIDIQ summary stats.

#Some of these limitations in the previous version of PIDIQ  are:
#1) Inclusion of non-plant/background pixels (perlite, well/plate plastic
#boundaries) after spectral filtering, which in some instances results in
#spurious PIDIQ scores, especially for wells with absent/small plants.

#2) Overgrowth of plants from neighbouring crops/wells and off-center growth of
#plants within a well, causing over- and under- estimation of plant sizes,
#respectively, and potential artifacts in reliably calculating diseased plant tissue.


#The approach relies on spectral filtering and image segmentation algorithms to
#distinguish plants from plate background and also individual plants from one
#another.

#First, spectral green/yellow filtering is applied to extract foreground pixels,
#which are next assigned to segmentation groups based on:

#1) neighbouring foreground pixels (bw), and 2) watershed segmentation around
#seed foreground-pixels determined using distance from nearest background
#pixels.

#- Any background pixels which pass can be filtered out if they belong to
#segmentation groups with < the supplied parameterL <gsize_cut> pixels.

#bw segmentation is first used to assign segmented pixel groups to defined crop/well
#regions in the image (proportion of pixels in crop > 75% or maximum # pixels, and other criteria...).

#If any bw groups are merged across multiple wells, then watershed groups are used to
#separate different crops and make final pixel-crop assignments.

#To do this, watershed groups found within a merged bw group are searched,
#and distinguished into those which are uniquely assigned to a crop, or
#overlap multiple crops. For the latter, these can either be reassigned to a crop
#based on 3 nearest-neighbouring watershed groups (if they are all uniquely assigned to
#the same crop), otherwise they are omitted remove them. With the final segmentation group
#assignments made, their corresponding pixel-crop assignments are updated accordingly.

#Crops/wells without any assignable segmentation groups are excluded from PIDIQ
#analysis, and indicated as NA's in the final summary stats table.


#' Calculate PIDIQ statistics from a plant-immunity elicitation experiment image.
#'
#' @param img_file Path to plant/plate spray image.
#' @param output_dir Path to output directory for storing processed images, if test_plt == TRUE (default: './test_pidiq2').
#' @param gsize_cut Minimum segmentation group size threshold for distinguishing plant from
#' spurious background pixels which pass filtering, if not supplied, then
#' a segmentation group size distribution plot will be generated to help you
#' to interactively select the value (default: NULL).
#' @param test_plt Plot preliminary segmentation group images (bw and watershed) (default: FALSE).
#' @param col_lab For multi-plant images: Supplied labels corresponding to plate/flat column indices (default: NULL).
#' @param row_lab For multi-plant images: Supplied labels corresponding to plate/flat row indices (default: NULL).
#' @param filter Name of filtering preset corresponding to the 'name' column of the pidiq_presets.csv file,
#' which is assumed to be found in the current working directory where the script is called.
#' @param msg Print progress-messages to console (default: FALSE).
#'
#' @return A data.frame with PIDIQ statistics calculated from the given plant image file.
#' @export
#'
#' @examples
#' \dontrun{pidiq_calc2()}
#'
pidiq_calc2 <- function(img_file,
                        output_dir = './test_pidiq2',
                        gsize_cut = NULL,
                        test_plt = FALSE,
                        col_lab = NULL,
                        row_lab = NULL,
                        filter = NULL,
                        msg = FALSE){


  ##Initialization of PIDIQ analysis run:

  #Set up a test directory for storing output:
  if(length(list.dirs('./test_pidiq2')) == 0) system(paste('mkdir', output_dir))


  #Check if image contains single or multiple plants:
  if(is.null(row_lab)){
    n_r <- 1
    n_c <- 1
    well_lab <- NULL
  }
  else{
    #Get the number of rows and columns from row_lab and col_lab:
    n_r <- length(row_lab)
    n_c <- length(col_lab)

    #Assign plate wells/crop areas to their corresponding letter-digit codes (row-wise):
    well_lab <- paste0(rep(row_lab, each = n_c), col_lab)
  }


  #Print an opening message:
  message('PIDIQv2.2:')

  if(n_r * n_c == 1) message('Single Image.') else message(paste('Multi-Well Image:', n_r, 'x', n_c))

  if(msg) message(rep('-', 25))


  #Preliminary step, make sure that spectral filtering parameters are provided:
  message(paste('Loading Spectral Filtering Defaults:', './pidiq_presets.csv'))

  filt_params <- list()

  if(file.exists('./pidiq_presets.csv')){
    #Load the presets file:
    filt_presets  <- utils::read.csv('./pidiq_presets.csv', header = TRUE, sep = '\t')
    #print(filt_presets)


    #Check if the parameter name exists in the presets file:
    if(!is.null(filter)){
      if(!filter %in% filt_presets$name)
        stop(paste0('"', filter, '"', ' preset not found.'), call. = FALSE)

      #If it does, then filter the presets table (note, all preset names must be unique!)

      #And split into 'healthy' and 'diseased' parameters and save to the filt_params list:
      filt_params <- split(subset(filt_presets, name == filter), ~ type)

      filt_params[[1]] <- filt_params[[1]][,-(1:2)]
      filt_params[[2]] <- filt_params[[2]][,-(1:2)]

      #Remember to normalize filtering parameters by 255 (if they haven't been already....):
      if(any(filt_params[[1]] > 1) | any(filt_params[[2]] > 1)){
        filt_params[[1]] <- filt_params[[1]] / 255
        filt_params[[2]] <- filt_params[[2]] / 255
      }

    }
    else{
      stop(paste('Please specify a spectral filtering preset.\n',
                 'Either one of:', paste(unique(filt_presets$name), collapse = ', ')),
           call. = FALSE)
    }

  }
  else{
    stop('...Spectral Filtering Presets File Cannot Be Found!', call. = FALSE)

  }

  #print(filt_params)


  #Step 1 Image Loading:

  #Check if the image path has been supplied or exists:
  if(!file.exists(img_file))
     stop(paste('Image File:', img_file, 'Cannot Be Found.'))

  message(paste('1) Loading Image:', img_file))


  #Load the image:
  img <- magick::image_read(img_file)

  #If the image width or height exceed 1000 pixels,
  #rescale to make it faster to work with:
  if(max(magick::image_info(img)[c('width', 'height')]) > 1000) img <- magick::image_scale(img, "1000")


  #Rasterize the image into a table: (x,y) of pixel coords and corresponding hex
  #colours:
  img2 <- magick::image_raster(img)

  #Convert the hex colour codes of the rasterized image to rgb and hsv/b:
  img2 <- cbind.data.frame(img2,
                           grDevices::col2rgb(img2$col) %>%
                             t() %>%
                             as.data.frame(),
                           grDevices::rgb2hsv(grDevices::col2rgb(img2$col)) %>%
                             t() %>%
                             as.data.frame())



  #####
  #Pixel Spectral Filtering:
  #####

  #Assign pixels which pass hsv cutoffs to 'yellow' and 'green' categories:
  if(msg) message('2) Performing Pixel Green/Yellow Filtering')


  t <- with(img2,
            ifelse(h >= filt_params[['healthy']][1, 'h_min'] & h <= filt_params[['healthy']][1, 'h_max'] &
                     s >= filt_params[['healthy']][1, 's_min'] & s <= filt_params[['healthy']][1, 's_max'] &
                     v >= filt_params[['healthy']][1, 'v_min'] & v <= filt_params[['healthy']][1, 'v_max'],
                   'green',
                   ifelse(h >= filt_params[['diseased']][1, 'h_min'] & h <= filt_params[['diseased']][1, 'h_max'] &
                            s >= filt_params[['diseased']][1, 's_min'] & s <= filt_params[['diseased']][1, 's_max'] &
                            v >= filt_params[['diseased']][1, 'v_min'] & v <= filt_params[['diseased']][1, 'v_max'],
                          'yellow',
                          NA
                          )
                   ))

  #Use the green/yellow spectral filtering to generate a binary image for input
  #into watershed and bw segmentation:


  #####
  #Segmentation:
  #####
  if(msg) message('3) Performing Image Segmentation')

  #Convert the spectral filtered image pixel array (t) back into a matrix, and convert it to the 'Image' class
  #for use with EBImage:

  #Note that x and y coordinates of a raster image correspond to the rows/columns of the
  #Image() matrix object, which are the image height and width, respectively. It is also important to note
  #that raster images are sorted column-wise:
  img3 <- EBImage::Image(matrix(ifelse(!is.na(t), TRUE, FALSE), nrow = max(img2$x), ncol = max(img2$y), byrow = FALSE))


  #Watershed Segmentation:
  nmask = EBImage::watershed( EBImage::distmap(img3),
                     tolerance = 1, ext = 1)

  #BW - Neighbouring Foreground Pixel Segmentation:
  bw <- EBImage::bwlabel(img3)



  #####
  #Detecting potential noise with segmentation group size histograms:
  #####

  #Calculate the sizes (n pixels) assigned per watershed group:
  grp_size <- data.frame(grp = matrix(nmask, ncol = 1),
                         npix = 1) %>%
    subset(grp != 0) %>%
    stats::aggregate(npix ~ grp, data = ., sum) %>%
    with(., array(.$npix, dimnames = list(.$grp)))


  #Same for bw groups:
  grp2_size <- data.frame(grp2 = matrix(bw, ncol = 1),
                          npix = 1) %>%
    subset(grp2 != 0) %>%
    stats::aggregate(npix ~ grp2, data = ., sum)

  grp2_size <- with(grp2_size, array(npix, dimnames = list(grp2)))



  #Add the spectral filtering pixel green/yellow assignments and segmentation
  #groups to the rasterized image data.frame:
  df <- data.frame(img2,
                   cat = matrix(t, ncol = 1),
                   #Watershed groups:
                   grp = matrix(nmask, ncol = 1),
                   #BW groups:
                   grp2 = matrix(bw, ncol = 1),
                   #Crop assignment:
                   crop = 0)


  #####
  #Check if group size cutoff filter has been set:
  #####

  if(is.null(gsize_cut)){
    message('###Seg. Group Size Threshold has not been set.')

    #Check distribution of segmentation group sizes:
    #watershed (npix) by bw (ordered by npix).

    g2_s <- table(df$grp2)
    g2_o <- array(1:length(g2_s), dimnames = list(names(g2_s)[order(g2_s)]))

    grp_size_dist <- df %>%
      subset(grp != 0 | grp2 != 0, select = c(grp, grp2)) %>%
      transform(npix = 1) %>%
      stats::aggregate(npix ~ grp2 + grp, data = ., length) %>%
      transform(g2_s = as.numeric(g2_s[as.character(grp2)])) %>%
      transform(ord = as.numeric(g2_o[as.character(grp2)])) %>%
      .[order(.$ord),]

    message('###Generating Seg. Group Size Distribution.')
    with(grp_size_dist,
         {
           plot(1:length(npix), sort(npix, decreasing = TRUE),
                col = grDevices::rainbow(max(grp2))[grp2],
                main = 'Segmentation Group Size Distribution',
                log = 'xy',
                cex = 0.5,
                pch = 16)
           #If the image contains multiple plants,
           #Add a vertical line roughly indicating the number of plants:
           if(n_r * n_c > 1)
            abline(v = n_r * n_c)
         }
    )


    while(is.null(gsize_cut)){

      opt <- readline('Set Segmentation Group Size Cutoff (Input Value Or Click On Plot = P/p)?\n')

      if(toupper(opt) == 'P'){
        gsize_cut <- graphics::locator(n = 1)
        gsize_cut <- floor(gsize_cut$y)

        message(paste('gsize_cut set to:', gsize_cut))

      }
      else if(is.na(as.numeric(opt)) | as.numeric(opt) < 0){
        print(as.numeric(opt))
        message('Value must be numeric and greater than or equal to 0.')
        gsize_cut <- NULL

        #stop()
      }
      else{
        gsize_cut <- as.numeric(opt)
      }
    }
  }


  #Here is where background noise (non-plant pixels passing spectral filtering)
  #can be excluded by assigning bw and watershed segmentation groups falling
  #under a specified size to background (= 0).

  #Note that the segmentation group size threshold for filtering will depend on
  #the overall resolution of the image....

  df <- df %>%
    transform(grp2 = ifelse(grp2 %in% as.numeric(names(grp2_size[grp2_size <= gsize_cut])), 0, grp2),
              grp = ifelse(grp %in% as.numeric(names(grp_size[grp_size <= gsize_cut])), 0, grp)) %>%
    transform(cat = ifelse(grp == 0 | grp2 == 0, NA, cat)#,
              # col = ifelse(grp2 %in% as.numeric(names(grp2_size[grp2_size <= gsize_cut])) |
              #              grp %in% as.numeric(names(grp_size[grp_size <= gsize_cut])),
              #              'white', col)
              )

  #After this step, single-image PIDIQ is finished.

  ##Perform Single or Multi-Well PIDIQ##
  if(is.null(row_lab) | n_c * n_r == 1){
    pidiq_single(img_file, output_dir, msg, df)
  }
  else{
    pidiq_multi(img_file, output_dir, well_lab, n_r, n_c, msg, test_plt, df, grp_size, grp2_size)
  }

}


#Multi-Well PIDIQ:
pidiq_multi <- function(img_file, output_dir, well_lab, n_r, n_c, msg, test_plt, df, grp_size, grp2_size){

  #####
  #Subdivide the image into individual crops and wells:
  #####

  if(msg) message('4) Determining Image Crop Regions')

  #Get average crop/well pixel height and width:
  df_row_l = floor(max(df$y) / n_r)
  df_col_l = floor(max(df$x) / n_c)

  #Get breaks between the crop/wells:
  df_row_br <- seq(0, max(df$y), df_row_l)
  df_row_br[length(df_row_br)] <- max(df$y)

  df_col_br <- seq(0, max(df$x), df_col_l)
  df_col_br[length(df_col_br)] <-max(df$x)


  #####
  #Get the dimensions of each crop:
  #####
  #y-axis/row well limits (invert these afterward for consistency with visualizations?):
  r_lims <- sapply(1:(length(df_row_br) - 1), function(i) df_row_br[i:(i+1)])

  #x-axis/column well limits:
  c_lims <- sapply(1:(length(df_col_br) - 1), function(i) df_col_br[i:(i+1)])


  #To make the ordering of the crops (mappings) consistent with their labels for
  #the well/circles (see df_c_cat), extract limits row-wise:

  #Define extraction indices for the c_lims and r_lims matrices:
  #Note: Wells are numbered row-wise:
  r_ind <- rep(1:n_r, each = n_c)
  c_ind <- rep(1:n_c, length.out = n_c * n_r)

  crop_dim <- data.frame(well = 1:(n_r * n_c),
                         x_min = c_lims[1, c_ind] + 1,
                         x_max = c_lims[2, c_ind],
                         y_min = r_lims[1, r_ind] + 1,
                         y_max = r_lims[2, r_ind])


  #####
  #Assign pixels to crops:
  #####
  #Initialize a matrix with the same dimensions of the original
  #image, and map entries according to well #, then flatten into an 1-d array.
  #Note: matrix is flattened column-wise, so to match up the
  #corresponding crop indices with x-y image pixel raster
  #coordinates in df, the well-matrix indices dimensions have to be reversed,
  #i.e. rows = width/x, columns = height/y

  #Note: x = cols / width, y = rows / height.
  m <- matrix(0, max(df$x), max(df$y))

  for(i in 1:nrow(crop_dim)){
    m[crop_dim$x_min[i] : crop_dim$x_max[i],
      crop_dim$y_min[i] : crop_dim$y_max[i]] <- crop_dim$well[i]
  }

  df$crop <- array(matrix(m, ncol = 1))


  #Data for calculating well/circle crops based on the grids:
  #*Will be used for well-radius group assignment and plotting.
  df_c_dat <- data.frame(well = 1:(n_c * n_r),
                         a = rep(df_col_br[-length(df_col_br)] + diff(df_col_br)/2, times = length(df_row_br) - 1),
                         b = rep(df_row_br[-length(df_row_br)] + diff(df_row_br)/2, each = length(df_col_br) - 1),
                         r = ifelse(min(diff(df_row_br)) < min(diff(df_col_br)), min(diff(df_row_br))/2, min(diff(df_col_br))/2)
  )


  #Find midpoints of segmentation groups (x_mid, y_mid) by finding the median of the ranges of x and y
  #pixel coordinates which delimit each group.

  grp_mids <- df %>%
    subset(grp != 0, select = c(grp, x, y)) %>%
    #Note, the groupings for x and y are stored as separate matrices, have to
    #reformat the final data.frame:
    stats::aggregate(. ~ grp, data = .,
              function(x) return(c(npix = length(x),
                                   min = min(x),
                                   max = max(x),
                                   mid = stats::median(c(min(x), max(x)))
                                   )
              )
    )

  grp_mids <- data.frame(grp_mids$grp,
                         x=grp_mids$x[,-4],
                         y = grp_mids$y[,-c(1,4)],
                         mid = cbind(grp_mids$x[,4], grp_mids$y[,4]))

  colnames(grp_mids) <- c('grp', 'npix', 'xmin', 'xmax', 'ymin', 'ymax', 'x_mid', 'y_mid')



  ####
  #Plot the segmentation groups (pre-crop assignment):
  #Put this into a code block, and assign a parameter to make plotting optional.
  ####

  if(test_plt){

    if(msg) message('# Generating Test Plots')

    watershed_file = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_seggrp_orig_watershed.png')
    bw_file = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_seggrp_orig_bw.png')

    if(msg) message(paste('\t###', watershed_file, "\n"), paste('\t###', bw_file))

    with(df, #%>%
           #Filter out watershed segmentation groups with <= 200 pixels?
           #(These will appear as 'white' regions on the watershed group plot):
           #transform(cat = ifelse(grp %in% as.numeric(names(grp_size[grp_size <= 200])), NA, cat)),
         {
           #Initialize matrices for raster-plotting:

           #Palettes for groups:
           set.seed(1234)
           grp_col <- c(NA, sample(grDevices::rainbow(max(df$grp))))

           grp2_col <- c(NA, sample(grDevices::rainbow(max(df$grp2))))

           #Watershed segmentation group crop assignments:
           img_ws <- matrix(ifelse(!is.na(cat), grp_col[grp + 1], col),
                            nrow = max(y), ncol = max(x), byrow = TRUE)

           #BW segmentation group crop assignments:
           img_bw <- matrix(ifelse(!is.na(cat), grp2_col[grp2 + 1], col),
                            nrow = max(y), ncol = max(x), byrow = TRUE)


           ####
           #Plot images:
           ####

           #1) Watershed Segmentation Groups:

           #Open a file for plotting:
           png(watershed_file, width = 10, height = 10, units = 'in', res = 72)

           #Initialize the plot:
           par(bty = 'n', mar = c(0,0,0,0), xpd = NA)
           plot(c(min(x), max(x)), c(min(y), max(y)), type = 'n', xaxt = 'n', yaxt = 'n', axes = 0)

           #Add spectral-filtered green/yellow pixels:
           rasterImage(img_ws,
                       xleft = min(x), ybottom = min(y), xright = max(x), ytop = max(y),
                       interpolate = FALSE)

           #Plot group midpoints:
           with(grp_mids,
                {
                  points(x_mid, max(df$y) - y_mid, col = 'white', pch = 16, cex = 1)
                  points(x_mid, max(df$y) - y_mid, col = 'purple', pch = 16, cex = 0.6)
                })

           #Close the file:
           dev.off()


           #2) BW Segmentation Groups:

           #Open a new file for plotting:
           #outfile = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_seggroups.png')
           png(bw_file, width = 10, height = 10, units = 'in', res = 72)

           #Initialize the plot (un-comment the following lines to produce a separate plot):
           par(bty = 'n', mar = c(0,0,0,0), xpd = NA)
           plot(c(min(x), max(x)), c(min(y), max(y)), type = 'n', xaxt = 'n', yaxt = 'n', axes = 0)

           #Add segmentation groups (merged bw and watershed groups with updated crop assignments):
           rasterImage(img_bw,
                       xleft = min(x), ybottom = min(y), xright = max(x), ytop = max(y),
                       interpolate = FALSE)

           #Close the file:
           dev.off()

         })
    #####
  }



  ###
  #Assignment of spectral filtered pixels to specific crops/plants:
  ###

  #First check if bw groups can be identified which are unique to different crop regions,
  #I.e. if > 75% of their pixels fall into a respsective crop, or by the crop where the
  #maximum number of the group's pixels can be found.

  if(msg) message('5) Identifying Segmentation Groups Unique to Crops:')

  #Calculate the Sum the number of foreground pixels in each crop area:
  crop_size <- df %>%
    #Do not include background pixels:
    subset(grp2 != 0, select = crop) %>%
    transform(n = 1) %>%
    stats::aggregate(n ~ crop, data = ., sum)

  crop_size <- with(crop_size, array(n, dimnames = list(crop)))

  grp2_crop <- df %>%
    #Remove background bw groups, not passing spectral filtering:
    subset(grp2 != 0, select = c(grp2, crop)) %>%
    cbind.data.frame(., n_pix = 1) %>%
    stats::aggregate(n_pix ~ crop + grp2, data = ., sum) %>%
    #Filter out bw groups with <= 200 pixels (Could do it as a separate step above):
    #subset(n_pix > 200) %>%
    #Similar to the above, but group by crops instead:
    split(., as.factor(.$crop)) %>%
    #Identify the groups which are exclusively associated with the crop (> 50 % of pixels),
    #Or have the greatest overall coverage of the crop:
    lapply(., function(x) return(cbind.data.frame(x,
                                                  prop_crop = x$n_pix / crop_size[as.character(x$crop)],
                                                  prop_grp2 = x$n_pix / grp2_size[as.character(x$grp2)]))) %>%
    lapply(., function(x) return(x[x$prop_grp2 >= 0.75 | x$prop_crop >= 0.5, ])) %>%
    #lapply(., function(x) return(x[which.max(x$n_pix),])) %>%
    do.call(what = 'rbind.data.frame', args = .)


  grp2_summary <- grp2_crop %>%
    #Summarize the number of crops in which each group is found:
    by(., .$grp2, function(.) {
      d <- data.frame(grp2 = unique(.$grp2),
                      n_crop = length(unique(.$crop)))
      d$crops = list(unique(.$crop))
      return(d)
    }) %>%
    do.call('rbind.data.frame', .)


  #Extract groups uniquely assigned to crops:
  grp2_uni <- grp2_summary %>%
    subset(n_crop == 1) %>%
    transform(crops = unlist(crops)) %>%
    .[order(.$crops),]

  #Note that sometimes an over-grown group might be misassigned as the dominant one in a
  #crop, because the actual plant in the crop has a smaller area of plant tissue present. This
  #can be corrected downstream when examining watershed groups, by first making
  #unique crop assignments based on the majority crop pixel proportion of their encompassing
  #bw groups, and finding the segmentation group to split on.


  if(msg) message(paste('\t###', nrow(grp2_uni), 'Groups Assigned To', length(unique(grp2_uni$crop)), 'Crops.'))



  #Next, check if any bw groups span multiple crops.
  #If they do, set the groups aside for reassignment by
  #checking the crop localizations of watershed groups
  #encompassed by them.

  if(msg) message('6) Identifying Segmentation Groups Merged Across Crops:')


  #Extract groups assigned/spanning multiple crops:
  grp2_ol <- grp2_summary %>%
    subset(n_crop != 1)


  if(msg) message(paste('\t###', nrow(grp2_ol), 'Segmentation Groups Merged Across', length(unlist(grp2_ol$crops)), 'Crops.'))


  #Initialize an array for storing updated watershed group crop assignments
  #for merged bw groups:
  grp_crop_map = NA


  ###Check if any merged segmentation groups exist:
  if(nrow(grp2_ol) != 0){

    if(msg) message('7) Splitting Merged Segementation Groups')


    #Now, for each bw group overlapping multiple crops,
    #partition into distinct crops based on their constituent
    #watershed groups.

    #Examine its corresponding watershed segmentation groups,
    #Identify what proportion of each group's pixels belongs to
    #a given crop area.

    #Instances to check for:
    ##Watershed groups with pixels fall across multiple crop regions:
    ##1) These could overlap or extend into the borders of a neighbouring crop that is
    ##not included by the parent bw segmentation group, or;
    ##2) Lead in to a neighoubring crop included by the bw segmentation group.

    #After identifying these suspect watershed groups,
    #Reassign them based on nearest-neighbouring groups
    #that are uniquely assigned to a particular crop.

    #Use these watershed groups to then update the cropping assignments
    #of pixels in the image.

    #Extract  bw ('parent') segmentation groups with > 75%? of pixels
    #assigned to a given crop, and transfer these to ambiguously assigned
    #watershed group;
    #Note, also make sure that the bw group is not the predominant within the crop,
    #prop_crop < 50%.
    uni_grp2 <- with(grp2_crop,
                     {
                       i <- prop_grp2 >= 0.8 & prop_crop < 0.5
                       array(crop[i], dimnames = list(grp2[i]))
                     })

    #Do the equivalent for watershed groups:
    uni_grp <- df %>%
      subset(grp != 0, select = c(grp, grp2, crop)) %>%
      transform(n_pix_crop = 1) %>%
      #Tally the number of pixels of each watershed group by
      #each distinct crop they fall into:
      stats::aggregate(n_pix_crop ~ grp + crop, data = ., FUN = sum) %>%
      #Now calculate the total proportion of watershed group pixels
      #found in each crop region:
      transform(grp_prop = n_pix_crop / grp_size[as.character(grp)]) %>%
      #Keep all watershed groups which have >= 80% pixels in a specific crop-region,
      #and return as a grp-crop mapping array:
      subset(grp_prop >= 0.8) %>%
      with(.,
           array(crop, dimnames = list(grp)))



    grp_ord <- df %>%
      subset(grp2 %in% grp2_ol$grp2, select = c(grp2, grp, crop)) %>%
      #Reassign any watershed groups which belong to a bw group with
      # > 90% pixels assigned to a given crop:
      transform(crop = ifelse(!is.na(uni_grp2[as.character(grp2)]),
                              as.numeric(uni_grp2[as.character(grp2)]),
                              crop)) %>%
      #(Could also check the proportion of a watershed group's pixels
      #across crops, and assign based on a proportional cutoff) -
      #But this results in problems with significantly overgrown wells...
      # transform(crop = ifelse(!is.na(uni_grp[as.character(grp)]),
      #                         as.numeric(uni_grp[as.character(grp)]),
      #                         crop)) %>%
      #Reduce dataframe down to the number of unique crop assignments
      #made for watershed groups:
      unique() %>%
      #Transform the data.frame to wide-format:
      by(., .$grp, function(.) {
        d <- data.frame(grp = unique(.$grp),
                        grp2 = unique(.$grp2),
                        n_crop = length(unique(.$crop)))
        d$crops = list(unique(.$crop))
        return(d)
      }
      ) %>%
      do.call('rbind.data.frame', .) %>%
      #Add mid-points of groups, for downstream visualization:
      merge(., grp_mids, by = 'grp') %>%
      #Check each watershed segmentation group against the grp2_ol data.frame to
      #identify whether they fall between/merge bw crops (>= 2) or fall outside of the
      #crops (0) (sel == 0), Otherwise they fall inside the crop (sel == 1), Note,
      #make sure to only check the crops included by the grp2 (bw-group)
      #assignment of each watershed segmentation group!:
      cbind.data.frame(.,
                       sel = sapply(1:nrow(.), function(i){
                         j <- which(unlist(.$crops[i]) %in% unlist(grp2_ol$crops[grp2_ol$grp2 == .$grp2[i]]))
                         ifelse(length(j) >= 2 | length(j) == 0, 0, 1)
                       })
      ) %>%
      #Sort by group x and y coordinate midpoints to
      #enable checking of neighbouring groups below:
      .[order(.$x_mid, .$y_mid),]


    #For the calculated midpoints of the watershed groups,
    #assign them to their corresponding BW groups:
    #
    #Identify groups which haven't been uniquely assigned to
    #a specific crop, or if they belong to a crop different from the
    #crop assigned to the respective bw group:
    #Reassign them based on the crop assignments of the nearest groups:


    if(msg) message('8) Reassigning Segmentation Groups to Crops')

    g_check <- grp_ord %>%
      subset(n_crop >= 2 | sel == 0, select = c(grp, grp2, sel))


    crop_new <- NULL

    if(nrow(g_check) != 0){

      for(g in g_check$grp){

        g_i <- grp_ord %>%
          subset(grp == g)

        #Identify the nearest watershed groups to the current one:
        check <- grp_ord %>%
          subset(#npix > 200 &
            grp2 %in% g_i$grp2) %>%
          transform(d_x = g_i$x_mid - x_mid,
                    d_y = g_i$y_mid - y_mid) %>%
          transform(d = sqrt(d_x^2 + d_y^2)) %>%
          subset(n_crop == 1 &
                   #Exclude neighbouring watershed groups
                   #that have to be reassigned:
                   !grp %in% g_check$grp &
                   #And only check segmentation groups
                   #which belong to major-crops merged by the
                   #parent BW group where they are present:
                   crops %in% (grp2_ol %>%
                                 subset(grp2 %in% g_i$grp2) %>%
                                 .$crops %>%
                                 unlist())) %>%
          transform(crops = unlist(crops)) %>%
          .[order(.$d),]

        #Check the crop assignments of the three nearest watershed groups,
        #If they are all the same, then assign it to the suspect group:
        if(length(unique(check$crops[1:3])) == 1){
          crop_new <- append(crop_new, check$crops[1])
        }
        else{
          crop_new <- append(crop_new, NA)
        }

      }
    }

    #Make updated group assignments to the remaining watershed groups:
    crop_update <- rbind.data.frame(data.frame(grp = g_check$grp, crop_new),
                                    #Add watershed groups delimited by bw groups which overlap
                                    #multiple wells but uniquely belong to a single crop:
                                    grp_ord %>%
                                      subset(sel == 1 & n_crop == 1) %>%
                                      transform(crop_new = unlist(crops)) %>%
                                      subset(select = c(grp, crop_new))
    )


    #Create updated crop-assignment mapping array for watershed groups:
    grp_crop_map <- with(crop_update, array(crop_new, dimnames = list(grp)))

  }


  #Create updated crop-assignment mapping array for bw groups:
  grp2_crop_map <- with(grp2_uni %>%
                          #Remember to remove any merged bw groups if they exist:
                          subset(!grp2 %in% grp2_ol$grp2,
                                 select = c(grp2, crops)),
                        array(crops, dimnames = list(grp2)))



  #Map updated crop assignments to segmented pixel groups from the original
  #plate image:

  if(msg) message('9) Generating Final Pixel Crop Assignments')

  plt <- df %>%
    transform(crop_new.x = grp2_crop_map[as.character(grp2)],
              crop_new.y = grp_crop_map[as.character(grp)]) %>%
    transform(crop_new = ifelse(!is.na(crop_new.x), crop_new.x,
                                ifelse(!is.na(crop_new.y), crop_new.y, NA))) %>%
    transform(crop_col = ifelse(!is.na(crop_new),
                                grDevices::rainbow(n_r * n_c)[crop_new],
                                col)) %>%
    transform(cat_col = ifelse(!is.na(crop_new),
                               ifelse(cat == 'yellow', 'red', cat), col)) #%>%
  # transform(crop_col = as.numeric(str_split(well_lab[crop_new], '', simplify = TRUE)[,2]),
  #           crop_row = str_split(well_lab[crop_new], '', simplify = TRUE)[,1])



  ###
  #The plots:
  ###

  if(msg) message('10) Generating Output Images:')

  #Output image for green/yellow pixel assignments:
  outfile1 = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_painted.png')

  #Output image for final segmentation group assignments:
  outfile2 = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_seggrp_final.png')

  if(msg) message('\t###', outfile1)
  if(msg) message('\t###', outfile2)


  #The plot, Base-R version:
  #Visual check:
  with(plt,
       {
         #Initialize matrices for raster-plotting:

         #Spectral Filtered Pixels, sing pidiq group assignments (green == 'green', yellow = 'red'):
         img_cat <- matrix(cat_col,
                           nrow = max(y), ncol = max(x), byrow = TRUE)

         #Segmentation group crop assignments:
         img_seggrp <- matrix(crop_col,
                              nrow = max(y), ncol = max(x), byrow = TRUE)


         ####
         #Plot images:
         ####

         #1) Plot Green/Yellow pixels:

         #Open a file for plotting:
         #outfile = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_painted.png')
         png(outfile1, width = 10, height = 10, units = 'in', res = 72)

         #Initialize the plot:
         par(bty = 'n', mar = c(0,0,0,0), xpd = NA)
         plot(c(min(x), max(x)), c(min(y), max(y)), type = 'n', xaxt = 'n', yaxt = 'n', axes = 0)

         #Add spectral-filtered green/yellow pixels:
         rasterImage(img_cat,
                     xleft = min(x), ybottom = min(y), xright = max(x), ytop = max(y),
                     interpolate = FALSE)

         #Close the file:
         dev.off()


         #2) Plot the 'final' segmentation groups:

         #Open a new file for plotting:
         #outfile = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_seggroups.png')
         png(outfile2, width = 10, height = 10, units = 'in', res = 72)

         #Initialize the plot (un-comment the following lines to produce a separate plot):
         par(bty = 'n', mar = c(0,0,0,0), xpd = NA)
         plot(c(min(x), max(x)), c(min(y), max(y)), type = 'n', xaxt = 'n', yaxt = 'n', axes = 0)

         #Add segmentation groups (merged bw and watershed groups with updated crop assignments):
         rasterImage(img_seggrp,
                     xleft = min(x), ybottom = min(y), xright = max(x), ytop = max(y),
                     interpolate = FALSE)

         #Annotate the final crop assignments of the segmentation groups:
         crp_lab <- plt %>%
           subset(!is.na(crop_new), select = c(crop_new, x, y)) %>%
           stats::aggregate(cbind(x, y) ~ crop_new, data = ., stats::median) %>%
           transform(y = max(plt$y) - y)

         with(crp_lab,
              text(x, y, paste0(crop_new, ' : ', well_lab[crop_new]), col = 'white'))

         #Close the file:
         dev.off()

       })



  #Output PIDIQ summary statistics for each crop/well/plant in the plate:

  if(msg) message('11) Calcuating PIDIQ Summary Stats')


  #Also calculate the average area of each crop to determine roughly whether
  #plants are overgrown or undergrown on screened plates:
  crop_area <- floor(max(df$x) * max(df$y) / n_r * n_c)

  #Instead of square-cropped area, how about circular?
  well_area <- floor(pi * (min(max(df$x) / n_c, max(df$y) / n_r) / 2 ) ^2)


  #Make new crop assignments to the spectral filtered pixels,
  #and calculate pidiq summary stats:
  res <- plt %>%
    #Extract all of the segmentation group pixels assigned to different crops:
    subset(!is.na(crop_new),
           select = c(crop_new, cat)) %>%
    #Summarize the number of green and yellow pixels assigned to each crop:
    transform(npix = 1) %>%
    stats::aggregate(npix ~ crop_new + cat, data = ., sum) %>%
    #Pivot data.frame to wide format so that the number of green and yellow pixels
    #per crop are assigned to individual columns:
    stats::reshape(idvar = 'crop_new', timevar = 'cat', direction = 'wide') %>%
    #Calculate PIDIQ summary statistics:
    transform(npix.yellow = ifelse(is.na(npix.yellow), 0, npix.yellow),
              npix.green = ifelse(is.na(npix.green), 0, npix.green)) %>%
    transform(Well = well_lab[crop_new],
              Well_Sum = npix.yellow + npix.green,
              Well_Prop = (npix.yellow + npix.green) / well_area,
              Well_RawYellow = npix.yellow / (npix.yellow + npix.green)) %>%
    transform(Well_ArcYellow = asin(sqrt(Well_RawYellow)) / asin(1)) %>%
    #Add any potentially missing Wells:
    merge(data.frame(Well = well_lab), all.y = TRUE) %>%
    #Add column for processed input file name:
    cbind.data.frame(data.frame(File = img_file), .) %>%
    #Remove the crop_new column:
    subset(select = -crop_new)

  #Rename the green and yellow pixel sum columns:
  colnames(res)[3:4] <- c('Well_GreenArea', 'Well_YellowedArea')

  if(msg) message('-----')
  if(msg) message('Done.')

  #Output results table:
  return(res)
}




pidiq_single <- function(img_file, output_dir, msg, df){
  #Map updated crop assignments to segmented pixel groups from the original
  #plate image:

  if(msg) message('4) Generating Output Images:')

  #Generate final plot image by removing any pixels which passed spectral
  #filtering but did not pass segmentation group size thresholding:
  plt <- df %>%
    transform(cat_col = ifelse(grp != 0 | grp2 != 0,
                               ifelse(cat == 'yellow', 'red', cat), col))


  ###
  #The plot:
  ###

  #Output image for green/yellow pixel assignments:
  outfile1 = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_painted.png')

  if(msg) message('\t###', outfile1)

  #The plot, Base-R version:
  #Visual check:
  with(plt,
       {
         #Initialize matrices for raster-plotting:

         #Spectral Filtered Pixels, sing pidiq group assignments (green == 'green', yellow = 'red'):
         img_cat <- matrix(cat_col,
                           nrow = max(y), ncol = max(x), byrow = TRUE)

         ####
         #Plot image:
         ####

         #1) Plot Green/Yellow pixels:

         #Open a file for plotting:
         #outfile = paste0(output_dir, '/', sub('.*[/]', '', img_file), '_painted.png')
         png(outfile1, width = 10, height = 10, units = 'in', res = 72)

         #Initialize the plot:
         par(bty = 'n', mar = c(0,0,0,0), xpd = NA)
         plot(c(min(x), max(x)), c(min(y), max(y)), type = 'n', xaxt = 'n', yaxt = 'n', axes = 0)

         #Add spectral-filtered green/yellow pixels:
         rasterImage(img_cat,
                     xleft = min(x), ybottom = min(y), xright = max(x), ytop = max(y),
                     interpolate = FALSE)

         #Close the file:
         dev.off()

       })


  #Output PIDIQ summary statistics for each crop/well/plant in the plate:

  if(msg) message('5) Calcuating PIDIQ Summary Stats')

  #Calculate pidiq summary stats:
  res <- plt %>%
    #Extract all of the segmentation group pixels assigned to different crops:
    subset(!is.na(cat),
         select = cat) %>%
    #Summarize the number of green and yellow pixels assigned to each crop:
    transform(File = img_file, npix = 1) %>%
    stats::aggregate(npix ~ cat + File, data = ., sum) %>%
    #Pivot data.frame to wide format so that the number of green and yellow pixels
    #per crop are assigned to individual columns:
    stats::reshape(idvar = 'File', timevar = 'cat', direction = 'wide') %>%
    #Calculate PIDIQ summary statistics:
    transform(npix.yellow = ifelse(is.na(npix.yellow), 0, npix.yellow),
              npix.green = ifelse(is.na(npix.green), 0, npix.green)) %>%
    transform(Sum = npix.yellow + npix.green,
              RawYellow = npix.yellow / (npix.yellow + npix.green)) %>%
    transform(ArcYellow = asin(sqrt(RawYellow)) / asin(1))

  #Rename the green and yellow pixel sum columns:
  colnames(res)[1:2] <- c('Well_GreenArea', 'Well_YellowedArea')

  if(msg) message('-----')
  if(msg) message('Done.')

  #Output results table:
  return(res)
}

