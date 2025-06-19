# load libraries (you might need to install them first - use install.packages("xxx") command line)
#install.packages("imager")
#install.packages("pbapply")
#install.packages("data.table")
#install.packages("readxl")
#install.packages("EBImage")
#install.packages("BiocManager")
#devtools::install_github('omarwagih/gitter', force=T)

library("imager")
library("pbapply")
library("data.table")
library("readxl")
library("gitter")


# clear workspace
rm(list=ls())

# SET UP FILES and FOLDER
#########################
# folders need to be set up as provided ie: background folder, unwashed folder, washed folder and analysis folder
# the plate image files need to be in the washed and unwashed folders
# the background images should be in the background folder and systematically named similarly as washed and unwashed images
# image files should be name as follow:
	# for unwashed plates: BHET25mM_24h_plate_01.JPG
	# for washed plates: BHET25mM_24h_washed_plate_01.JPG
	# for background plates: BHET25mM_24h_background_plate_01.JPG

# plate numbering needs to be the same for background, wash and unwashed images

# SET UP MAIN FOLDER
####################
# ENTER PATH to main folder containing the analysis, unwashed and washed folders
	
	folder0 = "/Users/tania/Desktop/PETase\ project\ 2/Tania_imaging/2025-06-06"

			# THE FOLLOWING ARE AUTOMATICALLY DEFINED FROM folder0	
				# folder containing all the analysis results
				analysis.folder = paste0(folder0, "/analysis")
		
				# folder with background images
				folder00 = paste0(folder0, "/background")

				# folder with unwashed images  
				folder1 = paste0(folder0, "/unwashed")
		
				# folder with washed plates
				folder2 = paste0(folder0, "/washed")
		
				# template image for gridding with gitter (OPTIONAL)
				# template = paste0(folder1, "/TEMPLATE.JPG")


# PLATE map: if there is a plate map, the file should in an excel format and in the analysis folder
###################################################################################################
	# the map needs to contains 5 columns name as follow and in the following order: p, c, r, ORF (gene systematic name e.g YAL001C), gene (gene standard name)
	# VERY IMPORTANT: the map needs to be organized by ascending plate and rows and columns: p = 1, 1, 1, 1, 1, r = 1, 1, 1, 1, 1 and c = 1, 2, 3, 4, 5
	# if in doubt: see example map file in the analysis folder

	# provide path name to strain map in the code below
	
		setwd(analysis.folder)
		
		strain.map = as.data.frame(read_excel("plate_map.xlsx", sheet = 1))
		strain.map = strain.map[order(strain.map$p, strain.map$r, strain.map$c),]
		
		setwd(folder0)

	
# SET UP ARRAY SIZE
###################
# number of colonies on array (96, 384 or 1536) and array factor

	array.factor = 2			# 1 for 96, 2 for 384, 4 for 1536
	
				ncol.array = 12 * array.factor
				nrow.array = 8 * array.factor
	

##################### the code below does not need to be changed #####################


# QUANTIFY COLONY SIZE USING GITTER AND IDENTIFY GRID
#####################################################

# MODIFY GITTER FUNCTION SO THAT THE .DAT FILE CONTAINS THE GRID COORDINATES
	# run the following command to open the script editing window
	trace(gitter, edit=T) 

	# replace line 232 with the following - make sure to remove the # sign or the line won't run
    # results = results[,]
	# save the changes in the editing window and close

	setwd(analysis.folder)
		
	gitter.batch(folder1, plate.format=c(nrow.array, ncol.array), verbose = "p", grid.save = folder1, dat.save = folder1)

	# template = "path to ref image" - define this if a reference image is needed to identify colonies properly
	# gitter.batch(folder1, template, plate.format=c(nrow.array, ncol.array), verbose = "p", grid.save = folder1, dat.save = folder1)	# code for gitter using a template image for colony identification

	colonySize.ls = dir(folder1)[grep('.dat', dir(folder1))]

	# cbind all colony size data and annotate with plate names

		# generate plate IDs and colony coordinates
		
		plate = dir(folder1)
		plate = plate[grep(".JPG$", plate)]
		plate = plate[-grep("gridded", plate)]
		plate = sort(rep(plate, nrow.array * ncol.array))

		colNames = c("plate", "size", "circularity", "flags", "x", "y", "xl", "xr", "yt", "yb")

		setwd(folder1)

		# cbind plate map and colony size
		
		colonySize = cbind(plate , do.call(rbind,  lapply(colonySize.ls, function(i) read.table(i, sep = "\t"))))
		colnames(colonySize) = c("plate", "r", "c", "size", "circularity", "flags", "x", "y", "xl", "xr", "yt", "yb")

		setwd(analysis.folder)


# QUANTIFY HALO	
###############

	# generate list of background images
	backgdImages = dir(folder00)[grep('.JPG$', dir(folder00))]

	# generate list of halo images
	haloImages = dir(folder2)[grep('.JPG$', dir(folder2))]

	# generate list of mask images
	maskImages = dir(folder1)[grep('gridded', dir(folder1))]
	
	# generate list of colony coordinate files
	colonyCoord = dir(folder1)[grep('.dat', dir(folder1))]
	
	# run quantification pipeline
	source(paste0(analysis.folder, "/halo_quantif_function_v3.R"))

	# write data in a csv file
	write.csv(Halo.data, "Halo.data.csv", row.names=F)


# COMBINE DATA WITH PLATE MAP 
#############################

		setwd(analysis.folder)

		Halo.data2 = cbind(strain.map[,], Halo.data)
		colnames(Halo.data2) = c(colnames(strain.map) , colnames(Halo.data))
		
		write.csv(Halo.data2, "halo.data.med.csv", row.names=F)








