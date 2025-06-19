# function for quantifying halo 

Halo.data = do.call(rbind, pblapply(1:length(haloImages), function(i){

		setwd(folder00)
		BackgdImage = load.image(backgdImages[i])

		setwd(folder2)
		HaloImage = load.image(haloImages[i])
		
		setwd(folder1)
		MaskImage = load.image(maskImages[i])
		
		colony.coord = read.table(colonyCoord[i], sep = "\t")
		
		# make df with binary values for mask image
			MaskImage.df = as.data.frame(grayscale(MaskImage))
			MaskImage.df[,"value"] = ifelse(MaskImage.df[,"value"] < 0.6, 0, 1)
		
		# convert HaloImage and BackgdImage to grayscale and make dataframe 

		Gray.HaloImage.df = as.data.frame(grayscale(HaloImage))
		Gray.BackgdImage.df = as.data.frame(grayscale(BackgdImage))

		# correct Gray.HaloImage.df pixel intensity at each position by the pixel intensity in the Gray.BackgdImage.df
		
		Gray.HaloImage.df[,] = cbind(Gray.HaloImage.df[,1:2], value = Gray.HaloImage.df["value"] - Gray.BackgdImage.df[,"value"])
		
		# cbind values halo image and mask
		# convert to data table
		# the reason for converting into data.table is because it is much faster to subset in the loop (lapply) below

		Halo.Mask.dt = as.data.table(cbind(Gray.HaloImage.df, mask= MaskImage.df[,"value"]))
			
			# compute intensity per colony
			# this involve removing 0 values because they shouldn't be part of the colony intensity
			# the lapply below compute mean and median intensity for each colony - using horizontal compute (across columns first then across rows)
		
			image.name = haloImages[i]
			
			colony.intensity = do.call(rbind, pblapply(1:nrow(colony.coord), function(j){
				
				left.col.index = colony.coord[j,8]
				right.col.index = colony.coord[j,9]
				top.row.index = colony.coord[j,10]
				bottom.row.index = colony.coord[j,11]
			
				colony.grid = Halo.Mask.dt[x %in% left.col.index:right.col.index & y %in% top.row.index:bottom.row.index]
								
				#x.expand = c(colony.grid[mask == 1][,x]-ring.px, colony.grid[mask == 1][,x]+ring.px)
				#y.expand = c(colony.grid[mask == 1][,y]-ring.px, colony.grid[mask == 1][,y]+ring.px)
				
				# measure intensity of the colony grid 
				grid.value = colony.grid[,median(value)] 
				
				# measure intensity of the halo under colony
				colony.value = colony.grid[,median(value[mask == 1])]
											
				as.data.frame(cbind(image.name, grid.value, colony.value))
	
					}))

			}))

setwd(analysis.folder)

Halo.data[,"grid.value"] = as.numeric(Halo.data[,"grid.value"])
Halo.data[,"colony.value"] = as.numeric(Halo.data[,"colony.value"])

Halo.data = cbind(image.name = colonySize[,1], colonySize[2:ncol(colonySize)], Halo.data[,2:ncol(Halo.data)])





