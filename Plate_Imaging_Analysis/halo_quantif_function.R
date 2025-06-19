# function for quantifying halo 

Halo.data = do.call(rbind, pblapply(1:length(haloImages), function(i){

		setwd(folder00)
		BackgdImage = load.image(backgdImages[i])

		setwd(folder2)
		HaloImage = load.image(haloImages[i])
		
		setwd(folder1)
		MaskImage = load.image(maskImages[i])
		
		# make df with binary values for mask image
			MaskImage.df = as.data.frame(grayscale(MaskImage))
			MaskImage.df[,"value"] = ifelse(MaskImage.df[,"value"] < 0.6, 0, 1)
		
		# convert HaloImage and BackgdImage to grayscale and make dataframe 

		Gray.HaloImage.df = as.data.frame(grayscale(HaloImage))
		Gray.BackgdImage.df = as.data.frame(grayscale(BackgdImage))

		# cbind values from background and halo image and mask
		# convert to data table
		# the reason for converting into data.table is because it is much faster to subset in the loop (lapply) below

		Halo.Mask.dt = as.data.table(cbind(Gray.BackgdImage.df, halo.value = Gray.HaloImage.df[,"value"], mask= MaskImage.df[,"value"]))
				
			# compute intensity per colony
			# this involve removing 0 values because they shouldn't be part of the colony intensity
			# the lapply below compute mean and median intensity for each colony - using horizontal compute (across columns first then across rows)
		
			image.name = haloImages[i]
			
			colony.intensity = do.call(rbind, pblapply(1:nrow(colony.px.index), function(j){
				
				left.col.index = colony.px.index[j,1]
				right.col.index = colony.px.index[j,1] + width.col
				top.row.index = colony.px.index[j,2]
				bottom.row.index = colony.px.index[j,2] + width.col
			
				colony.grid = Halo.Mask.dt[x %in% left.col.index:right.col.index & y %in% top.row.index:bottom.row.index]
								
				x.expand = c(colony.grid[mask == 1][,x]-ring.px, colony.grid[mask == 1][,x]+ring.px)
				y.expand = c(colony.grid[mask == 1][,y]-ring.px, colony.grid[mask == 1][,y]+ring.px)

				backd = subset(colony.grid,	x %in% left.col.index:(left.col.index + band.px) & 
											y %in% top.row.index:(top.row.index + band.px) |
											x %in% (right.col.index - band.px):right.col.index & 
											y %in% (bottom.row.index - band.px):bottom.row.index)
											
				perim = subset(colony.grid, x %in% x.expand & y %in% y.expand & mask == 0)

				grid.value = colony.grid[,median(halo.value)] 
				backgd.value = backd[,median(value)]
				colony.value = colony.grid[,median(halo.value[mask == 1])]
				perim.value = perim[,median(halo.value)]
				
				norm.colony.value = colony.value - backgd.value
				norm.perim.value = perim.value - backgd.value
							
				as.data.frame(cbind(image.name, grid.value, backgd.value, colony.value, perim.value, norm.perim.value, norm.colony.value))
	
					}))

			}))

setwd(analysis.folder)

Halo.data[,"grid.value"] = as.numeric(Halo.data[,"grid.value"])
Halo.data[,"backgd.value"] = as.numeric(Halo.data[,"backgd.value"])
Halo.data[,"colony.value"] = as.numeric(Halo.data[,"colony.value"])
Halo.data[,"perim.value"] = as.numeric(Halo.data[,"perim.value"])
Halo.data[,"norm.perim.value"] = as.numeric(Halo.data[,"norm.perim.value"])
Halo.data[,"norm.colony.value"] = as.numeric(Halo.data[,"norm.colony.value"])

Halo.data = cbind(image.name = Halo.data[,1], map.complete, Halo.data[,2:ncol(Halo.data)])





