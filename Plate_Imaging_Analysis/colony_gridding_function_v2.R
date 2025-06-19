colonyImages = dir(folder1)[grep('.JPG$', dir(folder1))]
colonyImages = colonyImages[-grep('gridded', colonyImages)]	


lapply(1:i, function(i){
	
	setwd(folder1)
	qc.image = paste0(folder1, "/", colonyImages[i])

	white.grid = as.data.frame(grayscale(load.image(qc.image)))
	colnames(white.grid) = c("x", "y", "value")
	white.grid[,"value"] = 0
	
	white.grid[white.grid $x %in% colony.px.index[,1] ,"value"] = 1
	white.grid[white.grid $y %in% colony.px.index[,2],"value"] = 1
	
	setwd(analysis.folder)
	
	dir.create("gridQC")
	
	setwd(paste0(analysis.folder, "/gridQC"))
	
	save.image(plot(add(list(grayscale(load.image(qc.image)), as.cimg(white.grid)))), paste0("gridQC_", colonyImages[i], ".jpg"))
	
	dev.off()
	
})
