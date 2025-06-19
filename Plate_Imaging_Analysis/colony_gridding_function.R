	setwd(folder1)
	
	qc.image = dir(folder1)[1]

	white.grid = as.data.frame(grayscale(load.image(qc.image)))
	colnames(white.grid) = c("x", "y", "value")
	white.grid[,"value"] = 0
	
	white.grid[white.grid $x %in% colony.px.index[,1] ,"value"] = 1
	white.grid[white.grid $y %in% colony.px.index[,2],"value"] = 1
	
	setwd(analysis.folder)

	save.image(plot(add(list(grayscale(load.image(paste0(folder1, "/", qc.image))), as.cimg(white.grid)))), "gridQC.jpg")
	dev.off()