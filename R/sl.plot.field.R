sl.plot.field <-
function (plot.init.res,num,lon.v,lat.v,col.fill="colbar",col.border="colbar",colbar=sl.colbar.redgreyblue_256,colbar.breaks=NA,colbar.breaks.log=FALSE,border.lwd=0.01,border.lty=1) {
	
	Npoly = nrow(lon.v)
	colbar.res = NULL
	if (col.fill == "colbar" || col.border == "colbar") {
		colbar.res = sl.num2colbar(num,colbar,colbar.breaks,colbar.breaks.log)
		col.ind = colbar.res$colour.index
	}
	for (np in 1:Npoly) {
		cb.fill = col.fill
		cb.border = col.border
		if (col.fill == "colbar") {cb.fill = colbar[[col.ind[np]]]}
		if (col.border == "colbar") {cb.border = colbar[[col.ind[np]]]}
		slplotpolygon(plot.init.res,lon.v[np,],lat.v[np,],colfill=cb.fill,colborder=cb.border,borderlwd=border.lwd,borderlty=border.lty)
	}
	
	return(colbar.res)
}
