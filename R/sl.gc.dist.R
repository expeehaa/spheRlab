sl.gc.dist <-
function(lon,lat,Rsphere=1,sequential=TRUE) {
	
	if (sequential) {
	
		L = length(lon) - 1
		pointdist = rep(NA,L)
		for (l in 1:L) {
			alpha = lon[l] + 90
			beta = 90 - lat[l]
			gamma = 0
			rot.res = sl.rot(lon[l+1],lat[l+1],alpha,beta,gamma)
			pointdist[l] = (90 - rot.res$lat) * pi * Rsphere / 180
		}
		
	} else {
		
		L = length(lon) - 1
		alpha = lon[1] + 90
		beta = 90 - lat[1]
		gamma = 0
		rot.res = sl.rot(lon[2:(L+1)],lat[2:(L+1)],alpha,beta,gamma)
		pointdist = (90 - rot.res$lat) * pi * Rsphere / 180
		
	}
	
	return(pointdist)
	
}
