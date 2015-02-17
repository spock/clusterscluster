total_clusters = function(fname) {
	# get directory name where to run the script
	dname = dirname(file.path(fname))
	allname = file.path(dname, 'all_clusters.csv')
	# read in and return the df from all_clusters.csv file
	all_clusters = read.csv(file = allname, header = TRUE, sep = "\t")
	return(all_clusters)
}