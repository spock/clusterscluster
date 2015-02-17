# Sourcing this script creates a plot relating (ratio1 AND ratio2) joint threshold
# to the percentage of unique clustes obtained at different values of these two.
# Data is expected in the data60_nodup_noknown dataframe,
# complete list of clusters in all60_nodup_noknown.

diagnostic_ratio1_ratio2 = function() {
	# Source required functions.
	source('get_nonunique.R')

	# This was 3247 as of 2014-10-22.
	total_clusters = nrow(all60_nodup_noknown)

	# All resulting values (effectors and uniques) must be in [0, 1] range.

	# ratio1/ratio2 effect on %% unique.
	ratioa_unique = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	iter_range = seq(0, 1, by = 0.01)
	for (ratioa in iter_range) {
		cat('ratio1/ratio2', ratioa, '\n')
		filtered = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								 ratioc = ratioa, identc = 0, genec = 0,
								 pabs = 0, joint_identc = TRUE)
		percentage = (total_clusters - filtered$count) / total_clusters
		newrow = c(ratioa, percentage)
		cat('Percent unique:', percentage, '\n')
		ratioa_unique = rbind(ratioa_unique, newrow)
	}
	plot(ratioa_unique, xlab = 'minimal joint ratio1/ratio2 score',
		 ylab = '%% unique', col = 'coral', pch = 1, ylim = c(0,1))
}

diagnostic_ratio1_ratio2()
