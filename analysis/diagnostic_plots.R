# Sourcing this script creates 4 overlayed plots relating different thresholds
# to the percentage of unique clustes obtained with those thresholds.
# Data is expected in the data60_nodup_noknown dataframe,
# complete list of clusters in all60_nodup_noknown.

diagnostic_plots = function(point_estimate = FALSE) {
	# point_estimate: demand % unique estimate with that specific threshold value.
	# more specifically, this means looking for (0, 1; 1, 2; ...) threshold ranges.

	# Source required functions.
	source('get_nonunique.R')

	# This was 3247 as of 2014-10-22.
	total_clusters = nrow(all60_nodup_noknown)

	# All resulting values (effectors and uniques) must be in [0, 1] range.

	# Average identity threshold effect on %% unique.
	ident_unique = data.frame(identity_over = numeric(0), percent_unique = numeric(0))
	if (point_estimate) {
		iter_range = 0:99
	}
	else {
		iter_range = 0:100
	}
	for (ident in iter_range) {
		cat('Identity', ident, '\n')
		filtered = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								 ratioc = 0, identc = ident, genec = 0,
								 pabs = 0, spot = point_estimate)
		percentage = (total_clusters - filtered$count) / total_clusters
		newrow = c(ident/100, percentage)
		cat('Percent unique:', percentage, '\n')
		ident_unique = rbind(ident_unique, newrow)
	}
	plot(ident_unique, xlab='minimal average protein identity', ylab='%% unique',
		 col = 'blue', pch=1, ylim=c(0,1))

	# Synteny effect on %% unique, 3 graphs: Pearson, Spearman, Kendall
	synt_uniqueP = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	synt_uniqueS = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	synt_uniqueK = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	if (point_estimate) {
		iter_range = seq(0, 0.99, by = 0.01)
	}
	else {
		iter_range = seq(0, 1, by = 0.01)
	}
	for (synt in iter_range) {
		cat('Synteny', synt, '\n')
		# P
		filteredP = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								  ratioc = 0, identc = 0, genec = 0,
								  pabs = synt, spot = point_estimate)
		percentageP = (total_clusters - filteredP$count) / total_clusters
		newrowP = c(synt, percentageP)
		cat('Percent unique, P:', percentageP, '\n')
		synt_uniqueP = rbind(synt_uniqueP, newrowP)
		# S
		filteredS = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								  ratioc = 0, identc = 0, genec = 0,
								  sabs = synt, spot = point_estimate)
		percentageS = (total_clusters - filteredS$count) / total_clusters
		newrowS = c(synt, percentageS)
		cat('Percent unique, S:', percentageS, '\n')
		synt_uniqueS = rbind(synt_uniqueS, newrowS)
		# K
		filteredK = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								  ratioc = 0, identc = 0, genec = 0,
								  kabs = synt, spot = point_estimate)
		percentageK = (total_clusters - filteredK$count) / total_clusters
		newrowK = c(synt, percentageK)
		cat('Percent unique, K:', percentageK, '\n')
		synt_uniqueK = rbind(synt_uniqueK, newrowK)
	}
	#par(new = TRUE)
	plot(synt_uniqueP, xlab='minimal synteny, Pearson', ylab='%% unique',
		 col = 'green', pch = 2, ylim=c(0,1))
	plot(synt_uniqueS, xlab='minimal synteny, Spearman', ylab='%% unique',
		 col = 'green', pch = 2, ylim=c(0,1))
	plot(synt_uniqueK, xlab='minimal synteny, Kendall', ylab='%% unique',
		 col = 'green', pch = 2, ylim=c(0,1))

	# Score/ratioa effect on %% unique.
	ratioa_unique = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	if (point_estimate) {
		iter_range = seq(0, 0.99, by = 0.01)
	}
	else {
		iter_range = seq(0, 1, by = 0.01)
	}
	for (ratioa in iter_range) {
		cat('Score/ratioa', ratioa, '\n')
		filtered = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								 ratioc = ratioa, identc = 0, genec = 0,
								 pabs = 0, spot = point_estimate)
		percentage = (total_clusters - filtered$count) / total_clusters
		newrow = c(ratioa, percentage)
		cat('Percent unique:', percentage, '\n')
		ratioa_unique = rbind(ratioa_unique, newrow)
	}
	#par(new = TRUE)
	plot(ratioa_unique, xlab='minimal score', ylab='%% unique',
		 col = 'red', pch=3, ylim=c(0,1))

	# Absolute similar gene count effect on %% unique.
	similar_unique = data.frame(similar_over = integer(0), percent_unique = numeric(0))
	#maximal_similar = max(data60_nodup_noknown$similar_genes_count)
	maximal_similar = max(data60_nodup_noknown$similar_genes_count) / 2
	if (point_estimate) {
		iter_range = 1:(maximal_similar - 1)
	}
	else {
		iter_range = 1:maximal_similar
	}
	for (similar in iter_range) {
		cat('Similar genes', similar, '\n')
		filtered = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								 ratioc = 0, identc = 0, genec = similar,
								 pabs = 0, spot = point_estimate)
		percentage = (total_clusters - filtered$count) / total_clusters
		#newrow = c(similar/maximal_similar, percentage)
		newrow = c(similar, percentage)
		cat('Percent unique:', percentage, '\n')
		similar_unique = rbind(similar_unique, newrow)
	}
	cat('Maximal similar genes:', maximal_similar, '\n')
	#par(new = TRUE)
	plot(similar_unique, xlab='similar genes count', ylab='%% unique',
		 col = 'purple', pch=4, ylim=c(0,1))

 	# genes count ratio (genes1_count/genes2_count) effect on %% unique
	# this requires modification of get_nonunique and was not done; likely similar to ratioa

 	# plot all on the same canvas
 	#plot(ident_unique, synt_unique, ratioa_unique)
}

diagnostic_plots()
#diagnostic_plots(point_estimate = TRUE)
