get_nonunique = function(data, ignore_intra = T, ratioc = 0.5,
						 identc = 50, genec = 1,
						 # 3 measures of synteny: Pearson, Spearman, Kendall
						 pabs = 0, sabs = 0, kabs = 0,
						 spot = FALSE, joint_identc = FALSE) {
	# Print the count of non-unique clusters with given arguments.
	# ignore_intra: self-explanatory
	# ratioc: minimal allowed ratioa
	# identc: minimal allowed avg_protein_identity
	# genec: minimal allowed similar_genes_count
	# pabs: minimal allowed absolute value of P
	# spot: find links matching given argument + 1 / 0.01 (point estimate)
	# joint_identc: f TRUE, both ratio1 and ratio2 must be > identc; ratioa and spot are ignored
	# Return a list with 3 components:
	# - count, the count of non-unique clusters;
	# - list, the list of non-unique clusters, as pasted 'genomeIDcluster' strings
	# - df, filtered dataframe
	cat('starting with', nrow(data), "rows\n")
	if (ignore_intra) {
		# delete rows which have is_intra == 1
		# with this option, only inter-species clusters can be similar;
		# without it (ignore_intra=False), clusters within a single genome
		# may also be considered similar
		selector = -which(data$is_intra == 1)
		if (length(selector) > 0) {
			#cat("is_intra", selector, "\n")
			data = data[selector, ]
		}
	}
	cat('filtered is_intra:', nrow(data), "\n")
	# filter by ratioa[], ratio1, ratio2]
	if (joint_identc) {
		selector = -which(data$ratio1 < ratioc & data$ratio2 < ratioc)
	}
	else if (spot & ratioc > 0) {
		selector = -which(data$ratioa < ratioc | data$ratioa > ratioc + 0.01)
	}
	else {
		selector = -which(data$ratioa < ratioc)
	}
	if (length(selector) > 0) {
		#cat("ratioa", selector, "\n")
		data = data[selector, ]
	}
	cat('filtered by ratioa:', nrow(data), "\n")
	# filter by avg_protein_identity
	if (spot & identc > 0) {
		selector = -which(data$avg_protein_identity < identc | data$avg_protein_identity > identc + 1)
	}
	else {
		selector = -which(data$avg_protein_identity < identc)
	}
	if (length(selector) > 0) {
		#cat("avg_protein_identity", selector, "\n")
		data = data[selector, ]
	}
	cat('filtered by identity:', nrow(data), "\n")
	# filter by similar_genes_count
	if (spot & genec > 0) {
		selector = -which(data$similar_genes_count < genec | data$similar_genes_count > genec + 1)
	}
	else {
		selector = -which(data$similar_genes_count < genec)
	}
	if (length(selector) > 0) {
		#cat("similar_genes_count", selector, "\n")
		data = data[selector, ]
	}
	cat('filtered by similar genes:', nrow(data), "\n")
	# Synteny: Pearson
	if (spot & pabs > 0) {
		selector = -which(abs(data$P) < pabs | abs(data$P) > pabs + 0.01)
	}
	else {
		selector = -which(abs(data$P) < pabs)
	}
	if (length(selector) > 0) {
		#cat("minimal absolute P", selector, "\n")
		data = data[selector, ]
	}
	cat('filtered by minimal P:', nrow(data), "\n")
	# Synteny: Spearman
	if (spot & sabs > 0) {
		selector = -which(abs(data$S) < sabs | abs(data$P) > sabs + 0.01)
	}
	else {
		selector = -which(abs(data$S) < sabs)
	}
	if (length(selector) > 0) {
		data = data[selector, ]
	}
	cat('filtered by minimal S:', nrow(data), "\n")
	# Synteny: Kendall
	if (spot & kabs > 0) {
		selector = -which(abs(data$K) < kabs | abs(data$K) > kabs + 0.01)
	}
	else {
		selector = -which(abs(data$K) < kabs)
	}
	if (length(selector) > 0) {
		data = data[selector, ]
	}
	cat('filtered by minimal K:', nrow(data), "\n")
	# additional possible filters: log-ratio of size1_kb/size2_kb, ...
	# Summarize/count non-unique clusters.
	rows = nrow(data)
	# Code below is extremely slow.
# 	final = list()
# 	for (i in 1:rows) {
# 		final = c(final, data[i, 'clusterid1'])
# 		final = c(final, data[i, 'clusterid2'])
# 	}
# 	final = unique(final)
	final = unique(c(data[, 'clusterid1'], data[, 'clusterid2']))
	return(list(count=length(final), list=as.list(final), df=data))
}

get_unique = function(all_clusters, nonunique) {
	# Subtract nonunique clusters from all clusters to find unique clusters.
	# all_clusters: list of merged cluster IDs of all clusters (all_clusters$id)
	# nonunique: list of merged IDs of all clusters (count_nonunique$list)
	return(setdiff(all_clusters, nonunique))
}

ratio_unique = function(uniques, all) {
	return(length(uniques)/length(all))
}

ratio_nonunique = function(uniques, all) {
	return(1 - length(uniques)/length(all))
}