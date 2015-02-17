lookback_unique = function(data, all, threshold = 0.5, n.perm = 100) {
	# data is e.g. data60_nodup_noknown.
	# all is e.g. all60_nodup_noknown
	# n.perm is the number of genome order permutations to perform.
	# Examine genomes one by one.
	# The first genome would have 100% unique clusters, which we add to a variable.
	# The second may have some non-unique, so we only add new clusters to the variable.
	# At the end we have a vector along the number of genomes (e.g. 203), where
	# each cell has a cumulative sum of novel clusters from the current and
	# all the previous genomes.
	# Plot n.perm lines on the same plot (or maybe make a scatterplot?)
	# to show the bend of this "cumulative novel" line.
	genomes = unique(as.character(all$genome_ID))
	# Prepare 'data' subset with intra = 0 and aratio >= threshold
	d = data[which(data$is_intra == 0 & data$ratioa >= threshold), ]
	# Factors and their levels interfere with further processing, must convert to simple vectors.
	d$genome1_ID = as.character(d$genome1_ID)
	d$genome2_ID = as.character(d$genome2_ID)
	d$species1 = as.character(d$species1)
	d$species2 = as.character(d$species2)
	d$type1 = as.character(d$type1)
	d$type2 = as.character(d$type2)
	cat('Similar clusters dataframe has', nrow(d), "rows.\n")
	#cat ('xlim', length(genomes), 'type', typeof(length(genomes)), "\n")
	#cat ('ylim', nrow(all), 'type', typeof(nrow(all)), "\n")
	# setup empty plot with axis
	plot(1,
		 type = 'n',
		 xlim = c(0, length(genomes)),
		 ylim = c(0, nrow(all)),
		 xlab = 'Genomes, cumulative',
		 ylab = 'Novel clusters, cumulative',
		 main = 'Cumulative increase of novel clusters with new genomes')

	for (iteration in 1:n.perm) {
		# Prepare the final, empty "plotting" matrix.
		results = matrix(0, nrow = 2, ncol = length(genomes))
		# name matrix rows, for cumulative novel and seen values
		rownames(results) = c('novel', 'seen')
		# Dataframe of the clusters we had already seen, as a subset of 'd'.
		seen = data.frame()
		# randomize genome order
		genomes = sample(genomes)
		for (genome_counter in 1:length(genomes)) {
			genome = genomes[genome_counter]
			clusters = all[which(all$genome_ID == genome), 'cluster']
			clusters = paste0(genome, clusters)
			#cat('genome', genome, 'has', length(clusters), "clusters\n")
			if (nrow(seen) > 0) {
				# character vector of the clusters we had seen before
				seen_clusters = unique(c(seen$clusterid1,
										 seen$clusterid2))
				# simply count how many we had seen
				seen = sum(clusters %in% seen_clusters)
			}
			else {
				#cat('this is one of the first genomes, so no seen clusters yet\n')
				seen_clusters = as.character()
				seen = 0
			}
			# count 'novel' clusters
			novel = length(clusters) - seen
			#cat("\tof these,", novel, 'are novel and', seen, "were seen before\n")
			# save the numbers to a matrix to plot later
			if (genome_counter == 1) {
				# the 1st cell, special case
				results['novel', genome_counter] = novel
				results['seen',  genome_counter] = seen
			}
			else {
				results['novel', genome_counter] = results['novel', (genome_counter - 1)] + novel
				results['seen',  genome_counter] = results['seen',  (genome_counter - 1)] + seen
			}
			# update 'seen' to include current 'clusters'
			merged_seen = c(seen_clusters, clusters)
			seen = d[(d$clusterid1 %in% merged_seen
					  | d$clusterid2 %in% merged_seen),]
			# likely a faster alternative, but somehow does not work - overwrites every time:
			#subset = d[which(d$clusterid1 %in% clusters
			#				 | d$clusterid2 %in% clusters),]
			#seen = rbind(seen, subset)
			#cat('After rbinding, seen has', nrow(seen), "rows.\n")
		}
		# now plot this single line, or collect data for later scatterplotting
		par(new = TRUE)
		plot(results['novel', ],
			 xlim = c(0, length(genomes)),
			 ylim = c(0, nrow(all)),
			 type = 'l',
			 axes = FALSE,
			 xlab = '', ylab = '')
	}
}

#lookback_unique(data60_nodup_noknown, all60_nodup_noknown, 0.5, n.perm = 500)
#lookback_unique(data60_nodup_noknown, all60_nodup_noknown, 0.5, n.perm = 100)
#lookback_unique(data60_nodup_noknown, all60_nodup_noknown, 0.8, n.perm = 100)
#lookback_unique(data60_nodup_noknown, all60_nodup_noknown, 0.2, n.perm = 100)