# copy-based on lookback_unique.R; does not produce plot(s), instead
# generates a dataframe with all the permutations, with 1 column per genome,
# and 1 row per permutation (e.g. 100x203 DF for 100 permutations and 203 genomes).

lookback_unique_data = function(data, all, threshold = 0.5, n.perm = 100) {
	# data is e.g. data60_nodup_noknown.
	# all is e.g. all60_nodup_noknown
	# n.perm is the number of genome order permutations to perform.
	# Examine genomes one by one.
	# The first genome would have 100% unique clusters, which we add to a variable.
	# The second may have some non-unique, so we only add new clusters to the variable.
	# At the end we have a vector along the number of genomes (e.g. 203), where
	# each cell has a cumulative sum of novel clusters from the current and
	# all the previous genomes.

	genomes = unique(as.character(all$genome_ID))
	# Prepare 'data' subset with intra = 0 and aratio >= threshold
	d = data[which(data$is_intra == 0 & data$ratioa >= threshold), ]
	# Factors and their levels interfere with further processing: convert to simple vectors.
	d$genome1_ID = as.character(d$genome1_ID)
	d$genome2_ID = as.character(d$genome2_ID)
	d$species1 = as.character(d$species1)
	d$species2 = as.character(d$species2)
	d$type1 = as.character(d$type1)
	d$type2 = as.character(d$type2)
	cat('Similar clusters dataframe has', nrow(d), "rows.\n")
	#cat ('xlim', length(genomes), 'type', typeof(length(genomes)), "\n")
	#cat ('ylim', nrow(all), 'type', typeof(nrow(all)), "\n")

	# Build column names for the dataframe and each row.
	# Empty column names vector.
	col_names = character()
	for (genome_counter in 1:length(genomes)) {
		col_names = c(col_names, paste0('G', genome_counter))
	}

	# dataframe to hold all the data
	permuts = data.frame(matrix(NA, nrow = n.perm, ncol = length(genomes)),
						 row.names = 1:n.perm)
	# There was no real benefit in colnames(permuts) = col_names
	colnames(permuts) = 1:length(genomes)

	for (iteration in 1:n.perm) {
		# Prepare the final row for this iteration.
		row = integer(length = length(genomes))
		names(row) = col_names

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
				# the first element, also paste0('G', genome_counter)
				row[genome_counter] = novel
			}
			else {
				# following elements, sum of previous + current novel,
				# also paste0('G', genome_counter - 1)
				row[genome_counter] = row[(genome_counter - 1)] + novel
			}

			# update 'seen' to include current 'clusters'
			merged_seen = c(seen_clusters, clusters)
			seen = d[(d$clusterid1 %in% merged_seen
					  | d$clusterid2 %in% merged_seen),]
		}
		# collect data for later plotting
		permuts[iteration, ] = row
	}
	return(permuts)
}

#lookback_unique_data(data60_nodup_noknown, all60_nodup_noknown) # threshold 0.5, 100 permutations
#lookback_unique_data(data60_nodup_noknown, all60_nodup_noknown, 0.8, n.perm = 100)
#lookback_unique_data(data60_nodup_noknown, all60_nodup_noknown, 0.2, n.perm = 100)
