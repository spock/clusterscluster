genomes_to_groups = function(clusters2gids, all, uniq = TRUE) {
	# clusters2gids: a list, mapping each cluster ID
	# (concatenation of genome ID and cluster number)
	# to an arbitrary integer gid (group ID),
	# e.g. as returned by group_clusters(data, cutoff);
	# all: all clusters dataframe (e.g. all60_nodup_noknown,
	# with id/genome_ID/cluster columns), used to extract genome-cluster relations.
	# uniq: should repeated GIDs be removed (default:yes).
	# Return a list mapping each genome ID to a vector of integer group IDs (gids).
	genomes2gids = list() # final result
	genomes = unique(all$genome_ID) # all 203 genome IDs
	total_clusters = 0
	for (genome in genomes) {
		#cat('genome:', genome, "\n")
		groups = integer() # empty vector to store gids for this genome
		# vector of cluster numbers of this genome
		clusters = all[which(all$genome_ID == genome), 'cluster']
		total_clusters = total_clusters + length(clusters)
		for (cluster in clusters) {
			clusterid = paste0(genome, cluster)
			#cat('clusterid:', clusterid, "\n")
			#cat('GID for the cluster:', clusters2gids[[clusterid]], "\n")
			groups = c(groups, clusters2gids[[clusterid]])
		}
		# Save the resulting vector.
		#cat("groups:\n", groups, "\n");
		#cat("unique groups:\n", unique(groups), "\n");
		#cat('genomes2gids[genome] for', genome, 'is', genomes2gids[[genome]], "\n")
		if (uniq) {
			genomes2gids[[genome]] = unique(groups)
		}
		else {
			genomes2gids[[genome]] = groups
		}
	}
	cat('Processed', length(genomes), 'genomes and', total_clusters, "clusters.\n")
	return(genomes2gids)
}

groups_to_genomes = function(genomes2gids) {
	# Given genomes2gids, a list mapping every genome ID to an integer vector
	# of group IDs (gids),
	# return inverse mapping of gids to character vectors of genome IDs in the group.
	gids2genomes = list() # final result
	# Iterate over the list, populate the final result.
	for (genome in names(genomes2gids)) {
		# integers can't be list keys?..
		gids = as.character(genomes2gids[[genome]])
		# gids is a vector of IDs, iterate them all
		for (gid in gids) {
			if (gid %in% names(gids2genomes)) {
				#cat("\tgid", gid, "already exists, appending\n")
				gids2genomes[[gid]] = c(gids2genomes[[gid]], genome)
			}
			else {
				gids2genomes[[gid]] = genome
			}
		}
	}
	cat('Mapped', length(gids2genomes), "gids.\n")
	return(gids2genomes)
}

count_groups_per_genome = function(genomes2gids) {
	# Helper function: count how many groups of clusters each genome has.
	as.numeric(sapply(genomes2gids, length))
}

count_genomes_per_group = function(gids2genomes) {
	as.numeric(sapply(gids2genomes, length))
}
