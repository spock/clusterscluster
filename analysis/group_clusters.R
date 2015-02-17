group_clusters = function (data, all, cutoff) {
	# Given 'data' like data60_nodup_noknown, cutoff by the 'ratioa' column,
	# and 'all' like all60_nodup_noknown (including clusters
	# NOT IN 'data' at all, with columns id, genome_ID, cluster and type),
	# return a list mapping cluster IDs to the assigned groups.
	# Creates singletons (groups with only one cluster).
	# Somewhat duplicates what the get_nonunique() already does...

	clusters2gids = list()
	subset_cutoff = data[which(data$ratioa >= cutoff), ]
	# Order descending by clusterid1/custerid2: necessary for properly avoiding
	# full graph traversal and still finding all the clusters belonging to a
	# single group.
	subset_cutoff = subset_cutoff[order(subset_cutoff$clusterid1,
										subset_cutoff$clusterid2), ]

	# We must create separate (singleton) groups for "unique" clusters.
	cat('There are', nrow(all), "total clusters.\n")
	# "non-unique", or similar, clusters, are all which have links with ratioa > 0.5
	similar = unique(c(subset_cutoff[, 'clusterid1'], subset_cutoff[, 'clusterid2']))
	cat(length(similar), "of them are similar.\n")
	# then "unique" are all the remaining clusters
	singletons = setdiff(all$id, similar)
	cat(length(singletons), "of them are unique.\n")

	# counter for group ID (incremented before use)
	gid = 0
	for (i in 1:nrow(subset_cutoff)) {
		# Check if cluster1 or cluster2 are already in clusters2gids:
		first  = data[i, 'clusterid1'] %in% names(clusters2gids)
		second = data[i, 'clusterid2'] %in% names(clusters2gids)
		if (first && second) { # Both are in groups already, do nothing.
			next
		}
		# if yes for any single of them: add the 2nd to the same group
		if (first) {
			clusters2gids[data[i, 'clusterid2']] = clusters2gids[[data[i, 'clusterid1']]]
			next
		}
		if (second) {
			clusters2gids[data[i, 'clusterid1']] = clusters2gids[[data[i, 'clusterid2']]]
			next
		}
		# Only get here if none of the clusters is assigned to a group.
		# Simply create a new group.
		# Increment group ID.
		gid = gid + 1
		clusters2gids[data[i, 'clusterid1']] = gid
		clusters2gids[data[i, 'clusterid2']] = gid
	}
	cat('Created', gid, "groups for similar clusters.\n")

	for (clusterid in singletons) {
		gid = gid + 1
		clusters2gids[clusterid] = gid
	}
	cat('Created a total of', gid, "groups (similar+singletons) of clusters.\n")
	return(clusters2gids)
}

groups_to_clusters = function(clusters2gids) {
	# Given clusters2gids, a list mapping every cluster ID to an integer group ID (gid),
	# return inverse mapping of gids to character vectors of cluster IDs in the group.
	gids2clusters = list() # final result
	# Iterate over the list, populate the final result.
	for (clusterid in names(clusters2gids)) {
		gid = as.character(clusters2gids[[clusterid]])
		#cat("Found group", gid, "\n")
		#lapply(list, function(.element) {...process each element of the list....})
		if (gid %in% names(gids2clusters)) {
			#cat("\tgid", gid, "already exists, appending\n")
			gids2clusters[[gid]] = c(gids2clusters[[gid]], clusterid)
		}
		else {
			gids2clusters[[gid]] = clusterid
		}
	}
	return(gids2clusters)
}

count_clusters_per_group = function(gids2clusters) {
	# Given gids2clusters, a list mapping group IDs (gids) to character vectors of
	# cluster IDs, return a simple vector of the number of clusters per each group.
	# Length of the vector equals the number of groups,
	# sum of the vector elements equals total number of clusters.
	as.numeric(sapply(gids2clusters, length))
}
