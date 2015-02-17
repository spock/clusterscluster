make_pan.matrix = function(g2gids) {
	# g2gids: genomes2gids_full, a list mapping each genome ID to a vector of
	# group IDs (gids), *including* repeated groups, e.g.
	# genomes2gids_05_strepto_full = genomes_to_groups(clusters2gids_05_strepto,
	#													 all_strepto_clusters,
	#													 uniq = FALSE)
	# Return pan.matrix, usable by micropan.
	# Number of genomes/rows in the matrix.
	rows = length(g2gids)
	# Total number of unique groups/columns in the matrix.
	cols = length(count_genomes_per_group(groups_to_genomes(g2gids)))
	# Initialize a matrix.
	pan_matrix = matrix(0, nrow=rows, ncol=cols,
						dimnames=list(names(g2gids), # row names
									  as.character(1:cols))) # col names
	# Iterate g2gids, populate matrix.
	for (genome in names(g2gids)) {
		#cat('genome is', genome, "\n")
		gids = as.character(g2gids[[genome]])
		# gids is a vector of IDs, iterate them all
		for (gid in gids) {
			#cat('gid is', gid, "\n")
			#cat('cell value before is', pan_matrix[genome, gid], "\n")
			pan_matrix[genome, gid] = pan_matrix[genome, gid] + 1
			#cat('cell value after is', pan_matrix[genome, gid], "\n")
		}
	}
	return(pan_matrix)
}

pan_matrix_05_strepto = make_pan.matrix(genomes2gids_05_strepto_full)
pan_matrix_05 = make_pan.matrix(genomes2gids_05_full)