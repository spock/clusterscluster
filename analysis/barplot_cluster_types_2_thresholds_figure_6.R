# for 2 thresholds (1 and 2), create a single barplot with 3 bars:
# All, Non-unique at threshold 1, Non-unique at threshold 2.

# (all60_nodup_noknown,
#  nonunique60_05$list,
#  unique60_05,
#  nonunique60_05_10_05$list,
#  unique60_05_10_05)
prepare_2_barplot_data = function(all, nonunilist1, uni1, nonunilist2, uni2) {
	# all: all60_nodup_noknown (with at least id and type columns)
	# nonunilist1,2: nonunique60_05_10_05$list <- list!
	# uni1,2: unique60_05_10_05.
	# Return 3xcluster_types matrix.

	# make create_types_list() available
	source('create_dfs_with_cluster_types.R')
	# create_types_list returns a list with 2 components,
	# $unique_table and $nonunique_table, which summarize cluster counts per type
	clustertypes1 = create_types_list(all, nonunilist1, uni1)
	clustertypes2 = create_types_list(all, nonunilist2, uni2)

	# Cluster types I want to have in the graph.
	known_clusters = c('Bacteriocin', 'Butyrolactone', 'Composite', 'Ectoine',
					   'Lantipeptide', 'Melanin', 'NRPS', 'Other', 'Siderophore',
					   'T1PKS', 'T2PKS', 'T3PKS', 'Terpene')
	# Original letters case, except for COMPOSITE, which is preserved to maintain same index.
	known_small = c('bacteriocin', 'butyrolactone', 'COMPOSITE', 'ectoine',
					'lantipeptide', 'melanin', 'nrps', 'other', 'siderophore',
					't1pks', 't2pks', 't3pks', 'terpene')

	# Initialize the final matrix.
	types = matrix(0, nrow = 3, ncol = length(known_clusters))
	rownames(types) = c('All', 'Threshold 1', 'Threshold 2')
	colnames(types) = known_clusters

	# Assemble 'clustertypes' object with all the data
	clustertypes = list()
	# First, add all_table to clustertypes.
	clustertypes$all_table = table(all$type)
	clustertypes$threshold1 = clustertypes1[['nonunique_table']]
	clustertypes$threshold2 = clustertypes2[['nonunique_table']]

	# Iterate all items of clustertypes, assigning counts to where they
	# belong in the initialized empty matrix.
	for (iter_table in c('all_table', 'threshold1', 'threshold2')) {
		if (iter_table == 'all_table') {
			current_row = 'All'
		}
		else if (iter_table == 'threshold1') {
			current_row = 'Threshold 1'
		}
		else {
			current_row = 'Threshold 2'
		}
		for (current_cluster in names(clustertypes[[iter_table]])) {
			# Is this cluster in known_small?
			if (current_cluster %in% known_small) {
				# Yes! Assign count to the appropriate matrix cell, simply by the index.
				proper_name = known_clusters[which(known_small == current_cluster)]
				#cat('Cluster', current_cluster, 'is known, proper name is', proper_name, '.\n')
				types[current_row, proper_name] = types[current_row, proper_name] + clustertypes[[iter_table]][[current_cluster]]
			}
			# Maybe it has a dash in the name, making it composite?
			else if ( grepl('-', current_cluster, fixed = T) ) {
				# There is at least 1 dash in the cluster name, assign to Composite.
				types[current_row, 'Composite'] = types[current_row, 'Composite'] + clustertypes[[iter_table]][[current_cluster]]
				#cat('Cluster', current_cluster, 'is Composite.\n')
			}
			# Must be something Other...
			else {
				types[current_row, 'Other'] = types[current_row, 'Other'] + clustertypes[[iter_table]][[current_cluster]]
				#cat('Cluster', current_cluster, 'is Other.\n')
			}
		}
	}

	print(types)
	return(types)
}

# make the barplot function available
source('barplot_cluster_types_2.R')

# nonunique60_05$ list, df: filtered by score 0.5, threshold 1
# nonunique60_05_10_05: threshold 2; types60_05_10_05 = create_types_list(all60_nodup_noknown, nonunique60_05_10_05$list, unique60_05_10_05) - cluster types for the 2nd threshold
all_and_2_thresholds = prepare_2_barplot_data(all60_nodup_noknown,
											  nonunique60_05$list,
											  unique60_05,
											  nonunique60_05_10_05$list,
											  unique60_05_10_05)

#postscript("Figure 6.eps", width = 7, height = 5.25,
#		   horizontal = FALSE, onefile = FALSE)
#pdf("Figure 6.pdf", width = 7, height = 5.25, useDingbats = FALSE)
ppi = 300
png("Figure 6.png", width = 7*ppi, height = 5.25*ppi, res = ppi)
barplot_cluster_types_2(all_and_2_thresholds,
						title = 'Cluster types, all and non-unique for 2 thresholds')
dev.off()
