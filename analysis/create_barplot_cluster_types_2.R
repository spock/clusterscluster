# plot 3 cluster types graphs for 3 thresholds

prepare_barplot_data = function(all, nonunilist, uni) {
	# all: all60_nodup_noknown (with at least id and type columns)
	# nonunilist: nonunique60_05_10_05$list <- list!
	# uni: unique60_05_10_05
	# return "types" matrix expected by barplot_cluster_types_2()
	# make create_types_list() available
	source('create_dfs_with_cluster_types.R')
	# create_types_list returns a list with 2 components,
	# $unique_table and $nonunique_table, which summarize cluster counts per type
	clustertypes = create_types_list(all, nonunilist, uni)

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
	rownames(types) = c('All', 'Unique', 'Non-unique')
	colnames(types) = known_clusters

	# Add all_table to clustertypes.
	clustertypes$all_table = table(all$type)

	# Iterate unique_table and nonunique_table, assigning counts to where they
	# belong in the initialized empty matrix.
	for (iter_table in c('all_table', 'unique_table', 'nonunique_table')) {
		if (iter_table == 'all_table') {
			current_row = 'All'
		}
		else if (iter_table == 'unique_table') {
			current_row = 'Unique'
		}
		else {
			current_row = 'Non-unique'
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

	# Count sums, for verification.
	total_unique = sum(types['Unique', ])
	total_nonunique = sum(types['Non-unique', ])
	total = total_unique + total_nonunique
	cat('Total unique:', total_unique, '\n')
	cat('Total non-unique:', total_nonunique, '\n')
	cat('Total:', total, '\n')
	cat('Total from the all-list:', length(all$id), '\n')

# 	# Convert all values into percentages.
# 	types['All', ] = round(100 * types['All', ] / sum(types['All', ]))
# 	types['Unique', ] = round(100 * types['Unique', ] / sum(types['Unique', ]))
# 	types['Non-unique', ] = round(100 * types['Non-unique', ] / sum(types['Non-unique', ]))
	print(types)
	return(types)
}


# make the barplot function available
source('barplot_cluster_types_2.R')


# input datasets

# nonunique60_05$ list, df: filtered by score 0.5
barplot_cluster_types_2(prepare_barplot_data(all60_nodup_noknown, nonunique60_05$list, unique60_05),
					  title = 'Cluster types, threshold 1')
# nonunique60_05_10_05: threshold 2; types60_05_10_05 = create_types_list(all60_nodup_noknown, nonunique60_05_10_05$list, unique60_05_10_05) - cluster types for the 2nd threshold
barplot_cluster_types_2(prepare_barplot_data(all60_nodup_noknown,
										   nonunique60_05_10_05$list,
										   unique60_05_10_05),
					  title = 'Cluster types, threshold 2')
# nonunique60_05_10_05_65_90: threshold 3
barplot_cluster_types_2(prepare_barplot_data(all60_nodup_noknown,
										   nonunique60_05_10_05_65_90$list,
										   unique60_05_10_05_65_90),
					  title = 'Cluster types, threshold 3')
