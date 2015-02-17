# Prepare data for ggiNEXT():
# - two thresholds, 0.2 and 0.5;
# - two data types, abundance and incidence.

rarefaction = function() {
	# Threshold 0.5.
	clusters2gids_05 = group_clusters(data60_nodup_noknown,
									  all60_nodup_noknown, 0.5)
	gids2clusters_05 = groups_to_clusters(clusters2gids_05)
	clusters_per_group_05 = count_clusters_per_group(gids2clusters_05)

	genomes2gids_05 = genomes_to_groups(clusters2gids_05, all60_nodup_noknown)
	gids2genomes_05 = groups_to_genomes(genomes2gids_05)
	genomes_per_group_05 = count_genomes_per_group(gids2genomes_05)

	# Threshold 0.2.
	clusters2gids_02 = group_clusters(data60_nodup_noknown,
									  all60_nodup_noknown, 0.2)
	gids2clusters_02 = groups_to_clusters(clusters2gids_02)
	clusters_per_group_02 = count_clusters_per_group(gids2clusters_02)

	genomes2gids_02 = genomes_to_groups(clusters2gids_02, all60_nodup_noknown)
	gids2genomes_02 = groups_to_genomes(genomes2gids_02)
	genomes_per_group_02 = count_genomes_per_group(gids2genomes_02)

	# Threshold 0.5, Streptomyces and Kutzneria only.
	# We need to subset data60_nodup_noknown to only have Streptomyces rows.
	strepto60_nodup_noknown = data60_nodup_noknown[which(data60_nodup_noknown$genome1_ID %in% streptomyces_genomes$genomeID & data60_nodup_noknown$genome2_ID %in% streptomyces_genomes$genomeID), ]

	# we need a different list of all clusters here, build it first.
	all_strepto_clusters = all60_nodup_noknown[all60_nodup_noknown$genome_ID %in% streptomyces_genomes$genomeID, ]

	clusters2gids_05_strepto = group_clusters(strepto60_nodup_noknown,
									  		  all_strepto_clusters, 0.5)
	gids2clusters_05_strepto = groups_to_clusters(clusters2gids_05_strepto)
	clusters_per_group_05_strepto = count_clusters_per_group(gids2clusters_05_strepto)
	genomes2gids_05_strepto = genomes_to_groups(clusters2gids_05_strepto, all_strepto_clusters)
	gids2genomes_05_strepto = groups_to_genomes(genomes2gids_05_strepto)
	genomes_per_group_05_strepto = count_genomes_per_group(gids2genomes_05_strepto)
}

rarefaction()