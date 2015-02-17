per_genome = data.frame()
#per_genome = data.frame(row.names = all60_nodup_noknown$genome_ID)

rows = nrow(all60_nodup_noknown)
cat('processing', rows, "rows\n")

for (row in 1:rows) {
	#cat('row', row, "\n")
	genome = as.character(all60_nodup_noknown[row, 'genome_ID'])
	#cat('genome', genome, "\n")
	cid = as.character(all60_nodup_noknown[row, 'id'])
	#cat('cid', cid, "\n")
	# counter 1
	if (genome %in% rownames(per_genome)) {
		per_genome[genome, 'total'] = per_genome[genome, 'total'] + 1
	}
	else {
		per_genome[genome, 'total'] = 1
		per_genome[genome, 'nonunique'] = 0
		per_genome[genome, 'unique'] = 0
	}
	# counter 2
	if (cid %in% nonunique60_05$list) {
		# nonunique60_05 = get_nonunique(data60_nodup_noknown)
		per_genome[genome, 'nonunique'] = per_genome[genome, 'nonunique'] + 1
	}
	else if (cid %in% unique60_05) {
		# unique60_05 = get_unique(all60_nodup_noknown, nonunique60_05$list)
		per_genome[genome, 'unique'] = per_genome[genome, 'unique'] + 1
	}
	else {
		cat('Cluster', cid, 'neither unique nor similar to others?')
	}
}
# add the "unique percentage" column
per_genome$perc_uniq = per_genome$unique / per_genome$total
# finally, sort DESC by perc_uniq and uniq
per_genome = per_genome[with(per_genome, order(-unique, -perc_uniq, -nonunique)), ]
