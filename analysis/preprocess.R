preprocess = function(fname) {
	# read csv file
	raw = read.csv(file = fname, header = TRUE, sep = "\t")
	# add 3 extra columns
	raw$ratio1 = raw$similar_genes_count / raw$genes1_count
	raw$ratio2 = raw$similar_genes_count / raw$genes2_count
	raw$ratioa = (raw$ratio1 + raw$ratio2) / 2
	# create single-value cluster ID columns
	raw$clusterid1 = paste0(raw$genome1_ID, raw$cluster1)
	raw$clusterid2 = paste0(raw$genome2_ID, raw$cluster2)
	# sort: DESC by ratioa, ratio1, ratio2, size1_kb, size2_kb,
	# avg_protein_identity, similar_genes_count, P, K, S
	raw = raw[with(raw, order(-ratioa, -ratio1, -ratio2, -size1_kb, -size2_kb,
							  -avg_protein_identity, -similar_genes_count,
							  -P, -K, -S)), ]
	# remove numeric row names; may later replace with some IDs?
	row.names(raw) = NULL
	return(raw)
}