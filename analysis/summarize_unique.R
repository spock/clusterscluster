summarize_unique = function(df) {
	rows = nrow(df)
	final = list()
	for (i in 1:rows) {
		final = c(final, paste0(df[i, 'genome1_ID'], df[i, 'cluster1']))
		final = c(final, paste0(df[i, 'genome2_ID'], df[i, 'cluster2']))
		#paste(df[i, c('genome1_ID', 'genome2_ID')], df[i, c('cluster1', 'cluster2')])
	}
	return(unique(final))
}