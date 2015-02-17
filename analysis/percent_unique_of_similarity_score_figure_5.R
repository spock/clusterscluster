# Data is expected in the data60_nodup_noknown dataframe,
# complete list of clusters in all60_nodup_noknown.
# Build a graph of %% unique clusters at different values of similarity
# score threshold. Condition is similarity_score >= threshold.

unique_similarity_plot = function() {
	# Source required functions.
	source('get_nonunique.R')

	# This was 3247 as of 2014-10-22.
	total_clusters = nrow(all60_nodup_noknown)

	# All resulting values (effectors and uniques) must be in [0, 1] range.

	# Score/ratioa effect on %% unique.
	ratioa_unique = data.frame(synt_over = numeric(0), percent_unique = numeric(0))
	iter_range = seq(0, 1, by = 0.01)
	for (ratioa in iter_range) {
		cat('Similarity score (ratioa)', ratioa, '\n')
		filtered = get_nonunique(data60_nodup_noknown, ignore_intra = T,
								 ratioc = ratioa, identc = 0, genec = 0,
								 pabs = 0)
		percentage = (total_clusters - filtered$count) / total_clusters * 100
		newrow = c(ratioa, percentage)
		cat('Percent unique:', percentage, '\n')
		ratioa_unique = rbind(ratioa_unique, newrow)
	}
	plot(ratioa_unique,
		 xlab='Similarity score threshold',
		 ylab='% unique clusters',
		 col = 'coral',
		 type = 'p', # s, p, l, b
		 lty = 6, # 1-6
		 pch = 1,
		 #ylim=c(0,1),
		 cex = 1, # points size
		 cex.axis = 1.3, lwd = 2.2, cex.lab = 1.3, las = 3)
	grid()
}



#postscript("Figure 5.eps", width = 7, height = 5.25,
#		   horizontal = FALSE, onefile = FALSE)
#pdf("Figure 5.pdf", width = 7, height = 5.25, useDingbats = FALSE)
ppi = 300
png("Figure 5.png", width = 7*ppi, height = 5.25*ppi, res = ppi)
unique_similarity_plot()
dev.off()
