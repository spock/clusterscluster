# Note: edit the Y-lim before use.

# Create a 'beside' barplot from 'types' matrix with 1+ rows and 2+ columns.
# There are only 3 colours in the code.

barplot_cluster_types_2 = function(types, title = '') {
	# types: matrix with 1 row per data sequence
	# title: plot title ("main")
	
	op = par(mar = c(6.1, 4.1, 2, 0.1) + 0.1)
	mp = barplot(types, beside = T, density = c(25, 25), angle = c(45, 135),
				 col=c('aquamarine3', 'coral', 'lightcyan'), main = title,
				 ylab = 'Cluster count', #legend.text = TRUE
				 cex.axis = 1.3, lwd = 1.5, cex.lab = 1.3, las = 3,
				 ylim = c(0, 610)) # 30 for %%, 610 for absolute
	# Draw bar values above the bars
	text(mp + 0.4, types + 5, labels = format(types, digits = 1), pos = 3,
		 cex = .7, srt = 90)
	legend(x = "top", legend = rownames(types), #pch = 15,
		   density = c(25, 25), angle = c(45, 135), horiz = TRUE, cex = 1,
		   fill = c('aquamarine3', 'coral', 'lightcyan'),
		   bty = "n", inset = 0.025)
	par(op)
}
