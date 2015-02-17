create_types_list = function(all_name, nonunique_list, unique_list) {
  # all_name: df of all clusters, should have 'id' column
  # nonunique_list: a get-nonunique list, namely the $list component of it
  # unique_list: complementary to the above

  unique_types = data.frame()
  nonunique_types = data.frame()

  for (id in all_name$id) {
  	if (id %in% nonunique_list) {
  		nonunique_types[id, 'type'] = all_name[which(all_name$id == id), 'type']
  	}
  	else if (id %in% unique_list) {
  		unique_types[id, 'type'] = all_name[which(all_name$id == id), 'type']
  	}
  	else {
  		cat('Neither unique nor non-unqiue:', id, '\n')
  	}
  }

  unique_table = as.numeric(table(unique_types))
  names(unique_table) = names(table(unique_types))

  nonunique_table = as.numeric(table(nonunique_types))
  names(nonunique_table) = names(table(nonunique_types))

  result = list(unique_types=unique_types, nonunique_types=nonunique_types,
                unique_table=unique_table, nonunique_table=nonunique_table)
  return(result)
}