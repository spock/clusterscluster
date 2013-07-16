# this is a heavily simplified short segment from antismash2/antismash/generic_modules/hmm_detection/__init__.py


def create_rules_dict():
    "Create a cluster rules dictionary from the cluster rules file"
    rulesdict = {}
    first = True
    for line in open("cluster_rules.txt", "r"):
        # skip the first line with the legend
        if first:
            first = False
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        key = parts.pop(0)
        rules = parts.pop(0)
        cutoff = int(parts.pop(0)) * 1000
        extension = int(parts.pop(0)) * 1000
        rulesdict[key] = (rules, cutoff, extension)
    return rulesdict


def get_extension_size_by_cluster_type(clustertype, rulesdict):
    extension = ''
#    print 'clustertype', clustertype
    if "-" in clustertype:
        extension = max([rulesdict[value][2] for value in clustertype.split("-")])
    else:
#        print 'rulesdict[clustertype]', rulesdict[clustertype]
        extension = rulesdict[clustertype][2]
    return extension
