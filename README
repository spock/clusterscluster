Code contributions (patches) of all types are welcome.

Project status
==============
Frozen due to lack of human resources for further development.
Cooperation/collaboration may unfreeze the project, send your inquiries.

Purpose of the software
=======================
Annotate secondary metabolite biosynthesis clusters in supplied genomes with antismash2,
then compare clusters to find similar and unique ones.  
Analyze cluster pairs table, find genomes with the highest absolute unique clusters count.  
Examine changes to unique clusters count when changing thresholds.  
Plot clusters by type. Plot cumulative unique clusters scatterplot, using genome order resampling.  
Perform rarefaction analysis.  
For a fairly detailed example, please see the [Automating Assessment of the Undiscovered Biosynthetic Potential of Actinobacteria](http://biorxiv.org/content/early/2016/01/07/036087) preprint.

Dependencies
============
* Python 2.7+
    * SciPy
* antismash2 (possibly but not necessarily patched, see below)
    * antismash3 should be fully compatible, but wasn't tested yet
* usearch from http://www.drive5.com/usearch/
* InParanoid 4.1 from http://software.sbc.su.se/cgi-bin/request.cgi?project=inparanoid :
    * free for private study, education or nonprofit research only
    * Perl
        * XML::Parser
    * BLAST 2.2.x
* QuickParanoid from http://pl.postech.ac.kr/QuickParanoid/ (patched, included)
* R: per-analysis packages
    * iNEXT
    * micropan

Installation
============

antismash2 + optional patch
---------------------------
After you install antismash2, run_antismash.py must be added to system's PATH.

cluster.py can use a version of antismash which avoids extending clusters using  
the rules in cluster_rules.txt. Here's how you can add --no-extensions option to antismash2:

- install antismash2 from git to a location you can write to, add run_antismash
  to PATH;
- apply the supplied 0001-added-support-for-no-extensions-option.patch using
  `git am` to add support for the --no-extensions option (patch was made against
  commit 04532e8aa05ff1039acba5bd11613e929c409eca)
- resolve any possible conflicts while patching.

InParanoid 4.1 + patch
----------------------
InParanoid 4.1 has a restrictive license, so it is not present in the repository.  
It should be installed into custom_inparanoid directory.  
There is a patch custom_inparanoid/inparanoid.pl.patch with some improvements, clean ups,  
and changes necessary for InParanoid to be used by cluster.py.

QuickParanoid (patched, included)
---------------------------------
QuickParanoid (already patched) is included in the quickparanoid directory of the repository.  
Patch against the original downloaded version is also available as quickparanoid/qp.patch.

Dependencies
------------
It should be possible to implement --skip-antismash and feed already-processed genomes into cluster.py.
Otherwise, antismash2 is required.  
SciPy.stats dependency is used for gene order correlation calculations.
This dependency can be removed as soon as an alternative (more sensitive) gene synteny method is implemented.  
InParanoid and QuickParanoid are only necessary if you intend using genome-wide orthology calculations
to estimate the number of orthologous links between clusters. Use --skip-orthology to skip both.

Usage
=====
$ cluster.py --help
usage: cluster.py [-h] [-d] [-q] [-V] [--trim] [--fulldp] [--highmem]
                  [--cutoff CUTOFF] [--skip-putative] [--skip-orthology]
                  [--no-name-problems] [--no-tree-problems]
                  [--emulate-inparanoid] [--prefix PREFIX] [--project PROJECT]
                  [--force] [--no-extensions] [--threshold THRESHOLD]
                  [--from-file FROM_FILE]
                  [path [path ...]]

positional arguments:
  path                  paths to GenBank files with genomes to analyze

optional arguments:
  -h, --help            show this help message and exit
  -d, --debug           set verbosity level to debug [default: False]
  -q, --quiet           report only warnings and errors [default: False]
  -V, --version         show program's version number and exit
  --trim                trim away antismash2 cluster extensions [default:
                        False]
  --fulldp              use full dynamic programming solution in usearch
                        alignment (much slower!) [default: False]
  --highmem             assume huge RAM: all GenBanks are loaded early and
                        kept in RAM (faster processing) [default: False]
  --cutoff CUTOFF       protein identity cut-off when aligning with usearch
                        [default: 0.6]
  --skip-putative       exclude putative clusters from the analysis [default:
                        False]
  --skip-orthology      do not run any orthology analysis [default: False]
  --no-name-problems    only use ortho-clusters which do not have diff.names
                        tree_conflict problems [default: False]
  --no-tree-problems    only use ortho-clusters which do not have
                        [diff.names/diff.numbers] tree_conflict problems
                        [default: False]
  --emulate-inparanoid  only print generated inparanoid commands, do not run;
                        exit after inparanoid blasting; suppress some normal
                        output [default: False]
  --prefix PREFIX       output CSV files prefix [default: out]
  --project PROJECT     put all the project files into this directory
                        [default: cluster_project]
  --force               insist on re-using existing project directory (this
                        will re-use existing intermediate files) [default:
                        False]
  --no-extensions       pass --no-extensions option to the modified antismash2
                        (see README for details) [default: False]
  --threshold THRESHOLD
                        cluster links with weight below this one will be
                        discarded [default: 0.0]
  --from-file FROM_FILE
                        read paths to GenBank files (one per line) from the
                        provided file

License
=======
This software is currently available under the GNU Affero General Public License Version 3 (AGPL-3.0).  
Other (commercial, exclusive) licenses are possible.

TODO list
=========
- check the source code for many more FIXME and TODO tags
- lib/gb2fasta.py is symlinked from another (open) repository; should use git submodules instead
- possible optimization for ortholgoy code: only run InParanoid on the translations of genes
  inside clusters, not entire genomes (i.e. "clustome vs clustome" instead of "genome vs genome")
- add mechanisms for discerning highly similar NRPK/PKS clusters,
  e.g. by the biosynthetic domains composition + domains order
- cleanup (remove or move into external files) obsolete code
- check/add file function explanations
- provide a sample analysis R session script/history file