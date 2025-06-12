Simple tangle traverser, generates "not terrible" traversal through repetitive genomic tangle that somehow matches coverage and the read alignment.

For help run ./tangle_traverer.py --help

Requires pulp, ahocorasick, networkx, statistics, logging   python libraries.
## Usage

```bash
./tangle_traverser.py --gfa <gfa_file> --alignment <alignment_file> --output <output_directory> [options]
```

### Required Arguments:
- `--gfa`: Path to the GFA file with the graph structure
- `--alignment`: Path to a file with GraphAligner alignment

OR

- `--verkko-output` - HiFi graph ,coverage (ONT) and ONT alignments from verkko would be used.

- Tangle should be specified with either one internal node (`--tangle-node utig4-267`) or a file with complete list of internal tangle nodes one by line (`--tangle-file nodes.list`)
- `--output`: Output directory for all result files (will be created if it doesn't exist)

Be sure that you use the same graph for all the files (gfa, alignment, coverage, tangle and border nodes) -  HiFi graph (or --verkko-output) will not work with tangle nodes provided with respect to the final ONT resolved (utig4- in verkko case) graph.

### Example:
```bash
./tangle_traverser.py --gfa assembly.gfa --alignment reads.gaf --output results_dir --tangle-node utig4-267 --quality-threshold 20
```

### Verkko's final graph coverage fix
In verkko up to (and including )v2.2.1 coverage of the short nodes in tangles in final graph (assembly.homopolymer-compressed.gfa) is deeply flawed. To get the updated coverage file we suggest to run additional scripts

`./verkko_coverage_fix/utig4_to_utig1.py <assembly_folder> > utig42utig1.gaf`

`./verkko_coverage_fix/utig4_coverage_updater.py utig42utig1.gaf <assembly_folder>/assembly.homopolymer-compressed.noseq.gfa <assembly_folder>/2-processGraph/unitig-unrolled-hifi-resolved.ont-coverage.csv > utig4_upt.ont-coverage.csv`

and then pass `utig4_upt.ont-coverage.csv` as `--coverage-file` in main script.

Alternatively you can find how utig4- nodes match to the utig1- graph in utig42utig1.gaf and run tangle_traverser.py on the same tangle in hifi-only graph.


UNDER CONSTRUCTION
