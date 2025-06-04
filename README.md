Simple tangle traverser, generates "not terrible" traversal through repetitive genomic tangle that somehow matches coverage and the read alignment
For help run ./tangle_traverer.py --help

Requires pulp, ahocorasick, networkx, statistics, logging   python libraries.
## Usage

```bash
./tangle_traverser.py --gfa <gfa_file> --alignment <alignment_file> --output <output_directory> [options]
```

### Required Arguments:
- `--gfa`: Path to the GFA file with the graph structure
- `--alignment`: Path to a file with GraphAligner alignment
- `--output`: Output directory for all result files (will be created if it doesn't exist)
OR
- `--verkko-output` - HiFi graph and ONT alignments from verkko would be used.
- Tangle should be specified with either one internal node (`--tangle-node utig4-267`) or a file with complete list of internal tangle nodes one by line (`--tangle-file nodes.list`)
### Example:
```bash
./tangle_traverser.py --gfa assembly.gfa --alignment reads.gaf --output results_dir --tangle-node utig4-267 --quality-threshold 20
```

UNDER CONSTRUCTION
