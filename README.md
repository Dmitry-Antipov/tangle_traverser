TTT stands for Trivial Tangle Traverser. This tool generates "not terrible" traversals through repetitive genomic tangles that somehow matches coverage and the read alignment.

For help run ./TTT.py --help

Requires python &ge; 3.7 and dataclasses, pulp, ahocorasick, networkx, statistics, and logging python libraries.

<a href="https://docs.google.com/presentation/d/1TLjG0tR3o-Un8pnRawm0XEvTZY1hqpC4o6zhUthawAw/edit?usp=sharing">Slides explaining algorithmic details </a>

UNDER CONSTRUCTION!

## Usage

```bash
./TTT.py --gfa <gfa_file> --alignment <alignment_file> --output <output_directory> [options]
```
<details>
<summary><b>Will TTT help with this gap in my scaffold?</b></summary>

 
Generally there are three main reasons for gaps in a scaffold:
* Lack of coverage
  
  TTT searches for the "best" path in the assembly graph that traverses the gap. If there's no path because of the coverage gap --- nothing can be done.
  <p>
  <img width="400" height="400" alt="gap" src="https://github.com/user-attachments/assets/695419c0-a7fb-4728-8454-0433bfe66433" />  
    
  <em> Scaffold &lt;utig4-1497[N100000N:scaffold]&lt;utig4-340 --- nothing can be done </em>
</p>
  
* Long homozygous nodes
  
  Such gaps happen because of the read length being shorter than homozygous nodes. 
  Typical structure looks like a sequence of "bubbles" of similar length, interlaced with long homozygous nodes.
  TTT can be run on such tangles. But usually if those structures left unresolved in the assembly graph (especially if homozygous nodes are longer than ~100kbp homopolymer-compressed) then there's just no information in the read alignments helping to traverse this region, and thus it will be essentially a random guess.
  <p>
  <img width="400" height="400" alt="diploid_simple_tangle" src="https://github.com/user-attachments/assets/650052c6-5f53-43fa-bdfa-464c8a5d6fdb" />
    
  <em> Scaffolds &lt;utig4-1225&lt;utig4-1224[N5000N:ambig_bubble]&gt;utig4-1511&lt;utig4-1513 and &lt;utig4-1226&lt;utig4-1224[N5000N:ambig_bubble]&gt;utig4-1511&lt;utig4-1512. Because of long homozygous nodes utig4-1224 and utig4-1511 there's just no long reads connecting utig4-1228/utig4-1227 with utig4-1225/utig4-1226 or utig4-1512/utig4-1513 </em>
</p>

* Complex repeats
  
  TTT was designed for such cases. However, there are still limitations --- there can be no more than 2 haplotypes in the tangle (so rDNA tangles connecting multiple chromosomes are usually unresolvable), and for two tangle cases you should provide pairs of in- and out- nodes for each of the haplotypes.

  <p>
  <img width="400" height="400" alt="haploid tangle" src="https://github.com/user-attachments/assets/6df23394-811a-49a7-8606-993d2b8f1e89" />  

  <em>Gap caused by repeat array </em>
  </p>

  <p>
  <img width="400" height="400" alt="diploid tangle" src="https://github.com/user-attachments/assets/fc8418dd-1391-4032-ac62-0f9881c9a08c" />
  
  <em>Gap caused by large duplication of homozygous region, present inonly  one of the haplotypes</em>
  </p>

</details>

### Required Arguments:
- `--gfa`: Path to the GFA file with the graph structure
- `--alignment`: Path to a file with GraphAligner alignment

OR

- `--verkko-output` - HiFi graph ,coverage (ONT) and ONT alignments from verkko would be used.

- Tangle should be specified with either one internal node (`--tangle-node utig4-267`) or a file with complete list of internal tangle nodes one by line (`--tangle-file nodes.list`)
- `--output`: Output directory for all result files (will be created if it doesn't exist)

Be sure that you use the same graph for all the files (gfa, alignment, coverage, tangle and border nodes) -  HiFi graph (or --verkko-output) will not work with tangle nodes provided with respect to the final ONT resolved (utig4- in verkko case) graph.

Tangle traversal does no scaffolding. So, when running on tangles with two genomic paths you should provide incoming and outgoing boundary node pairs with --boundary-nodes boundary_file.tsv

File should be tab separated, in format:

`incoming_hap1_node  outgoing_hap1_node`

`incoming_hap2_node  outgoing_hap2_node`

Tangle_traverer does not support tangles with more than 2 traversing paths (i.e. most of the rDNA tangles)



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


