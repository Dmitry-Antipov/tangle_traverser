#!/usr/bin/env bash
test=helo1
rm tmp/*
./tangle_traverser.py --gfa tests/$test/graph.noseq.gfa --coverage tests/$test/ont-coverage.csv --boundary-nodes tests/$test/border.ids --alignment tests/$test/alignments.gaf --outdir tmp --tangle-file tests/$test/tangle.ids   --num-initial-paths 3
diff tmp/traversal.gaf tests/$test/expected_traversal.gaf
python3 diff_compare.py tmp/traversal.gaf tests/$test/expected_traversal.gaf
