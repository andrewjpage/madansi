# madansi
Reorders and orientates contigs from a fasta file using a pan genome graph.
## Installation

## Example usage
Run the following to get usage:
`madansi -h`

An example:
```
madansi 'madansi/tests/data/pan_genome_reference.fa' 'madansi/tests/data/assembly_4_sequences.fa' 'madansi/tests/data/graph_5_nodes.dot'
			--output_fasta_file '/tmp/ordered_contigs.fa'
			--output_directory '/tmp/all_files'

```
## Tests
Run `python setup.py test`
## Known Issues