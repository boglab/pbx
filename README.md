# PBX
A fully automated pipeline for TAL effector RVD sequence determination from raw PacBio data

## Requirements
* [AMOS 3.1.0](http://amos.sourceforge.net/wiki/index.php/AMOS)
* [BLAST command line tools 2.2.28+](http://www.ncbi.nlm.nih.gov/books/NBK279671/)
* [GNU Parallel 20130722+](http://www.gnu.org/software/parallel/)
* [MUMmer 3.23](http://mummer.sourceforge.net/)
* [Python 2.6/2.7](https://www.python.org/)
* [SMRTAnalysis 2.3](http://www.pacb.com/devnet/)

## Usage

To run the pipeline edit config.ini then do:

```python pbx.py config.ini```

config.ini has only 3 required parameters:

* `smrtanalysis_path` is the full path to SMRTAnalysis
* `mummer_path` is the full path to MUMmer
* `raw_reads_path` is the full path to a folder containing raw PacBio reads in .bas.h5 and .bax.h5 format. 

Other parameters of possible interest:

* `results_path` is the directory the results will be stored in
* `tale_seqs_file_whitelisting` is the sequence file in tale_seqs/whitelisting that raw reads will be aligned to identify TALE-containing reads. As Xanthomonas TALEs are all highly similar in sequence, the default set from Xoc will identify nearly all reads even in data sets from other Xanthomonas species.
* `tale_seqs_file_export`
  Automated TALE sequence extraction from assembled reads relies on identifying conserved TALE N-terminal and C-terminal coding regions.
  This is the file in tale_seqs/exporter containing known terminal sequences that will be used.
* `tale_seqs_file_export_boundaries`
  RVD sequence determination breaks apart repeat regions and identifies RVDs based on conserved boundary residues.
  This is the file in tale_seqs/exporter containing known boundaries that will be used.

## Results

After the pipeline is finished running, the determined RVD sequences will be at `results_path/resequencing/unique_tale_seqs.txt` and `results_path/combine_resequenced_tals/unique_tale_seqs.txt`.

These files should be interpreted as discussed in Booher et al.

If the number of identified TALEs seems low it may be that your library insert size was too small to produce a useful number of long reads at the 16 kbp threshold. Run the pipeline again using a lower value for `min_seed_read_length` such as 12000 or 10000.
