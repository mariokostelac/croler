Comparison tests
================

This directory contains comparions tests between
[minimus](http://amos.sourceforge.net/wiki/index.php/Minimus)
(one of AMOS
assemblers) and `croler`.

Each directory represents particular pipeline

- **o** stands for `croler overlap | minimus layout | minimus consensus`
- **ol** stands for `croler overlap | croler layout | minimus consensus`
- **olc** stands for `croler overlap | croler layout | croler consensus`

Running
-------
Each directory contains `run.sh` script that will run particular
pipeline against given file (*reads in .afg format*).

It will also trigger minimus on the same input to make comparison
easier.

Example:
```bash
cd tests/o/influenza-A/
../run.sh ../../test_data/influenza-A/influenza-A.afg
```
will run **o** pipeline and **minimus** on the reads of influenza-A
genome. Running such a script will produce files **contigs.we** and
**contigs.minimus**.

Comparing these two could be useful for observing how good is croler
overlapper.

**NOTE**: `run.sh` assumes that amos and croler bin directories are part
of `$PATH` var.

Results
-------
All results are calculated and written with [v1.0RC1]
(https://github.com/mariokostelac/croler/archive/v1.0RC1.zip), all default params.

#TODO
