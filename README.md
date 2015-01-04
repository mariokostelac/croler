Mile assembly
=============

Instalacija
-----------

    git submodule update --init --recursive
    make

Rikvajrmenti
------------
(molim nekog da prevede ovo, ne mogu se sjetit hrvatske rijeci)
- gcc (4.7.2 na debianu fino kompajlira)
- make
- ~~zlib (dev verzija - treba nam source)~~ (otkad smo se maknuli od citanje fasta/fastq, ovo vise ne trebamo)

Pokretanje
----------
    ./run.sh <reads.afg>
