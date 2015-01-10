# croler

*croler* (please, all letters lowercased) is genome OLC *de novo*
assembler written on [FER](http://www.fer.hr), under leadership of [Mile
Sikic](http://complex.zesoi.fer.hr/msikic.html).

As most OLC assemblers do, *croler* consists of three different phases:
- [overlap](pipeline/qpid/README.md)
- [layout](pipeline/brahle_assembly/README.md)
- [consensus](pipeline/msa/README.md)

*[for details follow the links above].*

Requirements
------------
- gcc (4.7.2 on debian will do the job)
- make

Installation
-----------
Change the directory where the project is cloned and make sure that
following commands flow through your command line

    git submodule update --init --recursive
    make

When build process finishes (successful, eventually), the `bin`
directory will contain symlinks to binaries of different stages:
`croler_overlap`, `croler_layout` and `croler_consensus`.
We think that names are verbose enough to not explain it which binary
refers to which stage.

Running
-------
*croler* reads input in `.afg`
([AMOS](http://amos.sourceforge.net/wiki/index.php/Message_Types)) format.

It can be run multiple ways.
The simplest way is to use `run.sh` script like

    ./run.sh <reads.afg>
The example runs the all phases, pipelined, with default parameters.
Resulting contigs will be saved in `{reads}_consensus.fasta` file.
Feel free to explore other side-effect files (intermediate results),
such as `{reads}_overlaps.afg`, `{reads}_layout.afg` and various `.dot`
graphs created by layout phase.

*{reads} refers to a prefix of reads filename*

**Note**: It is the simplest way to run the assembly, but running each
phase separetely is more flexible and allows adjusting parameters to
your test data. We really encourage you to explore the running options
of each phase.

Results
-------
Jump to [tests](tests/README.md).

Contributors
------------
- Bruno Rahle
- Anton Grbin
- Mario Kostelac (mario.kostelac@gmail.com)
- Marko Culinovic (marko.culinovic@gmail.com).

Do not hesitate to drop an email to Marko and Mario, if you have any
questions.
