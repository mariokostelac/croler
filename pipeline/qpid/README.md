# qpid

*qpid* is a part of *croler*, serving the purpose of seeking for good
overlaps between given sequences. It is shamelessly inspired by *minimus
hash-overlap*.

*hash-overlap* does the job pretty much, but it is slow and  code
quality is far below today's standards so we created our implementation
of same algorithm (that's how qpid is born); with code quality and
efficiency in mind.

## Installation
hit `make` and that's it.

Also, it can be built recursively by triggering a global build for
*croler* (explained in *croler*'s instructions).

## Usage
Usage is pretty simple, just hit run `./bin/overlap` and you'll get

    SYNOPSIS
./bin/overlap [options] <reads.afg>

    DESCRIPTION
    The   following options are available:
    -t number of threads
    -a alignignment band radius
    -w offset wiggle
    -r merge radius
    -bande maximum error rate
    -o output file

For explanation of each argument meaning, proceed to
[algorithm](#algorithm) section (after explanation of the core algorithm).

## Input/output formats
*qpid* reads sequence reads in *.afg* format and outputs overlaps in the
same manner, in *.afg*.

Details about that format can be found on [AMOS
wiki](http://sourceforge.net/apps/mediawiki/amos/index.php?title=Message_Types#Overlap_t_:_Universal_t)

## Algorithm
In the first phase it uses
[minimizers](http://bioinformatics.oxfordjournals.org/content/20/18/3363.long)
for getting the information which pairs of reads could be good
candidates for calculating overlaps. This step is really necessary;
otherwise, we would have to calculate `n * n` overlaps, where `n` is the
size of reads set. Imagine assembling something that's broken on 10^6
reads (which is not much). We would never be finished.

After filling our database with minimizers, we use them to approximate
which parts of a read could align to which parts of some other read
(*region*).
It really matters if the real overlap looks like 
```
---------->
  ---------->
```
or
```
---------->
      ---------->
```
or
```
      ---------->
---------->
```

After we decide which parts (*regions*) are good to be aligned, we use
that information to create bands for *banded overlap*. We are not going
into details of banded overlap since it is usually explained on
bioinfromatics courses. One of the resources for learning could be
[lectures](http://www.site.uottawa.ca/~lucia/courses/5126-10/lecturenotes/03-05SequenceSimilarity.pdf)
from University of Ottawa.

Also, it is important to know that *qpid* returns the best read for each
pair, if error below `error_rate` parameter.

### What about all that parameters?
After we find bounds, we like to add algorithm some more space for
searching so we make these bounds a bit looser and extend them by
`alignment band radius`; both sides.

During transforming set of *regions* for a pair `(x, y)`, it could
happen that two *regions* do not agree about bands. If differences
between them are less than `merge_radius` we make merge them like they
fully agree. This step just improves the efficiency since we calculate
less *banded overlaps* for `(x, y)` pair.

`offset_wiggle` is very similiar as `merge_radius`, but it occurs more
often because it tries to merge two *regions* right after they are
created. Therefore, it has to be more strict than `merge_radius`.

After we find the best overlap for a pair `(x, y)`, it has to have an
`error_rate` below the defined value. Contrary to *hash-overlap*, *qpid*
approximates the error_rate in order to simplify *banded overlap*
calculation.

### Still not clear?
Create an issue or drop me an email.

## Found a bug, typo or anything wrong?
Make a pull request or raise an issue.
