# Layout phase

In layout phase of assembly, we tend to find as many contigs as possible from overlaps given by *qpid*, and we prefere if these contigs length is as long as possible. It is inspired by *SGA* ([https://github.com/jts/sga](https://github.com/jts/sga)) and uses String graph for reads and overlaps storage. Also, it uses the same methods as *SGA* for graph simplification, i.e. trimming and bubble popping.

## Installation
hit `make` and that's it.

Also, it can be built recursively by triggering a global build for
*croler* (explained in *croler*'s instructions).

## Usage
Usage is pretty simple, just hit run `./bin/croler_layout` and you'll get

```
SYNOPSIS
./bin/croler_layout [options] <reads.afg> <overlaps.afg>
```


```
DESCRIPTION
The following options are available
    -t   number of trimming rounds
    -b   number of bubble popping rounds
    -h   trimming read length threshold
    -n   maximum number of nodes during bfs in bubble popping
    -d   maximum walk sequence length in bubble
    -w   maximum number of walks in bubble
    -a   maximum diff between aligned bubble walk sequences
```

For explanation of each argument meaning, proceed to
[algorithm](#algorithm) section (after explanation of the core algorithm).

## Input/output formats
*Layout phase* reads sequence reads and overlaps in *.afg* format and outputs contigs in the
same manner, in *.afg*.

Details about that format can be found on [AMOS
wiki](http://sourceforge.net/apps/mediawiki/amos/index.php?title=Message_Types#Overlap_t_:_Universal_t)

## Algorithm

First, you need to know what types of overlaps you're handling, so check [overlap types](http://sourceforge.net/p/amos/mailman/message/19965222/). As in *minimus*, *croler_layout* handles only overlaps of type normal and innie.

Second thing you need to know, what is *String graph*. You can read section 2.3. from paper [Efficient construction of an assembly string graph using the FM-index](http://bioinformatics.oxfordjournals.org/content/26/12/i367.full.pdf+html) by Simpson and Darbin, but i will also try to put it simple here.

Let's suppose you have 3 reads (example from mentioned paper):
```
R1 ACATACGATACA
R2    TACGATACAGTT
R3       GATACAGTTGCA
```

How do you build string graph from reads? You store read data into vertices and
build bidirectional egde with two different labels marking directions from first read to second and reversed.
The simplest way to understand labels is this one - substring of second read you need to concatenate to first read so that you get merged sequence of two reads.
Example for overlap between reads R1 and R2 you'll build edge with labels ``GTT`` for direction one to two and ```ACA``` for direction two to one. Still not clear? - see picture below or read mentioned paper.

So, how does algorithm works? Foremost, we remove reads which are containment (see [overlap types](http://sourceforge.net/p/amos/mailman/message/19965222/)) and reads which are unnecessary because of transitive overlaps.

Transitive overlaps example:

![alt text](https://github.com/mariokostelac/croler/blob/master/images/transitive_graph.png "Transitive overlaps")

As you see in picture above, there exists overlap between reads (R1,R2), (R2,R3) and (R1,R3). All data which is necessary is stored in reads R1 and R3 and their overlap so we can remove read R2 and overlaps (R1,R2) and (R2,R3).

Next step is removing tips and disconnected vertices. Tips are reads which have overlaps in only one directions, i.e. only prefix or only suffix of read is part of overlaps. They emerge due to errors at the edges of reads. Disconnected vertices are reads which doesn't have any overlaps. This process is called *trimming*.

Tip example:
```
R4    TACGATACAGTA
R1 ACATACGATACA
R2    TACGATACAGTT
R3       GATACAGTTGCA
R5          ACAGTTGCACCT
```
And also, lets suppose that reads R1, R2, R3 and R5 have other overlaps with different directions which aren't listed here.

![alt text](https://github.com/mariokostelac/croler/blob/master/images/trimming_graph.png "Tip removal")


But, we want to simplify graph even more. We consider tho walks in graph redundant if they start and end at the same vertex(read) and contain similar sequences. This kind of structure is called *bubble* and emerges due to internal read errors or to nearby tips connecting.

Bubble example:

![alt text](https://github.com/mariokostelac/croler/blob/master/images/bubbles_graph.png "Bubble example")

Bubbles are found with [BFS algorithm](http://en.wikipedia.org/wiki/Breadth-first_search) starting in every vertex of graph. Vertices in bubble must not have edges to vertices outside the bubble. Also, direction of last edge in all walks must be same. Afterwards, sequences representing walks in bubble are extracted and aligned using [edlib](https://github.com/Martinsos/edlib). If sequences are similar enough, walk with greatest coverage is preserved. Vertices and edges of other walks are deleted from graph. Process of removing bubbles from graph is called *bubble popping*.

### What about all that parameters?
There are many options/parameters for *layout phase*. Options ```-t``` and ```-b``` are used to specify how many times you want to start trimming and bubble popping algorithms.

In trimming, you must specify read length threshold(```-h``` option), because you don't want to remove reads which represent a relatively larger part of genome even if they have overlaps in only one direction.

In bubble popping, there are several parameters.

First one is ```-n``` for maximum number of nodes created during BFS. Node in BFS is simply wrapper for vertex, therefore two nodes can wrap same vertex. Second one is ```-d``` for maximum sequence length extracted from walk in bubble. These two parameters are used to constrain search because you don't want to iterate through complete set of reads every time you're searching for bubble.

Another parameters for bubble poping are:
- ```-w``` if you want to limit number of different walks in bubble and
- ```-a``` for setting maximum difference percentage in score obtained by aligning two walk sequences. Difference percentage is calculated as *alignment_score / sequence_length*

### Still not clear?
Create an issue or drop me an email.

## Found a bug, typo or anything wrong?
Make a pull request or raise an issue.
