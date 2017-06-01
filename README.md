# lattice-word-length-distribution

This tool computes the distribution of the length of the transcription of lattices.
That is, what is the probability that the transcription length (in terms of words) is 1, 2, 3, 4, etc.

Lattices must be in Kaldi's format. The program will output one line for each lattice, with the ID
of the lattice, followed by the distribution of the lenghts (sorted in decreasing order). 

```
lattice1 length log-prob ; length log-prob ; ...
lattice2 length log-prob ; length log-prob ; ...
```

You can find a full example [here](egs/README.md).

### Options:

- `--acoustic-scale`    : Scaling factor for acoustic likelihoods in the lattices.
- `--graph-scale`       : Scaling factor for graph probabilities in the lattices.
- `--insertion-penalty` : Add this penalty to the lattice arcs with non-epsilon output label (typically, equivalent to word insertion penalty).
- `--nbest`             : Limit the distribution to this number of n-best lengths.
```
