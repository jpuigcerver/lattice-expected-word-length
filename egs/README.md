# Example for lattice-word-length-distribution

Suppose that we have some lattice from we wish to find the distribution
of the transcript length. That is, the probability that the length of the
transcript (in words) is 1, 2, 3, 4, etc.

```
lat1
0	1	1	0.35667494393873237891,0.0,1_1_1
0	2	0	1.60943791243410037460,0.0,
0	2	2	2.30258509299404568401,0.0,2_2
1	3	2	0.69314718055994530941,0.0,2_2_2_2
1	4	2	0.69314718055994530941,0.0,2_2_2
2	4	1	0.0,0.0,1_1_1_1
3	4	3	0.22314355131420975576,0.0,3_3
3	1.60943791243410037460,0.0,
4	0.0,0.0,
```

This lattice has the following paths with the associated probability:

- P(`a`)  = 0.2
- P(`ab`)  = 0.35 + 0.07 = 0.42
- P(`ba`)  = 0.1
- P(`abc`) = 0.28

Thus, the distribution of the length is:

- P(L = 1) = P(`a`) = 0.2
- P(L = 2) = P(`ab`) + P(`ba`) = 0.42 + 0.1 = 0.52
- P(L = 3) = P(`abc`) = 0.28

This is exactly what the program will compute:

```
$ ./lattice-word-length-distribution ark:egs/lattice1.txt
lat1 2 -0.6539416313 ; 3 -1.272688508 ; 1 -1.609786749 ;
```

Where the log-probabilities are shown for each length. The lengths are sorted
from the most likely to the least. You can choose to show only the
n-best lengths with the `--nbest`.

Of course, the acoustic scale, language model scale and word insertion penalty
values will affect the distribution. By default no scaling or penalty is
applied, but you can change this with the appropriate option.
