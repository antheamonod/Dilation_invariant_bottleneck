## Dilation-invariant bottleneck dissimilarity

DI bottleneck dissimilarity is defined as 

$$
D(A\|B) = \min_{c\ge 0}d_{bo}(cA,B)
$$

Note that $D(A||B)$ is asymmetric. Changing the variable order yields different dissimilarities. In comparison of persistence diagrams, put query at the second place to ensure everything is measured at the same scale. i.e.

$$
\arg\min_{\rm data_i} D(\rm data_i\| query)
$$

$\tt{DI\_dissimilarity.py}$ computes the DI dissimilarity.

## Dilation-invariant bottleneck distance

DI bottleneck distance is defined as 

$$
\bar{D}([A],[B]) = \bar{d}_{bo}(\log(A),\log(B))
$$

where $\bar{d}_{bo}$ is the shift-invariant distance defined by [Don Sheehy etc.](http://donsheehy.net/research/cavanna18computing.pdf). There are two algorithms to compute this distance:

- Naive Searching: fast but may not be accurate (precision depends on number of partitions). The python code is $\tt{DI\_distance.py}$
- Don's algorithm: accurate (may put as ground truth) but very slow. **Be careful NOT** to put full dataset into this algorithm. It takes a week to run mammal_euclidean_* data. The main function is $\tt{Don\_distance.py}$. Other dependent functions (written by Don) are 
  - $\tt{main\_algorithm.py}$,
  - $\tt{event\_queue.py}$,
  - $\tt{bipartite\_matching.py}$,
  - $\tt{plane_util.py}$

## Datasets

datasets are listed with $\tt{*.npz}$.
