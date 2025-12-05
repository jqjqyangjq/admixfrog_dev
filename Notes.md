# Original implementation

Assume the error rates are the same for all types of substitutions, then for the probability of the data (observe a number of alt alleles) we have:

$$
P(O \mid N, p) = \binom{N}{O} p^{O} (1-p)^{N-O}
$$

$N$ denotes total number of observed alleles, and $O$ denotes the number of observed alt alleles.  
$p$ denotes the probability of observing an alt allele given any allele. Therefore we further have:

$$
g = [0, 1, 2]
$$

$$
p = c \cdot P_{\text{cont}} + (1-c) \cdot \frac{g}{2}
$$

$$
p = p(1-\text{error}) + (1-p)\text{error}
$$

Note that `error` accounts for the uniform error rate. This is implemented in the original functions in *admixfrog*.



# Different error rates

Assume the error rates are different for each observed allele, as substitution rates are indeed different  
(deamination / sequencing errors). Then, for all the reads/alleles, we have the total likelihood of the data  
(observe a number of alt alleles) at this position:

$$
L[g] = \prod_i P(\text{read}_i \mid G = g)
$$

For a given read $\text{read}_i$, we have:

$$
P(\text{read}_i) =
\begin{cases}
p   & \text{if read}_i = \text{alt} \\
1-p & \text{if read}_i = \text{ref}
\end{cases}
$$

And $p$ is calculated in a similar way as above, but the error rate is now substitution-dependent  
(ref base, alt base, strand of the allele, observed allele):

$$
p = p(1 - \text{error}_{\text{current}}) + (1-p)\,\text{error}_{\text{current}}
$$
