# Original implementation

Assume the error rates are the same for all types of substitutions, then for the probability of the data (given a set of genotype, to observe a number of alt alleles from reads) we have:

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

Note that `error` accounts for the uniform error rate. `p` is then used to calcuapte the emission probability as above. This is implemented in the original functions in *admixfrog*.



# Different error rates

Assume the error rates are different for each observed allele, as substitution rates are indeed different  
(deamination / sequencing errors). Then, for all the reads/alleles, we have the total likelihood of the data  
(given a set of genotype, to observe a number of alt alleles from reads) at this position:

$$
L[g] = \prod_i P(\text{read}_i \mid G = g)
$$

For a given read $\text{read}_i$, if the allele seen is alt, we have:

$$
g = [0, 1, 2]
$$

$$
p = c \cdot P_{\text{cont}} + (1-c) \cdot \frac{g}{2}
$$

$$
p = p(1-\text{error}) + (1-p)\text{error}
$$ 

Considering that the allele seen here can be ref or alt, we have:

$$
P(\text{read}_i) =
\begin{cases}
p   & \text{if read}_i = \text{alt} \\
1-p & \text{if read}_i = \text{ref}
\end{cases}
$$

The error now depends on ref base, alt base, strand of the allele, observed allele.  
Typically,
for a forward strand T allele at a c/t position, the error will range between 40% to 2% for non-UDG-SS lib.  
for most others, I (arbitrarily) set to 0.002 (can be any other baseline).  
e.g. For a C(ref)-T(alt) sites, if there are 3 ref and 2 alts (one +, one -) observed  
the error rate matrix can be [0.002, 0.002, 0.002, 0.002, f(C-T at this pos)]
