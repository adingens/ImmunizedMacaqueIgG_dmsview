# Mutational antigenic profiling of purified IgG from immunized macaques

Adam S. Dingens, Jesse D. Bloom, In collaboration with Chris Cottrell and Andrew Ward

We are performing mutational antigenic profiling of sera from macaques vaccinated with BG505 trimers, using the BG505.T332N mutant Env libraries, first described and characterized in [Haddox, Dingens et al 2018](https://elifesciences.org/articles/34420). In some cases, the trimer has 241 or 289 holes glycosylated and the RM19R, a non-NAb targeting the highly immunodominant peptidic epitope on the exposed base of the trimer, was pre-incubated with trimer to mask this epitope. Note that the RM19R epitope is not present on live virus used in mutational antigenic profiling.

[This page](https://jbloomlab.github.io/dms_tools2/diffsel.html) documents the differential selection statsitics we use to analyze these data.


The site-metrics (dot plot) include:

- **positive diffsel**: The sum of all positive differential selection values at a site. This gives a sense to the total amount of escape/selective pressure at each site.
- **negative diffsel**: The sum of all negative differential selection values at a site. This gives a sense for mutations that are depleted, rather than enriched, during serum selection relative to a non-selected control library. It is intriguing that many of these potential serum sensitizing mutations cluster and are consistent across sera.
- **max diffsel**: The value of the largest effect mutation (largest mutation differential selection) at each site. This gives a sense of the maximal effect of any mutation at each site.
- **min diffsel**: The value of the smallest effect mutation (smalled mutation differential selection) at each site.
- **abs diffsel**: The absolute value of all mutation differential selection values at each site.

The mutation-metrics (logoplot) include

- **diffsel**: All mutation differential selection values, including negative values (which we oftentimes do not analyze rigorously for technical reasons), are plotted.



## Preliminary interpretation
I do not know which animals were in the experimental or ocntrol groups. However, there was a qualitative differentiaion between  groups of samples. I am not very familiar with study details, and I encourage you to explore these data beyond the crude interpetations outlined below. For example, some sera appear to target the glycan hole slightly more than others, but I did not comment on this directly as it is not totally clear from this data alone. 

CG41 and 33311 both strongly target the C3/V5 epitope, and also appear to target the V3 epitope.

33433 and 34945 both target the C3/V5 epitope, but likely also target a number of other epitopes such as V3, V1/V2, others (where signal:noise is not completely resolved). I would guess these are relatively more polyclonal, but this is speculative based on this data alone.

MA255, 33203 give very poor signal:noise, suggesting these neutralizing antibody responses could be relatively more polyclonal. 

Note that this is only the first set of data; most data is from single replicates, while few are averaged across two replicates (33311, CG41). It does *not* appear that the cleaner data for these samples if due to averaging noise, as single replicate data from these samples are quite clean. 