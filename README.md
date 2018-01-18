# Multilevel-Models-Monte-Carlo-Simulation

### What is this?

This is the script from my dissertation. It is written in SAS. It was my first step into creating original code - written under very patient tutelage from my dissertation chair.
While I did my then-best to make it readable, it is not optimized to run quickly. Knowing what I know now, I could have - for example - programmed it in Python and made it run much faster using things like list comprehensions. 
But, while it's not pretty it got the job done, and was a terrific introduction into thinking about statistical programming.

### What does it do (in _English_)?



### What does it do?

The script has 3 main steps:

1. Simulate data from various populations.
1. Analyze that data.
1. Synthesize the results from those analyses.


To accomplish this, it relies on 3 respective user-defined macros(i.e. functions):

2. `<mixedmod>`
2. `<mixedanalyze>`
2. `<mixedresults>`

mixedmod(condition, reps, clustdist, subs, nclust, ICC, L1effect, L2effect, L12effect):

	condition: Takes an integer (in the dissertation this ranged between 1-6144) and is the identifier for a specific population's parameters.

	reps: Takes a positive integer. The number of data sets to simulate under a population's parameters. In the dissertation this was uniformly set at 1000 to minimize sampling error.

	clustdist: Takes either 1 (for a uniform distribution of cluster sizes) or 2 (for a truncated Poisson distribution of cluster sizes).

	subs: Takes a positive integer. Returns an equivalent number of rows (i.e. subjects) PER CLUSTER for clustdist = 1, or some slightly-randomized number of rows around this value for clustdist = 2. Ranged between 5 and 20 in the dissertation.

	nclust: A positive integer. Specifies the number of clusters (i.e. subgroups) within the dataset. Thus, the total number of rows in a given dataset is roughly equal to subs * clustdist.

	ICC: A positive float value. In practice between 0 and 1; (0 and 0.8 in the dissertation) although no out-of-range error is raised for inputs exceeding this value. The inter-class correlation coefficient, which is a measure of the proportion of variance due to clustering between rows belonging to the same subgroup, or statistical non-independence. The higher the ICC, the more values of a row are similar to values of other rows in the same cluster. This occurence is the primary motivation for using mixed-model analysis.

	L1effect: The level 1 (i.e. row-level) effect size. A positive float value. In the dissertation this was one of 0, 0.2, 0.5, and 0.8 to simulate no, small, medium and large effects (cohen's D) respectively. However, there is no range error raised for values outside of this range, and such values are sometimes observed in certain sciences.

	L2effect: The level 2 (i.e. cluster-level) effect size. Takes the same range of values as L1effect. This effect is the impact at the subgroup/cluster level in a dataset.
	
	L12effect: The interaction effect between the L1 and the L2 effects. A positive float value in the dissertation set between 0 and 0.8. This argument quantifies how much statistical interaction there is between the row and the subgroup effects. Note it is possible to have null effects for the L1 and L2 effect, but still see large interaction effects, though this is rare in practice.

Example usage: %mixedmod(condition=1 , reps=10, clustdist=1, subs=5, nclust=25, ICC=.3, L1effect=0.5, L2effect=0.2, L12effect=0):

This would return 10 datasets labeled with condition label 1. Each dataset would contain 25 uniformly distributed clusters with 5 rows per cluster (total number of rows: 125). 
Roughly 30% in the variance on the outcome (Y) variable would be due to statistical non-independence. A medium-sized effect would be simulated across all datasets at the row level, a small effect at the cluster-level, but with no interaction between them.  

mixedanalyze(start, stop, type):

	start: Defines which condition's datasets to begin analysis on.

	stop: Defines which conditions's analyses to end analysis on. In conjunction with the start argument, the stop argument defines the overall batch size for analysis.

	type: Integer between 1 and 4. All models use SAS's `<proc mixed>` procedure, with the /random line omitted for regression-only models. Null models are necessary for computing the empirical ICC values. Regression-only models, are necessary to compute how much bias was produced in estimates by ignoring statistical non-independence.

		3.  