# Multilevel-Models-Monte-Carlo-Simulation

### What is this?

This is the script from my dissertation. It is written in SAS. It was my first step into creating original code - written under very patient tutelage from my dissertation chair.
While I did my then-best to make it readable, it is not optimized to run quickly. Knowing what I know now, I could have - for example - programmed it in Python and made it run much faster using things like list comprehensions. 
But, while it's not pretty it got the job done, and was a terrific introduction into thinking about statistical programming.

### What does it do (in _English_)?



### What does it do?
#
The script has 3 main steps:

1. Simulate data from various populations.
1. Analyze that data using some variant of SAS's `proc mixed` function.
1. Synthesize the results from those analyses.


To accomplish this, it relies on 3 respective user-defined macros(i.e. functions):

1. `mixedmod`
1. `mixedanalyze`
1. `mixedresults`

`%mixedmod(condition, reps, clustdist, subs, nclust, ICC, L1effect, L2effect, L12effect)`

	condition: Takes an integer (in the dissertation this ranged between 1-6144) and is the identifier for a specific population's parameters.

	reps: Takes a positive integer. The number of data sets to simulate under a population's parameters. In the dissertation this was uniformly set at 1000 to minimize sampling error.

	clustdist: Takes either 1 (for a uniform distribution of cluster sizes) or 2 (for a truncated Poisson distribution of cluster sizes).

	subs: Takes a positive integer. The m value in mixed model notation. Returns an equivalent number of rows (i.e. subjects) PER CLUSTER for clustdist = 1, or some slightly-randomized number of rows around this value for clustdist = 2. Ranged between 5 and 20 in the dissertation.

	nclust: A positive integer. The k term in mixed model notation. Specifies the number of clusters (i.e. subgroups) within the dataset. Thus, the total number of rows in a given dataset is roughly equal to subs * clustdist.

	ICC: A positive float value. Denoted with rho in mixed model notation. In practice between 0 and 1; (0 and 0.8 in the dissertation) although no out-of-range error is raised for inputs exceeding this value. The inter-class correlation coefficient, which is a measure of the proportion of variance in Y due to clustering between rows belonging to the same subgroup, or statistical non-independence. The higher the ICC, the more values of a row are similar to values of other rows in the same cluster. This occurence is the primary motivation for using mixed-model analysis.

	L1effect: The level 1 (i.e. row-level) effect size. A positive float value. In the dissertation this was one of 0, 0.2, 0.5, and 0.8 to simulate no, small, medium and large effects (cohen's D) respectively. However, there is no range error raised for values outside of this range, and such values are sometimes observed in certain sciences.

	L2effect: The level 2 (i.e. cluster-level) effect size. Takes the same range of values as L1effect. This effect is the impact at the subgroup/cluster level in a dataset.
	
	L12effect: The interaction effect between the L1 and the L2 effects. A positive float value in the dissertation set between 0 and 0.8. This argument quantifies how much statistical interaction there is between the row and the subgroup effects. Note it is possible to have null effects for the L1 and L2 effect, but still see large interaction effects, though this is rare in practice.

Example usage: `%mixedmod(condition=1 , reps=10, clustdist=1, subs=5, nclust=25, ICC=.3, L1effect=0.5, L2effect=0.2, L12effect=0)`

This would return 10 datasets labeled with condition label 1. Each dataset would contain 25 uniformly distributed clusters with 5 rows per cluster (total number of rows: 125). 
Roughly 30% in the variance on the outcome (Y) variable would be due to statistical non-independence. A medium-sized effect would be simulated across all datasets at the row level, a small effect at the cluster-level, but with no interaction between them.  

`%mixedanalyze(start, stop, type)`

Returns analysis on datasets created by mixedmod, using SAS's `proc mixed` procedure. Outputs single files across start-stop batches. See mixedresults below.

	start: Defines which condition's datasets to begin analysis on.

	stop: Defines which conditions's analyses to end analysis on. In conjunction with the start argument, the stop argument defines the overall batch size for analysis.

	type: Integer between 1 and 4. Null models feature no input variables, and are necessary for computing the empirical ICC values. Regression-only models, are necessary to compute how much bias was produced in estimates by ignoring statistical non-independence. Regression models under the proc mixed procedure are specified by omitting the /random line.

		1 Requests null, regression, and mixed-model analysis for all datasets in the batch. Outputs and saves results to file specified under libname, instead of the SAS Output Delivery System (ODS) GUI. 
			The option used during the dissertation analysis.

		2 As above, except instead of saving output to file, output is instead rendered using SAS's ODS. Using this option with a "large" batch of data may cause an OOM situation and cause SAS to crash or freeze.

		3 Requests only [full] mixed-model analyses for all datasets in the batch. Renders using ODS. Useful for testing code.

		4 Requests only regression and [full] mixed-model analyses for all datasets in batch. Renders using ODS. Useful for quickly observing differences in estimates between methods.

Example usage: `%mixedanalyze(start = 64, stop = 128, type = 1)`

This would request null, regression, and mixed-model analysis using the 3 input (i.e. predictor) variables specified in 'mixedmod` above. The results would be output and saved to folders located at "C:\SAS_MC\Test\" by default. SAS produces several difference output tables. The current script requests the saving of NFixed, NCov, Fixed, Random, Coverge, Interate, NObs, FitStats, Type3, CovParms, and Clust output tables.
Note that in this case the call is requesting analysis on 128-64 = 64 conditions. Assuming 1,000 datasets per condition by convention, this is 64,000 datasets (with subs * nclust rows per dataset) for analysis.
For further information on `proc mixed` in SAS and its associated output, [go here.](http://support.sas.com/documentation/cdl/en/statug/66859/HTML/default/viewer.htm#statug_mixed_overview.htm)

`mixedresults(batch, type)`

Returns estimates of the ICCs and design effects averaged across datasets in each condition. 

	batch: Identifies the conditions in the batches from the mixedanalyze call previously. Specified with 'a_to_b' string, where a is the value of 'start' and b is the value of 'stop' in the mixedanalyze macro above.

	type: Specifies which aggregated estimations to perform. Currently only type = 0 is available.

Example usage: `%mixedresults(batch = 64_to_128, type = 0)`

This call would return estimates for the 64_to_128 file created by `mixedanalyze` above.

The output - delivered in Excel format - specifies the conditions and the characteristics of the population of that condition, as specified in `mixedmod` above. 
More importantly, it also contains:
1. The empirical (i.e. observed) ICC, averaged across all datasets in that condition  
1. The empirical design effect for the overall model, which may be different than the expected/calculated design effect, due to other factors in the population.
1. The empirical design effect contributed by each term in the model (including the intercept). 

The actual results from all this (and a more meaningful walkthrough) might potentially be available if I can post a draft of my disseration here. This may be included in a future update.
