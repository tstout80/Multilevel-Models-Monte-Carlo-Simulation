/*********************************************/
/* Stout Mixed Model
/*********************************************/

options symbolgen;

%macro mixedmod(condition, reps, clustdist, subs, nclust, ICC, L1effect, L2effect, L12effect);

proc printto log = "C:\SAS_MC\Logs_1\condition_&condition..log";
run;

/* Uniform Distribution */


%do i = 1 %to &reps;
%if &clustdist = 1 %then %do;

data one;
seed = -1;

i = &i;
nclust = &nclust;
subs = &subs;
ICC = &ICC;
L1effect = &L1effect;
L2effect = &L2effect;
L12effect = &L12effect;
intercept = 50;
reps = &reps;
clustdist = &clustdist;

do clust = 1 to &nclust;
                indvarnum = (25*ICC);
                indvarden = (1-ICC);
                intvar =  sqrt(indvarnum/indvarden)*rannor(-1);
                SD = sqrt(16 + 25);
                b01 = L2effect*SD;
                x2 = rannor(123);
                do indiv = 1 to &subs;
                        err = sqrt(25)*rannor(-1);
                        b10 = L1effect*SD;
                        b11 = L12effect*SD;
                        x1 = rannor(-1);
                        x3 = x1*x2;
                        y = intercept + b10*(x1) + b01*(x2) + b11*(x3) + intvar + err;
                        output;
keep reps y i x1 x2 x3 clustdist clust subs nclust indiv ICC L1effect L2effect L12effect b10 b01 b11 err indvarnum indvarden intvar SD;
                end;
end;
run;

%end;


/* Binomial distribution */

%if &clustdist = 2 %then %do;

data one;
seed = -1;

i = &i;
nclust = &nclust;
subs = &subs;
ICC = &ICC;
L1effect = &L1effect;
L2effect = &L2effect;
L12effect = &L12effect;
intercept = 50;
reps = &reps;
clustdist = &clustdist;

do clust = 1 to &nclust;
                indvarnum = (25*ICC);
                indvarden = (1-ICC);
                intvar =  sqrt(indvarnum/indvarden)*rannor(-1);
                SD = sqrt(16 + 25);
                b01 = L2effect*SD;
                x2 = rannor(123);
                do indiv = 1 to ranbin (-1, 1000000, %sysevalf(&subs/1000000));
                        err = sqrt(25)*rannor(-1);
                        b10 = L1effect*SD;
                        b11 = L12effect*SD;
                        x1 = rannor(-1);
                        x3 = x1*x2;
                        y = intercept + b10*(x1) + b01*(x2) + b11*(x3) + intvar + err;
                        output;
keep reps y i x1 x2 x3 clustdist clust subs nclust indiv ICC L1effect L2effect L12effect b10 b01 b11 err indvarnum indvarden intvar SD;
                end;
end;
run;

%end;

data out_data;
set one;
file "C:\SAS_MC\Data_1\Condition_&condition..txt" mod;
put reps y i x1 x2 x3 clustdist clust subs nclust indiv ICC L1effect L2effect L12effect b10 b01 b11 err indvarnum indvarden intvar SD;
run;

%end;
%mend mixedmod;

%mixedmod(      condition=      1      ,      reps=10,      clustdist=1,      subs=5,      nclust=25,      ICC=.1,      L1effect=0,      L2effect=0,      L12effect=0)
%mixedmod(      condition=      2      ,      reps=10,      clustdist=1,      subs=5,      nclust=25,      ICC=.3,      L1effect=0,      L2effect=0,      L12effect=.2)
%mixedmod(      condition=      3      ,      reps=10,      clustdist=1,      subs=5,      nclust=50,      ICC=.3,    L1effect=.2,      L2effect=0,      L12effect=.5)
%mixedmod(      condition=      4      ,      reps=10,      clustdist=1,    subs=10,      nclust=50,      ICC=.1,      L1effect=0,      L2effect=0,      L12effect=.8)
%mixedmod(      condition=      5      ,      reps=10,      clustdist=1,    subs=10,      nclust=50,      ICC=.3,      L1effect=0,      L2effect=.2,      L12effect=0)
%mixedmod(      condition=      6      ,      reps=10,      clustdist=2,      subs=5,      nclust=25,      ICC=.1,     L1effect=0,      L2effect=.2,      L12effect=.2)
%mixedmod(      condition=      7      ,      reps=10,      clustdist=2,     subs=5,      nclust=25,      ICC=.3,     L1effect=0,      L2effect=.2,      L12effect=.5)
%mixedmod(      condition=      8      ,      reps=10,      clustdist=2,     subs=10,      nclust=50,      ICC=.3,     L1effect=0,      L2effect=.2,      L12effect=.8)
%mixedmod(      condition=      9      ,      reps=10,      clustdist=2,     subs=10,      nclust=25,      ICC=.3,     L1effect=0,      L2effect=.2,      L12effect=.8)
%mixedmod(      condition=      10  ,         reps=10,      clustdist=2,      subs=5,      nclust=50,      ICC=.3,     L1effect=0,      L2effect=.2,      L12effect=.8)


/************************************************************************************************************************************/


%macro mixedanalyze(start, stop, type);

proc printto log="C:\SAS_MC\Logs_1\Acondition_&start._to_&stop..log";
run;

%do i=&start %to &stop;
/* Pull in marked data files and make a working file from them, renaming them with their respective condition no. */
data Cond_&i;
infile "C:\SAS_MC\Data_1\Condition_&i..txt";
input reps y repno x1 x2 x3 clustdist clust subs nclust indiv ICC L1effect L2effect L12effect b10 b01 b11 err indvarnum indvarden intvar
 SD;
run;

%end;
data analyze;
set Cond_&start - Cond_&stop;

/*sort 'analyze' by each replication and each cluster */
proc sort data= analyze;
by clustdist subs nclust ICC L1effect L2effect L12effect repno clust indiv;
run;

*Get the mean x1 for each clust;
*output these to a new SAS data set called clust_means;
proc means data=analyze noprint;
by clustdist subs nclust ICC L1effect L2effect L12effect repno clust;
var x1;
output out=clust_means mean(x1)=clustmean;
run;


*Open the clust_means data set and keep only the relevant variables;
data clust_means (keep=clustdist subs nclust ICC L1effect L2effect L12effect repno clust clustmean);
set clust_means;
run;

proc means data=analyze noprint;
var x2;
*adding BY statement to perform calculation for each repetition;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
output out=grand_mean mean(x2)=grandmean;
run;

data grand_mean (keep = clustdist subs nclust ICC L1effect L2effect L12effect repno grandmean mergevar);
set grand_mean;
mergevar=1;
run;

*to add grand_mean to everyone;
data analyze;
set analyze;
mergevar=1;
run;

data analyze2 (drop=mergevar);
merge analyze grand_mean;
by clustdist subs nclust ICC L1effect L2effect L12effect repno mergevar;
run;

*Merge the clust_means data set with analyze2 and create centered vars;

data analyze3;
merge analyze2  clust_means;
by clustdist subs nclust ICC L1effect L2effect L12effect repno clust;
x1_CWC=x1-clustmean;
x2_GMC=x2-grandmean;
x3_Cen=x1_CWC*x2_GMC;
run;

/*MYOUTPUT only */

%if &type=1 %then %do;

libname NFixed "C:\SAS_MC\Test\NFixed";
libname NCov "C:\SAS_MC\Test\NCov";
libname Fixed "C:\SAS_MC\Test\Fixed";
libname Random "C:\SAS_MC\Test\Random";
libname Converge "C:\SAS_MC\Test\Converge";
libname Interate "C:\SAS_MC\Test\IntHist";
libname NObs "C:\SAS_MC\Test\NObs";
libname FitStats "C:\SAS_MC\Test\FitStats";
libname Type3 "C:\SAS_MC\Test\Type3";
libname CovParms "C:\SAS_MC\Test\CovParms";
libname Clust "C:\SAS_MC\Test\Clust";


/* Null model */
ods select none;
ods output SolutionF = NFixed.NFixed_&start._to_&stop.;
ods output CovParms = NCov.Cov_&start._to_&stop.;
proc mixed data=analyze3 covtest method=ML;
class clust;
model y= /s ddfm=satterthwaite;
random int /s type=un subject=clust g;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

/* Full model */

ods output SolutionF = Fixed.Fixed_&start._to_&stop.;
ods output SolutionR = Random.Random_&start._to_&stop.;
ods output ConvergenceStatus = Converge.Convergence_&start._to_&stop.;
ods output IterHistory = Interate.Iteration_&start._to_&stop.;
ods output NObs = NObs.NoObservations_&start._to_&stop.;
ods output FitStatistics = FitStats.FitStat_&start._to_&stop.;
ods output Tests3 = Type3.Type3_&start._to_&stop.;
ods output CovParms = CovParms.Cov_&start._to_&stop.;
ods output ClassLevels = Clust.Clust_&start._to_&stop.;
proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
random int /subject = clust g s type = un;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
ods output close;
ods select all;
quit;

/* Regression model */

ods output SolutionF = Fixed.RFixed_&start._to_&stop.;
*ods output ConvergenceStatus = Converge.RConvergence_&start._to_&stop.;
*ods output IterHistory = Interate.RIteration_&start._to_&stop.;
ods output NObs = NObs.RNoObservations_&start._to_&stop.;
ods output FitStatistics = FitStats.RFitStat_&start._to_&stop.;
ods output Tests3 = Type3.RType3_&start._to_&stop.;
ods output CovParms = CovParms.RCov_&start._to_&stop.;
ods output ClassLevels = Clust.RClust_&start._to_&stop.;
proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
ods output close;
ods select all;
quit;

%end;

/* ODS output only ALL MODELS */

%if &type=2 %then %do;

proc mixed data=analyze3 covtest method=ML;
class clust;
model y= /s ddfm=satterthwaite;
random int /s type=un subject=clust g;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
random int /subject = clust g s type = un;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
quit;


%end;

/* ODS output for Mixed Model only */
%if &type=3 %then %do;
proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
random int /subject = clust g s type = un;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
quit;
%end;

/* ODS output for Mixed and Regression models */

%if &type=4 %then %do;
proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
random int /subject = clust g s type = un;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
quit;
%end;

/*Proc GENMOD REG and MIXED comparison; ODS */

%if &type=5 %then %do;

proc mixed data = analyze3 covtest method=ML;
class clust;
model y = x1_CWC x2_GMC x3_Cen /s corrb ddfm=sat;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc genmod data = analyze3;
model y = x1_CWC x2_GMC x3_Cen;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc reg data=analyze3;
model y = x1_CWC x2_GMC x3_Cen;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;
quit;


%end;

%mend mixedanalyze;

%mixedanalyze (start=1, stop=5, type=1)
%mixedanalyze (start=6, stop=10, type=1)

/********************************************************************************************/
/*Stout Dissertation Results Syntax */


/* Results file */

%macro mixedresults(batch, type);

proc printto log="C:\SAS_MC\Test\Batch_&batch..log";
run;

%if &type=0 %then %do;
/*Transpose data: Create columns from rows for the CovParm estimates*/

proc transpose data="C:\SAS_MC\Test\NCov\cov_&batch." out=Results1_0 name=Estimate_column;
var estimate;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
id CovParm;
run;

/*Transpose data: Create columns from rows for the Estimates of SE from mixed and regression estimates*/
proc transpose data="C:\SAS_MC\Test\Fixed\fixed_&batch." out=Results1_1 name=Estimate_column;
var stderr;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
id effect;
run;

data Results1_1;
set Results1_1;
IntSE_MLM=Intercept;
x1SE_MLM=x1_CWC;
x2SE_MLM=x2_GMC;
x3SE_MLM=x3_Cen;
drop Intercept x1_CWC x2_GMC x3_Cen Estimate_column _Label_;
run;

/**/

proc transpose data="C:\SAS_MC\Test\Fixed\rfixed_&batch" out=Results1_2 name=Estimate_column;
var stderr;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
id effect;
run;

data Results1_2;
set Results1_2;
IntSE_OLS=Intercept;
x1SE_OLS=x1_CWC;
x2SE_OLS=x2_GMC;
x3SE_OLS=x3_Cen;
drop Intercept x1_CWC x2_GMC x3_Cen Estimate_column _Label_;
run;

/*Transpose data: Create columns from rows for MLM p-values of x1 and x3; Type I error*/
proc transpose data="C:\SAS_MC\Test\Fixed\fixed_&batch." out=Results1_3 name=Estimate_column;
var probt;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
id effect;
run;

data Results1_3;
set Results1_3;
x1M_pVal=x1_CWC;
x3M_pVal=x3_Cen;
drop Intercept x1_CWC x2_GMC x3_Cen Estimate_column _Label_;
run;

data Results1_3;
set Results1_3;
if L1Effect=0 and x1M_pVal <.05 then
 do;
 x1M_Typ1=1;
/* x1M_Typ1Cnt+1; */
 end;
else
 do;
 x1M_Typ1=0;
 end;
if L12Effect=0 and x3M_pVal <.05 then
 do;
 x3M_Typ1=1;
/* x3M_Typ1Cnt+1; */
 end;
else
 do;
 x3M_Typ1=0;
 end;
/*Type II */
if L1Effect>0 and x1M_pVal >.05 then
 do;
 x1M_Typ2=1;
 end;
else
 do;
 x1M_Typ2=0;
 end;
if L12Effect>0 and x3M_pVal >.05 then
 do;
 x3M_Typ2=1;
 end;
else
 do;
 x3M_Typ2=0;
 end;
drop x1M_pVal x3M_pVal;
run;


/*Transpose data: Create columns from rows for SRS p-values of x1 and x3; Type I error*/
proc transpose data="C:\SAS_MC\Test\Fixed\rfixed_&batch." out=Results1_4 name=Estimate_column;
var probt;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
id effect;
run;

data Results1_4;
set Results1_4;
x1R_pVal=x1_CWC;
x3R_pVal=x3_Cen;
drop Intercept x1_CWC x2_GMC x3_Cen Estimate_column _Label_;
run;

data Results1_4;
set Results1_4;
if L1Effect=0 and x1R_pVal <.05 then
 do;
 x1R_Typ1=1;
/* x1R_Typ1Cnt+1;*/
 end;
else
 do;
 x1R_Typ1=0;
 end;
if L12Effect=0 and x3R_pVal <.05 then
 do;
 x3R_Typ1=1;
/* x3R_Typ1Cnt+1;*/
 end;
else
 do;
 x3R_Typ1=0;
 end;
/*Type II */
if L1Effect>0 and x1R_pVal >.05 then
 do;
 x1R_Typ2=1;
 end;
else
 do;
 x1R_Typ2=0;
 end;
if L12Effect>0 and x3R_pVal >.05 then
 do;
 x3R_Typ2=1;
 end;
else
 do;
 x3R_Typ2=0;
 end;
drop x1R_pVal x3R_pVal;
run;

/*Merge the transposed files, create effect variances from SEs*/
data Results2_0;
merge Results1_0 Results1_1 Results1_2 Results1_3 Results1_4;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
IntVar_MLM=IntSE_MLM**2;
x1Var_MLM=x1SE_MLM**2;
x2Var_MLM=x2SE_MLM**2;
x3Var_MLM=x3SE_MLM**2;
IntVar_OLS=IntSE_OLS**2;
x1Var_OLS=x1SE_OLS**2;
x2Var_OLS=x2SE_OLS**2;
x3Var_OLS=x3SE_OLS**2;
run;


/*Calculate ratio of mixed var to regression var for each effect*/
/*Calculate the Total Variance, Observed ICC, and resulting Observed Design Effect in each replication*/
data Results3;
set Results2_0;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
TotVar=UN_1_1_+Residual;
ObsICC=UN_1_1_/TotVar;
ObsDE=1+(subs-1)*ObsICC;
Int_DE=IntVar_MLM/IntVar_OLS;
x1_DE=x1Var_MLM/x1Var_OLS;
x2_DE=x2Var_MLM/x2Var_OLS;
x3_DE=x3Var_MLM/x2Var_OLS;
drop Estimate_column;
run;

/*Cleanup*/
data Results3;
set Results3;
drop IntSE_MLM IntSE_OLS x1SE_MLM x1SE_OLS x2SE_MLM x2SE_OLS x3SE_MLM x3SE_OLS IntVar_MLM IntVar_OLS x1Var_MLM x1Var_OLS
x2Var_MLM x2Var_OLS x3Var_MLM x3Var_OLS UN_1_1_ Residual TotVar;
run;


/*Create Omnibus dataset with parameters for all datasets in all conditions, for ANOVA */
data Final_out;
set Results3;
file "C:\SAS_MC\Test\FinalOut\Final_Output.txt" mod;
put clustdist subs nclust ICC L1effect L2effect L12effect repno ObsICC ObsDE Int_DE x1_DE x2_DE x3_DE x1M_Typ1 x1R_Typ1 x3M_Typ1 x3R_Typ1;
run;


/*Average each CONDITION, export to internal/external SAS files and Excel*/
/*Calculate average observed ICC/DE in each condition*/
/*Creates internal SAS results tables for averages */
proc means data=Results3 noprint;
by clustdist subs nclust ICC L1effect L2effect L12effect;
var ObsICC ObsDE Int_DE x1_DE x2_DE x3_DE;
output out=batch_&batch._means
mean(ObsICC)=AvgObsICC mean(ObsDE)=AvgModelDE
mean(Int_DE)=Avg_Int_DE
mean(x1_DE)=Avg_x1_DE
mean(x2_DE)=Avg_x2_DE
mean(x3_DE)=Avg_x3_DE
mean(x1R_Typ1)=TypeI_x1R
mean(x1M_Typ1)=TypeI_x1M
mean(x3R_Typ1)=TypeI_x3R
mean(x3M_Typ1)=TypeI_x3M
mean(x1R_Typ2)=TypeII_x1R
mean(x1M_Typ2)=TypeII_x1M
mean(x3R_Typ2)=TypeII_x3R
mean(x3M_Typ2)=TypeII_x3M;

/*Creates external SAS results tables for averages */
output out="C:\SAS_MC\Test\MergedResults\SAS\Batch_&Batch._results" mean(ObsICC)=AvgObsICC mean(ObsDE)=AvgModelDE
mean(Int_DE)=Avg_Int_DE
mean(x1_DE)=Avg_x1_DE
mean(x2_DE)=Avg_x2_DE
mean(x3_DE)=Avg_x3_DE
mean(x1R_Typ1)=TypeI_x1R
mean(x1M_Typ1)=TypeI_x1M
mean(x3R_Typ1)=TypeI_x3R
mean(x3M_Typ1)=TypeI_x3M
mean(x1R_Typ2)=TypeII_x1R
mean(x1M_Typ2)=TypeII_x1M
mean(x3R_Typ2)=TypeII_x3R
mean(x3M_Typ2)=TypeII_x3M;
run;

data batch_&batch._means;
set batch_&batch._means;
drop _TYPE_ _FREQ_;
run;

/*Creates external Excel results tables for averages */

proc export data=batch_&batch._means outfile= "C:\SAS_MC\Test\MergedResults\Excel\Batch_&Batch..xlsx" dbms=xlsx replace;
run;

%end;

%mend;

%mixedresults (batch= 1_to_5,type=0)
%mixedresults (batch= 6_to_10,type=0)

/*
%mixedresults (batch= 1_to_256,type=0)
%mixedresults (batch= 257_to_512,type=0)
%mixedresults (batch= 513_to_768,type=0)
%mixedresults (batch= 769_to_1024,type=0)
%mixedresults (batch= 1025_to_1280,type=0)
%mixedresults (batch= 1281_to_1536,type=0)
%mixedresults (batch= 1537_to_1792,type=0)
%mixedresults (batch= 1793_to_2048,type=0)
%mixedresults (batch= 2049_to_2304,type=0)
%mixedresults (batch= 2305_to_2560,type=0)
%mixedresults (batch= 2561_to_2816,type=0)
%mixedresults (batch= 2817_to_3072,type=0)
%mixedresults (batch= 3073_to_3328,type=0)
%mixedresults (batch= 3329_to_3584,type=0)
%mixedresults (batch= 3585_to_3840,type=0)
%mixedresults (batch= 3841_to_4096,type=0)
%mixedresults (batch= 4097_to_4352,type=0)
%mixedresults (batch= 4353_to_4608,type=0)
%mixedresults (batch= 4609_to_4864,type=0)
%mixedresults (batch= 4865_to_5120,type=0)
%mixedresults (batch= 5121_to_5376,type=0)
%mixedresults (batch= 5377_to_5632,type=0)
%mixedresults (batch= 5633_to_5888,type=0)
%mixedresults (batch= 5889_to_6144,type=0)
*/

/* Merge internal SAS files by Distribution and m */
/*
data Uniform_5;
merge Batch_1_to_256_means Batch_257_to_512_means Batch_513_to_768_means Batch_769_to_1024_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

data Uniform_10;
merge Batch_1025_to_1280_means Batch_1281_to_1536_means Batch_1537_to_1792_means Batch_1793_to_2048_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

data Uniform_20;
merge Batch_2049_to_2304_means Batch_2305_to_2560_means Batch_2561_to_2816_means Batch_2817_to_3072_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

data Binomial_5;
merge Batch_3073_to_3328_means Batch_3329_to_3584_means Batch_3585_to_3840_means Batch_3841_to_4096_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

data Binomial_10;
merge Batch_4097_to_4352_means Batch_4353_to_4608_means Batch_4609_to_4864_means Batch_4865_to_5120_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

data Binomial_20;
merge Batch_5121_to_5376_means Batch_5377_to_5632_means Batch_5633_to_5888_means Batch_5889_to_6144_means;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;


data FINAL;
merge Uniform_5 Uniform_10 Uniform_20 Binomial_5 Binomial_10 Binomial_20;
by clustdist subs nclust ICC L1effect L2effect L12effect;
run;

proc export data=Uniform_5 outfile= "C:\SAS_MC\MergedResults\Excel\Uniform_05.xlsx" dbms=xlsx replace;
run;
proc export data=Uniform_10 outfile= "C:\SAS_MC\MergedResults\Excel\Uniform_10.xlsx" dbms=xlsx replace;
run;
proc export data=Uniform_20 outfile= "C:\SAS_MC\MergedResults\Excel\Uniform_20.xlsx" dbms=xlsx replace;
run;
proc export data=Binomial_5 outfile= "C:\SAS_MC\MergedResults\Excel\Binomial_05.xlsx" dbms=xlsx replace;
run;
proc export data=Binomial_10 outfile= "C:\SAS_MC\MergedResults\Excel\Binomial_10.xlsx" dbms=xlsx replace;
run;
proc export data=Binomial_20 outfile= "C:\SAS_MC\MergedResults\Excel\Binomial_20.xlsx" dbms=xlsx replace;
run;

*/

/* Pull in Omnibus and make a working file for it */

proc printto log = "C:\SAS_MC\Logs_1\ANOVA.log";
run;

data ANOVA;
infile "C:\SAS_MC\Test\FinalOut\Final_Output.txt";
input clustdist subs nclust ICC L1effect L2effect L12effect repno ObsICC ObsDE Int_DE x1_DE x2_DE x3_DE x1M_Typ1 x1R_Typ1 x3M_Typ1 x3R_Typ1;
run;

/*sort 'ANOVA' by each replication and each cluster */
proc sort data= ANOVA;
by clustdist subs nclust ICC L1effect L2effect L12effect repno;
run;

proc glm data=ANOVA;
class clustdist subs nclust ICC L1effect L2effect L12effect;
model ObsICC = clustdist subs nclust ICC L1effect L2effect L12effect/ ss1 ss2 ss3 ss4 effectsize alpha=.05;
run;

proc glm data=ANOVA;
class clustdist subs nclust ICC L1effect L2effect L12effect;
model ObsICC = clustdist subs nclust ICC L1effect L2effect L12effect
L12effect*ICC / ss1 ss2 ss3 ss4 effectsize alpha=.05;
run;

proc anova data=ANOVA;
class clustdist subs nclust ICC L1effect L2effect L12effect;
model ObsDE = clustdist subs nclust ICC L1effect L2effect L12effect / nouni;
run;

proc genmod data=ANOVA;
class clustdist subs nclust ICC L1effect L2effect L12effect;
model ObsICC = clustdist subs nclust ICC L1effect L2effect L12effect;
run;
