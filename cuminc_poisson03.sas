********************************cuminc_poisson************************;

*Find the cumulative incidence from a model with piecewice constant rates.

Input:

data:		data must contain variables describing the following:
				the case variable: the number of cases.
				the pyrs variable: the number of person years
				the strata variable: indicates what kind of stade the cases account for.
				the time variable
				the categorical covariates
				the continuous covarriates
			
data_cov	a dataset with the values of the covariates for which the cumulative incidences 
			are to be estiamted. This dataset must hold the following variables strata, the 
			categorical covariates, the continuous covarriates. Their must be as many observations as there is 
			failure types. The strata variable must take the values 1,2,3,...., and for each 
			strata value the other variables must take the values that one is interesset in.
case:		the number of cases.
pyrs:		the number of person years.
strata:		the variable in the dataset indicating the deferent types of cases in the analysis.
			This must be numeric and take the values 1,2,3...etc.
categorical:	The categorical covariate from the dataset to be used in the model. 
			If more then 1 categorical covariate, the names schould be seperatet with a space.
continuous:		The continuous covariates from the dataset to be used in the model. 
			If more then 1 continuous covariate, the names schould be seperatet with a space.
time:		The time variable in the dataset, which gives the starting point of every interval.
tau_1:		The value of the beginning og the first time interval.
tau_2:		The value of the end of the last time interval.
step:		the value of the interval between to estiamtes of the cumulative incidences. 
			every endpoint in every subinterval must bee devidable with step.
model:		The model to be analysed in a poisson regression using the genmod procedure. 
			This model curresponds to the right hand side of the equal sign in the genmod 
			proceure.
*Output:

*The Macro return a dataset 'cuminc' containing the following:

*CI01-CI0nstrata:	The estimates of the cumulative incidences for each strata. 
*var01-var0nstrata:	The variance of the estimates of the cumulative incidences for each strata.
*Futher more 'cuminc' containts 
*the time variable, indicating the intervals with length given in the step varible.


*see E:\config\makro\cuminc\eks.sas for an exampel;

%macro cuminc_poisson(data,data_cov,case,pyrs,strata,categorical,continuous,time,tau_1,tau_2,step,model);

*estimate in the model;

proc sort data=&data out=_code;
by &strata &time &categorical &continuous;
run;

data _code;
set _code;
logpyrs=log(&pyrs);
run;

proc genmod data=_code;
class &categorical &time &strata;
model &case = &model/dist=poisson link=log offset=logpyrs noint covb type3;
ods output ParameterEstimates=_parm;
ods output covB=_cov;
run;

data _parm;
set _parm;
if parameter in ('Intercept','Scale') then delete;
run;

data _null_;
set _parm nobs=nobs;
call symput ('nobs',nobs);
run;

%let nobs=%eval(&nobs-1+1);

*make the design matrix;

proc sort data=_code out=_code(keep = &strata &time &categorical) nodupkey;
by &strata &time &categorical;
run;

proc sort data=&data_cov out=_code_t(keep= &strata &continuous);
by &strata;
run;

data _code;
merge _code_t _code;
by &strata;
retain _obs 0;
_obs=_obs+1;
run;

proc glmmod data=_code outdesign=_design;
class  &categorical &time &strata;
model _obs = &model / noint;
ods output ClassLevels=_classlevels;
run;

proc sort data=_design;
by _obs;
run;

proc sort data=_code;
by _obs;
run;

data _design;
merge _design _code;
by _obs;
run;

proc sort data=_design;
by &strata &categorical &continuous;
run;

proc sort data=&data_cov;
by &strata &categorical &continuous;
run;

data _design;
merge _design &data_cov(in=a);
by &strata &categorical &continuous;
if a;
*keep col1-col&nobs;
run;

data _code;
set _design;
retain obs 0;
obs=obs+1;
keep &strata &time obs;
run;

data _design;
set _design;
keep col1-col&nobs;
run;

*transforming the estimates and the covariance matrix;

proc iml;
use _parm;
read all var {df estimate} into est;
close _parm;
use _design;
read all var _num_ into design;
close _design;

df=est[,1];
est=est[,2];
rate=design*est;

do i=1 to nrow(rate);
	rate[i,1]=exp(rate[i,1]);
end;

rd1=j(nrow(design),ncol(design),0);

do i=1 to nrow(design);
	v=rate[i,1];
	rd1[i,]=v*design[i,];
end;
rd=rd1[,1];
do i=2 to ncol(rd1);
	if df[i,1]^=0 then rd=rd||rd1[,i];
end;

use _cov;
read all var _num_ into cov;
close _cov;
cov1=rd*cov*t(rd);

obs=1:nrow(rate);
rate=t(obs)||rate;
varnames={'obs' 'rate'};
create _rate from rate [colname=varnames];
append from rate;
create _cov from cov1;
append from cov1;
quit;
run;

*prepare the rates for calculations.;

proc sort data=_code;
by obs;
run;

proc sort data=_rate;
by obs;
run;

data _rate;
merge _code _rate;
by obs;
run;

data _null_;
set _classlevels;
if class="&strata" then call symput('nstrata',levels);
run;

%let nstrata=%eval(&nstrata+1-1);

proc sort data=_rate out=_prepare(keep = &strata &time rate);
by &strata descending &time ;
run;

data _prepare;
set _prepare;
by &strata;
if first.&strata=1 then gem=&tau_2;
t_n=gem;
gem=&time;
retain gem;
drop gem;
run;

proc sort data=_prepare;
by &time &strata;
run;

data _prepare;
set _prepare;
by &time;
array r{*} rate1-rate&nstrata;
array t{*} time_n1-time_n&nstrata;
array tn{*} time_nr1-time_nr&nstrata (&nstrata*0);
array ts{*} time_s1-time_s&nstrata;
r{&strata}=rate;
t{&strata}=t_n;
tn{&strata}=tn{&strata}+1;
ts{&strata}=&time;
retain rate1-rate&nstrata time_n1-time_n&nstrata time_s1-time_s&nstrata;
if last.&time=1 then do;
	time_n=time_n1;
	do i=1 to &nstrata;
		if time_n>t{i} then time_n=t{i};
	end;
	call symput('nafl',sum(of time_nr1-time_nr&nstrata));
end;
else delete;
drop &strata rate t_n i;
run;

%let nafl=%eval(&nafl+1-1);

%do i=1 %to &nstrata;
	data _null_;
	set _prepare end=end;
	if end=1 then call symput("ntime&i",time_nr&i);
	run;

	%let ntime&i=%eval(&&ntime&i+1-1);
	
	data _time_nr&i (index=(time_nr&i));
	set _prepare;
	by time_nr&i;
	if first.time_nr&i;
	keep &time time_nr&i;
	run;
%end;

%do i=1 %to &nstrata;

	data _survival (index=(&time));
	set _prepare;
	rate=rate&i;
	survrate=sum( of rate1-rate&nstrata);
	t=lag(&time);
	s=lag(time_n);
	su=lag(survrate);
	if &time=&tau_1 then suml=0;
	if &time>&tau_1 then suml=suml+su*(s-t);
	retain suml;
	drop rate1-rate&nstrata s t su;
	run;

*calculate the cumulative incidences;

	data _cuminc;
	set _survival;
	if _n_=1 then do;
		_t=&tau_1;
		cuminc=0;
		sum1=0;
		sum2=0;
		output;
		_t=_t+&step;
	end;

	do while(_t<time_n);

		ci_til=(rate/survrate)*exp(-sum2)*(1-exp(-survrate*(_t-&time)));
		cuminc=sum1+ci_til;
		output;
		_t=_t+&step;
	end;

	sum1=sum1+(rate/survrate)*exp(-sum2)*(1-exp(-survrate*(time_n-&time)));
	sum2=sum2+survrate*(time_n-&time);

	retain sum1 sum2 _t;
	rename _t=&time;
	rename cuminc=CumInc0&i;
	keep _t cuminc;
	run;

*calculate the derivative of the cumulative incidence;

	data _afl;
	set _survival;
	array afl{&nafl} afl1-afl&nafl (&nafl*0);
	array t_st{&nstrata} time_s1-time_s&nstrata;
	array t_nr{&nstrata} time_nr1-time_nr&nstrata;
	array t_n{&nstrata} time_n1-time_n&nstrata;
	array ntime{&nstrata} ntime1-ntime&nstrata (%do n=1 %to &nstrata; &&ntime&n %end;);
	do tid=&time to time_n-&step by &step;
		index=0;
		%do w=1 %to &nstrata %by 1; *iterere over antallet af strata;

			if &w>1 then index=sum(of ntime1-ntime&w)-ntime&w;*opdatere index, således at der tages højde for forskellige antal rater under forskelleg strata;

			do v=1 to t_nr{&w} by 1;*iterater over antallet af rater i strata w;
				if &time=t_st{&w} and v=t_nr{&w} then do;
					afl{index+v}=exp(-suml)*(
								(survrate*(&i=&w)-rate)/(survrate*survrate)*(1-
													exp(-survrate*(tid-&time))) 
								+ rate/survrate*exp(-survrate*(tid-&time))*(tid-&time));
				end;
				if v=t_nr{&w} and t_st{&w}<&time then do;
					i=t_st{&w};
					stop=&time;
					afl_temp=0;
					swmw=t_st{&w};
					&time=i;
					set _survival key=&time /unique;
					do while(i<stop);
						afl_temp=afl_temp+
							exp(-suml)*(-1*(&time-swmw))*(rate/survrate)*(
							1-exp(-survrate*(time_n-&time)))
							+exp(-suml)*((survrate*(&i=&w)-rate)/(survrate*survrate)*(
							1-exp(-survrate*(time_n-&time)))
							+rate/survrate*exp(-survrate*(time_n-&time))*(time_n-&time));
						i=time_n;
						&time=i;
						set _survival key=&time /unique;
					end;	
					afl{index+v}=afl_temp+
						exp(-suml)*(-1*(&time-swmw))*(rate/survrate)*(
						1-exp(-survrate*(tid-&time)))
						+exp(-suml)*((survrate*(&i=&w)-rate)/(survrate*survrate)*(
						1-exp(-survrate*(tid-&time)))
						+rate/survrate*exp(-survrate*(tid-&time))*(tid-&time));
				end;
				if v<t_nr{&w} then do;
					stop=&time;
					afl_temp=0;
					time_nr&w=v;
					set _time_nr&w key=time_nr&w /unique;
					set _survival key=&time /unique;
					swmw=t_st{&w};
					swmw1=t_n{&w};
					do while(&time<swmw1);
						afl_temp=afl_temp+
							exp(-suml)*(-1*(&time-swmw))*rate/survrate*(
							1-exp(-survrate*(time_n-&time)))
							+ exp(-suml)*((survrate*(&i=&w)-rate)/(survrate*survrate)*(
								1-exp(-survrate*(time_n-&time)))
								+rate/survrate*exp(-survrate*(time_n-&time))*(time_n-&time));
						&time=time_n;
						set _survival key=&time /unique;
					end;
					do while(&time<stop);
						afl_temp=afl_temp+
							exp(-suml)*(-1*(swmw1-swmw))*rate/survrate*(
							1-exp(-survrate*(time_n-&time)));
						&time=time_n;
						set _survival key=&time /unique;
					end;
					afl{index+v}=afl_temp+
						exp(-suml)*(-1*(swmw1-swmw))*rate/survrate*(
							1-exp(-survrate*(tid-&time)));
				end;
			end;
		%end;
		output;
	end;
	rename tid=&time;
	keep tid afl1-afl&nafl;
	run;

*calculate the variances of the cumulative incidences;

	proc iml;

	v='afl1';
	if &nafl>1 then do;
		%do n=2 %to &nafl;
			v=v//"afl&n";
		%end;
	end;

	use _afl;
	read all var v into afl;
	read all var {&time} into t;
	close _afl;
	
	use _cov;
	read all var _num_ into cov;
	close _cov;

	m=afl*cov;
	v=j(nrow(t),1,0);
	do i=1 to nrow(t);
		v[i,1]=m[i,]*t(afl[i,]);
	end;

	v1=t||v;
	varnames={"&time" "var0&i"};
	create _var from v1 [colname=varnames];
	append from v1;
	quit;
	
*gather in one dataset;

	data _var;
	set _var;
	&time=round(&time,&step);
	run;

	proc sort data=_var;*indeholder varianserne;
	by &time;
	run;
	/*
	data _var;
	set _var;
	i=i+1;
	retain i 0;
	run;	
	*/

	data _cuminc;
	set _cuminc;
	&time=round(&time,&step);
	run;

	proc sort data=_cuminc;
	by &time;
	run;
	/*
	data _cuminc;
	set _cuminc;
	i=i+1;
	retain i 0;
	run;

	proc sort data=_var;
	by i;
	run;

	proc sort data=_cuminc;
	by i;
	run;
	*/
	data _cuminc;
	merge _cuminc _var;
	by &time;
	*drop i;
	run;

%if &i=1 %then %do;
		data cuminc;
		set _cuminc;
		if CumInc0&i<=0 then do;
			CumInc0&i=0; 
			var0&i=0;
		end;
		run;
	%end;
	%else %do;

		proc sort data=_cuminc;
		by &time;
		run;
		
		proc sort data=cuminc;
		by &time;
		run;

		data cuminc;
		merge cuminc _cuminc;
		by &time;
		if CumInc0&i<=0 then do;
			CumInc0&i=0; 
			var0&i=0;
		end;
		run;
	%end;
%end; *fra 'do i=1 to nstrata' løkken;
%mend;

