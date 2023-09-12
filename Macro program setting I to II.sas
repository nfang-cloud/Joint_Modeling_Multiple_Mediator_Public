***********************************************************************************
*************************Setting I***************************************************;

%macro jointmodel1(Data=);

%do ii=1 %to 200;
title "Iteration &ii.";

PROC IMPORT OUT= WORK.sim2 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data quant1_r;
infile "&InD.\quant_r.txt";
input qr1_min qr1_2 qr1_4 qr1_6 qr1_max aa;
run;

data quant2_r;
infile "&InD.\quant_r.txt";
input qr2_min qr2_2 qr2_4 qr2_6 qr2_max aa;
run;

data quant_d;
infile "&InD.\quant_d.txt";
input qd_min qd_2 qd_4 qd_6 qd_max aa;
run;

data four;
set sim2;
aa=1;
run;

proc sort data=four;
by id stoptime;
run;

data four2;
merge four quant1_r quant2_r quant_d;
by aa;
event_1=0;
event_2=0;
IF event_Tp=1 then event_1=1;
IF event_Tp=2 then event_2=1;
run;

* Calculate the number of recurrent events as a mediator for death;
data four3;
set four2;
by id;
retain last_stop nevent nevent_1 add_1 nevent_2 add_2;
if first.id then do;
		nevent=0;
        start=0;
		stop=stoptime;
		last_stop=stoptime;
        nevent_1=0;
        add_1=event_1;
        nevent_2=0;
        add_2=event_2;

end;
else do;
		nevent=nevent+1;
		start=last_stop;
		stop=stoptime;
		last_stop=stoptime;
        nevent_1=nevent_1 + 1*add_1;
		add_1=event_1;
		nevent_2=nevent_2 + 1*add_2;
		add_2=event_2;

end;

drop add_1 add_2;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quant1_r {5} qr1_min qr1_2 qr1_4 qr1_6 qr1_max;
array quant2_r {5} qr2_min qr2_2 qr2_4 qr2_6 qr2_max;
array quant_d {5} qd_min qd_2 qd_4 qd_6 qd_max;

array dur1_r {4} dur1_r1-dur1_r4;
array dur2_r {4} dur2_r1-dur2_r4;
array dur_d {4} dur_d1-dur_d4;


array event1_r {4} event1_r1-event1_r4;
array event2_r {4} event2_r1-event2_r4;
array event_d {4} event_d1-event_d4;


array median1_r {4} median1_r1-median1_r4;
array median2_r {4} median2_r1-median2_r4;
array median_d {4} median_d1-median_d4;

do i=1 to 4;
    dur1_r{i}=0;
	dur2_r{i}=0;
	event1_r{i}=0;
	event2_r{i}=0;
	dur_d{i}=0;
	
	event_d{i}=0;
	median1_r{i}=0;
	median2_r{i}=0;
	median_d{i}=0;
end;

* For recurrent event type 1;
if event=1 and event_Tp=1 then do;
	do i=2 to 5;
		if stoptime<=quant1_r{i} then do;
			event1_r{i-1}=1;
			i=5;
		end;
	end;
end;
* For recurrent event type 2;
if event=1 and event_Tp=2 then do;
	do i=2 to 5;
		if stoptime<=quant2_r{i} then do;
			event2_r{i-1}=1;
			i=5;
		end;
	end;
end;


if event in (0,2) then do;
/* If death or censored observation */
	do i=2 to 5;
		if stoptime<=quant1_r{i} then do;
			dur1_r{i-1}=stoptime-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1} + dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
			i=5;
		end;
		else do;
			dur1_r{i-1}=quant1_r{i}-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
		end;
	end;

do i=2 to 5;
		if stoptime<=quant2_r{i} then do;
			dur2_r{i-1}=stoptime-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
			i=5;
		end;
		else do;
			dur2_r{i-1}=quant2_r{i}-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
		end;
	end;

	do i=2 to 5;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=5;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;
* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;

data five2;
set five;
array quant {5} qd_min qd_2 qd_4 qd_6 qd_max;
array dur {4} dur1-dur4;
last_start=start;
do i=1 to 4;
	dur{i}=0;
    
end;

do i=2 to 5;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
           
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
        
		i=5;
	end;
end;
run;	


proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log1_r1=-0.43 log1_r2=-0.69 log1_r3=-0.80 log1_r4=-1.0
      log2_r1=-0.69 log2_r2=-1.0 log2_r3=-1.2 log2_r4=-1.4
	  log_h1=-3.4 log_h2=-2.5 log_h3=-2.3 log_h4=-0.5
       
	  betax1=0.2 betax2=0.15
	  betaz1=0.35 betaz2=0.4
      deltax1=0.7 deltax2=0.84
	  etax=0.15 etaz=0.35
	  delta1=0.7 etam1=0.18 etam2=0.15
	  ;

array dur {4} dur1-dur4;
array baseh {4} log_h1-log_h4;

base_haz1_r=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2 + exp(log1_r3) * event1_r3
            + exp(log1_r4) * event1_r4;
base_haz2_r=exp(log2_r1) * event2_r1 + exp(log2_r2) * event2_r2 + exp(log2_r3) * event2_r3
           + exp(log2_r4) * event2_r4;
cum_base_haz1_r=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2 + exp(log1_r3) * dur1_r3
               + exp(log1_r4) * dur1_r4;
cum_base_haz2_r=exp(log2_r1) * dur2_r1 + exp(log2_r2) * dur2_r2 + exp(log2_r3) * dur2_r3
                + exp(log2_r4) * dur2_r4;

/* different beta**/
/* also different coefficient for vi, only keep one vi*/
mu1_1= betaz1 * X1 + betax1 * X2 + deltax1*vi;
mu2_1= betaz2 * X1 + betax2 * X2 + deltax2*vi;
loglik1_1=-exp(mu1_1) * cum_base_haz1_r;
loglik2_1=-exp(mu2_1) * cum_base_haz2_r;

* Note here the mediator: nevent changes over time, so we need to calculate it for each interval;
* Each record is an interval (start, stop) for nevent;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};

sum2=0;
do k=1 to 4;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam1 * nevent_1 + etam2 * nevent_2);
  
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;;/* for death event */

loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3 + exp(log_h4) * event_d4 ;
	mu4= etaz * X1 + etax * X2 + etam1 * nevent_1  + etam2 * nevent_2 + delta1 * vi;	
end;

if (event=1 and event_Tp=1) then loglik= log(base_haz1_r) + mu1_1 +loglik0; /*log likelihood for recurrent lung event */
if (event=1 and event_Tp=2) then loglik= log(base_haz2_r) + mu2_1 +loglik0; /*log likelihood for recurrent other event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1_1 + loglik2_1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1_1 + loglik2_1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  1) subject=id;
ods output ParameterEstimates=est1&ii. FitStatistics=fit1&ii. CorrMatParmEst=corr1&ii. CovMatParmEst=cov1&ii.; 
run;

%end;

%mend;

%jointmodel1(Data=sim1_multi);

%macro outcov1();
%do ii=1 %to 200;
data _null_;
	set cov1&ii;
	file "&out.\covImulti&ii..txt";
	put parameter log1_r1 log1_r2 log1_r3 log1_r4
      log2_r1 log2_r2 log2_r3 log2_r4
	  log_h1 log_h2 log_h3 log_h4
       
	  betax1 betax2
	  betaz1 betaz2
      deltax1 deltax2
	  etax etaz
	  delta1 etam1 etam2
	  ;
	run;
%end;
%mend;


%macro outest1();
%do ii=1 %to 200;
data _null_;
	set est1&ii;
	file "&out.\estImulti&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outest1();
%outcov1();



%macro jointmodel2(Data=);

%do ii=1 %to 200;
title "Iteration &ii.";

PROC IMPORT OUT= WORK.sim2 
            DATAFILE= "&InD.\&Data.&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data quant1_r;
infile "&InD.\quant_r.txt";
input qr1_min qr1_2 qr1_4 qr1_6 qr1_max aa;
run;

data quant2_r;
infile "&InD.\quant_r.txt";
input qr2_min qr2_2 qr2_4 qr2_6 qr2_max aa;
run;

data quant_d;
infile "&InD.\quant_d.txt";
input qd_min qd_2 qd_4 qd_6 qd_max aa;
run;

data four;
set sim2;
aa=1;
run;

proc sort data=four;
by id stoptime;
run;

data four2;
merge four quant1_r quant2_r quant_d;
by aa;
event_1=0;
event_2=0;
IF event_Tp=1 then event_1=1;
IF event_Tp=2 then event_2=1;
run;

* Calculate the number of recurrent events as a mediator for death;
data four3;
set four2;
by id;
retain last_stop nevent nevent_1 add_1 nevent_2 add_2;
if first.id then do;
		nevent=0;
        start=0;
		stop=stoptime;
		last_stop=stoptime;
        nevent_1=0;
        add_1=event_1;
        nevent_2=0;
        add_2=event_2;

end;
else do;
		nevent=nevent+1;
		start=last_stop;
		stop=stoptime;
		last_stop=stoptime;
        nevent_1=nevent_1 + 1*add_1;
		add_1=event_1;
		nevent_2=nevent_2 + 1*add_2;
		add_2=event_2;

end;

drop add_1 add_2;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quant1_r {5} qr1_min qr1_2 qr1_4 qr1_6 qr1_max;
array quant2_r {5} qr2_min qr2_2 qr2_4 qr2_6 qr2_max;
array quant_d {5} qd_min qd_2 qd_4 qd_6 qd_max;

array dur1_r {4} dur1_r1-dur1_r4;
array dur2_r {4} dur2_r1-dur2_r4;
array dur_d {4} dur_d1-dur_d4;


array event1_r {4} event1_r1-event1_r4;
array event2_r {4} event2_r1-event2_r4;
array event_d {4} event_d1-event_d4;


array median1_r {4} median1_r1-median1_r4;
array median2_r {4} median2_r1-median2_r4;
array median_d {4} median_d1-median_d4;

do i=1 to 4;
    dur1_r{i}=0;
	dur2_r{i}=0;
	event1_r{i}=0;
	event2_r{i}=0;
	dur_d{i}=0;
	
	event_d{i}=0;
	median1_r{i}=0;
	median2_r{i}=0;
	median_d{i}=0;
end;

* For recurrent event type 1;
if event=1 and event_Tp=1 then do;
	do i=2 to 5;
		if stoptime<=quant1_r{i} then do;
			event1_r{i-1}=1;
			i=5;
		end;
	end;
end;
* For recurrent event type 2;
if event=1 and event_Tp=2 then do;
	do i=2 to 5;
		if stoptime<=quant2_r{i} then do;
			event2_r{i-1}=1;
			i=5;
		end;
	end;
end;


if event in (0,2) then do;
/* If death or censored observation */
	do i=2 to 5;
		if stoptime<=quant1_r{i} then do;
			dur1_r{i-1}=stoptime-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1} + dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
			i=5;
		end;
		else do;
			dur1_r{i-1}=quant1_r{i}-quant1_r{i-1};
			median1_r{i-1}=quant1_r{i-1}+dur1_r{i-1}/2; /* Get the median of each interval */
			lastmed1_r=median1_r{i-1};
		end;
	end;

do i=2 to 5;
		if stoptime<=quant2_r{i} then do;
			dur2_r{i-1}=stoptime-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
			i=5;
		end;
		else do;
			dur2_r{i-1}=quant2_r{i}-quant2_r{i-1};
			median2_r{i-1}=quant2_r{i-1}+dur2_r{i-1}/2; /* Get the median of each interval */
			lastmed2_r=median2_r{i-1};
		end;
	end;

	do i=2 to 5;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
			i=5;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
			median_d{i-1}=quant_d{i-1}+dur_d{i-1}/2; /* Get the median of each interval */
			lastmed_d=median_d{i-1};
		end;
	end;
end;

run;
* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;

data five2;
set five;
array quant {5} qd_min qd_2 qd_4 qd_6 qd_max;
array dur {4} dur1-dur4;
last_start=start;
do i=1 to 4;
	dur{i}=0;
    
end;

do i=2 to 5;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
           
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
        
		i=5;
	end;
end;
run;	


proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms log1_r1=-0.43 log1_r2=-0.69 log1_r3=-0.80 log1_r4=-1.0
      log2_r1=-0.69 log2_r2=-1.0 log2_r3=-1.2 log2_r4=-1.4
	  log_h1=-3.4 log_h2=-2.5 log_h3=-2.3 log_h4=-0.5
       
	  betax1=0.2 betax2=0.15
	  betaz1=0.35 betaz2=0.4
      deltax1=0.7 deltax2=0.84
	  etax=0.15 etaz=0.35
	  delta1=0.7 etam1=0.18 etam2=0.15
	  ;

array dur {4} dur1-dur4;
array baseh {4} log_h1-log_h4;

base_haz1_r=exp(log1_r1) * event1_r1 + exp(log1_r2) * event1_r2 + exp(log1_r3) * event1_r3
            + exp(log1_r4) * event1_r4;
base_haz2_r=exp(log2_r1) * event2_r1 + exp(log2_r2) * event2_r2 + exp(log2_r3) * event2_r3
           + exp(log2_r4) * event2_r4;
cum_base_haz1_r=exp(log1_r1) * dur1_r1 + exp(log1_r2) * dur1_r2 + exp(log1_r3) * dur1_r3
               + exp(log1_r4) * dur1_r4;
cum_base_haz2_r=exp(log2_r1) * dur2_r1 + exp(log2_r2) * dur2_r2 + exp(log2_r3) * dur2_r3
                + exp(log2_r4) * dur2_r4;

/* different beta**/
/* also different coefficient for vi, only keep one vi*/
mu1_1= betaz1 * X1 + betax1 * X2 + deltax1*vi;
mu2_1= betaz2 * X1 + betax2 * X2 + deltax2*vi;
loglik1_1=-exp(mu1_1) * cum_base_haz1_r;
loglik2_1=-exp(mu2_1) * cum_base_haz2_r;

* Note here the mediator: nevent changes over time, so we need to calculate it for each interval;
* Each record is an interval (start, stop) for nevent;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};

sum2=0;
do k=1 to 4;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etam1 * nevent_1 + etam2 * nevent_2);
  
end;

mu2= etaz * X1 + etax * X2 + delta1 * vi;;/* for death event */

loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3 + exp(log_h4) * event_d4 ;
	mu4= etaz * X1 + etax * X2 + etam1 * nevent_1  + etam2 * nevent_2 + delta1 * vi;	
end;

if (event=1 and event_Tp=1) then loglik= log(base_haz1_r) + mu1_1 +loglik0; /*log likelihood for recurrent lung event */
if (event=1 and event_Tp=2) then loglik= log(base_haz2_r) + mu2_1 +loglik0; /*log likelihood for recurrent other event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglik1_1 + loglik2_1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglik1_1 + loglik2_1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  1) subject=id;
ods output ParameterEstimates=est2&ii. FitStatistics=fit2&ii. CorrMatParmEst=corr2&ii. CovMatParmEst=cov2&ii.; 
run;

%end;

%mend;
%jointmodel2(Data=sim2_multi);


%macro outcov2();
%do ii=1 %to 200;
data _null_;
	set cov2&ii;
	file "&out.\covIImulti&ii..txt";
	put parameter log1_r1 log1_r2 log1_r3 log1_r4
      log2_r1 log2_r2 log2_r3 log2_r4
	  log_h1 log_h2 log_h3 log_h4
       
	  betax1 betax2
	  betaz1 betaz2
      deltax1 deltax2
	  etax etaz
	  delta1 etam1 etam2
	  ;
	run;
%end;
%mend;


%macro outest2();
%do ii=1 %to 200;
data _null_;
	set est2&ii;
	file "&out.\estIImulti&ii..txt";
	put Parameter Estimate;
	run;
%end;
%mend;

%outest2();
%outcov2();


