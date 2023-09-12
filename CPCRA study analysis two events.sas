*Import the dataset from lib ddiddc;
data ddiddc;
set ddiddc.ddiddc;
trt=randgrp-1;
gender=gender-1;
hemobl=hemobl-12;
cd4bl=cd4bl/100;
id=seq;

run;


data one;
set ddiddc;
array event_all{20} MAC PCP1 PCP2 PCP3 TB OMYC HIST 
CANE1 CANE2 CANE3 CMV WAST 
KSV ADC CRYC LYMP PML TOXO CRYS HZ1;
array time_all{20} T2MAC T2PCP1 T2PCP2 T2PCP3 T2TB T2OMYC T2HIST 
T2CANE1 T2CANE2 T2CANE3 T2CMV T2WAST 
T2KSV T2ADC T2CRYC T2LYMP T2PML T2TOXO T2CRYS T2HZ1;
aa=1;

do i=1 to 7;
    
	event_L=event_all{i};
	event_O=0;
	stoptime=time_all{i}/30;
	IF event_all{i}=. THEN event_O=.;
	output;
end;

do i=8 to 20;
    
	event_O=event_all{i};
	event_L=0;
	stoptime=time_all{i}/30;
	IF event_all{i}=. THEN event_L=.;
	output;
end;
label event_O="Other OI event" 
event_L="Lung related OI event";


run;


data Lung_two;
set one;
if event_L=1;
aa=1;
run;
proc sort;
by id stoptime;
run;

proc univariate data=Lung_two noprint;
var stoptime; 
output out=quantL_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qrL; 
run;
data quantL_r;
set quantL_r;
aa=1;
run;


data Other_two;
set one;
if event_O=1;
aa=1;
run;

proc sort;
by id stoptime;
run;

proc univariate data=Other_two noprint;
var stoptime; 
output out=quantO_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qrO; 
run;
data quantO_r;
set quantO_r;
aa=1;
run;


* The dataset for death event;
data three;
set ddiddc;
by seq;
if first.seq;
stoptime=t2death/30;
event=death*2;		* Set event=2 for death;
aa=1;
run;

* Get the quantiels for death time;
proc univariate data=three noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
where event=2;
run;
data quant_d;
set quant_d;
aa=1;
run;


data four;
set Lung_two other_two three;
IF event_L=1 THEN event=3; *lung related recurrence;
IF event_O=1 THEN event=4;* not lung related recurrence;
run;

proc sort data=four;
by id stoptime;
run;
* Merge data with the quantiles;

data four2;
merge four quantL_r quantO_r quant_d;
by aa;

run;


* Calculate the number of recurrent events as a mediator for death;
data four3;
set four2;
by id;
retain last_stop nevent nevent_L add_L nevent_O add_O;
if first.id then do;
		nevent=0;
        start=0;
		stop=stoptime;
		last_stop=stoptime;
        nevent_L=0;
        add_L=event_L;
        nevent_O=0;
        add_O=event_O;

end;
else do;
		nevent=nevent+1;
		start=last_stop;
		stop=stoptime;
		last_stop=stoptime;
        nevent_L=nevent_L + 1*add_L;
		add_L=event_L;
		nevent_O=nevent_O + 1*add_O;
		add_O=event_O;

end;

drop add_L add_O;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four3;
array quantL_r {6} qrL0 qrL20 qrL40 qrL60 qrL80 qrL100;
array quantO_r {6} qrO0 qrO20 qrO40 qrO60 qrO80 qrO100;
array quant_d {6} qd0 qd20 qd40 qd60 qd80 qd100;

array durL_r {5} durL_r1-durL_r5;
array durO_r {5} durO_r1-durO_r5;
array dur_d {5} dur_d1-dur_d5;

array eventL_r {5} eventL_r1-eventL_r5;
array eventO_r {5} eventO_r1-eventO_r5;
array event_d {5} event_d1-event_d5;

do i=1 to 5;
    durL_r{i}=0;
	durO_r{i}=0;
	eventL_r{i}=0;
	eventO_r{i}=0;
	dur_d{i}=0;
	event_d{i}=0;

end;

* For lung recurrent event;
if event=3 then do;
	do i=2 to 6;
		if stoptime<=quantL_r{i} then do;
			eventL_r{i-1}=1;
			i=6;
		end;
	end;
end;
if event=4 then do;
	do i=2 to 6;
		if stoptime<=quantO_r{i} then do;
			eventO_r{i-1}=1;
			i=6;
		end;
	end;
end;


if event in (0,2) then do;

	do i=2 to 6;
		if stoptime<=quantL_r{i} then do;
			durL_r{i-1}=stoptime-quantL_r{i-1};
			
			i=6;
		end;
		else do;
			durL_r{i-1}=quantL_r{i}-quantL_r{i-1};
			
		end;
	end;

do i=2 to 6;
		if stoptime<=quantO_r{i} then do;
			durO_r{i-1}=stoptime-quantO_r{i-1};
			
			i=6;
		end;
		else do;
			durO_r{i-1}=quantO_r{i}-quantO_r{i-1};
		
		end;
	end;


	do i=2 to 6;
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=stoptime-quant_d{i-1};
	
			i=6;
		end;
		else do;
			dur_d{i-1}=quant_d{i}-quant_d{i-1};
		
		end;
	end;
end;

run;
* For each of mediator "nevent", we need to calculate the duration in each quantile interval of death time;


data five2;
set five;
array quant {6} qd0 qd20 qd40 qd60 qd80 qd100;
array dur {5} dur1-dur5;
last_start=start;
do i=1 to 5;
	dur{i}=0;
    
end;

do i=2 to 6;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
           
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
        
		i=6;
	end;
end;
run;	


* Model I in Biometrics paper;


proc nlmixed data=five2 qpoints=5 MAXITER=5000 corr cov;

parms logL_r1=-3.5 logL_r2=-2.4 logL_r3=-1.9 logL_r4=-1.5 logL_r5=-1.1
      
      logO_r1=-3.2 logO_r2=-2.1 logO_r3=-1.6 logO_r4=-1.2 logO_r5=-0.74
      
	  log_h1=-3.7 log_h2=-2.5 log_h3=-1.9 log_h4=-1.5 log_h5=-0.94
       
	  betaL_x1=0.5 betaL_x2=-0.4 betaL_x3=.12 betaL_x4=-.06 betaL_x5=-0.56 betaL_z=0.19 /*lung*/
	  betaO_x1=0.34 betaO_x2=-.75 betaO_x3=-.18 betaO_x4=-.14 betaO_x5=-0.43 betaO_z=-0.15 /*other*/
      deltaO=0.7304
	  etax1=0.63 etax2=-0.04 etax3=-0.18 etax4=-.31 etax5=-0.73 etaz=-0.28 
	  etamL=0.01 etamO=-0.124 delta1=1.589 
      log_varc=-0.809
;

array dur {5} dur1-dur5;

array baseh {5} log_h1-log_h5;


base_hazL_r=exp(logL_r1) * eventL_r1 + exp(logL_r2) * eventL_r2 + exp(logL_r3) * eventL_r3
            + exp(logL_r4) * eventL_r4 + exp(logL_r5) * eventL_r5 
            ;

base_hazO_r=exp(logO_r1) * eventO_r1 + exp(logO_r2) * eventO_r2 + exp(logO_r3) * eventO_r3
           + exp(logO_r4) * eventO_r4 + exp(logO_r5) * eventO_r5 
;


cum_base_hazL_r=exp(logL_r1) * durL_r1 + exp(logL_r2) * durL_r2 + exp(logL_r3) * durL_r3
               + exp(logL_r4) * durL_r4 + exp(logL_r5) * durL_r5 

;

cum_base_hazO_r=exp(logO_r1) * durO_r1 + exp(logO_r2) * durO_r2 + exp(logO_r3) * durO_r3
                + exp(logO_r4) * durO_r4 + exp(logO_r5) * durO_r5 

;


/* different beta**/
/* also different coefficient for vi, only keep one vi*/
muL_1= betaL_x1 * prevoi + betaL_x2 * gender + betaL_z*trt  + betaL_x3*stratum + betaL_x4*hemobl 
       + betaL_x5* cd4bl + vi;			/* for recurrent event */

muO_1= betaO_x1 * prevoi + betaO_x2 * gender + betaO_z*trt  + betaO_x3*stratum + betaO_x4*hemobl 
       + betaO_x5* cd4bl + deltaO*vi;			/* for recurrent event */

loglikL1=-exp(muL_1) * cum_base_hazL_r;
loglikO1=-exp(muO_1) * cum_base_hazO_r;

* Note here the mediator: nevent changes over time, so we need to calculate it for each interval;
* Each record is an interval (start, stop) for nevent;
* dur{k} is the duration of nevent in each of the death quantiles, where nevent is a constant in each dur{k};

sum2=0;
do k=1 to 5;
	/* cumulative baseline hazard for time dependent measure */
	sum2=sum2 + exp(baseh{k}) * dur{k} * exp(etamL * nevent_L + etamO * nevent_O);
  
end;

mu2= etax1 * prevoi + etax2* gender + etaz*trt + etax3*stratum + etax4 * hemobl + etax5 * cd4bl + delta1 * vi;/* for death event */

loglik0=-exp(mu2) * sum2;

if event=2 then do; 				/* for death event */
	base_haz_d=exp(log_h1) * event_d1 + exp(log_h2) * event_d2 + exp(log_h3) * event_d3 + exp(log_h4) * event_d4 
               + exp(log_h5) * event_d5 
 ;
	mu4= etax1 * prevoi + etax2* gender + etaz*trt + etax3*stratum + etax4 * hemobl 
         + etax5 * cd4bl + etamL * nevent_L  + etamO * nevent_O + delta1 * vi;	
end;

if event=3 then loglik= log(base_hazL_r) + muL_1 +loglik0; /*log likelihood for recurrent lung event */
if event=4 then loglik= log(base_hazO_r) + muO_1 +loglik0; /*log likelihood for recurrent other event */
if event=2 then loglik=loglik0 + log(base_haz_d) + mu4 + loglikL1 + loglikO1;	/*log likelihood for death */
if event=0 then loglik=loglik0 + loglikL1 + loglikO1;							/*log likelihood for censoring */

model id ~ general(loglik);
random vi ~ normal(0,  exp(log_varc)) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1 CorrMatParmEst=corr1 CovMatParmEst=cov1; 
run;




data _null_;
	set cov1;
	file "&out.\covI_OL.txt";
	put parameter 
		logL_r1 logL_r2 logL_r3 logL_r4 logL_r5 
        logO_r1 logO_r2 logO_r3 logO_r4 logO_r5 
		log_h1 log_h2 log_h3 log_h4 log_h5
		betaL_x1 betaL_x2 betaL_z betaL_x3 betaL_x4 betaL_x5
	    betaO_x1 betaO_x2 betaO_z betaO_x3 betaO_x4 betaO_x5
        deltaO
	    etax1 etax2 etaz etax3 etax4 etax5
	    delta1 etamL etamO
		log_varc;

	run;

data _null_;
	set est1;
	file "&out.\estI_OL.txt";
	put Parameter Estimate;
	run;


ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\estimation OL.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=est1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "Estimation with different quant time for lung and other (5 jump points)";

	run;

ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;


ODS RESULTS OFF;
ODS LISTING CLOSE;
ODS EXCEL file="&out.\covariance OL.xlsx"
    options (start_at="B1" tab_color="red" absolute_row_height="15" embedded_titles="yes" embedded_footnotes="yes");
	proc report data=cov1 nowindows style(Header)=[background=white foreground=black];
	ods Excel;
	column _all_;

	title j=L h=2.5 "covariance with different quant time for lung and other (5 jump points)";

	run;

ODS Excel CLOSE;
ODS LISTING;
ODS RESULTS ON;

