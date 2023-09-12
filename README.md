# multiple_mediator
**************************************************************************************************************************
*****************************************Programs for the Study*********************************************************

*******************Part I Simulation Setting I to Setting II*************************
1. Generate Datasets for setting I to II. R
    - Used to generate 200 datasets for each setting

2. Macro program setting I to II. sas
   - Estimate the model parameters under settings I to II

3. Summary the fitting results from setting I to II.R
   -Summary the estimation of model parameters under settings I to II
   - Resutls shown in Table S1-S2



*******************Part II Self Defined Function*************************
1. Lambdainv.R
    -Function used to find vector T such that Lambda(T)=s for a vector input s

2. myprod.R
    -Fucntion used to product of two step function

3. mysimRec.R
   -Function using Cinlar's inversion Method to generate Non-homogeneous Poisson process

4. datagenT01_M_multi.R
   - Function used to generate 2 types of recurrent events

5. datagenT01_M_multi5.R
   -Function used to generate 5 types of recurrent events

6. datagen_X.R
   -Function used to generate datasets with simulated X and Z

7. mymed.R
   -Function used to compute the NDE/NIE/NIE_k for simulation settings

8. mymed_T01_2.R
   -Function used to compute the NDE/NIE/NIE_k for two recurrent events

9. mymed_T01_5.R
   -Function used to compute the NDE/NIE/NIE_k for five recurrent events

10. S.R
    - Used to compute the survival functions

 

*******************Part III Estimate NDE/NIE Under Simulation Setting I to II*************************

1. Folder "NDE NIE from Simulation I"
    *SimulationI_parallel.R
              - Function used to estimate NDE/NIE/NIE_k under simulation setting I
    *SimulationI_parallel_boot.R
             -Function of the bootstrap for NDE/NIE/NIE_k under simulation setting I
     *SimulationI_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting I


2. Folder "NDE NIE from Simulation II"
    *SimulationII_parallel.R
              - Function used to estimate NDE/NIE/NIE_k under simulation setting II
    *SimulationII_parallel_boot.R
             -Function of the bootstrap for NDE/NIE/NIE_k under simulation setting II
    *SimulationII_summary.R
             -Function used to summary and estimate bias, SD, Mese and CR for NDE/NIE under simulation setting II





*******************Part IV Real Data Analysis*************************

1. CPCRA study analysis five events.sas
    -Estimates the parameters under simulation setting I with five OIs for CPCRA study

2. CPCRA study analysis two events.sas
    -Estimates the parameters under simulation setting I with two OIs for CPCRA study

3. Folder "CPCRA NDE NIE Plot"

3.1 Subfolder "Five events"
   * NDE_NIE_TRT_5.R 
      -Estimate NDE/NIE/NIE-PCP/NIE-MAC/NIE-CMV/NIE-WAST/NIE-TOXO for treatment

   * NDE_NIE_CI_TRT_5.R
     -Boostrap of NDE/NIE/NIE-PCP/NIE-MAC/NIE-CMV/NIE-WAST/NIE-TOXO for treatment

   * Summary Five Events.R
    -Summary the results and plot the estimates of NDE/NIE/TE/NIE-PCP/NIE-MAC/NIE-CMV/NIE-WAST/NIE-TOXO with bootstraped 95% CI for treatment and CD4

  * NDE_NIE_CD4_5.R 
      -Estimate NDE/NIE/NIE-PCP/NIE-MAC/NIE-CMV/NIE-WAST/NIE-TOXO for CD4

   * NDE_NIE_CI_CD4_5.R
     -Boostrap of NDE/NIE/NIE-PCP/NIE-MAC/NIE-CMV/NIE-WAST/NIE-TOXO for CD4

3.2 Subfolder "Two events"
   * NDE_NIE_TRT_2.R 
      -Estimate NDE/NIE/NIE-Lung/NIE-Other for treatment

   * NDE_NIE_CI_TRT_2.R
     -Boostrap of NDE/NIE/NIE-Lung/NIE-Other for treatment

   * Summary Two Events.R
    -Summary the results and plot the estimates of NDE/NIE/TE/NIE-Lung/NIE-Other with bootstraped 95% CI for treatment and CD4

  * NDE_NIE_CD4_2.R 
      -Estimate NDE/NIE/NIE-Lung/NIE-Other for CD4

   * NDE_NIE_CI_CD4_2.R
     -Boostrap of NDE/NIE/NIE-Lung/NIE-Other for CD4


 
