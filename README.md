# SFU: Surface-Free Utility-Based Design for Dose Optimization in Cancer Drug-Combination Trials

R codes to implement surface-free utility-based design for dose optimization in drug-combination trials.

# Description

Precision oncology has demonstrated the potential of drug combinations in effectively enhancing anti-tumor efficiency and controlling disease progression. Nonetheless, dose optimization in early-phase drug-combination trials presents various challenges and is considerably more complex than single-agent dose optimization. To address this, we propose a surface-free design for exploring the optimal doses of combination therapy within the phase I-II framework. Rather than relying on parametric models to define the shape of toxicity and efficacy surfaces, our approach centers on characterizing dose-toxicity and dose-efficacy relationships between adjacent dose combinations using surface-free models. The proposed design encompasses a run-in phase, facilitating a swift exploration of the dose space, followed by a main phase where the dose-finding rule relies on the proposed surface-free model. Through extensive simulation studies, we have thoroughly examined the operating characteristics of this innovative design. Our method exhibits desirable operating characteristics across a wide range of dose-toxicity and dose-efficacy relationships.

# Functions

The repository includes two functions:

- sfu_4by4.R: The R code includes the function ```sfu_4by4()``` to obtain the operating characteristics of the SFU design for 4*4 scenarios by simulating trials.
  
  ```rscript
  sfu_4by4(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
           assessment.window, accrual.rate, u11, u00)
  ```

- sfu_tite_4by4.R: The R code includes the function `sfu_tite_4by4()` to obtain the operating characteristics of the SFU design for 4*4 scenarios with late-onset toxicities and efficacies by simulating trials.
  
  ```rscript
  sfu_tite_4by4(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
                assessment.window, accrual.rate, u11, u00, maxpen)
  ```

- sfu_3by5.R: The R code includes the function `sfu_3by5()` to obtain the operating characteristics of the SFU design for 3*5scenarios by simulating trials.
  
  ```rscript
  sfu_3by5(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
           assessment.window, accrual.rate, u11, u00)
  ```

- sfu_tite_3by4.R: The R code includes the function `sfu_tite_3by5()` to obtain the operating characteristics of the SFU design for 3*5 scenarios with late-onset toxicities and efficacies by simulating trials.
  
  ```rscript
  sfu_tite_3by5(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
                assessment.window, accrual.rate, u11, u00, maxpen)
  ```

# Inputs

- `ntrial`: The total number of trials to be simulated.

- `p.true.tox`: A matrix containing the true toxicity probabilities of the investigational dose combinations.

- `p.true.eff`: A matrix containing the true efficacy probabilities of the investigational dose combinations.

- `target.tox`: The upper limit of toxicity rate, e.g., `target.tox <- 0.35`.

- `target.eff`: The lower limit of efficacy rate, e.g., `target.eff <- 0.2`.

- `samplesize`: The maximum sample size.

- `cohortsize`: The cohort size.

- `assessment.window`: The assessment window for toxicity and efficacy.

- `accrual.rate`: The accrual rate of patients (patients per month).

- `u11`: The specified utility value for (toxicity=1, efficacy=1).

- `u00`: The specified utility value for (toxicity=0, efficacy=0).

- `maxpen`: The maximum proportion of patients who have pending DLT or efficacy outcomes at the current dose.

# Outputs

- `sfu_4by4()` , `sfu_tite_4by4()`, `sfu_3by5()`, and `sfu_tite_3by5()`will return the operating characteristics of the SFU design as a data frame, including:
  
  ```
   (1) the true utility of each dose combination (u.true);  
   (2) selection percentage at each dose combination (selpercent);  
   (2) the number of patients treated at each dose combination (n.pts);
   (2) the number of toxicities treated at each dose combination (n.tox);  
   (3) the number of efficacies treated at each dose combination (n.eff);
   (5) the average number of patients (totaln);    
   (4) the average number of toxicities (totaltox);  
   (4) the average number of efficacies (totaleff);  
   (6) the average trial duration (duration).
  ```

# Example

We consider using the `sfu_4by4()` as an illustration.

- Suppose the upper limit of toxicity rate is 35%. The lower limit of efficacy rate is 20%. A maximum of 51 patients will be recruited in cohorts of 3. Suppose 4 dose levels of agent A and 4 of agent B are considered. We can use the following code to simulate the scenario 1 in Table 1.

```rscript
p.true.tox <- matrix(c(0.12,0.16,0.18,0.20,  0.14,0.18,0.20,0.40,  0.18,0.20,0.42,0.55,  0.20,0.40,0.52,0.60),nrow = 4, byrow = TRUE)
p.true.eff <- matrix(c(0.20,0.25,0.30,0.50,  0.25,0.30,0.50,0.60,  0.30,0.50,0.60,0.66,  0.40,0.60,0.65,0.70),nrow = 4, byrow = TRUE)

sfu_4by4(rseed=1, ntrial=10, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, assessment.window, accrual.rate, u11, u00)

  -----------------------output------------------------

$u.true
     [,1] [,2] [,3] [,4]
[1,] 47.2 48.6 50.8 62.0
[2,] 49.4 50.8 62.0 60.0
[3,] 50.8 62.0 59.2 57.6
[4,] 56.0 60.0 58.2 58.0

$selpercent
     [,1] [,2] [,3] [,4]
[1,]    0    0    0   20
[2,]    0    0    0   20
[3,]    0   10   20   10
[4,]   10   10    0    0

$n.pts
     [,1] [,2] [,3] [,4]
[1,]  3.3  2.1  2.1  6.6
[2,]  1.2  1.8  3.0  3.9
[3,]  1.2  5.7  6.3  3.6
[4,]  1.5  7.2  1.5  0.0

$n.tox
     [,1] [,2] [,3] [,4]
[1,]  0.4  0.1  0.1  1.1
[2,]  0.1  0.2  0.8  1.4
[3,]  0.3  0.8  3.2  1.7
[4,]  0.1  2.4  0.8  0.0

$n.eff
     [,1] [,2] [,3] [,4]
[1,]  0.4  0.7  0.6  3.6
[2,]  0.2  0.5  1.6  2.2
[3,]  0.5  2.3  4.4  2.8
[4,]  0.6  4.8  0.9  0.0

$totaltox
[1] 13.5

$totaleff
[1] 26.1

$totaln
[1] 51

$duration
[1] 45.22
```

# Authors and Reference

* Jingyi Zhang, Nolan A. Wages, and Ruitao Lin
* Zhang, J., Wages, N. A., and Lin, R. (2023) “SFU: Surface-Free Utility-Based Design for Dose Optimization in Cancer Drug-Combination Trials.” Revision submitted to journal.
