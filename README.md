# Cost-effectiveness of ledipasvir/sofosbuvir (LS) vs pegylated interferon/ribavirin (PR) for chronic hepatitis C: a Markov model in R

A Markov cohort cost-effectiveness analysis in base R.
It re-implements a model originally built in Excel as coursework, including
deterministic and probabilistic sensitivity analyses.
This is an educational exercise, not a published analysis.

## Model overview

- **Cohort:** 1,000 patients with chronic hepatitis C, starting at age 40
- **Strategies:** PR vs LS (the new regimen)
- **States:** SVR (cured), Moderate disease, Cirrhosis, Liver failure, Death
- **Structure:** treatment success (SVR) sets the initial split of the cohort. Non-responders
  progress **irreversibly** Moderate → Cirrhosis → Liver failure. Death is absorbing, from liver failure or natural causes (age-specific) in any state. SVR patients are cured **for good**, with no re-infection.
- **Time:** 1-year cycles, 70 cycles (ages 41 to 110), lifetime horizon
- **Discounting:** 3.5% per year for **both** costs and outcomes
- **Outcomes:** costs, QALYs, ICER, net monetary benefit (NMB)
- **Willingness-to-pay (WTP):** 20%, 39% (base) and 49% of GDP per capita
  (GDP per capita 2023: 77,350.76 EGP, so 15,470 / 30,167 / 37,902 EGP per QALY)
- **Currency:** EGP
- **Conventions inherited from the Excel model:** cycle 1 is the initial allocation (decision tree); treatment cost is paid once at cycle 0 and not discounted.

## Uncertainty analyses

| Input | One-way (tornado) | Threshold analysis | PSA |
|---|---|---|---|
| SVR rates, transition probabilities, utilities, state costs (12 inputs with a standard error) | mean ± 1.96 × SE | not used | drawn from distributions |
| Discount rate | scenario, 2% to 6% (cost and QALY rates varied together) | not used | fixed |
| LS price | not varied | maximum price solved for | fixed |
| PR price, natural mortality | fixed | fixed | fixed |

**PSA distributions:** beta for probabilities and utilities, gamma for costs, with parameters
derived from each mean and SE by the method of moments. 1,000 Monte Carlo iterations (**`set.seed(123)`**).

**Threshold analysis:** NMB is linear in the LS price, so the maximum price at which NMB ≥ 0
is NMB(price = 0) / cohort size.

## Results

### Base case (per patient, discounted)

| | Cost (EGP) | QALYs |
|---|---|---|
| PR | 33,868 | 11.299 |
| LS | 48,994 | 12.520 |
| Incremental | 15,126 | 1.221 |

**ICER: 12,387 EGP per QALY gained.** Cohort NMB at the 39% GDP threshold: 21,711,933 EGP.



### Deterministic sensitivity analysis

[Summarise: which parameters drive the NMB most (e.g. SVR utility, discount rate, Moderate utility),and whether LS stays cost-effective across all ranges. The NMB is negative only when
the SVR utility is at its lower bound.]

![Tornado diagram](figures/tornado.png)

**Maximum LS price per patient (NMB = 0):** 
The current LS price is 44,000 EGP per patient.

### Probabilistic sensitivity analysis



![Cost-effectiveness plane](figures/ce_plane.png)
![CEAC](figures/ceac.png)

## Limitations

- Utilities are drawn independently, so the SVR utility can fall below the Moderate utility
  (in roughly a quarter of PSA draws), which is clinically implausible. This produces some draws with negative incremental QALYs and keeps the CEAC below 100% at high WTP. The same applies to the lower bound in the one-way analysis. 
- Discount rates for costs and outcomes are varied jointly, not separately.
- Natural mortality and treatment prices are fixed in the PSA.



## Author

Min Kang
