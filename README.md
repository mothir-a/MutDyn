
# MutDyn: Mutation Rate Dynamics Analysis Program

## Script Name
**Script:** mutation_rate_dynamics.py

## Development Date
**Date:** 20th Jan 2023

## Author
**Author:** Motohiro Akashi

## About
Mutations play a key role in the evolution of organisms and viruses. Traditional evolutionary simulations have often overlooked the dynamics of mutation rates, focusing instead on traits and adaptations. The disparity mutagenesis model highlights the conservation of duplicated genetic information through unequal mutation rates, theoretically avoiding extinction due to intracellular damage. This model contrasts with the parity model, which assumes equal mutation rates under a constant rate. This program, MutDyn (Mutation Rate Dynamics Analysis Program), extends the investigation by analyzing the evolution of two mutation rates, which can vary between 0 and 1, using a simple simulation model. The simulation imposed no constraints on the mutation rates, emphasizing their inherent functionality.

## Change Log
- README.md was made (26th Jul 2024).
- Fix error caused by ambiguous truth value of a pandas Series (pandas 2.x) (20th Nov 20).

## Development Environment
- **conda version:** 23.5.0
- **conda-build version:** 3.25.0
- **python version:** 3.8.16.final.0
- **jupyter:** v1.0.0
- **pandas:** v1.5.3
- **matplotlib:** v3.5.3
- **numpy:** v1.21.2

## Dependencies
- **python version:** > 3.8
- **pandas:** v1.5.3
- **matplotlib:** v3.5.3
- **numpy:** v1.21.2


## Installation Instructions

If necessary, plealse install the libraries.

i.e. biopython
```zsh
$ conda install conda-forge::pandas
```

## Input Values Requirements
- Before calculation, please input the following values located at Line 581~591.
- initial_mut_rate_1                    : Mutation rate 1
- initial_mut_rate_2                    : Mutation rate 2
- initial_population_size               : Initial population size
- population_max                        : Max population size. Under 5000 is recommended.
- output_mid_data_interval              : Output generation intervals of mid data.
- reputation_num                        : Max generation.
- sharp_fall                            : Generation in which sharp fall occur (Punctuated equilibrium), 'off' or int.
- sharp_fall_val = 'off'                : Sharp fall score, 'off' or int. Lower than 99 is recommended.
- mut_dif_selection                     : Selection by log10 difference of mut_rates, 'off' or int.
- mut_fix                               : Fix mutation rate, 'on' or 'off'.
- fig_dpi                               : Summary fig dpi setting.

## Output Description
- test_YYYY_MMD_HHMM_SS                                  : Directory including result files
  - initial_setting_YYYY_MMD_HHMM_SS.csv                 : Initial setting data.
  - mutation_rate_dynamics_genXXX_survived.csv           : Survived population data at generation XXX.
  - mutation_rate_dynamics_genXXX.csv                    : Total population data at generation XXX.
  - mutation_rate_dynamics_summary_YYYY_MMD_HHMM_SS.csv  : Summary data per generation.
  - summary_df_fig_YYYY_MMD_HHMM_SS.png                  : Figure file. Fitness score, Fidelity difference (FD) and average mutation rate per generation.

## Usage
Input the initial parameters mentioned above.

After that, execute the following code:

```zsh
$ conda activate XXXX           # (optional) activate your environment
$ python mutation_rate_dynamics.py  
```
**Note 1:** Replace `XXXX` with the name of your conda environment if necessary.
**Note 2:** If you execute the code without parameter setting, the calculate will be performed with following values.

```py
                initial_mut_rate_1 = 0.1              # mut rate 1
                initial_mut_rate_2 = 0.1              # mut rate 2
                initial_population_size = 10          # initial population size
                population_max = 50                   # under 5000 is recommended
                output_mid_data_interval = 100        # output mid df per generation
                reputation_num = 1000                 # generation
                sharp_fall = 'off'                    # generation in which sharp fall occur, 'off' or int
                sharp_fall_val = 'off'                # sharp fall score, 'off' or int, lower than 99 is recommended
                mut_dif_selection = 'off'             # selection by log10 difference of mut rates, 'off' or int
                mut_fix = 'off'                       # fix mutation rate, 'on' or 'off'
                fig_dpi = 200                         # summary fig dpi setting
```

## License
- GNU GENERAL PUBLIC LICENSE Version 3.

## Contact Information
- motohiro-akashi[at]st.seikei.ac.jp (M.A.)

**Note:** You can reach out for support or questions related to the script.

## References
- _in submission_
