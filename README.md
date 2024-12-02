# psilocybin_heirarchy
A repo for the analysis of the effects of psilocybin on the heirarchical structure of fMRI data in collaboration with Deco lab

## Brief data description
fMRI data in 28 participants scanned before and up to 4 times after receivied psilocybin
Subjects were scanned on two scanners with different TRs, so analyses are performed on each individually.

## Analyses

Model-free analyses

  1: Normalised modularity using BCT toolbox

  2: INSIDEOUT heirarchicality measure

  3: Assymetry of the forward time-lagged functional connectivity matrix


Model-based analyses

  Model 1 : Generative effective connectivity trained on FC matrix only

  1: Modularity

  2: Trophic score

  Model 2 : Generative effective connectivity trained on FC matrix and INSIDEOUT matrix

  1: Modularity

  2: Trophic score

## Statistics
Each psilocybin scan is accompanied by a plasma drug level measure (psi_conc). The association between each metric and plasma drug level is estimated using a random slope model. i.e., 
fitlme(df, 'metric ~ psi_conc + Age + Sex + (1|ID)', 'FitMethod', 'ML');

## Hypotheses
Broadly, psilocybin level will be associated with decreased heirarchicality in the brain.
Specifically: 

  1: Decreased modularity

  2: Decreased global INSIDEOUT score

  3: Decreased assymetry

  4: Narrowing of range of trophic scores
