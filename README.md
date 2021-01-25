# Customised-DRAGON2_beta
Customised version of the DRAGON2_beta code (https://github.com/cosmicrays/DRAGON2-Beta_version). It contains important updates in the diffusion routines.

The main change is that the diffusion coefficient is now calculated as function of rigidity, instead of as function of kinetic energy per nucleon, as in the standard code. This allows the implementation of spectral features depending of rigidity, like breaks in the diffusion coefficient. 

Three different functions are predefined for rigidity breaks in the diffusion coefficient, that are called with the <BreakDiff/> command belowe the Diffusion type instruction in the xml input file:

  -A break at high rigidity value is used by default. The parametrisation of the diffusion coefficient implemented is the same as in Genolini et. al, PRL 119, 241101 (2017).
  The break position (rho_b), spectral index after the break (delta_H) and smoothing parameter (s_H) can be set simply using:
    <delta_H value="" />
    <rho_H value="" />
    <s_H value="" />
  where the value must be given inside the "". 

  - A low energy break is also available, by adding to the xml input file the <LowE_Break/> instruction below the <BreakDiff/> one. The parameters can be set using:
    <delta_L value="" />
    <rho_L value="" />
    <s_L value="" />

  - A doubly broken diffusion coefficient is also implemented, following Weinrich et al, A&A 639, A131 (2020). It needs the instruction <Double_Break/> below the <BreakDiff/>.         The parameters to be set are delta_L, delta_H, rho_L, rho_H, s_L, s_H, specified as before. 

These routines can be called in the case of 2D model for the Galaxy structure, for the 3D model only the high energy break is implemented at the moment. These routines have been used, among other works, in Luque et. al, arXiv:2101.01547 (https://arxiv.org/abs/2101.01547), Luque et al., 2020 (https://iopscience.iop.org/article/10.1088/1742-6596/1690/1/012010) and Luque et al., IN PREPARISON. 


In addition, there are other routines implemented for an inhomogeneous diffusion coefficient (i.e. Diff. coefficient dependent on position in the Galaxy). These are based either in two_zones models (Diff. coefficient calculated different in the halo and the galactic disk), based on Fornieri et al, arXiv:2011.09197 (https://arxiv.org/abs/2011.09197) or a diffusion coefficient whose spectral index changes linearly with the radius of the Galaxy (as explored in Pothast et al, 2018 - https://iopscience.iop.org/article/10.1088/1475-7516/2018/10/045). For the former, example routines are left in the example folder. 
To implement the latter, the parametrisation of delta is delta(R) = delta_A*r + delta_B, where r is the galactic radius. 

To specify these parameters, the input xml file must contain, below the Diffusion type instruction:     
    <VariableDelta />
    <!TwoZone />
    <deltaB value="" />
    <deltaA value="" />


 ### Please, if you use this version of the code, cite as Pedro de la Torre Luque. (2021, January 25). tospines/Customised-DRAGON2_beta: First release of the customized DRAGON2_beta code (Version v0.1). Zenodo. http://doi.org/10.5281/zenodo.4461732
 
For any doubts about the implementation of these routines, please, email to pedro.delatorreluque@fysik.su.se
