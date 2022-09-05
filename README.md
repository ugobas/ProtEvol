# Prot_evol
Author  : Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM) ubastolla@cbm.csic.es

This repository has been substituted by https://github.com/ugobas/Prot_evol Please go there for the last version of the program Prot_evol that includes structure constrained substitution models of protein evolution.

The program Prot_evol computes site-specific amino acid substitution models where different protein sites evolve independently under selection on the thermodynamic stability of the native state against unfolding and misfolding.
The program implements two models: (1) The mean-field (MF) in which each amino acid site evolves in the mean-field produced by the other sites; (2) The wild-type (WT) model that considers the effect of mutations of the wild-type sequence at all possible positions.
Additionally, the program also allows simulating stability constrained protein evolution, without applying the approximation that protein sites evolve independently of each other.

Citations:
1) Jimenez MJ, Arenas M, Bastolla U (2018) Substitution rates predicted by stability-constrained models of protein evolution are not consistent with empirical data. Mol Biol Evol. 35: 743-755
2) Arenas M, Weber CC, Liberles DA, Bastolla U (2017) ProtASR: An Evolutionary Framework for Ancestral Protein Reconstruction with Selection on Folding Stability. Syst Biol. 66:1054-1064. 
3) Arenas M, Sánchez-Cobos A, Bastolla U (2015) Maximum-Likelihood Phylogenetic Inference with Selection on Protein Folding Stability. Mol Biol Evol. 32:2195-207.
4) Minning J, Porto M, Bastolla U (2013) Detecting selection for negative design in proteins through an improved model of the misfolded state. Proteins 81:1102-12.


==================================
GENERAL DESCRIPTION
==================================

The program Prot_Evol performs two kinds of computations.

(1) It computes site-specific and global amino acid frequencies F and
    exchangeability matrices E, which jointly allow computing global and
    site-specific amino acid substitution matrices Q=EF that can be used
    for phylogenetic inference.
    The amino acid frequencies are obtained imposing a global constraint on
    the average stability of the native state of the protein against both
    unfolding and misfolding, with the approximation that frequencies at
    different sites are independent and they have minimal differences with
    respect a global background distribution.
    Two approximations are used: mean-field (MF) in which the effect of a
    mutation on stability is estimated with respect to the mean field created
    by the other sites, and wild-type (WT) in which the effect of a mutation on
    stability is computed with respect to the wild-type sequence.
    The WT approach is expected to be more accurate at small evolutionary
    distances, and the MF approach is expected to be more accurate at large
    evolutionary distances.
    Both approaches make buried sites more hydrophobic and exposed sites more
    polar, while intermediate sites are most tolerant and more variable.
    The exchangeability matrices are computed under the Halpern-Bruno formula
    that expresses the fixation probabilities computed by Kimura.


(2) It simulates protein evolution subject to the constraint of selection
    on the folding stability of the native state against both unfolding and
    misfolding. It implements two selection models:
    
        (a) Neutral (stabilities above threshold have fitness 1 and below
	    threshold have fitness 0);
        (b) Based on the fitness function f=1/(1+exp(DG/T)) and on the fixation
	    probability computed by Kimura,
	    P_fix=(1-f(WT)/f(mut))/(1-(f(WT)/f(mut))^N),
            where N is the effective population size;

==================================
COMPILE AND RUN:
==================================

To compile and run the program execute the following commands:

    unzip Prot_evol.zip
    make
    ./Prot_evol # Without parameter file, it will print an help screen
    ./Prot_evol Input_Prot_evol.in

==================================
INPUT:
==================================

The parameter file Input_Prot_evol.in included in the package contains the
input parameters that can be modified and their explanation.

IMPORTANT: 
 - Modify the line PDB=... with the local path to your target PDB file
 - Modify the line FILE_STR=... with the local path to the file 
   structures.in, which is included in the package.

==================================
DETAILED DESCRIPTION
==================================

    The commands mentioned in this section can be modified by editing 
the attached input file Input_Prot_evol.in

A) MEAN_FIELD (MF) and WILD_TYPE (WT) models with independent sites.
    
    This computation is performed if MEANFIELD=1. 
    
    L independent mean-field amino acid distributions (P^MF,i)_a are 
computed by minimizing the Kullback-Leibler divergence with respect to 
the background distribution (P^mut)_a for fixed value of the average 
stability, ave(Delta G). This constraint it imposed through a Lagrange 
multiplier Lambda, whose value is optimized by imposing that the sequences 
in the PDB plus the optional input protein family has maximum likelihood
with respect to the model (option OPT_LAMBDA=1) or can be input through
the file (OPT_LAMBDA=0, LAMBDA=whatever).

    1) The BACKGROUND distribution can be computed in three different ways.

        1a) As the amino acid distribution observed across all protein 
    sites (command GET_FREQ=2, 19 degrees of freedom). If one of the 
    amino acids is not present, this distribution is combined with the 
    one at point (1b).
     
        1b) Fitted from a codon based mutation model, whose parameters 
    are the nucleotide frequencies pi_A, pi_T, pi_G, pi_C, the transition/
    tranversion ratio TT_RATIO and the enhancement of the mutation rate at
    CpG dinucleotides kCpG. These parameters are fitted imposing that the
    background distribution is as similar as possible to the observed one (1a)
    (GET_FREQ=1, 5 degrees of freedom). 
     
        1c) The parameters of the mutation model can be input manually 
    from the input file (GET_FREQ=0, 0 degrees of freedom).

    2) The program outputs the EXCHANGEABILITY matrices E^i_ab that allow
    computing the substitution rate between amino acids a and b at site i as
    (Q^i)_ab=(E^i)_ab*(P^MF,i)_b
    
        We impose that (E^i)_ab are symmetric, so detailed balance is 
    fulfilled and the stationary distributions coincide with the mean-field
    distributions. (E^i)_ab and (P^MF,i)_b fulfill the Halpern-Bruno formula
    that relates the fixation probability with the stationary distributions.
    Four kinds of exchangeability models are implemented. Three are based on
    an empirical substitution model that can be chosen with the command
    MATRIX=... (available options are WAG and JTT).

        2a) EXCHANGE=MUT The global exchangeability matrix is computed from
    the same mutation process that generates the background distribution.
        
        2b) EXCHANGE=EXCH The global exchangeability matrix is equal to the 
    empirical one (JTT or WAG, specified with MATRIX=WAG). The two former
    choices yield poor results.
    
        2c) EXCHANGE=RATE The exchangeability matrix is computed 
    imposing that the average amino acid substitution rate across 
    all sites is equal to the one generated by the empirical 
    rates, which is Q_ab=(E^emp)_ab*(f^emp)_b
        
        2d) EXCHANGE=FLUX (default) The exchangeability matrix is 
    computed imposing that the average flux of amino acids across all 
    sites is equal to the one generated by the empirical rates, which is 
    F_ab=(E^emp)_ab*(f^emp)_a*(f^emp)_b.

B) Simulations of protein evolution with selection on folding stability.
    
    If only the mean-field distributions are of interest, this 
computation can be switched off with the command TMAX=0.

    At each time step, one mutated sequence is generated, its folding 
stability against unfolding and misfolding DG_mut is estimated with the 
model of Minning, Porto and Bastolla (Proteins 2013, 81:1102-112), and 
the mutation is fixed or rejected with one of three possible SELECTION 
models:

    3a) Neutral: fixation if DG_mut < DG_thr, rejected otherwise.
    
    DG_thr= 0.95*DG_PDB
    Commands NEUTRAL=1 MEANFIELD=0

    3b) Fixation probability of the Moran model with very low mutation 
    rate, P_fix=(1-f_wt/f_mut)/[1-(f_wt/f_mut)^N], where f_wt and f_mut 
    are the wild type and mutant fitness, respectively, 
    f=exp(-DG/T)/[1+exp(-DG/T)], and N is the effective population size,
    which is input with the command NPOP=... (default is 100).
    
    Commands NEUTRAL=0 MEANFIELD=0

==================================
INPUT FILE
==================================
    
    A detailled documentation of the input parameters needed will be found on 
Input_Prot_evol.in file included in this package.
It includes three sections:
A) Input files (PDB file, chain identifier)
B) Thermodynamic model. Modifying the configurational entropy with respect to
unfolding SU1 and with respect to misfolding SC1 one type of stability will
prevail on the other one.
C) Selection model (recommended not to touch)
D) Evolution simulations. To activate them, set TMAX to a large number.
E) Mutation model to build the background distribution and the 
   exchangeability matrices
F) Output control.

==================================
OUTPUT FILES:
==================================
   
   There are several output files that can be classified in two 
different groups, depending on the meaning of their content.
    

A) Evolutionary model files. Two sets of files are produced, one for
   the wild-type model (_WT) and one for the mean-field model (_MF).

    a1) <PDB>_<parameters>_profiles.txt 
    Mean-field amino acid distributions and entropy (one site per line)

    a2) <PDB>_<parameters>_exchangeability_sites_<EXCH_MODEL>.txt
    Exchangeability matrix (one for each site)
    IMPORTANT: Due to their large size, site-specific exchangeability matrices
    can be disabled by setting PRINT_E=0

    a3) <PDB>_<parameters>_AA_profile_global.txt
    Average amino acid distribution across all sites
    
    a4) <PDB>_<parameters>_exchangeability_glob_<EXCH_MODEL>.txt
    Exchangeability matrix averaged over all sites.
    
    a5) <PDB>_<parameters>_summary.dat
    Performances of the mean-field model evaluated with the PDB 
    sequence. This is a tabular file with the following fields:
    
    Protein length; Fitted nucleotide frequencies; PDB code (last column);
    The following measures are given for various types of models: 
    fitted Lambda; KL divergence from background distribution;
    log-likelihood per site of the PDB sequence with respect to the model; 
    average hydrophobicity h;
    freezing temperature of the misfolded ensemble, Tf; 
    average free energy difference between native and misfolded state, DG; 
    	    Models:
	    mut: background distribution at all sites.
            nat: mean-field model constraining only the native energy.
            all: mean-field model constraining DG 
            wt:  sequence in the PDB.

    a6) <PDB>_<parameters>_likelihood.dat
    For each tested value of Lambda, DG, likelihood, hydrophobicity h,
    freezing temperature Tf and convergenge of the nat and all models
    are reported (the mut model coincides with Lambda=0). Conv=0 means 
    that the iterations did not converge.
    
    a7) <PDB>_<parameters>_rate_profile.dat 
    For each protein site we report:
    (1) the site identifier in the PDB,
    (2) the amino acid at the site in the PDB,
    (3) the secondary structure,
    (4) the number of native contacts,
    (5) the entropy of the predicted amino acid distribution,
    (6) the predicted average hydrophobicity.
    (7-9) the predicted substitution rate for the mut, exc and flux
    exchangeability models.

    Pearson correlation coefficients between these variables are 
    reported at the end of the file.
    
B) Simulations of evolution:

    b1) Contact_statistics_<L>.dat
    Statistics of alternative contact matrices of length L
    
    b2) <PDB>_<parameters>_stab.dat
    For each amino acid substitution, it reports (1) the mutated site,
    (2) the mutated amino acid, (3) the native energy, (4) DG, (5) fitness, 
    (6) number of synonymous substitutions since the previous aa substitution
    (7) number of attempted mutations.
    This file can be used to reconstruct the protein sequences generated
    by the simulated evolution.

    b3) <PDB>_<parameters>_ave.dat
    Every it_print (internal parameter) substitutions, it prints average
    and standard error of: (1,2) fitness, (3,4) native energy, 
    (5,6) DeltaG, (7) difference between amino acid entropy of mutation and 
    selection, (8) non_synonymos/synonymous subst. rate, (9) acceptance rate.

    b4) <PDB>_<parameters>_final.dat
    Same information as in stab file, but printed at the end of the 
    simulation.

