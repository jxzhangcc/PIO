# PIO
Principal Interacting Orbital

Requirement
---
This manual works for the following combinations.
- Python 2.7
- Numpy 1.14.1
- Gaussian 09 & 16
- NBO 6.0 and above are recommended but NBO 3.0 is also compatible (see below for detailed discussions).

For Python 3 users, please see https://github.com/jxzhangcc/PIO-py3.

Tutorial
---
1. Run quantum chemistry calculation
    
    Optmize the structure and calculate the electronic structure using standard quantum chemistry softwares. Gaussian09 is used for following demonstration.
    
    **Input**
    - *FILENAME.gjf*
        
        `
        "# ... pop=nboread
        ...
         $nbo archive file=FILENAME $end"
        `
    
    **Output**
    - *FILENAME.log*
  
        ordinary Gaussian output file
        
    - *FILENAME.FChk*
 
        FormCheck file is required for the current version of PIO analysis, for the purpose of visualizing orbital using GaussView. It can be generated via keyword *formcheck* specified in gjf file, or via running external formchk utility
    
    - *FILENAME.47*
        
        A FILENAME.47 file can be generated in Gaussian calculations via keyword *pop=nboread* and nbo keyword *$NBO archive $end* specified in gjf file. As the built-in NBO package in Gaussian 09 is NBO3.0, which is a quite old version of NBO package, an external NBO analysis with more recent version is suggested. 

2. Run an NBO analysis

    Run a NBO analysis on the system of interest to obtain the NAO basis and corresponding density matrix.
    
    ***NBO 6.0 or 7.0***
    
    For people who have bought individual NBO programs (such as NBO 6.0 or 7.0), modify the 47 file generated in the last step as following and run the NBO program.
    
    **Input**
    
    - *FILENAME.47*
        
        NBO input file, can be generated via running an internal NBO analysis in Gaussian
    
    The following NBO keywords should be specified in the NBO analysis:
    
        $NBO AONAO=W49 FNAO=W49 DMNAO=W49 SKIPBO $END
   
    **Running the NBO analysis**         
    
    Run the NBO analysis using the following command (replace FILENAME with your own file name, and note .47 should not be included):
    
        gennbo FILENAME
    
    **Output**
    - *FILENAME.nbo*
        
        ordinary NBO output file
    - *FILENAME.49*
        
        extra output required for PIO analysis, containing NAO coefficients, density matrix in NAO basis, and Fock matrix in NAO basis if available; only appear if the above mentioned NBO keywords have been specified properly
    
    ***NBO 3.0 built in Gaussian***
    For people who do not have an external NBO progam, the NBO 3.0 built in Gaussian is still compatible to complete PIO analysis. Just modify the Gaussian input file as following:
    
    **Input**
    - *FILENAME.gjf*
        
        `
        \# ... pop=nboread
        ...
         $nbo archive file=FILENAME AONAO=W49 FNAO=W49 DMNAO=W49 SKIPBO $end
        `
    
    The NBO 3.0 built in Gaussian will read in the keywords and do the same as external NBO programs.
    
    **Output**
    - *FILENAME.log*
    - *FILENAME.49*
    
    However, the results with NBO 3.0 should be dealed with cautious because it is known that NBO 3.0 may produce weird results such as counter-intuitive negative charge for transition metals in coordination complexes.
        
    ***For spin-polarized systems***
    No matter you use NBO 3.0 or 6.0, if you want to apply PIO analysis on spin-polarized systems (no matter restricted open-shell or unrestricted calculations), change the NBO keywords from "AONAO=W49 FNAO=W49 DMNAO=W49" to "AONAO=W33 FNAO=W61 DMNAO=W71" and run the following command to produce the required *FILENAME.49* file in correct format.
    
    **Command**
        
        `cat FILENAME.33 FILENAME.61 FILENAME.71 > FILENAME.49`

    ***Note***
    
    A bash script "genpionbos" is attached with the code that can perform all these steps for you by automatically modifying the keywords in .47 file and then generating the .49 file by calling the external NBO program. Make sure relevant paths are correctly set if you use this script. This script works for both open-shell and closed-shell species.

3. Run PIO analysis

    We now run PIO analysis with the following command. Note that FILENAME.49 must exist with the same name as FILENAME.FChk.
    
    **Command**
    
        `python PIO.py FILENAME.FChk`
    
    The program will then request for a fragmentation as input to specify two groups of atoms
    
    **Input**
    
        `1-5,8,13 6-7,9-12`
    
    Two groups of atom IDs should be input here separated by a space. Numbers in each group are separated by a comma. Hyphen is supported for sequential numbers. Atom numbering starts from 1. Complete fragmentation is always recommended (i.e. the specified two groups cover all the atoms present in the system). Incomplete fragmentation will lead to absence of mathematical elegance but is still meaningful if you really want to do it.

    **Output**
    - *FILENAME_pio.txt*
        
        PIO log file, containing basic information of the PIOs of the system subject to the input fragmentation
        
    - *FILENAME_pio.FChk*
        
        Gaussian FormCheck file containing PIOs labeled as in the txt file, could be visualized by GaussView and other compatable orbital visualization softwares
        
    - *FILENAME_pimo.FChk*
        
        A similar Gaussian FormCheck file containing PIMOs whose ordering is same as that of PIOs
        
    - *FILENAME_pio.raw* (discarded in latest PIO version)
        
        numpy-formatted raw data file for debugging or advanced user

Update at Jan 16, 2020
---
Principal Interacting Spin Orbital (PISO) analysis now available for spin-polarized systems. The usage is the same as above.

Related publication
---
Original method of PIO: doi.org/10.1002/chem.201801220

Extension to spin-polarized systems: doi.org/10.1039/D0CP00127A

A recent review on PIO: doi.org/10.1002/wcms.1469

Disclaimer
---
Copyright (c) 2018 jxzhangcc and fksheong

All rights reserved.

Redistribution and use in source and binary forms are permitted provided that the above copyright notice and this paragraph are duplicated in all such forms and that any documentation, advertising materials, and other materials related to such distribution and use acknowledge that the software was developed by the author. The name of the authors may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
