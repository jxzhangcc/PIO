# PIO
Principal Interacting Orbital

Requirement
---
- Gaussian 09
- NBO 6.0
- Python 2.7
- Numpy 1.13.1

The combination here has been well-tested. Other versions are not guaranteed to work but are welcomed to test.

Tutorial
---
1. Run quantum chemistry calculation
    
    Optmize the structure and calculate the electronic structure using standard quantum chemistry softwares. Gaussian09 is used for following demonstration.
    
    **Input**
    - *FILENAME.gjf*
    
    **Output**
    - *FILENAME.log*
        ordinary Gaussian output file
    - *FILENAME.FChk*
        FormCheck file is required for the current version of PIO analysis, for the purpose of visualizing orbital using GaussView. It can be generated via keyword *formcheck* specified in gjf file, or via running external formchk utility
    - *FILENAME.47*
        A FILENAME.47 file can be generated in Gaussian calculations via keyword *pop=nboread* and nbo keyword *$NBO archive $end* specified in gjf file. As the built-in NBO package in Gaussian 09 is NBO3.0, which is a quite old version of NBO package, an external NBO analysis with more recent version is suggested. 

2. Run an NBO analysis

    Run a NBO analysis on the system of interest. Only analysis done by an external NBO package has been tested and is introduced here.
    
    **Input**
    - *FILENAME.47*
        NBO input file, can be generated via running an internal NBO analysis in Gaussian
    
    The following NBO keywords should be specified in the NBO analysis:
    
        $NBO AONAO=W49 FNAO=W49 DMNAO=W49 SKIPBO $END
    
    **Output**
    - *FILENAME.nbo*
        ordinary NBO output file (if run by an external NBO package)
    - *FILENAME.49*
        extra output required for PIO analysis, containing NAO coefficients, density matrix in NAO basis, and Fock matrix in NAO basis if available

3. Run PIO analysis

    We now run PIO analysis with the following command. Note that FILENAME.49 must exist with the same name as FILENAME.FChk.
    
    **Command:**
    
        python PIO.py FILENAME.FChk
    
    The program will then request for a fragmentation as input to specify two groups of atoms
    
    **Input:**
    
        1,2,3,4 5-8
    
    Two groups of atom IDs should be input here separated by a space. Numbers in each group are separated by a comma. Hyphen is supported for sequential numbers. Atom numbering starts from 1.

    **Output**
    - *FILENAME_pio.txt*
        PIO log file, containing basic information of the PIOs of the system subject to the input fragmentation
    - *FILENAME_pio.FChk*
        GaussView FormCheck file, for orbital visualization using GaussView
    - *FILENAME_pio.raw*
        numpy-formatted raw data file, for debugging purpose

Related publication
---
*Coming soon*

TODO
---
- [ ] Upload the code

Contact
---
*jzhangbm@connect.ust.hk*

*fksheong@connect.ust.hk*

*chzlin@ust.hk*
