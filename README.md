linkTRIQS_2bH
=============
**Description:**
link between TRIQS CT-HYB code and a DMFT code for 2-band Hubbard model. 
Uses functions library *libqmc* that is included in different repository

**Files**
* Main scripts
 * xdmft.py - main code, runs in single-core mode, links solver to dmft code
 * cthyb_U.py - ct-hyb solver, solves the impurity problem and writes the outputs
 * params.py - reads parameter file *xdmft.in* using python *ConfigParser* module
* Helper scripts
 * sigmaoffdiag.py - puts a constant real off-diagonal part to the self-energy file
 * changesigma.py -  changes inverse temperature or number of Matsubara frequencies in self-energy file
 * changesgnoff.py - switches signs in first or second off-diagonal block in self-energy file, used when searching for broken symmetry phases

**Required input**
* *dmft.in*, *hamilt* (opt. *sigma*) files for dmft
* *xdmft.in* file for xdmft
* *Umm* file with interaction matrix
