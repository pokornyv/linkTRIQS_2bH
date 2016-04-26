linkTRIQS_2bH
=============
**Description:**
link between TRIQS CT-HYB code and a DMFT code for 2-band Hubbard model. 
Uses functions library _libqmc_ that is included in different repository

**Files**
* Main scripts
 * _xdmft.py_ - main code, runs in single-core mode, links solver to dmft code
 * _cthyb_U.py_ - ct-hyb solver, solves the impurity problem and writes the outputs
 * _params.py_ - reads parameter file _xdmft.in_ using python _ConfigParser_ module
* Helper scripts
 * _sigmaoffdiag.py_ - puts a constant real off-diagonal part to the self-energy file
 * _changesigma.py_ -  changes inverse temperature or number of Matsubara frequencies in self-energy file
 * _changesgnoff.py_ - switches signs in first or second off-diagonal block in self-energy file, used when searching for broken symmetry phases

**Required input**
* _dmft.in_, _hamilt_ (opt. _sigma_) files for dmft
* _xdmft.in_ file for xdmft
* _Umm.dat_ file with interaction matrix
