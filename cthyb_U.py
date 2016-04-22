# dmft calculation using triqs ct-hyb solver for a 2-band Hubbard model 
# with off-diagonal hybridizations
# ct-hyb solver for matrix form of Coulomb interaction
# Vladislav Pokorny; 2014-2015; pokornyv@fzu.cz

import scipy as sp
from numpy.random import randint
from time import time,ctime
from sys import argv,exit

import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import GfImFreq,BlockGf
from pytriqs.operators import *
from pytriqs.applications.impurity_solvers.cthyb import Solver

from libqmc import *
import params as p

NLoop   = int(argv[1])

J2   = p.J

if p.NBins > 1: p.measure_Gtau = True	# just in case so we don't lose statistics on G(tau)

# log file
if p.NBins == 1:
	logfname = 'xdmft.log'
	scflog   = 'scf.log'
else:
	logfname = 'xdmft_bin.log'
	scflog   = 'scf_bin.log'

MaxMats = OmegaN(p.NMats-1,p.beta) # Matsubaras begin with n = 0
FitMin  = int(0.8*p.NMats) # index of frequency, not frequency itself !
FitMax  = p.NMats

#GF = ReadGW(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,FitMin,FitMax,mpi.rank,p.offdiag,logfname)
#GF = FitTailGreen(GF,p.NBand,FitMin,FitMax,6)
GF = ReadGW2(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,FitMin,FitMax,p.NFit,mpi.rank,p.offdiag,1,logfname)

if p.offdiag: # off-diagonal solver, 1 block NBand x NBand
	gf_struct = { '0': p.bands_T }
else:       # diagonal solver, NBand blocks 1 x 1
	gf_struct = {}
	for band in p.bands_T: gf_struct[band] = [0]

if mpi.is_master_node() and p.NBins == 1:
	stars = '**********************************************************'
	PrintAndWrite('\n'+stars+'\n* CT_HYB solver START, iteration '+str(NLoop)+'.\n'+stars+'\n',logfname)

S = Solver(beta = p.beta, gf_struct = gf_struct, n_iw = p.NMats, n_l = p.NLeg, n_tau = p.NTau)

if p.offdiag: # off-diagonal solver, 1 block NBand x NBand
	for block, gw0 in S.G0_iw: gw0 << GF
else:       # diagonal solver, NBand blocks 1 x 1
	for block, gw0 in S.G0_iw: gw0 << GF[block,block]

# reading the Coulomb matrix from file ####################
U_M = sp.loadtxt('Umm.dat')

# model Hamiltonian and other operators ###################
h_int = Operator()
N_tot = Operator()

for b1 in range(p.NBand):
	N_tot = N_tot + n('0',p.bands_T[b1])
	for b2 in range(p.NBand):
		h_int += U_M[b1][b2]*n('0',p.bands_T[b1])*n('0',p.bands_T[b2])
h_int = h_int / 2.0

if p.SpinFlip: # add spin-flip part
	h_int += p.J*(c_dag('0','Au')*c_dag('0','Bd')*c('0','Ad')*c('0','Bu')+\
	            c_dag('0','Ad')*c_dag('0','Bu')*c('0','Au')*c('0','Bd'))

if p.PairHopping: # add pair hopping part
	h_int += J2*(c_dag('0','Au')*c_dag('0','Ad')*c('0','Bd')*c('0','Bu')+\
	             c_dag('0','Bu')*c_dag('0','Bd')*c('0','Ad')*c('0','Au'))

# solver parameters #######################################
p_D = {}
p_D['h_int']                       = h_int
p_D['partition_method']            = p.partition
p_D['quantum_numbers']             = [N_tot]
p_D['n_cycles']                    = p.NCycles
p_D['length_cycle']                = p.LengthCycle
p_D['n_warmup_cycles']             = p.NWarmup
p_D['random_name']                 = ''
p_D['random_seed']                 = 123*mpi.rank + 567*int(time())
p_D['max_time']                    = -1
p_D['measure_g_l']                 = p.measure_Gleg
p_D['measure_g_tau']               = p.measure_Gtau
p_D['measure_pert_order']          = p.measure_pert
p_D['measure_state_trace_contrib'] = p.measure_tc
p_D['move_shift']                  = True
p_D['move_double']                 = True
p_D['use_trace_estimator']         = False

mpi.barrier()

if mpi.is_master_node():
	if p.NBins == 1:
		#WriteGw(NLoop,NBand,NMats,p.beta,bands_T,S.G0_iw,offdiag,'triqs_files/g0_iw',True,logfname)
		[Gzero_tail1_D,Gzero_tail2_D,Gzero_tail3_D,Gzero_tail4_D] = TailCoeffs(S.G0_iw,p.bands_T,p.offdiag)
		PrintAndWrite('\n',logfname)	
		PrintAndWrite('G0(iw) tail fit: 1 / iw (-)',logfname)
		WriteMatrix(Gzero_tail1_D,p.bands_T,'D',p.offdiag,logfname)
		PrintAndWrite('G0(iw) tail fit: 1 / iw^2 (-) (local impurity levels)',logfname)
		WriteMatrix(Gzero_tail2_D,p.bands_T,'D',p.offdiag,logfname)
		PrintAndWrite('G0(iw) tail fit: 1 / iw^3 (+)',logfname)
		WriteMatrix(Gzero_tail3_D,p.bands_T,'D',p.offdiag,logfname)
		PrintAndWrite('G0(iw) tail fit: 1 / iw^4 (+)',logfname)
		WriteMatrix(Gzero_tail4_D,p.bands_T,'D',p.offdiag,logfname)

	if p.measure_Gleg: PrintAndWrite('Measuring G(L).',logfname)
	if p.measure_Gtau: PrintAndWrite('Measuring G(tau).',logfname)

# run the solver ####################################################
t = time()
for num_bin in range(p.PrevBins,p.PrevBins+p.NBins):
	p_D['random_name'] = ''
	p_D['random_seed'] = randint(0,128)*mpi.rank + 567*int(time())
	S.solve(**p_D)
	mpi.barrier()
	if mpi.is_master_node():
		# processing output from solver
		# processing the Legendre GF ####################################
		if p.measure_Gleg:
			[Gl_iw,Sl_iw,Sl_w0_D,Sl_hf_D,nl_D,ntot_leg] = \
			ProcessTriqsOut(S.G0_iw,S.G_l,p.offdiag,p.NBand,FitMin,FitMax,6,p.bands_T,p.equiv_F,'leg',p.SymmS,p.SymmG,logfname)
		# processing the imaginary time GF ##############################
		if p.measure_Gtau:
			[Gtau_iw,Stau_iw,Stau_w0_D,Stau_hf_D,ntau_D,ntot_tau] = \
			ProcessTriqsOut(S.G0_iw,S.G_tau,p.offdiag,p.NBand,FitMin,FitMax,6,p.bands_T,p.equiv_F,'tau',p.SymmS,p.SymmG,logfname)
		if p.NBins > 1:
			PrintAndWrite('Binning mode, bin '+str(num_bin+1)+'/'+str(p.PrevBins+p.NBins),logfname)
			PrintAndWrite('{0: 3d}\t{1: .3f}'.format(num_bin+1,float(S.average_sign)),'bin_sgn.dat')
			if p.measure_Gleg:
				# write parameters ##########################################
				AppendGnuplotFile(num_bin+1,p.bands_T,nl_D,p.offdiag,'bin_n_leg.dat')
				AppendGnuplotFile(num_bin+1,p.bands_T,Sl_hf_D,p.offdiag,'bin_SEinf_leg.dat')
				# write data files ##########################################
				WriteGleg(num_bin+1,p.NBand,p.NLeg,p.beta,p.bands_T,S.G_l,p.offdiag,'bin_files/gl_bin',logfname)
				WriteGw(num_bin+1,p.NBand,p.NMats,p.beta,p.bands_T,Sl_iw,p.offdiag,'bin_files/sl_iw_bin',True,logfname)
				#WriteGw(num_bin+1,NBand,NMats,beta,p.bands_T,Gl_iw,p.offdiag,'bin_files/gl_iw_bin',True,logfname)
			if p.measure_Gtau:
				# write parameters ##########################################
				AppendGnuplotFile(num_bin+1,p.bands_T,ntau_D,p.offdiag,'bin_n_tau.dat')
				AppendGnuplotFile(num_bin+1,p.bands_T,Stau_hf_D,p.offdiag,'bin_SEinf_tau.dat')
				# write data files ##########################################
				WriteGtau(num_bin+1,p.NBand,p.NTau,p.beta,p.bands_T,S.G_tau,p.offdiag,'bin_files/gtau_bin',logfname)
				WriteGw(num_bin+1,p.NBand,p.NMats,p.beta,p.bands_T,Stau_iw,p.offdiag,'bin_files/stau_iw_bin',True,logfname)
run_time = sp.around((time()-t)/60.0,2)

if mpi.is_master_node():
	PrintAndWrite('Stopping ct-hyb solver after {0: .2f} minutes = {1: .2f} hours.'.format(run_time,run_time/60.0),logfname)
	PrintAndWrite('Number of steps per core @ cycle length: {0: .2e} @{1: 3d}\n'.format(p.NCycles,p.LengthCycle),logfname)
	if p.NBins==1:
		PrintAndWrite('Average sign: {0: .3f}\n'.format(float(S.average_sign)),logfname)
		PrintAndWrite(':TIME-CTHYB:\t{0: .3f}'.format(float(run_time)),scflog)	
		PrintAndWrite(':AVGS:\t{0: .3f}'.format(float(S.average_sign)),scflog)
		# processing the atomic GF ######################################
		if p.measure_Gat:
			Gat_iw = S.G0_iw.copy()
			Gat_iw.zero()
			#nat_D = {}
			if p.offdiag:
				for band1 in p.bands_T:
					for band2 in p.bands_T: 
						Gat_iw['0'] = Fourier(S.atomic_gf['0'])
						#nat_D[band1,band2]  = Gat_iw['0'][band1,band2].total_density()
			else:
				for band in p.bands_T:
					Gat_iw[band] = Fourier(S.atomic_gf[band])
			WriteGw(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,Gat_iw,p.offdiag,'triqs_files/gatom',True,logfname)

		# writing the output ############################################
		# write the atomic Hamiltonian eigenvalues to file ##############
		WriteEig(S.eigensystems,p.NBand,p.bands_T,'triqs_files/eig'+str(NLoop)+'.dat')
		#print S.state_trace_contribs
		PrintAndWrite('{0: 3d}\t{1: .3f}'.format(NLoop,float(S.average_sign)),'sgn.dat')
		if p.measure_Gleg: 
			# write parameters ##########################################
			AppendGnuplotFile(NLoop,p.bands_T,nl_D,p.offdiag,'n_leg.dat')
			AppendGnuplotFile(NLoop,p.bands_T,Sl_w0_D,p.offdiag,'SE0_leg.dat')
			AppendGnuplotFile(NLoop,p.bands_T,Sl_hf_D,p.offdiag,'SEinf_leg.dat')
			# write data files ##########################################
			WriteGw(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,Sl_iw,p.offdiag,'sigma',False,logfname)
			WriteGw(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,Gl_iw,p.offdiag,'triqs_files/gl_iw',True,logfname)
			WriteGleg(NLoop,p.NBand,p.NLeg,p.beta,p.bands_T,S.G_l,p.offdiag,'triqs_files/gl',logfname)
			#WriteGtau(NLoop,NBand,NTau,beta,bands_T,S.Delta_tau,offdiag,'triqs_files/deltatau',logfname)
			PrintAndWrite('Results of measuring G(L), loop: {0: 3d}:'.format(NLoop),logfname)
			PrintAndWrite('Occupation numbers from G(L):',logfname)
			WriteMatrix(nl_D,p.bands_T,'D',p.offdiag,logfname)
			PrintAndWrite('Total occupation number: {0: .5f}\n'.format(float(ntot_leg)),logfname)
			PrintAndWrite('{0: 3d}\t{1: .6f}'.format(NLoop,ntot_leg),'ntot_leg.dat')
			PrintAndWrite('Self - energy at iw0 from G(L):',logfname)
			WriteMatrix(Sl_w0_D,p.bands_T,'D',p.offdiag,logfname)
			PrintAndWrite('Self - energy asymptotics from G(L):',logfname)
			WriteMatrix(Sl_hf_D,p.bands_T,'D',p.offdiag,logfname)

			# calculate order parameters
			[mu,ndmft] = WriteMuAndN(logfname)
			mA = nl_D['Au','Au'] - nl_D['Ad','Ad']
			mB = nl_D['Bu','Bu'] - nl_D['Bd','Bd']
			print 'NLoop NTot,  Ndmft,  mu      <AuBd>, <AdBu>, mag   Ndiff'
			PrintAndWrite('{0: 3d}\t{1: .4f}\t{2: .4f}\t{3: .4f}\t{4: .4f}\t{5: .4f}\t{6: .4f}\t{7: .4f}'\
			.format(NLoop,ntot_leg,ndmft,mu,float(nl_D['Au','Bd']),float(nl_D['Ad','Bu'])\
			,mA+mB,float(sp.fabs(ntot_leg-ndmft))),'order.dat')

		if p.measure_Gtau:
			# write parameters ##########################################
			AppendGnuplotFile(NLoop,p.bands_T,ntau_D,p.offdiag,'n_tau.dat')
			AppendGnuplotFile(NLoop,p.bands_T,Stau_w0_D,p.offdiag,'SE0_tau.dat')
			AppendGnuplotFile(NLoop,p.bands_T,Stau_hf_D,p.offdiag,'SEinf_tau.dat')
			# write data files ##########################################
			#WriteGw(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,Stau_iw,p.offdiag,'sigma',False,logfname)
			WriteGw(NLoop,p.NBand,p.NMats,p.beta,p.bands_T,Gtau_iw,p.offdiag,'triqs_files/gtau_iw',True,logfname)
			WriteGtau(NLoop,p.NBand,p.NTau,p.beta,p.bands_T,S.G_tau,p.offdiag,'triqs_files/gtau',logfname)
			PrintAndWrite('Results of measuring G(tau), loop: {0: 3d}:'.format(NLoop),logfname)
			PrintAndWrite('Occupation numbers from G(tau):',logfname)
			WriteMatrix(ntau_D,p.bands_T,'D',p.offdiag,logfname)
			PrintAndWrite('Total occupation number: {0: .5f}\n'.format(float(ntot_tau)),logfname)
			PrintAndWrite('{0: 3d}\t{1: .6f}'.format(NLoop,ntot_tau),'ntot_tau.dat')
			PrintAndWrite('Self - energy at iw0 from G(tau):',logfname)
			WriteMatrix(Stau_w0_D,p.bands_T,'D',p.offdiag,logfname)
			PrintAndWrite('Self - energy asymptotics from G(tau):',logfname)
			WriteMatrix(Stau_hf_D,p.bands_T,'D',p.offdiag,logfname)
	PrintAndWrite(argv[0]+' done, '+ctime(),logfname)

