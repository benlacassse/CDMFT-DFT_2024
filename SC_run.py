#!/usr/bin/env python

import sys
import os
import argparse
import importlib
import numpy as np
import pyqcm
from pyqcm import qcm
from pyqcm.cdmft import CDMFT, frequency_grid
from KH_import import KH_import, SuperCluster_model

#--------------------------------------------------------------------------------------
def to_supercluster_matrix(T, nC, ns, 
                           is_hyb=False, 
                           diag_only=False, 
                           off_diag_only=False):
    """
    Convert a hopping matrix from the cluster to the supercluster basis

    PARAMETERS:
    -----------
        T : np.array, list
            Hopping matrix in the cluster basis. If not hybridization, 
            should be a list of np.arrays, one for each cluster.
        nC : int
            Number of clusters
        ns : int
            Number of sites in each cluster
        is_hyb : bool
            Flag to indicate if the hopping matrix is a hybridization
            function (default: False). Transformation will differ 
            slightly in this case. 
        diag_only : bool
            Flag to indicate that only the diagonal elements of the 
            hopping matrix should be transformed (default: False).
            Other blocks are put to zero.
        off_diag_only : bool
            Flag to indicate that only the off-diagonal elements of the 
            hopping matrix should be transformed (default: False).
            Diagonal blocks are put to zero.

    RETURNS:
    --------
        T_out : np.array
            Hopping matrix in the supercluster basis

    """
    nb = ns*2 # size of each cluster block
    T_out = np.zeros((nC*nb,nC*nb), dtype=complex)

    if is_hyb:
        for iclus in range(nC):
            for jclus in range(nC):
                # up left block
                T_out[iclus*nb:iclus*nb+ns,
                      jclus*nb:jclus*nb+ns] = T[iclus*ns:(iclus+1)*ns,
                                                jclus*ns:(jclus+1)*ns]
                # bottom right block
                T_out[iclus*nb+ns:iclus*nb+2*ns,
                      jclus*nb+ns:jclus*nb+2*ns] = -np.conjugate(
                                                    T[iclus*ns:(iclus+1)*ns,
                                                    jclus*ns:(jclus+1)*ns]
                                                    ).T     
    else:
        mult_diag = 1
        mult_off_diag = 1
        if diag_only:
            mult_off_diag = 0
        elif off_diag_only:
            mult_diag = 0

        for iclus in range(nC):
            # up left block
            T_out[iclus*nb:iclus*nb+ns,
                  iclus*nb:iclus*nb+ns] = T[iclus][:ns,:ns]*mult_diag
            # bottom left block
            T_out[iclus*nb+ns:iclus*nb+2*ns,
                  iclus*nb:iclus*nb+ns] = T[iclus][ns:,:ns]*mult_off_diag
            # up right block
            T_out[iclus*nb:iclus*nb+ns,
                  iclus*nb+ns:iclus*nb+2*ns] = T[iclus][:ns,ns:]*mult_off_diag
            # bottom right block
            T_out[iclus*nb+ns:iclus*nb+2*ns,
                  iclus*nb+ns:iclus*nb+2*ns] = T[iclus][ns:,ns:]*mult_diag

    return T_out

def set_Hyb(I, ret=False):
    """
    Compute the CDMFT host function from the DFT+CDMFT k-dependent 
    hybridization function

    PARAMETERS:
    -----------
        I : pyqcm.model_instance
            Instance of the cluster model
        ret : bool
            Flag to return the host function instead of setting it
            Mostly for debug purposes (default: False)

    """
    # initialize matrices
    d = model.dimGF
    nk = np.sum(K_hyb.K[:,0])
    host = np.empty((w_grid.shape[0], d, d), dtype=complex)
    self = np.empty((w_grid.shape[0], d, d), dtype=complex)
    T_nambu = np.zeros((d,d), dtype=complex)
    T_norm = np.zeros((d,d), dtype=complex)
    # construct the local hopping matrix, we separate into
    # the "nambu" (anonalous) and "normal" parts
    T = I.cluster_hopping_matrix()
    T_nambu = to_supercluster_matrix([I.cluster_hopping_matrix(clus=iclus) for iclus in range(K_hyb.NC)], 
                                     K_hyb.NC, K_hyb.L, off_diag_only=True)
    T_norm = to_supercluster_matrix([I.cluster_hopping_matrix(clus=iclus) for iclus in range(K_hyb.NC)], 
                                     K_hyb.NC, K_hyb.L, diag_only=True)
    # for each cluster, print some information
    for iclus in range(K_hyb.NC):
        print("\n------ Gf & Cluster averages for cluster nbr "+str(iclus+1)
              +"/"+str(K_hyb.NC)+" ------")
        print(I.Green_function_average(clus=iclus).real)
        print(I.cluster_averages(clus=iclus))
    for i,w in enumerate(w_grid):
        # hybkw already contains normal hopping matrix
        # hence no T_norm here
        self[i,:,:] = to_supercluster_matrix([I.cluster_self_energy(w*1j,clus=iclus) for iclus in range(K_hyb.NC)],
                                             K_hyb.NC, K_hyb.L)
        g1 = np.eye(d)*w*1j - self[i,:,:] - T_nambu
        gproj = np.zeros((d,d), dtype=complex)
        for ik,k in enumerate(K_hyb.K):
            h = K_hyb.H[i,:,:,ik]
            h_nambu = to_supercluster_matrix(h, K_hyb.NC, K_hyb.L, is_hyb=True)
            g = g1 - h_nambu 
            g = np.linalg.inv(g)
            gproj += g*K_hyb.K[ik,0]
        gproj /= nk

        # we put back normal part of hopping matrix here for the right 
        # definition of the impurity hybridization function
        host[i,:,:] = np.linalg.inv(gproj) - g1 + T_norm
    if ret:
        return -host
    
    # save self-energy in the dedicated file
    with open(self_file, 'a') as f:
        np.savetxt(f, self.real.reshape(-1,d*d), fmt='%1.8e')
        f.write("\n")

    # set the host function for each cluster. Ignore inter-cluster
    # blocks
    print("\n------ Setting Host functions & proceed ------")
    for iclus in range(K_hyb.NC):
        qcm.set_CDMFT_host(I.label, 
                           w_grid, 
                           host[:,iclus*K_hyb.L*2:(iclus+1)*K_hyb.L*2,
                                iclus*K_hyb.L*2:(iclus+1)*K_hyb.L*2], 
                           iclus, False)

#--------------------------------------------------------------------------------------
class SymmetryError(Exception):
    """Exception raised for errors in the symmetry conditions."""
    def __init__(self, message="SymmetryError: 'with-ansym' requires 'with-nsym' to be enabled."):
        self.message = message
        super().__init__(self.message)

#--------------------------------------------------------------------------------------
if __name__ == '__main__':
    """
    Manages the calculation of the superconducting order parameter
    as a post-processing step after a converged DFT+CDMFT calculation. 
    The hybridization of the correlated degrees of freedom is fixed
    to the hybkw provided in the input.def file. Hence superconductivity
    only affects the self-energy, but not the uncorrelated charge. 

    OPTIONAL ARGUMENTS:
    -------------------
        -i, --input-file [filename]    input file (default: input.def)
        -b, --bath-model [filename]    bath model file without .py 
                                       extension (default: SC_BathModel)
        -p, --poke [value]             value of the pairing field to 
                                       poke the system in an first cdmft
                                       run (default: 0.1)     
        --no-sc                        flag to disable superconductivity
                                       (default: False)
        --beta [beta]                  target inverse temperature for the 
                                       cdmft frequency grid
                                       (default: 50)
        --wc [wc]                      cutoff frequency for the cdmft
                                       frequency grid (default: 2)
        --accur [accur]                accuracy for the cdmft convergence
                                       criterion (default: 1e-4)
        -m, --maxiter [maxiter]        maximum number of iterations for the
                                        cdmft convergence (default: 32)
        --conv [conv]                  type of iteration method for the 
                                       cdmft convergence (default: 'fixed_point')
        # --with-seb                     Flag to include the on-site anomalous
        #                                bath terms usually called seb (default: False)
        --with-nsym                     Introduce the x/y normal symmetry in the 
                                        hopping matrix  
        --with-ansym                    Introduce the anormal anti-symmetry x2-y2 
                                        on the superconductivity, nsym must be True 
                                        (default:False)
    """

    # check arguments provided
    parser = argparse.ArgumentParser(prog='SC_run',
                                        usage='%(prog)s [options]')
    # optional arguments
    parser.add_argument('-i', '--input-file', type = str,
                        default = 'input.def',
                        help = 'input file (default: input.def)')
    parser.add_argument('-b', '--bath-model', type = str,
                        default = 'SC_BathModel',
                        help = ('bath model file without .py' 
                                + ' extension (default: SC_BathModel)')
                        )
    parser.add_argument('-p', '--poke', type = float, default = 0.1,
                        help = ('value of the pairing field to poke the system'
                                + 'in an frist cdmft run (default: 0.1).'
                                + 'If zero and sc is enabled, the program will'
                                + 'run a sc cdmft run without first poke')
                        ) 
    parser.add_argument('--no-sc', action = 'store_true', default = False,
                        help = ('flag to disable superconductivity' 
                               + ' (default: False)')
                        )
    parser.add_argument('--beta', type = float, default = 50,
                        help = ('target inverse temperature for the cdmft'
                                + ' frequency grid (default: 50)')
                        )
    parser.add_argument('--wc', type = float, default = 2,
                        help = ('cutoff frequency for the cdmft'
                                + ' frequency grid (default: 2)')
                        )
    parser.add_argument('--accur', type = float, default = 1e-4,
                        help = ('accuracy for the cdmft convergence'
                                + ' criterion (default: 1e-4)')
                        )
    parser.add_argument('-m','--maxiter', type = int, default = 32,
                        help = ('maximum number of iterations for the'
                                + ' cdmft convergence (default: 32)')
                        )
    parser.add_argument('-a', '--alpha', type = float, default = 0.0,
                        help = ('Damping factor'
                                + ' (default: 0.0)')
                        )
    parser.add_argument('--conv', type = str, default = 'fixed_point',
                        help = ('type of iteration method for the'
                                + ' cdmft convergence (default: fixed_point)')
                        )
    # parser.add_argument('--with-seb', action = 'store_true', default = False,
    #                     help=('Flag to include the on-site anomalous'
    #                           + ' bath terms usually called seb (default: False)')
    #                     )
    parser.add_argument('--with-nsym', action = 'store_true', default = False,
                        help=('Flag to introduce x/y symmetry on the hopping'
                              + ' matrix in the normal state (default: False)')
                        )
    parser.add_argument('--with-ansym', action = 'store_true', default = False,
                        help=('Flag to introduce x2-y2 symmetry on the anormal'
                              + ' state with superconductivity, nsym must be true'
                              + ' (default: False)')
                        )      
    args = parser.parse_args()

    # assign arguments to variables
    input_file = args.input_file
    SC = not args.no_sc
    beta = args.beta
    wc = args.wc
    bath_model = args.bath_model
    accur = args.accur
    maxiter = args.maxiter
    alpha = args.alpha
    conv  = args.conv
    # wseb = args.with_seb
    nsym = args.with_nsym
    ansym = args.with_ansym

    # Raise error if ansym is True without activating nsym
    if not nsym and ansym:
        raise SymmetryError()
    
    # set up some global parameters
    pyqcm.set_global_parameter('Hamiltonian_format', 'E')
    pyqcm.set_global_parameter('parallel_sectors')
    pyqcm.set_global_parameter('Ground_state_method', 'P')
    pyqcm.set_global_parameter('Ground_state_init_last')

    # create the frequency grid
    w_grid = frequency_grid(beta=beta, wc=wc).wr

    # read input from DFT+DMFT preparation
    K_hyb = KH_import(input_file, w_grid, nsym)

    
    # Construct the lattice model
    # SuperC = SuperCluster_model(K_hyb, wseb)
    SuperC = SuperCluster_model(K_hyb, nsym, ansym)
    model = SuperC.get_model()

    # set target
    model.set_target_sectors(['R0:S0' for i in range(K_hyb.NC)])
    
    # set pairing field to poke value
    model.set_parameter('D',args.poke)

    # declare the variational parameters
    # first the normal ones
    varia = []
    varia += [bp.split("=")[0].strip() 
                  for bp in SuperC.bath_params_N.strip().splitlines()
                  if "*" not in bp.split("=")[1]]
    # if superconductivity is enabled, add the anomalous ones
    if SC:
        # we check in bath_params_SC which parameters are independent
        # the dependent are signaled by * in the bath_params_SC string
        varia += [bp.split("=")[0].strip() 
                  for bp in SuperC.bath_params_SC.strip().splitlines()
                  if "*" not in bp.split("=")[1]]
    # if not SC, put all anomalous parameters to zero
    else:
        for bp in SuperC.bath_params_SC.strip().splitlines():
            if "*" not in bp.split("=")[1]:
                model.set_parameter(bp.split("=")[0].strip(),0.0)
    
    # print the model for the user to check
    model.print_model('SuperCluster_model.tsv')

    if (SC and np.abs(args.poke) > 1e-8):
        # define file for saving self-energies at each iteration
        # empty if it exists
        self_file='self_energies_poke.dat'
        if os.path.exists(self_file):
            os.remove(self_file)
        # run the CDMFT with poke
        print("\n##### Running CDMFT with poke #####")
        sol = CDMFT(model, varia, beta=beta, wc=wc, host_function=set_Hyb, 
              method='Nelder-Mead', convergence='self-energy', accur=accur, 
              max_value=1e6, depth=1, accur_dist = 1e-8, file='cdmft_poke.tsv',
              maxiter=maxiter, alpha=alpha, iteration = conv)
    
    # turn off poke and re-converge
    # define file for saving self-energies at each iteration
    self_file='self_energies.dat'
    if os.path.exists(self_file):
        os.remove(self_file)
    print("\n##### Running CDMFT without poke #####")
    model.set_parameter('D',1e-9)
    sol = CDMFT(model, varia, beta=beta, wc=wc, host_function=set_Hyb, 
        method='Nelder-Mead', convergence='self-energy', accur=accur, 
        max_value=1e6, depth=1, accur_dist = 1e-8, file='cdmft.tsv',
        maxiter=maxiter, alpha=alpha, iteration = conv)

