import os
import numpy as np
import matplotlib.pyplot as plt
import pyqcm
# import cluster model class
from SC_BathModel import ClusterModel

from scipy.interpolate import splrep, BSpline



class interpolated_array:

    def __init__(self, W, _M):
        assert W.shape[0] == _M.shape[0], \
                ("the lengths of the arrays do not match" 
                 + " in 'interpolated_matrix'")
                
        self.M = np.empty(_M.shape[1], dtype=object)
        for i in range(self.M.shape[0]):
            self.M[i] = splrep(W, _M[:,i], s=0)
        
    def eval(self, w):
        Z = np.empty(self.M.shape[0])
        for i in range(self.M.shape[0]):
            Z[i] = BSpline(*self.M[i])(w)
        return Z
 

#-----------------------------------------------------------------------------------------
class KH_import:

    """
    Manages the import of the data generated for the SC calculation, 
    especially the host hybridization function. 
    """

    def __init__(self, filename, w_grid, nsym):
        """
        :param str filename: name of text file containing the data
        :param ndarray w_grid: array of frequencies along the imaginary axis (real)
        :param nsym: boolean, if True we add the presence of x/y 
                        symmetry in hoping matrix
        """
        
        self.nsym = nsym
        self.w_grid = w_grid
        self.nw = len(self.w_grid)
        # extract lines from the input file. 
        with open(filename, 'r') as F:
            self.lines = F.readlines()
        self.nlines = len(self.lines)
        self.current_line = 0

        # reading the model names:  self.model_name (array)
        # number of names provides the number of clusters: self.NC
        self._read_model_names()

        # reading the model parameter string: self.parameter_string
        self._read_param_str()

        # reading pair information: will set 
        # self.L : number of sites in each cluster
        # self.pair : sigind matrix, sets which component of host hybridization are equivalent
        # self.NP : number of non-equivalent pairs
        # self.nterm : number of terms in the uppper diag part of host hybridization
        self._read_pair()

        # reading the cluster hopping terms (impurity T matrix): we discard 
        # the inter-cluster terms as they only enter via 
        # the host hybridization function
        # sets: self.T (dict containing self.NC arrays)
        self._read_clusterT()

        # reading the interaction terms : one interaction matrix per custer
        # sets: self.V (dict containing self.NC arrays)
        self._read_interaction()

        # reading the wavevectors and the frequencies, 
        # sets: self.K, self.NK, self.W, self.NW
        self._read_kpt_and_freq()
        

        # finally read the host hybridization function
        # and construct self.H
        self._read_host_hyb()

        print("Input file "+filename+" read successfully.")
        #......................................................................
        # reading done.


    def find_label(self, s, lines, start, end):
        """
        Searches the line containg the label s, starting from line start
        and ending at line end.
        """
        label_found = False
        for i in range(start, end):
            if s in lines[i]:
                current_line = i
                label_found = True
                break
        if label_found == False: 
            raise ValueError("label '" + s + "' not found in input file")
        return current_line

    def _read_model_names(self):
        """
        Read the model names from the lines. Number of 
        model names provides the number of clusters.
        """
        self.current_line = self.find_label("model name: ", self.lines, 
                                       self.current_line, self.nlines)
        # loop over following lines, if "\n" is found, break
        self.NC = 0
        self.model_name = []
        for i in range(self.current_line+1, self.nlines):
            if self.lines[i] == '\n':
                self.current_line = i
                break
            else:
                self.NC += 1
                self.model_name.append(self.lines[i].strip().split()[0])
    
    def _read_param_str(self):
        """
        Read the bath parameters for each cluster in lines. We have 
        to change the labels _1 to _i where i is the index of the 
        cluster. 
        """
        p_str = []
        self.current_line = self.find_label("model parameters:", self.lines, 
                                       self.current_line, self.nlines)
        _iclus = 1 if (self.NC == 1) else 0
        for i in range(self.current_line+1, self.nlines):
            if self.lines[i] == '\n':
                self.current_line = i+1
                break
            elif "super-cluster" in self.lines[i]:
                _iclus +=1
            else:
                # the replacement here is necessary to differentiate
                # bath operators of each cluster within PyQCM
                if self.nsym:
                    # introduction of the symmetry between x and y 
                    if "5_1" in self.lines[i]:
                        p_str.append(self.lines[i][:4]+str(_iclus)+'=1*'+self.lines[i][:2]+'3_'+str(_iclus)+'\n')
                    elif "6_1" in self.lines[i]:
                        p_str.append(self.lines[i][:4]+str(_iclus)+'=1*'+self.lines[i][:2]+'4_'+str(_iclus)+'\n')
                    else:
                        p_str.append(self.lines[i].replace("_1", "_"+str(_iclus)))
                else:
                    p_str.append(self.lines[i].replace("_1", "_"+str(_iclus)))
        
        self.parameter_string = ["" for i in range(_iclus)]
        for i in range(_iclus):
            for p in p_str:
                if "_"+str(i+1) in p:
                    self.parameter_string[i] += p

    def _read_pair(self):
        """
        Read the pair information from the lines. 
        """

        # nbr of pairs
        self.current_line = self.find_label("site pairs", self.lines, 
                                       self.current_line, self.nlines)
        self.NP = int(self.lines[self.current_line].strip().split()[2])

        # then read the pair matrix
        self.current_line = self.find_label("pair index", self.lines, 
                                       self.current_line, self.nlines)
        self.pair = []
        for i in range(self.current_line+1, self.nlines):
            if self.lines[i] == '\n':
                self.current_line = i
                break
            else:
                self.pair.append(list(map(int, self.lines[i].strip().split())))
        self.pair = np.array(self.pair)

        # number of sites in each cluster
        self.L = int(len(self.pair)/self.NC)

        # number of terms in the uppper diagonal part of host hybridization
        self.nterm = int(0.5*np.shape(self.pair)[0]*(np.shape(self.pair)[0]+1))
        
    def _read_clusterT(self):
        """
        Read the cluster hopping terms (impurity T matrix) from the lines.
        self.T is a dict containing an array for each cluster. We discard
        the inter-cluster terms as they only enter via the host hybridization. 
        """
        self.current_line = self.find_label("hopping terms", self.lines, 
                                       self.current_line, self.nlines)
        self.T = {}
        # loop over the clusters
        for iclus in range(self.NC):
            # for each cluster, we have to find the associated pair indices
            clus_pindex = np.unique(self.pair[iclus*self.L:(iclus+1)*self.L,
                                              iclus*self.L:(iclus+1)*self.L]
                                    )
            self.T[iclus+1] = np.zeros(len(clus_pindex))
            for ind in clus_pindex:
                self.T[iclus+1][ind-np.min(clus_pindex)] = float(self.lines[self.current_line+ind].strip().split()[0])
        self.current_line += self.NP+1
    
    def _read_interaction(self):
        """
        Read the interaction terms from the lines. self.V is a dict
        containing an array for each cluster. 
        """
        self.current_line = self.find_label("interaction terms", self.lines, 
                                       self.current_line, self.nlines)
        self.V = {}
        # loop over the lines, each "super-cluster" line indicates a new cluster
        _iclus = 1 if (self.NC == 1) else 0
        for i in range(self.current_line+1, self.nlines):
            if self.lines[i] == '\n':
                self.current_line = i+1
                break
            elif "super-cluster" in self.lines[i]:
                _iclus +=1
                # initialize the array for the cluster
                self.V[_iclus] = []
            else:
                self.V[_iclus].append(float(self.lines[i].strip().split()[0]))

        # tranform all list in dict to arrays
        for iclus in range(1, self.NC+1):
            self.V[iclus] = np.array(self.V[iclus]) 

    def _read_kpt_and_freq(self):
        """
        Read the wavevectors and the frequencies from the lines. 
        """
        self.current_line = self.find_label("wavevectors", self.lines, 
                                       self.current_line, self.nlines)
        self.NK = int(self.lines[self.current_line].strip().split()[1])
        self.K = np.loadtxt(self.lines[self.current_line+1 : self.current_line+1+self.NK], 
                            delimiter='\t')
        self.current_line += 1 + self.NK
        self.K[:,1:] /= 2*np.pi # divide by 2pi to conform to the pyqcm norm

        # reading the frequencies
        self.current_line = self.find_label("frequencies", self.lines, 
                                       self.current_line, self.nlines)
        self.NW = int(self.lines[self.current_line].strip().split()[1])
        self.W = np.loadtxt(self.lines[self.current_line+1 : self.current_line+1+self.NW], 
                       delimiter='\t')
        self.current_line += 1 + self.NW

    def _read_host_hyb(self):
        """
        Reads the host hybridization function from the lines.
        """
        # Real part
        self.current_line = self.find_label("Hybridization function (real)", self.lines, 
                                       self.current_line, self.nlines)
        Hr = np.loadtxt(self.lines[self.current_line+1 : self.current_line+1+self.NW], 
                        delimiter='\t', usecols=range(self.nterm*self.NK))
        Hr_inter = interpolated_array(self.W, Hr)
        self.current_line += 1 + self.NW

        # Imaginary part
        self.current_line = self.find_label("Hybridization function (imag)", self.lines, 
                                       self.current_line, self.nlines)
        Hi = np.loadtxt(self.lines[self.current_line+1 : self.current_line+1+self.NW], 
                        delimiter='\t', usecols=range(self.nterm*self.NK))
        Hi_inter = interpolated_array(self.W, Hi)
        self.current_line += 1 + self.NW

        # construction of the Hybridization function matrix
        def foldhyb(g):
            gu = np.empty((self.L*self.NC, self.L*self.NC))
            cnt = 0
            for iorb1 in range(self.L*self.NC):
                for iorb2 in range(iorb1,self.L*self.NC):
                    # reconstruct from upper triangle
                    gu[iorb1,iorb2] = g[cnt]
                    gu[iorb2,iorb1] = gu[iorb1,iorb2]
                    cnt += 1
            return gu

        self.H = np.empty((self.nw, self.L*self.NC, self.L*self.NC, self.NK), dtype=complex)
        for i in range(self.nw):
            gr = Hr_inter.eval(self.w_grid[i])
            gi = Hi_inter.eval(self.w_grid[i])
            gr = gr.reshape((self.nterm, self.NK))
            gi = gi.reshape((self.nterm, self.NK))
            for k in range(self.NK):
                self.H[i,:,:,k] = foldhyb(gr[:,k]) + foldhyb(gi[:,k])*1j

        ## here we enforce symmetry of the hybridization function
        ## to follow the provided self.pair, this is to avoid
        ## instabilities especially when SC is emerging
        M = np.sum(np.sum(self.H,axis=3),axis=0)/(len(self.w_grid)*self.NK)
        M_sym = np.zeros(np.shape(M), dtype=complex)
        # average over the symmetry equivalent sites
        for p in np.unique(self.pair):
            M_sym[self.pair==p] = np.mean(M[self.pair==p])
        
        # Simple symmetrization : we assume that A*M = M_sym, hence 
        # A = M_sym*M^-1
        A = np.dot(M_sym, np.linalg.inv(M))
        for i in range(self.nw):
            for k in range(self.NK):
                self.H[i,:,:,k] = np.dot(A, self.H[i,:,:,k])

#-----------------------------------------------------------------------------------------
class SuperCluster_model:
    """
    Manages the creation of the PyQCM lattice model, using the imput 
    from KH_import object. Is adapted to a single-cluster problem (i.e.
    single layer), or multiple-cluster problem (e.g. tri-layer).
    """

    def __init__(self, KH, nsym, ansym):
        """
        Constructs the PyQCM lattice model. 
        :param KH: KH_import object containing the input information
        # :param wseb: boolean, if True we add the on-site anomalous
        #                 terms to the model
        :param nsym: boolean, if True we add the presence of x/y 
                        symmetry in hoping matrix
        :param ansym: boolean, if True we add the presence of x2-y2 
                        symmetry on the superconductivity
        """
        # Construct the cluster objects
        # sets : self.clusters (cluster objects), 
        #        self.CM (cluster model objects and SC bath parameter strings)
        self._construct_clusters(KH, nsym, ansym)

        # Set up the lattice model
        # sets: self.model (PyQCM lattice model object)
        self._setup_lattice_model(KH)

    
    def get_model(self):
        return self.model

    def _construct_clusters(self, KH, nsym, ansym):
        """
        Constructs the set of clusters. We check if some clusters
        are equivalent (they share the same KH.pair matrix). The equivalent
        ones are constructed such that only one impurity problem will
        be solved in the DMFT loop.
        """

        # start by extracting the pair matrix for each cluster
        self.c_pairs = {}
        for iclus in range(KH.NC):
            self.c_pairs[iclus+1] = KH.pair[iclus*KH.L:(iclus+1)*KH.L, 
                                       iclus*KH.L:(iclus+1)*KH.L]

        self.CM = {}                # contains the cluster model objets
        previous = []               # list of pair matrices of non-equivalent clusters
        previous_idx = []           # Also need to register index corresponding to first occurence
        self.which_cm = []               # index to which CM the cluster belongs
        # loop over clusters to create the cluster models
        # There is one cluster model per non-equivalent cluster
        n_cm = 1
        for iclus in range(KH.NC):
            exists = False
            # check if the pair matrix is already in the list of previous
            for i in range(len(previous)):
                if np.array_equal(self.c_pairs[iclus+1], previous[i]):
                    exists = True
                    self.which_cm.append(self.which_cm[previous_idx[i]])
            # if not, create a new cluster model
            if not exists:
                self.which_cm.append(n_cm)
                previous.append(self.c_pairs[iclus+1])
                previous_idx.append(iclus)
                # self.CM[n_cm]= ClusterModel(str(n_cm),wseb)
                self.CM[n_cm]= ClusterModel(str(n_cm), nsym, ansym)
                n_cm += 1
        
        print("\n-------------- Creation of cluster objects --------------")
        # now we can construct the cluster objects
        self.clusters = np.empty(KH.NC, dtype=object)
        # we need a registery of the created clusters
        created={}
        for iclus in range(KH.NC):
            # set the Cu sites, we increase z by 1 for each
            # We don't care about the ordering along z
            # effective distance between planes is taken care
            # of by the host hybridization function
            Cu_sites = ((0, 0, iclus), (1, 0, iclus), 
                        (0, 1, iclus), (1, 1, iclus))

            # if the cluster is equivalent to another already created
            # then we construct the cluster object by feeding the 
            # equivalent cluster object and not model
            # if tuple alrady exists, then create from cluster object
            if (self.which_cm[iclus]) in created.keys():
                indx = created[self.which_cm[iclus]]
                print("Creating cluster "+str(iclus+1)+"/"+str(KH.NC)
                      +" from equivalent cluster nbr",indx+1)
                self.clusters[iclus] = pyqcm.cluster(self.clusters[indx],
                                                 Cu_sites)
            # if not existing, create from model
            else:
                print("Creating cluster "+str(iclus+1)+"/"+str(KH.NC)
                      +" from cluster model")
                self.clusters[iclus] = pyqcm.cluster(self.CM[self.which_cm[iclus]].CM, 
                                                 Cu_sites)
                created[self.which_cm[iclus]] = iclus
        print("---------------------------------------------------------")

    def _setup_lattice_model(self,KH):
        """
        Creates the PyQCM lattice model object from the cluster
        objects constructed in _construct_clusters.
        """   

        # create object 
        self.model = pyqcm.lattice_model("SuperCluster_lattice", 
                                         self.clusters, 
                                         ((2,0,0), (0,2,0), (0,0,KH.NC)), 
                                         ((1,0,0), (0,1,0), (0,0,1)))
        # set some operators
        self.model.anomalous_operator('D', [1,0,0], 1)
        self.model.anomalous_operator('D', [0,1,0],-1)

        # construction of the hopping matrix and interaction matrix
        def fold(g,L,c_pairs):
            gu = np.empty((L, L))
            for i in range(L):
                for j in range(L):
                    # be careful with indexing with c_pair: 
                    # should always start by zero hence I should
                    # substract the minimum value
                    gu[i,j] = g[c_pairs[i,j]-np.min(c_pairs)]
            return gu
        # Each cluster has its own T and V matrix
        for iclus in range(KH.NC):
            Tmat = fold(KH.T[iclus+1], KH.L, self.c_pairs[iclus+1])
            Vmat = fold(KH.V[iclus+1], KH.L, self.c_pairs[iclus+1])
            T_elem = []
            V_elem = []
            no = (self.model.clus[iclus].cluster_model.n_sites 
                  + self.model.clus[iclus].cluster_model.n_bath)
            for x in range(KH.L):
                 for y in range(x, KH.L):
                    if Tmat[x,y] != 0: 
                        T_elem.append((x+1, y+1, Tmat[x,y]))
                        T_elem.append((x+1+no, y+1+no, Tmat[x,y]))
                    if Vmat[x,y] != 0: 
                        # n+n-: full up-down block
                        V_elem.append((x+1, y+1+no, Vmat[x,y]))
                        if y>x:
                            V_elem.append((y+1, x+1+no, Vmat[y,x]))
            self.model.clus[iclus].cluster_model.new_operator('T', 'one-body', T_elem)
            self.model.clus[iclus].cluster_model.new_operator('V', 'interaction', V_elem)

        # set all starting parameters, they may be re-written in SC_run.py
        # by convention the chemical potential is 
        # absorbed in the self-energy
        TV_params = """
        mu = 1e-9
        D = 1e-9
        """
        # concatenate the bath parameters of the non-equivalent clusters
        self.bath_params_SC = ""
        self.bath_params_N = ""
        last = 0
        # first cluster is always non-equivalent
        TV_params += "T_"+str(1)+" = 1.0\n"
        TV_params += "V_"+str(1)+" = 1.0\n"
        self.bath_params_SC += self.CM[self.which_cm[0]].bath_params_SC
        self.bath_params_N += KH.parameter_string[0]
        for iclus in range(1,KH.NC):
            # check if this cluster is equivalent to another
            # if not, we can set parameters
            if self.which_cm[iclus] != self.which_cm[last]:
                last = iclus
                # set the parameters
                TV_params += "T_"+str(iclus+1)+" = 1.0\n"
                TV_params += "V_"+str(iclus+1)+" = 1.0\n"
                self.bath_params_SC += self.CM[self.which_cm[iclus]].bath_params_SC.replace("_1", "_"+str(iclus+1))
                self.bath_params_N += KH.parameter_string[iclus]
            # equivalent, will set parameters to be equal to the non-equivalent
            else:
                TV_params += "T_"+str(iclus+1)+" = 1*T_"+str(last+1)+"\n"
                TV_params += "V_"+str(iclus+1)+" = 1*V_"+str(last+1)+"\n"
                for p in self.CM[self.which_cm[iclus]].bath_params_SC.split("\n"):
                    if p.strip():
                        self.bath_params_SC += (p.strip().split('=')[0].replace("_1", "_"+str(iclus+1))
                                           +"=1*"
                                           +p.strip().split('=')[0]
                                           +"\n")
                for p in KH.parameter_string[iclus].split("\n"):
                    if p.strip():
                        self.bath_params_N += (p.strip().split('=')[0]
                                          +"=1*"
                                          +p.strip().split('=')[0].replace("_"+str(iclus+1), 
                                                                      "_"+str(last+1))
                                          +"\n")
        # remove empty lines to avoid conflict in SC_run.py
        self.bath_params_N = os.linesep.join([s for s in self.bath_params_N.splitlines() if s]) 
        self.bath_params_SC = os.linesep.join([s for s in self.bath_params_SC.splitlines() if s])  
        # set all model parameters
        self.model.set_parameters(os.linesep.join([TV_params,
                                                   self.bath_params_N,
                                                   self.bath_params_SC]))
        print(os.linesep.join([TV_params,
                                                   self.bath_params_N,
                                                   self.bath_params_SC]))

#
##########################################################################################
