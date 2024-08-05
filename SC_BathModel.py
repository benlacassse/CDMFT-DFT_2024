import numpy as np
import pyqcm

################### ROUTINES ######################

# Class for representation
class rep:
    def __init__(self, _C, _signs, _name):
        self.C = _C
        self.signs = _signs
        self.name = _name


class ClusterModel:
    """
    Regroup cluster model generating functions
    """
    def __init__(self, suffix, wseb):
        # construct cluster operators and model
        # sets: self.CM
        self.cluster_2x2('A1A1B1B1B2B2A2A2',suffix,wseb)

        # instructions for bath parameters
        # non-independent parameters are signaled by * in the string
        # if wseb we add the seb parameters
        if wseb:
            self.bath_params_SC="""
db1_1 =  0.0001
db2_1 =  0.0001
db3_1 =  0.0001
db4_1 =  0.0001
db5_1 =  1*db3_1
db6_1 =  1*db4_1
db7_1 = 0.0001
db8_1 = 0.0001
seb1_1 =  0.0001
seb2_1 =  0.0001
seb3_1 =  0.0001
seb4_1 =  0.0001
seb5_1 =  -1*seb3_1
seb6_1 =  -1*seb4_1
seb7_1 = 0.0001
seb8_1 = 0.0001
        """
        else:
            # if not, only the db
            self.bath_params_SC="""
db1_1 =  0.0001
db2_1 =  0.0001
db3_1 =  0.0001
db4_1 =  0.0001
db5_1 =  1*db3_1
db6_1 =  1*db4_1
db7_1 = 0.0001
db8_1 = 0.0001
        """

    # function to create the cluster model
    def cluster_2x2(self, reps_str, suffix, wseb):
        """
        Contruct the 2x2 cluster bath model
        with the given irreps
        """

        # The 4 authorized irreps
        A1 = rep((0,0), (1, 1, 1, 1), 'A1')
        B1 = rep((2,0), (1, 1,-1,-1), 'B1')
        B2 = rep((0,2), (1,-1, 1,-1), 'B2')
        A2 = rep((2,2), (1,-1,-1, 1), 'A2')
        REPS = {'A1': A1, 'A2': A2, 'B1': B1, 'B2': B2}

        L = len(reps_str)//2
        reps = []
        for i in range(L):
            try:
                reps += [REPS[reps_str[2*i:2*i+2]]]
            except:
                raise('Error in reading the irreps in cluster_2x2')

        # number of sites, bath orbitals
        ns = 4
        nb = len(reps)
        no = ns + nb
        g1 = [3, 4, 1, 2]
        g2 = [2, 1, 4, 3]
        for i in range(nb):
            g1 += [reps[i].C[0]]
            g2 += [reps[i].C[1]]
        self.CM = pyqcm.cluster_model(ns, nb, generators=(g1, g2), 
                                      bath_irrep=True, 
                                      name='2x2_C2v_'+reps_str+'_'+suffix)

        varia_N = []
        varia_SC = []
        # bath energies
        for i in range(nb):
            name = 'eb'+str(i+1)
            lab = i+4+1
            self.CM.new_operator(name, 'one-body', (
                (lab, lab, 1.0),
                (lab + no, lab + no, 1.0)
            ))
            varia_N.append(name)

            if wseb:
                name = 'seb'+str(i+1)
                lab = i+4+1
                self.CM.new_operator(name, 'anomalous', [(lab, lab+no, 1.0),])
                varia_SC.append(name)

            elem = []
            for j in range(ns):
                elem.append((j+1, i+ns+1, reps[i].signs[j]))
                elem.append((j+1+no, i+ns+no+1, reps[i].signs[j]))
            name = 'tb'+str(i+1)
            self.CM.new_operator(name, 'one-body', elem)

            elem = []
            for j in range(ns):
                elem.append((j+1, i+ns+no+1, reps[i].signs[j]))
                elem.append((i+ns+1, j+1+no, reps[i].signs[j]))
            name = 'db'+str(i+1)
            self.CM.new_operator(name, 'anomalous', elem)








