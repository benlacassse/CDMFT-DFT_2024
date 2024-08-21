import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HBCCO = False
CCOC = True
separated = True
unified = False
layer_3 = True
layer_4 = False
layer_5 = False
IP0 = False
IP1 = True
OP = False
# for i in range(1,4):
#     x = pd.read_csv('CCOC_Dop.csv').iloc[1].values[2]


if HBCCO:
    
    if separated:
        if layer_3:
            print()
        if layer_4:
            print()
        if layer_5:
            print()
    if unified:
        print()


if CCOC:
    df = pd.read_csv('CCOC_Dop.csv', index_col = 'Couche')
    dop_dict = {}
    for row in ['3l_C', '4l_C', '5l_C']:
        row_dict = {}
        if row == '5l_C':
            row_dict['IP0'] =  np.array(df.loc[row].dropna().tolist()[:5])
        row_dict['IP1'] = np.array(df.loc[row].dropna().tolist()[-10:-5])
        row_dict['OP'] = np.array(df.loc[row].dropna().tolist()[-5:])
        dop_dict[row] = row_dict

    df = pd.read_csv('CCOC_S_Msc.csv', index_col = 'Couche')
    sm_dict = {}
    for row in ['3l_C', '4l_C', '5l_C']:
        row_dict = {}
        if row == '5l_C':
            row_dict['IP0'] =  abs(np.array(df.loc[row].dropna().tolist()[:5]))
        row_dict['IP1'] = abs(np.array(df.loc[row].dropna().tolist()[-10:-5]))
        row_dict['OP'] = abs(np.array(df.loc[row].dropna().tolist()[-5:]))
        sm_dict[row] = row_dict

    if separated:
        if layer_3:
            if IP1:
                x = dop_dict['3l_C']['IP1']
                y = sm_dict['3l_C']['IP1']
                plt.plot(x,y, marker = '^', linestyle ='dashed', label = 'IP1 (seperated)', color = 'goldenrod')

            if OP:
                x = dop_dict['3l_C']['OP']
        if layer_4:
            print()
        if layer_5:
            print()
    if unified:
        print()



plt.show()










# eff_dop_i = np.array([2.8, 3.4677168, 5.65847639, 8.26775846, 8.54520343])
# eff_dop_e = np.array([2.3, 5.428927, 8.1015256, 10.14835008, 11.28090702])

# coef_i = abs(np.array([0.000103306, -0.0009976704265392772, -0.080911916, -0.09090267, -0.092865159]))
# coef_e = abs(np.array([0.000195955, -6.731532374127056e-05, -0.089689408, -0.091366085, -0.090797264]))

# old_coef_i = abs(np.array([0.00098891, 0., 0.07209275, 0.0878354, 0.09039749]))
# old_coef_e = abs(np.array([0.00282208, 0., 0.0853819, 0.09052051, 0.08978711]))

# inside = True
# outsie = True

# if inside:
#     plt.plot(eff_dop_i, coef_i, marker = '^', label = 'IP', color = 'goldenrod')
#     plt.plot(eff_dop_i, old_coef_i, marker = '^', linestyle = 'dashed', label = 'IP (seperated)', color = 'goldenrod')

# if outsie:
#     plt.plot(eff_dop_e, coef_e, marker = 'o', label = 'OP', color = 'midnightblue')
#     plt.plot(eff_dop_e, old_coef_e, marker = 'o', linestyle = 'dashed', label = 'OP (seperated)', color = 'midnightblue')

# plt.xlabel('%'+ ' Eff. Doping')
# plt.ylabel('$|m_{SC}|$')
# plt.legend()

# x = np.array([0, 0.001, 0.01, 0.1])
# externe = abs(np.array([-0.000183774, -0.000346736,-0.030791829, -0.032079848]))
# interne = abs(np.array([-1.08976E-05, -0.000127063,-0.077712362, -0.07823838]))
# plt.plot(x, interne, marker = '^', label = 'IP', color = 'goldenrod')
# plt.plot(x, externe, marker = 'o', label = 'OP', color = 'midnightblue')
# plt.xscale('log')
# plt.xlabel('Poke')
# plt.ylabel('$|m_{SC}|$')
# plt.legend()
# plt.show()