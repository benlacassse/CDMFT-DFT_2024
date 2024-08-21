import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def specific_n(n_couche, material, data = ['w_prox'], poke = 0.001):
    marker_dict = {'IP0':'s', 'IP1':'^', 'OP':'o'}
    line_dict = {'w_prox':'solid', 'wo_prox':'dashed'}
    color_dict = {'IP0':'lime', 'IP1':'goldenrod', 'OP':'navy'}
    label_dict = {'w_prox':'', 'wo_prox': ' (seperated)'}
    dop_csv = str(material)+'_Dop.csv'
    df = pd.read_csv(dop_csv, index_col = 'Couche')
    dop_dict = {}
    if n_couche == 5:
        dop_dict['IP0'] =  np.array(df.loc[n_couche].dropna().tolist()[:5])
    dop_dict['IP1'] = np.array(df.loc[n_couche].dropna().tolist()[-10:-5])
    dop_dict['OP'] = np.array(df.loc[n_couche].dropna().tolist()[-5:])

    for presence in data:
        if presence == 'wo_prox':
            df = pd.read_csv(material + '_wo_prox.csv', index_col = 'Couche')
        else:
            df = pd.read_csv(material + '_w_prox_'+str(poke)+'.csv', index_col = 'Couche')
        m_dict = {}
        if n_couche == 5:
            m_dict['IP0'] =  abs(np.array(df.loc[n_couche].dropna().tolist()[:5]))
        m_dict['IP1'] = abs(np.array(df.loc[n_couche].dropna().tolist()[-10:-5]))
        m_dict['OP'] = abs(np.array(df.loc[n_couche].dropna().tolist()[-5:]))

        for key in dop_dict.keys():
            plt.plot(dop_dict[key], m_dict[key], marker = marker_dict[key], label = key + label_dict[presence], linestyle = line_dict[presence], color = color_dict[key])
    plt.title('$n$'+material+', $n=$'+str(n_couche)+', poke='+str(poke))
    plt.xlabel('%'+ ' Eff. Doping')
    plt.ylabel('$|m_{SC}|$')
    return 0


specific_n(3, 'CCOC', data = ['wo_prox', 'w_prox'], poke = 0.001)

plt.legend()
plt.show()














# HBCCO = False
# CCOC = True
# separated = True
# unified = False
# layer_3 = True
# layer_4 = False
# layer_5 = False
# IP0 = False
# IP1 = True
# OP = False
# # for i in range(1,4):
# #     x = pd.read_csv('CCOC_Dop.csv').iloc[1].values[2]


# if HBCCO:
    
#     if separated:
#         if layer_3:
#             print()
#         if layer_4:
#             print()
#         if layer_5:
#             print()
#     if unified:
#         print()


# if CCOC:
#     df = pd.read_csv('CCOC_Dop.csv', index_col = 'Couche')
#     dop_dict = {}
#     for row in ['3l_C', '4l_C', '5l_C']:
#         row_dict = {}
#         if row == '5l_C':
#             row_dict['IP0'] =  np.array(df.loc[row].dropna().tolist()[:5])
#         row_dict['IP1'] = np.array(df.loc[row].dropna().tolist()[-10:-5])
#         row_dict['OP'] = np.array(df.loc[row].dropna().tolist()[-5:])
#         dop_dict[row] = row_dict

#     df = pd.read_csv('CCOC_S_Msc.csv', index_col = 'Couche')
#     sm_dict = {}
#     for row in ['3l_C', '4l_C', '5l_C']:
#         row_dict = {}
#         if row == '5l_C':
#             row_dict['IP0'] =  abs(np.array(df.loc[row].dropna().tolist()[:5]))
#         row_dict['IP1'] = abs(np.array(df.loc[row].dropna().tolist()[-10:-5]))
#         row_dict['OP'] = abs(np.array(df.loc[row].dropna().tolist()[-5:]))
#         sm_dict[row] = row_dict

#     if separated:
#         if layer_3:
#             if IP1:
#                 x = dop_dict['3l_C']['IP1']
#                 y = sm_dict['3l_C']['IP1']
#                 plt.plot(x,y, marker = '^', linestyle ='dashed', label = 'IP1 (seperated)', color = 'goldenrod')

#             if OP:
#                 x = dop_dict['3l_C']['OP']
#         if layer_4:
#             print()
#         if layer_5:
#             print()
#     if unified:
#         print()

# plt.xlabel('%'+ ' Eff. Doping')
# plt.ylabel('$|m_{SC}|$')
# plt.legend
# plt.show()







