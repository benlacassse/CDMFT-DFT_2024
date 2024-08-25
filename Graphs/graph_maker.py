import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def specific_n(n_couche, material, data = ['w_prox'], poke = 0.001, save = False):
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
    plt.legend()
    if save:
        # Create the folder if it doesn't exist
        folder_name = 'saved graphs'
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        
        # Construct the file name
        file_name = f"{n_couche}_{material}_{poke}_{'_'.join(data)}.png"
        file_path = os.path.join(folder_name, file_name)
        
        # Save the plot
        plt.savefig(file_path)
        print(f"Graph saved as {file_path}")
        plt.close()
    else:
        plt.show()

    return 0


specific_n(4, 'HBCCO', data = ['wo_prox', 'w_prox'], poke = 0.1, save=True)
specific_n(4, 'HBCCO', data = ['wo_prox'], poke = 0.1, save=True)
specific_n(4, 'HBCCO', data = ['w_prox'], poke = 0.1, save=True)