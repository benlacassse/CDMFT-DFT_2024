import matplotlib.pyplot as plt
import numpy as np

eff_dop_i = np.array([2.8, 3.4677168, 5.65847639, 8.26775846, 8.54520343])
eff_dop_e = np.array([2.3, 5.428927, 8.1015256, 10.14835008, 11.28090702])

coef_i = abs(np.array([0.000103306, -0.0009976704265392772, -0.080911916, -0.09090267, -0.092865159]))
coef_e = abs(np.array([0.000195955, -6.731532374127056e-05, -0.089689408, -0.091366085, -0.090797264]))

old_coef_i = abs(np.array([0.00098891, 0., 0.07209275, 0.0878354, 0.09039749]))
old_coef_e = abs(np.array([0.00282208, 0., 0.0853819, 0.09052051, 0.08978711]))

inside = True
outsie = True

if inside:
    plt.plot(eff_dop_i, coef_i, marker = '^', label = 'IP', color = 'goldenrod')
    plt.plot(eff_dop_i, old_coef_i, marker = '^', linestyle = 'dashed', label = 'IP (seperated)', color = 'goldenrod')

if outsie:
    plt.plot(eff_dop_e, coef_e, marker = 'o', label = 'OP', color = 'midnightblue')
    plt.plot(eff_dop_e, old_coef_e, marker = 'o', linestyle = 'dashed', label = 'OP (seperated)', color = 'midnightblue')

plt.xlabel('%'+ ' Eff. Doping')
plt.ylabel('$|m_{SC}|$')
plt.legend()
plt.show()