"""
Output related functions
"""
import os
import pickle
import numpy as np

def write_data_to_txt_file(savedir, solu_fL2, solu_fL2_inner, solu_T, Lacc_over_LEdd, QL2Qadv_over_Qrad, logM1dotmin, logM1dotmax, NM1dot, logamin, logamax, M2_in_Msun, q, Na):
    """
    Function to write the data to txt files. Old format
    """

    # write the data into txt files [old formats]
    name_list = ['fL2out', 'fL2in', 'T', 'LaccLEdd', 'QQoQ']
    data_list = [np.log10(solu_fL2), np.log10(solu_fL2_inner), np.log10(solu_T),
                 np.log10(Lacc_over_LEdd), np.log10(QL2Qadv_over_Qrad)]
    for n in range(len(name_list)):
        name = name_list[n] + 'M%.2f_q%.3f' % (M2_in_Msun, q)
        savedata = data_list[n]
        with open(savedir+name+'.txt', 'w') as f:
            f.write('logM1dotmin\tlogM1dotmax\tNlogM1dot\t%.5f\t%.5f\t%.5f\n'
                    % (logM1dotmin, logM1dotmax, NM1dot))
            f.write('logamin\tlogamax\tNa\t%.5f\t%.5f\t%.5f\n'
                    % (logamin, logamax, Na))
            for i in range(NM1dot):
                f.write('\n')
                for j in range(Na):
                    if j == 0:
                        f.write('%.8f' % savedata[i, j])
                    else:
                        f.write('\t%.8f' % savedata[i, j])

def write_data_to_pickle_file(savedir, savename, logM1dotgrid, logagrid, solu_fL2, solu_fL2_inner, solu_the, solu_T, Lacc_over_LEdd, QL2Qadv_over_Qrad, solu_tau, solu_rho):
    """
    Function to write the data to a pickle file
    """

    # write the data into a pickle file
    info = ['info', 'logM1dotgrid', 'logagrid', 'logfL2out', 'logfL2in', 'logtheta(H/R)',
            'logT(outer_disk_temperature)', 'logLaccLEdd(luminosity)', 'logQQoQ(radiative_efficiency)',
            'logtau(vertical_optical_depth)', 'logrho(outer_disk_density)']
    data_all = [info, logM1dotgrid, logagrid, np.log10(solu_fL2), np.log10(solu_fL2_inner),
                np.log10(solu_the), np.log10(solu_T), np.log10(Lacc_over_LEdd),
                np.log10(QL2Qadv_over_Qrad), np.log10(solu_tau), np.log10(solu_rho)]

    #
    full_path_output_file = os.path.join(savedir, savename)
    with open(full_path_output_file, 'wb') as f:
        pickle.dump(data_all, f)
        print('data saved at: {}'.format(full_path_output_file))
