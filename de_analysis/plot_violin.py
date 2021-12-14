# violin plot
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
from de_analysis import as_numpy
import time
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, \
            VPacker

# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#      '01_DE_clusters/Lung_cl0_cl12/data/Lung_'
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  '02_DE_COPD_vs_Control/Lung_cl012/data/Lung_'
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Blood/02_DE_COPD_vs_control/Blood_0123456/data/Blood_'
#
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  '02_DE_COPD_vs_Control/Lung_cl012_try_effectsize/plots/'
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Blood/02_DE_COPD_vs_control/Blood_0123456/plots/uniqueOur/'
#
# patients_group1 = ['Pat133', 'Pat135', 'Pat139', 'Pat141']#'Pat135','Pat111',
# patients_group2 = ['Pat137', 'Pat155', 'Pat233', 'Pat175']
#
#
# genes_oi=['CORO1A', 'TXNIP', 'FOSB', 'CXCL8', 'ITM2B', 'PYGL', 'LIMD2', 'PLAUR',
#        'BID', 'RPS13', 'SLK', 'RPL36AL', 'OAZ2', 'SF3B6', 'UBE2D1', 'MEGF9',
#        'PTMA', 'RPL35A', 'ISG20', 'RGS18', 'RPLP1', 'CLIC1', 'RNF24', 'ACSL1',
#        'NINJ1', 'TMEM59', 'PAK1', 'LIMK2', 'CALM2', 'CASP8', 'IFIT2', 'NDUFB1',
#        'UBA52', 'NCOA1', 'BTG1', 'CMTM2', 'WASF2']
# genes_oi = ['IRF1','EFHD2','CAP1','FTL','S100A6','PLEK','SRGN','B2M',
#             'IFNAR1','PPP1R15A']
#
#
# ct = '0vs12'
# ct='0123456'
#
# title1 = 'Control'
# title2 = 'COPD'
# # title1 = patient_control[0] + ' Cluster0'
# # title2 = patient_copd[0] + ' Cluster12'
# # title1 = ' Cluster0'
# # title2 = ' Cluster12'
#
# fs = 13  # fontsize
# figsize = (7,4.5)


def group_data_into_list(all_cells_path,
                         ct,
                         patient_control,
                         genes_oi,
                         dissstat,
                         read_pd = True,
                         strcontain='Pt'
                         ):
    '''
    get data to plot
    Parameters
    ----------
    patient_control: list of patient names
    genes_oi: list of genes of interest (e.g. genes which should be plotted)
    dissstat: str, status of disease / Group name

    :return:
    '''


    len_arr = np.zeros(len(patient_control))
    nr_zero = np.zeros(len(patient_control))
    for i_pat in range(0, len(patient_control)):
        directory = all_cells_path + ct + '_'
        if read_pd:
            ending = '.tsv'
        else:
            ending = '.npy'

        if os.path.exists(directory + patient_control[i_pat] + ending):
            # get indices of genes of interest
            index_array = np.zeros(1)
            index_array2 = np.zeros(1)
            if read_pd:
                t1 = time.time()
                genes = pd.read_csv(
                    directory + patient_control[i_pat] +'.tsv', sep="\t", index_col=None, usecols=[0])
                header_col = pd.read_csv(
                    directory + patient_control[i_pat] +'.tsv', sep="\t", index_col=0, nrows=1, header=None)
                t2 = time.time()
                print('read header and index: ' + str(t2-t1))
                # get header column
                # header_col.iloc[0, 0] = 'index'
                header_col_list = header_col.iloc[0]


                # for i in range(0, len(genes_oi)):
                index_array = genes[genes['Unnamed: 0'] == genes_oi].index.values
                # index_array = genes[
                #     genes['Unnamed: 0.1'] == genes_oi].index.values
                # index_array = genes[
                #     genes['index'] == genes_oi].index.values
            else:
                t1 = time.time()
                genes = np.load(directory + patient_control[i_pat] +
                                '_indx.npy', allow_pickle=True)
                header_col_list2= np.load(directory + patient_control[i_pat] + '_col.npy', allow_pickle=True)
                t2 = time.time()
                print('read header and index: ' + str(t2 - t1))

                # for i in range(0, len(genes_oi)):
                index_array = np.where(genes == genes_oi)[0][0]

            # transform float to integer
            index_array = pd.to_numeric(index_array, downcast='signed')
            # index_array.sort()

            # read only certain rows
            rows_to_keep = index_array + 1
            if read_pd:
                all_cells = pd.read_csv(
                    directory + patient_control[i_pat] +'.tsv', sep="\t",
                    index_col=0, header=None, skiprows=lambda x: x not in rows_to_keep)
                all_cells.columns = header_col_list
                #TODO
                # strcontain = 'Pt'
                # strcontain = 'cell'
                d = all_cells.iloc[
                    0, all_cells.columns.str.contains(strcontain)].values
                # title_print = all_cells['Unnamed: 0'].iloc[0]
                title_print = all_cells.index.values[0]

            else:
                all_cells_full = as_numpy.read_numpy_to_df(directory + patient_control[i_pat])
                all_cells = all_cells_full.iloc[index_array,:]
                d = all_cells.values
                # if np.shape(all_cells)[0]>1:    # Dataframe
                title_print = all_cells.name
                # else: # Series
                #     title_print = all_cells.name

            # get patient name
            # pat_name = all_cells.columns[0]
            # pat_name_split = pat_name.split("_")
            # pat_name_legend = pat_name_split[0]

            # mean
            if i_pat == 0:
                mean_pat = np.zeros(len(patient_control))
            mean_pat[i_pat] = np.mean(d)



            d_nonzero = d[np.nonzero(d)]

            len_arr[i_pat] = len(d_nonzero)
            print(i_pat)
            print(len(d_nonzero))
            print('zero counts:')
            print(len(d)-len(d_nonzero))
            nr_zero[i_pat]= (len(d)-len(d_nonzero))*100/len(d)

            if i_pat == 0:
                d_nonzero0_li = [np.array(d[np.nonzero(d)].tolist())]
                d_withzero_li = [np.array(d.tolist())]
            else:
                d_nonzero = np.array(d[np.nonzero(d)].tolist())
                d_withzero = np.array(d.tolist())
                d_nonzero0_li.append(d_nonzero)
                d_withzero_li.append(d_withzero)

            pd_expr_val = pd.DataFrame(np.transpose(d_withzero_li[i_pat]),columns=['cells'])
            pd_patnames = pd.DataFrame(np.empty([np.shape(d_withzero_li[i_pat])[0],2]),columns=['Patient','Group'])
            pd_patnames['Patient'] = patient_control[i_pat]
            pd_patnames['Group'] = dissstat
            if i_pat == 0:
                pd_expr_pat0 = pd.concat([pd_expr_val,pd_patnames],axis=1)
            else:
                pd_expr_pat1 = pd.concat([pd_expr_val, pd_patnames], axis=1)
                pd_expr_pat0 = pd.concat([pd_expr_pat0,pd_expr_pat1],axis=0)
        else:
            print(directory + patient_control[i_pat] + ending + ' not exist')
    max_width = np.max(len_arr)
    percent_width_arr = 1/max_width * len_arr

    return d_nonzero0_li, d_withzero_li, percent_width_arr,title_print, \
           pd_expr_pat0, nr_zero, mean_pat


def plot_violin(all_cells_path,
                ct,
                where_to_save,
                patient_control,
                patient_copd,
                genes_oi,
                title1='Control',
                title2='COPD',
                ylabel='Expression values',
                specify_title = None,
                xticklabels = None,
                legend=True,
                fs=15,
                figsize=(12, 7.5),
                read_pd = True,
                subtitle = '',
                strcontain = '',
                swarmplot = False,
                plot_percentage_zero=True,
                palette=['#CD5C5C','#85ABDD']):
    '''
    Plotting function violin plots (cell expression per gene) for
    multiple patients: 1 group versus another group
    :param all_cells_path: str
        path to data files (including path + fileprename), e.g.
        '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
                 '02_DE_COPD_vs_Control/Lung_cl012/data/Lung_'
    :param where_to_save: str
        path to where the plots should be saved
    :param patient_control: list of strings
        list of patient names of group 1
    :param patient_copd: list of strings
        list of patient names of group 2
    :param genes_oi: list of strings
        names of genes which should be plotted,
        e.g., ['gene1'], or ['gene1, 'gene2',..]
    :param ct: string
        name of celltype (called in filename)
    :param title1: optional: string, default = 'Control'
        name of group 1
    :param title2: optional: string, default = 'COPD'
        name of group 2
    :param specify_title: string
        You can specify your title of interest, if nothing is given, the gene
        name is printed
    :param fs: int or float
        fontsize, default = 13
    :param figsize: (float, float)
        size of the figure, default = (7, 4.5)
    :param plot_percentage_zero: boolean
        if True: Plot the grey bars that show the percentage of zero counts
    :param palette: list of str
        colors of the two groups
    :return:
        saves plots
    '''



    legend_g1 = title1
    legend_g2 = title2
    withsubplot = False
    withzerocounts = True
    #      ####################################################
    # get data
    for i_gene in range(0,len(genes_oi)):
        d_nonzero_list,\
        d_withzero_li, \
        percent_width_arr1,\
        title_print, pd_expr_patcontrol, nr_zero1, mean_pat_g1 = \
            group_data_into_list(all_cells_path, ct, patient_control,
                                 genes_oi[i_gene], dissstat=legend_g1 ,
                                 read_pd=read_pd, strcontain=strcontain)


        d_nonzero2_li,\
        d_withzero2_li, \
        percent_width_arr2,\
        title_print, \
        pd_expr_patcopd, nr_zero2, mean_pat_g2 = \
            group_data_into_list(all_cells_path, ct, patient_copd,
                                 genes_oi[i_gene], dissstat=legend_g2,
                                 read_pd=read_pd, strcontain=strcontain)

        # concat control and copd
        pd_expr = pd.concat([pd_expr_patcontrol, pd_expr_patcopd], axis=0)
        # number of zero counts
        nr_zero_ges = np.append(nr_zero1, nr_zero2)
        #
        # fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(18,8))
        num_subplot = len(genes_oi)
        num_row = int(np.round(np.sqrt(num_subplot)))
        num_col = int(np.ceil(num_subplot / num_row))

        # prepare figures
        if withsubplot:
            if i_gene == 0:
                fig, axes = plt.subplots(num_row, num_col, squeeze=False, figsize=figsize)
                axes = dict(zip(genes_oi, axes.flat))
                if num_row > 1 or num_col > 1:
                    title3 = 'subplots_'+ str(len(genes_oi))
                else:
                    title3 = str(len(genes_oi))
                ax = [genes_oi[i_gene]]
        else:
            fig, axes = plt.subplots(squeeze=False, figsize=figsize)
            title3 = str(len(genes_oi))
            ax = 0

    #   ##################################################
        # plot violinplot
        # ax = plt.figure(figsize=figsize)
        t1 = time.time()
        # ax = sns.violinplot(x="Patient", y="cells", data=pd_expr,
        #                                             palette="Set2", split=False, scale="count",
        #                                             inner="stick", cut=0, hue='Group', dodge=False,
        #                                             bw=0.1)

        #palette = ['#b2df8a', '#F08080']
        ax = sns.violinplot(x="Patient", y="cells", data=pd_expr,
                                                    palette=palette, split=False, scale="count",
                                                    cut=0, hue='Group', dodge=False,
                                                    bw=0.2)
        x_pos = np.arange(len(nr_zero_ges))
        colorbars = 'white'
        t2 = time.time()
        print('time plot violin: ' + str(t2-t1))

        # plot mean
        mean_pat = np.concatenate([mean_pat_g1, mean_pat_g2])
        ax.plot(x_pos, mean_pat, 'k',label='mean')


        # plot bars for percentage of zero counts -------------------------
        if plot_percentage_zero:
            t1 = time.time()
            ax.bar(x_pos, nr_zero_ges / 50, width=0.3, align='edge', color=colorbars,
                   edgecolor='k', label='% zero counts',linewidth=1.4, alpha=0.2)

            # colors1 = np.repeat('mediumseagreen', len(patient_control))
            # colors2 = np.repeat('sandybrown', len(patient_copd))
            colors1 = np.repeat('#6CA62A', len(patient_control))
            colors2 = np.repeat('#D25B18', len(patient_copd))
            colorsss = np.concatenate((colors1, colors2), axis=None)


            lw1 = np.repeat(1, len(patient_control))
            lw2 = np.repeat(0.9, len(patient_copd))
            lw = np.concatenate((lw1, lw2), axis=None)
            linsty1 = np.repeat('--', len(patient_control))
            linsty2 = np.repeat('-.', len(patient_copd))
            linsty = np.concatenate((linsty1, linsty2), axis=None)


            for i_l in range(0,len(nr_zero_ges)):
                ax.axhline(y=nr_zero_ges[i_l]/50, xmin=0, xmax=1, color=colorsss[i_l], linestyle=linsty[i_l],linewidth=lw[i_l])

            # ax2 = plt.bar(range(0,6),nr_zero_ges/50,color=colorsss)
            t2 = time.time()
            print('time plot bars: ' + str(t2 - t1))

        # plot swarmplot on top of violinplot for dots ---------------------
        if swarmplot:
            ax = sns.swarmplot(x="Patient", y="cells", color="k", size=5, data=pd_expr,
                          ax=ax);

        t1 = time.time()
        # legend
        if legend:
            patches, labels = ax.get_legend_handles_labels()
            # ax.legend([patches[0], patches[1], patches[2]],
            #           [title1, title2, '% zero counts'], loc='best')

            ax.legend([patches[0], patches[1],patches[2]],
                  ['mean', title1, title2, '% zero counts'], loc='upper left',
                      fontsize=fs,frameon=False)
        else:
            ax.legend_.remove()


        if xticklabels:
            ax.set_xticklabels(xticklabels)
        else:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

        if specify_title:
            ax.set_title(specify_title, fontsize=fs)
        else:
             ax.set_title(title_print + ' - ' + subtitle, fontsize=fs)
        ax.set_xlabel('Patients', fontsize=fs)
        # ax.set_ylabel('')
        ax.set_ylabel(ylabel, fontsize=fs)
        # axes.set_xticks(range(1, len(patient_control)+1))
        # axes[0].set_xticklabels(np.append('0',patient_control))


        ax.tick_params(labelsize=fs)
        # axes[0].set_yticklabels(fontsize=fs)
        ylim0 = ax.get_ylim()[1]
        ax.set_ylim([0, ylim0])
        xti = ax.get_xticks()
        xlim0 = np.max(xti)
        ax.set_xticks(xti)
        ax.set_xlim([-0.5,xlim0+0.5])


        if plot_percentage_zero:
            # add 2nd y-axis, grey
            xv = ax.get_xlim()[1]
            colorlines = 'darkgrey'
            ax.axvline(x=xv, ymin=0, ymax=1, color=colorlines, zorder=2, linewidth=3)
            ytick = ax.get_yticks()
            # for i_y in range(1, len(ytick) - 1):
            i_y = 2
            ax.hlines(y=i_y - 0.005, xmin=xv-0.1, xmax=xv, color=colorlines,
                      zorder=1)
            ax.text(xv + 0.1, i_y - 0.05, '100', color=colorlines,
                    fontsize=fs)
            ax.axhline(y=i_y, xmin=0, xmax=1, color=colorlines, linestyle=':',linewidth=0.7)
            i_y = 1
            ax.hlines(y=i_y - 0.005, xmin=xv - 0.1, xmax=xv, color=colorlines,
                      zorder=1)
            ax.text(xv + 0.1, i_y - 0.05, '50', color=colorlines,
                    fontsize=fs)
            # ax.axhline(y=i_y, xmin=0, xmax=1, color=colorlines, linestyle=':',
            #            linewidth=0.7)
            #
            # add ylabel

        #
            # ybox1 = TextArea("Expression values / ",
            #                  textprops=dict(color="k", size=fs, rotation=90, ha='left',
            #                                 va='bottom'))

            ybox2 = TextArea(" Percentage zero counts ",
                             textprops=dict(color=colorlines, size=fs+1, rotation=90,
                                            ha='left',
                                            va='bottom'))
            ybox = VPacker(children=[ybox2], align="center", pad=0, sep=5)
            anchored_ybox = AnchoredOffsetbox(loc=10, child=ybox, pad=0.,
                                              frameon=False,
                                              bbox_to_anchor=(1.1, 0.55),
                                              bbox_transform=ax.transAxes,
                                              borderpad=0.)#(1.04, 0.55)
            ax.add_artist(anchored_ybox)

        t2 = time.time()
        print('doig rest plot preparations: ' + str(t2-t1))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        # t1 = time.time()
        if swarmplot:
            nametosave2 = '_catswarmplot_withzero_COPD_Control'
        else:
            nametosave2 = '_catplot_withzero_COPD_Control'

        t1 = time.time()
        fig.tight_layout()
        plt.savefig(where_to_save + title_print +'_'+ title3 +'_cl' + ct + '_'+
                    nametosave2 + '.png',bbox_inches='tight')
        plt.savefig(
            where_to_save + title_print + '_' + title3 + '_cl' + ct + '_' +
            nametosave2 + '.svg', bbox_inches='tight')
        t2 = time.time()
        print('time saving png: ' + str(t2 - t1))

        print('finished')
        print(i_gene)




    return

## example:
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Blood/' \
#      '01_between_clusters/0_vs_rest/all_pat/data/Blood_'
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/Blood/' \
#      '01_between_clusters/0_vs_rest/plots/'
#
# # patients_group1 = ['Pat111','Pat135','Pat172']
# # patients_group2 = ['Pat133','Pat139','Pat141']
# cl1 = '0'
# cl2  = '123456'
# ct = cl1 + 'vs' + cl2
# patients_group1 = ['Pat133_'+cl1,'Pat135_'+cl1,'Pat137_'+cl1,'Pat139_'+cl1,'Pat141_'+cl1,
#         'Pat155_'+cl1,'Pat175_'+cl1,'Pat233_'+cl1]
# patients_group2= ['Pat133_'+cl2,'Pat135_'+cl2,'Pat137_'+cl2,'Pat139_'+cl2,'Pat141_'+cl2,
#        'Pat155_'+cl2,'Pat175_'+cl2,'Pat233_'+cl2]
#
# title1 = 'cluster '+ cl1
# title2 = 'cluster '+ cl2
#
# genes_oi = ['S100A6','FOS','STAT1','NEAT1']
#
# plot_violin(all_cells_path,ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title1, title2=title2)
##
# all_cells_path = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'compare_edgeR_myData_control/data/AllCells_'
# where_to_save = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'compare_edgeR_myData_control/de_results/plots/'
# ct = '0'
# # patients_group1 = ['Pat111','Pat135','Pat172']
# # patients_group2 = ['Pat133','Pat139','Pat141']
# patients_group1 = ['Pat111']
# patients_group2 = ['Pat133','Pat139','Pat141','Pat135','Pat172']
# genes_oi = ['FTL']
#
# plot_violin(all_cells_path,ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title1, title2=title2)
##
# all_cells_path = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'sc_simulation_w_muscat/kang/pandasDF/logcounts/de10_ss4;1/data/de10_ss4;1_'
#
# where_to_save = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'sc_simulation_w_muscat/kang/pandasDF/logcounts/de10_ss4;1/plots/'
# ct = 'cl1'
# # patients_group1 = ['Pat111','Pat135','Pat172']
# # patients_group2 = ['Pat133','Pat139','Pat141']
# patients_group1 = ['Pt1','Pt2','Pt3']
# patients_group2 = ['Pt4','Pt5','Pt6']
# title1 = 'Group1'
# title2 = 'Group2'
# genes_oi = ['gene3825','gene5304','gene3930','gene2899','gene1474','gene4635',
#             'gene2769','gene1277','gene4560']
#
# plot_violin(all_cells_path,ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title1, title2=title2)
##
# example:
# all_cells_path = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'sc_simulation_w_muscat/kang/pandasDF/logcounts/de10;1/data/simdata_cl'
# where_to_save = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#                  'sc_simulation_w_muscat/kang/pandasDF/logcounts/de10;1/de_results/figures/'
#
# patients_group1 = ['Pt1','Pt2','Pt3']
# patients_group2 = ['Pt4','Pt5','Pt6']
# ct = '1'
#
#
# title_group1 = 'Group 1'
# title_group2 = 'Group 2'
#
# genes_oi = ['gene1382']
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2)

##
#
# all_cells_path = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#      'own_simulation_example/different_distributions/cl01/data/simdata_cl'
# # all_cells_path = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
# #      'own_simulation_example/same_distributions/cl0/data/simdata_cl'
# where_to_save = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
#      'own_simulation_example/different_distributions/cl01/de_results/figures/'
# # where_to_save = '/home/erika/Documents/Projects/Evaluation_DE_method/' \
# #      'own_simulation_example/same_distributions/cl0/de_results/figures/'
#
# patients_group1 = ['pat_G1_0','pat_G1_1','pat_G1_2','pat_G1_3','pat_G1_4','pat_G1_5']
# patients_group2 = ['pat_G2_0','pat_G2_1','pat_G2_2','pat_G2_3','pat_G2_4','pat_G2_5',
#                    'pat_G2_6','pat_G2_7','pat_G2_8']
#
# ct = '01'
#
#
# title_group1 = 'Group 1'
# title_group2 = 'Group 2'
# goi = ['65','10','88','1','54']
# for i_g in range(0,5):#9):
#     # genes_oi = ['gene'+ str(i_g)]
#     genes_oi = ['gene' + goi[i_g]]
#
#     strcontain = 'cell' # prename in data columns 'Pt'
#
#     # p_val1 = ['0.003','0.011','0.023','0.0079','0.0039','0.0039','0.0039',
#     #           '0.0039','0.0079','0.0079']#7
#     # p_val1 = ['0.09','0.13','0.11','0.08','0.07','0.11','0.19','0.00',
#     #           '0.02','0.01']#8
#     # p_val1 = ['0.20','0.33','0.03','0.00','0.22','0.04','0.02','0.00','0.00','0.36']#9
#     p_val1 = ['0.0004','0.0026','0.0018','0.0038','0.2452','0.0008','0.0328',
#               '0.0156','0.0010','0.0030']#00
#     p_val1 = ['0.11','0.48','0.25','0.12','0.37','0.29','0.05','0.34','0.25','0.08']#0same
#     p_val1 = ['0.486','0.478','0.348','0.293','0.290']
#     p_val2 = ['0.0099','0.117','0.25','0.099','0.12']
#     # p_val2 = ['0.98','0.86','0.96','0.83','0.56','0.76','0.88','0.99',
#     #           '0.36','0.87']#7
#     # p_val2 = ['0.99','0.86','0.89','0.99','0.98','0.94','0.96','0.86',
#     #           '0.95','0.93']#8
#     # p_val2 = ['0.70','0.06','0.69','0.76','0.41','0.93','0.51','0.16','0.26','0.03']#9
#     # p_val2 = ['0.42','0.93','0.26','0.36','0.08','0.89','0.08','0.05','0.25','0.55']#00
#     # p_val2 = ['0.25','0.35','0.65','0.15','0.94','0.73','0.35','0.78','0.99','0.35']#0same
#
#     subtitle = 'p-val_our = '+ str(p_val1[i_g]) + \
#                ', p-val_edgeR = ' + str(p_val2[i_g])
#     plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#                 genes_oi, title1=title_group1, title2=title_group2, subtitle= subtitle,
#                 strcontain=strcontain,swarmplot = False)
#
##
# all_cells_path = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#            'benchmarking_edgeR_vs_our/data/allCells_Filtered_025genes_'
# where_to_save = '/home/erika/PycharmProjects/Kevin_GeneSetEnrichmentAnalysis/' \
#            'benchmarking_edgeR_vs_our/plots/'
#
#
# patients_group1 = ['Pat135','Pat133', 'Pat139', 'Pat111', 'Pat141', 'Pat172']
# patients_group2 = ['Pat137',  'Pat142', 'Pat115', 'Pat145',  'Pat155', 'Pat175',
#             'Pat149', 'Pat173','Pat162']
#
# ct = '1'
#
# title_group1 = 'control'
# title_group2 = 'COPD'
# strcontain = 'Pt'
#
# genes_oi = ['FTL']
#
# # subtitle = 'p-val our=0.04, p-val edgeR=0.85'
# subtitle = 'p-val our=0.33, p-val edgeR=0.02'
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2, subtitle=subtitle,
#             strcontain=strcontain)
##
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Oktober20/02_DE_COPD_vs_control/ct0/data/normalized201020_'
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Oktober20/02_DE_COPD_vs_control/ct0/'
#
# patients_group1 = ['Pat133', 'Pat135', 'Pat139', 'Pat141']
# patients_group2 = ['Pat137', 'Pat155', 'Pat233', 'Pat175']
# ct = '0'
#
# title_group1 = 'control'
# title_group2 = 'COPD'
# strcontain = 'Pt'
#
# genes_oi = ['RGS2','CFLAR','NAMPT']
#
# # subtitle = 'p-val our=0.04, p-val edgeR=0.85'
# subtitle = ''
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2, subtitle=subtitle,
#             strcontain=strcontain)
##
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Dec20/Lung/02_DE_COPD_vs_control/lung_feb_2021_012/' \
#                  'de_results/allCells_Filtered_0.05genes_'
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Dec20/Lung/02_DE_COPD_vs_control/lung_feb_2021_012/' \
#                  'de_results/plots/'
#
# # patients_group1 = ['Pat133_3', 'Pat137_3', 'Pat139_3', 'Pat141_3',
# #                        'Pat142_3','Pat149_3']
# # patients_group2 = ['Pat133_0124', 'Pat135_0124', 'Pat139_0124',
# #                        'Pat141_0124','Pat142_0124', 'Pat190_0124',
# #                        'Pat192_0124','Pat137_0124', 'Pat155_0124',
# #                        'Pat233_0124', 'Pat175_0124','Pat149_0124','Pat173_0124']
# patients_group1 = ['Pat111','Pat133', 'Pat135', 'Pat139', 'Pat141','Pat172']
# patients_group2 = ['Pat115','Pat137','Pat142', 'Pat145', 'Pat155', 'Pat162', 'Pat175']
# ct = '012'
#
# title_group1 = 'control'
# title_group2 = 'COPD'
# strcontain = 'Pt'
#
# genes_oi = ['FKBP1A','PTPN2','RPL10A']
#
# # subtitle = 'p-val our=0.04, p-val edgeR=0.85'
# subtitle = ''
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2, subtitle=subtitle,
#             strcontain=strcontain)
##
# all_cells_path = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Protein_mar21/01_between_clusters/0vs1234/de_results/' \
#                  'allCells_Filtered_0.05genes_'
# where_to_save = '/home/erika/Documents/Projects/Theo_Maerz_2020/Analysis/' \
#                  'Protein_mar21/01_between_clusters/0vs1234/de_results/plots/'
#
# ct1 = '0'
# ct2 = '1234'
# # patients_group1 = ['Pat133_3', 'Pat137_3', 'Pat139_3', 'Pat141_3',
# #                        'Pat142_3','Pat149_3']
# # patients_group2 = ['Pat133_0124', 'Pat135_0124', 'Pat139_0124',
# #                        'Pat141_0124','Pat142_0124', 'Pat190_0124',
# #                        'Pat192_0124','Pat137_0124', 'Pat155_0124',
# #                        'Pat233_0124', 'Pat175_0124','Pat149_0124','Pat173_0124']
# # patients_group1 = ['Pat111','Pat133', 'Pat135', 'Pat139', 'Pat141','Pat172']
# # patients_group2 = ['Pat115','Pat137','Pat142', 'Pat145', 'Pat155', 'Pat162', 'Pat175']
# # patients_group1 = ['Pat4', 'Pat5', 'Pat6']
# # patients_group2 = ['Pat1', 'Pat2', 'Pat3']
# patients_group1 = ['Pat1'+'_'+ct1, 'Pat2'+'_'+ct1, 'Pat3'+'_'+ct1,
#                    'Pat4'+'_'+ct1, 'Pat5'+'_'+ct1, 'Pat6'+'_'+ct1]
# patients_group2 = ['Pat1'+'_'+ct2, 'Pat2'+'_'+ct2, 'Pat3'+'_'+ct2,
#                    'Pat4'+'_'+ct2, 'Pat5'+'_'+ct2, 'Pat6'+'_'+ct2]
#
# ct = '0vs1234'
#
# title_group1 = 'control'
# title_group2 = 'COPD'
# strcontain = 'Pt'
#
# genes_oi = ["AB-CD44:G44-26",
# "AB-CD15",
# "AB-CD294",
# "AB-CD66",
# "AB-CD62L:DREG-56"
# ]
#
# # subtitle = 'p-val our=0.04, p-val edgeR=0.85'
# subtitle = ''
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2, subtitle=subtitle,
#             strcontain=strcontain)

# all_cells_path = "/home/erika/Documents/Projects/Evaluation_DE_method/" \
#                  "sc_simulation_w_muscat/kang/pandasDF/counts/db10;1/" \
#                  "data/db10;1_"
# where_to_save = "/home/erika/Documents/Projects/Evaluation_DE_method/" \
#                  "sc_simulation_w_muscat/kang/pandasDF/counts/db10;1/"
#
#
# patients_group1 = ['Ptctrl1015.A', 'Ptctrl1016.A', 'Ptctrl1244.A']
# patients_group2 = ['Ptctrl1015.B', 'Ptctrl1016.B', 'Ptctrl1244.B']
#
# ct = 'CD4 T cells'
#
# title_group1 = 'A'
# title_group2 = 'B'
# strcontain = ''
#
# genes_oi = ['gene737']
#
# # subtitle = 'p-val our=0.04, p-val edgeR=0.85'
# # subtitle = 'p-val our=0.33, p-val edgeR=0.02'
#
# plot_violin(all_cells_path, ct, where_to_save, patients_group1, patients_group2,
#             genes_oi, title1=title_group1, title2=title_group2,
#             strcontain=strcontain)#subtitle=subtitle,