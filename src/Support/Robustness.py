import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Support.Additional import clean_names
import gender_guesser.detector as gender
gendet = gender.Detector()


def traits_robustness(EFO_NoCom):
    """A function to save space in the main notebook which plots the same
    plots the same traits figure as the main body but drops rows which
    have multiple traits in
    """
    plt.style.use('seaborn-ticks')
    plt.rcParams["font.family"] = "Helvetica"
    Trait_Count = EFO_NoCom.groupby('EFO term')['EFO term'].count()
    Trait_Count = Trait_Count.to_frame()
    Trait_Count = Trait_Count.rename(columns={'EFO term': '# Studies'})
    Trait_Count.reset_index(inplace=True)
    Trait_AncSum = EFO_NoCom.groupby('EFO term')['NUMBER OF INDIVDUALS'].sum()
    Trait_AncSum = Trait_AncSum.to_frame()
    Trait_AncSum = Trait_AncSum.rename(columns={'NUMBER OF INDIVDUALS':
                                                'Total Sample'}).reset_index()
    Trait_AssSum = EFO_NoCom.groupby('EFO term')['ASSOCIATION COUNT'].sum()
    Trait_AssSum = Trait_AssSum.to_frame()
    Trait_AssSum = Trait_AssSum.rename(columns={'ASSOCIATION COUNT':
                                                'Total Assoc.'}).reset_index()
    Trait_plot = pd.merge(Trait_AssSum,
                          Trait_AncSum, how='left', on='EFO term')
    Trait_plot = pd.merge(Trait_plot,
                          Trait_Count, how='left', on='EFO term')
    Trait_plot['Studies'] = (
        Trait_plot['# Studies'] / Trait_plot['# Studies'].sum())
    Trait_plot['Studies'] = Trait_plot['Studies'] * 100
    Trait_plot['Associations'] = (
        Trait_plot['Total Assoc.'] / Trait_plot['Total Assoc.'].sum())
    Trait_plot['Associations'] = Trait_plot['Associations'] * 100
    Trait_plot['Sample'] = (
        Trait_plot['Total Sample'] / Trait_plot['Total Sample'].sum())
    Trait_plot['Sample'] = Trait_plot['Sample'] * 100
    Trait_plot = Trait_plot.set_index('EFO term')
    Trait_plot = Trait_plot.sort_values(by='Studies', ascending=False)[0:17]
    Trait_plot = Trait_plot[['Sample', 'Studies', 'Associations']]

    Par_Count = EFO_NoCom.groupby('Parent term')['Parent term'].count()
    Par_Count = Par_Count.to_frame()
    Par_Count = Par_Count.rename(columns={'Parent term': '# Studies'})
    Par_Count = Par_Count.reset_index()
    Par_AncSum = EFO_NoCom.groupby('Parent term')['NUMBER OF INDIVDUALS'].sum()
    Par_AncSum = Par_AncSum.to_frame().rename(columns={'NUMBER OF INDIVDUALS':
                                                       'Total Sample'})
    Par_AncSum = Par_AncSum.reset_index()
    Par_AssSum = EFO_NoCom.groupby('Parent term')['ASSOCIATION COUNT'].sum()
    Par_AssSum = Par_AssSum.to_frame().rename(columns={'ASSOCIATION COUNT':
                                                       'Total Assoc.'})
    Par_AssSum = Par_AssSum.reset_index()
    Par_plot = pd.merge(Par_AssSum, Par_AncSum, how='left', on='Parent term')
    Par_plot = pd.merge(Par_plot, Par_Count, how='left', on='Parent term')
    Par_plot['Studies'] = Par_plot['# Studies'] / Par_plot['# Studies'].sum()
    Par_plot['Studies'] = Par_plot['Studies'] * 100
    Par_plot['Associations'] = (
        Par_plot['Total Assoc.'] / Par_plot['Total Assoc.'].sum()) * 100
    Par_plot['Sample'] = (
        Par_plot['Total Sample'] / Par_plot['Total Sample'].sum()) * 100
    Par_plot = Par_plot.set_index('Parent term')
    Par_plot = Par_plot.sort_values(by='Studies', ascending=False)
    Par_plot = Par_plot[['Sample', 'Studies', 'Associations']]

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(18, 7.5))
    fig_a = Trait_plot.sort_values(by='Studies', ascending=True)
    fig_a = fig_a.plot(kind='barh', ax=axes[0], legend=False,
                       colors=['goldenrod', 'steelblue', 'tomato'],
                       edgecolor='k', width=0.66)
    fig_a.set_ylabel("")
    fig_a.set_xlabel("Percent of Mapped Traits (%)", fontsize=12)
    fig_b = Par_plot.sort_values(by='Studies', ascending=True)
    fig_b = fig_b.plot(kind='barh', ax=axes[1], legend=True,
                       colors=['goldenrod', 'steelblue', 'tomato'],
                       edgecolor='k', width=0.66)
    fig_b.set_ylabel("")
    fig_b.set_xlabel("Percent of Parent Traits (%)", fontsize=12)
    fig_b.legend(prop={'size': 12}, frameon=True, loc='bottom right',
                 bbox_to_anchor=(1, 0.2), edgecolor='k')
    sns.despine()


def gender_robustness(AuthorMaster, EFO_NoCom):
    """A function to save space in the main notebook which plots the same
    plots the same gender figure as the main body but drops rows which
    have multiple traits in
    """
    AuthorMaster['CleanForename'] = AuthorMaster['FORENAME'].map(
        lambda x: clean_names(x))
    AuthorMaster['CleanGender'] = AuthorMaster['CleanForename'].map(
        lambda x: gendet.get_gender(x))
    AuthorMaster['MaleFemale'] = AuthorMaster['CleanGender'].str.replace(
        'mostly_', '')
    AuthorMaster['isfemale'] = np.where(
        AuthorMaster['MaleFemale'] == 'female', 1, 0)
    AuthorMenorWomen = AuthorMaster[(AuthorMaster['MaleFemale'] == 'male') |
                                    (AuthorMaster['MaleFemale'] == 'female')]
    AuthorMaster_merged = pd.merge(AuthorMenorWomen,
                                   EFO_NoCom, how='left',
                                   on='PUBMEDID')
    AuthorMaster_EFO = AuthorMaster_merged.groupby(
        ['EFO term'])['isfemale'].mean().to_frame()
    AuthorMaster_Parent = AuthorMaster_merged.groupby(
        ['Parent term'])['isfemale'].mean().to_frame()

    countperEFO = EFO_NoCom.groupby(
        ['EFO term'])['PUBMEDID'].count().sort_values(ascending=False)[0:17]
    countperEFO = countperEFO.index.tolist()
    meanfemacc = AuthorMaster_merged.groupby(
        ['STUDY ACCESSION'])['isfemale'].mean()
    meanfemacc = meanfemacc.to_frame().reset_index()
    meanfemacc_withparent = pd.merge(EFO_NoCom[['STUDY ACCESSION',
                                                'Parent term']],
                                     meanfemacc, on='STUDY ACCESSION',
                                     how='left')
    meanfemacc_withparent = meanfemacc_withparent[
        ~meanfemacc_withparent['isfemale'].isnull()]

    meanfemacc_wEFO = pd.merge(EFO_NoCom[['STUDY ACCESSION', 'EFO term']],
                               meanfemacc, on='STUDY ACCESSION', how='left')
    meanfemacc_wEFO = meanfemacc_wEFO[meanfemacc_wEFO['EFO term'].isin(
        countperEFO)]
    meanfemacc_wEFO = meanfemacc_wEFO[~meanfemacc_wEFO['isfemale'].isnull()]

    AuthorMaster_EFO.reset_index(inplace=True)
    ranks = AuthorMaster_EFO[AuthorMaster_EFO['EFO term'].isin(countperEFO)]
    ranks = ranks.sort_values(by='isfemale', ascending=False)['EFO term']
    ranks = ranks.tolist()
    sns.set(font_scale=1.5, font="Times New Roman", style='ticks')

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=False, figsize=(18, 7.5))

    p = sns.boxplot(data=meanfemacc_wEFO,
                    y='EFO term', ax=ax1, order=ranks,
                    x='isfemale', whis=1.5, palette="vlag", fliersize=0)
    sns.swarmplot(y='EFO term', x='isfemale', data=meanfemacc_wEFO,
                  order=ranks, size=2,  color="white", edgecolor="black",
                  linewidth=0.5, ax=ax1)
    p.set_ylabel('')
    p.set_xlim(0, 1.05)
    p.set(axis_bgcolor='white')
    p.set_xlabel('Female/Male Ratio Per Paper Across Mapped Terms')
    plt.setp(p.spines.values(), color='k')
    plt.setp([p.get_xticklines(), p.get_yticklines()], color='k')
    sns.despine(bottom=False, left=False, offset=10, ax=ax1)

    AuthorMaster_Parent.reset_index(inplace=True)
    ranks = AuthorMaster_Parent.sort_values(by='isfemale', ascending=False)[
        'Parent term'].tolist()
    p = sns.boxplot(data=meanfemacc_withparent,
                    y='Parent term', ax=ax2, order=ranks,
                    x='isfemale', whis=1.5, palette="vlag", fliersize=0)
    sns.swarmplot(y='Parent term', x='isfemale', data=meanfemacc_withparent,
                  ax=ax2, order=ranks, size=1,  color="white",
                  edgecolor="black", linewidth=0.5)
    p.set_ylabel('')
    p.set_xlabel('Female/Male Ratio Per Paper Across Parent Terms')
    p.set_xlim(0, 1.05)
    sns.despine(bottom=False, left=False, offset=10, ax=ax2)
    sns.set(font_scale=1, font="Times New Roman")

    holdstring = 'Parent terms out of dataframe from the figure above: '
    for index, row in AuthorMaster_Parent.iterrows():
        holdstring = holdstring + \
            row['Parent term'].title() + \
            ' (' + str(round(row['isfemale'], 3) * 100) + '%), '
    print('\n' + holdstring[:-2])

    holdstring = 'Trait terms out of dataframe from the figure above: '
    tempdf = AuthorMaster_EFO[AuthorMaster_EFO['EFO term'].isin(countperEFO)]
    tempdf = tempdf.sort_values(by='isfemale', ascending=False)
    for index, row in tempdf.iterrows():
        holdstring = holdstring + \
            row['EFO term'].title() + \
            ' (' + str(round(row['isfemale'], 3) * 100) + '%), '
    print('\n' + holdstring[:-2])


def funder_robustness(EFO_NoCom, FunderInfo, Cat_Ancestry):
    """A function to save space in the main notebook which plots the same
    plots the same funder heatmap as the main body but drops rows which
    have multiple traits in
    """
    EFO_NoCom['Parent term'] = EFO_NoCom['Parent term'].str.replace(
        'Disorder', '')
    FInf_Par = pd.merge(FunderInfo, EFO_NoCom, left_on="PUBMEDID",
                        right_on="PUBMEDID", how="left")
    df_index = FunderInfo.groupby(['Agency'])['Agency'].count()
    df_index = df_index.sort_values(ascending=False)
    cols = FInf_Par.groupby(['Parent term'])['Parent term'].count()
    cols = cols.sort_values(ascending=False)
    fun_par = pd.DataFrame(index=df_index.index.values[0:20].tolist(),
                           columns=cols.index.values.tolist())

    for par in cols.index.values.tolist():
        for fun in df_index.index.values[0:20].tolist():
            fun_par[par][fun] = len(FInf_Par[(FInf_Par['Agency'] == fun) & (
                FInf_Par['Parent term'] == par)])

    FInf_Anc = pd.merge(FunderInfo, Cat_Ancestry, left_on="PUBMEDID",
                        right_on="PUBMEDID", how="left")

    cols = FInf_Anc.groupby(['Broader'])['Broader'].count()
    cols = cols.sort_values(ascending=False)

    fun_anc = pd.DataFrame(index=df_index.index.values[0:20].tolist(),
                           columns=cols.index.values.tolist())

    for anc in cols.index.values.tolist():
        for fun in df_index.index.values[0:20].tolist():
            fun_anc[anc][fun] = len(FInf_Anc[(FInf_Anc['Agency'] == fun) & (
                FInf_Anc['Broader'] == anc)])

    sns.set(font_scale=8, font="Times New Roman", style='white')
    f, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharex=False,
                                 sharey=True, figsize=(18, 9),
                                 gridspec_kw={'width_ratios': [.25, .7],
                                              'wspace': 0, 'hspace': 0})
    gg = sns.heatmap(fun_anc.astype(float), ax=ax1, fmt=".0f", annot=True,
                     cmap="Oranges", xticklabels=True, yticklabels=True,
                     linewidth=2, robust=True, cbar=False,
                     annot_kws={"size": 8})
    gg.tick_params(axis='both', which='major', labelsize=9)
    hh = sns.heatmap(fun_par.astype(float), ax=ax2, fmt=".0f", annot=True,
                     cmap='Blues', xticklabels=True, yticklabels=True,
                     linewidth=2, robust=True, cbar=False,
                     annot_kws={"size": 8})
    hh.tick_params(axis='both', which='major', labelsize=9)
    plt.gcf()
    plt.setp(ax2.get_yticklabels(), visible=False)
