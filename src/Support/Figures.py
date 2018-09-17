from wordcloud import WordCloud, STOPWORDS
from PIL import Image
import os
import seaborn as sns
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from Support.Additional import get_grey_colour
from Support.LoadData import make_timely
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
plt.style.use('seaborn-ticks')
plt.rcParams["font.family"] = "Helvetica"
mpl.rcParams.update(mpl.rcParamsDefault)


def wordcloud_figure(abstract_count, output_file):
    """ Make the double helix word cloud: just a bit of fun."""
    words_array = []
    with open(abstract_count, 'r', errors='replace') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['word'].lower() not in STOPWORDS:
                if len(row['word']) > 3:
                    words_array.append(
                        (row['word'].upper(), float(row['num_words'])))
    mask = Image.new('RGBA', (8000, 4467))
    icon = Image.open(os.path.abspath(
                      os.path.join('__file__', '../..', 'Data', 'Support',
                                   'doublehelix_mask.png'))).convert('RGBA')
    mask.paste(icon, icon)
    mask = np.array(mask)
    wc = WordCloud(background_color='white', max_words=1250, mask=mask,
                   max_font_size=5000)
    wc.generate_from_frequencies(dict(words_array))
    wc.recolor(color_func=get_grey_colour)
    wc.to_file(output_file)
    plt.figure(figsize=(16, 8))
    plt.imshow(wc, interpolation='bilinear')
    plt.axis('off')
    plt.show()


def gwas_growth(output_file, Cat_Studies, Cat_Ancestry,
                Cat_Ancestry_groupedbyN):
    """ Plot the growth of GWAS over time (Figure 1)"""
    plt.style.use('seaborn-ticks')
    plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['axes.linewidth'] = 0.75
    yearlist = []
    yearquarterlist = []
    for year in range(2007, 2018):
        yearlist.append(str(year))
        for quarter in ['Q1', 'Q2', 'Q3', 'Q4']:
            yearquarterlist.append(str(year) + quarter)
    variables = ['N ≤ 5,000', '5,001 ≤ N ≤ 50,000', '50,001 ≤ N ≤ 100,000',
                 '100,001 ≤ N', 'N', 'Associations', 'Journals Printing GWAS',
                 '# Diseases Studied']
    df_years, df_quarters = make_timely(variables,
                                        yearlist,
                                        yearquarterlist,
                                        Cat_Studies,
                                        Cat_Ancestry,
                                        Cat_Ancestry_groupedbyN)
    plt.figure(figsize=(15, 10))
    axA = plt.subplot(2, 1, 1)
    ax0variables = ['N ≤ 5,000', '5,001 ≤ N ≤ 50,000',
                    '50,001 ≤ N ≤ 100,000', '100,001 ≤ N']
    ax0 = df_quarters[ax0variables].plot(kind='bar', stacked=True, ax=axA,
                                         color=['#e41a1c', '#377eb8',
                                                '#4daf4a', '#ff7f00'],
                                         alpha=0.6, edgecolor='k')
    sns.despine(top=True, right=True, ax=ax0)
    ax0.set_ylabel('Number of Study Accessions', fontsize=12)
    ax0.tick_params(labelsize=10)
    ax0.legend(fontsize=12, loc='upper left')
    axB = plt.subplot(2, 2, 3)
    ax1a = df_years[['Associations']].plot(ax=axB, color='#e41a1c', alpha=0.75,
                                           rot=90, marker='o', linewidth=1.5,
                                           markersize=8,
                                           label='Associations Discovered',
                                           markeredgecolor='k',
                                           markeredgewidth=0.5)
    ax1b = axB.twinx()
    ax1b.plot(df_years[['N']], color='#377eb8', marker='s', markersize=7,
              linewidth=1.5, markeredgecolor='k', markeredgewidth=0.5)
    ax1a.set_ylabel('Number of Associations Discovered', fontsize=12)
    ax1b.set_ylabel('Number of Study Participants Analyzed', fontsize=12)
    ax1b.grid(False)
    axB.plot(0, 0, '-r', color='#377eb8', marker='s', markersize=7,
             markeredgecolor='k', markeredgewidth=0.5)
    axB.legend(['Associations (left)', 'Participants (right)'],
               fontsize=12, loc='upper left')
    ax1a.tick_params(labelsize=10)
    ax1b.tick_params(labelsize=10)
    plt.axis('tight')
    axC = plt.subplot(2, 2, 4)
    axtest = axC.twinx()
    ax_2a = df_years[['Journals Printing GWAS']].plot(kind='bar',
                                                      ax=axC,
                                                      position=1,
                                                      color='#377eb8',
                                                      legend=False,
                                                      width=0.35,
                                                      alpha=0.75,
                                                      edgecolor='k')
    ax_2b = df_years[['# Diseases Studied']].plot(kind='bar',
                                                  ax=axtest,
                                                  position=0,
                                                  color='#ff7f00',
                                                  width=0.35,
                                                  legend=False,
                                                  alpha=0.75,
                                                  edgecolor='k')
    ax_2a.set_ylabel('Unique Number of Journals Publishing GWAS', fontsize=12)
    ax_2b.set_ylabel('Unique Number of Diseases Studied', fontsize=12)
    ax_2b.grid(False)
    axC.plot(np.nan, '#377eb8', linewidth=4)
    axC.plot(np.nan, '#ff7f00', linewidth=4)
    axC.legend(['Journals (left)', 'Diseases (right)'],
               fontsize=12, loc='upper left')
    ax_2a.margins(1, 0.5)
    ax_2a.tick_params(labelsize=10)
    ax_2b.tick_params(labelsize=10)
    plt.axis('tight')
    plt.tick_params(axis='both', which='minor', labelsize=10)
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')


def choropleth_map(df, input_series, cmap, output_path):
    """ fairly generic function to make a choropleth map
    feed in either 'N' or 'ParticipationPerPerson': Population adjusted is
    just for robustness"""

    cm = plt.get_cmap(cmap)
    df['scheme'] = [
        cm(df[input_series][i] / df[input_series].max()) for i in df.index]
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, facecolor='#deeff5', frame_on=False)
    m = Basemap(lon_0=0, projection='robin', resolution='i')
    m.drawmapboundary(color='k', linewidth=0.075)
    m.drawcountries(color='k', linewidth=0.025)
    m.drawmeridians(np.arange(0, 360, 60), labels=[False, False, False, True],
                    color='#bdbdbd', dashes=[6, 6], linewidth=0.1, fontsize=10)
    m.drawparallels(np.arange(-90, 90, 30), labels=[True, False, False, False],
                    color='#bdbdbd', dashes=[6, 6], linewidth=0.1, fontsize=10)
    m.readshapefile(os.path.abspath(os.path.join('__file__',
                                                 '../..',
                                                 'data',
                                                 'ShapeFiles',
                                                 'ne_10m_admin_0_countries')),
                    'units', color='#444444', linewidth=.075)
    for info, shape in zip(m.units_info, m.units):
        country = info['NAME_CIAWF']
        if country not in df.index:
            color = 'w'
        else:
            color = df.loc[country]['scheme']
        patches = [Polygon(np.array(shape), True)]
        pc = PatchCollection(patches)
        pc.set_facecolor(color)
        ax.add_collection(pc)
    norm = Normalize()
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cm)
    if input_series == 'N':
        mapper.set_array(df[input_series] / 1000000)
    else:
        mapper.set_array(df[input_series])
    clb = plt.colorbar(mapper, shrink=0.75)
    clb.ax.tick_params(labelsize=10)
    if input_series == 'N':
        clb.ax.set_title('Number of\n People (m)', y=1.02, fontsize=11)
    elif input_series == 'ParticipationPerPerson':
        clb.ax.set_title('Participations\nPer Recruitment',
                         y=1.02, fontsize=11)
    plt.savefig(output_path, bbox_inches='tight')
    plt.show()


def plot_heatmap(funder_ancestry, funder_parent, output_path):
    ''' Build the funder heatmaps '''
    sns.set(font_scale=1, font='Arial', style='white')
    f, (ax1, ax2) = plt.subplots(nrows=1,
                                 ncols=2,
                                 sharex=False,
                                 sharey=True,
                                 figsize=(18, 11),
                                 gridspec_kw={'width_ratios': [.25, .7],
                                              'wspace': .05, 'hspace': 0})
    gg = sns.heatmap(funder_ancestry.astype(float),
                     ax=ax1,
                     fmt='.0f',
                     annot=True,
                     cmap='Oranges',
                     xticklabels=True,
                     yticklabels=True,
                     linewidth=0.1,
                     linecolor='k',
                     robust=True,
                     cbar=False,
                     annot_kws={'size': 10})
    gg.tick_params(axis='both', which='major', labelsize=12)
    gg.set_xlabel('Ancestry', fontsize=12)
    hh = sns.heatmap(funder_parent.astype(float),
                     ax=ax2,
                     fmt='.0f',
                     annot=True,
                     cmap='Blues',
                     xticklabels=True,
                     yticklabels=True,
                     linewidth=0.1,
                     linecolor='k',
                     robust=True,
                     cbar=False,
                     annot_kws={'size': 11})
    hh.tick_params(axis='both', which='major', labelsize=12)
    hh.set_xlabel('Broad EFO Category', fontsize=12)
    plt.gcf()
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.tight_layout(h_pad=5)
    plt.savefig(output_path, bbox_inches='tight')


def plot_bubbles(output_path, Cat_Ancestry,
                 Broad_Ancestral_NoNR, countriesdict):
    """ This makes the Broader Ancestry bubble plot (Figure 2) """
    fig = plt.figure(figsize=(12, 6), dpi=800)
    ax = fig.add_subplot(1, 1, 1)
    for obs in Cat_Ancestry.index:
        for key, value in countriesdict.items():
            if Cat_Ancestry['Broader'][obs].strip() == key:
                ax.plot_date(x=pd.to_datetime(Cat_Ancestry['Dates'][obs]),
                             y=Cat_Ancestry['N'][obs],
                             color=value,
                             marker='.',
                             label='the data',
                             alpha=0.6,
                             markersize=Cat_Ancestry['N'][obs] / 6500)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=24))
    ax.set_xlim(pd.Timestamp('2007-01-01'), pd.Timestamp('2018-12-31'))
    ax.set_ylim(0, Cat_Ancestry['N'].max()+50000)
    ax.yaxis.tick_right()
    ax.set_ylabel('Individual Sample Sizes')
    ax.yaxis.set_label_position("right")
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    handles, labels = ax.get_legend_handles_labels()
    inset_fig = inset_axes(ax, width=3.25,
                           height=1.5, loc=3,
                           bbox_to_anchor=(0.275, .6),
                           bbox_transform=ax.figure.transFigure)
    barplot = Broad_Ancestral_NoNR[['% Discovery',
                                    '% Replication']].transpose()
    barplot = barplot.plot(kind='bar', colors=['#4daf4a', '#984ea3',
                                               '#e41a1c', '#377eb8',
                                               '#ff7f00', '#ffff33'],
                           alpha=0.6, edgecolor='k', linewidth=1,
                           ax=inset_fig, width=0.9)
    barplot.set_yticklabels('')
    barplot.set_ylabel('')
    barplot.set_xlabel('')
    barplot.set_ylim(0, 120)
    barplot.get_yaxis().set_visible(False)
    for p in barplot.patches:
        barplot.annotate(str(round(p.get_height(), 2)) + '%',
                         (p.get_x() + p.get_width() / 2., p.get_height() + 20),
                         ha='center', va='center', rotation=90, fontsize=8)
    plt.legend(loc='center right', bbox_to_anchor=(0.025, 0.675),
               facecolor='white', edgecolor='k', frameon=False,
               prop={'size': 9})
    sns.despine(top=True, left=True, right=False, ax=ax, offset=10)
    sns.despine(top=True, right=True, left=True, offset=5, ax=inset_fig)
    plt.savefig(output_path, bbox_inches='tight')


def boxswarm_plot(dataframe, ranks, output_path):
    """ plot the box swarm plots: these are relegated mostly to the notebook
    and dont feature in the paper at present """
    fig = plt.figure(figsize=(14, 7), dpi=800)
    ax1 = fig.add_subplot(1, 1, 1)
    p = sns.boxplot(data=dataframe,
                    y='EFO term',
                    ax=ax1,
                    order=ranks,
                    x='isfemale',
                    whis=1.5,
                    palette="vlag",
                    fliersize=0)
    sns.swarmplot(y='EFO term',
                  x='isfemale',
                  data=dataframe,
                  order=ranks,
                  size=2.25,
                  color="white",
                  edgecolor="black",
                  linewidth=0.5,
                  ax=ax1)
    p.set_ylabel('')
    p.set_xlim(0, 1.05)
    p.set_xlabel('Female/Male Ratio Per Paper Across Mapped Terms',
                 fontsize=12)
    p.tick_params(labelsize=12)
    plt.setp(p.spines.values(), color='k')
    plt.setp([p.get_xticklines(), p.get_yticklines()], color='k')
    sns.despine(bottom=False, left=False, offset=10, ax=ax1)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
