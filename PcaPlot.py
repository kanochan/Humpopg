
import re
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from functools import reduce
from itertools import filterfalse
from collections import Counter
from adjustText import adjust_text
from Bio.Cluster import pca
from DataFunc import hap_made, get_mutation, get_haplogroup, get_initial, remove_star, remove_bracket, index_num
from HapFileIO import data_path

class PcaPlot(object):

    def __init__(self, population_file):
        self.population_data = pd.read_table(population_file, sep='\t', header=0, encoding='utf-8')
        self.logger = logging.getLogger()

    def hap_read(self, hap_data, extra_file, hap_type, filter, classify):
        hap_df = pd.merge(hap_data, self.population_data, on='SampleID', how='inner').iloc[:, 0:5]
        key_info = pd.read_csv(data_path + 'ISOGG_keymutation_index.csv', sep=',', header=0, encoding='utf-8')
        hap_df = hap_made(hap_df, hap_type, filter)

        key_info['mutation'] = key_info['mutation'].map(remove_bracket)
        if classify == 'rough':
            hap_df['Haplogroup'] = hap_df['Haplogroup'].map(lambda x: x[0])
            hap_df['Mutation'] = [key_info.at[key_info['haplogroup'].tolist().index(i), 'mutation'] for i in hap_df['Haplogroup']]
        elif classify == 'smart':
            initial = set(hap_df['Haplogroup'].map(lambda x: x[0]))
            common_list = []
            for i in initial:
                initial_group = sorted(list(filterfalse(lambda x: x[0] != i, hap_df['Haplogroup'])), key=lambda x: len(x))
                while True:
                    most_common = Counter(initial_group).most_common(1)[0][0]
                    initial_group = [most_common if x.startswith(most_common) and len(x) > len(most_common) else x for x in initial_group]
                    hap_num = len(most_common)-1
                    for j in range(len(most_common)-1):
                        hap = most_common[0:hap_num]
                        if hap in initial_group:
                            initial_group = [hap if x == most_common else x for x in initial_group]
                            most_common = hap
                        hap_num-=1
                    common_list.append(most_common)
                    initial_group = list(filterfalse(lambda x: x.startswith(most_common), initial_group))
                    if len(initial_group) == 0:
                        break
            common_list = sorted(common_list, key=lambda x: len(x), reverse=True)
            for i in range(len(hap_df['Haplogroup'])):
                for j in common_list:
                    if hap_df.at[i, 'Haplogroup'].startswith(j):
                        hap_df.at[i, 'Haplogroup'] = j
            hap_df['Mutation'] = [key_info.at[key_info['haplogroup'].tolist().index(i), 'mutation'] for i in hap_df['Haplogroup']]
        elif classify == 'accurate':
            pass

        hap_df['Haplogroup'] = hap_df['Haplogroup'].map(get_initial)
        hap_df['Haplogroup'] = hap_df['Haplogroup']+ '-' +hap_df['Mutation']
        del hap_df['Mutation']

        if extra_file:
            hap_extra = pd.read_csv(extra_file, sep='\t', header=0, encoding='utf-8')
            if hap_extra.iloc[:, 0:2].tolist == ['SampleID', 'Haplogroup']:
                hap_df = pd.merge(hap_df, hap_extra.iloc[:, 0:2], how='outer')
            else:
                self.logger.warning('Extra file format error, not to input.')

        return hap_df

    def pca_plot(self, hap_df, pca_output_file, freq=False):
        if self.population_data.columns[0:2].tolist() == ['SampleID', 'Population']:
            hap_df = pd.merge(hap_df, self.population_data, on='SampleID').drop('SampleID', axis=1)
            num_list = []
            num = 0
            for i in sorted(list(set(self.population_data['Population'])), key=self.population_data.drop('SampleID', axis=1).drop_duplicates()['Population'].tolist().index):
                num_list.append((num, num + self.population_data['Population'].tolist().count(i)))
                num = num + self.population_data['Population'].tolist().count(i)
            self.population_data = self.population_data.drop('SampleID', axis=1).drop_duplicates()
            if hap_df.columns.size > 2:
                hap_df = hap_df.iloc[:, :2]
        elif self.population_data.columns[0:2].tolist() == ['Population', 'Number']:
            population_df = pd.DataFrame(columns=['Population', 'Index'])
            num_list = []
            num = 0
            for i in range(self.population_data.index.size):
                population, index = self.population_data.loc[i, :].tolist()[:2]
                index = index.split(', ')
                index = list(reduce(lambda x,y: x+y, list(map(index_num, index))))
                population = population*len(index)
                num_list.append((num, num + len(index)))
                num = num + len(index)
                population_df_sub = pd.DataFrame({'Population':population, 'Index':index})
                population_df = pd.merge(population_df, population_df_sub, how='outer')

            hap_df['Index'] = hap_df.index
            hap_df = pd.merge(hap_df, population_df, on='Index', how='outer').drop(['SampleID', 'Index'], axis=1)
            hap_df = hap_df.sort_values(by='Population')
        else:
            self.logger.error('Population file format error, please check your file.')
            sys.exit()

        population_freq = pd.DataFrame(index=sorted(list(set(hap_df['Population'])), key=self.population_data['Population'].tolist().index),
                                       columns=sorted(list(set(hap_df['Haplogroup']))),
                                       dtype=np.float)
        hap_df = pd.Series(hap_df['Haplogroup'].tolist(), index=[hap_df['Population']])

        for i in population_freq.index:
            hap_list = hap_df[i].value_counts().index.tolist()
            freq_list = hap_df[i].value_counts().tolist()
            sum_num = sum(freq_list)
            for j in range(len(hap_list)):
                population_freq.at[i, hap_list[j]] = freq_list[j]/sum_num
        population_freq = population_freq .fillna(0)
        if freq:
            population_freq.to_csv(pca_output_file.rstrip('png')+'freq', index=True, sep='\t', encoding='utf-8')

        matrix = population_freq.as_matrix()
        columnmean, coordinates, components, eigenvalues = pca(matrix)
        color_list = ['slategrey']*self.population_data.index.size
        shape_list = ['o']*self.population_data.index.size
        if 'Color' in self.population_data.columns:
            color_list = self.population_data['Color'].tolist()
        if 'Shape' in self.population_data.columns:
            shape_list = self.population_data['Shape'].tolist()
        if 'Label' in self.population_data.columns:
            label_list = self.population_data['Label'].tolist()

        self.logger.info('Start to plot PCA.')
        figure, ax = plt.subplots(figsize = (12,8))
        plt.xlabel('PC1  %.2f%%'% (eigenvalues[0]/sum(eigenvalues)*100))
        plt.ylabel('PC2  %.2f%%'% (eigenvalues[1]/sum(eigenvalues)*100))

        if 'label_list' in vars():
            plot_list = []
            for i in range(self.population_data.index.size):
                plot = list(zip(color_list, shape_list, label_list))[i]
                if plot in plot_list:
                    plt.plot(coordinates[i, 0], coordinates[i, 1], shape_list[i], color=color_list[i], markersize=5, alpha=0.7)
                else:
                    plt.plot(coordinates[i, 0], coordinates[i, 1], shape_list[i], color=color_list[i], label=label_list[i], markersize=5, alpha=0.7)
                plot_list.append(plot)
            plt.legend()
        else:
            for i in range(self.population_data.index.size):
                plt.plot(coordinates[num_list[i][0]:num_list[i][1], 0], coordinates[num_list[i][0]:num_list[i][1], 1], shape_list[i], color=color_list[i], markersize=5, alpha=0.7)

        #plt.plot(coordinates[num_list[1][0]:num_list[1][1], 0], coordinates[num_list[1][0]:num_list[1][1], 1], shape_list[0], color=color_list[2], label=label_list[0], markersize=5, alpha=0.7)

        if population_freq.index.size < 40:
            texts = [plt.text(x=coordinates[i, 0], y=coordinates[i, 1] ,s=population_freq.index.tolist()[i]) for i in range(len(matrix))]
            adjust_text(texts)
        plt.savefig(pca_output_file, dpi=100)
        self.logger.info('PCA plot finished.')
        print('-' * 80 + '\n')
