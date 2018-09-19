
import os
import re
import sys
import time
import logging
import pandas as pd
import FilterMatch

from functools import reduce
from itertools import filterfalse
from collections import Counter
from DataFunc import hap_determine, get_final_hap
from HapFileIO import get_reference_index, data_path, tree_info

def get_matched_info(arguments, input_file):
    file_type, notes, missing, build = arguments.input[0], arguments.notes, arguments.missing, str(arguments.build)
    reference = get_reference_index(data_path, notes)
    if file_type == 'vcf':
        vcf = FilterMatch.VcfProcess(input_file)
        head_num = vcf.open_vcf()
        vcf_info, ind_info = vcf.filter_female(head_num, missing)
        (info_matched,
        ind_hap_dict,
        ind_pos_dict,
        ind_keyhap_dict,
        nohap_ind,
        key_hap_num) = vcf.match_info(reference, vcf_info, ind_info, 'pos_'+build)
    elif file_type == 'inp':
        inp = FilterMatch.InpProcess(input_file)
        head_num = inp.open_inp()
        vcf_info, ind_info = inp.filter_female(head_num, missing)
        (info_matched,
        ind_hap_dict,
        ind_pos_dict,
        ind_keyhap_dict,
        nohap_ind,
        key_hap_num) = inp.match_info(reference, vcf_info, ind_info, 'pos_'+build)

    return (info_matched,
            ind_hap_dict,
            ind_pos_dict,
            ind_keyhap_dict,
            nohap_ind,
            key_hap_num)

class Trackbuild(object):
    """Match Vcf information to reference index."""

    def __init__(self, key_hap_num, data_path):
        self.key_hap_num = key_hap_num
        self.data_path = data_path
        self.logger = logging.getLogger()

    # count frequency of haplogroup of each individual

    @staticmethod
    def count_freq(ind_dict):

        freq_info = {}
        for key in ind_dict.keys():
            freq = Counter(ind_dict[key])
            freq = list(freq.items())
            freq_info[key] = freq

        return freq_info

    # dispose frequency information
    @staticmethod
    def num_inter(freq_info_ind, freq_info_sum):
        for key in freq_info_ind.keys():
            hap_ind, freq_ind = zip(*freq_info_ind[key])
            hap_sum, freq_sum = zip(*freq_info_sum[key])
            freq_table_ind = pd.DataFrame({'hap_ind': hap_ind, 'freq_ind': freq_ind})
            freq_table_sum = pd.DataFrame({'hap_sum': hap_sum, 'freq_sum': freq_sum})
            freq_table = pd.merge(freq_table_ind, freq_table_sum, left_on='hap_ind', right_on='hap_sum',
                                  how='inner').drop('hap_sum', axis=1)
            freq_info_ind[key] = list(zip(freq_table['hap_ind'], freq_table['freq_ind'], freq_table['freq_sum']))

        return freq_info_ind

    # search final haplogroup
    def search_final_hap(self, filter, freq_info_ind, ind_keyhap_dict):
        hap_info = {}
        key_hap_info = {}
        track_rate = 0.8
        common_trunk = ['A0T', 'A1', 'A1b', 'BT', 'CT', 'CF', 'DE', 'GHIJK',
                        'HIJK', 'IJK', 'IJ', 'LT', 'P1', 'NO', 'NO1']

        for key, value in freq_info_ind.items():
            key_hap = '.'
            final_hap = ''
            key_hap_list = sorted(ind_keyhap_dict[key])
            hap_list = list(list(zip(*value))[0])
            if len(key_hap_list) > 0:
                length = 0
                for hap in filterfalse(lambda x: x.startswith('K') or x in common_trunk, key_hap_list):
                    length2 = len(hap)
                    if length2 == 1:
                        continue
                    if length2 >= length:
                        exist_num = 0
                        hap_num = length2 - 1
                        for i in hap:
                            hap_str = hap[0:hap_num]
                            if hap_str in key_hap_list:
                                exist_num += 1
                            hap_num -= 1
                        if exist_num == length2 - 1:
                            length = length2
                            final_hap = hap
                            continue
                        if length2 >= 5 and (exist_num + 1)/length2 > track_rate:
                            length = length2
                            final_hap = hap
                            continue
                        if hap.startswith('A') and (exist_num + 1)/length2 > (track_rate - 0.2):
                            length = length2
                            final_hap = hap

                if len(final_hap) > 1:
                    key_hap = final_hap
                    length = len(final_hap)
                    flag_num = list(list(zip(*value))[1])[hap_list.index(final_hap)]
                    for hap, n, m in filterfalse(lambda x: len(x[0]) != len(final_hap)+1, value):
                        if not hap.startswith(final_hap):
                            continue
                        if hap in common_trunk:
                            continue
                        if m > 5:
                            if not n/m > 0.85:
                                continue
                        elif m <= 5:
                            if not (m-n) == 1:
                                continue
                        if not len(hap) >= length:
                            continue
                        else:
                            length2 = len(hap)
                        if ((length2 > length) or
                                (length2 == length and n > flag_num and hap in key_hap_list)):
                            final_hap = hap
                            flag_num = n

            if len(key_hap_list) == 0:
                self.logger.warning('No key mutation in sample %s.' % (key))
                if filter:
                    self.logger.info('Sample %s filtered.' % (key))
                    continue
                else:
                    final_hap, key_hap = hap_determine(value, final_hap, hap_list, key_hap_list, track_rate, common_trunk)
            elif len(key_hap_list) > 0 and final_hap == '':
                final_hap, key_hap = hap_determine(value, final_hap, hap_list, key_hap_list, track_rate, common_trunk, key_hap)

            if final_hap == '':
                length = 0
                if len(list(filterfalse(lambda x: x in common_trunk, key_hap_list))) > 1:
                    last_hap_list = list(filterfalse(lambda x: x in common_trunk, key_hap_list))
                    for hap in last_hap_list:
                        length2 = len(hap)
                        if length2 >= length:
                            final_hap = hap
                            length = length2
                            key_hap = '*' + final_hap
                else:
                    hap_rate = 0
                    last_hap_list = sorted(value, key=lambda value: value[0])
                    last_hap_list = sorted(last_hap_list, key=lambda last_hap_list: len(value[0]))
                    for hap, n, m in last_hap_list:
                        length2 = len(hap)
                        if len(hap) > 1 and hap[-1].isupper():
                            continue
                        if length2 >= length and n / m >= hap_rate:
                            final_hap = hap
                            hap_rate = n / m
                            length = length2

            hap_info[key] = final_hap
            key_hap_info[key] = key_hap

        return hap_info, key_hap_info

    # match mutation information to final haplogroup
    def search_mut(self, hap_info, key_hap_info, info_matched):
        for ind, hap in hap_info.items():
            temp_info = info_matched[info_matched['haplogroup'].map(lambda x: x == hap)]
            index_list = temp_info.index.tolist()
            for i in index_list:
                if temp_info.at[i, ind] == 'matched':
                    if i < self.key_hap_num:
                        hap_info[ind] = [hap, temp_info.at[i, 'mutation']]
                    else:
                        hap_info[ind] = [hap, '*' + temp_info.at[i, 'mutation']]
                    break
                else:
                    continue

        for ind, key_hap in key_hap_info.items():
            key_hap_info[ind] = [key_hap, '.']
            if key_hap == '.':
                if '*' not in hap_info[ind][1]:
                    key_hap_info[ind] = [hap_info[ind][0], hap_info[ind][1]]
                else:
                    continue
            elif key_hap.startswith('*'):
                key_hap_index = info_matched.loc[:self.key_hap_num, ]['haplogroup'].tolist().index(key_hap[1:])
                key_hap_info[ind] = [key_hap[1:], info_matched.loc[:self.key_hap_num, ].at[key_hap_index, 'mutation']]
            else:
                key_hap_index = info_matched.loc[:self.key_hap_num, ]['haplogroup'].tolist().index(key_hap)
                key_hap_info[ind] = [key_hap, info_matched.loc[:self.key_hap_num, ].at[key_hap_index, 'mutation']]

        self.logger.info('Mutation matched.')
        return hap_info, key_hap_info

    def search_population(self, info_matched, key_hap_info):
        pop_info = pd.read_csv(self.data_path + 'ISOGG_population_index.csv', sep=',', header=0, encoding='utf-8')
        key_mut_info = pd.merge(pop_info, info_matched.loc[:self.key_hap_num, ], on='mutation', how='inner')
        key_mut_info = key_mut_info.iloc[:, :3]

        for key, values in key_hap_info.items():
            if key_mut_info.index.size == 0:
                break
            key_hap_info[key] = [values[0], values[1], '.']
            if values[0] == '.':
                continue
            if values[0].startswith('*'):
                continue
            else:
                hap_num = len(values[0])
                for i in values[0]:
                    hap = values[0][0:hap_num]
                    if hap in key_mut_info['haplogroup_x'].tolist():
                        pop_num = key_mut_info['haplogroup_x'].tolist().index(hap)
                        key_hap_info[key] = [values[0], values[1], key_mut_info.at[pop_num, 'population']]
                        break
                    else:
                        hap_num -= 1
        return key_hap_info

    # construct treetrack
    def tree_track(self, hap_info):
        for key, value in hap_info.items():
            track = ''
            for num_str in range(len(value[0])):
                temp = value[0][0:len(value[0]) - num_str]
                track = track + '->' + temp
            track = track[2:]
            # A is special
            if track[-1] == 'A':
                track_temp = track[:-3]
                track = '->' + track_temp + '->A0T->A00T->A000T->Y-Adam'
                pattern_track = re.compile(r'A00[abc]')
                if re.search(pattern_track, track):
                    track_temp = pattern_track.search(track).group()
                    track = '->' + track_temp + '->A00->A00T->A000T->Y-Adam'
                if 'A000b1' in track:
                    track = '->A000b1->A000b->A000->A000T->Y-Adam'
            hap_info[key][0] = track

        pattern = re.compile(r'[A-Za-z0-9]+')
        for key, value in hap_info.items():
            length = 0
            num = 0
            for i in tree_info:
                if i[-2] == '\t' and value[0][-1] == i[-1]:
                    length = i.count('-')
                    for j in list(reversed(tree_info)):
                        if j.count('-') + 1 == length and tree_info.index(j) < tree_info.index(i):
                            temp = pattern.search(j)
                            temp = temp.group()
                            hap_info[key][0] = hap_info[key][0] + '->' + temp
                            length = length - 1
                        else:
                            continue
                    hap_info[key][0] = '->' + hap_info[key][0] + '->' + 'Y-Adam'

        return hap_info

    # dispose all the information
    @staticmethod
    def tree_iter(hap_info, key_hap_info, freq_info_ind):

        for key, value in freq_info_ind.items():
            for hap, n, m in value:
                if hap in hap_info[key][0]:
                    string = '->' + hap + '->'
                    string_new = '->' + hap + '(' + str(n) + '/' + str(m) + ')->'
                    hap_info[key][0] = hap_info[key][0].replace(string, string_new)
            hap_info[key][0] = hap_info[key][0][2:]
            hap_info[key] = tuple(hap_info[key])

        ind_info = list(hap_info.keys())
        track_info, mut_info = zip(*list(hap_info.values()))
        key_hap_info, key_mut_info, population = zip(*list(key_hap_info.values()))
        info_gathered = [ind_info, key_mut_info, key_hap_info, mut_info, track_info, population]

        return info_gathered

    # write file
    def write_info(self, start, info_gathered, output_file, nohap_ind):
        print('Start writing……')
        hap_data = pd.DataFrame({'SampleID': info_gathered[0],
                                 'Key_mutation': info_gathered[1],
                                 'Final_mutation': info_gathered[3],
                                 'Key_hap': info_gathered[2],
                                 'Final_hap': list(map(get_final_hap, info_gathered[4])),
                                 'Treetrack': info_gathered[4],
                                 'Reference_population': info_gathered[5]},
                                 columns=['SampleID',
                                          'Key_mutation',
                                          'Final_mutation',
                                          'Key_hap',
                                          'Final_hap',
                                          'Treetrack',
                                          'Reference_population'])
        hap_data.to_csv(output_file, sep='\t', index=False, encoding='utf-8')

        if len(nohap_ind) != 0:
            nohap = reduce(lambda x, y: x + '\t' + y, nohap_ind)
            self.logger.info('%d individuals exists 0 pos that could match to reference info.' % len(nohap_ind))

        interval = time.clock() - start
        if interval > 60:
            minute = int(interval/60)
            second = interval - minute * 60
            self.logger.info('Run time: %dmin%.2fs' % (minute, second))
        else:
            self.logger.info('Run time: %.2fs' % (interval))
        print('-' * 80 + '\n')

        return hap_data
