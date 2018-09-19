
import re
import sys
import gzip
import time
import logging
import itertools
import numpy as np
import pandas as pd

from collections import deque
from DataFunc import missing_count

class VcfProcess(object):
    """Read vcf file and filter female individuals."""

    def __init__(self, input_vcf):
        self.input_vcf = input_vcf
        self.logger = logging.getLogger()

    # open vcf_file and filter female individuals
    def open_vcf(self):
        if self.input_vcf.endswith('.vcf'):
            vcf_file = open(self.input_vcf, 'r+')
        elif self.input_vcf.endswith('vcf.gz'):
            vcf_file = gzip.open(self.input_vcf, 'rt')

        for i in itertools.count(start=0, step=1):
            vcf_line = vcf_file.readline()
            if vcf_line.startswith('#CHROM'):
                head_num = i
                break

        return head_num

    # filter female in vcf file
    def filter_female(self, head_num, missing):
        vcf_info = pd.read_table(self.input_vcf, header=head_num, sep='\t', encoding='utf-8', dtype='object')

        chr_type = len(set(vcf_info['#CHROM'].tolist()))
        if chr_type > 1:
            vcf_info = vcf_info[~vcf_info['#CHROM'].map(lambda x: 'Y' in x)]

        # get the number of individuals
        ind_info = vcf_info.columns.tolist()[9:]
        ind_num = len(ind_info)
        self.logger.info('There are %d individuals in this vcf file.' % ind_num)

        # count missing rate and delete female individuals
        del_list = []
        female_num = 0
        ind_turn = 0
        pos_num = vcf_info.index.size

        try:
            for i in ind_info:
                missing_num = vcf_info[i].map(missing_count).tolist().count('missing')
                missing_rate = missing_num / pos_num
                percent = float(ind_turn + 1) * 100 / float(ind_num)
                sys.stdout.write('Filtering felamle individuals…… %.2f' % percent)
                sys.stdout.write('%\r')
                sys.stdout.flush()
                ind_turn += 1
                if missing_rate > missing:
                    female_num += 1
                    del_list.append(i)
            vcf_info.drop(del_list, inplace=True, axis=1)
        except TypeError:
            self.logger.error('Vcf format error.')
            sys.exit()

        ind_info = vcf_info.columns.tolist()[9:]
        left_num = vcf_info.columns.size - 9
        self.logger.info('%d female individuals filtered and %d individuals left.' % (female_num, left_num))

        return vcf_info, ind_info

    # match reference info to vcf info, and count haplogroup of each individual
    def match_info(self, ref_info, vcf_info, ind_info, pos_build):

        info_matched = pd.merge(ref_info, vcf_info,
                                left_on=pos_build,
                                right_on='POS',
                                how='inner',
                                sort=False).drop('POS', axis=1)
        key_hap_num = len(set(ref_info.loc[0:2859, ][pos_build]) & set(vcf_info['POS']))

        ind_hap_dict = {}
        ind_pos_dict = {}
        ind_keyhap_dict = {}
        nohap_ind = []
        ind_num = 0

        for i in ind_info:
            hap_list = deque([])
            pos_list = deque([])
            keyhap_list = deque([])
            for j in range(len(info_matched.loc[:, i])):
                if not info_matched.at[j, i][0] == '.':
                    pos_list.append(info_matched.at[j, 'haplogroup'])
                else:
                    continue

                if not ('ins' in info_matched.at[j, 'alt'] or 'del' in info_matched.at[j, 'alt']):
                    if info_matched.at[j, 'alt'] == info_matched.at[j, 'REF'] and info_matched.at[j, i][0] == '0':
                        hap_list.append(info_matched.at[j, 'haplogroup'])
                        info_matched.at[j, i] = 'matched'
                        if j < key_hap_num:
                            keyhap_list.append(info_matched.at[j, 'haplogroup'])
                    else:
                        alt_list = info_matched.at[j, 'ALT'].split(',')
                        if len(alt_list) > 1:
                            if info_matched.at[j, 'alt'] in alt_list:
                                alt_num = str(alt_list.index(info_matched.at[j, 'alt']) + 1)
                                if info_matched.at[j, i][0] == alt_num:
                                    hap_list.append(info_matched.at[j, 'haplogroup'])
                                    info_matched.at[j, i] = 'matched'
                                    if j < key_hap_num:
                                        keyhap_list.append(info_matched.at[j, 'haplogroup'])
                            else:
                                continue
                        elif len(alt_list) == 1:
                            if info_matched.at[j, 'alt'] == info_matched.at[j, 'ALT'] and info_matched.at[j, i][0] == '1':
                                hap_list.append(info_matched.at[j, 'haplogroup'])
                                info_matched.at[j, i] = 'matched'
                                if j < key_hap_num:
                                    keyhap_list.append(info_matched.at[j, 'haplogroup'])

                elif info_matched.at[j, 'alt'] == 'ins':
                    alt_list = info_matched.at[j, 'ALT'].split(',')
                    if len(alt_list) > 1:
                        ins_list = list(map(lambda x: len(info_matched.at[j, 'REF']) < len(x), alt_list))
                        ins_loc_list = np.argwhere(ins_list == True).shape
                        print((ins_list,ins_loc_list,info_matched.at[j, i][0]))
                        if True in ins_list and int(info_matched.at[j, i][0]) in ins_loc_list:
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])
                    elif len(alt_list) == 1:
                        if len(info_matched.at[j, 'ALT']) > len(info_matched.at[j, 'REF']) and info_matched.at[j, i][0] == '1':
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])

                elif re.search(r'ins \d+ bp', info_matched.at[j, 'alt']):
                    ins_bp = info_matched.at[j, 'alt'][4]
                    alt_list = info_matched.at[j, 'ALT'].split(',')
                    if len(alt_list) > 1:
                        ins_list = list(map(lambda x: str(len(x) - len(info_matched.at[j, 'REF'])), alt_list))
                        ins_loc_list = np.argwhere(ins_list == ins_bp).shape
                        if ins_bp in ins_list and int(info_matched.at[j, i][0]) in ins_loc_list:
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])
                    elif len(alt_list) == 1:
                        if len(info_matched.at[j, 'ALT']) - len(info_matched.at[j, 'REF']) == ins_bp and info_matched.at[j, i][0] == '1':
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])

                elif info_matched.at[j, 'alt'] == 'del':
                    alt_list = info_matched.at[j, 'ALT'].split(',')
                    if len(alt_list) > 1:
                        del_list = list(map(lambda x: len(info_matched.at[j, 'REF']) > len(x), alt_list))
                        del_loc_list = np.argwhere(del_list == True).shape
                        if True in del_list and int(info_matched.at[j, i][0]) in del_loc_list:
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])
                    elif len(alt_list) == 1:
                        if len(info_matched.at[j, 'ALT']) < len(info_matched.at[j, 'REF']) and info_matched.at[j, i][0] == '1':
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                            if j < key_hap_num:
                                keyhap_list.append(info_matched.at[j, 'haplogroup'])

                elif re.search(r'del -\d+ bp', info_matched.at[j, 'alt']):
                        del_bp = info_matched.at[j, 'alt'][5]
                        alt_list = info_matched.at[j, 'ALT'].split(',')
                        if len(alt_list) > 1:
                            del_list = list(map(lambda x: str(len(info_matched.at[j, 'REF']) - len(x)), alt_list))
                            del_loc_list = np.argwhere(del_list == del_bp).shape
                            print((del_list,del_loc_list,info_matched.at[j, i][0]))
                            if del_bp in del_list and int(info_matched.at[j, i][0]) in del_loc_list:
                                hap_list.append(info_matched.at[j, 'haplogroup'])
                                info_matched.at[j, i] = 'matched'
                                if j < key_hap_num:
                                    keyhap_list.append(info_matched.at[j, 'haplogroup'])
                        elif len(alt_list) == 1:
                            if len(info_matched.at[j, 'REF']) - len(info_matched.at[j, 'ALT']) == del_bp and info_matched.at[j, i][0] == '1':
                                hap_list.append(info_matched.at[j, 'haplogroup'])
                                info_matched.at[j, i] = 'matched'
                                if j < key_hap_num:
                                    keyhap_list.append(info_matched.at[j, 'haplogroup'])



            percent = float(ind_num + 1) * 100 / float(len(ind_info))
            sys.stdout.write('Matching vcf information to reference index…… %.2f' % percent)
            sys.stdout.write('%\r')
            sys.stdout.flush()
            ind_num+=1
            if len(hap_list) == 0:
                del info_matched[i]
                nohap_ind.append(i)
                continue
            else:
                ind_hap_dict[i] = hap_list
                ind_pos_dict[i] = pos_list
                ind_keyhap_dict[i] = keyhap_list

        time.sleep(0.05)
        self.logger.info('Matching vcf information to reference index is completed.')

        return info_matched, ind_hap_dict, ind_pos_dict, ind_keyhap_dict, nohap_ind, key_hap_num


class InpProcess(object):
    """Read vcf file and filter female individuals."""

    def __init__(self, input_inp):
        self.input_inp = input_inp
        self.logger = logging.getLogger()

    # open vcf_file and filter female individuals
    def open_inp(self):
        if self.input_inp.endswith('.inp'):
            inp_file = open(self.input_inp, 'r+')
        elif self.input_inp.endswith('inp.gz'):
            inp_file = gzip.open(self.input_inp, 'rt')

        blank_num = 0
        for i in itertools.count(start=0, step=1):
            inp_line = inp_file.readline()
            if inp_line == '\n':
                blank_num += 1
            if inp_line.startswith('dbSNP'):
                head_num = i - blank_num
                break

        return head_num

    # filter female in vcf file
    def filter_female(self, head_num, missing):
        inp_info = pd.read_table(self.input_inp, header=head_num, sep='\t', encoding='utf-8')
        pos_num = inp_info.index.size
        Chr, Pos, Allele = inp_info.columns.tolist()[1], inp_info.columns.tolist()[2], inp_info.columns.tolist()[4]

        chr_type = len(set(inp_info[Chr].tolist()))
        if chr_type > 1:
            inp_info = inp_info[~inp_info[Chr].map(lambda x: 'Y' in x)]

        # get the number of individuals
        ind_info = inp_info.columns.tolist()[5:]
        ind_num = len(ind_info)
        self.logger.info('There are %d individuals in this inp file.' % ind_num)

        # count missing rate and delete female individuals
        del_list = []
        female_num = 0
        ind_turn = 0

        try:
            for i in ind_info:
                missing_num = inp_info[i].tolist().count('U')
                missing_rate = missing_num / pos_num
                percent = float(ind_turn + 1) * 100 / float(ind_num)
                sys.stdout.write('Filtering felamle individuals…… %.2f' % percent)
                sys.stdout.write('%\r')
                sys.stdout.flush()
                ind_turn += 1
                if missing_rate > missing:
                    female_num += 1
                    del_list.append(i)
            inp_info.drop(del_list, inplace=True, axis=1)
        except:
            self.logger.error('Inp format error.')
            sys.exit()

        ind_info = inp_info.columns.tolist()[5:]
        left_num = inp_info.columns.size - 5
        self.logger.info('%d female individuals filtered and %d individuals left.' % (female_num, left_num))

        return inp_info, ind_info, Pos, Allele

    def match_info(self, ref_info, inp_info, ind_info, pos_build, Pos, Allele):

        info_matched = pd.merge(ref_info, inp_info, left_on=pos_build, right_on='Position', how='inner',
                                sort=False).drop('Position', axis=1)
        key_hap_num = len(set(ref_info.loc[0:2859, ][pos_build]) & set(inp_info[Pos]))

        ind_hap_dict = {}
        ind_pos_dict = {}
        ind_keyhap_dict = {}
        nohap_ind = []
        ind_num = 0
        for i in ind_info:
            hap_list = deque([])
            pos_list = deque([])
            keyhap_list = deque([])
            ind_alt_info = info_matched.loc[:, i].tolist()
            for j in range(len(ind_alt_info)):
                if not info_matched.at[j, i] == 'U':
                    pos_list.append(info_matched.at[j, 'haplogroup'])
                else:
                    continue
                fst_al, scd_al = info_matched.at[j, Allele].split('/')
                if info_matched.at[j, i] == 'H':
                    hap_list.append(info_matched.at[j, 'haplogroup'])
                    info_matched.at[j, i] = 'matched'
                    if j < key_hap_num:
                        keyhap_list.append(info_matched.at[j, 'haplogroup'])
                elif info_matched.at[j, 'alt'] == fst_al and info_matched.at[j, i] == 'A':
                    hap_list.append(info_matched.at[j, 'haplogroup'])
                    info_matched.at[j, i] = 'matched'
                    if j < key_hap_num:
                        keyhap_list.append(info_matched.at[j, 'haplogroup'])
                elif info_matched.at[j, 'alt'] == scd_al and info_matched.at[j, i] == 'B':
                    hap_list.append(info_matched.at[j, 'haplogroup'])
                    info_matched.at[j, i] = 'matched'
                    if j < key_hap_num:
                        keyhap_list.append(info_matched.at[j, 'haplogroup'])
                else:
                    continue
            percent = float(ind_num + 1) * 100 / float(len(ind_info))
            sys.stdout.write('Matching vcf information to reference index…… %.2f' % percent)
            sys.stdout.write('%\r')
            sys.stdout.flush()
            ind_num += 1
            if len(hap_list) == 0:
                del info_matched[i]
                nohap_ind.append(i)
                continue
            else:
                ind_hap_dict[i] = list(hap_list)
                ind_pos_dict[i] = list(pos_list)
                ind_keyhap_dict[i] = list(keyhap_list)

        time.sleep(0.05)
        self.logger.info('Matching inp information to reference index is completed.')

        return info_matched, ind_hap_dict, ind_pos_dict, ind_keyhap_dict, nohap_ind, key_hap_num
