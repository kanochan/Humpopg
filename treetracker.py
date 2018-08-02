import os
import gc
import re
import sys
import time
import gzip
import getopt
import pandas as pd
from functools import reduce
from collections import Counter

class CmdDispose(object):
    '''Set command of this script.'''
    def __init__(self):
        super(CmdDispose, self).__init__()
        self.start = time.clock()
        self.scriptpath = os.path.split(os.path.realpath(__file__))[0] + '/'
        self.workpath = os.getcwd() + '/'
        self.argv = sys.argv[1:]
        self.ncmd = False
        self.indelcmd = False
        self.mvalue = 0.4
        self.tvalue = 0.6
        self.bvalue = 37
        self.loginfo = ['Run Date: ' + time.asctime(time.localtime(time.time())),
                        'Filename: None',
                        'Not count haplogroup with (Notes).',
                        'Not remove the indel stastics.',
                        'Reference genome version: 37',
                        'Female missing rate: 0.4',
                        'Treetrack rate: 0.6']
        self.logfile = open(self.workpath + 'treetracker.log', 'w')

    def cmd_control(self):
        inputfile = ''
        if not '-h' in self.argv:
            if not '-i' in self.argv:
                print('\n' + '-'*80)
                print('Input command not found, please input \'-i filename\'.')
                print('-'*80 + '\n')
                sys.exit()
        try:
            opts, args = getopt.getopt(self.argv, 'b:hi:m:nt:', ['bvalue=', 'ifile=', 'mvalue=', 'tvalue=', 'rmindel'])
        except:
            print('\n' + '-'*80)
            print('Input format error, please try again.')
            print('-'*80 + '\n')
            sys.exit()
        for opt, arg in opts:
            if opt in ('-h', '--help'):
                self.cmd_h()
                sys.exit()
            elif opt in ('-i', '--ifile'):
                inputfile = arg
                self.cmd_i(inputfile)
            elif opt == '-n':
                self.cmd_n()
            elif opt == '--rmindel':
                self.cmd_rmindel()
            elif opt in ('-b', '--bvalue'):
                self.bvalue = arg
                self.cmd_b()
            elif opt in ('-m', '--mvalue'):
                self.mvalue = arg
                self.cmd_m()
            elif opt in ('-t', '--tvalue'):
                self.tvalue = arg
                self.cmd_t()

        for i in self.loginfo:
            print(i)
            self.logfile.write(i+'\n')
        self.logfile.close()

        return self.start, self.scriptpath, self.workpath, inputfile, self.ncmd, self.indelcmd, self.bvalue, self.mvalue, self.tvalue

    def cmd_h(self):
        print('\n' + '-'*80)
        print(  '  Y-haplogroup strain @ver-1.0 BY ChenHao\n'
              + '  -b build: Set the build version of reference genome.\n'
              + '  -i inputfile: input one or more vcf files.\n'
              + '  -h help: See the usage and command help for this program.\n'
              + '  -m missing: Set the gender missing rate of vcf files.\n'
              + '     gender missing rate is used to filter female individuals.\n'
              + '     it means the rate of missing point in vcf files.\n'
              + '     the defualt value is 0.4.\n'
              + '  -n notes: Set whether keep the haplogroup attatching (Notes).\n'
              + '     defualt not setting.\n'
              + '  -t treetrack: Set the treetrack rate of output file.\n'
              + '     treetrack rate is used to measure the acuracy of treetrack.\n'
              + '     it means the rate of haplogroup that could be traced in treetrack.\n'
              + '     the higher this rate, the more accurate the treetrack, but the trace may be shorter.\n'
              + '  -rmindel Remove indel stastics.\n'
              + '           defualt not setting')
        print('-'*80 + '\n')

    def cmd_i(self, inputfile):
        if os.path.isfile(inputfile):
            print('\n' + '-'*80)
            self.loginfo[1]= 'Filename: ' + os.path.split(inputfile)[1]
            inputfile = self.workpath + os.path.split(inputfile)[1]
            if not (inputfile.endswith('.vcf') or inputfile.endswith('.vcf.gz')):
                print(os.path.split(inputfile)[1] + ' is not a vcf file, please check it again.')
                self.logfile.write('Inputfile not a vcffile')
                self.logfile.close()
                print('-'*80 + '\n')
                sys.exit()
        elif os.path.isfile(self.workpath + inputfile):
            print('\n' + '-'*80)
            self.loginfo[1]= 'Filename: ' + os.path.split(inputfile)[1]
            inputfile = self.workpath + inputfile
            if not (inputfile.endswith('.vcf') or inputfile.endswith('.vcf.gz')):
                print(os.path.split(inputfile)[1] + ' is not a vcf file, please check it again.')
                self.logfile.write('Inputfile not a vcffile')
                self.logfile.close()
                print('-'*80 + '\n')
                sys.exit()
        else:
            print('\n' + '-'*80)
            print(inputfile + ' not found, please check your filename.')
            self.logfile.write('Inputfile not exists')
            self.logfile.close()
            print('-'*80 + '\n')
            sys.exit()

    def cmd_n(self):
        self.ncmd = True
        self.loginfo[2] = 'Count haplogroup with (Notes).'

    def cmd_rmindel(self):
        self.indelcmd = True
        self.loginfo[3] = 'Remove the indel stastics.'

    def cmd_b(self):
        try:
            self.bvalue = int(self.bvalue)
            if self.bvalue != 37 or self.bvalue != 38:
                print(self.bvalue + ' is not a correct value, please input the reference genome version number of 37 or 38.')
                self.logfile.write('Reference genome version number error.')
                self.logfile.close()
                print('-'*80 + '\n')
                sys.exit()
            else:
                self.loginfo[4] = 'Reference genome version: ' + str(self.bvalue)
        except:
            print(str(self.bvalue) + ' is not a correct value, please input the reference genome version number of 37 or 38.')
            self.logfile.write('Reference genome version number error.')
            self.logfile.close()
            print('-'*80 + '\n')
            sys.exit()

    def cmd_m(self):
        try:
            self.mvalue = float(self.mvalue)
            if self.mvalue < 0 or self.mvalue > 1:
                print(self.mvalue + ' is not a correct value, please input a value between 0-1.')
                self.logfile.write('Female missing rate number error.')
                self.logfile.close()
                print('-'*80 + '\n')
                sys.exit()
            else:
                self.loginfo[5] = 'Female missing rate: ' + str(self.mvalue)
        except:
            print(str(self.mvalue) + ' is not a correct value, please input a value between 0-1.')
            self.logfile.write('Female missing rate number error.')
            self.logfile.close()
            print('-'*80 + '\n')
            sys.exit()

    def cmd_t(self):
        try:
            self.tvalue = float(self.tvalue)
            if self.tvalue < 0 or self.tvalue > 1:
                print(self.tvalue + ' is not a correct value, please input a value between 0-1.')
                self.logfile.write('Treetrack rate number error.')
                self.logfile.close()
                print('-'*80 + '\n')
                sys.exit()
            else:
                self.loginfo[6] = 'Treetrack rate: ' + str(self.tvalue)
        except:
            print(str(self.tvalue) + ' is not a correct value, please input a value between 0-1.')
            self.logfile.write('Treetrack rate number error.')
            self.logfile.close()
            print('-'*80 + '\n')
            sys.exit()

class RefDispose(object):
    '''Read reference file and dispose stastics accoring to the options.'''
    def __init__(self, scriptpath):
        super(RefDispose, self).__init__()
        self.scriptpath = scriptpath

    #open reference index
    def open_ref(self, ncmd, indelcmd):
        ref_info = pd.read_table(self.scriptpath + 'ref_index.txt', header = 0, sep = '\t', encoding = 'utf-8')

        if indelcmd == True:
            ref_info = ref_info[~ref_info['haplogroup'].map(lambda x: len(x)>4)]
        else:
            pass

        def note_dispose(x):
            if '(Notes)' in x:
                x = x.replace(' (Notes)', '')
            return x
        if ncmd == True:
            ref_info['haplogroup'] = ref_info['haplogroup'].map(note_dispose)
        else:
            ref_info = ref_info[~ref_info['haplogroup'].map(lambda x: '(Notes)' in x)]

        def alt_dispose(x):
            if '->' in x:
                num = x.index('->') + 2
                x = x[num:]
            return x
        ref_info['alt'] = ref_info['alt'].map(alt_dispose)

        return ref_info

class VcfDispose(object):
    '''Read vcf file and filter female individuals.'''
    def __init__(self, workpath):
        super(VcfDispose, self).__init__()
        self.workpath = workpath

    #open vcf_file and filter female individuals
    def open_vcf(self, inputfile, mvalue):
        if inputfile.endswith('.vcf'):
            vcf_file = open(inputfile, 'r+')
        elif inputfile.endswith('vcf.gz'):
            vcf_file = gzip.open(inputfile, 'rt')
        vcf_info = vcf_file.readlines()
        head_num = 0
        for i in vcf_info:
            if i.startswith('#CHROM'):
                head_num = vcf_info.index(i)
                break
        vcf_file.close()

        vcf_info = pd.read_table(inputfile, header = head_num, sep = '\t', encoding = 'utf-8')
        ind_info = vcf_info.columns.tolist()[9:]

        #get the number of individuals
        ind_num = vcf_info.columns.size - 9
        print('There are %d individuals in this vcf file.'%ind_num)

        #count missing rate and delete female individuals

        del_list = []
        female_num = 0
        ind_num = 0
        pos_num = vcf_info.index.size
        def missing_count(x):
            if './.' in x:
                x = 'missing'
            return x
        try:
            for i in ind_info:
                missing_num = vcf_info[i].map(missing_count).tolist().count('missing')
                missing_rate = missing_num/pos_num
                percent = float(ind_num)*100/float(len(ind_info))
                sys.stdout.write('Filtering felamle individuals…… %.4f'%percent)
                sys.stdout.write('%\r')
                sys.stdout.flush()
                ind_num+=1
                if missing_rate > mvalue:
                    female_num+=1
                    del_list.append(i)
            vcf_info.drop(del_list, inplace=True, axis=1)
        except:
            print('Vcf format error.')

        ind_info = vcf_info.columns.tolist()[9:]
        left_num = vcf_info.columns.size - 9
        print('%d female individuals filtered and %d individuals left.'%(female_num, left_num))

        return vcf_info, ind_info

class InfoInter(object):
    '''Match reference file and vcf file.'''
    def __init__(self, scriptpath):
        super(InfoInter, self).__init__()
        self.scriptpath = scriptpath

    #match reference info to vcf info, and count haplogroup of each individual
    def match_info(self, ref_info, vcf_info, ind_info, bvalue):

        def indel_dispose(x):
            if '..' in x:
                num = x.index('.')
                x = x[:num]
            elif '-' in x:
                num = x.index('-')
                x = x[:num]
            x = int(x)
            return x
        build = 'pos_' + str(bvalue)
        ref_info[build] = ref_info[build].map(indel_dispose)

        info_matched = pd.merge(ref_info, vcf_info, left_on = build, right_on = 'POS', how = 'inner', sort = False)

        ind_hap_dict = {}
        ind_pos_dict = {}
        nohap_ind = []
        ind_num = 0
        for i in ind_info:
            hap_list = []
            pos_list = []
            ind_alt_info = info_matched.loc[ : , i].tolist()
            for j in range(len(ind_alt_info)):
                if not info_matched.at[j, i][0] == '.':
                    pos_list.append(info_matched.at[j, 'haplogroup'])
                else:
                    continue
                if info_matched.at[j, 'alt'] == info_matched.at[j, 'REF'] and info_matched.at[j, i][0] == '0':
                    hap_list.append(info_matched.at[j, 'haplogroup'])
                    info_matched.at[j, i] = 'matched'
                else:
                    if len(info_matched.at[j, 'ALT']) > 1:
                        alt_list = info_matched.at[j, 'ALT'].split(',')
                        if info_matched.at[j, 'alt'] in alt_list:
                            alt_num = str(alt_list.index(info_matched.at[j, 'alt']) + 1)
                            if info_matched.at[j, i][0] == alt_num:
                                hap_list.append(info_matched.at[j, 'haplogroup'])
                                info_matched.at[j, i] = 'matched'
                        else:
                            continue
                    elif len(info_matched.at[j, 'ALT']) == 1:
                        if info_matched.at[j, 'alt'] == info_matched.at[j, 'ALT'] and info_matched.at[j, i][0] == '1':
                            hap_list.append(info_matched.at[j, 'haplogroup'])
                            info_matched.at[j, i] = 'matched'
                        else:
                            continue
            percent = float(ind_num)*100/float(len(ind_info))
            sys.stdout.write('Matching vcf information to reference…… %.4f'%percent)
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

        print('Matching vcf information to reference is completed.')
        return info_matched, ind_hap_dict, ind_pos_dict, nohap_ind

    #count frequency of haplogroup of each individual
    def count_freq(self, ind_dict):

        freq_info = []
        for key in ind_dict.keys():
            arr = ind_dict[key]
            result2 = {}
            for j in set(arr):
                result2[j] = arr.count(j)
            freq = list(result2.items())
            freq.append(key)
            freq_info.append(freq)

        return freq_info

    #dispose frequency information
    def num_inter(self, freq_info_ind, freq_info_sum):
        for i in range(len(freq_info_ind)):
            for j in freq_info_sum[i][:-1]:
                a = 0
                while a < len(freq_info_ind[i])-1:
                    if j[0] == freq_info_ind[i][a][0]:
                        freq_info_ind[i][a] = (freq_info_ind[i][a][0], freq_info_ind[i][a][1],j[1])
                        break
                    else:
                        a+=1
                        continue

        return freq_info_ind

    #search final haplogroup
    def search_string(self, freq_info_ind, tvalue):
        hap_info = []
        a = 0

        while a < len(freq_info_ind):

            final_hap = ''
            hap_list = []
            for num_hap in range(len(freq_info_ind[a][:-1])):
                hap_list.append(freq_info_ind[a][num_hap][0])
            length = 0
            for i in range(len(freq_info_ind[a][:-1])):
                length2 = len(freq_info_ind[a][i][0])
                if freq_info_ind[a][i][0].startswith('A'):
                    hap_list.append('A')
                if ((length2 > length and freq_info_ind[a][i][0][0] in hap_list
                    and not (len(freq_info_ind[a][i][0]) > 1 and freq_info_ind[a][i][0][-1].isupper()))
                    or
                    (length2 == length and freq_info_ind[a][:-1][i][0][0] in hap_list
                    and freq_info_ind[a][i][1] > flag_num
                    and not (len(freq_info_ind[a][i][0]) > 1 and freq_info_ind[a][i][0][-1].isupper()))):
                    temp = freq_info_ind[a][i][0]
                    temp_num = len(temp)-1
                    exist_num = 0
                    for j in temp:
                        if temp_num > 1:
                            hap_str = temp[0:temp_num]
                            if hap_str in hap_list:
                                exist_num+=1
                            temp_num-=1
                    if len(temp) <= 2 and freq_info_ind[a][i][1]/freq_info_ind[a][i][2] > 0.85:
                        length = length2
                        final_hap = freq_info_ind[a][i]
                        flag_num = freq_info_ind[a][i][1]
                    elif len(temp) == 3 and exist_num == 1 and freq_info_ind[a][i][1]/freq_info_ind[a][i][2] > 0.85:
                        length = length2
                        final_hap = freq_info_ind[a][i]
                        flag_num = freq_info_ind[a][i][1]
                    elif len(temp) > 3:
                        track_rate = (exist_num+2)/len(temp)
                        if track_rate > tvalue and freq_info_ind[a][i][1]/freq_info_ind[a][i][2] > 0.85:
                            length = length2
                            final_hap = freq_info_ind[a][i]
                            flag_num = freq_info_ind[a][i][1]

            try:
                if len(final_hap[0]) == 1:
                    single_list= []
                    for hap, n, m in freq_info_ind[a][:-1]:
                        if len(hap) == 1 and n/m > 0.9:
                            single_list.append(hap)
                    if len(single_list) > 1:
                        single_list.sort(reverse = True)
                        hap_info.append([freq_info_ind[a][-1], single_list[0]])
                        a+=1
                        continue
            except:
                print('Treetrack rate is too high to trace haplogroup, please lower the rate.')
                print('-'*80 + '\n')
                sys.exit()

            hap_info.append([freq_info_ind[a][-1], final_hap[0]])
            a+=1

        return hap_info

    #match mutation information to final haplogroup
    def search_mut(self, hap_info, info_matched):
        mut_info = []
        for ind, hap in hap_info:
            temp_info = info_matched[info_matched['haplogroup'].map(lambda x: x == hap)]
            index_list = temp_info.index.tolist()
            for i in index_list:
                if temp_info.at[i, ind] == 'matched':
                    mut_info.append([ind, temp_info.at[i, 'mutation']])
                    break
                else:
                    continue

        return mut_info

class TreeWrite(object):
    '''build tree-track and write info to new file'''
    def __init__(self, scriptpath, hap_info, mut_info):
        super(TreeWrite, self).__init__()
        self.scriptpath = scriptpath
        self.hap_info = hap_info
        self.mut_info = mut_info

    #construct treetrack
    def tree_track(self):
        a = 0
        track_info = []
        while a < len(self.hap_info):
            track = ''
            for num_str in range(len(self.hap_info[a][1])):
                temp = self.hap_info[a][1][0:len(self.hap_info[a][1])-num_str]
                track = track + '->' + temp
            track = track[2:]
            #A型特殊，单独写
            if track[-1] == 'A':
                track_temp = track[:-3]
                track = '->' + track_temp + '->A0T->A00T->A000T->Y-Adam'
                if ('A00a' or 'A00b' or 'A00c') in track:
                    track_temp = track[:4]
                    track = '->' + track_temp + '->A00->A00T->A000T->Y-Adam'
                if 'A000b1' in track:
                    track = '->A000b1->A000b->A000->A000T->Y-Adam'
            track_info.append([self.hap_info[a][0],track])
            a+=1

        tree_file = open(self.scriptpath+ 'haplogroup_tree_2018.txt', 'r+')
        tree_info = tree_file.readlines()
        for i in range(len(tree_info)):
            tree_info[i] = tree_info[i].replace('\n', '')
        tree_file.close()
        tree_info_reversed = list(reversed(tree_info))
        a = 0
        pattern = re.compile(r'[A-Za-z0-9]+')

        while a < len(track_info):
            length = 0
            num = 0
            for j in tree_info:
                if (j[-2] == '\t' and track_info[a][1][-1] == j[-1]):
                    length = j.count('-')
                    for k in tree_info_reversed:
                        if k.count('-')+1 == length and tree_info.index(k) < tree_info.index(j):
                            temp = pattern.search(k)
                            temp = temp.group()
                            track_info[a][1] = track_info[a][1] + '->' + temp
                            length = length-1
                        else:
                            continue
                    track_info[a][1] = '->' + track_info[a][1] + '->' + 'Y-Adam'
            a+=1

        return track_info

    #dispose all the information
    def tree_iter(self, track_info, freq_info_ind):
        for i in range(len(track_info)):
            track_info[i].append(self.mut_info[i][1])

        for j in range(len(freq_info_ind)):
            b = 0
            while b < len(freq_info_ind[j][:-1]):
                if freq_info_ind[j][b][0] in track_info[j][1]:
                    string = '->' + freq_info_ind[j][b][0] + '->'
                    string_new = '->' + freq_info_ind[j][b][0]+ '(' + str(freq_info_ind[j][b][1])+ '/' + str(freq_info_ind[j][b][2])+ ')->'
                    track_info[j][1] = track_info[j][1].replace(string, string_new)
                b+=1
            track_info[j][1] = track_info[j][1][2:]

        for k in range(len(track_info)):
            info_num = track_info[k][1].count('>')
            if info_num == 0:
                track_info[k][1] = 'Could not paint tree-track. It could be caused by the sequencing quality.'
                track_info[k][2] = 'None'

        return track_info

    #write file
    def write_info(self, start, track_info_w, inputfile, nohap_ind):
        if inputfile.endswith('.vcf.gz'):
            inputfile = inputfile.rstrip('.gz')
        outputfile = inputfile.rstrip('vcf') + 'hap'
        print('Start writing……')
        write_file = open(outputfile, 'w')
        header = 'ind' + '\t' + 'mutation' + '\t' + 'final_hap' + '\t' + 'treetrack' + '\n'
        write_file.write(header)

        for i in range(len(track_info_w)):
            final_hap_index = track_info_w[i][1].index('(')
            final_hap = track_info_w[i][1][:final_hap_index]
            line = track_info_w[i][0] + '\t' + track_info_w[i][2] + '\t' + final_hap + '\t' + track_info_w[i][1] + '\n'
            write_file.write(line)

        if len(nohap_ind) != 0:
            write_file.write('\n')
            nohap = reduce(lambda x, y: x + '\t' + y, nohap_ind)
            nohap = nohap + '\n'
            write_file.write(nohap)
            print(str(len(nohap_ind)) + 'individuals exists 0 pos that could match to reference info.')
        print('Finished.')
        write_file.close()

        end = time.clock()
        interval = end - start
        if interval > 60:
            min = int(interval/60)
            second = interval - min*60
            print('Run time: %dmin%.2fs'%(min,second))
        else:
            print('Run time: %.2fs'%(end-start))
        print('-'*80 + '\n')

def main():

    try:
        C = CmdDispose()
        start, scriptpath, workpath, inputfile, ncmd, indelcmd, bvalue, mvalue, tvalue = C.cmd_control()

        R = RefDispose(scriptpath)
        ref_info = R.open_ref(ncmd, indelcmd)

        V = VcfDispose(workpath)
        vcf_info, ind_info = V.open_vcf(inputfile, mvalue)

        gc.collect()

        I = InfoInter(scriptpath)
        info_matched, ind_hap_dict, ind_pos_dict, nohap_ind = I.match_info(ref_info, vcf_info, ind_info, bvalue)
        freq_info_ind = I.count_freq(ind_hap_dict)
        freq_info_sum = I.count_freq(ind_pos_dict)
        freq_info_ind = I.num_inter(freq_info_ind, freq_info_sum)
        hap_info = I.search_string(freq_info_ind, tvalue)
        mut_info  = I.search_mut(hap_info, info_matched)

        T = TreeWrite(scriptpath, hap_info, mut_info)
        track_info = T.tree_track()
        track_info_w = T.tree_iter(track_info, freq_info_ind)
        T.write_info(start, track_info_w, inputfile, nohap_ind)

    except KeyboardInterrupt:
        print('Treetracker stopped.')
        print('-'*80 + '\n')
        sys.exit()

if __name__ == '__main__':
    main()
