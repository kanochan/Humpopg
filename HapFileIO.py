
import os
import sys
import time
import logging
import argparse
import pandas as pd

from DataFunc import remove_note_tag

start = time.clock()
data_path = os.path.split(os.path.realpath(__file__))[0] + '/data/'
work_path = os.getcwd() + '/'

tree_file = open(data_path + 'haplogroup_tree_2018.txt', 'r+')
tree_info = list(map(lambda x: x.replace('\n', ''), tree_file.readlines()))
tree_file.close()

def hap_parser():
    parser = argparse.ArgumentParser('hap', description='Y-haplogroup tracker @ver-1.0 BY ChenHao')
    parser.add_argument('hap',
                        help='Run to generate hap file and visualize haplogroup data.')
    parser.add_argument('-i', '--input', required=True, type=str, nargs=2,
                        help='input: Input vcf or inp file.')
    parser.add_argument('-o', '--output', required=False, type=str, nargs=1,
                        help='output: Output haplogroup file.')
    parser.add_argument('-b', '--build', required=False, type=int, nargs=1, default=37, choices=[37, 38],
                        help='build: Set the build version of reference genome.')
    parser.add_argument('-m', '--missing', required=False, type=float, nargs=1, default=0.4,
                        help='missing rate: Set missing rate to filter female samples, default is 0.4')
    parser.add_argument('-n', '--notes', required=False, action='store_true',
                        help='notes: Set whether keep the haplogroup attatching (Notes), defualt not setting.')
    parser.add_argument('-f', '--filter', required=False, action='store_true',
                        help='filter: filter samples without key mutation.')
    parser.add_argument('--pca', required=False, type=str, nargs='+',
                        help='pcaplot: Output a pca plot.')
    parser.add_argument('--phylo', required=False, type=str, nargs='*',
                        help='phylotreeplot: output a phylotree plot.')

    args = parser.parse_args()

    return args

def judge_inputfile_type(input_info, work_path):
    file_type, input_file = input_info
    if file_type == 'vcf':
        if os.path.isfile(input_file):
            print('\n' + '-' * 80)
            filename = os.path.split(input_file)[1]
            if input_file.endswith('.vcf') or input_file.endswith('.vcf.gz'):
                remove_indel = False
            else:
                print('%s is not a vcf file, please check it again.' % filename)
                print('-' * 80 + '\n')
                sys.exit()
            return filename, input_file, remove_indel
        elif os.path.isfile(work_path + input_file):
            print('\n' + '-' * 80)
            input_file = work_path + input_file
            filename = os.path.split(input_file)[1]
            if input_file.endswith('.vcf') or input_file.endswith('.vcf.gz'):
                remove_indel = False
            else:
                print('%s is not a vcf file, please check it again.' % filename)
                print('-' * 80 + '\n')
                sys.exit()
            return filename, input_file, remove_indel
        else:
            print('\n' + '-' * 80)
            print('%s not found, please check your filename.' % input_file)
            print('-' * 80 + '\n')
            sys.exit()
    elif file_type == 'inp':
        if os.path.isfile(input_file):
            print('\n' + '-' * 80)
            filename = os.path.split(input_file)[1]
            if input_file.endswith('.inp') or input_file.endswith('.inp.gz'):
                remove_indel = True
            else:
                print('%s is not a inp file, please check it again.' % filename)
                print('-' * 80 + '\n')
                sys.exit()
            return filename, input_file, remove_indel
        elif os.path.isfile(work_path + input_file):
            print('\n' + '-' * 80)
            input_file = work_path + input_file
            filename = os.path.split(input_file)[1]
            if input_file.endswith('.inp') or input_file.endswith('.inp.gz'):
                remove_indel = True
            else:
                print('%s is not a inp file, please check it again.' % filename)
                print('-' * 80 + '\n')
                sys.exit()
            return filename, input_file, remove_indel
        else:
            print('\n' + '-' * 80)
            print('%s not found, please check your filename.' % input_file)
            print('-' * 80 + '\n')
            sys.exit()
    else:
        print('Input command error, please input as format: -i [vcf|inp] [file]')
        sys.exit()

def set_outputfile(output_file, work_path):
    if os.path.exists(os.path.split(output_file)[0]):
        if not output_file.endswith('.hap'):
            output_file = output_file + '.hap'
        return output_file
    elif os.path.split(output_file)[0] == '':
        output_file = work_path + output_file
        if not output_file.endswith('.hap'):
            output_file = output_file + '.hap'
        return output_file
    else:
        print('%s is not a correct path, please check it again.' % output_file)
        print('-' * 80 + '\n')
        sys.exit()
    return output_file

def judge_missing(missing_value):
    if missing_value < 0 or missing_value > 1:
        print('Input command error, please input as format: -i [vcf|inp] [file]')
        sys.exit()


def set_log(input_file, filename, remove_indel, args_log):
    logger = logging.getLogger()
    logger.setLevel(level=logging.INFO)

    handler = logging.FileHandler(input_file + '.log', mode='w')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] - [%(levelname)s]: %(message)s')
    handler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    logger.addHandler(handler)
    logger.addHandler(console)

    log_info = ['Run Date: ' + time.asctime(time.localtime(time.time())),
                'Filename: %s' % filename,
                'Female missing rate: %.2f' % args_log.missing,
                'Not to count haplogroup with (Notes).',
                'Not to remove indel data',
                'Keep samples without key mutation.']

    if args_log.notes:
        log_info[2] = 'Count haplogroup with (Notes).'
    if remove_indel:
        log_info[3] = 'Remove indel data.'
    if args_log.filter:
        log_info[4] = 'Filter samples without key mutation.'

    for i in log_info:
        logger.info(i)


def get_reference_index(data_path, notes):
    ref_info = pd.read_csv(data_path + 'ISOGG_haplogroup_index.csv', header=0, sep=',', encoding='utf-8', dtype='object')

    if notes:
        ref_info['haplogroup'] = ref_info['haplogroup'].map(remove_note_tag)
    else:
        ref_info = ref_info[~ref_info['haplogroup'].map(lambda x: '(Notes)' in x)]

    return ref_info

def get_pca_file(pca_info, work_path):
    logger = logging.getLogger()
    classify_types = ['rough', 'smart', 'accurate']
    hap_types = ['key', 'final']
    classify = 'smart'
    hap_type = 'key'
    extra_file = None
    classify_flag, hap_type_flag, extra_info_flag = (0,0,0)
    if os.path.isfile(pca_info[0]) and pca_info[0].endswith('.info'):
        population_file = pca_info[0]
    elif os.path.isfile(work_path + pca_info[0]) and pca_info[0].endswith('.info'):
        population_file = work_path + pca_info[0]
    else:
        logger.error('%s is not a correct population file, please assure file format correct.' % pca_info[0])
        sys.exit()
    if len(pca_info) > 1:
        for i in pca_info[1:]:
            if i in classify_type:
                if classify_flag == 0:
                    classify = i
                    classify_flag+=1
                elif classify_flag != 0:
                    logger.warning('Input more than one haplogroup classifying arguments and only count the first.')
            elif i in hap_type:
                if hap_type_flag == 0:
                    hap_type = i
                    hap_type_flag+=1
                elif classify_flag != 0:
                    logger.warning('Input more than one haplogroup type arguments and only count the first.')
            elif os.path.isfile(i):
                if extra_info_flag == 0:
                    if i.endswith('.extra'):
                        extra_file = i
                        extra_info_flag+=1
                    else:
                        logger.warning('%s is not a correct extra file, please assure file format correct.' % i)
                elif extra_info_flag == 1:
                    logger.warning('Input more than one extra file and only count the first.')
            elif os.path.isfile(work_path + i):
                if extra_info_flag == 0:
                    if i.endswith('.extra'):
                        extra_file = i
                        extra_info_flag+=1
                    else:
                        logger.warning('%s is not a correct extra file, please assure file format correct.' % i)
                elif extra_info_flag == 1:
                    logger.warning('Input more than one extra file and only count the first.')
            else:
                logger.error('%s is not a correct arguments.' % i)
                sys.exit()

        return population_file, extra_file, classify, hap_type
