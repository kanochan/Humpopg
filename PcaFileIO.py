
import os
import sys
import time
import logging
import argparse
import pandas as pd

start = time.clock()

def pca_parser():
    parser = argparse.ArgumentParser('pca', description='Y-haplogroup tracker @ver-1.0 BY ChenHao')
    parser.add_argument('pca',
                        help='Run to plot pca plot.')
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+',
                        help='input: Input hap or extra file.')
    parser.add_argument('-o', '--output', required=False, type=str, nargs=1,  action='store',
                        help='output: Output haplogroup file.')
    parser.add_argument('-m', '--mode', required=False, type=str, nargs=1, default=['smart'], choices=['rough', 'smart', 'accurate'],
                        help='missing rate: Set missing rate to filter female samples, default is 0.4')
    parser.add_argument('-t', '--type', required=False, type=str, nargs=1, default='key', choices=['key', 'final'],
                        help='notes: Set whether keep the haplogroup attatching (Notes), defualt not setting.')
    parser.add_argument('-f', '--filter', required=False, action='store_true',
                        help='filter: filter samples without key mutation.')
    parser.add_argument('--info', required=True, type=str, nargs=1, action='store',
                        help='input: Input info file.')
    parser.add_argument('--freq', required=False, action='store_true',
                        help='pcaplot: Output a pca plot.')

    args = parser.parse_args()

    return args

def judge_pca_inputfile(input_info, work_path):
    pca_inputfile_extra = None
    if os.path.isfile(input_info[0]) and input_info[0].endswith('.hap'):
        print('\n' + '-' * 80)
        pca_inputfile = input_info[0]
        filename = os.path.split(pca_inputfile)[1]
    elif os.path.isfile(work_path + input_info[0]) and input_info[0].endswith('.hap'):
        print('\n' + '-' * 80)
        pca_inputfile = work_path + input_info[0]
        filename = os.path.split(pca_inputfile)[1]
    else:
        print('\n' + '-' * 80)
        print('%s is not a correct hap file, please assure file format correct.' % input_info[0])
        print('-' * 80 + '\n')
        sys.exit()
    pca_hap_df = pd.read_table(pca_inputfile, sep='\t', header=0, encoding='utf-8')
    if len(input_info) == 2:
        if os.path.isfile(input_info[1]) and input_info[1].endswith('.extra'):
            pca_inputfile_extra = input_info[1]
        elif os.path.isfile(work_path + input_info[1]) and input_info[1].endswith('.extra'):
            pca_inputfile_extra = work_path + input_info[1]
        else:
            print('%s is not a correct extra file, please assure file format correct.' % input_info[0])
    elif len(input_info) > 2:
        print('Please input no more than 2 file.')
    return pca_hap_df, pca_inputfile, pca_inputfile_extra, filename

def set_outputfile(output_file, work_path):
    if os.path.exists(os.path.split(output_file)[0]):
        if not output_file.endswith('.png'):
            output_file = output_file + '.png'
        return output_file
    elif os.path.split(output_file)[0] == '':
        output_file = work_path + output_file
        if not output_file.endswith('.png'):
            output_file = output_file + '.png'
        return output_file
    else:
        print('%s is not a correct path, please check it again.' % output_file)
        print('-' * 80 + '\n')
        sys.exit()
    return output_file

def set_log(input_file, filename, args_log):
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
                'Classifying mode: %s' % args_log.mode[0],
                'Haplogroup type: %s' % args_log.type,
                'Keep samples without key mutation.']

    if args_log.filter:
        log_info[4] = 'Filter samples without key mutation.'
    if args_log.freq:
        log_info.append('Output frequency data.')

    for i in log_info:
        logger.info(i)

def judge_info_file(info_file, work_path):
    if isinstance(info_file, list):
        info_file = info_file[0]
    if os.path.isfile(info_file) and info_file.endswith('.info'):
        population_file = info_file
    elif os.path.isfile(work_path + info_file) and info_file.endswith('.info'):
        population_file = info_file
    else:
        print('%s is not a correct population file, please assure file format correct.' % info_file)
        print('-' * 80 + '\n')
        sys.exit()

    return population_file
