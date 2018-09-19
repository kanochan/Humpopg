import os
import sys
import time
import logging
import argparse
import pandas as pd

def pca_parser():
    parser = argparse.ArgumentParser('pca', description='Y-haplogroup tracker @ver-1.0 BY ChenHao')
    parser.add_argument('phy',
                        help='Run to plot phylotree plot.')
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
