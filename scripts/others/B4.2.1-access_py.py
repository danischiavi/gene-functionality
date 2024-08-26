#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:16:25 2020

@author: bikash
"""

import sys
import os
import secrets
import string
import tempfile
from subprocess import run, PIPE, DEVNULL
import argparse
from distutils.spawn import find_executable

class Accessibility:
    
    '''Accessibility
    '''

    def __init__(self, seq, winsize=None, subseg=None):
        self.seq = seq
        self.winsize = winsize
        if self.winsize is None:
            self.winsize = len(self.seq)
        self.subseg = subseg
        if self.subseg is None:
            self.subseg = len(self.seq)
        self.all_values = None #Result from RNAplfold


    @staticmethod
    def accession_gen():
        '''Random accession numbers
        '''
        rand_string = ''.join(secrets.choice(string.ascii_uppercase + \
                                            string.digits) for _ in range(10))
        accession = '>' + rand_string + '\n'
        return accession, rand_string


    def run(self):
        '''
        Run RNAplfold
        '''
        #make temp directory for files
        tmp = os.path.join(tempfile.gettempdir(), 'plfold')
        try:
            os.makedirs(tmp)
        except FileExistsError:
            pass

        #Generate random accession for sequence
        seq_accession, rand_string = Accessibility.accession_gen()
        input_sequence = seq_accession + self.seq

        if self.subseg is None:
            position = len(self.seq)
        else:
            position = self.subseg

        #Run RNAplfold
        all_args = ['RNAplfold', '-W', str(self.winsize), '-u', \
                    str(position), '-O']
        run(all_args, stdout=PIPE, stderr=DEVNULL, input=input_sequence, \
            cwd=tmp, encoding='utf-8')

        out1 = '/' + rand_string + '_openen'
        out2 = '/' + rand_string + '_dp.ps'

        try:
            self.all_values = pd.read_csv(tmp+out1, sep='\t', skiprows=2, \
                                          header=None)

            open_en = self.all_values.loc[len(self.seq) - 1, position]
        except Exception as exp:
            print('THIS happened!\n' + str(exp) + '\n')
            open_en = None

        os.remove(tmp+out1)
        os.remove(tmp+out2)

        return open_en

def check_arg(args=None):
    '''arguments.
    '''
    parser = argparse.ArgumentParser(prog='Program',
                     description='Description',
                     epilog='Author')
    parser.add_argument('-v', '--version',
                    action='version',
                    version='%(prog)s ' + '1',
                    help="Show program's version number and exit.")
    parser.add_argument('-s', '--sequence',
                    type=str,
                    help='Input sequence',
                    required=True)
    parser.add_argument('-W', '--window_size',
                    type=int,
                    help='Window size',
                    required=False)
    parser.add_argument('-u', '--ulength',
                type=int,
                help='Compute the mean probability that regions of length' + \
                ' 1 to a given length are unpaired. ',
                required=False)
    results = parser.parse_args(args)
    return (results.sequence,
            results.window_size,
            results.ulength)


def main():
    '''Main func
    '''
    new_ = Accessibility(s, W, u)
    result = new_.run()
    print(result)

if __name__ == '__main__':
    s, W, u = check_arg(sys.argv[1:])

    try:
     import pandas as pd
    except ImportError:
        print('\nPandas not found!\nUse this:\n\npip3 install --user pandas\n')
        exit()
    if find_executable('RNAplfold') is None:
        print('\nInstall RNAplfold first!')
        exit()
    main()
