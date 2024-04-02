#!/usr/bin/env python
import os,sys
import pandas as pd

def samplesheet(sspath):
    '''
    Input:
    Sample sheet with following columns:
    SAMPLE	PATIENT	GENDER  EVENT	SID	CRAM	MD5
    
    Output:
    sample: dict=[SAMPLE:{PATIENT:***, GENDER:M/F/U, EVENT:TOD/CR, SID:***, CRAM:***, MD5:***}]
    '''
    df = pd.read_csv(sspath, sep='\t', index_col=0)
    patient = {}
    sample = {}
    tasks = {'germline':[],'somatic_ss':[], 'somatic_paired':[]}
    for s in df.index:
        tasks['germline'].append(s)
        p = df.loc[s,'PATIENT']
        e = df.loc[s,'EVENT']
        if p not in patient:
            patient[p] = {}
        sample[s] = df.loc[s].to_dict()
        patient[p][e] = s
    for p in patient:
        if 'TOD' in patient[p]:
            if 'CR' in patient[p]:
                tasks['somatic_paired'].append([patient[p]['TOD'], patient[p]['CR']])
            else:
                tasks['somatic_ss'].append(patient[p]['TOD'])
    return sample, patient, tasks

def test_samplesheet():
    sspath='%s/../sample_sheet.txt'%(os.path.dirname(os.path.abspath(__file__)))
    print(samplesheet(sspath))
    
if __name__ == '__main__':
    test_samplesheet()