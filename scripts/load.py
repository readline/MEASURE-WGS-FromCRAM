#!/usr/bin/env python
import os,sys
import pandas as pd

def samplesheet(sspath):
    '''
    Input:
    Sample sheet with following columns:
    SAMPLE	PATIENT	GENDER  EVENT	SID	CRAM	MD5
    
    Output:
    samples:  dict={SAMPLE:{PATIENT:***, GENDER:M/F/U, EVENT:TOD/CR, SID:***, CRAM:***, MD5:***}, ...}
    patients: dict={PATIENT:{EVENT:SAMPLE}, ...}
    tasks:    dict={
        germline:[SAMPLE, SAMPLE, ...], 
        somatic_ss:[SAMPLE, SAMPLE, ...], 
        somatic_paired:[[SAMPLE|TOD, SAMPLE|CR], ...]
    }
    '''
    df = pd.read_csv(sspath, sep='\t', index_col=0)
    patients = {}
    samples = {}
    tasks = {'germline':[],'somatic_ss':[], 'somatic_paired':[]}
    for s in df.index:
        tasks['germline'].append(s)
        p = df.loc[s,'PATIENT']
        e = df.loc[s,'EVENT']
        if p not in patients:
            patients[p] = {}
        samples[s] = df.loc[s].to_dict()
        patients[p][e] = s
    for p in patients:
        if 'TOD' in patients[p]:
            if 'CR' in patients[p]:
                tasks['somatic_paired'].append([patients[p]['TOD'], patients[p]['CR']])
            else:
                tasks['somatic_ss'].append(patients[p]['TOD'])
    return samples, patients, tasks

def test_samplesheet():
    sspath='%s/../sample_sheet.txt'%(os.path.dirname(os.path.abspath(__file__)))
    print(samplesheet(sspath))
    
if __name__ == '__main__':
    test_samplesheet()