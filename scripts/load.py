#!/usr/bin/env python
import os,sys
import pandas as pd

def samplesheet(sspath):
    '''
    Input:
    Sample sheet with following columns:
    SAMPLE	PATIENT	GENDER  EVENT	SID	RUN	CRAM	MD5
    
    Output:
    dic_patient_to_eventsamples:    dict={PATIENT:{EVENT:SAMPLE}, ...}
    dic_sample_to_runs:             dict={SAMPLE:[SAMPLE.RUN, SAMPLE.RUN, ...]}
    dic_run:                        dict={SAMPLE.RUN:{SAMPLE:***, PATIENT:***, GENDER:M/F/U, EVENT:TOD/CR, SID:***, CRAM:***, MD5:***}, }
    tasks:    dict={
        germline:[SAMPLE, SAMPLE, ...], 
        somatic_ss:[SAMPLE, SAMPLE, ...], 
        somatic_paired:[[SAMPLE|TOD, SAMPLE|CR], ...]
    }
    '''
    df = pd.read_csv(sspath, sep='\t')
    dic_patient_to_eventsamples = {}
    dic_sample_to_runs = {}
    dic_run = {}
    tasks = {'germline':[],'somatic_ss':[], 'somatic_paired':[]}
    for i in df.index:
        s = df.loc[i,'SAMPLE']
        r = df.loc[i,'RUN']
        sr = '%s.%s'%(s,r)
        tasks['germline'].append(s)
        p = df.loc[i,'PATIENT']
        e = df.loc[i,'EVENT']
        if p not in dic_patient_to_eventsamples:
            dic_patient_to_eventsamples[p] = {}
        if s not in dic_sample_to_runs:
            dic_sample_to_runs[s] = []
        dic_run[sr] = df.loc[i].to_dict()
        dic_patient_to_eventsamples[p][e] = s
        dic_sample_to_runs[s].append(sr)
    for p in dic_patient_to_eventsamples:
        if 'TOD' in dic_patient_to_eventsamples[p]:
            if 'CR' in dic_patient_to_eventsamples[p]:
                tasks['somatic_paired'].append([dic_patient_to_eventsamples[p]['TOD'], dic_patient_to_eventsamples[p]['CR']])
            else:
                tasks['somatic_ss'].append(dic_patient_to_eventsamples[p]['TOD'])
    return dic_patient_to_eventsamples, dic_sample_to_runs, dic_run, tasks

def test_samplesheet():
    sspath='%s/../sample_sheet.txt'%(os.path.dirname(os.path.abspath(__file__)))
    for i in samplesheet(sspath):
        print(i)
        print()
    
if __name__ == '__main__':
    test_samplesheet()