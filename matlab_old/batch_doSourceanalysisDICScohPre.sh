#!/bin/sh

ssh compute-0-0 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[1];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-20 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[2];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-2 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[3];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[5];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-4 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[6];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-5 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[8];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-6 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[9];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-17 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[10];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-8 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[11];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
#ssh compute-0-11 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[12];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[13];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[14];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-16 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[15];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-18 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[16];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-21 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[18];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-0-25 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[19];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
ssh compute-1-6 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[20];freqno=[16:31];batch_doSourceanalysisDICScohPre;exit" exit' &
