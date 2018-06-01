#!/bin/sh

#ssh compute-1-0 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[1:2];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-12 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[3];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-1 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[5:6];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[7:8];freqno=[3:5];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-11 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[9:10];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[11:12];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-7 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[13:14];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-8 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[15:16];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-9 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[17:18];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[19:20];freqno=[3:5 7:12];batch_doSourceanalysisDICSglm2;exit" exit' &
#alpha = [3:5];
#beta  = [7:12];
#gamma1 = [19:2:29];
#gamma2 = [35:2:45];

#ssh compute-1-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[18];freqno=[19:2:29];batch_doSourceanalysisDICSprepst3;exit" exit' &

#ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;freqno=[41:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=12;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=13;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=14;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=15;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=16;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=17;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=18;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=19;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
#ssh compute-1-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;freqno=[19:2:49];batch_doSourceanalysisDICSprepst3;exit" exit' &
