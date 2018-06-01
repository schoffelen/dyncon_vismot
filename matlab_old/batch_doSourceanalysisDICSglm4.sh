#!/bin/sh

ssh compute-0-1 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[1];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-2 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[2];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[3];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-5 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[5];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-6 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[6];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-7 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[7];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-8 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[8];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-9 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[9];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-10 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[10];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-11 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[11];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-13 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[12];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-14 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[13];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-15 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[14];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-16 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[15];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-17 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[16];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-18 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[17];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-19 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[18];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-20 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[19];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
#ssh compute-0-21 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[20];freqno=[3:5 7:12 19:2:29 35:2:45];batch_doSourceanalysisDICSglm4;exit" exit' &
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
