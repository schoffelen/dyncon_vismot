#!/bin/sh

#ssh compute-0-30 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[1];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-4 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[2];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-1 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[3];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-17 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[5];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-22 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[6];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
#ssh compute-0-25 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[7];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-23 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[8];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-31 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[9];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-0 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[10];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[11];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[12];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-21 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[13];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[14];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-7 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[15];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-26 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[16];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
#ssh compute-0-24 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[17];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-29 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[18];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-27 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[19];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
ssh compute-0-16 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[20];freqno=[1:17];batch_doSourceanalysisDICSpre2pow;exit" exit' &
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
