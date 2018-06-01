#!/bin/sh

ssh compute-1-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=1;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=2;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=3;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=5;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=6;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=7;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=8;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=9;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=10;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=11;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=12;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=13;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=14;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=15;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=16;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=17;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=18;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=19;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;freqno=[19:2:49];batch_doSourceanalysisDICSglm;exit" exit' &
