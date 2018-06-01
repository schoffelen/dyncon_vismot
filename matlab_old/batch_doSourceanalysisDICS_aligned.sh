#!/bin/sh

ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=1;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=2;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=3;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=5;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=6;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=7;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=8;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=9;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=10;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=11;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=12;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=13;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=14;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=15;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=16;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=17;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=18;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=19;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;freqno=[1:17 19:2:49];batch_doSourceanalysisDICS_aligned;exit" exit' &
