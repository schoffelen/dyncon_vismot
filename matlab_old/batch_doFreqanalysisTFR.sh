#!/bin/sh

ssh compute-0-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=1;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=2;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=3;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=5;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=6;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=7;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=8;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=9;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=10;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=11;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=12;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=13;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-6 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=14;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-6 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=15;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=16;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=17;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=18;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=19;batch_doFreqanalysisTFR;exit" exit' &
ssh compute-0-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;batch_doFreqanalysisTFR;exit" exit' &
