#!/bin/sh

ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=1;batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[2 3];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5 6];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[7 8];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[9 10];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11 12];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[13 14];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[15 16];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];batch_doFreqanalysisTFRperi;exit" exit' &
ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];batch_doFreqanalysisTFRperi;exit" exit' &
