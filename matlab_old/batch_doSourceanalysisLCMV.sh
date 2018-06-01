#!/bin/sh

ssh compute-1-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3 5];batch_doSourceanalysisLCMV;exit" exit' &
ssh compute-1-6 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[6:10];batch_doSourceanalysisLCMV;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:15];batch_doSourceanalysisLCMV;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[16:20];batch_doSourceanalysisLCMV;exit" exit' &
