#!/bin/sh

ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=5:6;batch_doRealignment;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=7:8;batch_doRealignment;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=9:10;batch_doRealignment;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=11:12;batch_doRealignment;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=13:14;batch_doRealignment;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=15:16;batch_doRealignment;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=17;batch_doRealignment;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=18;batch_doRealignment;exit" exit' &
ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=19;batch_doRealignment;exit" exit' &
ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=20;batch_doRealignment;exit" exit' &
