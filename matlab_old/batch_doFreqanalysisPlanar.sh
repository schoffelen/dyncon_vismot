#!/bin/sh

ssh compute-0-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=4;batch_doFreqanalysisPlanar;exit" exit' &
ssh compute-0-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:8];smoothing=4;batch_doFreqanalysisPlanar;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[9:14];smoothing=4;batch_doFreqanalysisPlanar;exit" exit' &
ssh compute-0-13 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[15:17];smoothing=4;batch_doFreqanalysisPlanar;exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[18:20];smoothing=4;batch_doFreqanalysisPlanar;exit" exit' &
#ssh compute-0-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=12;batch_doFreqanalysisPlanar;exit" exit' &
#ssh compute-0-14 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:8];smoothing=12;batch_doFreqanalysisPlanar;exit" exit' &
#ssh compute-0-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[9:14];smoothing=12;batch_doFreqanalysisPlanar;exit" exit' &
#ssh compute-0-6 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[15:17];smoothing=12;batch_doFreqanalysisPlanar;exit" exit' &
#ssh compute-0-16 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[18:20];smoothing=12;batch_doFreqanalysisPlanar;exit" exit' &
