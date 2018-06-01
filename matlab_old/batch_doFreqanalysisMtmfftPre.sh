#!/bin/sh

#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=0;foilim=[0 20];batch_doFreqanalysisMtmfftPre;exit" exit' &

#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-0-31 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-0-31 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-0-31 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
ssh compute-1-9 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=3.75;foilim=[8 40];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:2];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[3 5];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[6:7];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:9];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[10:11];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[12:13];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:15];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[16:17];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[18:19];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
#ssh compute-1-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[20];smoothing=7.5;foilim=[35 100];batch_doFreqanalysisMtmfftPre;exit" exit' &
