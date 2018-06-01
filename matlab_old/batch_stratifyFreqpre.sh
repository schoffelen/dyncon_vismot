#!/bin/sh

ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=12;batch_stratifyFreqpre;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=12;batch_stratifyFreqpre;exit" exit' &
