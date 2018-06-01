#!/bin/sh

ssh compute-0-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-13 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-17 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-0-16 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=4;foilim=[6 40];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-0 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:2];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-1 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[3 5];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-2 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[6:7];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-3 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:9];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[10:11];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[12:13];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:15];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-10 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[16:17];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[18:19];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
ssh compute-1-11 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[20];smoothing=12;foilim=[36 100];batch_doFreqanalysisBalancePrePostResplock;exit" exit' &
