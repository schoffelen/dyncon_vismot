#!/bin/sh

#ssh compute-0-18 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-0-12 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-0-13 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-0-14 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-0-15 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-0-16 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=4;foilim=[6 40];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[1:3];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-7 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[5:7];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[8:10];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-8 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[11:13];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[14:16];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-5 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[17 18];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
#ssh compute-1-4 'matlab -nodisplay -r "addpath /home/jan/projects/visuomotor/matlab;subjno=[19 20];smoothing=12;foilim=[36 100];batch_doFreqanalysis_aligned;exit" exit' &
