#!/bin/sh

#ssh compute-0-0 'matlab -nodisplay -r "ft_defaults;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(13)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-20 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(14)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-2 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(15)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(16)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-4 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(17)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-23 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(1)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-6 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(2)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-17 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(3)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-8 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(4)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-11 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(5)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(6)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(7)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-16 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(8)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-18 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(9)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-21 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(10)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-0-25 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(11)).name;collect4Dsmoothpow(sx);exit" exit' &
ssh compute-1-6 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjlist=[1:3 5:6 8:16 18:20];subjinfo;sx=SUBJ(subjlist(12)).name;collect4Dsmoothpow(sx);exit" exit' &
