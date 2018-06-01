#!/bin/sh

ssh compute-0-1 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[1];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-20 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[2];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-2 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[3];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-3 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[5];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-28 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[6];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-5 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[8];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-27 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[9];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-17 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[10];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-8 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[11];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-11 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[12];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-12 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[13];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-14 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[14];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-16 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[15];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-18 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[16];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-21 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[18];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-25 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[19];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
ssh compute-0-26 'matlab -nodisplay -r "fieldtripdefs;addpath /home/jan/projects/visuomotor/matlab;subjno=[20];freqno=[1:31];batch_doSourceanalysisDICSpowPstPrevious;exit" exit' &
