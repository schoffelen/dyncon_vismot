%this line loads in the .mat file
load('nameofdatafile');

%this line puts the data matrix in the variable dat
dat = left.trial;

%this line saves the variable dat, as an ascii file and the file will
%be called 'nameofnewdatafileleft' (of course the name can be whatever
save('nameofnewdatafileleft', 'dat', '-ascii');

%this step can then be repeated for the other conditions
dat = right.trial;
save('nameofnewdatafileright', 'dat', '-ascii');
    
%if you want to have access to the channel names
label = right.label;
save('nameoflabelfile', 'label', '-ascii');
