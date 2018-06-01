function datafile = vismot_subject2datafile(subject)

% helper function to extract the full-path to the data

if ischar(subject)
  subject = vismot_subjinfo(subject);
end

datafile = cell(numel(subject.runnames),1);
for k = 1:numel(subject.runnames)
  datafile{k} = fullfile(subject.rawpath,...
                      subject.name,...
                      subject.scanname,...
                      subject.sessionname,...
                      subject.runnames{k},...
                      subject.datafile);

end
 
