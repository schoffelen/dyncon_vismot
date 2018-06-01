subjinfo
if 0,
for k = [16:20]
  subject = SUBJ(k);
  tlckall = doTimelockanalysisPlanar(subject);
  for kk = 1:5
    cd('/analyse/4/Project0030/ascii/');
    fname = [subject.name,'condition',num2str(kk)];
    dat = tlckall(kk).avg;
    label = tlckall(kk).label;
    l = [];
    for mm = 1:length(label)
      l(mm) = str2num(label{mm}(2:end));
    end
    save([fname,'data'],'dat','-ascii');
    save([fname,'channels'],'l','-ascii');
  end
end

batch_analyzeRT
reactiontimes = trimm;
names = {SUBJ([1:3 5:20]).name}';
cd('/analyse/4/Project0030/ascii/');
save('reactiontimes','reactiontimes','-ascii');
names{2}(end+1)=' ';
names = char(names);
save('names','names','-ascii');
end

if 1,
  for slop = [1:3 5:20]
  subject = SUBJ(slop);
  [alltlck1, alltlck2] = doTimelockanalysisPlanar2(subject,0)
  for kk = 1:4
    cnd = [1:4];
    cd('/analyse/4/Project0030/ascii/withreplicates2');
    fname = [subject.name,'condition',num2str(cnd(kk))];
    siz = size(alltlck1(kk).trial);
    dat = reshape(permute(alltlck1(kk).trial,[3 1 2]),[siz(3) siz(1)*siz(2)])';
    label = alltlck1(kk).label;
    nrepl = siz(1);
    l = [];
    for mm = 1:length(label)
      l(mm) = str2num(label{mm}(2:end));
    end
    save([fname,'data'],'dat','-ascii');
    save([fname,'channels'],'l','-ascii');
    save([fname,'nrepl'], 'nrepl', '-ascii'); 
    %fname = [fname,'planar'];
    %siz = size(alltlck2(kk).trial);
    %dat = reshape(permute(alltlck2(kk).trial,[3 1 2]),[siz(3) siz(1)*siz(2)])';
    %label = alltlck2(kk).label;
    %l = [];
    %for mm = 1:length(label)
    %  labtok = tokenize(label{mm},'_');
    %  l(mm) = str2num(labtok{1}(2:end));
    %end
    %save([fname,'data'],'dat','-ascii');
    %save([fname,'channels'],'l','-ascii');
    %save([fname,'nrepl'], 'nrepl', '-ascii'); 
  end
  end
end

if 0,
  subject = SUBJ(13);
  [alltlck1, alltlck2] = doTimelockanalysisPlanar2(subject, 0)
  for kk = 1:4
    cnd = [1 2 3 4];
    cd('/analyse/4/Project0030/ascii/withreplicates');
    fname = [subject.name,'condition',num2str(cnd(kk))];
    siz = size(alltlck1(kk).trial);
    label = alltlck1(kk).label;
    nrepl = siz(1);
    l = [];
    for mm = 1:length(label)
      l(mm) = str2num(label{mm}(2:end));
    end
    save([fname,'channels'],'l','-ascii');
    save([fname,'nrepl'], 'nrepl', '-ascii'); 
    for mm = 1:nrepl
      mm
      dat = squeeze(alltlck1(kk).trial(mm,:,:));
      save([fname,'replicate',num2str(mm, '%03d')],'dat','-ascii');
    end
  end
end
