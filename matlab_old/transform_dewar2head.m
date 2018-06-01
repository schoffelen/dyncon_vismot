function [M] = transform_dewar2head(subject)

cd(subject.pathname);
cd('data');
%load([subject.name,'data']);
load([subject.name,'datatmp']);

nsmp = cellfun(@numel,data1.time);
nsmp = [nsmp cellfun(@numel,data2.time)];
nsmp = [nsmp cellfun(@numel,data3.time)];
nsmp = [nsmp cellfun(@numel,data4.time)];
nsmp = [nsmp cellfun(@numel,data5.time)];

%the following is based on the hard-coded assumption that the last 10
%channels contain headposition data, that is nas,lpa,rpa and rv of the fit
ntrl1 = length(data1.trial);
ntrl2 = length(data2.trial);
ntrl3 = length(data3.trial);
ntrl4 = length(data4.trial);
ntrl5 = length(data5.trial);
ntrl  = [ntrl1 ntrl2 ntrl3 ntrl4 ntrl5];
spnt  = zeros(9,sum(ntrl));
cnt   = 0;
for k=1:ntrl1
    cnt = cnt+1;
    spnt(:, cnt)= data1.trial{k}(end-9:end-1,1);
end
for k=1:ntrl2
    cnt = cnt+1;
    spnt(:, cnt)= data2.trial{k}(end-9:end-1,1);
end
for k=1:ntrl3
    cnt = cnt+1;
    spnt(:, cnt)= data3.trial{k}(end-9:end-1,1);
end
for k=1:ntrl4
    cnt = cnt+1;
    spnt(:, cnt)= data4.trial{k}(end-9:end-1,1);
end
for k=1:ntrl5
    cnt = cnt+1;
    spnt(:, cnt)= data5.trial{k}(end-9:end-1,1);
end
nsmp = double(nsmp);
spnt = double(spnt);
spnt = sum(spnt.*nsmp(ones(9,1),:),2)./sum(nsmp);

M = headcoordinates(spnt(1:3),spnt(4:6),spnt(7:9));
M(1:3,4) = M(1:3,4) * 100; %from m to cm

%---save the results
cd(subject.pathname);
cd('dewar2head_avg');
save([subject.name,'dewar2head_avg'], 'M');