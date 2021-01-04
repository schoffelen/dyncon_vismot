% replace the original names (e.g. AHE08 with subject numbers, e.g. sub01)
%% Define new subject names
list_orig =     [
  {'AHE08'}
  {'AMN20'}
  {'BKA01'}
  {'CBE22'}
  {'CDE04'}
  {'DBD12'}
  {'GAR12'}
  {'JGN27'}
  {'JOE22'}
  {'LCE09'}
  {'MAI27'}
  {'MHH14'}
  {'MME25'}
  {'PCL19'}
  {'RBE13'}
  {'SZQ02'}
  {'T001' }
  {'TMR04'}
  {'VIA12'}];
list_anon = cell(numel(list_orig),1);
for k=1:numel(list_orig)
  list_anon{k} = sprintf('sub%02d', k);
end
list_transform = [list_orig, list_anon];
save('list_transform.mat', 'list_transform')

list = list_anon;
save('list.mat', 'list')

%% Change filenames of folders
go = false; % change to true to change filenames

 % change all filenames in this folder
filelist = dir('/project/3011085.03/analysis/**');

% first change directory names
filelist = filelist([filelist.isdir]);  %only keep folders

% remove '.' and '..' directories
rem = [];
for k=1:numel(filelist)
  if strcmp(filelist(k).name, '.') || strcmp(filelist(k).name, '..')
    rem = [rem, k];
  end
end
filelist(rem) = [];

% loop through folders and replace subjectnames
status = zeros(numel(filelist),1);
for k=1:numel(filelist)
  fname_orig = fullfile(filelist(k).folder, filelist(k).name);
  if contains(fname_orig, list_orig)
    for j=1:numel(list)
      if contains(fname_orig, list_orig{j})
        fname_new = strrep(fname_orig, list_orig{j}, list_anon{j});
        if go, status(k) = movefile(fname_orig, fname_new); end
        break
      end
    end
  end
end
if any(status) % possibly because a directory contains double subjectnames, 
  % and these were changed in a previous iteration
  error('SOMETHING WENT HORRIBLY WRONG! CHECK WHICH FOLDERS WERE NOT CORRECTLY MOVED!')
end


%% Now continue with files
filelist = dir('/project/3011085.03/analysis/**');
filelist = filelist(~[filelist.isdir]);  %only keep files

% loop through folders and replace subjectnames
status = zeros(numel(filelist),1);
for k=1:numel(filelist)
  fname_orig = fullfile(filelist(k).folder, filelist(k).name);
  if contains(fname_orig, list_orig)
    for j=1:numel(list)
      if contains(fname_orig, list_orig{j})
        fname_new = strrep(fname_orig, list_orig{j}, list_anon{j});
        if go, status(k) = movefile(fname_orig, fname_new); end
        break
      end
    end
  end
end
if any(status)
  error('SOMETHING WENT HORRIBLY WRONG! CHECK WHICH FILES WERE NOT CORRECTLY MOVED!')
end





