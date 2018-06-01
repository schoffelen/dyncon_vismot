%subjlist = [1 2 3 5 6 8 9 10 11 12 13 14 15 16 18 19 20];
subjlist = [1 2 3 5 6 8 9 10 11 13 14 15 16 18 19 20];

for k = 1:numel(subjlist)
  subjname = SUBJ(subjlist(k)).name;
  load([subjname,'dcohFastSlow']);
  if k==1,
    ss = s;
  else
    for kk = 1:4
      ss{kk}.dcoh = s{kk}.dcoh + ss{kk}.dcoh; 
      ss{kk}.coh1 = s{kk}.coh1 + ss{kk}.coh1;
      ss{kk}.coh2 = s{kk}.coh2 + ss{kk}.coh2; 
    end
  end
end
