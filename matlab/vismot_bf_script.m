% this script performs spectral analysis of a given subject, using pre-computed
% data, and divides pre and post cue onset intervals per condition
if ~exist('latoi', 'var'), error('latency of interest (latoi) should be specified (pre/post)'); end
if ~exist('frequency', 'var'), error('frequency should be defined'); end
if ~exist('subjectname', 'var'), error('subjectname needs to be defined'); end
if ~exist('smoothing', 'var'), smoothing = []; end
if ~exist('prewhiten', 'var'), prewhiten = false; end
if ~exist('lambda', 'var'), lambda = []; end
if ~exist('nrand', 'var'), nrand = 0; end
if ~exist('stratifyflag', 'var'), stratifyflag = false; end
if isempty(smoothing)
  if frequency < 30, smoothing = 4;
  else, smoothing = 8; end
end
if ~exist('singletrialpow'); singletrialpow = false; end

subject = vismot_subjinfo(subjectname);

load(fullfile(subject.pathname,'grid',sprintf('%s_sourcemodel3d4mm',subject.name)),'sourcemodel');
[source, stat, filter] = vismot_bf(subject,'frequency',frequency,'sourcemodel',sourcemodel,'prewhiten',prewhiten, 'lambda', lambda, 'nrand', nrand, 'latoi', latoi, 'stratifyflag', stratifyflag);

%% post process stats
filename = fullfile(subject.pathname,'source', [subject.name,sprintf('_source3d4mm_%s_', latoi), num2str(frequency,'%03d')]);
if istrue(prewhiten), filename = [filename '_prewhitened']; end
if istrue(stratifyflag), filename = [filename, '_stratified']; end
if nrand>0, filename = [filename, '_resamp']; end

stat13 = stat.stat13;
stat42 = stat.stat42;
stat12 = stat.stat12;
stat43 = stat.stat43;
stat15 = stat.stat15;
stat25 = stat.stat25;
stat35 = stat.stat35;
stat45 = stat.stat45;
statCvsIC = stat.statCvsIC;
statCvsIC2 = stat.statCvsIC2;


% scrub the headmodel and grid from the output cfg
for k = 1:numel(source)
  source(k).cfg = removefields(source(k).cfg, {'grid' 'headmodel'});
  source(k).cfg.callinfo.usercfg = removefields(source(k).cfg.callinfo.usercfg, {'grid' 'headmodel'});
end

if nrand>0
  zx13 = rmfield(stat13, 'stat');
  zx13.stat = stat.zx13;
  zx42 = rmfield(stat42, 'stat');
  zx42.stat = stat.zx42;
  stat_resamp.zx13 = zx13;
  stat_resamp.zx42 = zx42;
  stat_resamp.statResp = rmfield(stat13, 'stat');
  tmpzx42 = hemiflip(stat_resamp.zx42, 'stat');
  stat_resamp.statResp.stat = (stat_resamp.zx13.stat + tmpzx42.stat)./2;
else
  stat_resamp=[];
end

% hemiflip right handed response/ right hemifield
if isfield(stat42, 'statsmooth')
  parameter = {'stat', 'statsmooth'};
else
  parameter = 'stat';
end
tmp42 = hemiflip(stat42, parameter);
tmp43 = hemiflip(stat43, parameter);
tmp25 = hemiflip(stat25, parameter);
tmp45 = hemiflip(stat45, parameter);

% treat as if everything is left handed response
statResp = stat13;
statResp.stat = (statResp.stat + tmp42.stat)/2;

% treat as if cue presented in left hemifield
statHemi = stat12;
statHemi.stat = (statHemi.stat + tmp43.stat)/2;

statCvsN       = stat15;
statCvsN.stat  = (statCvsN.stat + tmp45.stat)/2;
statICvsN      = stat35;
statICvsN.stat = (statICvsN.stat + tmp25.stat)/2;

%save(filename, 'source', 'stat13', 'stat42','stat12', 'stat43','stat15','stat25','stat35','stat45', 'statResp', 'statHemi', 'statCvsN', 'statICvsN', 'statCvsIC', 'smoothing');

stat = rmfield(stat13,'stat');
stat.stat13 = stat13.stat;
stat.stat42 = stat42.stat;
stat.stat12 = stat12.stat;
stat.stat43 = stat43.stat;
stat.stat15 = stat15.stat;
stat.stat25 = stat25.stat;
stat.stat35 = stat35.stat;
stat.stat45 = stat45.stat;
stat.statResp = statResp.stat;
stat.statHemi = statHemi.stat;
stat.statCvsN = statCvsN.stat;
stat.statICvsN = statICvsN.stat;
stat.statCvsIC = statCvsIC.stat;
stat.statCvsIC2 = statCvsIC2.stat;


save(filename, 'stat', 'smoothing', 'lambda', 'stat_resamp', 'nrand');

% compute single trial source power
if singletrialpow
  freq    = vismot_spectral(subject, 'foilim', [frequency frequency], 'toi', latoi, 'balance', true, 'smoothing', smoothing, 'output', 'fourier', 'prewhiten', prewhiten);
  F = zeros(numel(source(1).inside), numel(freq(1).label));
  F(source(1).inside,:) = cat(1, filter{:});
  
  for k=1:numel(freq)
    tmpfreq = freq(k);
    nrpt = numel(tmpfreq.cumtapcnt);
    ntap = tmpfreq.cumtapcnt(1);
    
    ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
    iy = 1:(nrpt*ntap);
    iz = ones(nrpt*ntap,1)./ntap;
    P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);
    
    tmpsource{k} = removefields(source(k), {'avg', 'cumtapcnt'});
    
    tmpsource{k}.pow = transpose((abs(F*transpose(tmpfreq.fourierspctrm)).^2)*P);
  end
  cfg=[];
  cfg.parameter = 'pow';
  sourcepow = ft_appendsource(cfg, tmpsource{:});
  sourcepow.dim = tmpsource{1}.dim;
  sourcepow.inside = tmpsource{1}.inside;
  
  filename = fullfile(subject.pathname,'pow', [subject.name,sprintf('_source3d4mm_%s_', latoi), num2str(frequency,'%03d')]);
  sourcepow = removefields(sourcepow, 'cfg');
  save(filename, 'sourcepow', 'filter');
end

