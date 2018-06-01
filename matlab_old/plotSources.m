mri = read_mri('/home/jan/matlab/mri/templateMRI.mnc');

cd /analyse/4/Project0030/source/4Dprepst/
load grandavg_alpha;
statC.dimord  = 'pos';
statC         = rmfield(statC, 'dim');
statC.prob(statC.outside) = 1;
st.dimord  = 'pos';
st         = rmfield(st, 'dim');
cfg.parameter = {'stat' 'prob' 'mask'};
cfg2.parameter = 'stat';
ialpha        = sourceinterpolate(cfg, statC, mri);
ialphaB       = sourceinterpolate(cfg2, st, mri);
load grandavg_beta;
statC.dimord  = 'pos';
statC         = rmfield(statC, 'dim');
statC.prob(statC.outside) = 1;
st.dimord  = 'pos';
st         = rmfield(st, 'dim');
ibeta         = sourceinterpolate(cfg, statC,mri);
ibetaB       = sourceinterpolate(cfg2, st, mri);
load grandavg_gamma1;
statC.dimord  = 'pos';
statC         = rmfield(statC, 'dim');
statC.prob(statC.outside) = 1;
st.dimord  = 'pos';
st         = rmfield(st, 'dim');
igamma1       = sourceinterpolate(cfg, statC,mri);
igamma1B       = sourceinterpolate(cfg2, st, mri);
load grandavg_gamma2;
statC.dimord  = 'pos';
statC         = rmfield(statC, 'dim');
statC.prob(statC.outside) = 1;
st.dimord  = 'pos';
st         = rmfield(st, 'dim');
igamma2       = sourceinterpolate(cfg, statC,mri);
igamma2B       = sourceinterpolate(cfg2, st, mri);

cd /analyse/4/Project0030/figures/sourcedata/4Dprepst/
cfgp = [];
cfgp.method = 'slice';
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'mask';
cfgp.nslices = 9;
cfgp.slicerange = [70 155];
cfgp.opacitylim = [0.5 1.2];
figure;sourceplot(cfgp,ialpha);
print(gcf,'-dpng','statsAlphaCongruency');
figure;sourceplot(cfgp,ibeta);
print(gcf,'-dpng','statsBetaCongruency');
figure;sourceplot(cfgp,igamma2);
print(gcf,'-dpng','statsGamma2Congruency');
cfgp.slicerange = [45 130];
figure;sourceplot(cfgp,igamma1);
print(gcf,'-dpng','statsGamma1Congruency');


cfgp.nslices = 16;
cfgp.slicerange = [45 155];
cfgp.maskparameter = 'stat2';
ialphaB.stat2 = abs(ialphaB.stat);
cfgp.opacitylim = [5 30];
figure;sourceplot(cfgp,ialphaB);
print(gcf,'-dpng','AlphaPstPre');
ibetaB.stat2 = abs(ibetaB.stat);
ibetaB.stat2(~ibetaB.inside) = 0;
figure;sourceplot(cfgp,ibetaB);
print(gcf,'-dpng','BetaPstPre');
igamma1B.stat2 = abs(igamma1B.stat);
cfgp.opacitylim = [3 15];
figure;sourceplot(cfgp,igamma1B);
print(gcf,'-dpng','Gamma1PstPre');
igamma2B.stat2 = abs(igamma2B.stat);
cfgp.opacitylim = [10 20]
figure;sourceplot(cfgp,igamma2B);
print(gcf,'-dpng','Gamma2PstPre');

cfgp2 = [];
cfgp2.method = 'ortho';
cfgp2.interactive = 'yes';
cfgp2.funparameter = 'stat';
cfgp2.maskparameter = 'mask';
cfgp2.opacitylim = [0.5 1.2];
cfgp2.atlas      = '/home/jan/matlab/fieldtrip/TTatlas+tlrc.BRIK';
cfgp2.inputcoord = 'mni';
figure;sourceplot(cfgp2,ialpha);

mri = read_mri('/home/jan/matlab/mri/templateMRI.mnc');
cd /analyse/1/Project0002/tmpProject0030/source;
cfg = [];
cfg.parameter = {'stat' 'prob' 'mask'};
load grandavg_012
i1 = ft_sourceinterpolateJM(cfg,stat,mri);
load grandavg_020
i2 = ft_sourceinterpolateJM(cfg,stat,mri);
load grandavg_048
i3 = ft_sourceinterpolateJM(cfg,stat,mri);
load grandavg_080
i4 = ft_sourceinterpolateJM(cfg,stat,mri);

cfgp = [];
cfgp.funparameter  = 'stat';
cfgp.maskparameter = 'prob';
cfgp.opacitymap    = 'rampdown';
cfgp.opacitylim    = [0 0.01];
%cfgp.method        = 'ortho';
%cfgp.interactive   = 'yes';
cfgp.method = 'slice';
cfgp.nslices = 16;
cfgp.slicerange = [30 160];
cfgp.colorbar = 'no';
i1.prob(~i1.inside)=1;
i2.prob(~i2.inside)=1;
i3.prob(~i3.inside)=1;
i4.prob(~i4.inside)=1;
figure;ft_sourceplotJM(cfgp, i1);
title('colorrange-25');
print(gcf,'-dpng','prepst012');
cfgp.opacitylim    = [0 0.01];
figure;ft_sourceplotJM(cfgp, i2);
title('colorrange-20');
print(gcf,'-dpng','prepst020');
cfgp.opacitylim    = [0 0.07];
figure;ft_sourceplotJM(cfgp, i3);
title('colorrange-13');
print(gcf,'-dpng','prepst048');
cfgp.opacitylim    = [0 0.0025];
figure;ft_sourceplotJM(cfgp, i4);
title('colorrange-18');
print(gcf,'-dpng','prepst080');
