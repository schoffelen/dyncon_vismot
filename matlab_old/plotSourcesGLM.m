mri = read_mri('/home/jan/matlab/mri/templateMRI.mnc');

cd /analyse/4/Project0030/source/4DprepstGLM/
load grandavg_alpha;
statC.dimord  = 'pos';
statC         = rmfield(statC, 'dim');
statC.prob(statC.outside) = 1;
statB.dimord  = 'pos';
statB         = rmfield(statB, 'dim');
statB.prob(statB.outside) = 1;
cfg.parameter = {'stat' 'prob' 'mask'};
ialpha        = sourceinterpolate(cfg, statC, mri);
ialphaB       = sourceinterpolate(cfg, statB, mri);
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

cd /analyse/4/Project0030/figures/sourcedata/4DprepstGLM/
cfgp = [];
cfgp.method = 'slice';
cfgp.funparameter = 'stat';
cfgp.maskparameter = 'mask';
cfgp.nslices = 9;
cfgp.slicerange = [70 155];
cfgp.opacitylim = [0.5 1.2];
figure;sourceplot(cfgp,ialpha);
print(gcf,'-dpng','statsAlphaGLMCongruency');
figure;sourceplot(cfgp,ialphaB);
print(gcf,'-dpng','statsAlphaGLMCongruencyB');
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
