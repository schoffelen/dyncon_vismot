function vismot_anatomy_coreg_qc(subject)

if nargin<1
  % do it for all subjects
  vismot_subjinfo;
  for k = [1:3 5:numel(SUBJ)]
    vismot_anatomy_coreg_qc(SUBJ(k));
  end
  return;
end

if ischar(subject)
  vismot_subjinfo;
  sel = find(strcmp({SUBJ(:).name}',subject));
  subject = SUBJ(sel);
end

if ~exist(fullfile(subject.pathname,'mri',[subject.name,'_transform_vox2bti.mat']),'file')
  return;
end
load(fullfile(subject.pathname,'mri',[subject.name,'_transform_vox2bti.mat']));


% quality check for the coregistration between headshape and mri
hs     = ft_convert_units(hs,    'mm');
hs_mri = ft_convert_units(hs_mri, 'mm');

 figure;
 subplot('position',[0.01 0.51 0.48 0.48]);hold on;
 ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
 ft_plot_mesh(hs_mri,'edgecolor','none','vertexcolor',hs_mri.distance); view(180,-90);
 plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
 subplot('position',[0.51 0.51 0.48 0.48]);hold on;
 ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
 ft_plot_mesh(hs_mri,'edgecolor','none','vertexcolor',hs_mri.distance); view(0,90);
 plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
 subplot('position',[0.01 0.01 0.48 0.48]);hold on;
 ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
 ft_plot_mesh(hs_mri,'edgecolor','none','vertexcolor',hs_mri.distance); view(90,0);
 plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
 subplot('position',[0.51 0.01 0.48 0.48]);hold on;
 ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
 ft_plot_mesh(hs_mri,'edgecolor','none','vertexcolor',hs_mri.distance);  view(0,0); colorbar('east');
 plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
 axis on;
 grid on;
 set(gcf,'color','w')


v = hs_mri.pnt;
f = hs_mri.tri;
[f,v]=reducepatch(f,v, 0.2);
hs_mri.pnt = v;
hs_mri.tri = f;
           
figure;
subplot('position',[0.01 0.51 0.48 0.48]);hold on;
ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
ft_plot_headshape(hs,'vertexsize',3); view(180,-90);
plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
subplot('position',[0.51 0.51 0.48 0.48]);hold on;
ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
ft_plot_headshape(hs,'vertexsize',3); view(0,90);
plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
subplot('position',[0.01 0.01 0.48 0.48]);hold on;
ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
ft_plot_headshape(hs,'vertexsize',3); view(90,0);
plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
subplot('position',[0.51 0.01 0.48 0.48]);hold on;
ft_plot_mesh(hs_mri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
ft_plot_headshape(hs,'vertexsize',3); view(0,0);
plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
axis on;
grid on;
set(gcf,'color','w')
