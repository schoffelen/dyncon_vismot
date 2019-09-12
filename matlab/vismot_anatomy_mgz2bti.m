function vismot_anatomy_mgz2bti(subject)

global ft_default;
ft_default.checksize = inf;

if nargin<1
  subject = vismot_subjinfo;
elseif ischar(subject)
  subject = vismot_subjinfo(subject);
end

for k = [1:numel(subject)]
  if exist(fullfile(subject(k).pathname,'mri',[subject(k).name,'_transform_vox2bti.mat']),'file')
%     continue;
  end
  
  fname = fullfile(subject(k).pathname,'mri',subject(k).name);
  mri   = ft_read_mri([fname,'.mgz']);
  
  cfg = [];
  cfg.coordsys = 'bti';
  mri = ft_volumerealign(cfg, mri);
  transform_interactive = mri.transform;
  
  hsfile = fullfile(subject(k).rawpath,subject(k).name,subject(k).scanname,subject(k).sessionname,subject(k).runnames{1},'hs_file');
          
  cfg           = [];
  cfg.method    = 'headshape';
  cfg.headshape.headshape = ft_read_headshape(hsfile);
  cfg.headshape.icp       = 1;
  cfg.headshape.interactive = 0;
  cfg.coordsys = 'bti';
  mri           = ft_volumerealign(cfg, mri);
  
  hs     = mri.cfg.headshape.headshape;
  hs_mri = mri.cfg.headshape.headshapemri;
  
  transform = mri.transform;
  mriname   = [fname,'.mgz'];
  save([fname,'_transform_vox2bti'],'transform','mriname','transform_interactive','hs','hs_mri');
  
  clear mri;
end
%         
%           % quality check for the coregistration between headshape and mri
%           headshape    = hcp_ensure_units(headshape,    'mm');
%           headshapemri = hcp_ensure_units(headshapemri, 'mm');
%           
%           figure;
%           subplot('position',[0.01 0.51 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
%           ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(180,-90);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.51 0.51 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
%           ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(0,90);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.01 0.01 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
%           ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance); view(90,0);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.51 0.01 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.7 0.7 0.7],'fidcolor','y','facealpha',0.3);
%           ft_plot_mesh(headshapemri,'edgecolor','none','vertexcolor',headshapemri.distance);  view(0,0); colorbar('east');
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           axis on;
%           grid on;
%           set(gcf,'color','w')
%           hcp_write_figure([outputprefix,'_headshape_distance.png'], gcf, 'resolution', 300); close all;
%           
%           v = headshapemri.pnt;
%           f = headshapemri.tri;
%           [f,v]=reducepatch(f,v, 0.2);
%           headshapemri.pnt = v;
%           headshapemri.tri = f;
%           
%           figure;
%           subplot('position',[0.01 0.51 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
%           ft_plot_headshape(headshape,'vertexsize',3); view(180,-90);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.51 0.51 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
%           ft_plot_headshape(headshape,'vertexsize',3); view(0,90);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.01 0.01 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
%           ft_plot_headshape(headshape,'vertexsize',3); view(90,0);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           subplot('position',[0.51 0.01 0.48 0.48]);hold on;
%           ft_plot_mesh(headshapemri,'edgecolor','none','facecolor',[0.5 0.6 0.8],'fidcolor','y','facealpha',0.3);
%           ft_plot_headshape(headshape,'vertexsize',3); view(0,0);
%           plot3([-130 130],[0 0],[0 0],'k');plot3([0 0],[-120 120],[0 0],'k');plot3([0 0],[0 0],[-100 150],'k');
%           axis on;
%           grid on;
%           set(gcf,'color','w')
%           hcp_write_figure([outputprefix,'_headshape.png'], gcf, 'resolution', 300); close all;
          