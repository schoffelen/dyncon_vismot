[x,y,alpha1,alpha2,beta1,beta2] = anovacode;

cnt = 1;
figure;volplotJM(alpha1.m0, 'montage');
title('omnibus alpha roi left');
print(gcf, '-dpsc2', 'Ar1M0');
fname{cnt} = 'Ar1M0.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.m0, 'montage');
title('omnibus alpha roi right');
print(gcf, '-dpsc2', 'Ar2M0');
fname{cnt} = 'Ar2M0.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.m1, 'montage');
title('main 1 alpha roi left');
print(gcf, '-dpsc2', 'Ar1M1');
fname{cnt} = 'Ar1M1.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.m1, 'montage');
title('main 1 alpha roi right');
print(gcf, '-dpsc2', 'Ar2M1');
fname{cnt} = 'Ar2M1.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.m2, 'montage');
title('main 2 alpha roi left');
print(gcf, '-dpsc2', 'Ar1M2');
fname{cnt} = 'Ar1M2.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.m2, 'montage');
title('main 2 alpha roi right');
print(gcf, '-dpsc2', 'Ar2M2');
fname{cnt} = 'Ar2M2.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.m3, 'montage');
title('main 2 alpha roi left');
print(gcf, '-dpsc2', 'Ar1M3');
fname{cnt} = 'Ar1M3.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.m3, 'montage');
title('main 2 alpha roi right');
print(gcf, '-dpsc2', 'Ar2M3');
fname{cnt} = 'Ar2M3.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.i12, 'montage');
title('interaction 12 alpha roi left');
print(gcf, '-dpsc2', 'Ar1I12');
fname{cnt} = 'Ar1I12.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.i12, 'montage');
title('interaction12 alpha roi right');
print(gcf, '-dpsc2', 'Ar2I12');
fname{cnt} = 'Ar2I12.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.i13, 'montage');
title('interaction 13 alpha roi left');
print(gcf, '-dpsc2', 'Ar1I13');
fname{cnt} = 'Ar1I13.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.i13, 'montage');
title('interaction13 alpha roi right');
print(gcf, '-dpsc2', 'Ar2I13');
fname{cnt} = 'Ar2I13.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.i23, 'montage');
title('interaction 23 alpha roi left');
print(gcf, '-dpsc2', 'Ar1I23');
fname{cnt} = 'Ar1I23.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.i23, 'montage');
title('interaction23 alpha roi right');
print(gcf, '-dpsc2', 'Ar2I23');
fname{cnt} = 'Ar2I23.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha1.i123, 'montage');
title('interaction 123 alpha roi left');
print(gcf, '-dpsc2', 'Ar1I123');
fname{cnt} = 'Ar1I123.ps';
close; cnt = cnt+1;
figure;volplotJM(alpha2.i123, 'montage');
title('interaction123 alpha roi right');
print(gcf, '-dpsc2', 'Ar2I123');
fname{cnt} = 'Ar2I123.ps';
close; cnt = cnt+1;

%beta
figure;volplotJM(beta1.m0, 'montage');
title('omnibus beta roi left');
print(gcf, '-dpsc2', 'Br1M0');
fname{cnt} = 'Br1M0.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.m0, 'montage');
title('omnibus beta roi right');
print(gcf, '-dpsc2', 'Br2M0');
fname{cnt} = 'Br2M0.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.m1, 'montage');
title('main 1 beta roi left');
print(gcf, '-dpsc2', 'Br1M1');
fname{cnt} = 'Br1M1.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.m1, 'montage');
title('main 1 beta roi right');
print(gcf, '-dpsc2', 'Br2M1');
fname{cnt} = 'Br2M1.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.m2, 'montage');
title('main 2 beta roi left');
print(gcf, '-dpsc2', 'Br1M2');
fname{cnt} = 'Br1M2.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.m2, 'montage');
title('main 2 beta roi right');
print(gcf, '-dpsc2', 'Br2M2');
fname{cnt} = 'Br2M2.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.m3, 'montage');
title('main 2 beta roi left');
print(gcf, '-dpsc2', 'Br1M3');
fname{cnt} = 'Br1M3.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.m3, 'montage');
title('main 2 beta roi right');
print(gcf, '-dpsc2', 'Br2M3');
fname{cnt} = 'Br2M3.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.i12, 'montage');
title('interaction 12 beta roi left');
print(gcf, '-dpsc2', 'Br1I12');
fname{cnt} = 'Br1I12.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.i12, 'montage');
title('interaction12 beta roi right');
print(gcf, '-dpsc2', 'Br2I12');
fname{cnt} = 'Br2I12.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.i13, 'montage');
title('interaction 13 beta roi left');
print(gcf, '-dpsc2', 'Br1I13');
fname{cnt} = 'Br1I13.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.i13, 'montage');
title('interaction13 beta roi right');
print(gcf, '-dpsc2', 'Br2I13');
fname{cnt} = 'Br2I13.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.i23, 'montage');
title('interaction 23 beta roi left');
print(gcf, '-dpsc2', 'Br1I23');
fname{cnt} = 'Br1I23.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.i23, 'montage');
title('interaction23 beta roi right');
print(gcf, '-dpsc2', 'Br2I23');
fname{cnt} = 'Br2I23.ps';
close; cnt = cnt+1;
figure;volplotJM(beta1.i123, 'montage');
title('interaction 123 beta roi left');
print(gcf, '-dpsc2', 'Br1I123');
fname{cnt} = 'Br1I123.ps';
close; cnt = cnt+1;
figure;volplotJM(beta2.i123, 'montage');
title('interaction123 beta roi right');
print(gcf, '-dpsc2', 'Br2I123');
fname{cnt} = 'Br2I123.ps';
close; cnt = cnt+1;

