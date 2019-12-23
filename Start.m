load distance.mat
load site.mat %区县名字
min=9999;
tic;
for j=1:1:10000
    dmat=distance;
    salesmen=2;
    min_tour=6500;
    pop_size=8;
    num_iter=20000;
    [opt_rte,opt_brk,min_dist]= mtspofs_ga(dmat,salesmen,min_tour,pop_size,num_iter,1,1);
    % 输出路径

    str1 = '甲：碑林区';
    str3=[min_dist,"单位:km",j];
    opt_brk_position=find(opt_rte==opt_brk);

    for i=1:1:opt_brk_position-1
        site_print = site(opt_rte(i));
        str1 = [str1,'->',site_print{1}];
    end

    str2 = '乙：碑林区';
    for i=opt_brk_position:1:106
        site_print = site(opt_rte(i));
        str2 = [str2,'->',site_print{1}];
    end

    disp(str1)
    disp(str2)
    disp(str3)
    
    if min>min_dist
        min=min_dist;
        min_index=j;
        str1_min=str1;
        str2_min=str2;
    end
end
str3=[min,min_index];
disp(str1_min)
disp(str2_min)
disp(str3)
toc