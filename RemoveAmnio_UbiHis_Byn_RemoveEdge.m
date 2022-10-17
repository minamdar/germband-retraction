clear all;
close all;
clc;
%% Reading all the nematic data

A = readtable('./Myosin/Nematic.csv');
%%
expt = 3;

if expt == 2
    BasicFileName = "./UbiHis/Images_Byn/Byn_";
    BasicTotalFileName = "./UbiHis/Images_All/Total_";
    saveFile = "./UbiHis/PIVImages/";
    pivbase = "./UbiHis/pivlab_all_2022_03_03/PIVlab_";
    ds = 0.573;
    dt = 0.98;
    
    Nstart = 1;
    %NImages = 80;
    NImages = 50;
    mean_int = [];
    pts = 1:5:329;
    yy = 1:200;
    int_th = 0.1;
else
    BasicFileName = "./UbiHis_1/Images_Byn/Byn_";
    BasicTotalFileName = "./UbiHis_1/Images_All/Total_";
    saveFile = "./UbiHis_1/PIVImages/";
    pivbase = "./UbiHis_1/pivlab_all_2022_03_03/PIVlab_";
    ds = 0.496;
    dt = 0.94;
    
    Nstart = 1;
    NImages = 144;
    mean_int = [];
    pts = 1:5:509;
    %pts = 509:-5:1
    yy = 1:162;
    int_th = 0.2;
end

%%

%pts = 1:5:512;


% data for ubihis_1
%{
Nstart = 1;
NImages = 144;
mean_int = [];
pts = 1:5:509;
%pts = 509:-5:1
yy = 1:162;
int_th = 0.2;
%}

recbox = [];
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
for i = Nstart:Nstart+NImages-1  
    file = char(BasicFileName + sprintf('%04d', i) + ".tif");
    img = imread(file);    
    imshow(img);
    %img = double(img);
    
    k = 1;
    for j = pts
       imsub = img(:, 1:j);
       mean_int(k) = mean(imsub(:));
       k = k + 1;
    end
    
    mean_int = mean_int - mean_int(1);
    mean_int = mean_int/max(mean_int);
    max_pts = max(pts(mean_int < int_th));
    
    if (isempty(max_pts))
        max_pts = recbox(i-1);
        recbox(i) = max_pts;
    else
        recbox(i) = max_pts;
    end
    
    xx = max_pts*ones(length(yy), 1);
    hold on
    plot(xx, yy, 'r', 'LineWidth', 1);
    %plot(mean_int)
    pause(0.2);
    
    hold off;
    %imname = char("./UbiHis/ImageBnd/Img_Bnd_" + sprintf('%04d', i) + ".tif");
    %imwrite(
end

close all;

%%
%Nstart = 1;
%NImages = 36;
Veldata = cell(NImages, 1);

vtime = [];
vytime = [];
vartime = [];
varytime = [];
for i = Nstart:Nstart+NImages-1 
    file = char(pivbase + sprintf('%04d', i) + ".txt");
    data = dlmread(file);
    
    x = data(:,1);
    y = data(:,2);
    vx = data(:,3);
    vy = data(:, 4);
  
   if (i >= 126)
       leftbdr = 410;
        X = x((x > recbox(i)) & (x < leftbdr)); % beyond 200 px we have blanck
        Y = y(x > recbox(i) & x < leftbdr);
        Vx = vx((x > recbox(i)) & (x < leftbdr));  
        Vy = vy(x > recbox(i) & x < leftbdr);
   else    
        X = x((x > recbox(i)) );
        Y = y(x > recbox(i) );
        Vx = vx((x > recbox(i)) );  
        Vy = vy(x > recbox(i) );
   end
    
    [min(X)]
    
    [Xg, Yg] = meshgrid(unique(X), unique(Y));
    Vxg = griddata(X, Y, Vx, Xg, Yg);
    Vyg = griddata(X, Y, Vy, Xg, Yg);
    
    % converting velocity from px/frm to microns/min
    Vxg = Vxg*ds/dt;
    Vyg = Vyg*ds/dt;
    
    Vxm = mean(Vxg, 1);
    Vym = mean(Vyg, 1);
    Veldata{i}{1} = Vxm;
    Veldata{i}{2} = unique(X');
    figure(2)
    %plot(Veldata{i}{2}, Veldata{i}{1});
    xlabel('x position');
    ylabel('vx (px/frame)');
    title(int2str(i));
    xlim([1, 500]);
    ylim([0, 13]);
    figname = char("./UbiHis_1/VelPlot_3/Plot_" + sprintf('%d', i) + ".png");
    %export_fig(figname);
    %close 
    %pause(0.1);
    vtime = [vtime; mean(Vxm)];
    vytime = [vytime; mean(Vym)];
    %vtime = [vtime; mean(Vx(:))];
    vartime = [vartime; std(Vx(:))];
    varytime = [vartime; std(Vy(:))];
    %figure(3)
    %quiver(Xg, Yg, Vxg, Vyg);
    %axis equal
    %pause(0.1);
    close;
    
    
    % plot the images and the corresponding piv
    file = char(BasicTotalFileName + sprintf('%04d', i) + ".tif");
    img = imread(file);  
    figure(2)
    imshow(img);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
    hold on
    quiver(X(:), Y(:), Vx(:), Vy(:), 'g', 'LineWidth', 1, 'Autoscale', 'off');
    axis equal
    %pause(0.1);
    hold off
    file = char(saveFile + sprintf('%04d', i) + ".png");
    export_fig(file, '-m2.5');
    
    close;
end



%%
hold on
time = dt*(Nstart:Nstart+NImages-1);
a = (vtime ~= 0);
%plot( Nstart:Nstart+NImages-1 , vtime, 'b*')
errorbar(time(a), vtime(a), vartime(a), 'r*');

%%
close all;
openfig('./Sliding_Speed/Figures_2022_03_03/CombinedPlot_NewPIV.fig');

hold on
if (expt == 2)
    lag = 81;
    a = (vtime ~= 0);
%time = time';
%time = 1:sum(a);
%time = time';
    dy = vartime;  % made-up error values
    %fill([time(a)'+ lag;flipud(time(a)'+ lag)],[vtime(a)-dy(a);flipud(vtime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
    line(time(a)+lag, vtime(a));
    hold on
    scatter(time(a)+lag, vtime(a), 'b*')
else
     a = (vtime ~= 0);
%time = time';
%time = 1:sum(a);
%time = time';
    dy = vartime;  % made-up error values
    %fill([time(a)';flipud(time(a)')],[vtime(a)-dy(a);flipud(vtime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
    %line(time(a), vtime(a));
    hold on
    scatter(time(a), vtime(a), '*')
end

%fill([time(a);flipud(time(a))],[vtime(a)-dy(a);flipud(vtime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
%line(time(a), vtime(a));

xlabel('time (min)');
ylabel('velocity (microns/min)');
%xlim([0, 135]);
%%


close all;
%openfig('./Sliding_Speed/Figures_2022_03_03/CombinedPlot_NewPIV.fig');

hold on
if (expt == 2)
    lag = 81;
    a = (vtime ~= 0);
%time = time';
%time = 1:sum(a);
%time = time';
    dy = varytime;  % made-up error values
    fill([time(a)'+ lag;flipud(time(a)'+ lag)],[vytime(a)-dy(a);flipud(vytime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
    line(time(a)+lag, vtime(a));
    hold on
    scatter(time(a)+lag, vtime(a), 'b*')
else
     a = (vtime ~= 0);
%time = time';
%time = 1:sum(a);
%time = time';
    dy = varytime;  % made-up error values
    fill([time(a)';flipud(time(a)')],[vytime(a)-dy(a);flipud(vytime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
    %line(time(a), vtime(a));
    hold on
    scatter(time(a), vytime(a), '*')
end

%fill([time(a);flipud(time(a))],[vtime(a)-dy(a);flipud(vtime(a)+dy(a))],[.9 0.9 0.9 ],'linestyle','none');
%line(time(a), vtime(a));

xlabel('time (min)');
ylabel('velocity (microns/min)');