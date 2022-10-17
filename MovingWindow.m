clear all;
close all;

clc;

% loading the list of directories
Dirs = cell(5, 1);
Dirs{1} = getenv('HOME')+"/Dropbox/GermbandRetraction/myosin_orientationJ/";
Dirs{2} = getenv('HOME')+"/Dropbox/GermbandRetraction/control_ECadherin_movies_for_PIV_and_shapes_analysis/";
Dirs{3} = getenv('HOME')+"/Dropbox/GermbandRetraction/control_nuclei_(byn in_red_all_in green)_labelled/";
Dirs{4} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_movies_ECadh_for_PIV/";
Dirs{5} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_nuclei/";
Dirs{6} = getenv('HOME')+"/Dropbox/GermbandRetraction/cable_ablation/";

% Loading images from that Directory
dir_i = 2; % corresponding to which number of Dirs above
j = 2; % corresponding to the Expt number inside this Dir
k = 2; % which set for the current experiment
% creating the Directory name which should be loaded
if (dir_i == 3)
    SubDir = "Expt-" + int2str(j) + "/Total";
else
    SubDir = "Expt-" + int2str(j);
end

pivfile = char(Dirs{dir_i} + SubDir + "/piv/PIVData.mat");
% The starting and the ending frame for piv
load(pivfile);

ExptData

% space from px to \mu m
dsp = cell(5, 2);
dsp{1}{1} = 0.248;
dsp{1}{2} = 0.207;
dsp{1}{3} = 0.248;
dsp{2}{2} = 0.248;
dsp{2}{5} = 0.262;
dsp{2}{6} = 0.262;

ds = dsp{dir_i}{j};

%%
dtime = cell(5, 2);
dtime{1}{1} = 0.442;
dtime{1}{2} = 0.456;
dtime{1}{3} = 0.165;
dtime{2}{2} = 1.37;
dtime{2}{5} = 0.575;
dtime{2}{6} = 0.575;

dt = dtime{dir_i}{j};
%%
% These two parameters have to changed 

if (dir_i == 2)
    
    FrmBgn = CadCtrl{j}{k}{1}(1);
    FrmEnd = CadCtrl{j}{k}{1}(2);

    w = CadCtrl{j}{k}{2}(1); % 
    h = CadCtrl{j}{k}{2}(2); % 

    ax = CadCtrl{j}{k}{3}(1);
    ay = CadCtrl{j}{k}{3}(2);
end

if (dir_i == 1)
    FrmBgn = MyoCtrl{j}{k}{1}(1);
    FrmEnd = MyoCtrl{j}{k}{1}(2);

    w = MyoCtrl{j}{k}{2}(1); % 
    h = MyoCtrl{j}{k}{2}(2); % 

    ax = MyoCtrl{j}{k}{3}(1);
    ay = MyoCtrl{j}{k}{3}(2);  
end
%%

% load ROI
%coods = [53, 68, 131, 128]; % this has to be loaded automatically
%coods = [123, 35, 124, 211];
%coods = [138, 12, 85, 122]; % This has to keep changing between experiments
%coods = [58, 63, 137, 246];
% parameters for ROI
%ax = coods(1); % top left x
%ay = coods(2); % top left y
%w = coods(3);
%h = coods(4);


%ExptDir = MainDir + "e6/Pt2_44-67/";
%ImDir = MainDir + "e6/Pt2_44-67/Img/";
ROIStrip =  "/ROI/ROI_Full/" + sprintf('%d', k) + "/";
FullROI =   "/ROI/FullROI_Full/" + sprintf('%d', k) + "/";

ROIDir = Dirs{dir_i} + SubDir + ROIStrip;
FullROIDir = Dirs{dir_i} + SubDir + FullROI;

if ~exist(ROIDir, 'dir')
       mkdir(ROIDir);
end


if ~exist(FullROIDir, 'dir')
       mkdir(FullROIDir);
end

% what kind of analysis: rectangle or ellipse. 
% we will predominantly use rectangular domain

typeanalysis = 'rectangle';
%%
switch typeanalysis
    case 'rectangle'

        Nframes = FrmEnd - FrmBgn + 1;
        Ux = zeros(Nframes, 1);
        Uy = zeros(Nframes, 1);
        Uvarx = zeros(Nframes, 1);
        Uvary = zeros(Nframes, 1);
       

        rectx = [0 w w 0 0] + ax;
        recty = [0 0 h h 0] + ay;
%%        
%%

        %imPartStack = zeros(h+1, w+1, Nframes);
        %j = FrmBgn;
        for i = FrmBgn:FrmEnd

           %figure(1)
           
           imfile = char(Dirs{dir_i} + SubDir + "/Original/Original_" + sprintf('%04d', i) + ".tif");
           im = imread(imfile);
    
           % check if the image is rgb and convert it to gray
           [sy, sx, sz] = size(im);
           
           if sz > 1
               im = rgb2gray(im);
           end

           imrangex = floor(ax)+1:floor(ax)+w+1;
           imrangey = floor(ay)+1:floor(ay)+h;

           impart = im(imrangey, imrangex);
           size(impart)
           %imPartStack(:, :, i) = impart;
           %imshow(impart);
           impartfile = char(ROIDir +  "ROI_" + sprintf('%04d', i) + ".tif");
           imwindowfile = char(FullROIDir  + "FullROI_"+ sprintf('%04d', i) + ".tif");
           %clear impart

           %uxi = ux(i).w;
           %uyi = uy(i).w;

           uxi = u_filt{i};
           uyi = v_filt{i};
           X = x{i}(1,:);
           Y = y{i}(:,1);

           xr = (X >= ax & X <= ax + w);
           yr = (Y >= ay & Y <= ay + h);

           vxrect = uxi(yr, xr);
           vyrect = uyi(yr, xr);

           [xrmesh, yrmesh] = meshgrid(X(xr), Y(yr));


           fh = figure(1)
           imshow(im);
           hold on;
           plot(rectx, recty, 'r-');
           rectangle('Position', [ax, ay, w, h], 'EdgeColor', 'r');
           axis equal;
           hold on
           quiver(xrmesh, yrmesh, 20*vxrect, 20*vyrect, 'color', 'red')  
           set(gcf,'Position',[0 0 511 225])
           %hold off
           %pause(0.5)
           frm = getframe( fh )
      
           imwrite(frm.cdata, imwindowfile);
           export_fig(imwindowfile, '-tif', '-m4');
           close


           figure(2)
           imshow(impart);
           %pause(0.1)
           close
           imwrite(impart, impartfile);
           %hold off;
           close;

           Ux(i) = mean(vxrect(:));
           Uy(i) = mean(vyrect(:));
           
           Uvarx(i) = std(vxrect(:));
           Uvary(i) = std(vyrect(:));

           %ax = ax + 10.74;
           ax = ax + Ux(i);
           ay = ay + Uy(i);
           %ax = ax + 3.1;
           rectx = rectx + Ux(i);
           recty = recty + Uy(i);

           if (ax + w > sx || ay + h > sy)
               i
               break;
           end
           %j = j + 1;
       %figure(3)

        end
%%
        figure(2)
        a = Ux ~= 0;
        %plot(Ux(a));
        %t = (FrmBgn-1 + linspace(1, sum(a)))*dt;
        t = ((1:sum(a)) + FrmBgn-1)*dt;
        %t = (FrmBgn:FrmEnd)*dt;
        errorbar(t, ds/dt*Ux(a), ds/dt*Uvarx(a));
        

        hold on;
        %plot(Uy(a));
        errorbar(t, ds/dt*Uy(a), ds/dt*Uvary(a));
        
        legend('<vx>','<vy>')
        xlabel('time (min)');
        ylabel('velocity (microns/min)');
        
        FileName = FullROIDir + "Velocity.tif";
        export_fig(FileName, '-tif', '-m4');
        
        hold on;
        xmean = mean(ds/dt*Ux(a));
        ymean = mean(ds/dt*Uy(a));
        
        plot(t, xmean*ones(length(t), 1), t, ymean*ones(length(t), 1))
    
    % ellipse case needs to be remodified
    case 'ellipse'

    % With ellipse fitting

    % Used for Maithreyi's Zoomed Cadherin Data

        Nframes = 30;
        start = 1;
        Ux = zeros(Nframes, 1);
        Uy = zeros(Nframes, 1);
        %ux = spaverf(vec2scal(vfil, 'vx'), 'xy');
        %uy = spaverf(vec2scal(vfil, 'vy'), 'xy');
        ExptDir = MainDir + "e6/Pt1_1-42/";
        ImDir = MainDir + "e6/Pt1_1-42/Img/";
        %creating rectangle moving it and plotting it
        % Outside
        % ax = 91;
        % ay = 36;
        % w = 102;
        % h = 44;

        %ax = 82;
        % ay = 44;
        % w = 98;
        % h = 36;

        % Middle
         ax = 80
         ay = 63
         w = 159
         h = 63

        % Inside
        % ax = 75;
        % ay = 111;
        % 
        % w = 80;
        % h = 70;

        % vertical 
        %ax = 162;
        %ay = 83;
        %w = 34;
        %h = 103;
        %w = 102;
        %h = 44;
        %h = 102

        %rectx = [0 w w 0 0] + ax;
        %recty = [0 0 h h 0]+ay;


        for i = start:Nframes

           %figure(1)
           imfile = ImDir+ "Images" + sprintf('%04d', i) + ".tif";
           im = imread(imfile);


           imrangex = floor(ax)-w/2:floor(ax)+w/2+1;
           imrangey = floor(ay)-h/2:floor(ay)+h/2+1;


           impart = im(imrangey, imrangex);
           size(impart)
           %imPartStack(:, :, i) = impart;
           %imshow(impart);
           impartfile = ExptDir + "Window_Ellipse/" + "Window_" + sprintf('%04d', i) + ".tif";

           %clear impart

           %uxi = ux(i).w;
           %uyi = uy(i).w;


           uxi = v(i).vx;
           uyi = v(i).vy;
           x = v(i).x;
           y = v(i).y;

           [xmesh, ymesh] = meshgrid(x, y);


            r = (xmesh-ax).^2/(w/2)^2 + (ymesh-ay).^2/(h/2)^2 < 1
            xr = xmesh(r);
            yr = ymesh(r);

            vxellipse = uxi(r);
            vyellipse = uyi(r);


           fh = figure(1)
           imshow(im);
           hold on;
           %plot(rectx, recty, 'r-');
           %rectangle('Position', [ax, ay, w, h], 'EdgeColor', 'r');
           drawellipse('Center',[ax,ay],'SemiAxes',[w/2,h/2],'StripeColor','r');


           axis equal;
           hold on
           quiver(xr, yr, 20*vxellipse, 20*vyellipse, 'color', 'red')  
           set(gcf,'Position',[0 0 511 225])
           %hold off
           pause(0.);
           frm = getframe( fh )
           imwindowfile = ExptDir + "ImWindow_Ellipse/" + "ImBox_" + sprintf('%04d', i) + ".tif";
           imwrite(frm.cdata, imwindowfile);
           %export_fig(imwindowfile, '-tif', '-m4');
           close

             s = size(impart);
            sx = s(2);
            sy = s(1);

            [xm, ym] = meshgrid(1:sx, 1:sy)
            impart = double(impart);
            cx = sx/2;
            cy = sy/2;

            r = (xm - cx).^2/(w/2)^2 + (ym - cy).^2/(h/2)^2 <= 1;

            newellipse = r.*impart;
            newellipse = uint8(newellipse)


           figure(2)
           %imshow(impart);
           imshow(newellipse);
           pause(0);
           close
           imwrite(newellipse, impartfile);
           %hold off;
           close;

           Ux(i) = mean(vxellipse(:));
           Uy(i) = mean(vyellipse(:));

           %ax = ax + Ux(i);
           ax = ax + 4.5;
           %ay = ay + Uy(i);
           %ax = ax + 3.1;
           %rectx = rectx + Ux(i);

           if (ax + w/2 > max(x))
               i
               break;
           end

           %figure(3)

        end

        figure(2)
        plot(Ux)
        hold on
        plot(Uy)
        plot(Uy)
        
end



