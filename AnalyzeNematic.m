clear all;
close all;


% loading the list of directories
Dirs = cell(5, 1);
Dirs{1} = getenv('HOME')+"/Dropbox/Research/Sudeepa/myosin_orientationJ/";
Dirs{2} = getenv('HOME')+"/Dropbox/Research/Sudeepa/control_ECadherin_movies_for_PIV_and_shapes_analysis/";
Dirs{3} = getenv('HOME')+"/Dropbox/GermbandRetraction/control_nuclei_(byn in_red_all_in green)_labelled/";
Dirs{4} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_movies_ECadh_for_PIV/";
Dirs{5} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_nuclei/";
Dirs{6} = getenv('HOME')+"/Dropbox/GermbandRetraction/cable_ablation/";



% Loading images from that Directory
dir_i = 1; % corresponding to which number of Dirs above
j = 1; % corresponding to the Expt number inside this Dir
k = 1; % which set for the current experiment
% creating the Directory name which should be loaded
if (dir_i == 3)
    SubDir = "Expt-" + int2str(j) + "/Total";
else
    SubDir = "Expt-" + int2str(j);
end

nematicfile = char(Dirs{dir_i} + SubDir + "/ROI/ROI_CircleFull_Mid/" + sprintf('%d', k) + "/"...
               +  "NematicInfo.dat");
            
%nematicfile = char(Dirs{dir_i} + SubDir + "/ROI/OrientationJ/" + sprintf('%d', k) + "/"...
%                +  "NematicInfo.dat");
% The starting and the ending frame for piv

%nematicfile = char(Dirs{dir_i} + SubDir + "/Original_Detach_Big/NematicInfo.dat");

NematicData = dlmread(nematicfile);
%NematicData = dlmread('Nematic_Info_Expt-1_Median.dat');
% Images directory and reading files from there
imDir = char(Dirs{dir_i} + SubDir + "/ROI/ROI_CircleFull_Mid/" + sprintf('%d', k) + "/");
%imDir = char(Dirs{dir_i} + SubDir + "/ROI/ROI_Full/" + sprintf('%d', k) + "/");
%imDir = char(Dirs{dir_i} + SubDir + "/Original_Detach_Big/");


% Create list of images inside specified directory
directory = imDir; % this is the directory from which the experiments are to be loaded
%directory=uigetdir; %directory containing the images you want to analyze
suffix='*.tif'; %*.bmp or *.tif or *.jpg or *.tiff or *.jpeg
direc = dir([directory,filesep,suffix]); filenames={};
[filenames{1:length(direc),1}] = deal(direc.name);
imfilenames = sortrows(filenames); %sort all image files

% time and space steps
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
frames = max(NematicData(:,3)) + 1;
seq = 1:length(NematicData(:, 3));
for i = 1:frames
   imfile = strcat(imDir, imfilenames{i});
   im = imread(imfile);
   a = (NematicData(:, 3) == i-1);
   x = NematicData(seq(a),1);
   y = NematicData(seq(a),2);
   
   xuniq = unique(x);
   dx = xuniq(2) - xuniq(1);
   dy = dx;
   %y = max(y) - y;
   px = NematicData(seq(a),4); % px 
   py = NematicData(seq(a),5); % py
   C = NematicData(seq(a),7).*NematicData(seq(a), 8);
   %C = NematicData(seq(a), 7); % coherence
   figure(1)
   imshow(im);
   set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
   pause(0.2);
   hold on
   quiver(x, y, C.*px, C.*py, 'y', 'ShowArrowHead', 'off', 'LineWidth', 2);
   quiver(x, y, -C.*px, -C.*py, 'y', 'ShowArrowHead', 'off', 'LineWidth', 2);
   %pause(0.5);
   comb_imfile = char(imDir + "Nematic_Window_" + sprintf('%04d', i) + ".tif");
   %export_fig(comb_imfile)
   saveas(gcf, imfile, 'tif');
   axis equal
   hold off
   close
end
%close all;

%% Now doing processing of the data

a = (NematicData(:, 3) == i-1);
x = NematicData(seq(a),1);
y = NematicData(seq(a),2);
[xmesh, ymesh] = meshgrid(unique(x), unique(y));
[ny, nx] = size(xmesh);

Coh = NematicData(:, 7);
Cave = zeros(ny, nx);
Cmean_time = zeros(ny,1);

% including color scheme
color = gray(frames);

figure(2)
for i = 1:frames
   imfile = strcat(imDir, imfilenames{i});
   im = imread(imfile);
   a = (NematicData(:, 3) == i-1);
   x = NematicData(seq(a),1);
   y = NematicData(seq(a),2);
   %y = max(y) - y;
   px = NematicData(seq(a),4); % cos phi
   py = NematicData(seq(a),5); 
   %C = NematicData(seq(a),7).*NematicData(seq(a), 8);
   C = NematicData(seq(a), 7);
   Qxx = C.*(2*px.^2 - 1);
   Qxy = C.*(2.*px.*py);
   
   [xmesh, ymesh] = meshgrid(unique(x), unique(y));
   Cgrid = griddata(x, y, C, xmesh, ymesh);
   pxgrid = griddata(x, y, px, xmesh, ymesh);
   Qxx = griddata(x, y, Qxx, xmesh, ymesh);
   Qxy = griddata(x, y, Qxy, xmesh, ymesh);
   Q = sqrt(Qxx.^2 + Qxy.^2)
   
   Px = Cgrid.*abs(pxgrid);
   Px = abs(pxgrid);
   Cmean = mean(Cgrid')';
   %Cmean = mean(Qxy')';
   %Cmean = mean(Q);
   %Cmean = mean(Qxx');
   Cmean_time = Cmean_time + Cmean;
   %Cmean = mean(Px')';
   Cave = Cave + Cgrid;
   yuniq = unique(y)
   %plot(ds*yuniq(1:end-1), Cmean(1:end-1), 'LineWidth', 0.1, 'Color', color(i, :));
   ymax = yuniq(end-1);
   ymin = yuniq(1);
   ymid = (ymax + ymin)/2;
   %plot(ds*(yuniq(1:end-1)-ymax), Cmean(1:end-1), 'LineWidth', 0.1, 'Color', color(i, :));
   plot(ds*(yuniq(1:end-1)-ymax), Cmean(1:end-1), 'LineWidth', 0.1, 'Color', [0.8, 0.8, 1.0000]);
   
   xlabel('y');
   ylabel('Q');
   %pause(0.3)
   hold on
   grid off
%    %S = surf(xmesh, ymesh, Cgrid);
%    %contourf(xmesh, ymesh, Cave)
%    figure(1)
%    imfile = strcat(imDir, imfilenames{i});
%    im = imread(imfile);
%    imshow(im);
%    hold on;
%    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
%    Nc = 5;
%    v = linspace(min(Coh(:)), max(Coh(:)), Nc)
%    [cc, hh] = contour(xmesh, ymesh, Cgrid, v, 'LineWidth', 1, 'Color', 'y');
%    %clabel(cc,hh, 'FontWeight','bold','Color','yellow');
%    view(2);
%    axis equal;
%    colorbar;
%    %show Colorbar;
%    %S.EdgeColor = 'none';
%    %S.FaceColor = 'interp';
%    pause(0.3);
%    hold off
end
%%
%close 
%figure(100)
hold on
%plot(ds*yuniq(1:end-1), Cmean_time(1:end-1)/frames, 'r', 'LineWidth', 3);

% for ecad start from the other end
% for myosin start from the middle
   
ymax = yuniq(end-1);
ymin = yuniq(1);
ymid = (ymax + ymin)/2;
qmin = min(Cmean_time(1:end-1)/frames);
qmax = max(Cmean_time(1:end-1)/frames);
plot(ds*(yuniq(1:end-1)-ymax), Cmean_time(1:end-1)/frames, 'b', 'LineWidth', 3);
xlabel('DV distance, y (\mum)');
ylabel('Anisotropy strength (Myosin)');
axis([-ds*ymax, ds, qmin-0.03, qmax + 0.03]); % for ecad-5
%axis([ds*(ymin - ymid)-dx/5, ds*(ymax - ymid)+dx/5, qmin-0.03, qmax + 0.03]); % for myo
%%
   %S = surf(xmesh, ymesh, Cave);
   %contourf(ds*xmesh, ds*ymesh(end:-1:1, :), Cave(end:-1:1, :)/frames);
   ymax = max(ymesh(:));
   ymin = min(ymesh(:));
   ymid = (ymax + ymin)/2;
   contourf(ds*xmesh, ds*(ymesh(end:-1:1, :)-ymid), Cave/frames);
   view(2);
   colorbar
   %show Colorbar;
   %S.EdgeColor = 'none';
   %S.FaceColor = 'interp';
   pause(0.1);
   xlabel('AP distance, x (\mum)');
   ylabel('DV distance, y (\mum)');
   axis equal;
   hold on
%% Now doing processing of the data for equivalent Nematic
frames = max(NematicData(:,3));
seq = 1:length(NematicData(:, 3));
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
for i = 1:frames
   imfile = strcat(imDir, imfilenames{i});
   im = imread(imfile);
   
   a = (NematicData(:, 3) == i);
   x = NematicData(seq(a),1);
   y = NematicData(seq(a),2);
   %y = max(y) - y;
   px = NematicData(seq(a),4);
   py = NematicData(seq(a),5);
   C = NematicData(seq(a),7);
   
   Px = C.*px;
   Py = C.*py;
   
   Pxmean = mean(Px(:));
   Pymean = mean(Py(:));
   
   Cx = mean(x(:));
   Cy = mean(y(:));
   
   imshow(im);
   hold on
   sc = 100
   quiver(Cx, Cy, sc*Pxmean, sc*Pymean, 'y', 'ShowArrowHead', 'off', 'LineWidth', 5);
   quiver(Cx, Cy, -sc*Pxmean, -sc*Pymean, 'y', 'ShowArrowHead', 'off', 'LineWidth', 5);
   pause(0.5);
   hold off
   
end
close all;