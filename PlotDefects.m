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
j = 3; % corresponding to the Expt number inside this Dir
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
%%
% In new version of OrientationJ the time starts from  0
NematicData(:, 3) = NematicData(:, 3);

% discretization dx and dy are obtained
dx = NematicData(2, 1)- NematicData(1,1);
dy = dx;

FrameBegin = 1; 
NFrames = max(NematicData(:,3) + 1); % +1 because time starts from zero


NdataCell = cell(NFrames,1); % create a cell to save data for every frame
PlusHalf = cell(NFrames, 1); % cell to save +1/2 defects
MinusHalf = cell(NFrames,1); % cell to save -1/2 defects

% fill up the cell array
for i = FrameBegin : (FrameBegin + NFrames-1)
    row = find( NematicData(:,3) == i-1 );
    NdataCell{i} = NematicData(row,1:end);
    YYY = NematicData(row, 2);
end

%%
frames = max(NematicData(:,3)) + 1;
seq = 1:length(NematicData(:, 3));

for i = 1:frames
    
   [PlusHalf{i} MinusHalf{i}] = Index(NdataCell{i}, dx, dy);
    
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
   scl = 5;
   quiver(x, y, scl*px, scl*py, 'y', 'ShowArrowHead', 'off', 'LineWidth', 2, 'AutoScale', 'off');
   quiver(x, y, -scl*px, -scl*py, 'y', 'ShowArrowHead', 'off', 'LineWidth', 2, 'AutoScale', 'off');
   hold on
   plot(PlusHalf{i}(:,1), PlusHalf{i}(:,2),'gs', 'MarkerFaceColor', 'c', 'MarkerSize', 12);
   plot(MinusHalf{i}(:,1), MinusHalf{i}(:,2),'mo', 'MarkerFaceColor', 'm',  'MarkerSize', 12);
   %pause(0.5);
   comb_imfile = char(imDir + "Nematic_Window_" + sprintf('%04d', i) + ".tif");
   %export_fig(comb_imfile)
   saveas(gcf, comb_imfile, 'tif');
   axis equal
   hold off
   close
end
%close all;


%{


% provide the cell array information one frame at a time to Index function
% maybe be faster numerically to do that in the previous loop
% doing this, however, to prevent mental clutter

PlHf = [];
MnHf = [];
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
for i = FrameBegin:FrameBegin+NFrames-1
%for i = FrameBegin:1 
   [PlusHalf{i} MinusHalf{i}] = Index(NdataCell{i}, dx, dy);
   
   ds = 0.248;
   PlusHalf{i} = ds*PlusHalf{i};
   MinusHalf{i} = ds*MinusHalf{i};
   
   Data = NdataCell{i};
   xx = Data(:,1)*ds;
   yy = Data(:,2)*ds;
   %yy = max(yy) - yy;
   nx = Data(:,4);
   ny = Data(:,5);
   C = Data(:,7);
   E = Data(:,8);
   %figure('units','normalized','outerposition',[0 0 1 1]); 
   
   %ImageFile = char("./Img_Cad_5/ROI_" + sprintf('%04d', i) + ".tif");
   %ImageFile = ImDir + imagefiles(i).name
   %ImageFile = char("Median_ROI_Expt-1.tif");
   ImageFile = char("Median_ROI_Expt-3.tif");
   %ImageFile = 'Expt-3_Detachment_Big_180.tif';

   hold off
   %im = imread(ImageFile); % load that image
   %K = imadjust(im,[0.3 0.7],[]);
   %imshow(K);
   %imshow(ImageFile); % plot the image
   hold on;
   % plot the nematic directors
   % scale = 1.5*C*dx/2;
   scale = 1.5*dx/14;
   
   [XX, YY] = meshgrid(unique(xx(:)), unique(yy(:)));
   Nx = griddata(xx, yy, nx, XX, YY);
   Ny = griddata(xx, yy, ny, XX, YY);
   Em = griddata(xx, yy, C, XX, YY);
   Cm = griddata(xx, yy, E, XX, YY);
   wgt = Em.*Cm;
   wgt = wgt./wgt;

   
   %% Minimize error
   error_pm = @(alpha) err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, Nx, Ny);
   x0 = 0;
   alpha = fminsearch(error_pm, x0);
   %XX = XX(2:end-1, 2:end-1);
   %YY = YY(2:end-1, 2:end-1);
   %Nx = Nx(2:end-1, 2:end-1);
   %Ny = Ny(2:end-1, 2:end-1);
   
   sh = 0; % shift in dx
   col_or = [0.8500, 0.3250, 0.0980];
   col_r  = [1, 0, 0];
   % -dx/2 and -dy/2 below is to match up with the actual image
   %q = quiver(XX - sh*dx/2, YY - sh*dy/2, scale.*Nx.*1, scale.*Ny.*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off'); % quiver plot for the nematic defect
   xuniq = unique(XX(:));
   yuniq = unique(YY(:));
   
   xmid = (max(xuniq) + min(xuniq))/2;
   ymid = (max(yuniq) + min(yuniq))/2;
   
   q = quiver(XX-xmid - sh*dx/2, YY(end:-1:1, :) -ymid - sh*dy/2, scale.*Nx.*1, -scale.*Ny.*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off'); % quiver plot for the nematic defect

   set(q, 'color', col_r, 'linewidth', 2.);
   %q = quiver(XX - sh*dx/2, YY - sh*dy/2, -scale.*Nx.*1, -scale.*Ny.*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off');
   q = quiver(XX-xmid - sh*dx/2, YY(end:-1:1, :) -ymid - sh*dy/2, -scale.*Nx.*1, scale.*Ny.*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off');

   set(q, 'color', col_r, 'linewidth', 2.);
   axis equal;
   hold on
   %title('Nematic Directors (red), +1/2 (green square) and -1/2 (blue circles) defects');
   %title('+1/2 (green), -1/2 (blue)');
   %xlabel('X');
   %ylabel('Y');
   %alpha = -0.3;
   %figure(2)
   phi = defect(PlusHalf, MinusHalf, XX, YY, alpha);
   %q = quiver(XX - sh*dx/2, YY - sh*dy/2, scale.*cos(phi)*1, scale.*sin(phi).*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off'); % quiver plot for the nematic defect
   q = quiver(XX-xmid - sh*dx/2, YY(end:-1:1, :)-ymid - sh*dy/2, scale.*cos(phi)*1, -scale.*sin(phi).*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off'); % quiver plot for the nematic defect
    set(q, 'color', 'b', 'linewidth', 1.5);
   %q = quiver(XX - sh*dx/2, YY - sh*dy/2, -scale.*cos(phi).*1, -scale.*sin(phi).*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off');
   q = quiver(XX-xmid - sh*dx/2, YY(end:-1:1, :)-ymid - sh*dy/2, -scale.*cos(phi)*1, scale.*sin(phi).*1, 'LineStyle', '-', 'Autoscale', 'off', 'ShowArrowhead', 'off'); % quiver plot for the nematic defect
   set(q, 'color', 'b', 'linewidth', 1.5);
   axis equal;
   
   hold on
   %plot(XX(:), YY(:), 'bx');
   
   set(gca,'FontSize',18)
   set(findall(gcf,'type','text'),'FontSize',5)
   %xlim([0 - dx, max(xx(:))+ dx]);
   %ylim([0 - dy, max(yy(:))+ dy]);
   
   xlim([0 - sh*dx - xmid, max(xx(:))-xmid + sh*dx]);
   xlim([-30, 30]);
   ylim([0 - sh*dy-ymid, max(yy(:))-ymid + sh*dy]);
   ylim([-30, 30]);
   box on;

   hold on
   % plot the +1/2 and -1/2 defects
   lightBlue = [91, 207, 244] / 255
  
   %plot(PlusHalf{i}(:,1)-sh*dx/2, (PlusHalf{i}(:,2))-sh*dy/2,'gs', 'MarkerFaceColor', 'c', 'MarkerSize', 12);
   plot(PlusHalf{i}(:,1)-xmid-sh*dx/2, (max(yuniq) + min(yuniq)-PlusHalf{i}(:,2))-ymid-sh*dy/2,'gs', 'MarkerFaceColor', 'c', 'MarkerSize', 12);
   %plot(MinusHalf{i}(:,1)-sh*dx/2, (MinusHalf{i}(:,2))-sh*dy/2, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 12);
   plot(MinusHalf{i}(:,1)-xmid-sh*dx/2, (max(yuniq) + min(yuniq) - MinusHalf{i}(:,2))-ymid-sh*dy/2, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 12);
   %plot(xx, yy, 'ko');
   %EpsName = strcat('./ControlExperimentsDavid/Expt-16/NematicDefects/Nematic', int2str(i), '.tif');
   
   %NematicImage = char(NematicDir + "defects_" + sprintf('%04d', i) + ".png");
   NematicImage = char("./Defects_Cad_5/" + "defects_" + sprintf('%04d', i) + ".tif");
   %export_fig(NematicImage, '-tif', '-m2.5');
   %saveas(gca, EpsName,'tif');
   %pause(0.1)
   hold off
  
   pause(1)
   %close
end
%%
er_a = err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, Nx, Ny);
%{
alpha = -0.3;
%er_a = err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, Nx, Ny);
error_pm(0.1)
%%
al = 0.1;
err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, Nx, Ny)
%%
   %% Minimize error
   error_pm = @(alpha) err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, Nx, Ny);
   x0 = 0;
   alpha = fminsearch(error_pm, x0);
%}
%%
function phi = defect(PlusHalf, MinusHalf, XX, YY, alpha)
phi = 0*XX;
Ph = PlusHalf{1};
Mn = MinusHalf{1};
for i = 1:length(Ph(:, 1))
    phi = phi + 0.5*atan2(YY - Ph(i, 2), XX - Ph(i, 1));
end
for i = 1:length(Mn(:, 1))
    phi = phi - 0.5*atan2(YY - Mn(i, 2), XX - Mn(i, 1));
end
phi = phi + alpha;
end

function er = err(PlusHalf, MinusHalf, XX, YY, alpha, wgt, nx, ny)

    phi = defect(PlusHalf, MinusHalf, XX, YY, alpha);
    px = cos(phi);
    py = sin(phi);
    errmat = wgt.*(nx.*px + ny.*py).^2;
    er = 1- sum(errmat(:))/sum(wgt(:));
end

%}

