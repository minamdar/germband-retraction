clear all;
close all;


%NematicFile = "Nematic_Info_Expt-5_Cad_end_23.dat";
NematicFile = "Nematic_Info_Expt-1_Median.dat";
%NematicFile = "Nematic_Info_Expt-3_Median.dat";
%NematicFile = "Nematic_Info_Expt-3_Median_2022_02_26.dat";
%NematicFile = "Nematic_Info_Cad_5.dat";
%ImageFile = 'Median_ROI_Expt-3.tif';
%ImageFile = 'Median_ROI_Expt-3_2022_02_26.tif';

Ndata = load(NematicFile); % load nematic data file
%%
% In new version of OrientationJ the time starts from  0
Ndata(:, 3) = Ndata(:, 3);

% discretization dx and dy are obtained
dx = Ndata(2, 1)- Ndata(1,1);
dy = dx;

FrameBegin = 1; 
NFrames = max(Ndata(:,3) + 1); % +1 because time starts from zero


NdataCell = cell(NFrames,1); % create a cell to save data for every frame
PlusHalf = cell(NFrames, 1); % cell to save +1/2 defects
MinusHalf = cell(NFrames,1); % cell to save -1/2 defects

% fill up the cell array
for i = FrameBegin : (FrameBegin + NFrames-1)
    row = find( Ndata(:,3) == i-1 );
    NdataCell{i} = Ndata(row,1:end);
    YYY = Ndata(row, 2);
end


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





