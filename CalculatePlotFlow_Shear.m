clear all;
close all;
%{
1: Myosin Data - Contains: Expt-1/2/3
2: Control eCadherin - Contains: Expt-1/2/3/4/5
%}

Dirs = cell(5, 1);

%Dirs{1} = getenv('HOME')+"/Dropbox/GermbandRetraction/myosin_orientationJ/";
Dirs{1} = getenv('HOME')+"/Dropbox/Research/Sudeepa/myosin_orientationJ/";
%Dirs{2} = getenv('HOME')+"/Dropbox/GermbandRetraction/control_ECadherin_movies_for_PIV_and_shapes_analysis/";
Dirs{2} = getenv('HOME')+"/Dropbox/Research/Sudeepa/control_ECadherin_movies_for_PIV_and_shapes_analysis/";
Dirs{3} = getenv('HOME')+"/Dropbox/GermbandRetraction/control_nuclei_(byn in_red_all_in green)_labelled/";
Dirs{4} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_movies_ECadh_for_PIV/";
Dirs{5} = getenv('HOME')+"/Dropbox/GermbandRetraction/SqhAA_nuclei/";
Dirs{6} = getenv('HOME')+"/Dropbox/GermbandRetraction/cable_ablation/";
Dirs{7} = "/Users/minamdar/Work/GitRepos/Gitlab/germband-retraction/Codes/TestCodes/";

% loading the directories

iex = 1; % which of the experimental conditions
jex = 3; % which particular experiment 
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_detach/";
DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_detach_big_16_new/";
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_detach_big/";
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_ROI_check/";
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_ROI/";
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/piv_ROI_New/";
%ImDir = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/ROI/ROI_Full/1/";
%ImDir = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/ROI/ROI_CircleFull_Mid/1/";
%ImDir = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/ROI/ROI_Full/2/";
%ImDir = Dirs{iexe} + "Expt-" + sprintf('%d', jex) + "/Original_Detach/";
ImDir = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/Original_Detach_Big/";
%DirName = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/PIV/";
%ImDir = Dirs{iex} + "Expt-" + sprintf('%d', jex) + "/SimpleShear/";

% space from px to \mu m
dsp = cell(5, 2);
dsp{1}{1} = 0.248;
dsp{1}{2} = 0.207;
dsp{1}{3} = 0.248;
dsp{2}{2} = 0.248;
dsp{2}{3} = 0.262;
dsp{2}{5} = 0.262;
dsp{2}{6} = 0.262;
dsp{7}{1} = 1; % test data
dsp{7}{2} = 1; % test data
dsp{7}{3} = 1; % test data
dsp{7}{4} = 1; % test data
ds = dsp{iex}{jex};

%%
dtime = cell(5, 2);
dtime{1}{1} = 0.442;
dtime{1}{2} = 0.456;
dtime{1}{3} = 0.165;
dtime{2}{2} = 1.37;
dtime{2}{3} = 0.626;
dtime{2}{5} = 0.575;
dtime{2}{6} = 0.575;
dtime{7}{1} = 1; % test data
dtime{7}{2} = 1; % test data
dtime{7}{3} = 1; % test data
dtime{7}{4} = 1; % test data

dt = dtime{iex}{jex};
%DirName = "./OpticalFlow/";

% Convert .mat file into bunch of 'PIVlab_xxxx.txt' files
%FrmEnd = Mat2txt(DirName, 'PIVData.mat');

%%
% opening all the files
FrmBgn = 1;
FrmEnd = 177;
%FrmEnd = 23;
kernsize = 0.;

vfil = CalcFlow(DirName, FrmBgn, FrmEnd, kernsize);
vid=VideoWriter(char(Dirs{1}{1}+ "Flow.avi"));
%showf(vfil, 'duxdx', 'ScaleArrow', 3, 'surf', 'Spacing', [1, 1],'savevideo', vid);
%showf(vfil, 'curl', 'ScaleArrow', 3, 'surf', 'Spacing', [1, 1]);
view(2);

% Getting derived quantities
vx = vec2scal(vfil, 'ux'); 
vy = vec2scal(vfil, 'uy');
curl = vec2scal(vfil, 'curl'); % curl
% getting the general gradient terms for getting pure shear rate or shear
% strain rate
duxdx = vec2scal(vfil, 'duxdx'); 
duxdy = vec2scal(vfil, 'duxdy');
duydx = vec2scal(vfil, 'duydx'); 
duydy = vec2scal(vfil, 'duydy'); 

% just get dimensions
len_y = length(duydy(1).y);
%%
% pure shear
nfrms = FrmEnd - FrmBgn + 1;
uxyavg = zeros(1, len_y);
u_mag_avg = zeros(1, len_y);
uxx_yy_avg = zeros(1, len_y);
div_avg = zeros(1, len_y); % for storing divergence

%uxyavg = zeros(1, 36);% expt-2 ecad
%%uxyavg = zeros(1, 22);% expt-5 ecad
%uxyavg = zeros(1, 20);% expt-6 ecad
%uxyavg = zeros(1, 30);% expt-6 ecad detach
%uxyavg = zeros(1, 38); % expt-3 myo detach
%uxyavg = zeros(1, 18); % expt-1 myo ROI
%uxyavg = zeros(1, 28);% expt-3 myosin


frames = FrmEnd;
Vxxmean = zeros(1, frames);
Vxymean = zeros(1, frames);
Vyymean = zeros(1, frames);
color = gray(length(vfil));

for i = 1:length(vfil)
    uxx = duxdx(i).w/dt;
    uxy = duxdy(i).w/dt;
    uyx = duydx(i).w/dt;
    uyy = duydy(i).w/dt;
    velx = vfil(i).vx/dt;
    vely = vfil(i).vy/dt;
    x =  duxdx(i).x;
    y = duydx(i).y;
    
    v_mag = sqrt((uxx - uyy).^2/4 + (uxy + uyx).^2/4);
    v_xy = 0.5*(uxy + uyx);
    v_xx_yy = 0.5*(uxx - uyy);
    div = (uxx + uyy);
    %shear = 0.5*(uxx - uyy);
    %shear = v;
    %shear = 0.5*(uxy + uyx);
    %v = 0.5*(uxy + uyx);
    %v = velx;
    %shear = uxy;
    Im  = char(ImDir + "ROI_" + sprintf('%04d', i) + ".tif");
 
   %{
    fh = figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.75]);
    imshow(Im);
    hold on;
    [xg, yg] = meshgrid(x, y);

    theta = 0.5*atan2((uxy + uyx)/2, (uxx - uyy)/2);
    scl = 5;
    quiver(xg, yg, scl*v'.*cos(theta'), scl*v'.*sin(theta'),'r', 'ShowArrowhead', 'off', 'LineWidth', 1);
    hold on
    quiver(xg, yg, -scl*v'.*cos(theta'), -scl*v'.*sin(theta'), 'r', 'ShowArrowhead', 'off', 'LineWidth', 1);
    %quiver(xg, yg, velx', vely', 'y', 'LineWidth', 1);
    axis equal 
    file = char("im_" + sprintf('%04d', i) + ".tif");
    frm = getframe( fh )
      
    %imwrite(frm.cdata, file);
 
    pause(0.5);
    hold off
    
    %}
    
    figure(9)
    ymid = (max(y) + min(y))/2.;
    ymax = max(y(1:end-1));
    meanv_mag = mean(v_mag)
    %plot((y(1:end-1)-ymax)*ds, meanv_mag(1:end-1) , 'LineWidth', 0.1, 'Color', color(i,:));
    %plot((y(1:end-1)-ymax)*ds, meanv_mag(1:end-1) , 'LineWidth', 0.1, 'Color', color(17,:));
    plot((y(1:end)-ymid)*ds, meanv_mag(1:end) , '-', 'LineWidth', 0.1, 'Color', [1.0000, 0.8, 0.8]);
    %xlabel('y (\mum)');
    %ylabel('pure shear magnitude (min^{-1})');
    %pause(1)
    hold on
    u_mag_avg = u_mag_avg + mean(v_mag);
    
    figure(10)
    ymid = (max(y) + min(y))/2.;
    meanv_xy = mean(v_xy);
    %plot((y)*ds, mean(v_xy), 'LineWidth', 0.05, 'Color', color(i, :));
    plot((y(1:end)-ymid)*ds, meanv_xy(1:end) , '-', 'LineWidth', 0.1, 'Color', [0.8, 0.8, 0.8]);
    %plot((y-ymid)*ds, mean(v_xy), 'LineWidth', 0.1, 'Color', [0, 0, 0] + i/frames);
    %xlabel('y (\mum)');
    %ylabel('0.5(v_{xy} + v_{yx}) (min^{-1})');
    %pause(1)
    hold on
    uxyavg = uxyavg + mean(v_xy);
    
    figure(11)
    ymid = (max(y) + min(y))/2.;
    meanv_xx_yy = mean(v_xx_yy);
    %plot((y-ymid)*ds, uxx_yy_mean, 'r', 'LineWidth', 2)
    %plot((y)*ds, mean(v_xx_yy), 'LineWidth', 0.1, 'Color', color(i, :));
    plot((y(1:end)-ymid)*ds, meanv_xx_yy(1:end) , '-', 'LineWidth', 0.1, 'Color', [0.8, 0.8, 0.8]);
    %xlabel('y (\mum)');
    %ylabel('0.5(v_{xx} - v_{yy}) (min^{-1})');
    %pause(1)
    hold on
    uxx_yy_avg = uxx_yy_avg + mean(v_xx_yy);
    
    
    figure(21)
    ymid = (max(y) + min(y))/2.;
    ymax = max(y(1:end-1));
    mean_div = mean(div);
    %plot((y-ymid)*ds, mean_div, 'r', 'LineWidth', 2)
    %plot((y)*ds, mean(v_xx_yy), 'LineWidth', 0.1, 'Color', color(i, :));
    plot((y(1:end-1)-ymax)*ds, mean_div(1:end-1) , '-', 'LineWidth', 0.1, 'Color', [0.8, 0.8, 0.8]);
    %xlabel('y (\mum)');
    %ylabel('0.5(v_{xx} - v_{yy}) (min^{-1})');
    %pause(1)
    hold on
    div_avg = div_avg + mean(div);
       
    Vxxmean(i) = mean(uxx(:));
    Vxymean(i) = mean(uxy(:) + uyx(:))*0.5;
    Vyymean(i) = mean(uyy(:));
    
end

uxx_yy_mean = uxx_yy_avg/length(vfil);
uxymean = uxyavg/length(vfil);
u_mag_mean = u_mag_avg/length(vfil);
div_mean = div_avg/length(vfil);

%%
figure(9)
ymax = max(y(1:end-1));
%axes('YAxisLOcation', 'right');
plot((y(1:end-1)-ymid)*ds, u_mag_mean(1:end-1), 'r-', 'LineWidth', 2);
plot((y(1:end)-ymid)*ds, u_mag_mean(1:end), 'r-', 'LineWidth', 2)
xlabel('DV distance $y({\rm \mu m})$', 'interpreter', 'latex', 'FontSize', 18);
ylabel('pure shear magnitude (min$^{-1}$)', 'interpreter', 'latex', 'FontSize', 18);
%
set(gca, 'ycolor', 'r')
%set(gca,'FontSize',18, 'FontName', 'Times')
hold off;

figure(10)
ymid = (max(y) + min(y))/2.;
%plot((y)*ds, uxymean, 'r', 'LineWidth', 2)
plot((y(1:end)-ymid)*ds, uxymean(1:end), 'r', 'LineWidth', 2)
xlabel('DV distance $y{\rm (\mu m)}$', 'interpreter', 'latex', 'FontSize', 18);
%ylabel('$x$');
%ylabel('off $\frac{x}{y}$');
ylabel('off diagonal shear $\left \langle \frac{v_{xy} + v_{yx}}{2} \right \rangle_{x}~({\rm min}^{-1})$', 'interpreter', 'latex', 'FontSize', 18);
%set(gca,'FontSize',18, 'FontName', 'Times')
hold off

figure(11)
ymid = (max(y) + min(y))/2.;
%plot((y)*ds, uxx_yy_mean, 'r', 'LineWidth', 2)
plot((y(1:end)-ymid)*ds, uxx_yy_mean(1:end), 'r', 'LineWidth', 2)
xlabel('DV distance $y({\rm \mu m})$', 'interpreter', 'latex', 'FontSize', 18);
ylabel('diagonal shear $\left \langle \frac{v_{xx} - v_{yy}}{2} \right \rangle_{x}~({\rm min}^{-1})$ ', 'interpreter', 'latex', 'FontSize', 18);
%set(gca,'FontSize',18, 'FontName', 'Times')
hold off


figure(21)
ymid = (max(y) + min(y))/2.;
ymax = max(y(1:end-1));
plot((y(1:end-1)-ymax)*ds, div_mean(1:end-1), 'r', 'LineWidth', 2)
%plot((y(1:end)-ymid)*ds, div_mean(1:end), 'r', 'LineWidth', 2)
xlabel('DV distance $y({\rm \mu m})$', 'interpreter', 'latex', 'FontSize', 18);
ylabel('divergence $\left \langle v_{xx} + v_{yy} \right \rangle_{x}~({\rm min}^{-1})$ ', 'interpreter', 'latex', 'FontSize', 18);
%set(gca,'FontSize',18, 'FontName', 'Times')
hold off

%%
%figure(12)
%plot(y, uxx_yy_mean, 'r', 'LineWidth', 2)
%xlabel('y');
%ylabel('diagonal shear');

%figure(2)
%hold on
%plot(y, uxymean, 'r','LineWidth', 3);
%xlabel('y');
%ylabel('$\left \langle \frac{v_{xx} - v_{yy}}{2} \right \rangle_{x, t}$', 'interpreter', 'latex', 'FontSize', 18);
%ylabel('pure shear magnitude');
%ylabel('v_{xy}')
%

%dt = 0.442;% Expt-1 Myosin
%dt = 1.37; % Expt-2 Ecad
%dt = 0.165; % Expt-3 Myosin
%dt = 0.575; % Expt-6 Ecad

time = (0:length(Vxxmean))*dt;

hold off
figure(5)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
plot(time, [0, 0.5*cumsum(Vxxmean-Vyymean)]*dt);
xlabel('time (min)', 'interpreter', 'latex', 'FontSize', 18)
%ylabel('$\left \langle \int_t \frac{v_{xx} - v_{yy}}{2} dt\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18);
ylabel('time cumulative diagonal shear $\left \langle \frac{v_{xx} - v_{yy}}{2}\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18)

hold off
figure(6)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
plot(time, [0, 1*cumsum(Vxymean)]*dt);
xlabel('time (min)', 'interpreter', 'latex', 'FontSize', 18)
ylabel('time cumulative off-diagonal shear $\left \langle \frac{v_{xy} + v_{yx}}{2}\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18)
%ylabel('$\left \langle \int_t \frac{v_{xy} + v_{yx}}{2} dt\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18);


hold off
figure(15)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
plot(time, [0, cumsum(Vxxmean + Vyymean)]*dt);
xlabel('time (min)', 'interpreter', 'latex', 'FontSize', 18)
ylabel('$\left \langle \int_t \frac{v_{xx} + v_{yy}}{1} dt\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18);
%%

hold off
figure(16)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
plot(time, [0,cumsum(Vxxmean)]*dt);
xlabel('time (min)')
ylabel('$\left \langle \int_t \frac{v_{xx}}{1} dt\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18);

hold off
figure(17)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
plot(time, [0, cumsum(Vyymean)]*dt);
xlabel('time (min)')
ylabel('$\left \langle \int_t \frac{v_{yy}}{1} dt\right \rangle_{x, y}$', 'interpreter', 'latex', 'FontSize', 18);

%figure(11)
%Vmean = sqrt((Vxxmean-Vyymean).^2/4 + (Vxymean).^2/4);
%plot(cumsum(Vmean));
%%
%a = CalcFlow(DirName, FrmBgn, FrmEnd, 0.*kernsize)
function vfil = CalcFlow(DirName, FrmBgn, FrmEnd, kernsize)
% Starting and ending frame

    dt = 1.; % in terms of hours
    scaling = 1.; % default scaling between px and \mu-m is unity
    scaling = 1.; % scaling 1 without pixels 

    j = 1;
    for i = FrmBgn:FrmEnd
        %FileName = strcat('xyuv_', sprintf('%04d', i), '.txt');
        FileName = strcat('PIVlab_', sprintf('%04d', i), '.txt');
        %FileName = strcat('VelData_', sprintf('%04d', i), '.txt');
        File = strcat(DirName, FileName);
        v(j) = loadpivlabtxt(File, scaling, dt);
        File
        %v(j) = loadvec(File);
        j
        j = j + 1;
    end

    vsfil = filterf(v, kernsize, 'gauss', 'valid'); % space filter
    Ntime = 1; % steps for filtering in time
    vfil = smoothf(vsfil, Ntime); % filter after both space and time in 'x' direction
    %vfil = v;
    %vrfil = rotatef(vfil, pi/2); % field rotated by 90 deg for usage while computing spatial correlation in 'y' direction
end





