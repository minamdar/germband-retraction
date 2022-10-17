% This function is to obtain the winding number of the nematic defect
function [PlusSource MinusSource] = Index(Frame, dx, dy)

D = Frame; % Data for a given frame

xx = D(:,1); % x position
yy = D(:,2); % y position
maxY = max(yy); % top co-ordinate of yy
%yy = maxY - yy; % this is to modify the up down symmtery
nx = D(:,4); % x component of nematic unit director
ny = D(:,5); % y component of nematic unit director
Coherency = D(:,7); % coherency of the image small value implies almost isotropic
E = D(:,8); % energy or the trace of the shape tensor


% if the nx and ny vector is not unit
n = sqrt(nx.^2 + ny.^2);
nx = nx./n;
ny = ny./n;


clear D; % get rid of the big D array

SF = 10.; % factor for plotting arrows in quiver

% define a grid as per the experimental data
xGrid = min(xx):dx:max(xx);
yGrid = min(yy):dy:max(yy);

[xMesh, yMesh] = meshgrid(xGrid, yGrid);

% reformat all the points on the grid points systematically
nxMesh = griddata(xx, yy, nx.*Coherency, xMesh, yMesh);
nyMesh = griddata(xx, yy, ny.*Coherency, xMesh, yMesh);
x1Mesh = griddata(xx, yy, nx, xMesh, yMesh);
x2Mesh = griddata(xx, yy, ny, xMesh, yMesh);
CohMesh = griddata(xx, yy, Coherency, xMesh, yMesh);

% maximum coherency
size(CohMesh);
CohMax = max(CohMesh(:));
CohMax;

% get the extent of the grid in the horizontal (x) and vertical (y)
% direction
[yMax xMax] = size(xMesh);
nContours = (xMax - 1)*(yMax - 1); 
Contours = cell(nContours,1);

% The for loop below to set up the connectivity of the contours 
% around which we will obtain the index
ind = 1;
for j = 1:xMax-1
    for i = 1:yMax-1           
         int1 = (j-1)*yMax + (i + 1); % first node of the contour
         int2 = (j-1)*yMax + i; % second node
         int3 = j*yMax + i; % third node
         int4 = j*yMax + (i + 1); % fourth node
         Contours{ind} = [int1 int2 int3 int4]; % This is the loop
         ind = ind + 1;
     end
end
% Now obtain the winding number for each loop
% Each loop corresponds to the smallest rectangle 

eZ4 = [0 0 0 0;0 0 0 0;1 1 1 1]; % unit vector in the 'z' direction 
Index = zeros(nContours,1); % Array to obtain the index for every loop



% This for-loop obtains the index for all the smallest rectangles   
    for i = 1:nContours
        Rect = Contours{i}; % Recalling the ith contour
        Vec1 = [x1Mesh(Rect);x2Mesh(Rect);zeros(1,4)];
        Vec2 = Vec1(:,[2 3 4 1]);
        Coh = [CohMesh(Rect)];
        CrossP = cross(Vec1', Vec2'); % The cross-product decides the sense
        DotP = dot(Vec1, Vec2); % The dot product gives the smallest angle
        WindMat = acos(abs(DotP')).*sign(CrossP(:,3)).*sign(DotP');% change in rotation
        Index(i) = 1/(2*pi)*sum(WindMat);% sum of all the changes 
        max(abs(WindMat));
        
        % check if the coherence is too small and if yes reject the defect
        loccoh = max(Coh(:));
       % if (loccoh < 0.01*CohMax)
       %     'rejected defect'
       %     Index(i) = 0;
       % end
        % if the maximum change in angle is more than 80 deg discard that
        %if(max(abs(WindMat)) > 90*pi/180. && abs(Index(i)) > 0.1)
        %    Index(i) = 0;
        %end
    end
    

% The Index will be +1/2 or -1/2 if the disclination defect is present
% The Index will be +1 or -1 if aster or vortex defect is present
% Else the Index will be zero


% Find + 1/2 defects
 clear a;
 size(Index)
 a = find(Index > 0.4);
 Index(a);
 PlusSource = zeros(length(a),2); % Storing center of source
 clear rect;
for i = 1:length(a)    
    rect = Contours{a(i)};
    %rect = [rect rect(1)];
    xMesh(rect)
    % The center of the rectangle containing the +1/2 defect
    PlusSource(i, 1) = mean(xMesh(rect));
    PlusSource(i, 2) = mean(yMesh(rect));
end

% Finding -1/2 defect
clear b;
clear rect;
b = find(Index < -0.4);
Index(b);
MinusSource = zeros(length(b),2);
for i = 1:length(b)
    hold on
    rect = Contours{b(i)};
    %rect = [rect rect(1)]
    % The center of the rectangle containing the -1/2 defect
    xMesh(rect);
    yMesh(rect);
    MinusSource(i, 1) = mean(xMesh(rect));
    MinusSource(i, 2) = mean(yMesh(rect)); 
end

end



