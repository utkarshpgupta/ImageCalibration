function [T_est,R_est,fxe,fye, oc, or]= calibrateCamera(imageFilename)
% Have you ever wondered how a photographer got a shot? This code computes
% the camera's original position and orientation from a single photo -- all
% you need to do is measure the real-world coordinates of points seen in
% the photograph and then click on those points in the photo. The camera's
% estimated position is T_est with orientation R_est.
% The code also determines the camera's focal length divided by the 
% scaling factor for the pixels in x and y to get fxe and fye, and the
% optical center (oc,or).
% 
% All you need is a photograph (imageFilename), and to update
% gather3Ddata.m with the coordinates of points in the 3D world visible in
% the photo.
%
% Instruction video: https://youtu.be/WEYwitb6dTo
% 
% This uses the method in Robot Modeling and Control by Spong, Vidgasgar
% and Hutchinson.  This methods is also known as Tsai's calibration method
% http://people.csail.mit.edu/bkph/articles/Tsai_Revisited.pdf Tsai, Roger
% Y. (1986) "An Efficient and Accurate Camera Calibration Technique for 3D
% Machine Vision", Proceedings of IEEE Conference on Computer Vision and
% Pattern Recognition, Miami Beach, FL, 1986, pp. 364-374.
% 
%  By Steban Soto, Cora Yuzhu, Aaron T. Becker in 2018, updated with mouse
%  handler by Shreyas Poyrekar in 2018.
% 
%  This code requires an image in the working directory named
%  'imageFilename.jpg' It also requires the function gather3Ddata.m and
%  findOrthocenter.m
%
% TODOs for assignment:
% 	1. Find an area with many orthogonal lines.  Physically add a world
% 	coordinate frame somewhere (I drew a right angle on a sheet of paper
% 	and called these x^W and y^W).
% 	
%   2. Measure the world frame coordinates of at least 20 points in your
%   image. (I chose the corners of rectangles). 
%
% 	3. Take a picture of your world, and save it as a .jpg.  YOU MUST BE
% 	VISIBLE IN YOUR IMAGE.
%
% 	4. Add your 3D world coordinates to gather3Ddata.m.
%
% 	5. Run calibrateCamera.m to learn the intrinsic and extrinsic
% 	parameters for your camera.
%


format compact

if nargin<1
    imageFilename = 'CameraCalP1.jpg';%'hw4pic.jpg';%'CameraCalP1.jpg';
end
close all
%% Load (or generate) pairs of camera pixel locations and 3D data
% http://people.csail.mit.edu/bkph/articles/Tsai_Revisited.pdf
% try to load (r,c) values if the file exists.
if exist([imageFilename(1:end-4),'Data.mat'],'file')
    load CameraCalP1Data.mat  rcVal  DataXYZ DataW 
else % otherwise have user click on them
    [rcVal,DataXYZ, DataW] = gather3Ddata(imageFilename);
    %rc values is x (points to right) and y (points down)
end
% try to load (or,oc) values if the file exists.
if exist([imageFilename(1:end-4),'Ortho.mat'],'file')
    load([imageFilename(1:end-4),'Ortho.mat'],'orthoc') 
else  % otherwise have user generate the orthocenter
    orthoc = findOrthocenter(imageFilename);
end

%% load all the data values
img = imread(imageFilename);
or = orthoc(2); %find these from orthocenter %x
oc = orthoc(1); %find these from orthocenter %y
% convert (r,c) values into r' and c' values
%  or = size(img,1)/2; oc = size(img,2)/2;
rP = -rcVal(:,2)+or;  %x
cP =  rcVal(:,1)-oc;  %y
xW = DataXYZ(:,1);
yW = DataXYZ(:,2);
zW = DataXYZ(:,3);

%% Solve for the missing parameters
% pg 384 -> A matrix
A = [rP.*xW, rP.*yW, rP.*zW, rP, -cP.*xW, -cP.*yW, -cP.*zW, -cP];
%A = A(10:20,:);  % worse results when you use a subset

% pg 384 -> solve for xbar, using SVD
% xbar = k*[r21 r22 r23 Ty alpha*r11 alpha*r12 alpha*r13 alpha*Tx]
[U,S,V] = svd(A);   %#ok<ASGLU>
xbar = V(:,end); %grab the last column of V to get the solution corresponding to the smallest singular value

% pg 384 -> solve for abs(k) and alpha
kAbs = sum(xbar(1:3).^2).^0.5;
alpha = (sum(xbar(5:7).^2).^0.5)/kAbs;

% pg 384 -> find sign for k
% Choose k such that r/(r11*x+r12*y+r13*z+Tx) > 0 NOTE: THIS IS DIFFERENT
% FROM BOOK.  SHOULD BE DIVIDE, and GREATER THAN
if sum(rP./( xbar(5).*xW +  xbar(6).*yW +  xbar(7).*zW +xbar(8)) ) > 0 
    k = -kAbs;
else
    k = kAbs;
end

% pg 384-385 -> Calculate R_est, col 3 calculated by cross product of
% col 1 and col 2
r1all = xbar(5:7)/(alpha*k);
r2all = xbar(1:3)/k;
r3all = cross(r1all,r2all);  % the 3rd row of any rotation matrix is the cross product of the first two rows.
R_est = [r1all,r2all,r3all]';
% RENORMALIZE   (rotation matrices should be orthogonal, but numeric errors may have have been introduced when you measured your points.)
R_est = orthogonalize(R_est);

% pg 384 -> use k and alpha to find Ty and Tx  
Ty = xbar(4)/k;  
Tx = xbar(8)/(alpha*k); 

% pg 385 -> find Tz and Fx
C = [rP, (R_est(1,1)*xW + R_est(1,2)*yW + R_est(1,3)*zW + Tx)];
d = -rP.*(R_est(3,1)*xW + R_est(3,2)*yW + R_est(3,3)*zW );

ybar = C\d;  %Solve C.ybar=d, ybar = [Tz,fx]

%maxErr = max(abs(C*ybar- d));
%display(maxErr)

Tz  = ybar(1);
fxe = ybar(2);

% pg 385 -> calculate fy estimate, using alpha = f_x / f_y
fye = fxe/alpha; 

% pg 381 -> calculate T_est, using T = -R_w^c * O_c^w  from 11.2.1
T_est =-R_est'*[Tx;Ty;Tz]; 

% Print out errors (TODO)
display([or,oc])

save([imageFilename(1:end-4),'Param'],'T_est','R_est','fxe','fye', 'oc', 'or','imageFilename')

% draw the camera in the 3D scene
% draw 3D world
f2 = figure(2); clf;
set(gca,'fontsize', 18);
set(f2,'name', '3D and Camera reconstruction');
for j = 1:size(DataW,1)/4
    recCoords = cell2mat(DataW(4*(j-1)+(1:4),1:3));
    label = DataW(4*(j-1)+1,4);
    
    label = strsplit(cell2mat(label),',');
    label = label(1);
    p=patch(recCoords(:,1),recCoords(:,2),recCoords(:,3),rand(1,1));
    set(p,'FaceColor','flat');
    hold on
    t=text( mean(recCoords(:,1)),mean(recCoords(:,2)),mean(recCoords(:,3)),label);
    set(t,'HorizontalAlignment','center')
end

% draw camera in 3D world
%  introduced in 2015a: cam = plotCamera('Location',T,'Orientation',R,'Opacity',0);
xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
axis equal tight
plotCamera('Location',T_est,'Orientation',R_est,'Opacity',0,'Size' ,10);
hold on
camc = R_est'*(-T_est);
plot3(camc(1),camc(2),camc(3),'o')

drawCoordframe(eye(3),20,'r','w') % world coordinate frame

drawCoordframeT(R_est',T_est,22,'m','c')% camera coordinate frame in the world frame

%plot3(-T(1),-T(2),-T(3),'o') %should be in center of camera

%% draw projected 3D coords onto the image
% to prove, display the hand-clicked coordinates and the 3D coordinates,
% converted to camera coordinates, and with the projective transformation
% applied.

%f1 = figure( 'units','normalized','outerposition',[0 0 1 1]);
f1 = figure();
set(f1,'name', 'Calibration check');
image(img,'Parent',gca);
set(gca,'fontsize', 18);
axis equal tight

a = axis;
% draw the principle point
line(oc,or,'Marker','*','color','white')
% ortho center line (?)
line( oc - fxe*[0,4]./[10,10],or - fye*[0,0]./[10,10], 'color','white','linewidth',2)
line( oc - fxe*[0,0]./[10,10], or - fye*[0,4]./[10,10],'color','white','linewidth',2)
axis(a)
for i = 1:size(DataXYZ,1)
    if mod(i,4)==0
        line(rcVal(i-[0,1,2,3,0],1),rcVal(i-[0,1,2,3,0],2),'color','red','linewidth',2)
        ptXYZw = DataXYZ(i-[0,1,2,3,0],:);
        % convert to camera coords
        ptXYZc=ptXYZw;
        for j = 1:5
            ptXYZc(j,:) = (R_est*(ptXYZw(j,:)'-T_est))';
        end
        % projection and covert to rc
        ptRC = [ oc - fxe*ptXYZc(:,2)./ptXYZc(:,3), or + fye*ptXYZc(:,1)./ptXYZc(:,3)];
        % draw the new line in a different color.  If everything is perfect, they will overlap.
        line(ptRC(4-[0,1,2,3,0],1),ptRC(4-[0,1,2,3,0],2),'color','blue','linewidth',2)
    end
end
% % (optional) Convert all points into camera coords:
% XYZinCam = (R_est*(DataXYZ'-T_est))';
% % Project points into 2D camera frame:
% RCreconPts =  [ oc - fxe*XYZinCam(:,2)./XYZinCam(:,3), or + fye*XYZinCam(:,1)./XYZinCam(:,3)];
% % Compare the standard deviation:
% std(rcVal -RCreconPts)

title({'Red are points clicked by user,', 'blue are reconstructed from 3D data and camera params'})
function drawCoordframe(R,s,color,superscript)
% draws coordinate frame rotated by R at the origin with arrow lengths s and
% color color.
z = zeros(3,1);
drawCoordframeT(R,z,s,color,superscript)

function drawCoordframeT(R,T,s,color,superscript)
% draws coordinate frame rotated by R at position T (in base frame with 
% arrow lengths s and color color).
if nargin < 4
    superscript = '';
end

drawArrow(R(:,1),T,s,color);
text(T(1) + 1.1*s*R(1,1),T(2) + 1.1*s*R(2,1),T(3) +1.1*s*R(3,1),['x^',superscript],'color',color);
drawArrow(R(:,2),T,s,color);
text(T(1) + 1.1*s*R(1,2),T(2) +1.1*s*R(2,2),T(3) +1.1*s*R(3,2),['y^',superscript],'color',color);
drawArrow(R(:,3),T,2*s,color);
text(T(1) + 2.3*s*R(1,3),T(2) +2.3*s*R(2,3),T(3) +2.3*s*R(3,3),['z^',superscript],'color',color);

function drawArrow(U,T,s,color)
%draws an arrow from position T in direction given by unit vector U that is
%s units long and color color.
l=line( T(1)+[0,s*U(1)],...
    T(2)+[0,s*U(2)],...
    T(3)+[0,s*U(3)]);
set(l,'color',color);
lm = line(T(1)+s*U(1),T(2)+s*U(2),T(3)+s*U(3), 'Marker','*');
set(lm,'color',color);

function  Rp = orthogonalize(Rest) % use SVD to ensure matrix is orthogonal
[Ur,Sr,Vr] = svd(Rest);
Rp = Ur*eye(size(Sr,1))*Vr';

