function orthoc = findOrthocenter( imageFilename)
% Finds the principle point of the camera by first finding a triangle whose
% vertices are defined by the vanishing point of 3 mutually
% othogonal sets of parallel lines.  The orthocenter of this triangle is
% the principle point of the camera.
%
%  Instruction video: https://youtu.be/WEYwitb6dTo
%
% by  Aaron T. Becker, 10/30/2018
% Updated by Nima Eskandari in 2018 nima.eskandari@outlook.com such
% that draws a line from the first clicked point to the mouse position 
% as you find the three sets of orthogonal lines.
if nargin<1
    imageFilename = 'CameraCalP1.jpg';%    1.0e+03 * [1.3198; 0.7816];
end
close all

%imshow(img);
f1 = figure( 'units','normalized','outerposition',[0 0 1 1]);
set(f1,'name', 'find orthocenter','MenuBar', 'none', 'ToolBar', 'none');
img = imread(imageFilename);
image(img,'Parent',gca);
set(gca,'fontsize', 18);
hold on
% % zoom in (optional
% title('zoom in: click twice to set the zone to be zoomed into')
% bb =ginput(2);
% axis([bb(1,1),bb(2,1),bb(1,2),bb(2,2)]);

%pause()
intxy1 = drawParallelLines(1);
intxy2 = drawParallelLines(2);
intxy3 = drawParallelLines(3);

%draw triangle
l1 = plot([intxy1(1),intxy2(1),intxy3(1),intxy1(1)],[intxy1(2),intxy2(2),intxy3(2),intxy1(2)],'.-');
set(l1,'color','k','linewidth',2);

% calculate the orthocenter (??)
orthoc = orthocenterTri(intxy1,intxy2,intxy3);
plot(orthoc(1),orthoc(2),'w*')


axis([min([0,intxy1(1),intxy2(1),intxy3(1)]),...
    max([size(img,1),intxy1(1),intxy2(1),intxy3(1)])...
    min([0,intxy1(2),intxy2(2),intxy3(2)]),...
    max([size(img,2),intxy1(2),intxy2(2),intxy3(2)])]);

%Save the ortho center!
save([imageFilename(1:end-4),'Ortho'], 'orthoc','imageFilename')
axis equal


    function orthoc = orthocenterTri(t1,t2,t3)
        %finds the orthocenter of a triangle with vertices at t1,t2,t3
        %slope of side 2-3:
        slope23 = [t2(2)-t3(2);t3(1)-t2(1)];
        slope13 = [t1(2)-t3(2);t3(1)-t1(1)];
        
        orthoc = LineIntersectionPoint(t1,t1+slope23,t2,t2+slope13);
        plot(t1(1)+[0,slope23(1)],t1(2)+[0,slope23(2)],'-')
        plot(t2(1)+[0,slope13(1)],t2(2)+[0,slope13(2)],'-')
        
    end

    function intxy = LineIntersectionPoint(a1,a2,b1,b2)
        %given two lines, each defined by two points a1 to a2 and b1 to b2,
        % computes the intersection point intxy or returns NaN
        denom = det([a1-a2,b1-b2]);
        if denom ~= 0
            intxy = ( det([a1, a2])*(b1-b2) - det([b1, b2])*(a1-a2)) / denom;
        else
            intxy = [NaN,NaN];
        end
    end

    function intxy = drawParallelLines(num)
        
        colors = ['r','b','g'];
        color = colors(num);
        strs = { 'step 1: find parallel lines';
            'step 2: find lines perpendicular to red lines';'step 3: find lines perpendicular to red and blue lines'}; 
        % query user to select two lines that are parallel in 3D world but have
        % a vanishing point in the image
        title({strs{num};'select two points on line 1'})
        pts = interactiveBoundingBoxNima(2, f1, [0, size(img,2)], [0, size(img,1)]);
        a1= pts(:,1);a2= pts(:,2);
        l1 = plot([a1(1),a2(1)],[a1(2),a2(2)],'.-');
        set(l1,'color',color,'linewidth',2);
        
        title({strs{num};'select two points on line 2 parallel to line 1'})
        pts = interactiveBoundingBoxNima(2, f1, [0, size(img,2)], [0, size(img,1)]);
        b1= pts(:,1);b2= pts(:,2);
        l2=plot([b1(1),b2(1)],[b1(2),b2(2)],'.-');
        set(l2,'color',color,'linewidth',2);
        
        intxy = LineIntersectionPoint(a1,a2,b1,b2);
        
        l3=plot([b1(1),intxy(1) a1(1)],[b1(2),intxy(2), a1(2)],'.-');
        set(l3,'color',color);
        
    end

end


function points = interactiveBoundingBoxNima(numberOfPolyonSides, figureToDrawOn, boundsX, boundsY)
% This file demonstrates the interactivty I want for 'gather3Ddata.m'
%
% Aaron Becker
S.numberOfPolygonSides = '';
S.fH = '';
S.aH = '';
S.maxX = 0;
S.minX = 0;
S.maxY = 0;
S.minY = 0;
S.title = '';

if nargin < 1
    S.numberOfPolygonSides = 2;
    S.fH = figure('menubar','none');
    im = imread( 'CameraCalP1.jpg' );
    S.aH = axes;
    S.iH = imshow( im ); hold on
    [height, width] = size(im);
    S.maxX = width;
    S.minX = 0;
    S.maxY = height;
    S.minY = 0;
else
    S.numberOfPolygonSides = numberOfPolyonSides;
    S.fH = figureToDrawOn;
    S.aH = gca;
    S.maxX = boundsX(2);
    S.minX = boundsX(1);
    S.maxY = boundsY(2);
    S.minY = boundsY(1); 
    S.title = S.aH.Title.String;
end

axis on


S.done = 0;
X = [];
Y = [];
S.draw = 0;
lineHandle = [];
patchHandle = [];
color = 'white';

%TODO: I'm unconvinced if cross hairs makes it easier or harder to click on
%points (Aaron Becker)
crosshair(1) = plot([0,0], [0,0], 'w', 'LineWidth', 1);
crosshair(2) = plot([0,0], [0,0], 'w', 'LineWidth', 1);
crosshair(3) = plot([0,0], [0,0], 'k', 'LineWidth', 0.5);
crosshair(4) = plot([0,0], [0,0], 'k', 'LineWidth', 0.5);

title([S.title ; '(Point ', num2str(1), '/' , num2str(S.numberOfPolygonSides), ')'])
set(S.fH, 'WindowButtonUpFcn', @startDragFcn);
set(S.fH, 'WindowButtonMotionFcn', @draggingFcn );
while S.done == 0
    waitforbuttonpress;
    pause(0.5)
end

    function startDragFcn(varargin)
        % reacts to each mouse click
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);
        [x,y] = check_axis(x, y);

        S.draw = S.draw+1;
        if S.draw == 1
            X = x;
            Y = y;
            lineHandle = plot([X,X], [Y,Y], color, 'LineWidth', 2);
        elseif S.draw < S.numberOfPolygonSides 
            X = [X,x];
            Y = [Y,y];
        else
            set(lineHandle, 'Xdata', [X,x,X(1)],'Ydata',[Y y Y(1)])
            delete(patchHandle)
            S.done = 1;
            set(S.fH, 'WindowButtonUpFcn','')
            set(S.fH, 'WindowButtonMotionFcn','')
            title(S.title);
            points = [[X,x];[Y,y]];
        end

        if S.draw == S.numberOfPolygonSides - 1
            patchHandle = patch(X,Y,'b','FaceAlpha',0.2);
        end
        if  S.draw ~= 0 && S.draw < S.numberOfPolygonSides
            title([S.title ; '(Point ', num2str(S.draw + 1), '/' , num2str(S.numberOfPolygonSides), ')'])
        end
    end

    function draggingFcn(varargin)
        % reacts to mouse movement
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);
        [x,y] = check_axis(x, y);
        
        if  S.draw > 0
            set(lineHandle, 'Xdata', [X,x],'Ydata',[Y y])
            if S.draw == S.numberOfPolygonSides - 1
                set(patchHandle, 'Xdata', [X,x],'Ydata',[Y y])
            end
        end
        if  S.draw < S.numberOfPolygonSides - 1
            %TODO: I'm unconvinced if cross hairs makes it easier or harder to click on points (Aaron Becker)
            a = axis;
            gap = 6;
            set(crosshair(1),'Xdata',[a(1),x-gap,NaN,x+gap,a(2)],'Ydata',y*[1,1,1,1,1]);
            set(crosshair(2),'Xdata',x*[1,1,1,1,1],'Ydata',[a(3),y-gap,NaN,y+gap,a(4)]);
            set(crosshair(3),'Xdata',[a(1),x-gap,NaN,x+gap,a(2)],'Ydata',y*[1,1,1,1,1]);
            set(crosshair(4),'Xdata',x*[1,1,1,1,1],'Ydata',[a(3),y-gap,NaN,y+gap,a(4)]);
        else
            set(crosshair(1),'Xdata',[0,0],'Ydata',[0,0]);
            set(crosshair(2),'Xdata',[0,0],'Ydata',[0,0]);
            set(crosshair(3),'Xdata',[0,0],'Ydata',[0,0]);
            set(crosshair(4),'Xdata',[0,0],'Ydata',[0,0]);
        end
        drawnow 
    end

    function [x,y] = check_axis(x, y)
        if (y > S.maxY)
            y = S.maxY;
        end
        if y < S.minY
            y = S.minY;
        end
        if x > S.maxX
            x = S.maxX;
        end
        if x < S.minX
            x = S.minX;
        end
    end

end
