%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   User clicks on points in the camera image that correspond to 3D 
%   coordinates. The pixel coordinates and 3D coordinates are returned 
%   as arrays.
% 
%   Instruction video: https://youtu.be/WEYwitb6dTo
%   Part of https://www.mathworks.com/matlabcentral/fileexchange/73079
%   "calibrate-camera-with-one-photo-linear-method"
%
%     For the image  'imageFilename', returns
%     rcVal:    (row,column) entries for each XYZ coordinate in DataXYZ
%     DataXYZ:  (X,Y,Z) entries for each (r,c) coordinate in rcVal
%     DataW:  the user must generate this before running the program by
%     measuring the 3D coordinates of multiple points in the image.  Each
%     line in DataW is x, y, z in world coordinates, a text description.
%     i.e.:  0,0,-.001, {'table, bottom left'};
%
%    Fig. 1 is the image, where the user must click on pixel coordinates.
%    Fig. 2 shows a 3D view of the world (rotate the view manually if
%    needed).  The 3D view is populated with patches and titles
%    corresponding to the names in DataW (must be entered by the user).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By  Aaron T. Becker, 10/25/2019, 
% Updated by Shreyas Poyrekar 10/25/2019 for mouse clicks.
% Updated by Utkarsh Gupta on 01/12/2024 improving the code readability
% by using arrays to store plots and added functionality 
% that lets user 'undo' their clicks.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rcVal,DataXYZ, DataW] = gather3Ddata(imageFilename)

if nargin<1
    imageFilename = 'CameraCalP1.jpg';
end

close all

DataW = [... % x, y, z in world coordinates, a text description. All distances are in cm. Designed to work with 4-sided shapes
    % TODO 1: for your picture, you need to generate these yourself
    0,0,-.001, {'table, bottom left'};
    100,0,-.001, {'table, bottom right'};
    100,240,-.001, {'table, top right'};
    0,240,-.001, {'table, top left'};

    0,0,0, {'paper 1, bottom left'};
    28,0,0, {'paper 1, bottom right'};
    28,22,0, {'paper 1, top right'};
    0,22,0, {'paper 1, top left'};

    51,37,2.5, {'book top, bottom left'};
    75,37,2.5, {'book top, bottom right'};
    74,56,2.5, {'book top, top right'};
    50,56,2.5, {'book top, top left'};
 
    7,58,0, {'paper 2, bottom left'};
    29,58,0, {'paper 2 bottom right'};
    29,86,0, {'paper 2, top right'};
    7,86,0, {'paper 2, top left'};

    78.5,56,0, {'paper 3, bottom left'};
    100,56,0, {'paper 3, bottom right'};
    100,84,0, {'paper 3, top right'};
    78.5,84,0, {'paper 3, top left'};
   
    78.5,108.5,0, {'paper 4, bottom left'};
    100,108.5,0, {'paper 4, bottom right'};
    100,136,0, {'paper 4, top right'};
    78.5,136,0, {'paper 4, top left'};
   
    32,83,30, {'cube top, bottom left'};
    63,83,30, {'cube top, bottom right'};
    63,113,30, {'cube top, top right'};
    32,113,30, {'cube top, top left'};
   
    32,83,0, {'cube front, bottom left'};
    63,83,0, {'cube front, bottom right'};
    63,83,30, {'cube front, top right'};
    32,83,30, {'cube front, top left'};
   
    63,83,0, {'cube right, bottom left'};
    63,113,0, {'cube right, bottom right'};
    63,113,30, {'cube right, top right'};
    63,83,30, {'cube right, top left'};
   
    74,150,0, {'tiny cube front, bottom left'};
    80.5,150,0, {'tiny cube front, bottom right'};
    80.5,150,6.5, {'tinycube front, top right'};
    74,150,6.5, {'tinycube front, top left'};
   
    13,190,0, {'pancake box front, bottom left'};
    30,190,0, {'pancake box front, bottom right'};
    30,190,26, {'pancake box front, top right'};
    13,190,26, {'pancake box front, top left'};
   
    56,192.5,0, {'tea box front, bottom left'};
    76,192.5,0, {'tea box front, bottom right'};
    76,192.5,29, {'tea box front, top right'};
    56,192.5,29, {'tea box front, top left'};
   
    -30,55,0, {'cake box front, bottom left'};
    -30,69,0, {'cake box front, bottom right'};
    -30,69,18, {'cake box front, top right'};
   -30,55,18, {'cake box front, top left'};
    ];

rcVal = [];

function setTitle(click_count)
        title(DataW(click_count+1,4)+" (click "+(click_count+1)+"/"+size(DataW,1)+")");
end

% draw 3D world
m = "Figure 2. Use the Rotate 3D button to move around in 3D Space";
S.fH1 = figure('Name',m,'NumberTitle','off'); clf;
set(gca,'fontsize', 18);
for j = 1:size(DataW,1)/4
    recCoords = cell2mat(DataW(4*(j-1)+(1:4),1:3));
    label = DataW(4*(j-1)+1,4); 
    label = strsplit(cell2mat(label),',');
    label = label(1);
    p=patch(recCoords(:,1),recCoords(:,2),recCoords(:,3),rand(1,1));
    set(p,'FaceColor','flat','FaceAlpha',0.5);
    t=text( mean(recCoords(:,1)),mean(recCoords(:,2)),mean(recCoords(:,3)),label);
    set(t,'HorizontalAlignment','center')
end
axis equal
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');

lineHandles = [];
redlineHandles =  [];
patchHandles = [];

click_count = 0;
polygon_flag = 0;

lineHandles = gobjects(size(DataW, 1), 1);  % initializing set number of objects based on no. of elements in DataW 
redlineHandles = gobjects(size(DataW, 1)/4, 1);
patchHandles = gobjects(size(DataW, 1)/4, 1);

color = 'g';

figTitle = "Figure 1. Gather Data: click on the points requested, press 'u' to undo.";
S.fH2 = figure( 'units','normalized','outerposition',[0 0 1 1],'MenuBar', 'none', 'ToolBar', 'none', 'Name',figTitle,'NumberTitle','off');
img = imread(imageFilename);
S.aH = axes;
S.iH = image(img,'Parent',gca); hold on
set(gca,'fontsize', 17);
setTitle(click_count);
size(img);

draglineHandles = plot(NaN,NaN, color, 'LineWidth', 2, 'LineStyle', '--');  % temporarily stores the lines that help guide a click

for k = 1:size(DataW,1)/4
    patchHandles(k) = patch(NaN, NaN, 'b', 'FaceAlpha', 0.2);
end

set(S.fH2, 'pointer','crosshair');
set(S.fH2, 'WindowButtonUpFcn', @startDragFcn);
set(S.fH2, 'WindowButtonMotionFcn', @draggingFcn);
set(S.fH2, 'KeyPressFcn', @keyPressedFcn);

while (click_count <= size(DataW,1))
    waitforbuttonpress
    if click_count + 1 == size(DataW,1)
        startDragFcn()  % save the last click
        break;  % Exit the loop if enough points have been collected
    end
end
close(S.fH2);  % clean up by closing both windows 
close(S.fH1);


    function checkPolygon() % checks if a polygon has been formed or not
        if mod(click_count,4) == 0
            polygon_flag = 1;
            set(patchHandles(ceil(click_count/4)), 'XData', [], 'YData', []);
        else 
            polygon_flag = 0;
        end
    end

    function plotLines()  % draws outlines around completed shapes and the shape in progress
        if polygon_flag == 1
            % bounds the polygon formed with red line.
            redlineHandles(click_count/4) = plot([rcVal(click_count-3:click_count,1)',rcVal(click_count-3,1)], [rcVal(click_count-3:click_count,2)',rcVal(click_count-3,2)], 'r', 'LineWidth', 2,'LineStyle','-');
        elseif mod(click_count,4) ~= 1
            % doesn't connect points 4th and 5th...
            lineHandles(click_count) = plot([rcVal(click_count-1:click_count, 1)], [rcVal(click_count-1:click_count, 2)], color, 'LineWidth', 2, 'LineStyle', '--');
        end
    end
    
    function plotPatches(x,y)  % plots a transparent 4-sided polygon every 4th click to aid in positioning
        if polygon_flag == 0  
            if mod(click_count,4) ~= 0
                % only plots guiding lines until the 3rd point of the
                % polygon (avoids connecting the 4th and 1st point of next
                % polygon)
                set(draglineHandles, 'Xdata', [rcVal(end,1)',x],'Ydata',[rcVal(end,2)', y])
            end
            if mod(click_count+1,4) == 0
                % creates a patch only after the 3rd, 7th, 11th,... point is clicked
                set(patchHandles((click_count+1)/4), 'Xdata', [rcVal(click_count-2:end,1)',x],'Ydata',[rcVal(click_count-2:end,2)' y])
            end
        end
    end
    
    function startDragFcn(varargin)
        % reacts to each mouse click, drawing a line connecting two clicked points 
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);

        click_count = click_count + 1;
        rcVal(click_count,:)=[x,y];
        
        if click_count < size(DataW,1)
            setTitle(click_count);
        end

        checkPolygon();
        plotLines();          
    end

    function draggingFcn(varargin)
        % reacts to mouse movement drawing a line to guide the user
        % also guides in clicking the 4th point of the polygon by creating a
        % patch
        pt = get(S.aH, 'CurrentPoint');
        x = pt(1,1);
        y = pt(1,2);
        if  click_count > 0
            plotPatches(x,y);
        end
    end
    
    function keyPressedFcn(~, event)
        keyPressed = event.Key;
        % Check if the pressed key is 'u'
        if strcmpi(keyPressed, 'u') && click_count>0
            % Actions for undo 
            rcVal = rcVal(1:end-1,:);

            delete(lineHandles(click_count));
            delete(draglineHandles);
            set(patchHandles(ceil(click_count/4)), 'XData', [], 'YData', []);
            
            checkPolygon();

            if polygon_flag == 1
                delete(redlineHandles(floor(click_count/4)));
                polygon_flag = 0;
            end

            lineHandles(click_count) = plot(NaN,NaN, color, 'LineWidth', 2, 'LineStyle', '--');
            draglineHandles = plot(NaN,NaN, color, 'LineWidth', 2, 'LineStyle', '--');
            click_count = click_count-1;
            setTitle(click_count);
        end
    end

DataXYZ = cell2mat(DataW(:,1:3));
save([imageFilename(1:end-4),'Data'], 'rcVal','DataXYZ', 'DataW','imageFilename')

end