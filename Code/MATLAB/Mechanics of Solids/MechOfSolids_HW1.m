clc;
clear;
[x,y1,y2] = createData();
h = animatedline;

% Iterate over the vectors and use addpoints, drawnow to animate line
for k = 1:length(x)
    addpoints(h,x(k),y1(k)); % Add points to the animated line 
    drawnow limitrate        % Use limitrate to speed up animation
end

function [x, y1, y2] = createData()
    x = linspace(0,4*pi,100000);
    y1 = sin(x);
    y2 = cos(x);
end