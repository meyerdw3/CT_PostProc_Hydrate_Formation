%This function takes in the three points selected at the edge of the sample
%in the CT slices and produces a circle that intersects all three points,
%defined by the radius and center coordinates. 

%This script written by Dr. David DiCarlo at the Petroleum and Geosystems
%Engineering Department at the University of Texas at Austin and was used
%with permission.

%Inputs
%   P1, P2, P3 = x, y coordinates of three points along the edge of the
%                sample

%Outputs:
%   x_center, y_center = x, y coordinates of the center of the circle.
%   radius = radius of the circle.

function [x_center y_center radius]=circle_eq (P1, P2, P3)
%Calculate slope of the perpendicular bisector between P1 and P2
m1=-1*(P1(1)-P2(1))/(P1(2)-P2(2));

%Calculate slope of the perpendicular bisector between P2 and P3
m2=-1*(P3(1)-P2(1))/(P3(2)-P2(2));

%Determine the circle center as the intersection point of the perpendicular 
%bisectors
x_center=(m1*(P1(1)+P2(1))/2-m2*(P2(1)+P3(1))/2+(P3(2)-P1(2))/2)/(m1-m2);
y_center=(P1(2)+P2(2))/2+m1*(x_center-(P1(1)+P2(1))/2);

%Calculate the radius of the circle
radius=sqrt((P1(1)-x_center)^2+(P1(2)-y_center)^2);
end


