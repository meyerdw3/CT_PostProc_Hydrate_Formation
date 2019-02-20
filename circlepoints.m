%This function calculates a series of points around the edge of a circle
%with a particular center point and radius. Values are returned as integers
%for easy display on a pixel image.

%This script provided by Dr. David DiCarlo at the Petroleum and Geosystems
%Engineering Department at the University of Texas at Austin and was used
%with permission.

%Inputs:
%   xCenter, yCenter = x, y coordinates of circle center
%   radius = circle radius

function pts=circlepoints(xCenter,yCenter,radius)

%Check to make sure all the inputs are provided
if(nargin<3)
   error('Too few arguements');
end

%Check to make sure the inputs are integers
if(rem(xCenter,1)~=0 | rem(yCenter,1)~=0 | rem(radius,1)~=0)
   warning('Increments are by whole numbers and using non-integers might not produce desired results');
end

x = 0;
y = radius;
p = 1 - radius;

pts=[];
pt=GetPoints(xCenter, yCenter, x, y);

%Calculate points on the circle at different distances away from the center
pts=[pts pt];
while(x<y)
   x=x+1;
   if(p < 0)
      p =p + (2*x + 1);
   else
      y=y-1;
      p = p+ (2*(x-y) + 1);
   end
   pt=GetPoints(xCenter, yCenter, x, y);
   
   %Append new points to total points list
   pts=[pts;pt];
end

%Sort the list of points
pts=sortrows(pts);

%Force the following while loop to run once
prevsz=length(pts)+1;

%Sort through the list and remove duplicate values
while(length(pts)~=prevsz)
   prevsz=length(pts);
   n=1;
   while(n<length(pts)),
      if(pts(n,:)==pts(n+1,:))
         pts(n,:)=[];
      end
      n=n+1;
   end
end

%This function produces a series of coordinates at various x, y offsets
%from the circle center
function pt=GetPoints(xCenter,yCenter,x,y)
pt(1,:)=[xCenter + x, yCenter + y];
pt(2,:)=[xCenter - x, yCenter + y];
pt(3,:)=[xCenter + x, yCenter - y];
pt(4,:)=[xCenter - x, yCenter - y];
pt(5,:)=[xCenter + y, yCenter + x];
pt(6,:)=[xCenter - y, yCenter + x];
pt(7,:)=[xCenter + y, yCenter - x];
pt(8,:)=[xCenter - y, yCenter - x];   
