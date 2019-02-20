%This function produces a square matrix of zeros with a circle of ones 
%circumscribed within the matrix. This is used in the main script as an
%initial start point of data masks. This script only accepted integer
%inputs.

%This script provided by Dr. David DiCarlo at the Petroleum and Geosystems
%Engineering Department at the University of Texas at Austin and was used
%with permission.

%Input:
%   n = diameter of the circle/width and height of the complete matrix

%Output:
%   y = matrix of dimension (n x n) with a circle circumscribed.

function y = imcircle(n)

%Check if the input is an integer
if rem(n,1) > 0, 
   disp(sprintf('n is not an integer and has been rounded to %1.0f',round(n)))
   n = round(n);
end

%Check if the input is less than 1
if n < 1     % invalid n
   error('n must be at least 1')

%If the input is 1 - 3, output a square matrix entirely filled with ones
elseif n < 4 % trivial n
   y = ones(n);

%If the input is an even number
elseif rem(n,2) == 0,
   
   DIAMETER = n;
   diameter = n-1;
   RADIUS = DIAMETER/2;
   radius = diameter/2;
   height_45 = round(radius/sqrt(2));
   width = zeros(1,RADIUS);
   semicircle = zeros(DIAMETER,RADIUS);   
   
   for i  = 1 : height_45
       upward = i - 0.5;
       sine = upward/radius;
       cosine = sqrt(1-sine^2);
       width(i) = ceil(cosine * radius);
   end

   array = width(1:height_45)-height_45;

   for j = max(array):-1:min(array)
       width(height_45 + j) = max(find(array == j));
   end

   if min(width) == 0
      index = find(width == 0);
      width(index) = round(mean([width(index-1) width(index+1)]));
   end

   width = [fliplr(width) width];

   for k  = 1 : DIAMETER
       semicircle(k,1:width(k)) = ones(1,width(k));
   end   

   y = [fliplr(semicircle) semicircle];

%If the input is an odd number
else
   
   DIAMETER = n;
   diameter = n-1;
   RADIUS = DIAMETER/2;
   radius = diameter/2;
   semicircle = zeros(DIAMETER,radius);
   height_45 = round(radius/sqrt(2) - 0.5);
   width = zeros(1,radius);

   for i  = 1 : height_45
       upward = i;
       sine = upward/radius;
       cosine = sqrt(1-sine^2);
       width(i) = ceil(cosine * radius - 0.5);
   end

   array = width(1:height_45) - height_45;

   for j = max(array):-1:min(array)
       width(height_45 + j) = max(find(array == j));
   end

   if min(width) == 0
      index = find(width == 0);
      width(index) = round(mean([width(index-1) width(index+1)]));
   end

   width = [fliplr(width) max(width) width];

   for k  = 1 : DIAMETER
       semicircle(k,1:width(k)) = ones(1,width(k));
   end   

   y = [fliplr(semicircle) ones(DIAMETER,1) semicircle];

end

