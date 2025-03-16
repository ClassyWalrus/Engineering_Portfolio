%A.2

alpha = 35;
beta = 23;

Left = cosd(alpha)*cosd(beta);
Right = (1/2)*(cosd(alpha-beta)+cosd(alpha+beta));

if Left == Right
   output = true;
else
    output = false;
end

output

%If output is equal to 1, then each side of the equation is equivalent. If
%output is equal to 0, then each side is not equivalent.


%A.5

RowVector = [14: -3: -10];
ColumVector = RowVector.'


%A.8

A = [2.5:1:7.5;
    42:-3.4:25; 
    15:-.4:13; 
    3:-1:-2]

va = A(2,3:1:6)

vb = A(2:1:4,5)

%A.10

x = [-1.2:0.8:3.6];
y = ((2.*x.^2 - 16.*x + 4).^2)./(x+15)

plot(x,y,"Color","k","Marker","*")
xlabel('X Values')
ylabel('Y Values')


%A.24

%TrueVolume = [];
for height = 0:.1:2.8
    
 Volume = Vfuel(height);
 %TrueVolume = [TrueVolume; {Volume}]
end

plot(Volume, height)

function [V] = Vfuel(h)
r = .60;
  if(h>=0) && (h<=r) && (h<=2)
          V = (pi*h^2)/3*(3*r-h);

  elseif (h>=0) &&(h>=r) && (h<=2)
          V = 2*pi*r^3+pi*r^2*h;
  
  elseif (h>=0) && (h>=2) && (h<=2.8)
          V =  2*pi*r^3+pi*r^2*1.8+pi*r^2*(h/3);
  else
  end

end

