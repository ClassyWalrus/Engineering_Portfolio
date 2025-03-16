%% PartA

h = 0.7;  % step size
x = 0:h:2.1;  % the range of x
y = zeros(size(x));  % allocate the result y
y(1) = 2;  % the initial y value
n = numel(y);  % the number of y values
% The loop to solve the DE
for i=1:n-1
    derivative = (x(i)^2)/y(i)
    y(i+1) = y(i)+h*derivative;
end

%% Part B

