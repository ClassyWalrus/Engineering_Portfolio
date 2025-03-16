%%Problem 3.5
disp('3.5 Solution:')
f = @(x) x^4-200;

imax = 5;
X_Newton = zeros(1,imax);
X_Newton(1) = (3*8^4+200)/(4*8^3);

for i = 2:imax
 X_Newton(i) = (3*X_Newton(i-1)^4+200)/(4*X_Newton(i-1)^3);
end


fprintf("3-5_X: %f\n", X_Newton);


%%Problem 3.7
disp('3.7 Solution:');
Xa = -4;
Xb = -5;
imax = 4;


Fun = @(x) 1.2*x^3+2*x^2-20*x-10;

for i = 1:imax
 FUNXb = Fun(Xb);
 Xi = Xb-FUNXb*(Xa-Xb)/(Fun(Xa)-FUNXb);
 Xs = Xi;
 Xa = Xb;
 Xb = Xi;
 fprintf("3-7_X: %f\n", Xs);
end


%3.5 with FZERO
disp('3.5 FZERO Solution:')

fnaut = 3;
f = @(x) x^4-200;

fz = fzero(f,fnaut);
disp(fz);



%3.7 with FZERO
disp('3.7 FZERO Solution:')

fnaut = -4;
f = @(x) 1.2*x^3+2*x^2-20*x-10;

fz = fzero(f,fnaut);
disp(fz);



%%Problem 3.24
disp('3.24 Solution: ')

Fun = @(x) x-2*exp(-x);

Xest = 20;
Xs = SteffensenRoot(Fun,Xest);
disp(Xs);

function [Xs] = SteffensenRoot(Fun,Xest)
Xi = Xest;

    for i = 1:100
     Fxi = Fun(Xi);
     Xii = Xi - ((Fxi^2)/(Fun(Xi+Fxi)-Fxi));
     Err = abs((Xii-Xi)/Xi);
     Xi = Xii;
     Xs = Xi;

     if Err < 10^-6
      break
     end

     if i == 100 && Err>10^-6
      disp('Answer not found in 100 iterations');
     break
     end

    end
end