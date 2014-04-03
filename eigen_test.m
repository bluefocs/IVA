%
% Code to test if the method behind the eigenvalue/eigenvector C code 
% actually works. 
%
%
clear all
close all
clc

 X = [-1+1j,7+1.5j;2-.5j,.7+3.2j];
%X = [-1, 7; 2, 2.1];

[v,d] = eig(X)

D = (real(X(1,1))*real(X(2,2))) - (imag(X(1,1))*imag(X(2,2))) - (real(X(1,2))*real(X(2,1))) + (imag(X(1,2)) * imag(X(2,1)));
D = D + j*((real(X(1,1))*imag(X(2,2))) + (imag(X(1,1))*real(X(2,2))) - (real(X(2,1))*imag(X(1,2))) - (imag(X(2,1)) * real(X(1,2))));%// Determinant (real part)

%	T.real = X[0].real + X[3].real; %//Trace real part
%	T.imag = X[0].imag + X[3].imag;
% Trace
T = trace(X);

% http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
% for the eigenvalues to be found need to express the bit under the square
% root as (stuff) + j(stuff) so that it can be rewritten sqrt(r)(cos(0.5*arg) + sin(0.5*arg))



realpart = (real(T)^2 - imag(T)^2)/4 - real(D);% Real part under the square root
imagpart = (0.5 * imag(T) * real(T)) - imag(D); % imaginary part under the square root
%arg = tan((imagpart) / (realpart));% Shouldn't this be arctan???
arg = atan((imagpart) / (realpart));

r = sqrt(realpart^2 + imagpart^2);


eigvals(1) = (real(T)/2) + sqrt(r)*cos(0.5*arg) + (1i)*((imag(T)/2) + sqrt(r)*sin(0.5*arg));
eigvals(2) = (real(T)/2) - sqrt(r)*cos(0.5*arg) + (1i)*((imag(T)/2) - sqrt(r)*sin(0.5*arg));


if((real(X(1,2))==0) && (real(X(2,1))==0) && (imag(X(1,2))==0) && (imag(X(2,1))==0))
    
    eigvecs = eye(2);
    
elseif((real(X(2,1))~=0) && (imag(X(2,1))~=0))% if c is not zero
    
    eigvecs(1,1) = real(eigvals(1)) - real(X(2,2)) + (1i)*(imag(eigvals(1)) - imag(X(2,2)));
    eigvecs(2,1) = X(2,1);%c
    eigvecs(1,2) = real(eigvals(2)) - real(X(2,2)) + (1i)*(imag(eigvals(2)) - imag(X(2,2)));
    eigvecs(2,2) = X(2,1);%c
    
else %// if b!=0
    eigvecs(1,1) = X(1,2);%b
    eigvecs(2,1) = real(eigvals(1)) - real(X(1,1)) + (1i)*(imag(eigvals(1)) - imag(X(1,1)));
    eigvecs(1,2) = X(1,2);%b
    eigvecs(2,2) = real(eigvals(2)) - real(X(1,1)) + (1i)*(imag(eigvals(2)) - imag(X(1,1)));

end


%//Normalise
% norm_fact = real(eigvecs(1,1))^4 + imag(eigvecs(1,1))^4 + real(eigvecs(2,1))^4 + imag(eigvecs(2,1))^4;  %//normalisation factor
% norm_fact = norm_fact^0.25; %// normalisation factor
% eigvecs(:,1) = eigvecs(:,1) / norm_fact;
eigvecs(:,1) = eigvecs(:,1) / sqrt(eigvecs(1,1)^2 + eigvecs(2,1)^2);

% norm_fact = real(eigvecs(1,2))^4 + imag(eigvecs(1,2))^4 + real(eigvecs(2,2))^4 + imag(eigvecs(2,2))^4;%//normalisation factor
% norm_fact = norm_fact^0.25;
% eigvecs(:,2) = eigvecs(:,2) / norm_fact;
eigvecs(:,2) = eigvecs(:,2) / sqrt(eigvecs(1,2)^2 + eigvecs(2,2)^2);


eigvecs
eigvals

disp(num2str(det(X)))
D