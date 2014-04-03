% Create_signals.m
%
% Creates an array of twiddle factors and saves it as header file (.h) for
% use with an FFT function in C
%
%
%
clc
close all 
clear all
fs = 8000;
s_num = 2;
fft_size = 1024; %  <--------- Change the length of the FFT HERE
len = fft_size*2;% Actual length of twiddles vector/array

delta = 2*pi/len;
for n=1:2:len	    % set up twiddle constants in w 		
  		w(n+1) = sin(delta*n);
  		w(n) = cos(delta*n);
end


%%%% Generate C code
coefflen=length(w);
fname = input('enter filename for coefficients ','s');
fid = fopen(fname,'wt');
fprintf(fid,'#ifndef SIGS_\n');
fprintf(fid,'#define SIGS_\n');
fprintf(fid,'// %s\n',fname);
fprintf(fid,'// this file was generated automatically using function create_signals.m\n',fname);
fprintf(fid,'// Created on %s\n',datestr(now));
fprintf(fid,'\n#define LEN %d\n',coefflen);
%fprintf(fid,'\n#define SNUM %d\n',s_num);
fprintf(fid,'\nfloat w[LEN] = {');
% j is used to count coefficients written to current line
% in output file
j=0;
% i is used to count through coefficients
for i=1:coefflen  
% if six coeffs have been written to current line
% then start new line
  if j>5    
    j=0; 
    fprintf(fid,'\n');  
  end  
% if this is the last coefficient then simply write
% its value to the current line
% else write coefficient value, followed by comma
  if i==coefflen
   fprintf(fid,'%2.4E',w(i));
  else
    fprintf(fid,'%2.4E,',w(i));  
    j=j+1;
  end
end

fprintf(fid,'\n};\n');
fprintf(fid,'\n\n#endif /*SIGS_*/');
fclose(fid);   
%----------- End of file part



% MANY * 2 format

return
fprintf(fid,'\nfloat w[LEN] = { {');
% j is used to count coefficients written to current line
% in output file
j=0;
% i is used to count through coefficients
for i=1:coefflen  
% if six coeffs have been written to current line
% then start new line
  if j>rowlength    
    j=0; 
    fprintf(fid,'},\n\t\t\t\t{');  
  end  
% if this is the last coefficient then simply write
% its value to the current line
% else write coefficient value, followed by comma
  %if i==coefflen
  % fprintf(fid,'%2.4E',w(i));
  if j==rowlength
      fprintf(fid,'%2.4E',w(i));  
    j=j+1;    
  else
    fprintf(fid,'%2.4E,',w(i));  
    j=j+1;
  end
end

fprintf(fid,'}\n};\n');
fprintf(fid,'\n\n#endif /*SIGS_*/');
fclose(fid);   