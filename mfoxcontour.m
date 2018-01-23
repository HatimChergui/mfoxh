function c = mfoxcontour(W, dim, index, A, varargin)
% Copyright 2018 Hatim Chergui, Mustapha Benjillali, 
% and Mohamed-Slim Alouini. Contact email <chergui@ieee.org>
% 
% ****************************************************************************************************************
% Citation: If you use this software or any (modified) part of it, please cite it as:
% Hatim Chergui, Mustapha Benjillali and Mohamed-Slim Alouini. (2018, January 22). 
% Multivariate Fox H-Function C/MEX Package: mfoxh (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.1157194
% ****************************************************************************************************************
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%  Help:
%  W  : control the width of the integration interval in [-i\infty +i\infty]
%  dim: stands for the diemension
%  
%  Example (2 dimensions):
%  index = [0 1 1 1 1 1]; % [0 m m1 n1 m2 n2 ...mM nM]
%  z = [1 2]; % [z1...zM]
%  A = [1.5 1.0 1.0; 2.0 1.0 1.0]; % [a1 alpha_1,1...alpha_1,M; ...; ap alpha_p,1...alpha_p,M]
%  B = [2.0 1.0 1.0]; % [b1 beta_1,1...beta_1,M;...; bq beta_q,1 ...beta_q,M]
%  A1 = [-1 1]; % [a1_1 alpha1_1;...;a1_p1 alpha1_p1]
%  B1 = [0 1]; % [b1_1 beta1_1;...;b1_q1 beta1_q1]
%  A2 = [-1 1];
%  B2 = [3 1];
%  c = mfoxcontour(10, 2, index, A, A1, B1, A2, B2);

Nvar    = length(varargin);
epsilon = 1/10;
f  = ones(1,dim);
Q  = -A(1:index(2),2:size(A,2));
b  = 1-A(1:index(2),1)-epsilon;
lb = [];
ub = [];

for i = 1 : Nvar/2
 Ai = cell2mat(varargin(2*(i-1)+1));
 Bi = cell2mat(varargin(2*(i-1)+2));  
 lb = [lb max((Ai(1:index(2*i),1)-1)./Ai(1:index(2*i),2))]+epsilon;
 ub = [ub min(Bi(1:index(2*i+1),1)./Bi(1:index(2*i+1),2))]-epsilon;    
end

options = optimoptions('linprog','Algorithm','interior-point');
out = linprog(f, Q, b, [], [], lb, ub, options);
c = [out.' - 1i * W; out.' + 1i * W] ;
    
