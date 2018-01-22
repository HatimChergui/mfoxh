 %*************************** 2-Dimensions ******************************* 
 index = [0 1 1 1 1 1]; % [0 m m1 n1 m2 n2 ...mM nM]
 z = [1 2]; % [z1...zM]
 A = [1.5 1.0 1.0; 2.0 1.0 1.0]; % [a1 alpha_1,1...alpha_1,M; ...; ap alpha_p,1...alpha_p,M]
 B = [2.0 1.0 1.0]; % [b1 beta_1,1...beta_1,M;...; bq beta_q,1 ...beta_q,M]
 A1 = [-1 1]; % [a1_1 alpha1_1;...;a1_p1 alpha1_p1]
 B1 = [0 1]; % [b1_1 beta1_1;...;b1_q1 beta1_q1]
 A2 = [-1 1];
 B2 = [3 1];
 %c  = [-0.5-10i 1.5-10i;-0.5+10i 1.5+10i];
 c = mfoxcontour(10, 2, index, A, A1, B1, A2, B2); 
 % Tolerence settings
 MaxFunEval = 2*1e5; % increase it to get more precision (especially for more than 2 dimensions)
 AbsTol     = 1e-4;
 tic;
 mfoxh(index, z, A, B, A1, B1, A2, B2, c, MaxFunEval, AbsTol)
 toc;

 %************************** 3-Dimensions ***************************** 
 index = [0 1 1 1 1 1 1 1];
 z = [1 2 0.5];
 A = [1.5 1.0 1.0 1; 2.0 1.0 1.0 1];
 B = [2.0 1.0 1.0 1];
 A1 = [-1 1];
 B1 = [0 1];
 A2 = [-1 1];
 B2 = [3 1];
 A3 = [0 1];
 B3 = [1 1];
 %c  = [-0.5-10i 1.5-10i 0.5-10i;-0.5+10i 1.5+10i 0.5+10i];
 c = mfoxcontour(10, 3, index, A, A1, B1, A2, B2, A3, B3)
% Tolerence settings
 MaxFunEval = 1e5; % increase it to get more precision (especially for more than 2 dimensions)
 AbsTol     = 1e-4;
 tic;
 mfoxh(index, z, A, B, A1, B1, A2, B2, A3, B3, c, MaxFunEval, AbsTol)
 toc;