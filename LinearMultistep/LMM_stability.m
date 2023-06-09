% Script to generate and plot the linear stability regions for LMM methods:
% {Adams-Bashforth, Adams-Moulton, Backwards Differentiation Formulas}.
%
% D.R. Reynolds
% Math 6321 @ SMU
% Fall 2016

% clear

% set LMM coefficients for each method
%   Adams-Bashforth
AB_q0_a = [1];
AB_q0_b = [0,1];
AB_q1_a = [1];
AB_q1_b = [0,3,-1]/2;
AB_q2_a = [1];
AB_q2_b = [0,23,-16,5]/12;
AB_q3_a = [1];
AB_q3_b = [0,55,-59,37,-9]/24;

%   Adams-Moulton
AM_q0_a = [1];
AM_q0_b = [1,0];
AM_q1_a = [1];
AM_q1_b = [1,1]/2;
AM_q2_a = [1];
AM_q2_b = [5,8,-1]/12;
AM_q3_a = [1];
AM_q3_b = [9,19,-5,1]/24;

%   BDF
BDF_p1_a = [1];
BDF_p1_b = [1];
BDF_p2_a = [4,-1]/3;
BDF_p2_b = [2/3];
BDF_p3_a = [18,-9,2]/11;
BDF_p3_b = [6/11];
BDF_p4_a = [48,-36,16,-3]/25;
BDF_p4_b = [12/25];
BDF_p5_a = [300,-300,200,-75,12]/137;
BDF_p5_b = [60/137];
BDF_p6_a = [360,-450,400,-225,72,-10]/147;
BDF_p6_b = [60/147];


% set the thetas resolution
nthetas = 500;
thetas = linspace(0,2*pi,nthetas);

% plot the A-B stability regions, one at a time
fprintf('\nAdams-Bashforth stability regions:\n');
figure(1)
fprintf('   q = 0:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AB_q0_a,AB_q0_b);
end
fill(x,y,'b'), axis([-6,2,-4,4]), grid on
hold on
pause
fprintf('   q = 1:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AB_q1_a,AB_q1_b);
end
fill(x,y,'r'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 2:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AB_q2_a,AB_q2_b);
end
fill(x,y,'k'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 3:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AB_q3_a,AB_q3_b);
end
fill(x,y,'g'), axis([-6,2,-4,4]), grid on
legend('q=0','q=1','q=2','q=3');
plot([-6,2],[0,0],'k--','LineWidth',1)
plot([0,0],[-4,4],'k--','LineWidth',1)
hold off
title('Adams-Bashforth Stability Regions');
pause


% plot the A-M stability regions, one at a time
fprintf('\nAdams-Moulton stability regions:\n');
fprintf('   q = 0:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AM_q0_a,AM_q0_b);
end
figure(2)
plot(x,y,'k-'), axis([-6,2,-4,4]), grid on
hold on
pause
fprintf('   q = 1:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AM_q1_a,AM_q1_b);
end
plot(x,y,'r-'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 2:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AM_q2_a,AM_q2_b);
end
fill(x,y,'b'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 3:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),AM_q3_a,AM_q3_b);
end
fill(x,y,'g'), axis([-6,2,-4,4]), grid on
legend('q=0','q=1','q=2','q=3');
plot([-6,2],[0,0],'k--','LineWidth',1)
hold off
title('Adams-Moulton Stability Regions');
pause


% plot the BDF stability regions, one at a time
fprintf('\nBDF stability regions:\n');
fprintf('   q = 0:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p1_a,BDF_p1_b);
end
figure(3)
plot(x,y,'k-'), axis([-6,2,-4,4]), grid on
hold on
pause
fprintf('   q = 1:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p2_a,BDF_p2_b);
end
plot(x,y,'r-'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 2:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p3_a,BDF_p3_b);
end
plot(x,y,'b-'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 3:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p4_a,BDF_p4_b);
end
plot(x,y,'g-'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 4:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p5_a,BDF_p5_b);
end
plot(x,y,'c-'), axis([-6,2,-4,4]), grid on
pause
fprintf('   q = 5:\n')
for i=1:nthetas
   [x(i),y(i)] = LMM_stab_eval(thetas(i),BDF_p6_a,BDF_p6_b);
end
plot(x,y,'m-'), axis([-6,2,-4,4]), grid on
legend('q=0','q=1','q=2','q=3','q=4','q=5');
plot([-6,2],[0,0],'k--','LineWidth',1)
plot([0,0],[-4,4],'k--','LineWidth',1)
hold off
title('BDF Stability Regions');

% end of script