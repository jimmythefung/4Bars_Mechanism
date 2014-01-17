clc
clear

% Enter angle here
 randAngle = [47   185   221   353    74    24];
 
% Prescribed positions
R1 = 1000*(0 - 0.045*i);
R2 = 1000*(0.0229 - 0.0288*i);
R3 = 1000*(0 - .055*i);

% Prescribed movements, del. del_i = R_i -R1
del1 = R1 - R1;
del2 = R2 - R1;
del3 = R3 - R1;

alpha2 = degtorad(randAngle(5));
alpha3 = degtorad(randAngle(6)); 

% Simulation
exit = 0;
Iterations = 0;
flagit = 0;
keptAngles = [];

% Left Dyad
beta2 = degtorad(randAngle(1));
beta3 = degtorad(randAngle(2)); %%%%%%%%% Modify angle input here
a11 = exp(i * beta2) - 1;
a12 = exp(i * alpha2) - 1;
a21 = exp(i * beta3) - 1;
a22 = exp(i * alpha3) - 1;
LHS = [a11, a12; a21, a22];
RHS = [del2; del3];
x1 = linsolve(LHS,RHS);

% Right Dyad
beta2 = degtorad(randAngle(3));
beta3 = degtorad(randAngle(4));
a11 = exp(i * beta2) - 1;
a12 = exp(i * alpha2) - 1;
a21 = exp(i * beta3) - 1;
a22 = exp(i * alpha3) - 1;
LHS = [a11, a12; a21, a22];
RHS = [del2; del3];
x2 = linsolve(LHS,RHS);

% Map linkage locations
A = R1 - x1(2);
A0 = A - x1(1);
B = R1 - x2(2);
B0 = B - x2(1);
AB = B-A;
A0B0 = B0 - A0;
Z1 = A0B0;
Z2 = x1(1);
Z3 = AB;
Z4 = x2(1);
Z5 = x1(2);
Z6 = x2(2);
barlengths = sort([abs(Z1) abs(Z2) abs(Z3) abs(Z4)]);

mA = abs(A);
mA0 = abs(A0);
mB = abs(B);
mB0 = abs(B0);
mAB = abs(AB);
mA0B0 = abs(A0B0);
mZ1 = abs(Z1);
mZ2 = abs(Z2);
mZ3 = abs(Z3);
mZ4 = abs(Z4);
mZ5 = abs(Z5);
mZ6 = abs(Z6);

% Plot of stationary 4-bars
close all
hold on
A0x = real(A0);
A0y = imag(A0);
B0x = real(B0);
B0y = imag(B0);
% figure
% axis([-10 20 -20 20]);

% A0 to A
Ax = A0x + real(Z2);
Ay = A0y + imag(Z2);
pl2=plot([A0x, Ax], [A0y, Ay], '-b^', 'LineWidth', 5); 

% A to B
ABx = Ax + real(Z3);
ABy = Ay + imag(Z3);
pl4=plot([Ax, ABx], [Ay, ABy], '-rO', 'LineWidth', 5); 

% B0 to B
Bx = B0x + real(Z4);
By = B0y + imag(Z4);
pl5=plot([B0x, Bx], [B0y, By], '-b^', 'LineWidth', 5); 

% A to P1
P1x = Ax + real(Z5);
P1y = Ay + imag(Z5);
pl6=plot([Ax, P1x], [Ay, P1y], '-rO', 'LineWidth', 5); 

% B to P1
pl7=plot([Bx, P1x], [By, P1y], '-rO', 'LineWidth', 5); 

% Treadmill Surface
pl8=plot([0, 100], [-55, -55], '-m^', 'LineWidth', 10); 

% % A0 to B0
% pl1=plot([A0x, B0x],[A0y, B0y],'-g*'); % Green
grid on;
xlabel('Horizontal mm');
ylabel('Vertical mm');


%%
% Solves Theta 4 numerically
z1 = abs(Z1);
z2 = abs(Z2);
z3 = abs(Z3);
z4 = abs(Z4);
Q2rr = [];
Q4 = [];
Q4r = [];    % Theta 4 in radian
Q4d = [];    % Theta 4 in degree
Q2d = [];    % Theta 2 in degree
syms t4;

Q2 = 0:20:360; % Theta 2 input angles /////////////////////////////

for i = 1:length(Q2)
    % Convert input angle theata2 to radian
    Q2r = degtorad(Q2(i));
    Q2rr(i) = Q2r;
    
    
    % Set up the equation to be solved by solver
    Aeq = 2*z2*z4*cos(Q2r) - 2*z1*z4;
    Beq = 2*z2*z4*sin(Q2r);
    Ceq = (z3)^2 - (z1)^2 - (z4)^2 - (z2)^2 + 2*z1*z2*cos(Q2r);
    %alpha = atan(Aeq/Beq);
    alpha = asin(Aeq/((Aeq^2 + Beq^2)^(1/2)));
    

    LHS = sin(alpha + t4);
    RHS = -1*Ceq/((Aeq^2 + Beq^2)^(1/2));
    
    % Solve for theata4
    symAns = solve(LHS == RHS, t4); % Solved symbolically
    numAns = vpa(symAns); % Take numerical value
    Q4(i) = numAns(1); 

    fprintf('Solving Theta4= %.2f rad, input angle Q2: %.2f rad\n', Q4(i),Q2r);

end

Q2 = degtorad(Q2);

% Solves Theta 3 numerically
syms t3
Q3 = [];
for i = 1:length(Q2)
    
    % Equation to be solved
    %LHSx = z2*cos(Q2(i)) + z3*cos(t3);
    %RHSx = z1 + z4*cos(Q4(i)); %Q4a or b
    %symAnsX = solve(LHSx == RHSx, t3); % Symbolic answer
    %numAnsX = vpa(symAnsX); % This gives complex solution
    
    LHSy = z2*sin(Q2(i)) + z3*sin(t3);
    RHSy = z4*sin(Q4(i));
    symAnsY = solve(LHSy == RHSy, t3); % Symbolic answer    
    numAnsY = vpa(symAnsY); % Take the numerical value in complex form

    Q3(i) = numAnsY(1);
    
    fprintf('Solving Theta3= %.2f, input angle Q4: %.2f rad\n', Q3(i),Q4(i));

end


%%
% Motion simulation for seat
A0x = real(A0);
A0y = imag(A0);
B0x = real(B0);
B0y = imag(B0);
Q2offset = angle(Z2);
Q3offset = angle(Z3);
Q4offset = angle(Z4);

Q2 = Q2 + Q2offset;
Q3 = Q3 - Q3(1) + Q3offset;
Q4 = Q4 - Q4(1) + Q4offset;

figure
axis([-100 100 -100 100]);
CTx = [];
CTy = [];

for i = 1:length(Q2)
    hold on
    % A0 to B0
    %pl1=plot([A0x, B0x],[A0y, B0y],'-g*'); % Green

    % A0 to A
    Ax = A0x + mZ2*cos(Q2);
    Ay = A0y + mZ2*sin(Q2);
    pZ2=plot([A0x, Ax(i)], [A0y, Ay(i)], '-b^', 'LineWidth', 5);

    % A to B
    ABx = Ax + mZ3*cos(Q3(i));
    ABy = Ay + mZ3*sin(Q3(i));
    pZ3=plot([Ax(i), ABx(i)], [Ay(i), ABy(i)], '-rO', 'LineWidth', 5);

    % B to B0
    pZ4=plot([B0x, ABx(i)], [B0y, ABy(i)], '-b^', 'LineWidth', 5);

    
    % A to coupler tip
     CTx(i) = Ax(i)+ mZ5*cos(angle(Z5)+Q3(i)-Q3(1));
     CTy(i) = Ay(i)+ mZ5*sin(angle(Z5)+Q3(i)-Q3(1));
     pZ5=plot([Ax(i), CTx(i)], [Ay(i), CTy(i)], '-rO', 'LineWidth', 5);   
    
     % coupler tip to B
     pZ6=plot([CTx(i), ABx(i)], [CTy(i), ABy(i)], '-rO', 'LineWidth', 5);


    
    % Treadmill Surface
    pTM=plot([0, 100], [-55, -55], '-m^', 'LineWidth', 10);
    
    hold off
    grid on;
    xlabel('distance cm');
    ylabel('height cm');
    
    % Two ways to animate: 1.getframe 2.pause
    F(i) = getframe;
    %pause(.3)
    
    if i ~= length(Q2)
        delete(pZ2);
        delete(pZ3);
        delete(pZ4);
        delete(pZ5);
        delete(pZ6);
        delete(pTM);
    end
    

end
movie(F,10,10)

%Output
A0
B0
dyad1=mZ2
dyad2=mZ4
coupler=mZ3
