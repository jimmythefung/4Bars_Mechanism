clear
clc

% Simulation
exit = 0;
Iterations = 0;
flagit = 0;
keptAngles = [];

while exit~=2
    clear
    exit = 0;
    
    % Prescribed positions
    R1 = 1000*(0 - 0.045*i);
    R2 = 1000*(0.0229 - 0.0288*i);
    R3 = 1000*(0.05 - .055*i);

    % Prescribed movements, del. del_i = R_i -R1
    del1 = R1 - R1;
    del2 = R2 - R1;
    del3 = R3 - R1;
    
    % Generate 6 random angles = [ beta2, beta3, beta2, beta3, alpha2,
    % alpha3]
    randAngle = round(360*rand(1,6));
    
    while (randAngle(1)>100) || (randAngle(1) > randAngle(2)) || (randAngle(3) > randAngle(4)) || (randAngle(6) > randAngle(5))
        randAngle = round(360*rand(1,6)); % Generate 6 random angle  
    end
    
    alpha2 = degtorad(randAngle(5));
    alpha3 = degtorad(randAngle(6));   
    
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

    if (mZ1 < 1) || (mZ2 < 1) || (mZ3 < 1) || (mZ4 < 1)
        flag = 1;
    else
        flag = 0;
    end
    
    if flag == 0
        % Check for real solutions only
        z1 = abs(Z1);
        z2 = abs(Z2);
        z3 = abs(Z3);
        z4 = abs(Z4);
        Q2 = -180:10:180;
        for j = 1:length(Q2)
            Q2r = degtorad(Q2(j));
            Aeq = 2*z2*z4*cos(Q2r) - 2*z1*z4;
            Beq = 2*z2*z4*sin(Q2r);
            Ceq = (z3)^2 - (z1)^2 - (z4)^2 - (z2)^2 + 2*z1*z2*sin(Q2r);
            RHS = -Ceq/((Aeq^2 + Beq^2)^(1/2));
            if RHS > 1
                flag = 1;
                %flagit = flagit + 1;
                break
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Define the angles requirements here
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check:
    % Grashoff conditions
    % shortest link is input dyad
    % linkage grounds located above tread mill
    
    S = barlengths(1);
    L = barlengths(4);
    P = barlengths(2);
    Q = barlengths(3);

    
    height = abs(abs(imag(A0)) - abs(imag(B0)));
    width = abs(abs(real(A0)) - abs(real(B0)));
    
    % Set constrants here
    minHeight = -.055*1000;
    minLeft = -.001*1000;
    LMax = abs(Z2+Z5);
    RMax = abs(Z4+Z6);
    
    if (flag == 0)                        &...
       (S + L < P + Q)                    &...% Grashoff condition
       ((S == abs(Z2)) || (S == abs(Z4))) &...% Shortest length is dyad
       ( imag(A0) - mZ2  > minHeight  )   &...% Left dyad doent hit Treadmill
       ( imag(B0) - mZ4  > minHeight  )   &...% Right dyad doesnt hit Treadmill
       ( imag(A0) - LMax > minHeight  )   &...% Coupler tip doesn't hit Trademill
       ( imag(B0) - RMax > minHeight  )   &...% Coupler tip doesn't hit Trademill
       S > 10 &...
       mZ4 < 100 &...
       mZ3 < 100 &...
       1

        % Check solution continuity
        exit = 2;    % assume solution continous
        Q4 = [];
        Q4r = [];    % Theta 4 in radian
        Q4d = [];    % Theta 4 in degree
        Q2d = [];    % Theta 2 in degree
        syms t4;
        Q2 = 0:20:360; % Theta 2 input angles
        for i = 1:length(Q2)
            % Convert input angle theata2 to radian
            Q2r = degtorad(Q2(i));


            % Set up the equation to be solved by solver
            Aeq = 2*mZ2*mZ4*cos(Q2r) - 2*mZ1*mZ4;
            Beq = 2*mZ2*mZ4*sin(Q2r);
            Ceq = (mZ3)^2 - (mZ1)^2 - (mZ4)^2 - (mZ2)^2 + 2*mZ1*mZ2*cos(Q2r);
            %alpha = atan(Aeq/Beq);
            alpha = asin(Aeq/((Aeq^2 + Beq^2)^(1/2)));


            LHS = sin(alpha + t4);
            RHS = -1*Ceq/((Aeq^2 + Beq^2)^(1/2));

            % Solve for theata4
            symAns = solve(LHS == RHS, t4); % Solved symbolically
            numAns = vpa(symAns); % Take numerical value
            Q4(i) = numAns(1);
            fprintf('Solving Theta4= %.2f rad, input angle Q2: %.2f rad\n', Q4(i),Q2r);

            if (i>1) & ( abs(Q4(i)-Q4(i-1))  > 1)
                exit = 0
                clearvars Q4 Q4r Q4d Q2d Aeq Beq Ceq alpha LHS RHS
                break
            end
      

        end

        
    end
    
    
    if exit == 2
        keptAngles = [randAngle];
        break
    end

end




% Output
A0
B0
dyad1=mZ2
dyad2=mZ4
coupler=mZ3
keptAngles

%plot
%close all
clf
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
pl2=plot([A0x, Ax], [A0y, Ay], '-bO', 'LineWidth', 5); 

% A to B
ABx = Ax + real(Z3);
ABy = Ay + imag(Z3);
pl4=plot([Ax, ABx], [Ay, ABy], '-rO', 'LineWidth', 5); 

% B0 to B
Bx = B0x + real(Z4);
By = B0y + imag(Z4);
pl5=plot([B0x, Bx], [B0y, By], '-bO', 'LineWidth', 5); 

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



