clc; clear; close all;

% Take Input From User:
% Take array of Resistors from user:
R = input('Enter values of n resistances in ohm, [R1...Rn] = ');
disp('==============================================================================');

% Check if Number of Resistor is less than three or having zero values,
%       or 
while length(R)<3 || sum(R==0) ~= 0 || ~isreal(R)
    disp('No. of Resistor must not be less than 3, each having non-zero real value!')
    disp('TRY AGAIN......')
    R = input('Enter values of n resistances in ohm, [R1...Rn] = ');
    disp('==============================================================================');
end

% Number of Resistors:
n = length(R);

% Take Capacitor from user:
C = input('Enter values of n capacitances in farads, [C1...Cn] = ');
disp('==============================================================================');

% Check if Number of cappacitors is equal to the No. of Resistors
% and the value of each is less non-zero and real.
while length(C) ~= n || sum(C==0) ~= 0 || ~isreal(R)
    disp('No. of Capacitances must be equal to No. of Resistances')
    disp('Each having non-zero real value!')
    disp('TRY AGAIN......')
    C = input('Enter values of n capacitances in farads, [C1...Cn] = ');
    disp('==============================================================================');
end


% Take Inductances from user:
L = input('Enter values of n-1 inductances in henries, [L1...Ln-1] = ');
disp('==============================================================================');

% Check if Number of Inductors is one less than the No. of Resistors & each 
%  having non-zero real value.
while length(L) ~= (n-1) || sum(L==0) ~= 0 || ~isreal(L)
    disp('No. of Inductor must be one less than to No. of Resistors')
    disp('Each having non-zero real value!')
    disp('TRY AGAIN......')
    L = input('Enter values of n capacitances in henries, [L1...Ln] = ');
    disp('==============================================================================');
end

% Take voltage Magnitudes from User
Vmag = input('Enter the amplitudes of the two voltage sources in volt, Vmag = ');
disp('==============================================================================');

% Check if No. of voltage magnitude sources are two  & each 
%  having non-zero real value.
while length(Vmag) ~= 2 || sum(Vmag==0) ~= 0 || ~isreal(Vmag)
    disp('There must be only two Voltage Magnitudes')
    disp('Each having non-zero real value!')
    disp('TRY AGAIN......')
    Vmag = input('Enter the amplitudes of the two voltage sources in volt, Vmag = ');
    disp('==============================================================================');
end

% Take voltage Phases from User
Vphase = input('Enter the phases of the two voltage sources in degrees, Vphase = ');
disp('==============================================================================');

% Check if No. of voltage magnitude sources are two  & each 
%  having real value.
while length(Vphase) ~= 2 || ~isreal(Vphase)
    disp('There must be only two Voltage Phases')
    disp('Each having real value!')
    disp('TRY AGAIN......')
    Vphase = input('Enter the phases of the two voltage sources in degrees, Vphase = ');
    disp('==============================================================================');
end
% Coverting Degrees in Radians:
Vphase = Vphase*pi/180;

% Take frequency from User
Freq = input('Enter the frequency in rad/s, Freq = ');
disp('==============================================================================');

% Take alpha from User
alpha = input('Enter the Alpha (between 0 and 2) to be used with the Successive Relaxation Method = ');
disp('==============================================================================');

% Computation:
% Evaluating Inductive & Capacitive Reactances:
XL = Freq*L;
XC = 1./(Freq*C);

% Impedance Matrix, Z:
Z = zeros(n);

% Kindly check word/pdf file for these variables, Equations (11, 12 & 13):
z_11 = R(1)-j*XC(1)+j*XL(1);
z_ii = j*XL(1:end-1) + R(2:end-1) - j*XC(2:end-1) + j*XL(2:end);
z_nn = R(end)-j*XC(end)+j*XL(end);

z_diag = [z_11 z_ii z_nn];

Z = [zeros(n-1, 1) eye(n-1).*(-j*XL); zeros(1, n)] + ... % Upper Diagonal 
    [zeros(1, n); eye(n-1).*(-j*XL) zeros(n-1, 1)] + ... % Lower Diagonal
    eye(n).*z_diag;                                      % Diagonal Elements

disp('The impedance matrix to solve the mesh current is:')
disp(Z)
disp('******************************************************************************');

% Coverting Voltage sources into complex numbers:
V = zeros(n,1);
temp = (Vmag.*(cos(Vphase) + j*sin(Vphase)));
V(1) = temp(1);
V(end) = temp(end);
disp('The voltage vector to solve the mesh current is:');
disp(V);
disp('******************************************************************************');

I = Z\V;
disp('Exact=');
disp(I)
disp('******************************************************************************');

% Checking Diagonal Domanance:
for i = 1:n
    row = Z(i, :);              % Extract a row matrix
    if abs(row(i)) < sum([row(1:i-1) row(i+1:end)])
        disp("The matrix is not strictly diagonally dominant at row " + string(i))
    end
end

% Successive Iteration M
er = 1;                     % Relative Error
tol = 0.001;                % Tolerance to the relative error
Xapp = zeros(n,1);          % Approximated Solution for Electric Current
Xp = Xapp;                  % Previous Iteration

D = eye(n).*diag(Z);        % Diagonal Impedance Matrix
b = V;                      % Vector b
L = tril(Z, -1);                % Lower Diagonal Matrix
U = conj(L');               % Upper Diagonal Matrix

while er > tol
    % Courtesy: https://en.wikipedia.org/wiki/Successive_over_relaxation
    Xapp = inv(D+alpha*L)*(alpha*b - (alpha*U+(alpha-1)*D)*Xapp);
    er = sum(abs((Xapp - Xp)./Xapp));
    temp = Xp;
    Xp = Xapp;
end

disp('Approximated=');
disp(Xapp)
disp('******************************************************************************');

disp('Relative Error')
ea = abs((Xapp - temp)./Xapp);  % Relative Error is given by: (current-previous)/current
disp(ea)
disp('******************************************************************************');

disp('Significant Digits:')
temp = ea;
while sum(temp < 1)
    temp = 10*temp;                 % Keep multiplying by 10 until the number become > 1
end
disp(round(temp));
disp('******************************************************************************');

disp('Displaying in polar form:')
disp('Magnitude in amperes:');
disp(abs(Xapp))
disp('Phase in degrees');
disp(180*angle(Xapp)/pi);
disp('******************************************************************************');
