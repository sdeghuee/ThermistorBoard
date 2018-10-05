%% Thermistor Temp Calculations
% Values and equations from TH-22 Datahseet 
A = 1.2870e-3;
B = 2.3585e-4;
C = 9.4346e-8;
R = 6e3;
T = 1/ (A +B * log(R) + C * (log(R))^3);
T = T - 273.15;

%% Thermistor R Calculations
clear
clc
% Values and equations from TH-22 Datahseet 
A = 1.2870e-3;
B = 2.3585e-4;
C = 9.4346e-8;
T = -80:.1:250;
T = T + 273.15;
alpha = ((A-(1./T))./C);
beta = sqrt(((B./(3*C))^3) + (alpha.^2/4));

R = exp((beta-(alpha./2)).^(1/3) - ((beta+(alpha./2)).^(1/3)));
subplot(2,1,1)
plot(T-273.15,R)
title('Resistance over Temp')
xlabel('Degrees (C)')
ylabel('Resistance (Ohms)')
% The diff of R will tell where in the temperature range the smallest
% change in R occurs
dR = diff(R);
subplot(2,1,2)
plot(T(2:end)-273.15,abs(dR));
title('Magnitude of dR/dT')
xlabel('Degrees (C)')
ylabel('dR/dT')
[dRmin, I] = min(abs(dR));
smallestSwing = T(I) - 273.15; %From this the worst case for ADC is at high temps
deltaRsmall = R(end-7) - R(end);
%% ADC Resolution Choice
% Selected manually to optimize power consumption while meeting spec
Rstatic = 500e3;
% Calculates smallest change in Voltage at the given specification
adcResreq = 1.8*(R(end-7) / (Rstatic + R(end-7)) - R(end)/(Rstatic + R(end)));
% Calculates smallest change in Voltage detectable 
adcRes = 1.8 / (2^24);
fprintf('ADC Resolution is good True/False: %d\n',adcResreq > adcRes)

%% ADC SCLK Choice
bitpch = 24;
fullsample = bitpch * 24 /8;
fclk = 27e6;
fdata = fclk / 2560;

%% Power Calculations
Vin = 28;
Vlogic = 1.8;
Vhigh = 5;
% All values are gotten from the respective datasheets
% Two assumptions are used: Iin to LDO = Iout
% Pin = Pout / Efficiency
% Microcontroller
Imicro = 2.5e-3;
% Thermocouples
Rmin = Rstatic + R(end);
Idivide = (Vlogic / Rmin)*24;
% ADCs Vsupply = 5V, 1V8
Iadc1V8 =  2.5e-3 + .035e-3;
Iadc5 = 9e-3;
%RS485 Vsupply = 5V
Irs485 = 28e-3;
%Level Shifter
Ipshift = 18e-6;
% Switching Supply Calculations
efficiencya = .75;
efficiencyp = .75;
% Total Power Calculation (active)
P1V8 = Vhigh*(24*Idivide + Imicro + 3*Iadc1V8 + 1*Ipshift)/efficiencya;
P5 = Vhigh*(Irs485 + 3*Iadc5 + 3*Ipshift)/efficiencya;
fprintf('Power Consumption(active) 1V8: %1.2f Watts\n',P1V8)
fprintf('Power Consumption(active) 5: %1.2f Watts\n',P5)
fprintf('Power Consumption(active) Total: %1.2f Watts\n',P1V8 + P5)
% Total Power Calculation (passive)
Pp1V8 = Vhigh*(24*Idivide + Imicro + 3*Iadc1V8 + 1*Ipshift)/efficiencyp;
Pp5 = Vhigh*(Irs485 + 3*Ipshift)/efficiencyp;
fprintf('Power Consumption(passive) 1V8: %1.2f Watts\n',Pp1V8)
fprintf('Power Consumption(passive) 5: %1.2f Watts\n',Pp5)
fprintf('Power Consumption(passive) Total: %1.2f Watts\n',Pp1V8 + Pp5)

