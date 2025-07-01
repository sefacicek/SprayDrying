%
%  Sefa Çiçek - Mechanical Engineering
% Yildiz Technical University
%
clc; clear; close all

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG muG dp mp0 vg rho_milk fat kG Diff R Tg0

%% ---------------------------------------------------------------------------------------------------------------------
% Physical Properties
%% ---------------------------------------------------------------------------------------------------------------------
Tp0 = 303;								        
Tg0 = 403;										
P = 1;											
CpL = 1;									    
CpG = 0.25;										
MWwater = 18;                                   
MWair = 29;                                     
R = 0.0821;                                     
dHev = 540;										
rhoL = 1000;					     			
rhoG = P/R/Tg0*MWair;                           
muG = 2.3*1e-5;                                 
kG = 8*1.0e-6;									
Diff = 1.8*1e-5;                               	
rho_milk  =  440;                               

%% ---------------------------------------------------------------------------------------------------------------------
% Antoine equation
%% ---------------------------------------------------------------------------------------------------------------------
A = 18.3036;
B = 3816.44;
C = -46.13;

%% ---------------------------------------------------------------------------------------------------------------------
% Process stream
%% ---------------------------------------------------------------------------------------------------------------------
Qmilk = 1900/3600;                            	
fat = 4.76/100;                           		
Qfat = Qmilk*fat;                     			
Ql = Qmilk - Qfat;                				
Win = Ql/Qfat;									
Wout = 0.005; 									
dp = 2*1e-4;     								
vp0 = 0.3;										
mp0 = rhoL * pi / 6 * (dp^3);         		    
mdry = mp0 / (1 + Win);                   		
mout = mdry * (1 + Wout);                  		
nAvp = Qmilk / mp0;                             

%% ---------------------------------------------------------------------------------------------------------------------
% Gas stream
%% ---------------------------------------------------------------------------------------------------------------------
Gdry = 20;										
D = 5.5;     									

%% ---------------------------------------------------------------------------------------------------------------------
% Momentum Balance
%% ---------------------------------------------------------------------------------------------------------------------
vg = Gdry / rhoG / 3.14 / D^2*4;                
vs = (rhoL-rhoG) * dp^2 * 9.81 / 18 / muG;      
vp = vg + vs;                                   
vs0 = vp0 - vg;                                 

%% ---------------------------------------------------------------------------------------------------------------------
% Mass and heat transfer
%% ---------------------------------------------------------------------------------------------------------------------
Re = rhoG * vs * dp / muG;                      
Pr = muG * CpG /kG; 							
Sc =  muG / rhoG / Diff;      			    	
Nu = 2 + 0.4 * Re^0.5 * Pr^(1/3);          		
Sh = 2 + 0.4 * Re^0.5 * Sc^(1/3);              	
h = Nu * kG / dp;                          		
Kc = Sh * Diff / dp;                   			
Kpw = Kc / R / Tg0 * MWwater;      				

%% ---------------------------------------------------------------------------------------------------------------------
% Integration
%% ---------------------------------------------------------------------------------------------------------------------
tspan = linspace(0,5.5,1001);
y0 = [mp0 0 Tp0 Tg0 vs0 0];

[t,y] = ode15s(@solver,tspan,y0,odeset("MaxStep",1e-4,"Refine",5));

%% ---------------------------------------------------------------------------------------------------------------------
% Results
%% ---------------------------------------------------------------------------------------------------------------------
m = y(:,1);
time = t(find(m<2e-10,1));
lenght = y(find(m<2e-10,1),6);

% Moisture ratio hesaplama
moisture = (y(:,1) - mdry) / mdry;   % kg water / kg fat

% Absolute humidity (mutlak nem) hesaplama
absolute_humidity = y(:,2) / Gdry;   % kg water / kg dry air

%% ---------------------------------------------------------------------------------------------------------------------
% Plots
%% ---------------------------------------------------------------------------------------------------------------------
cc = winter(6);
close all

figure

subplot(3,2,1)
plot(y(:,6), y(:,1),"LineWidth",1.4,"Color",cc(1,:))
hold on
grid on
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]')
ylabel('Particle mass [kg]')
ylim([min(y(:,1))/5,max(y(:,1))*1.1])
legend('Particle Mass','Useful Length')
hold off

subplot(3,2,2)
plot(y(:,6),y(:,3),"LineWidth",1.4,"Color",cc(1,:))
hold on
grid on
plot(y(:,6),y(:,4),"LineWidth",1.4,"Color",cc(4,:))
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]')
ylabel('Temperature [K]')
legend('Particle','Air','Useful Length')
ylim([min(y(:,3))/1.05,max(y(:,4))*1.05])
hold off

subplot(3,2,3)
plot(y(:,6),y(:,5),"LineWidth",1.4,"Color",cc(1,:))
hold on
grid on
plot(y(:,6), vg + y(:,5),"LineWidth",1.4,"Color",cc(4,:))
plot(y(:,6), vg * ones(1 * length(tspan)),"LineWidth",1.4,"Color",cc(6,:))
plot([lenght lenght],[-1000 1000],"LineStyle","-.","Color","r")
xlabel('Axial coordinate [m]')
ylabel('Velocity [m/s]')
legend('vs','vp','vg')
ylim([min(y(:,5))/5,max(vg + y(:,5))*1.1])
hold off

subplot(3,2,4)
plot(tspan,y(:,6),"LineWidth",1.4,"Color",cc(1,:))
hold on
plot([time time],[-1000 1000],"LineStyle","-.")
xlabel('Time [s]','FontWeight','bold','FontSize',14)
ylabel('Axial coordinate [m]','FontWeight','bold','FontSize',14)
legend('Axial Position','Useful Length')
ylim([min(y(:,6)),max(y(:,6))*1.1])
hold off

subplot(3,2,5) % Moisture ratio
plot(y(:,6), moisture, "LineWidth", 1.4, "Color", cc(1,:))
hold on
grid on
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]')
ylabel('Moisture ratio (kg water/kg fat)')
legend('Moisture Ratio', 'Useful Length')
ylim([min(moisture)/1.1, max(moisture)*1.1])
hold off

subplot(3,2,6) % Absolute humidity
plot(y(:,6), absolute_humidity, "LineWidth", 1.4, "Color", cc(4,:))
hold on
grid on
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]')
ylabel('Absolute Humidity (kg H2O/kg dry air)')
legend('Absolute Humidity', 'Useful Length')
ylim([min(absolute_humidity)/1.1, max(absolute_humidity)*1.1])
hold off

%% ---------------------------------------------------------------------------------------------------------------------
% Function
%% ---------------------------------------------------------------------------------------------------------------------
function spraydryer = solver(t,x)

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG muG dp mp0 vg rho_milk fat kG Diff R
    
    mp = x(1);
    Gvap = x(2);
    Tp = x(3);
    Tg = x(4);
    vs = x(5);
    z = x(6);
    
    mmin = mp0 * fat;
    mp = max(mmin, mp);
    
    x_H2O = (mp - mp0*fat) / mp;
    
    rho_avg = x_H2O * rhoL + (1 - x_H2O) * rho_milk;
    V = mp / rho_avg;
    dp = (6 * V / pi)^(1/3);
    Dmin = 83.87e-6;
    dp = max(Dmin, dp);
    
    rhoG = P / R / Tg * MWair; 
    Re = abs(rhoG * vs * dp / muG);                     
    Pr = muG * CpG / kG; 								
    Sc = muG / rhoG / Diff;      						
    Nu = 2 + 0.4 * Re^0.5 * Pr^(1/3);                  	
    Sh = 2 + 0.4 * Re^0.5 * Sc^(1/3);                  	
    f = (1 + 0.14*Re^0.7);                              
    h = Nu * kG / dp;                             		
    Kc = Sh * Diff/ dp;                      			
    Kpw = Kc / R / Tg * MWwater;         			    
    
    Pw = MWair / MWwater * Gvap * P / Gdry;
    P0 = exp(A - B / (Tp + C)) / 760;                   
    As = pi * dp^2;
    
    if mp > mmin + mmin/10000
        spraydryer(1) = Kpw * As * (Pw - P0);
    else
        spraydryer(1) = 0;
    end
    
    spraydryer(2) = -Kpw * As * (Pw-P0) * nAvp;
    spraydryer(3) = (h * As * (Tg - Tp) + spraydryer(1) * dHev) / mp / CpL;
    spraydryer(4) = -h * As * (Tg - Tp) * nAvp / (Gvap + Gdry) / CpG;
    spraydryer(5) = ((1 - rhoG / rho_avg) * 9.81 - 3 * f * vs * muG * pi * dp /mp - vs * spraydryer(1) / mp);
    spraydryer(6) = (vs+vg);
    
    spraydryer = spraydryer';
end
