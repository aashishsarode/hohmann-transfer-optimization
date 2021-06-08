% ------------------------------------------------------------------------------
%
% Copyright (c) 2021, Aashish Sarode
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Title: Trajectory Modelling of a satellite to GEO using Hohmann transfer 
%                           and Optimization using Particle Swarm Technique
% 
% ------------------------------------------------------------------------------

constastro;
    vtransa=readmatrix('Inputs.xlsx','Range','A2:A1001');   %  velocity of trans orbit at a in kmps
    vtransb=readmatrix('Inputs.xlsx','Range','B2:B1001');   %  velocity of trans orbit at b in kmps
    alt1=300;   % initial orbit altitude in km
    alt2=35786.2;   % fianl orbit altitude in km
    inc1=28.5;  % initial inclination angle in degrees
    inc2=0.0;   % fianl incilination angle in degrees
    rinit  = (re + alt1)/re;        % radius of perigee in earth radius (ER)
    rfinal = (re + alt2)/re;    % radius of apogee in earth radius (ER)
    fprintf(1,'from radius %11.5f to %11.5f \n', rinit*re, rfinal*re);
    einit  = 0.0;   % Initial Eccentricity
    efinal = 0.0;   % Final Eccentricity
    iinit = inc1/rad;   % Initial Inclination angle in radians
    ifinal=  inc2/rad;  % Final Inclination angle in radians
    fprintf(1,'from inclination %11.5f to %11.5f \n', iinit*rad, ifinal*rad);
    deltai = -ifinal + iinit;   % total change in inclination angle
    nuinit =  0.0/rad;  % true anomaly of first orbit    0 or pi rad
    nufinal=  180.0/rad;    % true anomaly of final orbit    0 or pi rad, opp of nuinit
    Vta = vtransa./velkmps; % Velocity of transfer orbit at a in er/tu
    Vtb= vtransb./velkmps;  % Velocity of transfer orbit at b in er/tu

    [deltai1, deltai2, deltava, deltavb, dttu, vinit, vfinal] = comman( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai,Vta,Vtb);

    fprintf(1,' deltava  %11.7f er/tu \n',deltava);
    fprintf(1,' deltavb  %11.7f er/tu \n',deltavb);
    fprintf(1,' deltai1  %11.7f  deltai2  %11.7f  \n',deltai1*rad,deltai2*rad);

    atran = (rinit+rfinal)*re/2.0;  % Semi major axis of the trasnfer orbit
    etran = re*(rfinal-rinit)/(2*atran);    % Eccentricty of the trasfer orbit
    ptran = atran*(1.0-etran^2);    % perigee of transfer orbit
    
    deltav = deltava + deltavb; % Summation of change in velocities of a and b in er/tu
    deltav1 = deltav.*velkmps;  % Summation of change in velocities of a and b in km/s
    
    T=readmatrix('Outputs.xlsx','Range','A2:A1001');
    s2=(1/deltai).*atan(sin(deltai)./(((vinit.*Vta)./(vfinal.*Vtb))+cos(deltai)));
    deltai11=s2*rad2deg(deltai);    % Inclination angle at Perigee (Degrees)
    deltai12=(1-s2)*rad2deg(deltai);    % Inclination angle at Apogee (Degrees)
    
    writematrix(deltav1,'Outputs.xlsx','Range','B2');
    writematrix(deltai11,'Outputs.xlsx','Range','C2');
    writematrix(deltai12,'Outputs.xlsx','Range','D2');
    
    %% Regression   
    X=[ones(size(deltav1)) T deltai11];
    [b,~,~,~,stats] = regress(deltav1,X);
    fprintf(1,' Coefficient of determination (R squared) = %11.7f \n',stats(1));

    %% PSO
    pso;
    fprintf(1,' Time Period (Seconds) = %11.7f \n',BestSol.Position(1));
    fprintf(1,' Inclination angle at Perigee (Degrees)  = %11.7f \n',BestSol.Position(2));
    fprintf(1,' Change in velocity at perigee (Km/s) = %11.7f \n',BestSol.Cost);

    %% Figures
    figure;
    plot(deltav1,deltai11,'LineWidth',2.0);
    xlabel('Change in Velocity (km/s)');
    ylabel('Inclination angle at Perigee (Degrees)');
    grid on;
    
    figure;
    plot(deltav1,deltai12,'LineWidth',2.0);
    xlabel('Change in Velocity (km/s)');
    ylabel('Inclination angle at Apogee (Degrees)');
    grid on;

    figure;
    plot(vtransa,deltai11,'LineWidth',2.0);
    xlabel('Change in velocity at Perigee (km/s)');
    ylabel('Inclination angle at Perigee (Degrees)');
    grid on;

    figure;
    plot(vtransb,deltai12,'LineWidth',2.0);
    xlabel('Change in velocity at Apogee (km/s)');
    ylabel('Inclination angle at Apogee (Degrees)');
    grid on;

%% create trajectory graphics
    % load orbital elements arrays, create state vectors and plot orbits
oevi(1) = rinit*re; %  oev(1) = semimajor axis (kilometers)
oevi(2) = 0.0;  %  oev(2) = orbital eccentricity (non-dimensional)  (0 <= eccentricity < 1)
oevi(3) = iinit;    %  oev(3) = orbital inclination (radians)   (0 <= inclination <= pi)
oevi(4) = 0.0;  %  oev(4) = argument of perigee (radians)   (0 <= argument of perigee <= 2 pi)
oevi(5) = 0.0;  %  oev(5) = right ascension of ascending node (radians)  (0 <= raan <= 2 pi)

% determine correct true anomaly (radians)

if (alt2 > alt1)
    
    oevi(6) = 0.0;  %  oev(6) = true anomaly (radians)  (0 <= true anomaly <= 2 pi)
    
else
    
    oevi(6) = 180.0/rad;
    
end

[ri, vi] = orb2eci(mu, oevi);

oevti(1) = atran;
oevti(2) = etran;
oevti(3) = deg2rad(BestSol.Position(2));
oevti(4) = 0.0;
oevti(5) = 0.0;

% determine correct true anomaly (radians)

if (alt2 > alt1)
    
    oevti(6) = 0.0;
    
else
    
    oevti(6) = 180.0/rad;
    
end

[rti, vti] = orb2eci(mu, oevti);

oevtf(1) = atran;
oevtf(2) = etran;
oevtf(3) = 0.0;
oevtf(4) = 0.0;
oevtf(5) = 0.0;

% determine correct true anomaly (radians)

if (alt2 > alt1)
    
    oevtf(6) = 180.0/rad;
    
else
    
    oevtf(6) = 0.0;
    
end

[rtf, vtf] = orb2eci(mu, oevtf);

oevf(1) = rfinal*re;
oevf(2) = 0.0;
oevf(3) = ifinal;
oevf(4) = 0.0;
oevf(5) = 0.0;

% determine correct true anomaly (radians)

if (alt2 > alt1)
    
    oevf(6) = 180.0/rad;
    
else
    
    oevf(6) = 0.0;
    
end

[rf, vf] = orb2eci(mu, oevf);

% compute orbital periods

period1 = 2.0 * pi * oevi(1) * sqrt(oevi(1) / mu);

period2 = 2.0 * pi * oevti(1) * sqrt(oevti(1) / mu);

period3 = 2.0 * pi * oevf(1) * sqrt(oevf(1) / mu);

deltat1 = period1 / 300;

simtime1 = -deltat1;

deltat2 = 0.5 * period2 / 300;

simtime2 = -deltat2;

deltat3 = period3 / 300;

simtime3 = -deltat3;

for i = 1:1:301
    
    simtime1 = simtime1 + deltat1;
    
    simtime2 = simtime2 + deltat2;
    
    simtime3 = simtime3 + deltat3;
    
    % compute initial orbit "normalized" position vector
    
    [rwrk, vwrk] = twobody2 (mu, simtime1, ri, vi);
    
    rp1_x(i) = rwrk(1) / re;
    
    rp1_y(i) = rwrk(2) / re;
    
    rp1_z(i) = rwrk(3) / re;
    
    % compute transfer orbit position vector
    
    [rwrk, vwrk] = twobody2 (mu, simtime2, rti, vti);
    
    rp2_x(i) = rwrk(1) / re;
    
    rp2_y(i) = rwrk(2) / re;
    
    rp2_z(i) = rwrk(3) / re;
    
    % compute final orbit position vector
    
    [rwrk, vwrk] = twobody2 (mu, simtime3, rf, vf);
    
    rp3_x(i) = rwrk(1) / re;
    
    rp3_y(i) = rwrk(2) / re;
    
    rp3_z(i) = rwrk(3) / re;
    
end

 figure(7);
 
 % create axes vectors
 
 xaxisx = [1 1.5];
 xaxisy = [0 0];
 xaxisz = [0 0];
 
 yaxisx = [0 0];
 yaxisy = [1 1.5];
 yaxisz = [0 0];
 
 zaxisx = [0 0];
 zaxisy = [0 0];
 zaxisz = [1 1.5];
 
 figure(7);
 
 hold on;
 
 grid on;
 
 %% plot earth
 
 [x, y, z] = sphere(24);
 h = surf(x, y, z);
 colormap([127/255 1 222/255]);
 set (h, 'edgecolor', [1 1 1]);
 
 %% plot coordinate system axes
 
 plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);
 plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);
 plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);
 
 % plot initial orbit
 
 plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);
 plot3(rp1_x(1), rp1_y(1), rp1_z(1), 'ob');
 
 % plot transfer orbit
 
 plot3(rp2_x, rp2_y, rp2_z, '-b', 'LineWidth', 1.5);
 plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');
 
 % plot final orbit
 
 plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);
 
 xlabel('X coordinate (ER)', 'FontSize', 12);
 ylabel('Y coordinate (ER)', 'FontSize', 12);
 zlabel('Z coordinate (ER)', 'FontSize', 12);
 title('Hohmann Transfer: Initial, Transfer and Final Orbits', 'FontSize', 16);
 
 axis equal;
 view(50, 20);
 rotate3d on;
 print -depsc -tiff -r300 hohmann1.eps