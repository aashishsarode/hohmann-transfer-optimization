  % ------------------------------------------------------------------------------
    %
    %                           procedure combined
    %
    %  this procedure calculates the delta v's for a hohmann transfer for either
    %    circle to circle, or ellipse to circle.
    %
    %  author        : david vallado                  719-573-2600   5 may  2012
    %
    %  inputs          description                    range / units
    %    rinit       - initial position magnitude     er
    %    rfinal      - final position magnitude       er
    %    einit       - eccentricity of first orbit
    %    efinal      - eccentricity of final orbit
    %    nuinit      - true anomaly of first orbit    0 or pi rad
    %    nufinal     - true anomaly of final orbit    0 or pi rad, opp of nuinit
    %
    %  outputs       :
    %    deltava     - change in velocity at point a  er / tu
    %    deltavb     - change in velocity at point b  er / tu
    %    dttu        - time of flight for the trans   tu
    %
    %  locals        :
    %    sme1        - mech energy of first orbit     er2 / tu
    %    sme2        - mech energy of transfer orbit  er2 / tu
    %    sme3        - mech energy of final orbit     er2 / tu
    %    vinit       - velocity of first orbit at a   er / tu
    %    vtransa     - velocity of trans orbit at a   er / tu
    %    vtransb     - velocity of trans orbit at b   er / tu
    %    vfinal      - velocity of final orbit at b   er / tu
    %    ainit       - semimajor axis of first orbit  er
    %    atrans      - semimajor axis of trans orbit  er
    %    afinal      - semimajor axis of final orbit  er
    %
    %  coupling      :
    %    none.
    %
    %  references    :
    %    vallado       2007, 352-359, ex 6-7
  % ------------------------------------------------------------------------------

  function [deltai1, deltai2,deltava, deltavb, dttu, vinit, vfinal] = comman( rinit, rfinal, einit, efinal, nuinit, nufinal, deltai,Vta,Vtb )

    show = 'y';
    % --------------------  initialize values   ------------------- }
    mu = 1.0; % canonical

    if show == 'y'
        fprintf(1,'rinit %11.7f ER %11.7f km rfinal %11.7f ER %11.7f km \n',rinit, rinit*6378.137 , rfinal, rfinal*6378.137  );
    end
    
    ainit  = (rinit * (1.0 + einit * cos(nuinit))) / (1.0 - einit * einit );
    atran  = (rinit + rfinal) * 0.5;
    afinal = (rfinal * (1.0 + efinal * cos(nufinal))) / (1.0 - efinal * efinal );
    sme1 = -mu / (2.0*ainit);
    sme2 = -mu / (2.0*atran);

    if show == 'y'
        fprintf(1,'ainit %11.7f ER %11.7f km \n',ainit, ainit*6378.1363 );
        fprintf(1,'atran %11.7f ER %11.7f km \n',atran, atran*6378.1363  );
        fprintf(1,'afinal %11.7f ER %11.7f km \n',afinal, afinal*6378.1363  );
    end

    % -----------------  find delta v at point a  ----------------- }
    vinit  = sqrt( 2.0*( mu/rinit + sme1 ) );
  
    %     fpa2a= atan( ( e2*sin(nu2a) ) / ( 1.0 + e2*cos(nu2a) ) );
    %     fpa1 = atan( ( einit*sin(nuinit) ) / ( 1.0 + einit*cos(nuinit) ) );
    %     deltava= sqrt( vtransa*vtransa + vinit*vinit - 2.0*vtransa*vinit* ...
    %                     ( sin(fpa2a)*sin(fpa1)+cos(fpa2a)*cos(fpa1)*cos(deltai)));

    % -----------------  find delta v at point b  ----------------- }
    vfinal = sqrt( mu/rfinal );  % assumes circular
    
    %     fpa2b= atan( ( e2*sin(nu2b) ) / ( 1.0 + e2*cos(nu2b) ) );
    %     fpa3 = atan( ( efinal*sin(nufinal) ) / ( 1.0 + efinal*cos(nufinal) ) );

    if show == 'y'
        vkmps = 7.905366149846074;
        fprintf(1,'vinit %11.7f er/tu %11.7f km/s vfinal %11.7f  %11.7f \n',vinit, vinit*vkmps, vfinal, vfinal*vkmps );
        fprintf(1,'vtransa %11.7f er/tu %11.7f km/s vtransb %11.7f  %11.7f \n',Vta, Vta*vkmps, Vtb, Vtb*vkmps );
    end

    %     deltavb= sqrt( vtransb*vtransb + vfinal*vfinal - 2.0*vtransb*vfinal* ...
    %( sin(fpa2b)*sin(fpa3)+cos(fpa2b)*cos(fpa3)*cos(deltai)));

    % -------------- find proportions of inclination change ---------------
    % ----------------- this is the approximate approach ------------------
    ratio = rfinal/rinit;
    s = 1.0/deltai * atan(sin(deltai)/(ratio^1.5 + cos(deltai) ) );
    if show == 'y'
        fprintf(1,' s %11.7f \n', s );
    end
    deltai1 = s*deltai;
    deltai2 = (1.0-s)*deltai;
 
    deltava= sqrt( vinit.^2  + Vta.^2 - 2.0*vinit.*Vta.*cos(deltai1) );
    deltavb= sqrt( vfinal.^2 + Vtb.^2 - 2.0*vfinal.*Vtb.*cos(deltai2) );

    dttu= pi * sqrt( (atran * atran * atran)/mu );
    
    % ----------------- figure orientation of the firings -----------------
    gam1 = acos( -(vinit.^2+deltava.^2-Vta.^2 ) / (2.0.*vinit.*deltava) );
    gam2 = acos( -(Vtb.^2+deltavb.^2-vfinal.^2 ) / (2.0.*Vtb.*deltavb) );
    