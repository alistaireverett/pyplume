#!/usr/bin/env python

import numpy as np

def melt(Y,T):
    a = lambda1*(GamT*c_w-GamS*c_i)

    b = GamS*c_i*(lambda1*Y[4]-lambda2-lambda3*T+iceTemp-(L/c_i))-GamT*c_w*(Y[3]-lambda2-lambda3*T)

    c = GamS*Y[4]*(c_i*(lambda2+lambda3*T-iceTemp)+L)

    Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5))
    Tb   = lambda1*Sb+lambda2+lambda3*T 
    mdot = GamS*(Cd**0.5)*Y[2]*(Y[4]-Sb)/Sb 

    return Sb, Tb, mdot

def melt2(speed,P,S,T):

    Aa = -gammaS * speed * cI * a + a * c0 * gammaT * speed

    Bb = -gammaS * speed * L + gammaS * speed * S * cI * a
    Bb = Bb - gammaS * speed * cI * (b + c*P) + gammaS * speed * cI * TI
    Bb = Bb - c0 * gammaT * speed * T + c0 * gammaT * speed * (b + c*P)

    Cc = gammaS * speed * S * L + gammaS * speed * S * cI * (b + c*P) + gammaS * speed * S * (-cI * TI)

    loc_Sb = (-Bb + (Bb**2 - 4.0*Aa*Cc)**0.5)/(2.0*Aa)
    if (loc_Sb < 0.0):
        loc_Sb = (-Bb - (Bb**2 - 4.0*Aa*Cc)**0.5)/(2.0*Aa)

    loc_Tb = a*loc_Sb + b + c*P
    loc_meltrate = (c0 * gammaT * speed * (T - loc_Tb))/(L + cI * (loc_Tb - TI))

    return loc_Sb, loc_Tb, loc_meltrate

if __name__=="__main__":
    c0    = 3994.
    cI     = 2009.
    L      = 334000.
    a = -0.0573
    b = 0.0832
    c = 0.000761
    gammaS   = 0.022
    gammaT   = 0.00062
    Cd     = 0.0025
    backgroundVel = 0.001
    TI  = -0.

    vel=0.1
    sal=30.
    temp=4.
    pressure=0.

    Sb, Tb, mdot = melt2(vel,pressure,sal,temp)
    print Sb, Tb, mdot
    assert Tb == a*Sb+b+c*pressure, "AAAARGH"
    print mdot*Sb, vel*gammaS*(sal-Sb)
    assert abs(mdot*Sb-vel*gammaS*(sal-Sb))<1e-5, "ooooh"
    print mdot*L + mdot*cI*(Tb-TI),c0*vel*gammaT*(temp-Tb)
    assert abs(mdot*L + mdot*cI*(Tb-TI)-c0*vel*gammaT*(temp-Tb))<1e-5, "ooooh"
    print Sb, Tb, mdot
