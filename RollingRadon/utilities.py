import numpy as np
from math import pi

def deg2rad(deg):
    r"""
    converts degrees to radians
    """
    return (pi/180.) * deg

def rad2deg(rad):
    r"""
    converts radians to degrees
    """
    return (180./pi) * rad

def pointdistance(p1, p2):
    r"""
    calculates the distance between two points p1 and p2
    """
    for i in range(len(p1)):
        distance = np.linalg.norm(p1-p2)
    return distance

def distance_vector(x,y):
    r"""
    distance vector given x,y coordinates
    """
    distance=np.empty(x.shape)
    distance[0]=0.0
    for i in range(0,len(x)-1):
        distance[i+1]=pointdistance(np.array([x[i],y[i]]),np.array([x[i+1],y[i+1]]))

    return distance

def polarstereo_fwd(ϕ,λ,ϕc,λ0,e=0.08181919,a=6378137.0):
    r"""
    convert lat lon to polar stereographic projection
    """
    ϕ=deg2rad(ϕ)
    ϕc=deg2rad(ϕc)
    λ=deg2rad(λ)
    λ0=deg2rad(λ0)

    t = np.tan(pi / 4 - ϕ / 2 ) / (1 - e * np.sin(ϕ)) / (1 + e * np.sin(ϕ))**(e / 2)
    tc = np.tan(pi / 4 - ϕc / 2) / ((1 - e * np.sin(ϕc)) / (1 + e * np.sin(ϕ_c)))**(e / 2)
    mc = np.cos(ϕc) / np.sqrt(1 - e**2 * (sin(ϕc))**2)
    ρ = a * mc * t / tc
    m = np.cos(ϕ) / np.sqrt(1 - e**2 * np.sin(ϕ)**2)
    x = pm * ρ * np.sin(λ - λ0)
    y = -pm * ρ * np.cos(λ - λ0)
    k = ρ / (a * m)
    return x,y

def regrid(xaxis,yaxis,data,nx,ny,interp_type='spline'):
    r"""
    reterpolates gridded data.
    """
    skipinterp = 0
    if min(xaxis.shape()) > 1:
        matflag = 1
    else:
        matflag = 0

    if len(nx) == 1:
        xstep = xaxis[1]-xaxis[0]
        ystep = yaxis[1]-yaxis[0]
        if nx == 1:
            ystep = ystep*cice/2

        if xstep > ystep:
            nx_temp = np.arange(xaxis[0],axis[-1],ystep)
            ny_temp = yaxis
        elif ystep > xstep:
            nx_temp = xaxis
            ny_temp = np.arange(yaxis[0],yaxis[-1],xstep)
        else:
            skipinterp = 1
            nx_temp = xaxis
            ny_temp = yaxis
        if nx == 1:
            f = ny
            ystep = (1/f)/20;
            ny_temp2 = np.arange(yaxis[0],yaxis[-1],ystep)
            nx_temp2 = np.arange(xaxis[0],ystep*cice/2,xaxis[-1])
            if len(nx_temp2) > len(nx_temp) & len(ny_temp2) > len(ny_temp):
                print('we made it')
            else:
                ny_temp = ny_temp2
                nx_temp = nx_temp2
        nx = nx_temp
        ny = ny_temp
    if skipinterp == 0:
        if matflag == 0:
            y_0, x_0 = np.meshgrid(yaxis,xaxis)
            interp_y, interp_x = np.meshgrid(ny,nx)
        else:
            y_0 = yaxis
            x_0 = xaxis
            interp_y = ny
            interp_x = nx
        A = griddedInterpolant(y_0,x_0,data,interp_type)
        output = A(interp_y,interp_x)
        if min(min(isnan(output))) == 1:
            A = griddedInterpolant(y_0,x_0,data,'cubic')
            output = A(interp_y,interp_x);
    else:
        output = data
    
    return output, nx, ny