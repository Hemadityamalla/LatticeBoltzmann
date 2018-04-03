# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 18:04:37 2018

@author: hemaditya
"""
import numpy as np
import matplotlib.pyplot as plt


def setBC(pop):
    for x in range(1,Nx+1):
        xm = x-1
        xp = x+1
        pop[0][x][2] = np.copy(pop[1][x][4])
        pop[0][xm][5] = np.copy(pop[1][x][7])
        pop[0][xp][6] = np.copy(pop[1][x][8])
        
        pop[Ny+1][x][4] = np.copy(pop[Ny][x][2])
        pop[Ny+1][xm][8] = np.copy(pop[Ny][x][6])
        pop[Ny+1][xp][7] = np.copy(pop[Ny][x][5])
        
    for y in range(1,Ny+1):
        pop[y][0] = np.copy(pop[y][Nx])
        pop[y][Nx+1] = np.copy(pop[y][1])
    
    return pop
            
def displace(pop):
    popTemp = np.zeros([Ny+2,Nx+2,9])
    for y in range(1,Ny+1):
        ym = y-1
        yp = y+1
        for x in range(1,Nx+1):
            xm = x-1
            xp = x+1
            
            popTemp[y][x][0] = np.copy(pop[y][x][0])
            popTemp[y][x][1] = np.copy(pop[y][xm][1])
            popTemp[y][x][5] = np.copy(pop[ym][xm][5])
            popTemp[y][x][2] = np.copy(pop[ym][x][2])
            popTemp[y][x][6] = np.copy(pop[ym][xp][6])
            popTemp[y][x][3] = np.copy(pop[y][xp][3])
            popTemp[y][x][7] = np.copy(pop[yp][xp][7])
            popTemp[y][x][4] = np.copy(pop[yp][x][4])
            popTemp[y][x][8] = np.copy(pop[yp][xm][8])
    for x in range(1,Nx+1):
        for y in range(1,Ny+1):
            for i in range(9):
                pop[y][x][i] = np.copy(popTemp[y][x][i])
    
    return pop
    
    
    
    
def collide(pop):
    feq = np.zeros(9)
    
    for y in range(1,Ny+1):
        for x in range(1,Nx+1):
            rho = np.sum(pop[y][x])
            u = np.dot(c[:,0],pop[y][x])/rho + tau*5e-5
            v = np.dot(c[:,1],pop[y][x])/rho
            usq = u**2
            vsq = v**2
            
            sumsq = (usq + vsq)/(2*csqr)
            sumsq2 = sumsq*(1.0 - csqr)/csqr
            u2 = usq/(2*csqr**2)
            v2 = vsq/(2*csqr**2)
            
            ui = u/csqr
            vi = v/csqr
            uv = ui*vi
            
            feq[0] = rho*w[0]*(1-sumsq)
            feq[1] = rho*w[1]*(1-sumsq + u2 + ui)
            feq[2] = rho*w[2]*(1-sumsq + v2 + vi) 
            feq[3] = rho*w[3]*(1-sumsq + u2 - ui)
            feq[4] = rho*w[4]*(1-sumsq + v2 - vi)
            feq[5] = rho*w[5]*(1+sumsq + ui + vi + uv)
            feq[6] = rho*w[6]*(1+sumsq - ui + vi - uv)
            feq[7] = rho*w[7]*(1+sumsq - ui - vi + uv)
            feq[8] = rho*w[8]*(1+sumsq + ui - vi - uv)
            
            for i in range(8):
                pop[y][x][i] = np.copy(pop[y][x][i] - (1.0/tau)*(pop[y][x][i] - feq[i]))
    
    return pop

                
Nx = 20
Ny = 50
t0 = 0
tFinal = 25
timeStep = 0
delt = 1


#D2Q9 parameters
q = 9
w = np.array([4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36])
c = np.array([np.array([0,0]),np.array([1,0]),np.array([0,1]),np.array([-1,0]),np.array([0,-1]),np.array([1,1]),np.array([-1,1]),np.array([-1,-1]),np.array([1,-1])])

#Initialization- extra cells have been initialized as ghost cells
f = np.zeros([Ny+2,Nx+2,q])
feq = np.zeros([Ny+2,Nx+2,q])
fNew = np.zeros([Ny+2,Nx+2,q])
tau = 0.666
csqr = 1./3#Speed of sound
nu = csqr*(tau - 0.5)#Kinematic viscosity

velocity = np.zeros([Ny+2,Nx+2,2])
vold = np.zeros([Ny+2,Nx+2,2])

fInit = np.zeros(9)
for i in range(8):
    fInit[i] = 1.0/9
for y in range(Ny+2):
    for x in range(Nx+2):
        f[y][x] = np.copy(fInit)
        
Max_Step = 10000
for time in range(Max_Step):
    f = np.copy(setBC(f))
    for y in range(1,Ny+1):
        for x in range(1,Nx+1):
            velocity[y][x][0] = np.dot(c[:,0],f[y][x])
            velocity[y][x][1] = np.dot(c[:,1],f[y][x])
            
    if (time%10 == 0):
        error = 0.0
        for y in range(1,Ny+1):
            for x in range(1,Nx+1):
                tmpx = velocity[y][x][0] - vold[y][x][0]
                tmpy = velocity[y][x][1] - vold[y][x][1]
                error = error + (tmpx**2 + tmpy**2)
        
        print('Error: ',error)
        if (error < 1e-11) and (time != 0):
            print ('Run Terminated')
            
    mass = 0.0
    for y in range(1,Ny+1):
        for x in range(1,Nx+1):
            mass = mass + np.sum(f[y][x])
            
    f = np.copy(displace(f))
    f = np.copy(setBC(f))
    f = np.copy(collide(f))
    
    for y in range(1,Ny+1):
        for x in range(1,Nx+1):
            vold[y][x][0] = np.copy(velocity[y][x][0])
            vold[y][x][1] = np.copy(velocity[y][x][1])