import math
import numpy as np
import csv

m_earth = 5.9722 * 10**24 #earth mass

earth_rad_eq = 6378.137 * 1000 #earth radius at equator 
earth_rad_pole = 6356.752 * 1000 #earth radius at pole
earth_rad_avg = 6371 * 1000

G = 6.6743 * 10**(-11) #gravitational constant
gama = 1.4
mol = 0.0289644 #molar mass of air
R = 8.31446261815324	

def set_initials(h):
    init = [] #P0,T0,L,base
    if(0<=h<11000):
        init.append(101325)
        init.append(288.15)
        init.append(0.0065)
        init.append(0)
    if(11000<=h<20000):
        init.append(22632.064)
        init.append(216.65)
        init.append(0)
        init.append(11000)
    if(20000<=h<32000):
        init.append(5474.88867)
        init.append(216.65)
        init.append(-0.001)
        init.append(20000)
    if(32000<=h<47000):
        init.append(868.018685)
        init.append(228.65)
        init.append(-0.0028)
        init.append(32000)
    if(47000<=h<51000):
        init.append(110.906306)
        init.append(270.65)
        init.append(0)
        init.append(47000)
    if(51000<=h<71000):
        init.append(66.9388731)
        init.append(270.65)
        init.append(0.0028)
        init.append(51000)
    if(71000<=h<84001):
        init.append(3.95642043)
        init.append(214.65)
        init.append(0.002)
        init.append(71000)
    return init

def geopotential_height(h):
    H = (h*earth_rad_avg) / (h+earth_rad_avg)
    return H
    
def temperature(h):
    init = set_initials(h)
    H = geopotential_height(h)
    T = init[1] - init[2] * (H-init[3])
    return T
    
def gravity(h=0,lat=45.5):

    IGF = 9.780327 * (1+0.0053024*math.sin(math.radians(lat))**2 - 0.0000058*2*math.sin(math.radians(lat))**2)
    FAC = -3.086 * 10**(-6) * h
    g = IGF + FAC

    return 9.80665 #g

def grav(h=0,lat=45.5):

    IGF = 9.780327 * (1+0.0053024*math.sin(math.radians(lat))**2 - 0.0000058*2*math.sin(math.radians(lat))**2)
    FAC = -3.086 * 10**(-6) * h
    g = IGF + FAC

    return g

def pressure(h=0):
    init = set_initials(h)
    P0 = init[0]
    T0 = init[1]
    L = init[2]
    b = init[3]
    H = geopotential_height(h)
    if(L != 0):
        P = P0 * ((1-((L*(H-b))/T0))**((gravity(h)*mol)/(R*L)))
    if(L == 0):
        P = P0 * np.exp((-mol*gravity(h)*(H-b))/(R*T0))
    return P
    

def density(h):
    rho = (pressure(h)*mol)/(R*temperature(h))
    return rho

def sound_speed(h):
    a = math.sqrt((gama*R*temperature(h))/mol)
    return a
    
fields = ["altitude","temperature","pressure","density","sound_speed","gravity"]

rows = []

s = [0,1000,3000,5000,10000,25000,50000,75000]
for i in range(0,84000):
    temp = []
    temp.append(str(i))
    temp.append(str(temperature(i)))
    temp.append(str(pressure(i)))
    temp.append(str(density(i)))
    temp.append(str(sound_speed(i)))
    temp.append(str(grav(i)))
    rows.append(temp)
    

fname = "atmosphere_standard.csv"

with open(fname, 'w',newline='') as csvfile: 
    csvwriter = csv.writer(csvfile) 
    csvwriter.writerow(fields)
    csvwriter.writerows(rows)
    
    
    



    
