# boussinesq function for line load

from math import pi
import numpy as np
import matplotlib.pyplot as plt
import math

def pll(q:float,z:float,x:float)->float:
    """ calculate the pressure a z, x below the foundation """
    """ for a line load """
    a = 2*q*z**3
    b=pi*(x**2+z**2)**2
    deltatP = a/b
    return deltatP

def strip_pressure(p:float, b:float, x:float, z:float)->float:
    """ calculate the pressure at a point x, z below the foundation """

    # calculate the angle from the edge of the foundation to the point being considered.
    x1 = b/2 + x
    x2 = x - b/2
    r1 = np.sqrt(x1**2 + z**2)
    r2 = np.sqrt(x2**2 + z**2)
    # calculate the pressure at a point
    theta1 = np.arctan(x1/z)
    theta2 = np.arctan(x2/z)
    s  = (p/pi) * ((theta1 - theta2) + np.sin(theta1)*np.cos(theta1) - np.sin(theta2)*np.cos(theta2))

    return s

def rect_pressure(q:float,b:float,l:float,z:float)->float:
    m = b/(2*z) # divided to split into rectangles
    n = l/(2*z) # divided to split into rectangles
    a = (2*m*n*math.sqrt(m**2+n**2+1))/(m**2+n**2+m**2 * n**2 +1)
    b = (m**2+n**2+2)/(m**2+n**2+1)
    c = np.arctan((2*m*n*math.sqrt(m**2+n**2+1))/(m**2+n**2-m**2 * n**2 +1))
    I2 = 1/(4*math.pi) * (a * b + c)
    return 4*I2 * q # pressure at the center of the rectangle

def plot_strip_pressure(q:float, b:float, x_max:float, z_max:float)->None:
    x = np.linspace (0, x_max, 100)
    z = np.linspace (0, -z_max, 100)
    X, Z = np.meshgrid(x, z)
    p = strip_pressure(q, b, X, Z)

    # plot the results
    intervals= np.linspace(-q, -5, 20)
    plt.contourf(x, z, p, intervals)
    plt.colorbar()
    plt.xlabel('x distance (m)')
    plt.ylabel('y distance (m)')
    plt.title('Half space plot of pressure (kPa) below a \nuniform strip contact pressure',fontsize=12)
    plt.grid()
    plt.show()

def plot_rect(q:float,b:float,l:float,z:float, g:f)->None:
    y = np.linspace(-.2,-z,1000)
    p = np.array([])
    ob = np.array([])
    for i in y:
        press = rect_pressure(q,b,l,i)
        p = np.append(p,press)
        ob = np.append(ob, -i*g)
    ratio = p/ob
    plt.figure(figsize=(10,6))
    plt.plot(p,y, label='Vertical Stress increase')
    plt.plot(ob,y,label='Overburden Pressure', ls='dotted')
    plt.plot(20/100*ob,y,label='20% Overburden Pressure', ls='dotted')
    plt.plot(10/100*ob,y,label='10% Overburden Pressure', ls='dotted')
    plt.plot(5/100*ob,y,label='5% Overburden Pressure', ls='dotted')
    plt.fill_between(p,y,0, alpha=0.20)
    plt.grid()
    plt.title('Vertical Pressure below the center of a rectangular flexible foundation')
    plt.xlabel('Pressure (kPa)')
    plt.ylabel('Depth (m)')
    plt.legend()
    plt.savefig('temp3.png')
    plt.show()
