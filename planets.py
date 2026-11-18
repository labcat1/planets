import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.patches as ptch
#2 вида блужданий (равнлм и гаус), сделать с 1 млн частиц
nofsteps = 100
m1 = 1.98892e30
m2 = 1
G = 6.6743e-11
x = np.empty((nofsteps, 2))
x[0, 1] =  14959787070.0
v = np.empty((nofsteps, 2))
v[0, 0] = np.sqrt(m1*G/x[0,1])
a = np.empty((nofsteps, 2))
dt = 365*24*60/12.0
size  = x[0, 1]*np.sqrt(2)

xH = np.empty((nofsteps, 2))
xH[0, 1] =  14959787070.0
vH = np.empty((nofsteps, 2))
vH[0, 0] = np.sqrt(m1*G/x[0,1])

xV = np.empty((nofsteps, 2))
xV[0, 1] =  14959787070.0
vV = np.empty((nofsteps, 2))
vV[0, 0] = np.sqrt(m1*G/x[0,1])
xV[1] = xV[0]+vV[0]*dt+np.float64(G)*xV[0]*m1*dt**2/((xV[0,0]*xV[0,0]+xV[0, 1]*xV[0, 1])**(3/2))/2

for n in range(1, nofsteps):
    #эйлер
    x1 = x[n-1] + v[n-1]*dt
    a1 = G*x1*m1/((x1[0]*x1[0]+x1[1]*x1[1])**(3/2))
    v1 = v[n-1] - a1*dt
    #print(v1.dtype)
    a[n] = np.float64(G)*x[n-1]*m1/((x[n-1,0]*x[n-1,0]+x[n-1,1]*x[n-1,1])**(3/2))
    x[n] = x[n-1] + (v1+v[n-1])*dt/2
    v[n] = v[n-1] - (a[n-1]+a1)*dt/2

    #хойн
    ap = np.float64(G)*xH[n-1]*m1/((xH[n-1,0]*xH[n-1,0]+xH[n-1, 1]*xH[n-1, 1])**(3/2))
    xp = xH[n-1] + vH[n-1]*dt
    vp = vH[n-1] + ap*dt
    aH = np.float64(G)*xp*m1/((xp[0]*xp[0]+ xp[1]*xp[1])**(3/2))
    xH[n] = xH[n-1] + (vp + vH[n-1])*dt/2
    vH[n] = vH[n-1] + (ap + aH)*dt/2
    
for n in range(2, nofsteps):
    #верле
    aV = np.float64(G)*xV[n-1]*m1/((xV[n-1,0]*xV[n-1,0]+xV[n-1, 1]*xV[n-1, 1])**(3/2))
    xV[n] = 2*xV[n-1] - xV[n-2]+ aV*dt**2
    vV[n] = (xV[n]-xV[n-1])/dt
    
E = m2*(v[:,0]*v[:,0]+v[:,1]*v[:,1])/2 - np.float64(G)*m1*m2/(np.sqrt(x[:,0]*x[:,0]+x[:, 1]*x[:, 1]))
EH = m2*(vH[:,0]*vH[:,0]+vH[:,1]*vH[:,1])/2 - np.float64(G)*m1*m2/(np.sqrt(xH[:,0]*xH[:,0]+xH[:, 1]*xH[:, 1]))
EV = np.float64(m2)*(vV[:,0]*vV[:,0]+vV[:,1]*vV[:,1])/2 - np.float64(G)*m1*m2/(np.sqrt(xV[:,0]*xV[:,0]+xV[:, 1]*xV[:, 1]))
t = np.arange(0, nofsteps, 1)

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(-size,size)
ax.set_ylim(-size,size)

graphH, = ax2.plot(t, EH)
graphE, = ax2.plot(t, E)
graphV, = ax2.plot(t, EV)
planet1 = ptch.Circle((0, 0), radius=696000000*3)
ax.add_patch(planet1)
planet2 = ptch.Circle(xH[0], radius=696000000)
ax.add_patch(planet2)


def animate(i):
    global xH, x, xV
    planet2.set_center(xH[i])
    return planet2,

#момент имп, эн для эйлераV, э-кромера, хойна и верле
ani = animation.FuncAnimation(
    fig, animate, cache_frame_data=False, interval=10)
plt.show()
