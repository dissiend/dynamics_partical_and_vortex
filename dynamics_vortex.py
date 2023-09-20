from numpy import*
import matplotlib.pyplot as plt
def increment(f, t, y, h): # Method Runge-Kutta Felberga
    k1 = h * f(t, y)
    k2 = h * f(t + (1 / 4) * h, y + (1 / 4) * k1)
    k3 = h * f(t + (3 / 8) * h, y + (3 / 32) * k1 + (9 / 32) * k2)
    k4 = h * f(t + (12 / 13) * h, y + (1932 / 2197) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3)
    k5 = h * f(t + h, y + (439 / 216) * k1 - 8 * k2 + (3680 / 513) * k3 - (845 / 4104) * k4)
    k6 = h * f(t + (1 / 2) * h, y - (8 / 27) * k1 + 2 * k2 - (3544 / 2565) * k3 + (1859 / 4104) * k4 - (11 / 40) * k5)
    return (16 / 135) * k1 + (6656 / 12825) * k3 + (28561 / 56430) * k4 - (9 / 50) * k5 + (2 / 55) * k6
def res(t,y): # A system of equations describing the dynamics of a vortex three-field and passive particles
    res = zeros(n * 2 + 6)
    res[6] = -(w[0] * (2 * y[7] - 2 * y[1])/ ((y[6] - y[0])**2 + (y[7] - y[1])**2) + w[1] * (2 * y[7] - 2 * y[3])/ ((y[6] - y[2])**2 + (y[7] - y[3])**2) + w[2] * (2 * y[7] - 2 * y[5])/ ((y[6] - y[4])**2 + (y[7] - y[5])**2))
    res[7] = (w[0] * (2 * y[6] - 2 * y[0])/ ((y[6] - y[0])**2 + (y[7] - y[1])**2) + w[1] * (2 * y[6] - 2 * y[2])/ ((y[6] - y[2])**2 + (y[7] - y[3])**2) + w[2] * (2 * y[6] - 2 * y[4])/ ((y[6] - y[4])**2 + (y[7] - y[5])**2))
    res[8] = -(w[0] * (2 * y[9] - 2 * y[1]) / ((y[8] - y[0]) ** 2 + (y[9] - y[1]) ** 2) + w[1] * (2 * y[9] - 2 * y[3]) / ((y[8] - y[2]) ** 2 + (y[9] - y[3]) ** 2) + w[2] * (2 * y[9] - 2 * y[5]) / ((y[8] - y[4]) ** 2 + (y[9] - y[5]) ** 2))
    res[9] = (w[0] * (2 * y[8] - 2 * y[0]) / ((y[8] - y[0]) ** 2 + (y[9] - y[1]) ** 2) + w[1] * (2 * y[8] - 2 * y[2]) / ((y[8] - y[2]) ** 2 + (y[9] - y[3]) ** 2) + w[2] * (2 * y[8] - 2 * y[4]) / ((y[8] - y[4]) ** 2 + (y[9] - y[5]) ** 2))
    res[10] = -(w[0] * (2 * y[11] - 2 * y[1]) / ((y[10] - y[0]) ** 2 + (y[11] - y[1]) ** 2) + w[1] * (2 * y[11] - 2 * y[3]) / ((y[10] - y[2]) ** 2 + (y[11] - y[3]) ** 2) + w[2] * (2 * y[11] - 2 * y[5]) / ((y[10] - y[4]) ** 2 + (y[11] - y[5]) ** 2))
    res[11] = (w[0] * (2 * y[10] - 2 * y[0]) / ((y[10] - y[0]) ** 2 + (y[11] - y[1]) ** 2) + w[1] * (2 * y[10] - 2 * y[2]) / ((y[10] - y[2]) ** 2 + (y[11] - y[3]) ** 2) + w[2] * (2 * y[10] - 2 * y[4]) / ((y[10] - y[4]) ** 2 + (y[11] - y[5]) ** 2))
    res[0] = (-w[0] * w[1] * (2 * y[1] - 2 * y[3]) / ((y[0] - y[2]) ** 2 + (y[1] - y[3]) ** 2) - w[0] * w[2] * (2 * y[1] - 2 * y[5]) / ((y[0] - y[4]) ** 2 + (y[1] - y[5]) ** 2)) / w[0]
    res[1] = -(-w[0] * w[1] * (2 * y[0] - 2 * y[2]) / ((y[0] - y[2]) ** 2 + (y[1] - y[3]) ** 2) - w[0] * w[2] * (2 * y[0] - 2 * y[4]) / ((y[0] - y[4]) ** 2 + (y[1] - y[5]) ** 2)) / w[0]
    res[2] = (-w[0] * w[1] * (-2 * y[1] + 2 * y[3]) / ((y[0] - y[2]) ** 2 + (y[1] - y[3]) ** 2) - w[1] * w[2] * (2 * y[3] - 2 * y[5]) / ((y[2] - y[4]) ** 2 + (y[3] - y[5]) ** 2)) / w[1]
    res[3] = -(-w[0] * w[1] * (-2 * y[0] + 2 * y[2]) / ((y[0] - y[2]) ** 2 + (y[1] - y[3]) ** 2) - w[1] * w[2] * (2 * y[2] - 2 * y[4]) / ((y[2] - y[4]) ** 2 + (y[3] - y[5]) ** 2)) / w[1]
    res[4] = (-w[0] * w[2] * (-2 * y[1] + 2 * y[5]) / ((y[0] - y[4]) ** 2 + (y[1] - y[5]) ** 2) - w[1] * w[2] * (-2 * y[3] + 2 * y[5]) / ((y[2] - y[4]) ** 2 + (y[3] - y[5]) ** 2)) / w[2]
    res[5] = -(-w[0] * w[2] * (-2 * y[0] + 2 * y[4]) / ((y[0] - y[4]) ** 2 + (y[1] - y[5]) ** 2) - w[1] * w[2] * (-2 * y[2] + 2 * y[4]) / ((y[2] - y[4]) ** 2 + (y[3] - y[5]) ** 2)) / w[2]
    return res
# the beginning of the countdown
t = 0
# intensity
w = [2,-1,-1]
n=3
# the end of the countdown
T =10
# initial conditions
d=1
e=0.4
# central vortex
xc=0.4
yc=e
# satellites
xstart1=0
ystart1=-d
xstart2=0
ystart2=d
# passive particles
xpart1 = -2.1
xpart2 = -2.3
xpart3 = -1.9
ypart = zeros(3)
ypart[0]=-1
ypart[1]=0
ypart[2]=1
yo = array([xc, yc, xstart1, ystart1, xstart2, ystart2, xpart1, ypart[0], xpart2, ypart[1], xpart3, ypart[2]])
h = 0.002  # step
yall = []
yall.append(yo)
while t < T:
    yo = yo + increment(res, t, yo, h)
    t = t + h
    yall.append(yo)

# Distribution of vortex values
x1 = array([i[0] for i in yall])
y1 = array([i[1] for i in yall])
x2 = array([i[2] for i in yall])
y2 = array([i[3] for i in yall])
x3 = array([i[4] for i in yall])
y3 = array([i[5] for i in yall])

# Distribution of particle values

xp1 = array([i[6] for i in yall])
yp1 = array([i[7] for i in yall])
xp2 = array([i[8] for i in yall])
yp2 = array([i[9] for i in yall])
xp3 = array([i[10] for i in yall])
yp3 = array([i[11] for i in yall])

# Visualization

plt.title('T=%i' %T)
plt.plot(x1,y1,'g',linewidth=0.5)
plt.plot(x1[len(x1) - 1],y1[len(y1) - 1],'gs')
plt.plot(x2,y2,'b',linewidth=0.5)
plt.plot(x2[len(x2) - 1],y2[len(y2) - 1],'b^')
plt.plot(x3,y3,'r',linewidth=0.5)
plt.plot(x3[len(x3) - 1],y3[len(y3) - 1],'r^')
plt.plot(xp1,yp1,'k',linewidth=1)
plt.plot(xp1[0],yp1[0],'ko')
plt.plot(xp1[len(xp1) - 1],yp1[len(yp1) - 1],'k*')
plt.plot(xp2,yp2,'k',linewidth=0.5)
plt.plot(xp2[0],yp2[0],'ko')
plt.plot(xp2[len(xp2) - 1],yp2[len(yp2) - 1],'k*')
plt.plot(xp3,yp3,'k',linewidth=0.3)
plt.plot(xp3[0],yp3[0],'ko')
plt.plot(xp3[len(xp3) - 1],yp3[len(yp3) - 1],'k*')
plt.axis('equal')
plt.show()