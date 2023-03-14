import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

p0 = 1 # *Па* = кг/(м*мс^2)
'''
ro1 = 2000 # кг/м^3
ro2 = 3000 # без отражения - работает!
'''

ro1 = 3000 # кг/м^3
ro2 = 2000 #

c1 = 3 # м/мс 
c2 = 2
L = 250 #м
L0 = 150 #м
n0 = 100
h = L0 / n0 # ~1/n0
N = int(L / h) #кол-во шагов по x ~n0
l = L0/2 #м длина волны
q = c1 / c2
y = 0.5 * (q-1)/(q+1) 
CFL2 = 0.5 - y
CFL1 = 0.5 + y
'''
CFL1 = 0.4
CFL2 = 0.4
'''
T = 50 #мс
tau = CFL1/c1*h #мс ~1/n0 
M = int(T/tau) + 1 #! #кол-во шагов по t ~n0
Z1 = ro1*c1
Z2 = ro2*c2
K1=Z1*c1
K2=Z2*c2

def p(x):
    if x <= l:
        #return 1
        return p0*np.sin(np.pi*x/l)**4
    else:
        return 0
    
def p_der(x):
    if x <= l:
        return (4*np.pi/l/ro1)*p0*(np.sin(np.pi*x/l)**3)*(np.cos(np.pi*x/l))
    else:
        return 0

p_curr = np.array([0] * N, 'float64')
for i in range(N):
    p_curr[i] = p(i*h)
p_next = np.array([0] * N, 'float64')

p_curr_lin = p_curr.copy()
p_next_lin = np.array([0] * N, 'float64')

v_curr = p_curr / Z1
v_next = np.array([0] * N, 'float64')

v_curr_lin = v_curr.copy()
v_next_lin = np.array([0] * N, 'float64')

p_der_curr = np.array([0] * N, 'float64')
for i in range(N):
    p_der_curr[i] = p_der(i*h)
p_der_next = np.array([0] * N, 'float64')
v_der_curr = Z1 * p_der_curr
v_der_next = np.array([0] * N, 'float64')

p_init = p_curr.copy()
v_init = v_curr.copy()                 

def min_max(a, b, c, d, x):
    y = a*x**3 + b*x**2 + c*x + d
    yl = -a*h**3 + b*h**2 + -c*h + d
    yr = d
    if(y>max(yl, yr)):
        return max(yl, yr)
    if(y<min(yl, yr)):
        return min(yl, yr)
    return y

def Transport(c_l, ro_l, c_r, ro_r, CFL_l, CFL_r):
    Zl = c_l*ro_l
    Kl = Zl*c_l
    Zr = c_r*ro_r
    Kr = Zr*c_r
    
    PosWl = 1/(2*Zl) * (p_curr[j-1] + v_curr[j-1]*Zl)
    PosW0 = 1/(2*Zl) * (p_curr[j] + v_curr[j]*Zl)
    NegW0 = 1/(2*Zr) * (-p_curr[j] + v_curr[j]*Zr)
    NegWr = 1/(2*Zr) * (-p_curr[j+1] + v_curr[j+1]*Zr)
    PosVl = 1/2 * (p_der_curr[j-1]/c_l + v_der_curr[j-1]/Kl)
    PosV0 = 1/2 * (p_der_curr[j]/c_l + v_der_curr[j]/Kl)
    NegV0 = 1/2 * (-p_der_curr[j]/c_r + v_der_curr[j]/Kr)
    NegVr = 1/2 * (-p_der_curr[j+1]/c_r + v_der_curr[j+1]/Kr)
    
    PosWl_lin = 1/(2*Zl) * (p_curr_lin[j-1] + v_curr_lin[j-1]*Zl)
    PosW0_lin = 1/(2*Zl) * (p_curr_lin[j] + v_curr_lin[j]*Zl)
    NegW0_lin = 1/(2*Zr) * (-p_curr_lin[j] + v_curr_lin[j]*Zr)
    NegWr_lin = 1/(2*Zr) * (-p_curr_lin[j+1] + v_curr_lin[j+1]*Zr)
    
    al = 2 * (PosWl + PosV0*h - PosW0)/h**3 + (PosVl-PosV0)/h**2
    bl = 3 * (PosWl + PosV0*h - PosW0)/h**2 + (PosVl-PosV0)/h
    cl = PosV0
    dl = PosW0
    ar = 2 * (NegW0 + NegVr*h - NegWr)/h**3 + (NegV0-NegVr)/h**2
    br = 3 * (NegW0 + NegVr*h - NegWr)/h**2 + (NegV0-NegVr)/h
    cr = NegVr
    dr = NegWr
    
    xl = -c_l*tau
    xr = -h + c_r*tau
    
    PosW = min_max(al, bl, cl, dl, xl)
    NegW = min_max(ar, br, cr, dr, xr)
    '''
    PosW = al*xl**3 + bl*xl**2 + cl*xl + dl
    NegW = ar*xr**3 + br*xr**2 + cr*xr + dr
    '''
    PosV = 3*al*xl**2 + 2*bl*xl + cl
    NegV = 3*ar*xr**2 + 2*br*xr + cr
    
    PosW_lin = PosW0_lin - CFL_l*(PosW0_lin - PosWl_lin)
    NegW_lin = NegW0_lin + CFL_r*(NegWr_lin - NegW0_lin)
    
    p_next_lin[j] = 2 * Zl*Zr/(Zl + Zr) * (PosW_lin - NegW_lin)
    v_next_lin[j] = 2/(Zl + Zr) * (Zl*PosW_lin + Zr*NegW_lin)
    
    p_next[j] = 2 * Zl*Zr/(Zl + Zr) * (PosW - NegW)
    v_next[j] = 2/(Zl + Zr) * (Zl*PosW + Zr*NegW)
    p_der_next[j] = 2*(Kl*PosV - Kr*NegV) / (Zl + Zr)
    v_der_next[j] = 2*Zl*Zr*(c_l*PosV + c_r*NegV) / (Zl + Zr)

for i in range(M):
    for j in range(1, n0):
        Transport(c1, ro1, c1, ro1, CFL1, CFL1)
    j = n0
    Transport(c1, ro1, c2, ro2, CFL1, CFL2)
    for j in range(n0+1, N-1):
        Transport(c2, ro2, c2, ro2, CFL2, CFL2)
    p_curr, p_next = p_next, p_curr
    v_curr, v_next = v_next, v_curr
    p_der_curr, p_der_next = p_der_next, p_der_curr
    v_der_curr, v_der_next = v_der_next, v_der_curr
    p_curr_lin, p_next_lin = p_next_lin, p_curr_lin
    v_curr_lin, v_next_lin = v_next_lin, v_curr_lin

#графики
x_ = np.linspace(0, L, N)
fs = 10 # fontsize
fig, ax = plt.subplots()
ax.plot(x_, p_init, linestyle='-', marker='', color='b', linewidth=1, label='t = 0 мс')
ax.plot(x_, p_curr, linestyle='-', marker='', color='r', linewidth=1, label='t = '+str(T)+' мс кубическое приближение')
ax.plot(x_, p_curr_lin, linestyle='--', marker='', color='g', linewidth=1, label='t = '+str(T)+' мс линейное приближение')
l_ = mlines.Line2D([L0, L0], [p_curr.min(), p0])
ax.add_line(l_)
plt.title('p(x,t)')
plt.xlabel('x, м')
plt.ylabel('p, 10^6 Па')
plt.legend()
plt.text(0.25*L0, 0.3*p0, 'с = '+str(c1*10**3)+' м/с', fontsize=fs)
plt.text(0.25*L0, 0.2*p0, 'ρ = '+str(ro1)+' кг/м^3', fontsize=fs)
plt.text(0.25*L0, 0.1*p0, 'CFL = '+str(round(CFL1, 2)), fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.3*p0, 'с = '+str(c2*10**3)+' м/с', fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.2*p0, 'ρ = '+str(ro2)+' кг/м^3', fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.1*p0, 'CFL = '+str(round(CFL2, 2)), fontsize=fs)
plt.grid()

fig, ax = plt.subplots()
ax.plot(x_, 10**3 * v_init, linestyle='-', marker='', color='b', linewidth=1, label='t = 0 мс')
ax.plot(x_, 10**3 * v_curr, linestyle='-', marker='', color='r', linewidth=1, label='t = '+str(T)+' мс кубическое приближение')
ax.plot(x_, 10**3 * v_curr_lin, linestyle='--', marker='', color='g', linewidth=1, label='t = '+str(T)+' мс линейное приближение')
l_ = mlines.Line2D([L0, L0], [0, v_curr.max()*10**3])
ax.add_line(l_)
plt.title('v(x,t)')
plt.xlabel('x, м')
plt.ylabel('v, м/c')
plt.legend()
plt.text(0.25*L0, 0.3*p0/Z1*10**3, 'с = '+str(c1*10**3)+' м/с', fontsize=fs)
plt.text(0.25*L0, 0.2*p0/Z1*10**3, 'ρ = '+str(ro1)+' кг/м^3', fontsize=fs)
plt.text(0.25*L0, 0.1*p0/Z1*10**3, 'CFL = '+str(round(CFL1, 2)), fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.3*p0/Z1*10**3, 'с = '+str(c2*10**3)+' м/с', fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.2*p0/Z1*10**3, 'ρ = '+str(ro2)+' кг/м^3', fontsize=fs)
plt.text(L0 + 0.25*(L - L0), 0.1*p0/Z1*10**3, 'CFL = '+str(round(CFL2, 2)), fontsize=fs)
plt.grid()
