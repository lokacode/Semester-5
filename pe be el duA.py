import numpy as np
from matplotlib import pyplot as plt

def qSimpson(Func, a, b, n):
    if (n % 2 == 0): n += 1 # increment n if even
    h = (b-a)/(n-1)
    s1 = s2 = 0e0
    for i in range(2,n-2,2): s1 += Func(a + i*h) # odd-index sum
    for i in range(1,n-1,2): s2 += Func(a + i*h) # even-index sum
    return (h/3)*(Func(a) + 4*s2 + 2*s1 + Func(b))

def percepatan(t):
    return 1.625

print(qSimpson(percepatan, 0, 10, 1000))

# area under the curve
x = np.linspace(-2, 12, 100)
x1 = np.linspace(0, 10, 100)
y = [percepatan(t) for t in x]
y1 = [percepatan(t) for t in x1]

# plot the function
plt.plot(x,y)
plt.fill_between(x1,y1,0,facecolor='cyan',alpha=0.5)
plt.xlabel('$t$ (s)')
plt.ylabel('$g$ (m/s$^2$)')
plt.grid(axis='x')
# x and y axis line
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
# remove down and left ticks
plt.gca().spines['bottom'].set_position(('data',0))
plt.gca().spines['left'].set_position(('data',0))
# remove up and right ticks
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("pe be el duA.png", dpi=300)
plt.show()