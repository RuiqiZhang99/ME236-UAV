import numpy as np
from pylab import *
import sympy as sy 
import seaborn as sns

sns.set(style='darkgrid', font='Times New Roman')

command = np.array([40, 80, 120, 160, 200, 240])
max_speed = np.array([12059, 16707, 21788, 20004, 23937, 27308]) * 2*np.pi/60
min_speed = np.array([9063, 13473, 16613, 19800, 22865, 25591]) * 2*np.pi/60
final_speed = np.array([10666, 13672, 16637, 19856, 22891, 25662]) * 2*np.pi/60
scatter(command, max_speed, color=sns.xkcd_rgb['windows blue'])
scatter(command, min_speed, color=sns.xkcd_rgb['amber'])
scatter(command, final_speed, color=sns.xkcd_rgb['grey'])
[a_max, b_max] = np.polyfit(command, max_speed, 1)
[a_min, b_min] = np.polyfit(command, min_speed, 1)
[a_fin, b_fin] = np.polyfit(command, final_speed, 1)
x=np.linspace(0, 300)
plt.plot(x, a_max*x+b_max, color=sns.xkcd_rgb['windows blue'], label='max speed')
plt.plot(x, a_min*x+b_min, color=sns.xkcd_rgb['amber'], label='min speed')
plt.plot(x, a_fin*x+b_fin, color=sns.xkcd_rgb['grey'], label='final speed')
text(200, 2e3,'y={:.2f}x+{:.2f}'.format(a_fin, b_fin))
xlabel('PWM commands')
ylabel('real motor speed (rad/s)')
title('PWM - Moter Speed')
plt.legend()
plt.show()