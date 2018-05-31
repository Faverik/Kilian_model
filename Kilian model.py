import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

"""
This program is designed to describe a certain set of quantities
(the area of particles) by the statistical model of Kilian.
To do this, create a data file containing a list of these values
and put the path to the file in "data_file" variable, using the \\ instead \
In case of wrong fitting try to change start params(a1-u2)
In fact, it does not work well with three or more ensembles,
as this requires fine-tuning the fitting bounds
In real systems for which it is created, more than two ensembles are not used
"""

# Parameters
ensemble = 2  # number of ensembles
s_bounds = ([-10, -7, -12, -7], [-7.5, -3.3, -7, -4])  # bounds for fitting,
# in case of more than 2 ensembles other bounds will be calculate
# s_bounds[0, 0] - lower bound of amplitude 1-st ensemble
# s_bounds[0, 1] - lower bound of  aggregation energy 1-st ensemble (rise-graph move down to the left)
# s_bounds[0, 2] - lower bound of  amplitude 2-nd ensemble
# s_bounds[0, 3] - lower bound of  aggregation energy 2-st ensemble
# s_bounds[1, :] - upper bounds
freq_count = 20  # number of columns in the histogram
pix_in_mkm = 1.826  # real pixel area on the photo
max_s = 25000  # maximum particle area in mkm**2
step = 10000  # smooth graphics
colon_space = 10000  # space between columns on hist
path = "C:\\Users\\nefae\\PycharmProjects\\Kilian_model\\"
data_file = "{0}data.txt".format(path)

# Constants
YB1 = 0  # distribution zero-point (common 0)
YB2 = 0  # distribution zero-point (common 0)
K = 2  # mean 2D (photo dimension)

# Declaration
cells = []
plot_y = []
text_p = ""
# Input data
f = open(data_file, "r")
for line in f.readlines():  # file reading in cells variable
    if line != "'\n'":
        cells.append(float(line[0:-1]))
f.close()
# Data preparation
number, area, dump = plt.hist(cells, freq_count, color="Gray")  # histogram calculation function
norm_number = number / sum(number)  # normalization
norm_area = area[0:-1] * pix_in_mkm  # targets area in mkm**2
colon = max_s / freq_count - colon_space  # width of columns on histogram
plt.close()
plt.show()
plt.figure()
plt.bar(norm_area, norm_number, width=colon, color="Yellow")  #
# print("number={0}, area={1}".format(np.round(norm_number, 3), np.round(norm_area, -2)))  # checking

# Fitting=============================================================
# Creation of additional initiating values
s_params = [(s_bounds[1][0]-s_bounds[0][0])/2+s_bounds[0][0],
            (s_bounds[1][1]-s_bounds[0][1])/2+s_bounds[0][1],
            (s_bounds[1][2]-s_bounds[0][2])/2+s_bounds[0][2],
            (s_bounds[1][3]-s_bounds[0][3])/2+s_bounds[0][3]]


def param_create(p, koef):  # make every next parameter of ensemble smaller
    return p - koef - 1  # each next parameter should be smaller than previous

if ensemble > 2:
    for m in range(ensemble-2):
        print(m)
        s_params.append(param_create(s_params[2], m))  # adding new start-point
        s_params.append(param_create(s_params[3], m))
        s_bounds[0].append(param_create(s_bounds[0][2], m))  # adding new bound
        s_bounds[0].append(param_create(s_bounds[0][3], m))
        s_bounds[1].append(param_create(s_bounds[1][2], m))
        s_bounds[1].append(param_create(s_bounds[1][3], m))

# Creating a sum fitting function with a variable number of members
def sum_func(x1_2_, i, args):
    if i <= -1:
        return 0
    else:
        return (1*10**args[i-1] * (x1_2_ - YB1) ** K) * (np.exp(-(x1_2_ - YB1) * 1*10**args[i])) + sum_func(x1_2_, (i - 2), args[:])


def func(x1_2_, *args):
    i = ensemble*2-1
    return sum_func(x1_2_, i, args[:]) + sum_func(x1_2_, (i - 2), args[:])

# Fittings process
try:
    popt, pcov = curve_fit(func, norm_area, norm_number, p0=s_params, bounds=s_bounds)
except:
    print("Perhaps you went out the bounds, check it")

# =====================================================================
# Plot
plot_x = np.linspace(0, max(norm_area), step)
plot_y.append(sum_func(plot_x, ensemble*2-1, popt))  # creates data for the entire function
print("params={0}".format(np.round(popt, 10)))
if ensemble > 1:
    for z in range(ensemble):
        k = 2 * z
        print("params={0}".format(np.round(popt[k:k+2], 10)))  # creates data for individual function members
        plot_y.append(sum_func(plot_x, 1, popt[k:k+2]))


color_t = ['k', 'r', 'b', 'y', 'g', 'm', 'c']

for t in range(len(plot_y)):
    plt.plot(plot_x, plot_y[t], color=color_t[t], label='{0}-ensemble'.format(t))
plt.legend(loc='upper right')
plt.xlabel(u'${\mu}m$²', fontsize=14)
plt.ylabel('share')
x_pos = max(plot_x) - (max(plot_x)/5)
y_pos = max(plot_y[0]) - 0.1
l = 1
for q in range(0, ensemble+1, 2):
    text_p += 'a{0}=1*10^{1} \nu{2}=1*10^{3} \n'.format(l, np.round(popt[q], 2), l, np.round(popt[q+1], 2))
    l += 1

plt.text(x_pos, y_pos, text_p)
plt.show()

l = 1
res_file = "{0}result.txt".format(path)
f = open(res_file, 'w')
for i in range(0, ensemble+1, 2):
    f.write("a{0}={1}\n".format(l, 1*10**popt[i]))
    f.write("u{0}={1}\n".format(l, 1*10**popt[i+1]))
    l += 1
f.close()
