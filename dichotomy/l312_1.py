import matplotlib.pyplot as plt
from tabulate import tabulate
from py_expression_eval import Parser
import timeit

class Function:
    def __init__(self, s):
        self.parser = Parser()
        self.s = s
        self.expr = self.parser.parse(s)

    def eval(self, x):
        return self.expr.evaluate({'x': x})


class FunctionPlotter:
    def __init__(self, f, a, b, step=0.01):
        self.func = f
        self.f = f.eval
        self.a = a
        self.b = b
        self.step = step

    def plot_function(self):
        xs = []
        ys = []
        x = self.a
        xs.append(x)
        ys.append(self.f(x))

        while x < self.b:
            x += self.step
            xs.append(x)
            ys.append(self.f(x))

        plt.plot(xs, ys, 'b', label='function')
        plt.ylabel('f(x)')
        plt.xlabel('x')
        plt.suptitle(self.func.s)

    def plot_minimum(self, min):
        plt.plot(min, self.f(min), 'D', label='minimum')

    def plot_maximum(self, max):
        plt.plot(max, self.f(max), 'D', label='maximum')

    def legend(self):
        plt.legend()

    def show(self):
        plt.show()

    def clear(self):
        plt.clf()


def dichotomy_max(f, a, b, eps, n, table_max):  # выборка где большее число
    start_time = timeit.default_timer()
    step = 0
    f_isp = 1
    table_max[step] = [step, a, b, f(a), f(b)]

    c = (a + b) / 2.0
    x1 = (a + c) / 2.0
    x2 = (c + b) / 2.0
    valc = f(c)

    while b - a > eps and step < n:
        step += 1
        f_isp += 2
        val1 = f(x1)
        val2 = f(x2)

        if val1 > valc and valc >= val2:  # f(x1) > F(c) >= f(x2)
            b = c
            c = x1
            valc = val1
        else:
            if val1 < valc and valc >= val2:  # f(x1) < F(c) >= f(x2)
                a = x1
                b = x2
                c = (a + b) / 2.0
            else:
                if val1 > valc and valc < val2 and val1 >= val2:
                    b = c
                    c = x1
                    valc = val1
                else:
                    a = c
                    c = x2
                    valc = val2
        table_max[step] = [step, a, b, f(a), f(b)]
        x1 = (a + c) / 2.0
        x2 = (c + b) / 2.0

    return [(b + a) / 2.0, step, timeit.default_timer() - start_time,f_isp]


def dichotomy_min(f, a, b, eps, n, table_min):
    start_time = timeit.default_timer()
    step = 0
    f_isp = 1
    table_min[step] = [step, a, b, f(a), f(b)]

    c = (a + b) / 2.0
    x1 = (a + c) / 2.0
    x2 = (c + b) / 2.0
    valc = f(c)

    while b - a > eps and step < n:
        step += 1
        f_isp += 2
        val1 = f(x1)
        val2 = f(x2)

        if val1 <= valc and valc < val2:
            b = c
            c = x1
            valc = val1
        else:
            if val1 > valc and valc <= val2:
                a = x1
                b = x2
                c = (a + b) / 2.0
            else:
                if val1 <= valc and valc > val2 and val1 <= val2:
                    b = c
                    c = x1
                    valc = val1
                else:
                    a = c
                    c = x2
                    valc = val2
        table_min[step] = [step, a, b, f(a), f(b)]
        x1 = (a + c) / 2.0
        x2 = (c + b) / 2.0

    return [(b + a) / 2.0, step, timeit.default_timer() - start_time,f_isp]

def menu():
    vf = int(input('f:1 simple 2     '))
    vg = int(input('figg 1=max 2=min  '))
    if vf == 1:
        var1 = ('sqrt(36-sin(x)^2)*sin(x)')
    else:
        var1 = ('((6-x)*x)')
    print(var1)
    f = Function(var1)
    a = float(input('a low  '))
    b = float(input('b up  '))
    max_steps = int(input('[max step]    n = '))
    eps = float(0.01)
    print(eps)

    table_min = {0: [0, 0, 0, 0, 0]}
    table_max = {0: [0, 0, 0, 0, 0]}
    min = []
    max = []
    min = dichotomy_min(f.eval, a, b, eps, max_steps, table_min)
    max = dichotomy_max(f.eval, a, b, eps, max_steps, table_max)
    print('find min fxy', f.eval(min[0]), 'dot', min[0], '_', min[1], 'step', 'time',min[2],'all df',min[3])
    print('find max fxy', f.eval(max[0]), 'dot', max[0], '_', max[1], 'step', 'time',max[2],'all df',max[3])

    print(tabulate(list(table_max.values()), headers=['step', 'a', 'b', 'f(a)', 'f(b)']))
    print(tabulate(list(table_min.values()), headers=['step', 'a', 'b', 'f(a)', 'f(b)']))

    plot = FunctionPlotter(f, a, b)
    plot.clear()
    plot.plot_function()
    plot.plot_minimum(min[0])
    plot.plot_maximum(max[0])

    if vg == 2:
        for i in range(table_min.__len__()-1):
            plt.vlines(table_min[i][1], table_min[i][3], table_min[i+1][3], linestyle="dashed")
            plt.vlines(table_min[i][2], table_min[i][3], table_min[i+1][3], linestyle="dashed")
            plt.vlines(table_min[i][1], table_min[i][4], table_min[i+1][4], linestyle="dashed")
            plt.vlines(table_min[i][2], table_min[i][4], table_min[i+1][4], linestyle="dashed")
            plt.hlines(table_min[i][3], table_min[i][1], table_min[i][2], linestyle="dashed")
            plt.hlines(table_min[i][4], table_min[i][1], table_min[i][2], linestyle="dashed")
            plt.hlines(table_min[i+1][3], table_min[i][1], table_min[i][2], linestyle="dashed")
            plt.hlines(table_min[i+1][4], table_min[i][1], table_min[i][2], linestyle="dashed")
    if vg == 1:
        for i in range(table_max.__len__()-1):
            plt.vlines(table_max[i][1], table_max[i][3], table_max[i+1][3], linestyle="dashed")
            plt.vlines(table_max[i][2], table_max[i][3], table_max[i+1][3], linestyle="dashed")
            plt.vlines(table_max[i][1], table_max[i][4], table_max[i+1][4], linestyle="dashed")
            plt.vlines(table_max[i][2], table_max[i][4], table_max[i+1][4], linestyle="dashed")
            plt.hlines(table_max[i][3], table_max[i][1], table_max[i][2], linestyle="dashed")
            plt.hlines(table_max[i][4], table_max[i][1], table_max[i][2], linestyle="dashed")
            plt.hlines(table_max[i+1][3], table_max[i][1], table_max[i][2], linestyle="dashed")
            plt.hlines(table_max[i+1][4], table_max[i][1], table_max[i][2], linestyle="dashed")
    plot.legend()
    plot.show()
if __name__ == '__main__':
    menu()
    while input('retry? [y/n]: ').lower().strip()[0] == 'y':
        print()
        menu()
