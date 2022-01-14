import matplotlib.pyplot as plt
import numpy as np
import sympy
from tabulate import tabulate
import timeit

class Function:
    def __init__(self, s):
        self.symbols = {'x1': sympy.Symbol('x_1'), 'x2': sympy.Symbol('x_2')}
        self.s = s
        self.func = sympy.parsing.sympy_parser.parse_expr(self.s, self.symbols)

    def primegrad(self,f):
        self.fprime_sym = [f.func.diff(f.symbols['x_']) for f.symbols['x_'] in
                           (f.symbols['x1'], f.symbols['x2'])]
        return self.fprime_sym

    def hessgrad(self,f):
        f.primegrad(f)
        self.matr=sympy.Matrix(f.fprime_sym)
        self.fhess_sym = [[f.func.diff(f.symbols['x1_'], f.symbols['x2_']) for f.symbols['x1_'] in
                           (f.symbols['x1'], f.symbols['x2'])] for f.symbols['x2_'] in (f.symbols['x1'], f.symbols['x2'])]
        return self.fhess_sym,self.matr

    def Hessian(self,f):
        f.hessgrad(f)
        self.hs=sympy.Matrix(f.fhess_sym)
        return self.hs

def newton_descent_iter(*, f,
                        f_grad,
                        f_hess, eps,  xstart,  max_steps, table_min):
    x0 = xstart
    xe = [eps,eps]
    iter_cnt = 0
    xb = [0,0]
    while True:
        if iter_cnt and iter_cnt % 100000 == 0:
            print("newton iterat: %d" % iter_cnt)
            print(x0)
        grad = f_grad(x0)
        hess = f_hess(x0)


        try:
            hess_inv = np.linalg.inv(hess)
        except np.linalg.LinAlgError:
            hess_inv = ([[0, 0],
                         [0, 0]])
            xb = ([0,0])
            xe = ([0, 0])
            table_min[iter_cnt] = [iter_cnt, grad, hess, hess_inv, x0, xe, xb]
            return tabulate(list(table_min.values()),
                            headers=['iter_cnt', 'grad', 'hess', 'hess_inv', 'x0', 'xe', 'xb'])


        grad_hessinv = np.matmul(grad, hess_inv)

        if abs(xe[0]) < eps or abs(xe[1]) < eps :
            return tabulate(list(table_min.values()),
                            headers=['iter_cnt', 'grad', 'hess', 'hess_inv', 'x0', 'xe', 'xb'])

        if iter_cnt == max_steps:
            return tabulate(list(table_min.values()),
                            headers=['iter_cnt', 'grad', 'hess', 'hess_inv', 'x0', 'xe', 'xb'])

        table_min[iter_cnt] = [iter_cnt, grad, hess, hess_inv, x0, xe, xb]
        xb = x0 + grad_hessinv
        xe=xb-x0
        x0 = [xb[0], xb[1]]
        iter_cnt += 1
        table_min[iter_cnt] = [iter_cnt, grad, hess, hess_inv, x0, xe, xb]

    return tabulate(list(table_min.values()), headers=['iter_cnt', 'grad', 'hess', 'hess_inv', 'x0','xe', 'xb'])


def newton_descent(*args, **kwargs) -> np.array:
    return (newton_descent_iter(*args, **kwargs))


def gyperballf_mark(var2,xstart,gyperballf):
    var2t=Function(var2)
    var2t=sympy.lambdify((var2t.symbols['x1'], var2t.symbols['x2']), var2t.func, 'numpy')
    var2t=var2t(xstart[0], xstart[1])

    gyperballfv1=0
    if var2t == 0:
        gyperballf = ('0')
    return gyperballf


class FunctionPlotter:
    def __init__(self,f_lmbda, xt, yt,var2):
       x = [-10.1 + 0.2 * i for i in range(100)]
       y = [-10.2 + 0.2 * i for i in range(100)]
       zt = list()

       gyperballf2 = ('1/(%s)') % var2
       ft = Function(gyperballf2)
       f_lmbdat = sympy.lambdify((ft.symbols['x1'], ft.symbols['x2']), ft.func, 'numpy')

       zt2 = list()

       for i in range(100):
          zt.append([])
          zt2.append([])
          for j in range(100):
              try:
                  f_lmbda(x[i], y[j])
                  zt[i].append(f_lmbda(x[i], y[j]))
              except ZeroDivisionError:
                zt[i].append(f_lmbda(x[i], y[j - 1]))

              try:
                  f_lmbdat(x[i], y[j])
                  zt2[i].append(f_lmbdat(x[i], y[j]))
              except ZeroDivisionError:
                  zt2[i].append(f_lmbdat(x[i], y[j - 1]))

       h = plt.contourf(x,y,zt)
       plt.contour(h, colors = 'k')
       plt.colorbar(h)
       plt.scatter(xt, yt, color = 'red', marker = 'x', label='minimum')
       plt.plot(xt, yt, color='red', marker='x', label='minimum')
       plt.xlim(-10.1,xt[len(xt)-1] + 5.5)
       plt.ylim(-10.2,yt[len(xt)-1] + 5)
       plt.show()

       plt.scatter(xt, yt, color='red', marker='x', label='minimum')
       plt.plot(xt, yt, color='red', marker='x', label='minimum')
       plt.xlim(-10.1, xt[len(xt) - 1] + 5.5)
       plt.ylim(-10.2, yt[len(xt) - 1] + 5)

       h2 = plt.contourf(x, y, zt2)
       plt.contour(h2, colors='black')
       plt.colorbar(h2)
       plt.show()


       h = plt.contourf(x,y,zt)
       plt.contour(h, colors = 'k')
       plt.colorbar(h)

       h2 = plt.contourf(x, y, zt2)
       plt.contour(h2, colors='black')
       plt.colorbar(h2)

       plt.scatter(xt, yt, color='red', marker='x', label='minimum')
       plt.plot(xt, yt, color='red', marker='x', label='minimum')
       plt.xlim(-10.1, xt[len(xt) - 1] + 5.5)
       plt.ylim(-10.2, yt[len(xt) - 1] + 5)




    def show(self):
        plt.show()

    def clear(self):
        plt.clf()

def ploter(table_min,f_lmbda,var2):
    xt=[0]*len(table_min)
    yt=[0]*len(table_min)
    for i in range(len(table_min)):
        xt[i]=((table_min[i])[6])[0]
        yt[i] = ((table_min[i])[6])[1]
    plot = FunctionPlotter(f_lmbda, xt, yt,var2)
    plot.show()

def bar_inter(xstart, eps, max_steps,var1,var2):

    gyperballf = ('1/(%s)') % var2

    gyperballf=gyperballf_mark(var2, xstart,gyperballf) # проверка на 0


    for zPower in range(0, 3):
        z = 100. * (0.1**(1 * zPower))
        z = -z
        print(z)

        var3 = ('((%s) + (%s*(%s)))') % (var1, z, gyperballf)  # sympy.ln

        print("z^keff=", -z)

        f = Function(var3)
        Function(var3).Hessian(f)


        f_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.func, 'numpy')
        fprime_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.fprime_sym, 'numpy')
        fhess_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.fhess_sym, 'numpy')

        def DelenieXY_X_Y(f):
            return lambda X: np.array(f(X[0], X[1]))

        f = DelenieXY_X_Y(f_lmbda)
        fprime = DelenieXY_X_Y(fprime_lmbda)
        fhess = DelenieXY_X_Y(fhess_lmbda)

        table_min = {0: [0, 0, 0, 0, 0, 0, 0]}


        Iterations =  newton_descent(f=f, f_grad=fprime, f_hess=fhess, eps=eps, xstart=xstart, max_steps=max_steps, table_min=table_min)
        print(Iterations)


        print('all iter', len(table_min), "solution fxy and hess",
              len(table_min))


        ploter(table_min,f_lmbda,var2)


def menu():
    var1 = '0.5*((x2)*(x1))'
    var2 = 'x1+x2-8'
    print(var1)
    print(var2)

    a = float(input('dot start x '))
    b = float(input('dot start y '))
    max_steps = int(input('[max step]    n = '))
    eps = float(0.001)
    print(eps)

    start_time = timeit.default_timer()

    print(timeit.default_timer() - start_time, 'time')

    bar_inter(xstart=np.array([a, b]), eps=eps, max_steps=max_steps,
                  var1=var1, var2=var2)

    print(timeit.default_timer() - start_time, 'time')


if __name__ == '__main__':
    menu()
    while input('retry? [y/n]: ').lower().strip()[0] == 'y':
        print()
        menu()