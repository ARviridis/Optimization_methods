from scipy import optimize
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

    def proizv(self,f):
        self.fprime_sym = [f.func.diff(f.symbols['x_']) for f.symbols['x_'] in
                           (f.symbols['x1'], f.symbols['x2'])]
        return self.fprime_sym

    def grad(self,f):
        f.proizv(f)
        self.matr=sympy.Matrix(f.fprime_sym)
        self.fhess_sym = [[f.func.diff(f.symbols['x1_'], f.symbols['x2_']) for f.symbols['x1_'] in
                           (f.symbols['x1'], f.symbols['x2'])] for f.symbols['x2_'] in (f.symbols['x1'], f.symbols['x2'])]
        return self.fhess_sym,self.matr

    def Hessian(self,f):
        f.grad(f)
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

        hess_inv = np.linalg.inv(hess)
        grad_hessinv = np.matmul(grad, hess_inv)

        if abs(xe[0]) < eps or abs(xe[1]) < eps:
            break
        if iter_cnt == max_steps:
            break
        table_min[iter_cnt] = [iter_cnt, grad, hess, hess_inv, x0, xe, xb]
        xb = x0 - grad_hessinv
        xe=xb-x0
        x0 = [xb[0], xb[1]]
        iter_cnt += 1
        table_min[iter_cnt] = [iter_cnt, grad, hess, hess_inv, x0, xe, xb]



def newton_descent(*args, **kwargs) -> np.array:
    return (newton_descent_iter(*args, **kwargs))


class FunctionPlotter:
    def __init__(self,f_lmbda, xt, yt):
       x = [-0.1 + 0.1 * i for i in range(100)]
       y = [-0.2 + 0.1 * i for i in range(100)]
       z = list()
       for i in range(100):
          z.append([])
          for j in range(100):
             z[i].append(f_lmbda(x[i], y[j]))
       h = plt.contourf(x,y,z)
       plt.contour(h, colors = 'k')
       plt.colorbar(h)
       plt.scatter(xt, yt, color = 'red', marker = 'x', label='minimum')
       plt.plot(xt, yt, color='red', marker='x', label='minimum')
       plt.xlim(-0.1,xt[len(xt)-1] + 0.5)
       plt.ylim(-0.2,yt[len(xt)-1] + 2)

    def show(self):
        plt.show()

    def clear(self):
        plt.clf()


def menu():
    var1 = ('sqrt(x1**2+2**2)+sqrt((x2-x1)**2+4**2)+sqrt((6-x2)**2+(3-4)**2)')
    print(var1)

    f = Function(var1)
    Function(var1).Hessian(f)

    a = float(input('dot start x '))
    b = float(input('dot start y '))
    max_steps = int(input('[max step]    n = '))
    eps = float(0.001)
    print(eps)

    start_time = timeit.default_timer()


    f_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.func, 'numpy')
    fprime_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.fprime_sym, 'numpy')
    fhess_lmbda = sympy.lambdify((f.symbols['x1'], f.symbols['x2']), f.fhess_sym, 'numpy')

    def DelenieXY_X_Y(f):
        return lambda X: np.array(f(X[0], X[1]))

    f = DelenieXY_X_Y(f_lmbda)
    fprime = DelenieXY_X_Y(fprime_lmbda)
    fhess = DelenieXY_X_Y(fhess_lmbda)

    table_min = {0: [0, 0, 0, 0, 0, 0, 0]}


    min = newton_descent(f=f, f_grad=fprime, f_hess=fhess, eps=eps, xstart=np.array([a, b]), max_steps=max_steps, table_min=table_min)
    print(tabulate(list(table_min.values()), headers=['iter_cnt', 'grad', 'hess', 'hess_inv', 'x0','xe', 'xb']))
    print('all iter', len(table_min), "solution fxy and hess", len(table_min) )
    print(timeit.default_timer() - start_time, 'time')
    # (newton_scipy)
    print('scipy')
    x_opt = optimize.fmin_ncg(f, (0, 0), fprime=fprime, fhess=fhess)
    x_opt

    xt=[0]*len(table_min)
    yt=[0]*len(table_min)
    for i in range(len(table_min)):
        xt[i]=((table_min[i])[6])[0]
        yt[i] = ((table_min[i])[6])[1]

    plot = FunctionPlotter(f_lmbda, xt, yt)
    plot.show()
if __name__ == '__main__':
    menu()
    while input('retry? [y/n]: ').lower().strip()[0] == 'y':
        print()
        menu()