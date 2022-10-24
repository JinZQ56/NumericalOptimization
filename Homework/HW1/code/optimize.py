import numpy as np
import matplotlib.pyplot as plt
import argparse
import time


class Rosenbrock:
    def __init__(self, dim: int):
        if dim % 2 == 0:
            self.n = dim
        else:
            raise Exception('dimension is false')

        self.fx = None
        self.dfx = None

    def f(self, x_vec):
        x_vec = np.array(x_vec)
        if x_vec.shape[0] != self.n:
            raise Exception('dimension is false')
        fx = 0
        for i in range(int(self.n/2)):
            i_1 = 2 * i
            i_2 = 2 * i + 1
            fx += 100 * (x_vec[i_1]**2 - x_vec[i_2])**2 + (x_vec[i_1] - 1)**2
        self.fx = fx
        return fx

    def df(self, x_vec):
        x_vec = np.array(x_vec)
        if x_vec.shape[0] != self.n:
            raise Exception('dimension is false')
        dfx = np.zeros(self.n)
        for i in range(int(self.n/2)):
            i_1 = 2 * i
            i_2 = 2 * i + 1
            dfx[i_1] = 400 * x_vec[i_1]**3 - 400 * x_vec[i_1] * x_vec[i_2] + 2 * x_vec[i_1] - 2
            dfx[i_2] = -200 * (x_vec[i_1]**2 - x_vec[i_2])
        self.dfx = dfx
        return dfx

    def ddf_2d(self, x_vec):
        x_vec = np.array(x_vec)
        if x_vec.shape[0] != 2:
            raise Exception('dimension is false')
        ddf_2d = np.zeros([2, 2])
        ddf_2d[0, 0] = 1200*x_vec[0]**2 - 400*x_vec[1] + 2
        ddf_2d[0, 1] = -400*x_vec[0]
        ddf_2d[1, 0] = -400*x_vec[0]
        ddf_2d[1, 1] = 200
        return ddf_2d

    @staticmethod
    def visualize2d():
        rb = Rosenbrock(2)
        x_1 = np.linspace(-2.0, 2.0, 400)
        x_2 = np.linspace(-1.0, 3.0, 400)
        X_1, X_2 = np.meshgrid(x_1, x_2)
        F = np.zeros(shape=X_1.shape)
        for i in range(F.shape[0]):
            for j in range(F.shape[1]):
                F[i, j] = rb.f([X_1[i, j], X_2[i, j]])
        fig, ax = plt.subplots(1, 1)
        levels = [0, 1, 10, 50, 100, 500, 1000, 10000]
        colors = []
        for i in range(10):
            colors.append("C{}".format(i))
        ct = ax.contourf(X_1, X_2, F, levels=levels, colors=colors)
        fig.colorbar(ct)
        ax.set_title('2D Rosenbrock Function')
        # plt.savefig('./2DRosenbrock.pdf')

        x = lssgd_dsz(rb, 0.001)
        # x = lssgd_armijo(rb, 1e-4)
        # x = newton(rb, 0.01)
        x = np.array(x)
        # print(x)
        ax.plot(x[:, 0], x[:, 1], color='white')
        ax.scatter(x[:, 0], x[:, 1], color='white', s=5, marker='*')

        plt.show()


def lssgd_dsz(func: Rosenbrock, c):
    x_s = np.zeros(func.n)
    x = np.array([-1.5, -0.5])
    s = time.time()
    iter_num = 0
    x_list = []
    for i in range(1, 20000):
        x_list.append(x)
        iter_num += 1
        d = -1 * func.df(x)
        if i/3000 == 0:
            c = c/0.1
        x = x + c * d
        if np.mean(np.abs(func.df(x))) < 1e-3:
            break
    e = time.time()
    print('*'*10, 'Linear-search Steepest Gradient Descent', '*'*10)
    print('constant: ', c)
    print('dimension of Rosenbrock function: ', func.n)
    print('start position: ', x_s)
    print('iteration number: ', iter_num)
    print('duration: ', e-s)
    print('final position: ', x)
    print('final gradient: ', func.df(x))
    print('minimum: ', func.f(x))

    return x_list


def lssgd_armijo(func: Rosenbrock, c):
    x_s = np.zeros(func.n)
    x = np.array([-1.5, -0.5])
    fx = func.f(x)
    iter_num = 0
    s = time.time()
    x_list = []
    for i in range(10000):
        x_list.append(x)
        iter_num += 1
        # print(i, ":", fx)
        d = -1 * func.df(x)
        t = 1
        while True:
            x_new = x + t * d
            dd = d @ d
            dd = -1 * c * t * dd
            if func.f(x_new) < func.f(x) + dd:
                break
            t = 0.5 * t
        x_old = x
        x = x + t * d
        # print(t)
        if np.mean(np.abs(func.df(x))) < 1e-3:
            break
        fx = func.f(x)
    e = time.time()

    print('*'*10, 'Linear-search Steepest Gradient Descent', '*'*10)
    print('constant: ', c)
    print('dimension of Rosenbrock function: ', func.n)
    print('start position: ', x_s)
    print('iteration number: ', iter_num)
    print('duration: ', e-s)
    print('final position: ', x)
    print('final gradient: ', func.df(x))
    print('minimum: ', func.f(x))

    return x_list


def newton(func: Rosenbrock, c):
    x = np.array([-1.5, -0.5])
    fx = func.f(x)
    iter_num = 0
    s = time.time()
    x_list = []
    while np.mean(np.abs(func.df(x))) > 1e-6:
        x_list.append(x)
        iter_num += 1
        g = func.df(x)
        h = func.ddf_2d(x)
        # m = h + c*np.ones(2)
        m = h
        d = np.linalg.solve(m, -1*g)
        t = 1
        while True:
            x_new = x + t * d
            if func.f(x_new) < func.f(x):
                break
            t = 0.5 * t
        x = x + t * d

    e = time.time()

    print('*'*10, 'Linear-search Steepest Gradient Descent', '*'*10)
    print('constant: ', c)
    print('dimension of Rosenbrock function: ', func.n)
    print('start position: ', [-2, -1])
    print('iteration number: ', iter_num)
    print('duration: ', e-s)
    print('final position: ', x)
    print('final gradient: ', func.df(x))
    print('minimum: ', func.f(x))

    return x_list


def quasi_newton(func: Rosenbrock, c):
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Linear-search Steepest Gradient Descent')
    # parser.add_argument('visualize2d', type=bool, help='Visualization of 2D Rosenbrock Function')
    parser.add_argument('dimension', type=int, help='Dimension of Rosenbrock Function')
    # parser.add_argument('constant', type=float, help='Constant in Algorithm between 0 to 1')
    args = parser.parse_args()

    Rosenbrock.visualize2d()

    # rb = Rosenbrock(args.dimension)
    # xl = lssgd_armijo(rb, 1e-3)



