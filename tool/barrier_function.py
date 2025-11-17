import itertools
import numpy as np
import matplotlib.pyplot as plt


def barrier_function(x, t):
    return (-t) * np.log(-x)


def barrier_gradient(x):
    return -1 / (1 + np.exp(x))


def constraint_one(x, y):
    return 2 * x - y - 2


def constraint_two(x, y):
    return -1 - x + y


def constraint_three(x, y):
    return -x - y + 1


def objective(x):
    return 2 * x[0] + 3 * x[1]


def show_log_function():
    plt.xlim(-1, 12)
    plt.xlabel("x")
    plt.ylabel("Barrier function : -t*log10(-1.0 * x)")
    # plt.ylim(-5, 5)
    x = np.linspace(0.01, 10, 100)
    for i in range(10):
        y = barrier_function(x, i)
        plt.plot(x, y, label=f"factor t={i}")
    plt.legend()
    plt.show()


def show_boundaries():
    x = np.linspace(-1, 3, 100)
    y1 = 2 * x - 2
    y2 = x + 1
    y3 = 1 - x
    
    y_o = 2/3 * x  + 1
    plt.plot(x, y1, label="2x - y - 2")
    plt.plot(x, y2, label="y - x - 1")
    plt.plot(x, y3, label="1 - x - y")

    plt.plot(x, y_o, c = 'r', label="2x + 3y")
    plt.xlim(-1, 3)
    plt.ylim(-1, 4)
    plt.legend()


def show_barrier_function():
    SIZE = 300
    x = np.linspace(-1, 4, SIZE, dtype=np.float64)
    y = np.linspace(-1, 4, SIZE, dtype=np.float64)
    X, Y = np.meshgrid(x, y)
    idc = np.zeros(X.shape, dtype=int)

    for i, j in itertools.product(range(SIZE), range(SIZE)):
        v = (
            constraint_one(X[i][j], Y[i][j]) >= 0
            or constraint_two(X[i][j], Y[i][j]) >= 0
            or constraint_three(X[i][j], Y[i][j]) >= 0
        )
        if v:
            idc[i][j] = 1
    # idc = np.flipud(idc)
    [horizon, vertical] = np.where(idc > 0)

    loc = horizon * SIZE + vertical
    X = np.delete(X, loc)
    Y = np.delete(Y, loc)
    t = 0.01
    Z = (
        2 * X
        + 3 * Y
        + barrier_function(constraint_one(X, Y), t)
        + barrier_function(constraint_two(X, Y), t)
        + barrier_function(constraint_three(X, Y), t)
    )
    print(X.shape, Y.shape, Z.shape)
    # Z = Z.reshape((len(X) // 2, len(Y) // 2))
    plt.scatter(X, Y, c = Z, cmap='inferno', s=5, alpha=0.5)
    plt.colorbar()
    # print(Z)
    # plt.contour(X, Y, Z, 20)
    plt.show()


if __name__ == "__main__":
    # show_log_function()
    show_boundaries()
    show_barrier_function()
