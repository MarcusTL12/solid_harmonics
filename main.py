import sympy as sp

x, y, z = sp.symbols("x y z")

r2 = x**2 + y**2 + z**2


def s(l, m):
    l = sp.Number(l)
    m = sp.Number(m)

    if m < -l or l < m:
        return sp.Number(0)

    if l == 0:
        return sp.Number(1)

    l = l - 1

    dl0 = sp.Number(1) if l == 0 else sp.Number(0)

    if abs(m) == l + 1:
        c = sp.sqrt(2**dl0 * (2 * l + 1) / (2 * l + 2))

        if m == l + 1:
            xy = x * s(l, l) - (1 - dl0) * y * s(l, -l)
        else:
            xy = y * s(l, l) + (1 - dl0) * x * s(l, -l)

        return sp.simplify(c * xy)
    else:
        n = (2 * l + 1) * z * s(l, m) - \
            sp.sqrt((l + m) * (l - m)) * r2 * s(l - 1, m)
        d = sp.sqrt((l + m + 1) * (l - m + 1))

        return sp.simplify(n / d)


print(s(2, -2))
