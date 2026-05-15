using SymPy

x, y, z = symbols("x y z")

r2 = x^2 + y^2 + z^2

function s(l, m)
    l = Sym(l)
    m = Sym(m)

    if !(-l <= m <= l)
        return Sym(0)
    end

    if l == 0
        return Sym(1)
    end

    l = l - 1

    δl0 = Sym(Int(l == 0))

    l = Sym(l)
    m = Sym(m)

    if abs(m) == l + 1
        c = √(2^δl0 * (2l + 1) / (2l + 2))

        xy = if m == l + 1
            x * s(l, l) - (1 - δl0) * y * s(l, -l)
        else
            y * s(l, l) + (1 - δl0) * x * s(l, -l)
        end

        simplify(c * xy)
    else
        n = (2l + 1) * z * s(l, m) - √((l + m) * (l - m)) * r2 * s(l - 1, m)
        d = √((l + m + 1) * (l - m + 1))

        simplify(n / d)
    end
end

function generate_cartesians_in_order(l)
    [(
        x^(l - lz - ly) * y^ly * z^lz,
        normalization(l - lz - ly) * normalization(ly) * normalization(lz)
    ) for lz in 0:l for ly in 0:(l-lz)]
end

function coeffs_in_order(l, m)
    ex = expand(s(l, m))
    [i => ex.coeff(cart) * √(Sym(n // normalization(l)))
     for (i, (cart, n)) in enumerate(generate_cartesians_in_order(l))
     if ex.coeff(cart) != 0]
end

function normalization(l)
    f = 2^l * x^(2l) / √Sym(π) * exp(-x^2)

    integrate(f, (x, -Inf, Inf))
end

function C(l, m, t, u, v)
    c = binomial(l, t) * binomial(l - t, abs(m) + t) *
        binomial(t, u) * binomial(abs(m), 2v + (m < 0))

    (iseven(t + v) ? c : -c) // 4^t
end

function NS(l, m)
    lf = factorial(l)

    f1 = factorial(big(l + m)) // lf
    f2 = factorial(l - m) // lf

    f3 = 1 // (2^(2 * abs(m) + (m == 0) - 1))

    √(Sym(f1 * f2 * f3))
end

function s_direct(l, m)
    ex = Sym(0)

    for t in 0:fld(l - abs(m), 2), u in 0:t, v in 0:fld(abs(m) - (m < 0), 2)
        ex += C(l, m, t, u, v) *
              x^(2t + abs(m) - 2u - 2v - (m < 0)) *
              y^(2u + 2v + (m < 0)) *
              z^(l - 2t - abs(m))
    end

    NS(l, m) * ex
end

function cart_to_index(lx, ly, lz)
    (lz * (2lx + 2ly + lz + 3) + 2ly) ÷ 2
end
