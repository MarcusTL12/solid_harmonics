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
    [x^(l - lz - ly) * y^ly * z^lz for lz in 0:l for ly in 0:(l-lz)]
end

function coeffs_in_order(l, m)
    ex = expand(s(l, m))
    [i => ex.coeff(cart)
     for (i, cart) in enumerate(generate_cartesians_in_order(l))
     if ex.coeff(cart) != 0]
end

function normalization(l)
    f = 2^l * x^(2l) / √Sym(π) * exp(-x^2)

    integrate(f, (x, -Inf, Inf))
end
