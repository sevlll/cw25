import sympy as sp

R, S, T = sp.symbols('R S T', real=True, positive=True)

cos_1t1 = 1 - 2 / (1 + T) ** 2
cos_1tr = 1 - 2 * R / ((R + T) * (1 + T))

cos_rts = (T ** 2 + T * (S + R) - S * R) / ((S + T) * (R + T))
cos_rst = (S ** 2 + S * (T + R) - R * T) / ((R + S) * (T + S))

cos_srt = (R ** 2 + R * (T + S) - S * T) / ((S + R) * (T + R))
cos_1rr = R / (1 + R)
cos_1rt = 1 - 2 * T / ((T + R) * (1 + R))

cos_111 = sp.Rational(1, 2)
cos_11t = 1 / (T + 1)
cos_r1t = (1 + (R + T) - R * T) / ((R + 1) * (T + 1))
cos_r1r = (1 + 2 * R - R * R) / ((R + 1) * (R + 1))


def make_sin(cos):
    return sp.powdenest(sp.factor(sp.simplify(sp.sqrt(1 - cos ** 2))), force=True)


sin_1t1 = make_sin(cos_1t1)
sin_1tr = make_sin(cos_1tr)

sin_rts = make_sin(cos_rts)

sin_rst = make_sin(cos_rst)
sin_1rr = make_sin(cos_1rr)

sin_1rt = make_sin(cos_1rt)
sin_111 = make_sin(cos_111)
sin_11t = make_sin(cos_11t)
sin_r1t = make_sin(cos_r1t)
sin_r1r = make_sin(cos_r1r)


def half(x):
    # returns cos(acos(x)/2)
    return sp.sqrt((1 + x) / 2)


def sum_of_cos(cos_alpha, cos_beta=1, cos_gamma=1):
    # returns cos(alpha + beta + gamma)
    sin_alpha = make_sin(cos_alpha)
    sin_beta = make_sin(cos_beta)
    sin_gamma = make_sin(cos_gamma)
    term1 = cos_alpha * cos_beta * cos_gamma
    term2 = sin_alpha * sin_beta * cos_gamma
    term3 = sin_beta * sin_gamma * cos_alpha
    term4 = sin_gamma * sin_alpha * cos_beta
    return term1 - term2 - term3 - term4


def dif_of_cos(cos_alpha, cos_beta=1):
    # returns cos(alpha - beta)
    sin_alpha = make_sin(cos_alpha)
    sin_beta = make_sin(cos_beta)
    return cos_alpha * cos_beta + sin_alpha * sin_beta


F1 = (sum_of_cos(cos_111, cos_11t, cos_r1t) + half(cos_r1r))
F1 = sp.simplify(F1)

FR = (sum_of_cos(cos_1rr, cos_1rt) + cos_srt)
FR = sp.simplify(FR)

FS = ((S ** 2 + R * S + T * S - R * T) / ((R + S) * (S + T)) - sp.Rational(1, 2))
FS = sp.simplify(FS)

FT = (sum_of_cos(cos_1tr, cos_rts) + half(cos_1t1))
FT = sp.simplify(FT)

print(F1)
print('\n')
print(FR)
print('\n')
print(FS)
print('\n')
print(FT)

