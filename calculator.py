"""
TILLIDON fx-991MS  -  Scientific Calculator  (Python + Tkinter)
Premium resizable canvas interface (F11 = full screen), layout traced from a real fx-991MS.

Run:   python fx991ms_premium_calculator.py          (optional:  ... 0.5   = window scale)
Linux: sudo apt install python3-tk   if Tkinter is missing.

Press  Help > Key guide  inside the app for every key, mode and shortcut.
"""

import ast
import cmath
import math
import operator
import random
import re
import tkinter as tk
from fractions import Fraction
from tkinter import ttk

# ======================================================================
#  Number formatting helpers
# ======================================================================


def plain(x, sig=10):
    if x == 0:
        return "0"
    dec = max(0, sig - 1 - math.floor(math.log10(abs(x))))
    s = f"{x:.{dec}f}"
    return s.rstrip("0").rstrip(".") if "." in s else s


def sci(x, p, strip=True):
    m, e = f"{x:.{p}e}".split("e")
    if strip and "." in m:
        m = m.rstrip("0").rstrip(".")
    return f"{m}×10^{int(e)}"


def eng(x):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x)) / 3) * 3)
    return f"{plain(x / 10 ** e)}×10^{e}"


def to_dms(x):
    sgn, x = ("-" if x < 0 else ""), abs(x)
    d = int(x)
    m = int((x - d) * 60)
    s = round(((x - d) * 60 - m) * 60, 4)
    if s >= 60:
        s, m = s - 60, m + 1
    if m >= 60:
        m, d = m - 60, d + 1
    return f"{sgn}{d}°{m}′{plain(s, 8)}″"


def lit(v):
    """number -> text that the calculator parser understands"""
    return "(" + repr(float(v)).replace("e", "ᴇ") + ")"


# ======================================================================
#  Calculus, linear algebra, base-n helpers
# ======================================================================


def integrate(g, a, b, tol=1e-10):
    def simp(fa, fm, fb, a, b):
        return (b - a) / 6 * (fa + 4 * fm + fb)

    def rec(a, b, fa, fm, fb, whole, tol, depth):
        m, lm, rm = (a + b) / 2, (3 * a + b) / 4, (a + 3 * b) / 4
        flm, frm = g(lm), g(rm)
        left, right = simp(fa, flm, fm, a, m), simp(fm, frm, fb, m, b)
        if depth <= 0 or abs(left + right - whole) <= 15 * tol:
            return left + right + (left + right - whole) / 15
        return (rec(a, m, fa, flm, fm, left, tol / 2, depth - 1) +
                rec(m, b, fm, frm, fb, right, tol / 2, depth - 1))

    if a == b:
        return 0.0
    total, n = 0.0, 8
    h = (b - a) / n
    for i in range(n):
        lo, hi = a + i * h, a + (i + 1) * h
        fa, fb, fm = g(lo), g(hi), g((lo + hi) / 2)
        total += rec(lo, hi, fa, fm, fb, simp(fa, fm, fb, lo, hi), tol / n, 25)
    return total


def derivative(g, a):
    h = 1e-3 * max(1.0, abs(a))
    return (-g(a + 2 * h) + 8 * g(a + h) - 8 * g(a - h) + g(a - 2 * h)) / (12 * h)


def solve_linear(A, b):
    n = len(A)
    M = [list(map(float, row)) + [float(b[i])] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-12:
            raise ValueError("Math ERROR (no unique solution)")
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def poly_roots(c):
    if c[0] == 0:
        raise ValueError("Coefficient a cannot be 0")
    if len(c) == 3:
        a, b, cc = c
        d = cmath.sqrt(b * b - 4 * a * cc)
        z = [(-b + d) / (2 * a), (-b - d) / (2 * a)]
    else:
        n, a = len(c) - 1, [x / c[0] for x in c]
        z = [(0.4 + 0.9j) ** k for k in range(n)]
        for _ in range(400):
            for i in range(n):
                p = 1
                for coef in a[1:]:
                    p = p * z[i] + coef
                d = 1
                for j in range(n):
                    if j != i:
                        d *= z[i] - z[j]
                z[i] -= p / d
    out = []
    for r in z:
        r = complex(round(r.real, 9), round(r.imag, 9))
        out.append(r.real if r.imag == 0 else r)
    return out


class Matrix:
    def __init__(self, rows):
        self.r = [[float(x) for x in row] for row in rows]

    n = property(lambda self: len(self.r))
    m = property(lambda self: len(self.r[0]))

    def _same(self, o):
        if not isinstance(o, Matrix) or (self.n, self.m) != (o.n, o.m):
            raise ValueError("Dimension ERROR")

    def __add__(self, o):
        self._same(o)
        return Matrix([[a + b for a, b in zip(x, y)] for x, y in zip(self.r, o.r)])

    def __sub__(self, o):
        self._same(o)
        return Matrix([[a - b for a, b in zip(x, y)] for x, y in zip(self.r, o.r)])

    def __mul__(self, o):
        if isinstance(o, Matrix):
            if self.m != o.n:
                raise ValueError("Dimension ERROR")
            return Matrix([[sum(a * b for a, b in zip(row, col)) for col in zip(*o.r)]
                           for row in self.r])
        return Matrix([[a * o for a in row] for row in self.r])

    __rmul__ = __mul__

    def __truediv__(self, o):
        return self * (1 / o)

    def __neg__(self):
        return self * -1

    def __pow__(self, k):
        if k == -1:
            return self.inv()
        if k != int(k) or k < 0 or self.n != self.m:
            raise ValueError("Math ERROR")
        res = Matrix([[1.0 if i == j else 0.0 for j in range(self.n)] for i in range(self.n)])
        for _ in range(int(k)):
            res = res * self
        return res

    def T(self):
        return Matrix(list(zip(*self.r)))

    def det(self):
        if self.n != self.m:
            raise ValueError("Dimension ERROR")
        M, d = [row[:] for row in self.r], 1.0
        for c in range(self.n):
            p = max(range(c, self.n), key=lambda r: abs(M[r][c]))
            if abs(M[p][c]) < 1e-14:
                return 0.0
            if p != c:
                M[c], M[p] = M[p], M[c]
                d = -d
            d *= M[c][c]
            for r in range(c + 1, self.n):
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        return d

    def inv(self):
        if self.n != self.m:
            raise ValueError("Dimension ERROR")
        cols = [solve_linear(self.r, [1.0 if i == j else 0.0 for i in range(self.n)])
                for j in range(self.n)]
        return Matrix(cols).T()

    def norm(self):
        return math.sqrt(sum(x * x for row in self.r for x in row))

    def text(self, f):
        return "\n".join("   ".join(f(x) for x in row) for row in self.r)


def flat(m):
    return [x for row in m.r for x in row]


def dot(a, b):
    u, v = flat(a), flat(b)
    if len(u) != len(v):
        raise ValueError("Dimension ERROR")
    return sum(x * y for x, y in zip(u, v))


def cross(a, b):
    u, v = flat(a), flat(b)
    if len(u) != 3 or len(v) != 3:
        raise ValueError("Dimension ERROR")
    return Matrix([[u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0]]])


BASE_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
            ast.FloorDiv: operator.floordiv, ast.BitAnd: operator.and_,
            ast.BitOr: operator.or_, ast.BitXor: operator.xor,
            ast.MatMult: lambda a, b: ~(a ^ b)}          # MatMult = xnor


def eval_base(text, base):
    s = text.upper().replace("×", "*").replace("÷", "/").replace("−", "-")
    s = re.sub(r"(?<!/)/(?!/)", "//", s)
    for w, o in (("XNOR", "@"), ("XOR", "^"), ("AND", "&"), ("OR", "|"),
                 ("NOT", "~"), ("NEG", "-")):
        s = re.sub(rf"\b{w}\b", o, s)
    s = re.sub(r"[0-9A-F]+", lambda m: str(int(m.group(), base)), s)

    def ev(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, int):
            return n.value
        if isinstance(n, ast.BinOp) and type(n.op) in BASE_OPS:
            return BASE_OPS[type(n.op)](ev(n.left), ev(n.right))
        if isinstance(n, ast.UnaryOp):
            f = {ast.USub: operator.neg, ast.UAdd: operator.pos, ast.Invert: operator.invert}
            if type(n.op) in f:
                return f[type(n.op)](ev(n.operand))
        raise ValueError("Syntax ERROR")

    v = ev(ast.parse(s, mode="eval").body)
    if not -2 ** 31 <= v < 2 ** 31:
        raise ValueError("Math ERROR")
    return v


def base_str(v, base):
    v = int(v)
    if base == 10:
        return str(v)
    u = v & 0xFFFFFFFF
    return {16: format(u, "X"), 2: format(u, "b"), 8: format(u, "o")}[base]


# ======================================================================
#  Statistics / regression
# ======================================================================


def sd_stats(rows):
    n = sum(w for _, w in rows)
    sx = sum(x * w for x, w in rows)
    sx2 = sum(x * x * w for x, w in rows)
    mean = sx / n
    ss = sum(w * (x - mean) ** 2 for x, w in rows)
    return {"n": n, "Σx": sx, "Σx²": sx2, "x̄ (mean)": mean, "σn": math.sqrt(ss / n),
            "σn-1": math.sqrt(ss / (n - 1)) if n > 1 else float("nan"),
            "min": min(x for x, _ in rows), "max": max(x for x, _ in rows)}


TRANS = {"Lin": (lambda x: x, lambda y: y), "Log": (math.log, lambda y: y),
         "Exp": (lambda x: x, math.log), "Pwr": (math.log, math.log),
         "Inv": (lambda x: 1 / x, lambda y: y)}


def regression(rows, kind):
    n = sum(w for *_, w in rows)
    S = lambda f: sum(w * f(x, y) for x, y, w in rows)
    res = {"n": n, "Σx": S(lambda x, y: x), "Σx²": S(lambda x, y: x * x),
           "Σy": S(lambda x, y: y), "Σy²": S(lambda x, y: y * y), "Σxy": S(lambda x, y: x * y)}
    mx, my = res["Σx"] / n, res["Σy"] / n
    vx, vy = S(lambda x, y: (x - mx) ** 2), S(lambda x, y: (y - my) ** 2)
    res.update({"x̄": mx, "ȳ": my, "σx": math.sqrt(vx / n), "σy": math.sqrt(vy / n)})
    if n > 1:
        res.update({"sx": math.sqrt(vx / (n - 1)), "sy": math.sqrt(vy / (n - 1))})
    if kind == "Quad":
        p = lambda k: S(lambda x, y: x ** k)
        q = lambda k: S(lambda x, y: y * x ** k)
        A, B, C = solve_linear([[n, p(1), p(2)], [p(1), p(2), p(3)], [p(2), p(3), p(4)]],
                               [q(0), q(1), q(2)])
        res.update({"A": A, "B": B, "C": C})
        yhat = lambda x: A + B * x + C * x * x

        def xhat(y):
            d = B * B - 4 * C * (A - y)
            if d < 0:
                raise ValueError("Math ERROR")
            return [(-B + math.sqrt(d)) / (2 * C), (-B - math.sqrt(d)) / (2 * C)]
        return res, yhat, xhat
    tx, ty = TRANS[kind]
    pts = [(tx(x), ty(y), w) for x, y, w in rows]
    ax = sum(w * a for a, b, w in pts) / n
    by = sum(w * b for a, b, w in pts) / n
    sxx = sum(w * (a - ax) ** 2 for a, b, w in pts)
    syy = sum(w * (b - by) ** 2 for a, b, w in pts)
    sxy = sum(w * (a - ax) * (b - by) for a, b, w in pts)
    b0 = sxy / sxx
    a0 = by - b0 * ax
    A = math.exp(a0) if kind in ("Exp", "Pwr") else a0
    B = b0
    res.update({"A": A, "B": B, "r": sxy / math.sqrt(sxx * syy)})
    yhat = {"Lin": lambda x: A + B * x, "Log": lambda x: A + B * math.log(x),
            "Exp": lambda x: A * math.exp(B * x), "Pwr": lambda x: A * x ** B,
            "Inv": lambda x: A + B / x}[kind]
    xh = {"Lin": lambda y: (y - A) / B, "Log": lambda y: math.exp((y - A) / B),
          "Exp": lambda y: math.log(y / A) / B, "Pwr": lambda y: (y / A) ** (1 / B),
          "Inv": lambda y: B / (y - A)}[kind]
    return res, yhat, lambda y: [xh(y)]


# ======================================================================
#  Constants and conversions (fx-991MS tables)
# ======================================================================

CONSTANTS = [
    ("Proton mass  mp (kg)", 1.67262192369e-27), ("Neutron mass  mn (kg)", 1.67492749804e-27),
    ("Electron mass  me (kg)", 9.1093837015e-31), ("Muon mass  mμ (kg)", 1.883531627e-28),
    ("Bohr radius  a0 (m)", 5.29177210903e-11), ("Planck constant  h (J·s)", 6.62607015e-34),
    ("Nuclear magneton  μN", 5.0507837461e-27), ("Bohr magneton  μB", 9.2740100783e-24),
    ("Reduced Planck  ħ", 1.054571817e-34), ("Fine-structure constant  α", 7.2973525693e-3),
    ("Classical electron radius  re", 2.8179403262e-15), ("Compton wavelength  λc", 2.42631023867e-12),
    ("Proton gyromagnetic ratio  γp", 2.6752218744e8), ("Proton Compton wavelength  λcp", 1.32140985539e-15),
    ("Neutron Compton wavelength  λcn", 1.31959090581e-15), ("Rydberg constant  R∞", 10973731.56816),
    ("Atomic mass unit  u (kg)", 1.6605390666e-27), ("Proton magnetic moment  μp", 1.41060679736e-26),
    ("Electron magnetic moment  μe", -9.2847647043e-24), ("Neutron magnetic moment  μn", -9.6623651e-27),
    ("Muon magnetic moment  μμ", -4.4904483e-26), ("Faraday constant  F (C/mol)", 96485.33212),
    ("Elementary charge  e (C)", 1.602176634e-19), ("Avogadro constant  NA", 6.02214076e23),
    ("Boltzmann constant  k (J/K)", 1.380649e-23), ("Molar volume ideal gas  Vm", 0.022413969545),
    ("Molar gas constant  R", 8.314462618), ("Speed of light  c0 (m/s)", 299792458),
    ("First radiation constant  c1", 3.741771852e-16), ("Second radiation constant  c2", 0.01438776877),
    ("Stefan-Boltzmann  σ", 5.670374419e-8), ("Electric constant  ε0", 8.8541878128e-12),
    ("Magnetic constant  μ0", 1.25663706212e-6), ("Magnetic flux quantum  φ0", 2.067833848e-15),
    ("Standard gravity  g (m/s²)", 9.80665), ("Conductance quantum  G0", 7.748091729e-5),
    ("Vacuum impedance  Z0 (Ω)", 376.730313668), ("Celsius temperature  t (K)", 273.15),
    ("Gravitational constant  G", 6.6743e-11), ("Standard atmosphere  atm (Pa)", 101325),
]


def lin(k):
    return (lambda x: x * k, lambda x: x / k)


CONVERSIONS = [
    ("in", "cm", *lin(2.54)), ("ft", "m", *lin(0.3048)), ("yd", "m", *lin(0.9144)),
    ("mile", "km", *lin(1.609344)), ("n mile", "m", *lin(1852)),
    ("acre", "m²", *lin(4046.8564224)), ("gal (US)", "L", *lin(3.785411784)),
    ("gal (UK)", "L", *lin(4.54609)), ("pc", "km", *lin(3.0856775814913673e13)),
    ("km/h", "m/s", *lin(1 / 3.6)), ("oz", "g", *lin(28.349523125)),
    ("lb", "kg", *lin(0.45359237)), ("atm", "Pa", *lin(101325)),
    ("mmHg", "Pa", *lin(133.322387415)), ("hp", "kW", *lin(0.745699872)),
    ("kgf/cm²", "Pa", *lin(98066.5)), ("kgf·m", "J", *lin(9.80665)),
    ("lbf/in²", "kPa", *lin(6.894757293168)),
    ("°F", "°C", lambda x: (x - 32) / 1.8, lambda x: x * 1.8 + 32),
    ("J", "cal", *lin(1 / 4.1855)),
]

# ======================================================================
#  Calculation engine
# ======================================================================

BIN_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
           ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod}
UN_OPS = {ast.UAdd: operator.pos, ast.USub: operator.neg}


class Engine:
    def __init__(self):
        self.angle = "DEG"
        self.complex_mode = False
        self.ans = 0
        self.vars = {k: 0 for k in "ABCDEFXYM"}
        self.pair = ""
        self.dm, self.dd = "norm1", 2            # display mode / digits
        self.frac, self.eng_on, self.polar = None, False, False

    # ---- angle helpers
    def _k(self):
        return {"DEG": math.pi / 180, "RAD": 1, "GRA": math.pi / 200}[self.angle]

    def _in(self, x):
        return x * self._k()

    def _out(self, x):
        return x / self._k()

    # ---- functions
    def functions(self):
        M = cmath if self.complex_mode else math
        cl = lambda v: v if isinstance(v, complex) or abs(v) > 1e-15 else 0.0

        def tan(x):
            a = self._in(x)
            if not self.complex_mode and abs(math.cos(a)) < 1e-15:
                raise ValueError("Math ERROR")
            return cl(M.tan(a))

        def pol(x, y):
            r, t = math.hypot(x, y), self._out(math.atan2(y, x))
            self.vars["X"], self.vars["Y"] = r, t
            self.pair = f"r = {self.fmt_real(r)}    θ = {self.fmt_real(t)}"
            return r

        def rec(r, t):
            x, y = r * math.cos(self._in(t)), r * math.sin(self._in(t))
            self.vars["X"], self.vars["Y"] = cl(x), cl(y)
            self.pair = f"x = {self.fmt_real(cl(x))}    y = {self.fmt_real(cl(y))}"
            return cl(x)

        def root(x, y):
            if x < 0 and y % 2:
                return -((-x) ** (1 / y))
            return x ** (1 / y)

        return {
            "sin": lambda x: cl(M.sin(self._in(x))), "cos": lambda x: cl(M.cos(self._in(x))),
            "tan": tan,
            "asin": lambda x: self._out(M.asin(x)), "acos": lambda x: self._out(M.acos(x)),
            "atan": lambda x: self._out(M.atan(x)),
            "sinh": M.sinh, "cosh": M.cosh, "tanh": M.tanh,
            "asinh": M.asinh, "acosh": M.acosh, "atanh": M.atanh,
            "ln": M.log, "log": M.log10, "exp": M.exp, "sqrt": M.sqrt,
            "cbrt": lambda x: math.copysign(abs(x) ** (1 / 3), x), "root": root,
            "abs": lambda x: x.norm() if isinstance(x, Matrix) else abs(x),
            "int": math.trunc, "intg": math.floor, "rnd": self._rnd,
            "ran": lambda: random.randint(0, 999) / 1000,
            "fact": lambda x: math.factorial(int(x)) if float(x).is_integer() and x >= 0
            else math.gamma(x + 1),
            "nPr": lambda n, r: math.perm(int(n), int(r)),
            "nCr": lambda n, r: math.comb(int(n), int(r)),
            "pol": pol, "rec": rec,
            "conjg": lambda z: complex(z).conjugate(), "arg": lambda z: self._out(cmath.phase(z)),
            "rep": lambda z: complex(z).real, "imp": lambda z: complex(z).imag,
            "cpolar": lambda r, t: cmath.rect(r, self._in(t)),
        }

    def _rnd(self, x):
        if self.dm == "fix":
            return round(x, self.dd)
        if self.dm == "sci":
            return float(f"{x:.{(self.dd or 10) - 1}e}")
        return float(f"{x:.10g}")

    # ---- parsing
    @staticmethod
    def preprocess(text):
        s = text.strip().replace("d/dx(", "deriv(").replace("∫(", "integ(").replace("ˣ√", "‡")
        s = re.sub(r"(\d+)⌟(\d+)⌟(\d+)", r"(\1+\2/\3)", s)
        s = re.sub(r"(\d+)⌟(\d+)", r"(\1/\2)", s)
        s = re.sub(r"(\d+(?:\.\d+)?)°(?:\s*(\d+(?:\.\d+)?)′)?(?:\s*(\d+(?:\.\d+)?)″)?",
                   lambda m: f"({m.group(1)}+{m.group(2) or 0}/60+{m.group(3) or 0}/3600)", s)
        s = re.sub(r"(\d+(?:\.\d+)?)∠([-−]?\d+(?:\.\d+)?)", r"cpolar(\1,\2)", s)
        for a, b in (("×", "*"), ("÷", "/"), ("−", "-"), ("^", "**"), ("π", "pi")):
            s = s.replace(a, b)
        s = re.sub(r"√\s*(\d+(?:\.\d+)?)", r"sqrt(\1)", s).replace("√", "sqrt")
        s = re.sub(r"ᴇ\s*([+-]?\d+)", r"*10**(\1)", s)
        opd = r"(\d+(?:\.\d+)?|\([^()]*\)|Ans|[A-FXYM])"
        s = re.sub(opd + "!", r"fact(\1)", s)
        s = re.sub(opd + r"\s*(nPr|nCr)\s*" + opd, r"\2(\1,\3)", s)
        s = re.sub(opd + r"\s*‡\s*" + opd, r"root(\2,\1)", s)
        s = s.replace("%", "/100")
        s = re.sub(r"(?<=[\d.)])\s*(?=[A-Za-z_(])", "*", s)
        s = re.sub(r"(?<=\))\s*(?=\d)", "*", s)
        return s

    def evaluate(self, text, extra_env=None, extra_fn=None):
        if ":" in text:
            *head, text = text.split(":")
            for part in head:
                self.ans = self.evaluate(part, extra_env, extra_fn)
        self.pair = ""
        tree = ast.parse(self.preprocess(text), mode="eval")
        env = {"pi": math.pi, "e": math.e, "i": 1j, "Ans": self.ans, **self.vars,
               **(extra_env or {})}
        fn = {**self.functions(), **(extra_fn or {})}
        v = self._ev(tree.body, env, fn)
        if isinstance(v, complex) and not self.complex_mode:
            raise ValueError("Math ERROR")
        return v

    def _ev(self, n, env, fn):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return n.value
        if isinstance(n, ast.Name) and n.id in env:
            return env[n.id]
        if isinstance(n, ast.BinOp) and type(n.op) in BIN_OPS:
            a, b = self._ev(n.left, env, fn), self._ev(n.right, env, fn)
            if isinstance(n.op, ast.Pow) and not isinstance(b, complex) and abs(b) > 10000:
                raise ValueError("Math ERROR")
            return BIN_OPS[type(n.op)](a, b)
        if isinstance(n, ast.UnaryOp) and type(n.op) in UN_OPS:
            return UN_OPS[type(n.op)](self._ev(n.operand, env, fn))
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name):
            name = n.func.id
            if name in ("integ", "deriv"):
                g = lambda t: self._ev(n.args[0], {**env, "X": t}, fn)
                lim = [self._ev(a, env, fn) for a in n.args[1:]]
                return integrate(g, *lim) if name == "integ" else derivative(g, *lim)
            if name in fn:
                return fn[name](*[self._ev(a, env, fn) for a in n.args])
        raise SyntaxError("Syntax ERROR")

    def solve(self, text, var="X", guess=1.0, env=None):
        l, _, r = text.partition("=")
        expr = f"({l})-({r or '0'})"

        def f(t):
            try:
                v = float(self.evaluate(expr, {**(env or {}), var: t}))
                return v if math.isfinite(v) else float("nan")
            except Exception:
                return float("nan")

        x = float(guess)
        for _ in range(100):                      # Newton
            fx = f(x)
            if math.isnan(fx):
                break
            if abs(fx) < 1e-13:
                return 0.0 if abs(x) < 1e-12 else x
            h = 1e-6 * max(1.0, abs(x))
            d = (f(x + h) - f(x - h)) / (2 * h)
            if not d or math.isnan(d):
                break
            nx = x - fx / d
            if abs(nx - x) < 1e-14 * max(1.0, abs(x)):
                return nx
            x = nx
        step = 0.25                               # fallback: scan for a sign change
        pts = sorted((guess + i * step for i in range(-4000, 4001)), key=lambda p: abs(p - guess))
        prev = None
        for p in pts:
            fp = f(p)
            if not math.isnan(fp) and abs(fp) < 1e-13:
                return p
        grid = sorted(pts)
        for lo, hi in zip(grid, grid[1:]):
            a_, b_ = f(lo), f(hi)
            if not (math.isnan(a_) or math.isnan(b_)) and a_ * b_ < 0:
                for _ in range(200):
                    mid = (lo + hi) / 2
                    if f(lo) * f(mid) <= 0:
                        hi = mid
                    else:
                        lo = mid
                return (lo + hi) / 2
        raise ValueError("Can't Solve")

    # ---- formatting
    def fmt_real(self, x):
        x = x.real if isinstance(x, complex) else x
        if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
            raise ValueError("Math ERROR")
        if self.frac and abs(x) < 1e10:
            f = Fraction(x).limit_denominator(10000)
            if abs(float(f) - x) <= 1e-10 * max(1, abs(x)):
                n, d = f.numerator, f.denominator
                if d == 1:
                    return str(n)
                if self.frac == "improper" or abs(n) < d:
                    return f"{n}/{d}"
                w, r = divmod(abs(n), d)
                return f"{'-' if n < 0 else ''}{w} {r}/{d}"
        if isinstance(x, int) and self.dm.startswith("norm") and abs(x) < 10 ** 10 \
                and not self.eng_on:
            return str(x)
        x = float(x)
        if self.eng_on:
            return eng(x)
        if self.dm == "fix":
            return f"{x:.{self.dd}f}"
        if self.dm == "sci":
            return sci(x, (self.dd or 10) - 1, strip=False)
        lo = 1e-2 if self.dm == "norm1" else 1e-9
        if x != 0 and (abs(x) < lo or abs(x) >= 1e10):
            return sci(x, 9)
        return plain(x)

    def fmt(self, v):
        if isinstance(v, Matrix):
            return v.text(self.fmt_real)
        if not self.complex_mode:
            if isinstance(v, complex):
                raise ValueError("Math ERROR")
            return self.fmt_real(v)
        z = complex(v)
        f = self.fmt_real
        if self.polar:
            return f"{f(abs(z))}∠{f(self._out(cmath.phase(z)))}"
        re_, im = z.real, z.imag
        if abs(im) < 1e-12 * max(1, abs(re_)):
            return f(re_)
        mag = f(abs(im))
        mag = "" if mag == "1" else mag
        if abs(re_) < 1e-12 * max(1, abs(im)):
            return ("-" if im < 0 else "") + mag + "i"
        return f(re_) + ("−" if im < 0 else "+") + mag + "i"


# ======================================================================
#  Dialog theme
# ======================================================================

BG, PANEL, FG = "#14161a", "#20232a", "#e8eaf0"
C_FN, C_EQ = "#2f6f9f", "#2ea043"


class Tool(tk.Toplevel):
    def __init__(self, app, title):
        super().__init__(app)
        self.app = app
        self.title(title)
        self.configure(bg=BG, padx=10, pady=10)

    def label(self, text):
        tk.Label(self, text=text, bg=BG, fg="#8b949e", justify="left").pack(anchor="w")

    def box(self, w, h, ro=False):
        t = tk.Text(self, width=w, height=h, bg=PANEL, fg=FG, insertbackground=FG, bd=0,
                    font=("Consolas", 11))
        t.pack(pady=4)
        if ro:
            t.config(state="disabled")
        return t

    def button(self, text, cmd, color=C_FN, parent=None):
        b = tk.Button(parent or self, text=text, command=cmd, bg=color, fg=FG, bd=0, padx=10, pady=4)
        b.pack(pady=4, **({"side": "left", "padx": 4} if parent else {}))
        return b

    def write(self, lines):
        self.out.config(state="normal")
        self.out.delete("1.0", "end")
        self.out.insert("end", "\n".join(lines))
        self.out.config(state="disabled")

    def f(self, x):
        try:
            return self.app.eng.fmt_real(x)
        except Exception:
            return "—"

    @staticmethod
    def num(text):
        return float(Engine().evaluate(text or "0"))


class EqnDialog(Tool):
    KINDS = ["Simultaneous – 2 unknowns", "Simultaneous – 3 unknowns",
             "Quadratic  aX² + bX + c = 0", "Cubic  aX³ + bX² + cX + d = 0"]

    def __init__(self, app):
        super().__init__(app, "EQN – equation solver")
        self.kind = ttk.Combobox(self, values=self.KINDS, state="readonly", width=34)
        self.kind.set(self.KINDS[0])
        self.kind.pack(fill="x")
        self.kind.bind("<<ComboboxSelected>>", lambda e: self.build())
        self.grid_f = tk.Frame(self, bg=BG)
        self.grid_f.pack(pady=6)
        self.button("Solve", self.solve, C_EQ)
        self.out = self.box(40, 6, True)
        self.build()

    def build(self):
        for w in self.grid_f.winfo_children():
            w.destroy()
        k = self.KINDS.index(self.kind.get())
        self.cells = []
        if k < 2:
            n, names = k + 2, "xyz"
            for r in range(n):
                row = []
                for c in range(n + 1):
                    e = tk.Entry(self.grid_f, width=6)
                    e.grid(row=r, column=2 * c, padx=2, pady=2)
                    row.append(e)
                    if c < n:
                        tk.Label(self.grid_f, bg=BG, fg=FG,
                                 text=f"{names[c]} {'+' if c < n - 1 else '='}").grid(row=r, column=2 * c + 1)
                self.cells.append(row)
        else:
            row = []
            for c, name in enumerate("abcd"[:k + 1]):
                tk.Label(self.grid_f, text=name, bg=BG, fg=FG).grid(row=0, column=c)
                e = tk.Entry(self.grid_f, width=7)
                e.grid(row=1, column=c, padx=2)
                row.append(e)
            self.cells = [row]

    def solve(self):
        try:
            k = self.KINDS.index(self.kind.get())
            v = [[self.num(e.get()) for e in row] for row in self.cells]
            if k < 2:
                sol = solve_linear([r[:-1] for r in v], [r[-1] for r in v])
                self.write([f"{n} = {self.f(x)}" for n, x in zip("xyz", sol)])
            else:
                roots = poly_roots(v[0])
                self.write([f"X{i + 1} = {self.cf(z)}" for i, z in enumerate(roots)])
        except Exception as ex:
            self.write([f"Error: {ex}"])

    def cf(self, z):
        if not isinstance(z, complex):
            return self.f(z)
        return f"{self.f(z.real)} {'−' if z.imag < 0 else '+'} {self.f(abs(z.imag))}i"


class MatDialog(Tool):
    def __init__(self, app, kind):
        super().__init__(app, "MAT – matrices" if kind == "MAT" else "VCT – vectors")
        self.label("Rows on separate lines, entries separated by spaces." if kind == "MAT"
                   else "Enter each vector as one line, e.g.  1 2 3")
        self.boxes = {}
        sample = {"MAT": ("1 2\n3 4", "5 6\n7 8", ""), "VCT": ("1 2 3", "4 5 6", "")}[kind]
        for name, s in zip("ABC", sample):
            self.label(f"{'Mat' if kind == 'MAT' else 'Vct'}{name}")
            self.boxes[name] = self.box(30, 3 if kind == "MAT" else 1)
            self.boxes[name].insert("1.0", s)
        self.label("Expression: A+B  A-B  A*B  2A  A^-1  A^3  det(A)  trn(A)  inv(A)"
                   if kind == "MAT" else "Expression: dot(A,B)  cross(A,B)  abs(A)  A+B  A-B  2A")
        self.expr = tk.Entry(self, font=("Consolas", 12), width=34)
        self.expr.insert(0, "A*B" if kind == "MAT" else "cross(A,B)")
        self.expr.pack(pady=4)
        self.expr.bind("<Return>", lambda e: self.calc())
        self.button("Calculate", self.calc, C_EQ)
        self.out = self.box(40, 6, True)

    def calc(self):
        try:
            env = {}
            for name, b in self.boxes.items():
                lines = [l for l in re.split(r"[\n;]", b.get("1.0", "end")) if l.strip()]
                if lines:
                    rows = [[self.num(t) for t in re.split(r"[,\s]+", l.strip()) if t] for l in lines]
                    if len({len(r) for r in rows}) != 1:
                        raise ValueError("Dimension ERROR")
                    env[name] = Matrix(rows)
            fn = {"det": lambda m: m.det(), "trn": lambda m: m.T(), "inv": lambda m: m.inv(),
                  "dot": dot, "cross": cross, "norm": lambda m: m.norm()}
            v = self.app.eng.evaluate(self.expr.get(), env, fn)
            self.write(self.app.eng.fmt(v).splitlines())
        except Exception as ex:
            self.write([f"Error: {ex}"])


class ConstDialog(Tool):
    def __init__(self, app):
        super().__init__(app, "CONST – scientific constants")
        self.label("Double-click a constant to insert it")
        self.lb = tk.Listbox(self, width=58, height=20, bg=PANEL, fg=FG, bd=0, font=("Consolas", 10))
        for i, (n, v) in enumerate(CONSTANTS, 1):
            self.lb.insert("end", f"{i:02d}  {n} = {v:.9g}")
        self.lb.pack()
        self.lb.bind("<Double-Button-1>", self.pick)

    def pick(self, _):
        sel = self.lb.curselection()
        if sel:
            self.app.put(lit(CONSTANTS[sel[0]][1]))
            self.destroy()


class ConvDialog(Tool):
    def __init__(self, app):
        super().__init__(app, "CONV – unit conversion")
        self.lb = tk.Listbox(self, width=30, height=20, bg=PANEL, fg=FG, bd=0,
                             font=("Consolas", 11), exportselection=False)
        for i, (a, b, *_) in enumerate(CONVERSIONS, 1):
            self.lb.insert("end", f"{i:02d}  {a}  ⇄  {b}")
        self.lb.pack()
        self.lb.selection_set(0)
        self.val = tk.Entry(self, font=("Consolas", 12))
        self.val.insert(0, app.var.get() or "1")
        self.val.pack(fill="x", pady=4)
        fr = tk.Frame(self, bg=BG)
        fr.pack()
        self.button("→ (left to right)", lambda: self.go(2), parent=fr)
        self.button("← (right to left)", lambda: self.go(3), parent=fr)
        self.res = tk.Label(self, text="", bg=BG, fg="#3fb950", font=("Consolas", 12, "bold"))
        self.res.pack()
        self.button("Insert result", self.insert, C_EQ)
        self.value = None

    def go(self, idx):
        try:
            c = CONVERSIONS[self.lb.curselection()[0]]
            self.value = c[idx](self.num(self.val.get()))
            unit = c[1] if idx == 2 else c[0]
            self.res.config(text=f"{self.app.eng.fmt_real(self.value)} {unit}")
        except Exception as ex:
            self.res.config(text=f"Error: {ex}")

    def insert(self):
        if self.value is not None:
            self.app.put(lit(self.value))
            self.destroy()


class SetupDialog(Tool):
    def __init__(self, app):
        super().__init__(app, "SETUP")
        e = app.eng
        self.ang, self.dm = tk.StringVar(value=e.angle), tk.StringVar(value=e.dm)
        self.dig = tk.IntVar(value=e.dd)
        self.frac = tk.StringVar(value=e.frac or "none")
        self.cpx = tk.StringVar(value="polar" if e.polar else "rect")
        for title, var, opts in (
                ("Angle unit", self.ang, [("Degree", "DEG"), ("Radian", "RAD"), ("Gradian", "GRA")]),
                ("Display", self.dm, [("Fix", "fix"), ("Sci", "sci"), ("Norm 1", "norm1"), ("Norm 2", "norm2")]),
                ("Fractions", self.frac, [("Off", "none"), ("a b/c", "mixed"), ("d/c", "improper")]),
                ("Complex format", self.cpx, [("a+bi", "rect"), ("r∠θ", "polar")])):
            self.label(title)
            fr = tk.Frame(self, bg=BG)
            fr.pack(anchor="w")
            for text, val in opts:
                tk.Radiobutton(fr, text=text, variable=var, value=val, bg=BG, fg=FG,
                               selectcolor=PANEL, activebackground=BG).pack(side="left")
        self.label("Digits for Fix / Sci (0–9)")
        tk.Spinbox(self, from_=0, to=9, textvariable=self.dig, width=4).pack(anchor="w")
        self.button("OK", self.ok, C_EQ)

    def ok(self):
        e = self.app.eng
        e.angle, e.dm, e.dd = self.ang.get(), self.dm.get(), self.dig.get()
        e.frac = None if self.frac.get() == "none" else self.frac.get()
        e.polar = self.cpx.get() == "polar"
        self.app.refresh_keys()
        self.app.redisplay()
        self.destroy()


def numtxt(v):
    return str(int(v)) if float(v).is_integer() else repr(float(v))


class StatDialog(Tool):
    TYPES = ["SD", "Lin", "Log", "Exp", "Pwr", "Inv", "Quad"]

    def __init__(self, app):
        super().__init__(app, "S-SUM / S-VAR  –  statistics")
        self.model, self.vals = None, []
        self.type = ttk.Combobox(self, values=self.TYPES, state="readonly", width=10)
        self.type.set(app.stat_type)
        self.type.pack(anchor="w")
        self.type.bind("<<ComboboxSelected>>", lambda e: self.hint())
        self.hl = tk.Label(self, bg=BG, fg="#8b949e", justify="left")
        self.hl.pack(anchor="w")
        self.data = self.box(34, 8)
        for r in app.stat_rows:
            vals = list(r[:-1]) + ([r[-1]] if r[-1] != 1 else [])
            self.data.insert("end", ", ".join(numtxt(v) for v in vals) + "\n")
        self.hint()
        self.button("Calculate", self.calc, C_EQ)
        self.lb = tk.Listbox(self, width=40, height=14, bg=PANEL, fg=FG, bd=0,
                             font=("Consolas", 11), highlightthickness=0)
        self.lb.pack()
        self.lb.bind("<Double-Button-1>", self.pick)
        fr = tk.Frame(self, bg=BG)
        fr.pack()
        self.ex = tk.Entry(fr, width=8)
        self.ex.pack(side="left")
        self.button("ŷ from x", lambda: self.est(True), parent=fr)
        self.button("x̂ from y", lambda: self.est(False), parent=fr)
        self.label("Double-click a result to insert it into the calculator")

    def hint(self):
        self.hl.config(text="x   or   x, freq" if self.type.get() == "SD" else "x, y   or   x, y, freq")

    def rows(self, kind):
        need = 1 if kind == "SD" else 2
        out = []
        for ln in self.data.get("1.0", "end").splitlines():
            if not ln.strip():
                continue
            v = [self.num(t) for t in re.split(r"[,;\s]+", ln.strip()) if t]
            if len(v) == need:
                v.append(1.0)
            if len(v) != need + 1:
                raise ValueError(f"Bad line: {ln}")
            out.append(tuple(v))
        if not out:
            raise ValueError("No data")
        return out

    def fill(self):
        self.lb.delete(0, "end")
        for k, v in self.vals:
            self.lb.insert("end", f"{k:>9} = {self.f(v)}")

    def calc(self):
        try:
            kind = self.type.get()
            rows = self.rows(kind)
            self.app.stat_rows, self.app.stat_type = rows, kind
            if kind == "SD":
                res, self.model = sd_stats(rows), None
            else:
                res, yh, xh = regression(rows, kind)
                self.model = (yh, xh)
            self.vals = list(res.items())
            self.fill()
        except Exception as ex:
            self.lb.delete(0, "end")
            self.lb.insert("end", f"Error: {ex}")

    def est(self, forward):
        if not self.model:
            return self.calc()
        try:
            v = self.num(self.ex.get())
            r = [self.model[0](v)] if forward else self.model[1](v)
            name = "ŷ" if forward else "x̂"
            self.vals += [(f"{name}{i + 1 if len(r) > 1 else ''}", x) for i, x in enumerate(r)]
            self.fill()
        except Exception as ex:
            self.lb.insert("end", f"Error: {ex}")

    def pick(self, _):
        sel = self.lb.curselection()
        if sel and sel[0] < len(self.vals):
            self.app.put(lit(self.vals[sel[0]][1]))


class DistrDialog(Tool):
    def __init__(self, app):
        super().__init__(app, "DISTR  –  normal distribution")
        self.ents = []
        for name, d in (("x", "1"), ("μ (mean)", "0"), ("σ (std dev)", "1")):
            fr = tk.Frame(self, bg=BG)
            fr.pack(anchor="w", pady=2)
            tk.Label(fr, text=name, bg=BG, fg=FG, width=12, anchor="w").pack(side="left")
            e = tk.Entry(fr, width=14)
            e.insert(0, d)
            e.pack(side="left")
            self.ents.append(e)
        fr = tk.Frame(self, bg=BG)
        fr.pack()
        self.button("Calculate", self.go, C_EQ, parent=fr)
        self.button("Use data mean / σn", self.use_data, parent=fr)
        self.lb = tk.Listbox(self, width=30, height=5, bg=PANEL, fg=FG, bd=0,
                             font=("Consolas", 11), highlightthickness=0)
        self.lb.pack(pady=4)
        self.lb.bind("<Double-Button-1>", self.pick)
        self.vals = []

    def use_data(self):
        rows = self.app.stat_rows
        if rows and len(rows[0]) == 2:
            s = sd_stats(rows)
            for e, v in zip(self.ents[1:], (s["x̄ (mean)"], s["σn"])):
                e.delete(0, "end")
                e.insert(0, numtxt(v))

    def go(self):
        try:
            x, mu, sg = (self.num(e.get()) for e in self.ents)
            t = (x - mu) / sg
            p = 0.5 * (1 + math.erf(t / math.sqrt(2)))
            self.vals = [("t", t), ("P(t)", p), ("Q(t)", abs(p - 0.5)), ("R(t)", 1 - p)]
            self.lb.delete(0, "end")
            for k, v in self.vals:
                self.lb.insert("end", f"{k:>5} = {self.f(v)}")
        except Exception as ex:
            self.lb.delete(0, "end")
            self.lb.insert("end", f"Error: {ex}")

    def pick(self, _):
        sel = self.lb.curselection()
        if sel and sel[0] < len(self.vals):
            self.app.put(lit(self.vals[sel[0]][1]))


class VarDialog(Tool):
    """CALC (evaluate with variable values) and SOLVE (find a variable)."""

    def __init__(self, app, expr, solve):
        super().__init__(app, "SOLVE" if solve else "CALC")
        self.expr, self.solve_mode = expr, solve
        names = sorted(set(re.findall(r"(?<![A-Za-z])[A-FXYM](?![A-Za-z])", expr)))
        if solve and not names:
            names = ["X"]
        self.label(expr)
        self.target = None
        if solve:
            self.label("Solve for:")
            self.target = ttk.Combobox(self, values=names, state="readonly", width=4)
            self.target.set("X" if "X" in names else names[0])
            self.target.pack(anchor="w")
            self.label("Initial values (for the solve variable: starting guess)")
        self.ents = {}
        for n in names:
            fr = tk.Frame(self, bg=BG)
            fr.pack(anchor="w", pady=2)
            tk.Label(fr, text=f"{n} =", bg=BG, fg=FG, width=4).pack(side="left")
            e = tk.Entry(fr, width=16)
            e.insert(0, numtxt(app.eng.vars[n]) if isinstance(app.eng.vars[n], (int, float)) else "0")
            e.pack(side="left")
            self.ents[n] = e
        self.button("Solve" if solve else "Calculate", self.go, C_EQ)
        self.res = tk.Label(self, text="", bg=BG, fg="#3fb950", font=("Consolas", 12, "bold"))
        self.res.pack()

    def go(self):
        eng = self.app.eng
        try:
            env = {n: self.num(e.get()) for n, e in self.ents.items()}
            if self.solve_mode:
                t = self.target.get()
                root = eng.solve(self.expr, t, env.pop(t, 1.0), env)
                env[t] = root
                msg = f"{t} = {eng.fmt_real(root)}"
            else:
                v = eng.evaluate(self.expr, env)
                msg = eng.fmt(v)
                eng.ans = v
            for n, v in env.items():
                eng.vars[n] = v
            self.res.config(text=msg)
            self.app.show(msg.split("=")[-1].strip() if self.solve_mode else msg,
                          sub=msg if self.solve_mode else "")
        except Exception as ex:
            self.res.config(text=f"Error: {ex}", fg="#f85149")


class HelpDialog(Tool):
    TEXT = (
        "KEY LAYOUT\n"
        "  SHIFT (yellow) and ALPHA (pink) select the printed labels.\n"
        "  Green labels work directly in BASE-N mode (DEC HEX BIN OCT, LOGIC).\n"
        "  Blue labels are mode functions: SHIFT+key (CMPLX: ∠ arg Abs Conjg ▶r∠θ ▶a+bi Re⇔Im),\n"
        "  SHIFT+4 MAT, SHIFT+5 VCT, SHIFT+1 S-SUM, SHIFT+2 S-VAR, SHIFT+3 DISTR,\n"
        "  M+ = DT in SD/REG mode (x;freq  or  x,y;freq), SHIFT+M+ = CL (clear data).\n\n"
        "CALC / SOLVE\n"
        "  Type an expression with variables (A–F, X, Y, M) and press CALC.\n"
        "  SOLVE: use ALPHA+CALC for '=', e.g.  X^2-2=0   then SHIFT+CALC.\n\n"
        "ENGINEERING PREFIXES (ALPHA+digit)\n"
        "  7 M  8 G  9 T  4 µ  5 m  6 k  1 f  2 p  3 n\n\n"
        "OTHER\n"
        "  ▲ ▼ replay earlier entries, ◀ ▶ move cursor, SHIFT+▲ copies the result.\n"
        "  SHIFT+DEL = INS (overwrite), SHIFT+AC = OFF, ON wakes the calculator.\n"
        "  ∫(f(X),a,b)   d/dx(f(X),a)   use the variable X.\n"
        "  MODE menu: COMP CMPLX SD REG BASE-N EQN MAT VCT and SETUP.")

    def __init__(self, app):
        super().__init__(app, "Key guide")
        t = self.box(88, 24)
        t.insert("1.0", self.TEXT)
        t.config(state="disabled")


# ======================================================================
#  Premium calculator window  (layout traced from the real fx-991MS)
# ======================================================================

WBG = "#0d0e11"
LCD, LCDT, LCDD = "#c6d0b3", "#1c2517", "#74826a"
KEY, KEY_HI = "#2f333b", "#3c414b"
COL = {"shift": "#f5c542", "alpha": "#ff7ab8", "green": "#6fd88a", "blue": "#69b7ff",
       "main": "#eef1f6"}
HYP = {"sin(": "sinh(", "cos(": "cosh(", "tan(": "tanh(",
       "asin(": "asinh(", "acos(": "acosh(", "atan(": "atanh("}
BASE_NAMES = {10: "DEC", 16: "HEX", 2: "BIN", 8: "OCT"}
DEL_TOKEN = re.compile(r"(?:a?(?:sin|cos|tan)h?|log|ln|exp|sqrt|cbrt|root|abs|int|intg|rnd|pol|"
                       r"rec|conjg|rep|imp|arg|integ|deriv|fact|not|neg|cpolar)\($|ran\(\)$|"
                       r"d/dx\($|\s(?:nCr|nPr|and|or|xor|xnor)\s$")
CMPLX, STAT = {"CMPLX"}, {"SD", "REG"}


def lerp(c1, c2, t):
    a = [int(c1[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(c2[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02x%02x%02x" % tuple(round(x + (y - x) * t) for x, y in zip(a, b))


DIM = {k: lerp(v, "#1a1c20", 0.62) for k, v in COL.items()}


def rrect(cv, x1, y1, x2, y2, r, **kw):
    pts = [x1 + r, y1, x2 - r, y1, x2, y1, x2, y1 + r, x2, y2 - r, x2, y2, x2 - r, y2,
           x1 + r, y2, x1, y2, x1, y2 - r, x1, y1 + r, x1, y1]
    return cv.create_polygon(pts, smooth=True, **kw)


def grad_rrect(cv, x1, y1, x2, y2, r, c1, c2):
    h = int(y2 - y1)
    for i in range(h):
        d = min(i, h - i)
        inset = 0 if d >= r else r - math.sqrt(max(0.0, r * r - (r - d) ** 2))
        cv.create_line(x1 + inset, y1 + i, x2 - inset, y1 + i, fill=lerp(c1, c2, i / h))


def k(label, act, sh=None, al=None, gr=None, mm=(), **kw):
    return dict(label=label, act=act, sh=sh, al=al, gr=gr, mm=list(mm), **kw)


DIG = {"size": 26}
ROWS = [  # (y, height, width, labels-below?, [(x, key), ...])
    (783, 46, 82, False, [
        (222, k("CALC", "@calc", ("SOLVE", "@solve"), ("=", "="), size=13)),
        (317, k("∫dx", "∫(", ("d/dx", "d/dx("), (":", ":"))),
        (602, k("x⁻¹", "^(-1)", ("x!", "!"), gr=("LOGIC", "@logic"))),
        (697, k("CONST", "@const", ("CONV", "@conv"), size=12))]),
    (865, 46, 82, False, [
        (222, k("a b/c", "⌟", ("d/c", "@dc"), size=13)),
        (317, k("√", "√(", ("³√", "cbrt("), size=18)),
        (412, k("x²", "^2", ("x³", "^3"), gr=("DEC", "@b10"))),
        (508, k("^", "^", ("ˣ√", "ˣ√"), gr=("HEX", "@b16"), size=18)),
        (602, k("log", "log(", ("10ˣ", "10^("), gr=("BIN", "@b2"))),
        (697, k("ln", "ln(", ("eˣ", "exp("), ("e", "e"), gr=("OCT", "@b8")))]),
    (945, 46, 82, False, [
        (222, k("(−)", "−", al=("A", "A"), mm=[("∠", "∠", CMPLX, True)])),
        (317, k("° ′ ″", "@dmsch", ("←", "@dms"), ("B", "B"), size=13)),
        (412, k("hyp", "@hyp", al=("C", "C"))),
        (508, k("sin", "sin(", ("sin⁻¹", "asin("), ("D", "D"))),
        (602, k("cos", "cos(", ("cos⁻¹", "acos("), ("E", "E"))),
        (697, k("tan", "tan(", ("tan⁻¹", "atan("), ("F", "F")))]),
    (1027, 46, 82, True, [
        (222, k("RCL", "@rcl", ("STO", "@sto"))),
        (317, k("ENG", "@eng", ("←", "@eng"), ("i", "i"))),
        (412, k("(", "(", mm=[("arg", "arg(", CMPLX, True)], size=18)),
        (508, k(")", ")", (",", ","), ("X", "X"), mm=[("Abs", "abs(", CMPLX, True)], size=18)),
        (602, k(",", ",", (";", ";"), ("Y", "Y"), mm=[("Conjg", "conjg(", CMPLX, True)], size=18)),
        (697, k("M+", "@mplus", ("M−", "@mminus"), ("M", "M"),
                mm=[("DT", "@dt", STAT, False), ("CL", "@cl", STAT, True)]))]),
    (1154, 60, 94, False, [
        (229, k("7", "7", al=("M", "ᴇ6"), **DIG)), (343, k("8", "8", al=("G", "ᴇ9"), **DIG)),
        (458, k("9", "9", al=("T", "ᴇ12"), **DIG)),
        (572, k("DEL", "@del", ("INS", "@ins"), fg="#ff4f9a", size=17)),
        (687, k("AC", "@ac", ("OFF", "@off"), fg="#ff4545", size=17))]),
    (1256, 60, 94, False, [
        (229, k("4", "4", al=("μ", "ᴇ-6"), mm=[("MAT", "@mat", None, True)], **DIG)),
        (343, k("5", "5", al=("m", "ᴇ-3"), mm=[("VCT", "@vct", None, True)], **DIG)),
        (458, k("6", "6", al=("k", "ᴇ3"), **DIG)),
        (572, k("×", "×", ("nPr", " nPr "), size=26)), (687, k("÷", "÷", ("nCr", " nCr "), size=26))]),
    (1357, 60, 94, False, [
        (229, k("1", "1", al=("f", "ᴇ-15"), mm=[("S-SUM", "@ssum", None, True)], **DIG)),
        (343, k("2", "2", al=("p", "ᴇ-12"), mm=[("S-VAR", "@svar", None, True)], **DIG)),
        (458, k("3", "3", al=("n", "ᴇ-9"), mm=[("DISTR", "@distr", None, True)], **DIG)),
        (572, k("+", "+", ("Pol(", "pol("), mm=[("▶r∠θ", "@polar", CMPLX, True)], size=26)),
        (687, k("−", "−", ("Rec(", "rec("), mm=[("▶a+bi", "@rect", CMPLX, True)], size=26))]),
    (1451, 60, 94, False, [
        (229, k("0", "0", ("Rnd", "rnd("), **DIG)), (343, k(".", ".", ("Ran#", "ran()"), **DIG)),
        (458, k("×10ˣ", "ᴇ", ("π", "π"), size=15)),
        (572, k("Ans", "Ans", ("DRG▶", "@drgmenu"), size=17)),
        (687, k("=", "@eq", al=("%", "%"), mm=[("Re⇔Im", "@reim", CMPLX, True)], size=26))]),
]


class App(tk.Tk):
    def __init__(self, scale=None):
        super().__init__()
        self.title("TILLIDON fx-991MS  ·  Scientific Calculator")
        self.configure(bg=WBG)
        self.minsize(520, 640)
        s0 = scale or max(0.42, min(0.66, (self.winfo_screenheight() - 170) / 1411))
        self.S, self.ox, self.oy, self._size = s0, 0.0, 0.0, None
        self._res, self._rd, self.fullscreen = ("0", None, ""), None, False
        self.eng = Engine()
        self.mode, self.base = "COMP", 10
        self.shift = self.alpha = self.hyp = self.fresh = self.overwrite = False
        self.on, self.pending = True, None
        self.expr_hist, self.hpos = [], 0
        self.stat_rows, self.stat_type = [], "SD"
        self.specs, self.labels, self._saved = {}, [], {}
        self.var = tk.StringVar()
        self.build_menu()
        self.build(s0)
        self.update_idletasks()
        self.redraw()
        self.set_mode("COMP")
        self.refresh()

    # ------------------------------------------------------------ geometry
    def X(self, x):
        return self.ox + x * self.S

    def Y(self, y):
        return self.oy + y * self.S

    def fs(self, n):
        return max(8, round(n * self.S / 0.62))

    def font(self, n, w="bold", fam="Segoe UI"):
        return (fam, -self.fs(n), w)

    # ------------------------------------------------------------ menu bar
    def build_menu(self):
        mb = tk.Menu(self)
        m = tk.Menu(mb, tearoff=0)
        for label, code in (("COMP", "COMP"), ("CMPLX", "CMPLX"), ("SD", "SD"), ("REG", "REG"),
                            ("BASE-N", "BASE"), ("EQN", "EQN"), ("MAT", "MAT"), ("VCT", "VCT")):
            m.add_command(label=label, command=lambda c=code: self.set_mode(c))
        m.add_separator()
        m.add_command(label="SETUP…", command=self.c_setup)
        mb.add_cascade(label="Mode", menu=m)
        t = tk.Menu(mb, tearoff=0)
        for label, cmd in (("Constants (CONST)", self.c_const), ("Conversions (CONV)", self.c_conv),
                           ("Statistics (S-SUM / S-VAR)", self.c_ssum), ("Distribution (DISTR)", self.c_distr),
                           ("Matrix (MAT)", self.c_mat), ("Vector (VCT)", self.c_vct),
                           ("Equations (EQN)", lambda: EqnDialog(self))):
            t.add_command(label=label, command=cmd)
        mb.add_cascade(label="Tools", menu=t)
        f = tk.Menu(mb, tearoff=0)
        for label, tok in (("Abs", "abs("), ("Int", "int("), ("Intg", "intg("), ("Rnd", "rnd("),
                           ("Ran#", "ran()"), ("Conjg", "conjg("), ("arg", "arg("), ("ReP", "rep("),
                           ("ImP", "imp("), ("mod  (a % b)", "%")):
            f.add_command(label=label, command=lambda x=tok: self.put(x))
        mb.add_cascade(label="Functions", menu=f)
        v = tk.Menu(mb, tearoff=0)
        v.add_command(label="Full screen  (F11)", command=self.toggle_full)
        v.add_command(label="Show / hide history", command=self.toggle_side)
        mb.add_cascade(label="View", menu=v)
        h = tk.Menu(mb, tearoff=0)
        h.add_command(label="Key guide", command=lambda: HelpDialog(self))
        mb.add_cascade(label="Help", menu=h)
        self.config(menu=mb)

    # ------------------------------------------------------------ drawing
    def build(self, s0):
        wrap = tk.Frame(self, bg=WBG)
        wrap.pack(fill="both", expand=True, padx=14, pady=10)
        self.build_side(wrap)
        self.cv = tk.Canvas(wrap, width=int(655 * s0 + 28), height=int(1411 * s0 + 28), bg=WBG,
                            highlightthickness=0)
        self.cv.pack(side="left", fill="both", expand=True)
        self.entry = tk.Entry(self.cv, textvariable=self.var, bd=0, highlightthickness=0,
                              relief="flat", bg=LCD, fg=LCDT, insertbackground=LCDT, justify="left")
        self.entry.bind("<Key>", self.on_key)
        self.entry.bind("<Up>", lambda e: self.c_up())
        self.entry.bind("<Down>", lambda e: self.c_down())
        self.bind("<Return>", lambda e: self._do("@eq"))
        self.bind("<KP_Enter>", lambda e: self._do("@eq"))
        self.bind("<Escape>", lambda e: self._do("@ac"))
        self.bind("<F11>", lambda e: self.toggle_full())
        self.var.trace_add("write", self.preview)
        self.cv.bind("<Configure>", self.on_resize)
        self.entry.focus_set()

    def toggle_full(self):
        self.fullscreen = not self.fullscreen
        self.attributes("-fullscreen", self.fullscreen)

    def on_resize(self, _):
        if self._rd:
            self.after_cancel(self._rd)
        self._rd = self.after(60, self.redraw)

    def redraw(self):
        """(Re)draw the whole calculator so it fills the window at any size."""
        self._rd = None
        cv = self.cv
        w, h = cv.winfo_width(), cv.winfo_height()
        if w < 50 or h < 50:
            w, h = int(cv["width"]), int(cv["height"])
        if (w, h) == self._size:
            return
        self._size = (w, h)
        self.S = max(0.3, min(1.6, (w - 28) / 655, (h - 28) / 1411))
        self.ox = (w - 655 * self.S) / 2 - 135 * self.S
        self.oy = (h - 1411 * self.S) / 2 - 125 * self.S
        cv.delete("all")
        self.specs, self.labels, self._saved = {}, [], {}
        self.draw_body()
        self.draw_lcd()
        for y, hh, ww, below, keys in ROWS:
            for x, spec in keys:
                self.draw_key(spec, x, y, ww, hh, below)
        self.draw_round_keys()
        self.draw_dpad()
        self.refresh()
        self._paint(*self._res)
        if not self.on:
            self.c_off()

    def draw_body(self):
        cv, S = self.cv, self.S
        x1, y1, x2, y2 = self.X(135), self.Y(125), self.X(790), self.Y(1536)
        r = 48 * S
        for d, c in ((12, "#08090b"), (8, "#0b0c0f"), (4, "#101216")):
            rrect(cv, x1 - d + 5, y1 - d + 9, x2 + d + 5, y2 + d + 9, r + d, fill=c, outline="")
        grad_rrect(cv, x1, y1, x2, y2, r, "#383c44", "#15171b")
        rrect(cv, x1, y1, x2, y2, r, fill="", outline="#4d525c")
        brand = ("Arial", -self.fs(46), "bold")
        cv.create_text(self.X(205) + 1.5, self.Y(188) + 2, text="TILLIDON", anchor="w",
                       fill="#050506", font=brand)
        cv.create_text(self.X(205), self.Y(188), text="TILLIDON", anchor="w", fill="#dde0e6",
                       font=brand)
        cv.create_text(self.X(208), self.Y(268), text="fx-991MS", anchor="w", fill="#aeb3bc",
                       font=("Segoe UI", -self.fs(24), "bold italic"))
        cv.create_line(self.X(205), self.Y(306), self.X(720), self.Y(306), fill="#41464f")

    def draw_lcd(self):
        cv, S = self.cv, self.S
        rrect(cv, self.X(176), self.Y(346), self.X(742), self.Y(540), 12 * S, fill="#08090b",
              outline="#2b2e35")
        sx1, sy1, sx2, sy2 = self.X(188), self.Y(358), self.X(730), self.Y(528)
        self.sx1, self.sy1, self.sx2, self.sy2 = sx1, sy1, sx2, sy2
        rrect(cv, sx1, sy1, sx2, sy2, 6 * S, fill=LCD, outline="#8a957a")
        cv.create_line(sx1 + 6, sy1 + 2, sx2 - 6, sy1 + 2, fill="#919c82", width=2)
        mono = ("Consolas", -self.fs(12), "bold")
        self.t_stat = cv.create_text(sx1 + 10, sy1 + 11, anchor="w", text="", font=mono, fill=LCDD)
        self.t_ind = cv.create_text(sx2 - 10, sy1 + 11, anchor="e", text="", font=mono, fill=LCDT)
        cv.tag_bind(self.t_ind, "<Button-1>", lambda e: self.cycle_angle())
        eh = int(self.fs(24) * 1.45)
        self.entry.config(font=("Consolas", -self.fs(24)))
        self.ent_win = cv.create_window(sx1 + 8, sy1 + 22, window=self.entry, anchor="nw",
                                        width=int(sx2 - sx1 - 16), height=eh)
        self.t_sub = cv.create_text(sx2 - 10, sy2 - self.fs(34) - 8, anchor="se", text="",
                                    font=("Consolas", -self.fs(12)), fill=LCDD)
        self.t_res = cv.create_text(sx2 - 10, sy2 - 4, anchor="se", text="0", fill=LCDT,
                                    font=("Consolas", -self.fs(34), "bold"))

    def bind_tag(self, tag):
        cv = self.cv
        cv.tag_bind(tag, "<ButtonPress-1>", lambda e, t=tag: self._press(t))
        cv.tag_bind(tag, "<ButtonRelease-1>", lambda e, t=tag: self._release(t, True))
        cv.tag_bind(tag, "<Enter>", lambda e: cv.config(cursor="hand2"))
        cv.tag_bind(tag, "<Leave>", lambda e, t=tag: (self._release(t, False), cv.config(cursor="")))

    def _press(self, tag):
        saved = {}
        for it in self.cv.find_withtag(tag):
            c = self.cv.itemcget(it, "fill")
            if isinstance(c, str) and len(c) == 7 and c.startswith("#"):
                saved[it] = c
                self.cv.itemconfig(it, fill=lerp(c, "#000000", 0.28))
        self._saved[tag] = saved

    def _release(self, tag, fire):
        saved = self._saved.pop(tag, None)
        if saved is None:
            return
        for it, c in saved.items():
            self.cv.itemconfig(it, fill=c)
        if fire:
            self._key(self.specs[tag])

    def place_labels(self, items, cx, y, limit, size=11):
        """Lay labels out side by side, shrinking the font until they fit in `limit` px."""
        if not items:
            return
        cv = self.cv
        for _ in range(4):
            ids = [cv.create_text(0, y, text=t, fill=COL[kd], font=self.font(size))
                   for t, kd, br in items]
            ws = []
            for it, (t, kd, br) in zip(ids, items):
                bb = cv.bbox(it)
                ws.append(bb[2] - bb[0] + (6 if br else 0))
            gap = max(8.0, self.fs(size) * 0.8)
            total = sum(ws) + gap * (len(items) - 1)
            if total <= limit or self.fs(size) <= 8:
                break
            for it in ids:
                cv.delete(it)
            size -= 1
        x = cx - total / 2
        for it, (t, kd, br), wd in zip(ids, items, ws):
            cv.coords(it, x + wd / 2, y)
            if kd != "main":
                self.labels.append((it, kd))
            if br:
                bb = cv.bbox(it)
                ln = cv.create_line(bb[0] - 2, bb[1] + 2, bb[0] - 2, bb[1] - 2, bb[2] + 2,
                                    bb[1] - 2, bb[2] + 2, bb[1] + 2, fill=COL[kd])
                self.labels.append((ln, kd))
            x += wd + gap

    def draw_labels(self, spec, cx, w, y1, y2, below):
        top = [(spec["sh"][0], "shift", False)] if spec["sh"] else []
        mm = [(m[0], "blue", True) for m in spec["mm"]]
        tail = []
        if spec["gr"]:
            tail.append((spec["gr"][0], "green", False))
        if spec["al"]:
            tail.append((spec["al"][0], "alpha", False))
        limit = w + 12 * self.S
        hy = self.fs(11) * 0.5
        if below:
            self.place_labels(top + tail, cx, y1 - hy - 7, limit)
            self.place_labels(mm, cx, y2 + hy + 6, limit)
        else:
            self.place_labels(top + mm + tail, cx, y1 - hy - 7, limit)

    def draw_key(self, spec, cx, cy, w, h, below=False):
        cv, S = self.cv, self.S
        tag = f"k{len(self.specs)}"
        self.specs[tag] = spec
        x1, x2 = self.X(cx - w / 2), self.X(cx + w / 2)
        y1, y2 = self.Y(cy - h / 2), self.Y(cy + h / 2)
        r = 11 * S
        rrect(cv, x1 + 1, y1 + 3, x2 + 1, y2 + 3, r, fill="#050506", outline="", tags=tag)
        rrect(cv, x1, y1, x2, y2, r, fill=KEY, outline="#0d0e11", tags=tag)
        rrect(cv, x1 + 2, y1 + 2, x2 - 2, y1 + (y2 - y1) * .52, r * .8, fill=KEY_HI, outline="",
              tags=tag)
        cv.create_text((x1 + x2) / 2, (y1 + y2) / 2 + 1, text=spec["label"],
                       fill=spec.get("fg", COL["main"]), font=self.font(spec.get("size", 15)),
                       tags=tag)
        self.bind_tag(tag)
        self.draw_labels(spec, (x1 + x2) / 2, x2 - x1, y1, y2, below)

    def silver(self, cx, cy, r, tag):
        cv = self.cv
        cv.create_oval(cx - r - 1, cy - r + 3, cx + r + 1, cy + r + 4, fill="#050506", outline="",
                       tags=tag)
        n = 9
        for i in range(n):
            t = i / (n - 1)
            rr = r * (1 - 0.78 * t)
            ox, oy = -r * 0.10 * t, -r * 0.14 * t
            cv.create_oval(cx + ox - rr, cy + oy - rr, cx + ox + rr, cy + oy + rr,
                           fill=lerp("#767c85", "#eceef2", t), outline="", tags=tag)

    def draw_round_keys(self):
        cv, S = self.cv, self.S
        specs = {"shift": k("SHIFT", "@shift"), "alpha": k("ALPHA", "@alpha"),
                 "mode": k("MODE", "@mode", ("CLR", "@clrmenu")), "on": k("ON", "@on")}
        for name, cx in (("shift", 220), ("alpha", 309), ("mode", 612), ("on", 697)):
            tag = f"k{len(self.specs)}"
            self.specs[tag] = specs[name]
            self.silver(self.X(cx), self.Y(665), 28 * S, tag)
            self.bind_tag(tag)
            if name in ("shift", "alpha"):
                ring = cv.create_oval(self.X(cx) - 31 * S, self.Y(665) - 31 * S,
                                      self.X(cx) + 31 * S, self.Y(665) + 31 * S,
                                      outline="", width=3)
                setattr(self, "ring_" + name, ring)
        for cx, items in ((220, [("SHIFT", "shift", False)]), (309, [("ALPHA", "alpha", False)]),
                          (604, [("MODE", "main", False), ("CLR", "shift", False)]),
                          (700, [("ON", "main", False)])):
            self.place_labels(items, self.X(cx), self.Y(621), 220 * S, 14)

    def draw_dpad(self):
        cv, S = self.cv, self.S
        cx, cy = 460, 712
        cv.create_oval(self.X(cx - 108), self.Y(cy - 108), self.X(cx + 108), self.Y(cy + 108),
                       fill="#1a1c20", outline="#33363d", width=2)
        dr = 1.5 * S
        for gy in range(-100, 101, 13):
            for gx in range(-100, 101, 15):
                x = gx + (7 if (gy // 13) % 2 else 0)
                if x * x + gy * gy < 96 * 96:
                    cv.create_oval(self.X(cx + x) - dr, self.Y(cy + gy) - dr,
                                   self.X(cx + x) + dr, self.Y(cy + gy) + dr, fill="#2b2e34", outline="")
        for dx, dy, kind, spec in ((0, -44, "up", k("▲", "@up", ("COPY", "@copy"))),
                                   (0, 46, "down", k("▼", "@down")),
                                   (-68, 0, "left", k("◀", "@left")),
                                   (68, 0, "right", k("▶", "@right"))):
            tag = f"k{len(self.specs)}"
            self.specs[tag] = spec
            kx, ky = cx + dx, cy + dy
            x1, x2, y1, y2 = self.X(kx - 32), self.X(kx + 32), self.Y(ky - 28), self.Y(ky + 28)
            rrect(cv, x1 + 1, y1 + 3, x2 + 1, y2 + 3, 12 * S, fill="#050506", outline="", tags=tag)
            rrect(cv, x1, y1, x2, y2, 12 * S, fill="#8d939c", outline="#585d66", tags=tag)
            rrect(cv, x1 + 2, y1 + 2, x2 - 2, y1 + (y2 - y1) * .5, 10 * S, fill="#d3d7dd",
                  outline="", tags=tag)
            tri = {"up": [(0, -9), (-11, 7), (11, 7)], "down": [(0, 9), (-11, -7), (11, -7)],
                   "left": [(-9, 0), (7, -11), (7, 11)], "right": [(9, 0), (-7, -11), (-7, 11)]}[kind]
            pts = [c for px, py in tri for c in (self.X(kx + px), self.Y(ky + py))]
            cv.create_polygon(pts, fill="#2a2d33", outline="", tags=tag)
            self.bind_tag(tag)

        cp = cv.create_text(self.X(460), self.Y(583), text="COPY", fill=COL["shift"],
                            font=self.font(14))
        self.labels.append((cp, "shift"))

    def build_side(self, wrap):
        self.side = tk.Frame(wrap, bg="#16181d", padx=10, pady=10)
        self.side.pack(side="right", fill="y", padx=(12, 0))
        tk.Label(self.side, text="HISTORY", bg="#16181d", fg="#9aa3b2",
                 font=("Segoe UI", 10, "bold")).pack(anchor="w")
        self.hist = tk.Listbox(self.side, width=30, bg="#1d2026", fg="#e6e9ef", bd=0,
                               font=("Consolas", 10), highlightthickness=0,
                               selectbackground="#3a3f4b", activestyle="none")
        self.hist.pack(fill="both", expand=True, pady=6)
        self.hist.bind("<Double-Button-1>", self.recall)
        tk.Button(self.side, text="Clear history", bg="#c9433f", fg="#fff", bd=0,
                  command=lambda: self.hist.delete(0, "end")).pack(fill="x")
        tk.Label(self.side, text="Double-click an entry to reuse it", bg="#16181d", fg="#6b7280",
                 font=("Segoe UI", 8)).pack(pady=(6, 0))

    def toggle_side(self):
        if self.side.winfo_ismapped():
            self.side.pack_forget()
        else:
            self.side.pack(side="right", fill="y", padx=(12, 0), before=self.cv)

    # ------------------------------------------------------------ display
    def refresh(self):
        cv, e = self.cv, self.eng
        act = self.shift or self.alpha
        for it, kind in self.labels:
            lit_ = (not act or (kind == "shift" and self.shift) or (kind == "alpha" and self.alpha)
                    or (kind == "blue" and self.shift))
            cv.itemconfig(it, fill=COL[kind] if lit_ else DIM[kind])
        cv.itemconfig(self.ring_shift, outline=COL["shift"] if self.shift else "")
        cv.itemconfig(self.ring_alpha, outline=COL["alpha"] if self.alpha else "")
        flags = [t for t, on in (("S", self.shift), ("A", self.alpha), ("HYP", self.hyp),
                                 (self.pending or "", bool(self.pending)),
                                 ("M", e.vars["M"] != 0), ("ENG", e.eng_on), ("FRAC", bool(e.frac)),
                                 ("OVR", self.overwrite)) if on]
        disp = {"fix": f"Fix{e.dd}", "sci": f"Sci{e.dd}", "norm1": "Norm1", "norm2": "Norm2"}[e.dm]
        mode = self.mode + (f" {BASE_NAMES[self.base]}" if self.mode == "BASE" else "")
        cv.itemconfig(self.t_stat, text="  ".join(flags))
        cv.itemconfig(self.t_ind, text=f"{e.angle[0]}  {mode}  {disp}")

    refresh_keys = refresh

    def _paint(self, text, color, sub):
        n = len(text)
        size = 34 if n <= 12 else max(15, int(34 * 12 / n))
        self.cv.itemconfig(self.t_res, text=text, fill=color or LCDT,
                           font=("Consolas", -self.fs(size), "bold"))
        self.cv.itemconfig(self.t_sub, text=sub)

    def show(self, text, color=None, sub=""):
        self._res = (text, color, sub)
        self._paint(*self._res)

    def redisplay(self):
        try:
            if self.mode == "BASE":
                self.show(base_str(self.eng.ans, self.base), sub=self.base_line(int(self.eng.ans)))
            else:
                self.show(self.eng.fmt(self.eng.ans))
        except Exception:
            pass

    @staticmethod
    def base_line(v):
        return "  ".join(f"{BASE_NAMES[b]} {base_str(v, b)}" for b in (10, 16, 2, 8))

    def preview(self, *_):
        s = self.var.get().strip()
        if not s or self.mode == "BASE" or re.search(r"integ|deriv|∫|d/dx|ran|pol\(|rec\(|=", s):
            return
        try:
            self.show(self.eng.fmt(self.eng.evaluate(s)), LCDD)
        except Exception:
            pass

    def error(self, ex):
        if isinstance(ex, ValueError) and "ERROR" in str(ex):
            msg = str(ex)
        elif isinstance(ex, (SyntaxError, KeyError)):
            msg = "Syntax ERROR"
        elif isinstance(ex, ZeroDivisionError):
            msg = "Math ERROR (÷0)"
        else:
            msg = "Math ERROR"
        self.show(msg, "#8a1c1c")

    # ------------------------------------------------------------ input
    def _key(self, spec):
        if not self.on and spec["act"] != "@on":
            return
        sh, al, act = self.shift, self.alpha, None
        for lab, a, modes, need in spec["mm"]:
            if (modes is None or self.mode in modes) and need == sh and not al:
                act = a
                break
        if act is None:
            if sh and spec["sh"]:
                act = spec["sh"][1]
            elif al and spec["al"]:
                act = spec["al"][1]
            elif not sh and not al and self.mode == "BASE" and spec["gr"]:
                act = spec["gr"][1]
            else:
                act = spec["act"]
        if act not in ("@shift", "@alpha"):
            self.shift = self.alpha = False
        self._do(act)
        self.refresh()
        self.entry.focus_set()

    def _do(self, a):
        if not a:
            return
        if a.startswith("@"):
            getattr(self, "c_" + a[1:])()
            self.refresh()
            return
        if self.hyp and a in HYP:
            a, self.hyp = HYP[a], False
        if self.pending and len(a) == 1 and a in "ABCDEFXYM":
            p, self.pending = self.pending, None
            return self.store(a) if p == "STO" else self.put(lit(self.eng.vars[a]))
        self.pending = None
        self.put(a)

    def put(self, tok):
        if not self.on:
            return
        if self.fresh:
            self.fresh = False
            prefix = ""
            if tok[0] in "+−×÷^!%":
                prefix = "Ans" if self.mode != "BASE" else base_str(self.eng.ans, self.base)
            self.entry.delete(0, "end")
            self.entry.insert(0, prefix)
            self.entry.icursor("end")
        elif self.overwrite:
            pos = self.entry.index("insert")
            self.entry.delete(pos, pos + len(tok))
        self.entry.insert("insert", tok)
        self.entry.focus_set()

    def on_key(self, e):
        ch = e.char
        if ch and ch.isprintable():
            self._do({"*": "×", "/": "÷", "-": "−"}.get(ch, ch))
            return "break"

    def recall(self, _):
        sel = self.hist.curselection()
        if sel:
            self.var.set(self.hist.get(sel[0]).split("  =  ")[0])
            self.entry.icursor("end")
            self.fresh = False

    def popup(self, menu):
        x, y = self.winfo_pointerxy()
        menu.tk_popup(x, y)

    # ------------------------------------------------------------ commands
    def c_shift(self):
        self.shift, self.alpha = not self.shift, False

    def c_alpha(self):
        self.alpha, self.shift = not self.alpha, False

    def c_hyp(self):
        self.hyp = not self.hyp

    def c_left(self):
        self.entry.icursor(max(0, self.entry.index("insert") - 1))

    def c_right(self):
        self.entry.icursor(self.entry.index("insert") + 1)

    def c_up(self):
        if self.expr_hist:
            self.hpos = min(self.hpos + 1, len(self.expr_hist))
            self.var.set(self.expr_hist[-self.hpos])
            self.entry.icursor("end")
            self.fresh = False

    def c_down(self):
        self.hpos = max(0, self.hpos - 1)
        self.var.set(self.expr_hist[-self.hpos] if self.hpos else "")
        self.fresh = False

    def c_copy(self):
        self.clipboard_clear()
        self.clipboard_append(self.cv.itemcget(self.t_res, "text"))

    def c_del(self):
        pos = self.entry.index("insert")
        if pos:
            m = DEL_TOKEN.search(self.entry.get()[:pos])
            self.entry.delete(pos - (len(m.group()) if m else 1), pos)

    def c_ins(self):
        self.overwrite = not self.overwrite

    def c_ac(self):
        self.var.set("")
        self.show("0")
        self.fresh = self.hyp = False
        self.pending = None

    def c_on(self):
        if not self.on:
            self.on = True
            for it in (self.t_stat, self.t_ind, self.t_sub, self.t_res):
                self.cv.itemconfigure(it, state="normal")
            self.cv.itemconfigure(self.ent_win, state="normal")
        self.c_ac()

    def c_off(self):
        self.on = False
        for it in (self.t_stat, self.t_ind, self.t_sub, self.t_res):
            self.cv.itemconfigure(it, state="hidden")
        self.cv.itemconfigure(self.ent_win, state="hidden")

    def c_mode(self):
        m = tk.Menu(self, tearoff=0)
        for n, (label, code) in enumerate((("COMP", "COMP"), ("CMPLX", "CMPLX"), ("SD", "SD"),
                                           ("REG", "REG"), ("BASE-N", "BASE"), ("EQN", "EQN"),
                                           ("MAT", "MAT"), ("VCT", "VCT")), 1):
            m.add_command(label=f"{n}  {label}", command=lambda c=code: self.set_mode(c))
        m.add_separator()
        m.add_command(label="SETUP…", command=self.c_setup)
        self.popup(m)

    def c_clrmenu(self):
        m = tk.Menu(self, tearoff=0)
        m.add_command(label="1  Scl  (clear statistics)", command=self.c_cl)
        m.add_command(label="2  Mode (reset modes & setup)", command=self.reset_modes)
        m.add_command(label="3  All  (reset everything)", command=self.reset_all)
        self.popup(m)

    def reset_modes(self):
        e = self.eng
        e.angle, e.dm, e.dd, e.frac, e.eng_on, e.polar = "DEG", "norm1", 2, None, False, False
        self.set_mode("COMP")

    def reset_all(self):
        self.eng = Engine()
        self.stat_rows.clear()
        self.hist.delete(0, "end")
        self.expr_hist.clear()
        self.set_mode("COMP")

    def set_mode(self, m):
        self.mode = m
        self.eng.complex_mode = m == "CMPLX"
        self.eng.polar = self.eng.polar and m == "CMPLX"
        self.c_ac()
        self.refresh()
        if m == "EQN":
            EqnDialog(self)
        elif m in ("MAT", "VCT"):
            MatDialog(self, m)
        elif m in ("SD", "REG"):
            self.show("Stat", sub="M+ = DT   SHIFT+1/2 = S-SUM / S-VAR")

    def cycle_angle(self):
        order = ["DEG", "RAD", "GRA"]
        self.eng.angle = order[(order.index(self.eng.angle) + 1) % 3]
        self.refresh()

    def c_setup(self):
        SetupDialog(self)

    def c_const(self):
        ConstDialog(self)

    def c_conv(self):
        ConvDialog(self)

    def c_mat(self):
        MatDialog(self, "MAT")

    def c_vct(self):
        MatDialog(self, "VCT")

    def c_ssum(self):
        StatDialog(self)

    c_svar = c_ssum

    def c_distr(self):
        DistrDialog(self)

    def c_calc(self):
        s = self.var.get().strip()
        if s:
            VarDialog(self, s, False)

    def c_solve(self):
        s = self.var.get().strip()
        if s:
            VarDialog(self, s, True)

    def c_logic(self):
        m = tk.Menu(self, tearoff=0)
        for label, tok in (("and", " and "), ("or", " or "), ("xor", " xor "), ("xnor", " xnor "),
                           ("Not", "not("), ("Neg", "neg(")):
            m.add_command(label=label, command=lambda t=tok: self.put(t))
        self.popup(m)

    def c_drgmenu(self):
        m = tk.Menu(self, tearoff=0)
        for label, f in (("° → ʳ", math.radians), ("° → ᵍ", lambda x: x * 10 / 9),
                         ("ʳ → °", math.degrees), ("ʳ → ᵍ", lambda x: x * 200 / math.pi),
                         ("ᵍ → °", lambda x: x * 0.9), ("ᵍ → ʳ", lambda x: x * math.pi / 200)):
            m.add_command(label=label, command=lambda f=f: self.drg(f))
        self.popup(m)

    def drg(self, f):
        try:
            self.eng.ans = f(float(self.entry_value()))
            self.show(self.eng.fmt(self.eng.ans))
        except Exception as ex:
            self.error(ex)

    def c_dmsch(self):
        pos = self.entry.index("insert")
        m = re.search(r"[\d.°′″]+$", self.entry.get()[:pos])
        seg = m.group() if m else ""
        self.put("°" if "°" not in seg else "′" if "′" not in seg else "″")

    def c_dms(self):
        try:
            self.show(to_dms(float(self.eng.ans)))
        except Exception as ex:
            self.error(ex)

    def c_dc(self):
        self.eng.frac = "improper" if self.eng.frac != "improper" else "mixed"
        self.redisplay()

    def c_eng(self):
        self.eng.eng_on = not self.eng.eng_on
        self.redisplay()

    def c_polar(self):
        self.eng.polar = True
        self.redisplay()

    def c_rect(self):
        self.eng.polar = False
        self.redisplay()

    def c_reim(self):
        z = complex(self.eng.ans)
        self.show(f"Re {self.eng.fmt_real(z.real)}", sub=f"Im {self.eng.fmt_real(z.imag)}")

    def _setbase(self, b):
        self.base = b
        self.redisplay()

    c_b10 = lambda self: self._setbase(10)
    c_b16 = lambda self: self._setbase(16)
    c_b2 = lambda self: self._setbase(2)
    c_b8 = lambda self: self._setbase(8)

    # memory ----------------------------------------------------------
    def c_sto(self):
        self.pending = "STO"

    def c_rcl(self):
        self.pending = "RCL"

    def entry_value(self):
        s = self.var.get().strip()
        return self.eng.evaluate(s) if s else self.eng.ans

    def store(self, name):
        try:
            self.eng.vars[name] = self.entry_value()
            self.show(self.eng.fmt(self.eng.vars[name]), sub=f"→ {name}")
        except Exception as ex:
            self.error(ex)

    def _mem(self, sign):
        try:
            self.eng.vars["M"] += sign * self.entry_value()
            self.show(self.eng.fmt(self.eng.vars["M"]), sub="M")
        except Exception as ex:
            self.error(ex)

    def c_mplus(self):
        self._mem(1)

    def c_mminus(self):
        self._mem(-1)

    def c_dt(self):
        s = self.var.get().strip()
        if not s:
            return
        try:
            body, _, fq = s.partition(";")
            freq = float(self.eng.evaluate(fq)) if fq.strip() else 1.0
            vals = [float(self.eng.evaluate(t)) for t in body.split(",")]
            if len(vals) != (1 if self.mode == "SD" else 2):
                raise SyntaxError("Syntax ERROR")
            self.stat_rows.append((*vals, freq))
            self.stat_type = "SD" if self.mode == "SD" else (self.stat_type if self.stat_type != "SD" else "Lin")
            self.var.set("")
            self.show(f"n = {self.eng.fmt_real(sum(r[-1] for r in self.stat_rows))}", sub="DT")
        except Exception as ex:
            self.error(ex)

    def c_cl(self):
        self.stat_rows.clear()
        self.show("Stat clear")

    # equals ----------------------------------------------------------
    def c_eq(self):
        expr = self.var.get().strip()
        if not expr:
            return
        try:
            if self.mode == "BASE":
                v = eval_base(expr, self.base)
                out, sub = base_str(v, self.base), self.base_line(v)
            else:
                v = self.eng.evaluate(expr)
                out, sub = self.eng.fmt(v), self.eng.pair
        except Exception as ex:
            return self.error(ex)
        if not isinstance(v, Matrix):
            self.eng.ans = v
        self.expr_hist.append(expr)
        self.hpos, self.fresh = 0, True
        self.hist.insert(0, f"{expr}  =  {out.splitlines()[0]}")
        self.show(out, sub=sub)


if __name__ == "__main__":
    import sys
    App(float(sys.argv[1]) if len(sys.argv) > 1 else None).mainloop()