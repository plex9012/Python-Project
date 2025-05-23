import tkinter



solutions = []
a1 = int(input("Input first parameter of quaternion algebra"))
b1 = int(input("Input second parameter of quaternion algebra"))
max_abs = int(input("Input size of search intervals"))

Q = QuaternionAlgebra(a1, b1)
order = Q.maximal_order().basis_matrix()


def Nrd(y):
    '''evaluates reduced norm of y'''
    return round(float(y[0][0] ^ 2 - a1 * y[1][0] ^ 2 - b1 * y[2][0] ^ 2 + a1 * b1 * y[3][0] ^ 2))


def phi(y):
    '''maps y from quaternion algebra to z from 2x2 matrices'''
    z = (
        y[0][0] + y[1][0] * sqrt(a1), y[2][0] + y[3][0] * sqrt(a1),
        b1 * y[2][0] - b1 * sqrt(a1) * y[3][0],
        y[0][0] - y[1][0] * sqrt(a1))
    return z


def psi(z):
    '''maps z from morphisms of upper-half plane to r from morphisms of unit circle'''
    m = matrix([[z[0], z[1]], [z[2], z[3]]])
    t = matrix([[1, 1j], [1j, 1]])
    r = t * m * t ^ (-1)
    return r[0][0], r[1][0], r[1][0], r[1][1]


def f(z, a, b, c, d):
    '''applies fractional-linear transformation to point z'''
    return (a * z + b) / (c * z + d)


def new_list(sol1, mx, my):
    '''applies fractional-linear transformation to fundamental domain'''
    n = 100
    v = float(mx) + 1j * float(my)
    v /= 10000
    check = 0
    while check == 0:
        try:
            for t in range(len(sol1)):
                z = sol1[t][2] + 1j * sol1[t][3]
                r = sol1[t][1]
                p = sol1[t][4]
                z1 = f(z + r, 1, v, v.conjugate(), 1)
                z2 = f(z - r, 1, v, v.conjugate(), 1)
                z3 = f(z + 1j * r, 1, v, v.conjugate(), 1)
                t1 = (1j * (z1 - z3) * ((z1 - z2).conjugate())).imag() / ((z3 - z2) * ((z1 - z2).conjugate())).imag()
                w = (z2 + z3) / 2 + 1j * t1 * (z2 - z3) / 2
                r1 = abs(w - z1)
                if abs(f(z, 1, v, v.conjugate(), 1) - w) > r1:
                    p = 1 - p
                sol1[t] = [sol1[0], r1, float(w.real()), float(w.imag()), p]
            check = 1
        except Exception:
            v += 1 / 10000
    return sol1


def click(event):
    '''saves coordinates of mouse'''
    global mouse
    mouse = [event.x, event.y]


def press(event):
    '''draws transformed fundamental domain'''
    global sol1
    global mouse
    sol1 = new_list(sol1, event.x - mouse[0], mouse[1] - event.y)
    list = []
    k = 4
    n = len(sol1)
    for i in range(n):
        list.append([k * sol1[i][2], k * sol1[i][3], k * sol1[i][1], sol1[i][4]])
    canvas.create_oval(2, 2, 722, 722, fill='white', outline='black')
    for q in list:
        if q[3] == 0:
            canvas.create_oval(361 + q[0] * 90 - q[2] * 90, 361 - q[1] * 90 - q[2] * 90, 361 + q[0] * 90 + q[2] * 90,
                               361 - q[1] * 90 + q[2] * 90, fill='white', outline='black')
        else:
            canvas.create_oval(361 + q[0] * 90 - q[2] * 90, 361 - q[1] * 90 - q[2] * 90, 361 + q[0] * 90 + q[2] * 90,
                               361 - q[1] * 90 + q[2] * 90, outline='black')
    for q in list:
        if q[3] == 0:
            canvas.create_oval(362 + q[0] * 90 - q[2] * 90, 362 - q[1] * 90 - q[2] * 90, 360 + q[0] * 90 + q[2] * 90,
                               360 - q[1] * 90 + q[2] * 90, fill='white', outline='white')
        else:
            canvas.create_oval(361 + q[0] * 90 - q[2] * 90 - 250, 361 - q[1] * 90 - q[2] * 90 - 250,
                               361 + q[0] * 90 + q[2] * 90 + 250,
                               361 - q[1] * 90 + q[2] * 90 + 250, outline='white', width=499)
    canvas.create_oval(2, 2, 722, 722, outline='black')
    canvas.create_oval(-148, -148, 872, 872, outline='white', width=299)


for x0 in range(-max_abs, max_abs + 1):
    for x1 in range(-max_abs, max_abs + 1):
        for x2 in range(-max_abs, max_abs + 1):
            for x3 in range(-max_abs, max_abs + 1):
                x = matrix([x0, x1, x2, x3])
                y = (x * order).transpose()
                if Nrd(y) == 1:
                    a, b, c, d = psi(phi(y))
                    if c != 0:
                        solutions.append([[a, b, c, d], 1 / abs(c), float((-d / c).real()), float((-d / c).imag()), 0])
s = set()
uniq = []
for i in solutions:
    rad = round(i[1], 3)
    real = round(i[2], 3)
    imag = round(i[3], 3)
    r = (rad, real, imag)
    if r not in s:
        s.add(r)
        uniq.append([i[0], rad, real, imag, 0])
sol = sorted(uniq, key=lambda x: -x[1])

for i in range(len(sol)):
    for j in range(len(sol)):
        if sol[j][1] + sqrt((sol[i][2] - sol[j][2]) ^ 2 + (sol[i][3] - sol[j][3]) ^ 2) - sol[i][1] + 0.05 < 0:
            sol[j][4] = 1
k = 4
n = min(100, len(sol))
list = []
sol1 = []
for i in range(n):
    if sol[i][4] == 0:
        sol1.append(sol[i])
n = len(sol1)
for i in range(n):
    list.append([k * sol1[i][2], k * sol1[i][3], k * sol1[i][1], sol1[i][4]])

main = tkinter.Tk()
main.title("Image")
mouse = []
canvas = tkinter.Canvas(main, width=724, height=724, bg="white")
canvas.pack()
label = tkinter.Label(main, text="Fundamental domain")
label.pack()
for q in list:
    if q[3] == 0:
        canvas.create_oval(361 + q[0] * 90 - q[2] * 90, 361 - q[1] * 90 - q[2] * 90,
                           361 + q[0] * 90 + q[2] * 90,
                           361 - q[1] * 90 + q[2] * 90, fill='white', outline='black')
    else:
        canvas.create_oval(361 + q[0] * 90 - q[2] * 90, 361 - q[1] * 90 - q[2] * 90,
                           361 + q[0] * 90 + q[2] * 90,
                           361 - q[1] * 90 + q[2] * 90, outline='black')
for q in list:
    if q[3] == 0:
        canvas.create_oval(362 + q[0] * 90 - q[2] * 90, 362 - q[1] * 90 - q[2] * 90,
                           360 + q[0] * 90 + q[2] * 90,
                           360 - q[1] * 90 + q[2] * 90, fill='white', outline='white')
    else:
        canvas.create_oval(361 + q[0] * 90 - q[2] * 90 - 250, 361 - q[1] * 90 - q[2] * 90 - 250,
                           361 + q[0] * 90 + q[2] * 90 + 250,
                           361 - q[1] * 90 + q[2] * 90 + 250, outline='white', width=499)
canvas.create_oval(2, 2, 722, 722, outline='black')
canvas.create_oval(-148, -148, 872, 872, outline='white', width=299)
canvas.bind('<Button-1>', click)
canvas.bind('<B1-Motion>', press)
main.mainloop()
