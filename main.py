from math import *
import numpy as np
from tkinter import *

root = Tk()
root.geometry('1024x620')
canvas = Canvas(root, width=1024, height=580, bg='white')
canvas.pack(side="bottom")

gus_x1 = np.array([0.1, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5, 2.9, 3.3, 3.7], float)
gus_y1 = np.array([-0.6575, -0.0875, 1.4905, 5.4205, 13.0465, 25.7125, 44.7625, 71.5405, 107.3905, 153.6565], float)

gus_x2 = np.array([0.25, 0.3, 0.65, 0.7, 0.85, 1.1, 1.5, 1.55, 1.6, 1.75], float)
gus_y2 = np.array([-0.4766, -0.4145, 0.3139, 0.4915, 1.1922, 3.0775, 8.6875, 9.6683, 10.72, 14.3266], float)

gus_x3 = np.array([0.1, 0.35, 0.7, 0.8, 0.95, 1.2, 1.4, 1.9, 2, 2.3], float)
gus_y3 = np.array([-0.6575, -0.3472, 0.4915, 0.928, 1.8256, 4.144, 6.928, 18.6655, 22, 34.3555], float)

gus_brige_viett_coef = np.array([3, -5.6, 1.43, 1.207], float)

my_x = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0], float)
my_y = np.array([0.17, -0.744, -0.992, -0.4, 1.2, 3.976, 8.096, 13.728, 21.04, 30.2], float)

brige_viett_coef = np.array([4, -5, -2, 2.98], float)

kurs_koef = np.array([25, 7.5, 1], float)
kurs_diap_min = 0
kurs_diap_max = 2
y0 = 15
fx = 14


def make_graphics(x, y, flag, color='red'):
    '''flag = 0 from array
       flag = 1 from coordinates'''
    x0 = 50  # =265
    k = 20  # 10

    # задание 1
    canvas.create_line(512, 0, 512, 580, fill='black')
    canvas.create_line(0, x0, 1024, x0, fill='black')
    i = 0
    while i < 1024:
        canvas.create_text(490, (x0 + i), text=i / k)
        i += 10
    i = 0
    if flag == 0:
        while i < len(x) - 1:
            # canvas.create_oval(x[i]*10+512,y[i]*10+265,x[i]*10+513,y[i]*10+266,width=0.1, fill='red')
            canvas.create_line(x[i] * k + 512, y[i] * k + x0, x[i + 1] * k + 512, y[i + 1] * k + x0, fill=color)
            i += 1
    elif flag == 1:
        canvas.create_oval(x * k + 512, y * k + x0, x * k + 513, y * k + x0 + 1, width=1, fill=color)


def newton_polynome_back_func_delta_y(y):
    del1_y = np.zeros(len(y) - 1)
    i = 0
    while i < len(y) - 1:
        del1_y[i] = y[i + 1] - y[i]
        i += 1
    return del1_y


def newton_polynome_back(x, y, fx, num_x=9):
    delta1_y = newton_polynome_back_func_delta_y(y)
    print('delta1_y', delta1_y)
    delta2_y = newton_polynome_back_func_delta_y(delta1_y)
    print('delta2_y', delta2_y)
    delta3_y = newton_polynome_back_func_delta_y(delta2_y)
    print('delta3_y', delta3_y)
    h = x[1] - x[0]
    print("h = ", h)
    q = (fx - x[num_x]) / h
    print('q = ', q)
    result_newton_polynome_back = y[num_x] + ((delta1_y[num_x - 1] * q) / 1) + (
            delta2_y[num_x - 2] * q * (q + 1)) / factorial(2) + (
                                          delta3_y[num_x - 3] * q * (q + 1) * q + 2) / factorial(3)
    print('Интерпполяционный полином Ньютона, при х = 1,7, у =', result_newton_polynome_back)


def lagrange_polynome(x, y, fx, num_x):
    '''num_x - получает первое число меньше заданного'''

    x0 = fx - x[num_x - 1]
    y0 = y[num_x - 1]
    x1 = fx - x[num_x]
    y1 = y[num_x]
    x2 = fx - x[num_x + 1]
    y2 = y[num_x + 1]
    x3 = fx - x[num_x + 2]
    y3 = y[num_x + 2]

    print('x0=', x[num_x - 1], ' y0=', y0, ' x1=', x[num_x], ' y1=', y1, ' x2=', x[num_x + 1], ' y2=', y2, ' x3=',
          x[num_x + 2], ' y3=', y3)

    result_lagrange_polynome = ((fx - x1) * (fx - x2) * (fx - x3) * y0) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) \
                               + ((fx - x0) * (fx - x2) * (fx - x3) * y1) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) \
                               + ((fx - x0) * (fx - x1) * (fx - x3) * y2) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) \
                               + ((fx - x0) * (fx - x1) * (fx - x2) * y3) / ((x3 - x0) * (x3 - x1) * (x3 - x2))
    print('Интерполяционный полином Лагранжа, при х=1,7 у=', result_lagrange_polynome)


def koef_smallest_sqr_1(x, y, st, m):
    '''m - размер матрицы'''
    st+=1
    koef_matr = np.zeros((st, st))
    #koef_matr[0, 0] = m
    sol_matr = np.zeros(st)
    temp_val = np.zeros(2 * st)

    i = 0
    while i < (2 * st):  # цикл присваиваетт коэффициенты
        j = 0
        while j < m:  # цикл суммирует значения
            temp_val[i] += x[j] ** i
            if i<st:
                sol_matr[i] += y[j] * (x[j] ** i)
            j += 1

        i += 1
    temp_val[0] = m
    i = 0
    while i<st:
        j = 0
        while j<st:
            koef_matr[i,j]=temp_val[j+i]
            j+=1
        i+=1
    #print(koef_matr)
    #print(sol_matr)
    sol_koef = np.linalg.solve(koef_matr, sol_matr)
    #print(sol_koef)
    return sol_koef


def smallest_sqr_1(x, y):
    '''массив х, массив у, кол-во эл
    коэффициенты расположены от минимальной степени х к максимальной'''
    choose_st = []
    for st in range(3,len(x)):
        sol_koef=koef_smallest_sqr_1(x, y, st, len(x))
        i=0
        temp_sol=0

        while i<st:
            temp_sol += sol_koef[i]*(x[i]**i)

            i+=1
        choose_st.append(abs(temp_sol-y[st]))
        #print(st, sol_koef)
    j=1

    st_min=choose_st[0]
    while j<len(choose_st):
        if choose_st[j]<st_min:
            st_min=choose_st[j]
            j_min=j+2
        j+=1

    sol_koef = koef_smallest_sqr_1(x, y,j_min, len(x))

    print('степень многочлена: ', j_min, '\n коэфициенты многочлена (расположены от минимальной степени х к максимальной):\n', sol_koef)
    return sol_koef, j_min




def brige_viett(coef):
    i = -10000
    while i < 10000:
        make_graphics(i, coef[0] * (i ** 3) + coef[1] * (i ** 2) + coef[2] * i + coef[3], 1)
        i += 1


def metod_Eyler(kurs_koef, kurs_diap_min, kurs_diap_max, y0, fx, h):
    x = kurs_diap_min
    i = 1
    kurs_y = []
    kurs_y.append(y0)
    kurs_x = []
    kurs_x.append(kurs_diap_min)

    while x <= kurs_diap_max:
        tmp_y = kurs_y[i - 1] + h * ((25 * cos(7.5 * x)) / (0.5 * x + 2))
        kurs_y.append(tmp_y)
        kurs_x.append(x)
        print('!',i,'\t', x,'\t', x+h,'\t', tmp_y,'\t', h,)
        i += 1
        x += h
    # print(kurs_y)
    # print(kurs_x)
    make_graphics(kurs_x, kurs_y, 0, 'red')
    # eiler_coef=smallest_sqr(kurs_x, kurs_y)

    eiler_coef, st = smallest_sqr_1(kurs_x, kurs_y)
    # print(eiler_coef)
    # return eiler_coef

    x = kurs_diap_min

    while x <= kurs_diap_max:
        i=0
        y=0
        while i<st:
            y += eiler_coef[i]*(x**i)
            i+=1

        #gr_y = eiler_coef[3] * (x ** 3) + eiler_coef[2] * (x ** 2) + eiler_coef[1] * x + eiler_coef[0]
        # gr_y=eiler_coef[4]*(x**4)+eiler_coef[3]*(x**3)+eiler_coef[2]*(x**2)+eiler_coef[1]*x+eiler_coef[0]
        canvas.create_oval(x * 20 + 512, y * 20 + 50, x * 20 + 513, y * 20 + 50 + 1, width=1, fill='green')

        x += h
    return eiler_coef


def proizvodnaya(arr):
    i = 0
    pr_arr = []
    while i < len(arr) - 1:
        pr_arr.append(arr[i] * (len(arr) - 1 - i))
        i += 1
    #print("производная ", pr_arr)

    return pr_arr


def give_x(coef, x):
    i = 0
    result = 0
    # print('coef', coef, '\nLENGTH',len(coef), 'X=', x)
    while i < len(coef):
        result += coef[i] * (x ** (len(coef) - 1 - i))
        #print(x, '=x : промежуточный', result)
        i += 1
    #print(len(coef), '=l : итоговый', result)
    return result


def find_root(metod_Eyler_coef, kurs_diap_min, kurs_diap_max, y0, fx, h):
    metod_Eyler_proizv = proizvodnaya(metod_Eyler_coef)
    i = kurs_diap_min
    xk = 0 - (give_x(metod_Eyler_coef, 0) / give_x(metod_Eyler_proizv, 0))
    xk1 = 5
    print(metod_Eyler_coef)
    print(metod_Eyler_proizv)
    while abs(xk - xk1) > 0.0001:
        xk1 = xk
        xk = xk - (give_x(metod_Eyler_coef, xk) / give_x(metod_Eyler_proizv, xk))

        ##i+=h
    print('\nрешение:', xk1)


# newton_polynome_back(my_x, my_y, 1.7, 8)
# lagrange_polynome(my_x, my_y, 1.7, 7)
#smallest_sqr_1(gus_x3, gus_y3)

# brige_viett(gus_brige_viett_coef)
metod_Eyler_coef = metod_Eyler(kurs_koef, kurs_diap_min, kurs_diap_max, y0, fx, 0.01)
find_root(metod_Eyler_coef, kurs_diap_min, kurs_diap_max, y0, fx, 0.1)


# make_graphics(my_x,my_y,0)


root.mainloop()
