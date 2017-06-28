# Jackie Tseng 2017.06.16
# Advanced Simplex Algorithm
# -*- coding:utf-8 -*-
from prettytable import PrettyTable
import numpy as np
import sys

class Simplex:
    def __init__(self):
        self.status = 0
        self.table = None
    
    def set_para(self, _A, _b, _c, _bases, _z = 0):
        self.A = _A
        self.b = _b
        self.c = _c
        self.z = _z
        self.bases = _bases
        self.m, self.n = self.A.shape
        self.c_b = np.zeros(len(_bases))
        self.C_B_inv = None

    def build_table(self):
        self.table = np.zeros((self.m + 1, self.n + 1))
        self.table[:-1, 1:]  = self.A
        self.table[:-1, :1]  = self.b.T
        self.table[ -1, 1:]  = self.c
        self.table[ -1,  0]  = self.z
        self.baseVars = self.bases

    def show(self):
        B_inv = self.table[0 : self.m, self.n - self.m + 1 : self.n + 1]
        self.C_B_inv = np.dot(self.c_b, B_inv)
        print ('------------Simplex Table------------')
        np.set_printoptions(precision=2)
        #print self.table
        table_viz = PrettyTable()
        table_viz.field_names = [' ', '  ','C_j', '   ']
        table_viz.add_row(['C_B', 'X_B', 'C_B_B^(-1)', 'P_0'])
        for i in range(self.m):
            if i != self.m:
                table_viz.add_row([str(self.c_b[i]), str('x') + str(self.baseVars[i] + 1), self.C_B_inv[i], self.table[i, 0]])
            else:
                table_viz.add_row([str(' '), str('x') + str(self.baseVars[i] + 1), self.table[i, 0]])
        table_viz.add_row([str(' '), str(' '), str('z'), -1 * self.table[self.m, 0]])
        index = 1
        while index < self.n + 1:
            table_viz.add_column(str(self.c[0, index - 1]), ([str('x') + str(index)] + [str(x) for x in self.table[:, index]]))
            index += 1
        print table_viz
        print ('------------Current Result-----------')
        print -1 * self.table[self.m, 0]
        print '\n'

    def is_optimal(self):
        for i in range(self.n):
            if i not in self.baseVars:
                if self.table[-1, i + 1] > 0:
                    return False
        self.status = 2
        return True

    def no_solution(self):
        for i in range(self.n):
            if i not in self.baseVars:
                if self.table[-1, i + 1] > 0:
                    no_solution_flag = True
                    for j in self.table[:-1, i + 1]:
                        if j > 0:
                            no_solution_flag = False
                    if no_solution_flag == True:
                        self.status = 1
                        return True
        return False

    def calculate(self):
        # Calculate swap-in variable
        max_value = 0
        inVar = None
        for i in range(self.n):
            if i not in self.baseVars:
                value = self.table[-1, i + 1]
                if value > max_value:
                    max_value = value
                    inVar = i

        # Calculate swap-out variable
        rates = []
        for nobaseVar in range(self.m):
            a = self.table[nobaseVar, inVar + 1]
            b = self.table[nobaseVar, 0]
            if a > 0:
                rate = b / a
                rates.append((rate, nobaseVar))
        outVar = min(rates)[1]

        # Calculate B^(-1)
        # Swap
        a = self.table[outVar, inVar + 1]
        self.table[outVar, :] /= a
        for i in range(self.m + 1):
            if i != outVar:
                self.table[i, :] -= self.table[outVar, :] * self.table[i, inVar + 1]
        self.baseVars[outVar] = inVar
        self.c_b[outVar] = self.c[0, inVar]

    def run(self):
        self.build_table()     
        while True:
            self.show()
            if self.is_optimal() or self.no_solution():
                break
            self.calculate()
        if self.status == 1:
            print ('Unbounded Solution')
        elif self.status == 2:
            print ('Optimal Solution = %.2f'%(-1 * self.table[self.m, 0]))
        else:
            print ('Infinite Optimal Solution')

def formulate(A, b, c, bases):
    total_item = A.shape[1]
    print(str('max z =')),
    index = 0
    while index < c.shape[1] - 1:
        print(str(c[0, index]) + 'x' + str(index + 1) + ' +'),
        index += 1
    print(str(c[0, index]) + 'x' + str(index + 1))
    for i in range(A.shape[0]):
        index = 0
        print(str(A[i, index]) + 'x' + str(index + 1)),
        index += 1
        flag = False
        while index < total_item - 1:
            if A[i, index] != 0:
                flag = True
                print('+ ' + str(A[i, index]) + 'x' + str(index + 1)),
            index += 1
        if A[i, index] != 0:
            if flag:
                print(str(A[i, index]) + 'x' + str(index + 1) + ' = ' + str(b[0, i]))
            else:
                print('+ ' + str(A[i, index]) + 'x' + str(index + 1) + ' = ' + str(b[0, i]))

        else:
            print('= ' + str(b[0, i]))
    print('xi >= 0, i ='),
    for i in range(total_item):
        if i != total_item - 1:
            print(str(i + 1) + ','),
        else:
            print(str(i + 1))
    print('\n')

if __name__ == "__main__":
    solution = Simplex()
    if sys.argv[1] == str(1):
        # Test case 1
        A = np.matrix([[1, 1, 2, 1, 0, 0],
                       [2, 0, 3, 0, 1, 0],
                       [2, 1, 3, 0, 0, 1]])
        b = np.matrix([[4, 5, 7]])
        c = np.matrix([[3, 2, 4, 0, 0, 0]])
        bases = [3, 4, 5]
    elif sys.argv[1] == str(2):
        # Test case 2
        A = np.matrix([[ 1, 1,  2, 1, 0, 0],
                       [ 1, 1, -1, 0, 1, 0],
                       [-1, 1,  1, 0, 0, 1]])
        b = np.matrix([[9, 2, 4]])
        c = np.matrix([[-1, -1, 4, 0, 0, 0]])
        bases = [3, 4, 5]
    else:
        # Test case 3
        A = np.matrix([[ 30, 20,  1, 0, 0],
                       [ 5,   1,  0, 1, 0],
                       [ 1,   0,  0, 0, 1]])
        b = np.matrix([[160, 15, 4]])
        c = np.matrix([[5, 2, 0, 0, 0]])
        bases = [2, 3, 4]
    formulate(A, b, c, bases)
    solution.set_para(A, b, c, bases)
    solution.run()
