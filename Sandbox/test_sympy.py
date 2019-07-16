from sympy import Symbol, Matrix, eye

A = Matrix([[0, 1, 0], [0, 0, 1], [-1, 0, 0]])
B = Matrix([[0], [0], [1]])
C = Matrix([[0, 1, 0]])

h = 0.1
G = eye(3) + h * A
H = h * B

print(G)
print(H)

zm1 = Symbol("z")

temp1 = eye(3) - zm1 * G
temp2 = temp1.inv()
temp3 = C * temp2 * zm1 * H
print(temp1)
print(temp2)
print(temp3)

