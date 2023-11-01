import timeit

# Given n+1 distinct nodes there exists a unique polynomial of degree n which interpolates each node
# test_set = [(1, 2), (2, 3), (4, 11), (6, 5)]


# for this set of nodes, there exists a polynomial such that:
# p(x0) = a0 + a1x0 + a2x0^2 + a3x0^3 = y0
# p(x1) = a0 + a1x1 + a2x1^2 + a3x1^3 = y0
# p(x2) = a0 + a1x2 + a2x2^2 + a3x2^3 = y0
# p(x3) = a0 + a1x3 + a2x3^2 + a3x3^3 = y0

# And so
# p(1) = a0 + a1 + a2 + a3 = 2
# p(1) = a0 + 2a1 + 4a2 + 8a3 = 3
# p(1) = a0 + 4a1 + 16a2 + 64a3 = 11
# p(1) = a0 + 6a1 + 36a2 + 216a3 = 5

# This is a *simulteneous equation* and so can be rewritten as:
# | 1  1  1  1   | | a0 |     | 2  |
# | 1  2  4  6   | | a1 | --- | 3  |
# | 1  4  16 64  | | a2 | --- | 11 |
# | 1  6  36 216 | | a3 |     | 5  |

# VA = Y where V is the vandermonde matrix of each value of x, A is the matrix of each a, and Y is the matrix of each y
# so A = V^-1 Y

# so just find the inverse of the matrix!

def get_cofactors(matrix):
    n = len(matrix)
    return [[(-1) ** ((n * i) + l) for l in range(n)] for i in range(n)]


def get_determinant(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    cofactored_row = [i * ((-1) ** index) for index, i in enumerate(matrix[0])]
    det = 0
    for index, cofactor in enumerate(cofactored_row):
        next_mat = []
        for row in matrix[1:]:
            new_row = row[:index] + row[index + 1:]
            next_mat.append(new_row)
        det += (cofactor * get_determinant(next_mat))
    return det


def get_cofactor_matrix(matrix):
    cofactors = get_cofactors(matrix)
    cofactor_mat = []
    for rowindex, row in enumerate(matrix):
        cofactor_row = []
        for itemindex, item in enumerate(row):
            temp_mat = matrix[:rowindex] + matrix[rowindex + 1:]
            minimat = []
            for temp_row in temp_mat:
                mini_row = temp_row[:itemindex] + temp_row[itemindex + 1:]
                minimat.append(mini_row)
            cofactor_row.append(get_determinant(minimat) * cofactors[rowindex][itemindex])
        cofactor_mat.append(cofactor_row)
    return cofactor_mat


def get_transposed_matrix(matrix):
    n = len(matrix)
    m = len(matrix[0])
    transposed_mat = [[0 for _ in range(n)] for _ in range(m)]
    for rowindex, row in enumerate(matrix):
        for itemindex, item in enumerate(row):
            transposed_mat[itemindex][rowindex] = item
    return transposed_mat


def get_inverse_matrix(matrix):
    determinant = get_determinant(matrix)
    adjugate = get_transposed_matrix(get_cofactor_matrix(matrix))
    return determinant, adjugate


def get_vandermonde_matrix(N, nodes):
    X = [i[0] for i in nodes]
    V = [[x ** n for n in range(N)] for x in X]
    return V


def get_y_matrix(nodes):
    Y = [[i[1]] for i in nodes]
    return Y


def matrix_multiply(N, M):
    tM = get_transposed_matrix(M)
    NM = []
    for Nrowindex, Nrow in enumerate(N):
        NMrow = []
        for tMrowindex, tMrow in enumerate(tM):
            NMrow.append(sum([Nrow[n] * tMrow[n] for n in range(len(Nrow))]))
        NM.append(NMrow)
    return NM


def constant_multiply(c, M):
    cM = [[c * item for item in row] for row in M]
    return cM


def interpolate(nodes):
    V = get_vandermonde_matrix(len(nodes), nodes)
    Y = get_y_matrix(test_set)
    DetV, AdjV = get_inverse_matrix(V)
    AdjVY = matrix_multiply(AdjV, Y)
    A = constant_multiply(1/DetV,AdjVY)
    return A


test_set = [(1, 2), (2, 3), (4, 11), (6, 5)]
A = interpolate(test_set)
out = "The polynomial that fits the points" + "".join([f" {i}," for i in test_set])
out = out[:-1] + " is y = " + "".join([f"{a[0]}x^{n} " for n,a in enumerate(A)])
print(out)
