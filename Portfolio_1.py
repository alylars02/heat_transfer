import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

def mega_function(guess, known):
    T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21 = guess
    To_inf, ho, Ti_inf, hi, k, dx = known
    
    # Row 1
    T1_eq = symmetric_plane_conv(T1, T7, T2, dx, k, ho, To_inf)
    T2_eq = plane_conv(T2, T8, T1, T3, dx, k, ho, To_inf)
    T3_eq = plane_conv(T3, T9, T2, T4, dx, k, ho, To_inf)
    T4_eq = plane_conv(T4, T10, T3, T5, dx, k, ho, To_inf)
    T5_eq = plane_conv(T5, T11, T4, T6, dx, k, ho, To_inf)
    T6_eq = symmetric_plane_conv(T6, T12, T5, dx, k, ho, To_inf)

    # Row 2
    T7_eq = symmetric_interior(T7, T8, T1, T13)
    T8_eq = interior(T8, T7, T9, T2, T14)
    T9_eq = interior(T9, T8, T10, T3, T15)
    T10_eq = interior(T10, T9, T11, T4, T16)
    T11_eq = interior(T11, T10, T12, T5, T17)
    T12_eq = symmetric_interior(T12, T11, T6, T18)

    # Row 3
    T13_eq = symmetric_interior(T13, T14, T7, T19)
    T14_eq = interior(T14, T13, T15, T8, T20)
    T15_eq = corner_conv(T15, T14, T9, T16, T21, dx, k, hi, Ti_inf)
    T16_eq = plane_conv(T16, T10, T15, T17, dx, k, hi, Ti_inf)
    T17_eq = plane_conv(T17, T11, T16, T18, dx, k, hi, Ti_inf)
    T18_eq = symmetric_plane_conv(T18, T12, T17, dx, k, hi, Ti_inf)

    #Row 4
    T19_eq = 2* T13 + 2* T20 - 4 * T19
    T20_eq = symmetric_interior(T20, T14, T19, T21)
    T21_eq = symmetric_plane_conv(T21, T20, T15, dx, k, hi, Ti_inf)

    return T1_eq, T2_eq, T3_eq, T4_eq, T5_eq, T6_eq, T7_eq, T8_eq, T9_eq, T10_eq, T11_eq, T12_eq, T13_eq, T14_eq, T15_eq, T16_eq, T17_eq, T18_eq, T19_eq, T20_eq, T21_eq

def interior(T_mn, T_L, T_R, T_U, T_D):
    # T_L: right node
    # T_R: right node
    # T_U: adjacent above node (up)
    # T_D: adjacent below node (down)
    return T_L + T_R + T_U + T_D - 4 * T_mn

def corner_conv(T_mn, T_in1, T_in2, T_air1, T_air2, dx, k, h, T_inf):
    # T_in: the internal adjacent nodes
    # T_air: the adjacent nodes open to air
    return 2*(T_in1 + T_in2) + (T_air1 + T_air2) + 2*h*dx/k*T_inf - 2*(3 + h*dx/k)*T_mn

def plane_conv(T_mn, T_in, T_air1, T_air2, dx, k, h, T_inf):
    # T_in: the internal adjacent nodes
    # T_air: the adjacent nodes open to air
    return (2*T_in + T_air1 + T_air2) + 2*h*dx/k*T_inf - 2*(h*dx/k + 2)*T_mn

def symmetric_interior(T_mn, T_S, T_pair1, T_pair2):
    # T_S: single adjacent node
    # T_pair: 1 of the paired adjacent nodes
    return T_S + T_S + T_pair1 + T_pair2 - 4 * T_mn

def symmetric_plane_conv(T_mn, T_in, T_air, dx, k, h, T_inf):
    # T_in: the internal adjacent nodes
    # T_air: the adjacent nodes open to air
    return (2*T_in + 2*T_air) + 2*h*dx/k*T_inf - 2*(h*dx/k + 2)*T_mn

def plot_matrix(matrix):
    # Create figure
    fig, ax = plt.subplots()

    # Display the matrix as an image with the red-orange colormap
    cax = ax.imshow(matrix, cmap='Reds', interpolation='nearest')

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)  # Remove border around plot

    # Remove color bar label
    cbar = plt.colorbar(cax)
    cbar.ax.set_ylabel("")

    # Annotate each cell with its value (smaller, non-bold text)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isnan(matrix[i, j]):  # Skip NaN values
                text = f"{matrix[i, j]:.1f}"  # Format to 1 decimal place
                ax.text(j, i, text, ha='center', va='center', color='grey', fontsize=10)

    # Set title
    plt.title("Temperture Field", fontsize=14)

    # Show plot
    plt.show()

#Knowns
To_inf = 1700
ho = 1000
Ti_inf = 400
hi = 200
k = 25
dx = .001
known = np.array([To_inf, ho, Ti_inf, hi, k, dx])

#Guesses
guess = np.array([2000]*21)


solution = root(mega_function, guess, args=(known))
T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21 = solution.x

matrix = np.array([[T1, T2, T3, T4, T5, T6], 
                   [T7, T8, T9, T10, T11, T12], 
                   [T13, T14, T15, T16, T17, T18], 
                   [T19, T20, T21, np.nan, np.nan, np.nan]])



plot_matrix(matrix)