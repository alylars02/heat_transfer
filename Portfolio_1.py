import numpy as np
import sympy as sp

one = 












def interior_node(T_mn, T_m1n, T_mn1, T_m_1n, T_mn_1):
    return T_m1n + T_mn1 + T_m_1n + T_mn_1 - 4 * T_mn

def internal_corner_convection(T_mn, T_m1n, T_mn1, h, dx, k, T_inf):
    return T_m1n + T_mn1 - 2 * T_mn + (h * dx / k) * (T_inf - T_mn)

def plane_surface_convection(T_mn, T_m1n, h, dx, k, T_inf):
    return T_m1n - T_mn + (h * dx / (2 * k)) * (T_inf - T_mn)

def external_corner_convection(T_mn, T_m1n, T_mn1, h, dx, k, T_inf):
    return T_m1n + T_mn1 - 2 * T_mn + (h * dx / k) * (T_inf - T_mn)

def plane_surface_heat_flux(T_mn, T_m1n, T_mn1, q_prime, dx, k):
    return T_m1n + T_mn1 - 2 * T_mn - (q_prime * dx / k)

# Example usage:
T_mn = 100  # Example temperature at node (m,n)
T_m1n, T_mn1, T_m_1n, T_mn_1 = 110, 105, 95, 90  # Neighboring node temperatures
h, dx, k, T_inf = 10, 0.01, 200, 300  # Convection parameters
q_prime = 500  # Heat flux

print("Interior Node Residual:", interior_node(T_mn, T_m1n, T_mn1, T_m_1n, T_mn_1))
print("Internal Corner with Convection Residual:", internal_corner_convection(T_mn, T_m1n, T_mn1, h, dx, k, T_inf))
print("Plane Surface with Convection Residual:", plane_surface_convection(T_mn, T_m1n, h, dx, k, T_inf))
print("External Corner with Convection Residual:", external_corner_convection(T_mn, T_m1n, T_mn1, h, dx, k, T_inf))
print("Plane Surface with Heat Flux Residual:", plane_surface_heat_flux(T_mn, T_m1n, T_mn1, q_prime, dx, k))
