"""
Generates setups for 2d shallow water MMS tests

"""
import sympy
from sympy import init_printing

init_printing()

# coordinates
x, y, z = sympy.symbols('xyz[0] xyz[1] xyz[2]')
z_tilde = sympy.symbols('z_tilde')
# domain lenght, x in [0, Lx], y in [0, Ly]
lx, ly = sympy.symbols('lx ly')
# depth scale
depth = sympy.symbols('depth', positive=True)
# coriolis scale
f0 = sympy.symbols('f0', positive=True)
# viscosity scale
nu0 = sympy.symbols('nu0', positive=True)
# constant salinity
salt0 = sympy.symbols('salt_const', positive=True)
# gravitational acceleration
g = sympy.symbols('g_grav')
# time
t = sympy.symbols('t', positive=True)
T = sympy.symbols('T', positive=True)

eos_alpha, eos_beta, eos_t0, eos_s0 = sympy.symbols('eos_alpha eos_beta eos_t0 eos_s0')
rho0 = sympy.symbols('rho_0', positive=True)
bath = depth
elev = 0

#temp = sympy.sin(0.2*sympy.pi*(3.0*x + 1.0*y)/lx)
temp = 5*sympy.cos((2*x + y)/lx)*sympy.cos((z/depth)) + 15

salt = salt0

# compute density and internal pressure gradient
rho = - eos_alpha*(temp - eos_t0) + eos_beta*(salt - eos_s0)
baroc_head = sympy.integrate(rho/rho0, (z, elev, z))
int_pg_x = -g*sympy.diff(baroc_head, x)  # NOTE why the minus sign? BUG
int_pg_y = -g*sympy.diff(baroc_head, y)

mom_source_x = int_pg_x
mom_source_y = int_pg_y


def print_expr(name, *expr):
    if len(expr) == 1:
        expr_str = str(expr[0])
    else:
        expr_str = 'as_vector((' + ', '.join([str(e) for e in expr]) + '))'
    print("out['{:}'] = {:}".format(name, expr_str))


print_expr('temperature_3d', temp)
print_expr('density_3d', rho)
print_expr('baroc_head_3d', baroc_head)
print_expr('int_pg_3d', int_pg_x, int_pg_y, 0)
print_expr('mom_source_3d', mom_source_x, mom_source_y, 0)
