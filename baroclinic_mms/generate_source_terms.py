"""
Generates setups for 2d shallow water MMS tests

"""
import sympy
import numbers
from sympy import init_printing
init_printing()

# coordinates
x, y, z = sympy.symbols('xyz[0] xyz[1] xyz[2]')
x_2d, y_2d = sympy.symbols('xy[0] xy[1]')
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
temp0 = sympy.symbols('temp_const', positive=True)
# gravitational acceleration
g = sympy.symbols('g_grav')
# time
t = sympy.symbols('t', positive=True)
T = sympy.symbols('T', positive=True)

eos_alpha, eos_beta, eos_t0, eos_s0 = sympy.symbols('eos_alpha eos_beta eos_t0 eos_s0')
rho0 = sympy.symbols('rho_0', positive=True)

# setup 1 -- temp only
#bath = depth
#elev = 0
#u = 0
#v = 0
#temp = 5*sympy.cos((2*x + y)/lx)*sympy.cos((z/depth)) + 15
#salt = salt0
#nu = nu0
#f = f0
#omit_int_pg = False

# setup 2 -- elevation only
#bath = depth
#elev = 5*sympy.cos((x + 5*y/2)/lx)
#u = 0
#v = 0
#temp = temp0
#salt = salt0
#nu = nu0
#f = f0
##omit_int_pg = True
#omit_int_pg = False  # setup2b

# setup 3 -- velocity only
#bath = depth
#elev = 0
#u = sympy.cos((2*x + y)/lx)
#v = 0
#temp = temp0
#salt = salt0
#nu = nu0
#f = f0
#omit_int_pg = True

# setup 4 -- velocity and temp
#bath = depth
#elev = 0
#u = sympy.cos((2*x + y)/lx)*sympy.cos(3*(z/depth))/2
#v = 0
#temp = 5*sympy.cos((x + 2*y)/lx)*sympy.cos((z/depth)) + 15
#salt = salt0
#nu = nu0
#f = f0
##omit_int_pg = True
#omit_int_pg = False  # setup4b

# setup 5 -- velocity and temp, symmetric
bath = depth
elev = 0
u = sympy.sin(2*sympy.pi*x/lx)*sympy.cos(3*(z/depth))/2
v = 0 # sympy.sin(sympy.pi*y/ly)*sympy.sin((z/depth))/3
temp = sympy.sin(sympy.pi*x/lx)*sympy.sin(sympy.pi*y/ly) + 15
salt = salt0
nu = nu0
f = f0
omit_int_pg = False


def split_velocity(u, v, eta, bath):
    u_2d = sympy.integrate(u, (z, -bath, eta))/(eta+bath)
    v_2d = sympy.integrate(v, (z, -bath, eta))/(eta+bath)
    u_3d = u - u_2d
    v_3d = v - v_2d
    return u_3d, v_3d, u_2d, v_2d

def evaluate_w(eta, u, v, bath):
    assert sympy.diff(bath, x) == 0 and sympy.diff(bath, y) == 0, 'bath source not implemented'
    return sympy.integrate(-(sympy.diff(u, x) + sympy.diff(v, y)), (z, -bath, z))


def evaluate_baroclinicity(elev, temp, salt):
    rho = - eos_alpha*(temp - eos_t0) + eos_beta*(salt - eos_s0)
    baroc_head = sympy.integrate(rho/rho0, (z, elev, z))
    return rho, baroc_head


def evaluate_tracer_source(eta, trac, u, v, w, bath, f, nu):
    return sympy.diff(trac, x)*u + sympy.diff(trac, y)*v + sympy.diff(trac, z)*w


def evaluate_mom_source(eta, baroc_head, u, v, w, bath, f, nu):
    int_pg_x = -g*sympy.diff(baroc_head, x)  # NOTE why the minus sign? BUG
    int_pg_y = -g*sympy.diff(baroc_head, y)
    adv_u = sympy.diff(u, x)*u + sympy.diff(u, y)*v + sympy.diff(u, z)*w
    adv_v = sympy.diff(v, x)*u + sympy.diff(v, y)*v + sympy.diff(v, z)*w
    res_u = adv_u
    res_v = adv_v
    if not omit_int_pg:
        res_u += int_pg_x
        res_v += int_pg_y
    return res_u, res_v, int_pg_x, int_pg_y


def evaluate_swe_source(eta, u, v, bath, f, nu, nonlin=True):
    # evaluate shallow water equations
    if nonlin:
        tot_depth = eta + bath
    else:
        tot_depth = bath
    div_hu = sympy.diff(tot_depth*u, x) + sympy.diff(tot_depth*v, y)
    res_elev = sympy.diff(eta, t) + div_hu
    u_x = sympy.diff(u, x)
    u_y = sympy.diff(u, y)
    v_x = sympy.diff(v, x)
    v_y = sympy.diff(v, y)
    if nonlin:
        adv_u = u*u_x + v*u_y
        adv_v = u*v_x + v*v_y
    else:
        adv_u = adv_v = 0
    cori_u = -f*v
    cori_v = f*u
    pg_u = g*sympy.diff(eta, x)
    pg_v = g*sympy.diff(eta, y)
    visc_u = -(2*sympy.diff(nu*sympy.diff(u, x), x) +
               sympy.diff(nu*sympy.diff(u, y), y) +
               sympy.diff(nu*sympy.diff(v, x), y))
    visc_v = -(2*sympy.diff(nu*sympy.diff(v, y), y) +
               sympy.diff(nu*sympy.diff(v, x), x) +
               sympy.diff(nu*sympy.diff(u, y), x))
    visc_u += -sympy.diff(tot_depth, x)/tot_depth * nu * 2 * sympy.diff(u, x)
    visc_v += -sympy.diff(tot_depth, y)/tot_depth * nu * 2 * sympy.diff(v, y)
    # NOTE in the coupled system 2d mom eq has no advection/diffusion terms
    res_u = sympy.diff(u, t) + cori_u + pg_u
    res_v = sympy.diff(v, t) + cori_v + pg_v
    return res_elev, res_u, res_v



u_3d, v_3d, u_2d, v_2d = split_velocity(u, v, elev, bath)
w = evaluate_w(elev, u, v, bath)
rho, baroc_head = evaluate_baroclinicity(elev, temp, salt)

mom_source_x, mom_source_y, int_pg_x, int_pg_y = evaluate_mom_source(elev, baroc_head, u, v, w, bath, f, nu)

vol_source_2d, mom_source_2d_x, mom_source_2d_y = evaluate_swe_source(elev, u_2d, v_2d, bath, f, nu, nonlin=True)
temp_source_3d = evaluate_tracer_source(elev, temp, u, v, w, bath, f, nu)


def expr2str(e):
    if isinstance(e, (numbers.Number, sympy.numbers.Zero)):
        return 'Constant({:})'.format(e)
    return str(e)


def print_expr(name, *expr):
    if len(expr) == 1:
        expr_str = expr2str(expr[0])
    else:
        if all(isinstance(e, numbers.Number) for e in expr):
            comp_str = ', '.join([str(e) for e in expr])
            expr_str = 'Constant((' + comp_str + '))'
        else:
            comp_str = ', '.join([expr2str(e) for e in expr])
            expr_str = 'as_vector((' + comp_str + '))'
    print("    out['{:}'] = {:}".format(name, expr_str))

def to_2d_coords(expr):
    if isinstance(expr, numbers.Number):
        return expr
    return expr.subs(x, x_2d).subs(y, y_2d)

print_expr('elev_2d', to_2d_coords(elev))
print_expr('uv_full_3d', u, v, 0)
print_expr('uv_2d', to_2d_coords(u_2d), to_2d_coords(v_2d))
print_expr('uv_dav_3d', u_2d, v_2d, 0)
print_expr('uv_3d', u_3d, v_3d, 0)
print_expr('w_3d', 0, 0, w)
print_expr('temp_3d', temp)
print_expr('density_3d', rho)
print_expr('baroc_head_3d', baroc_head)
print_expr('int_pg_3d', int_pg_x, int_pg_y, 0)

print_expr('vol_source_2d', to_2d_coords(vol_source_2d))
print_expr('mom_source_2d', to_2d_coords(mom_source_2d_x), to_2d_coords(mom_source_2d_y))
print_expr('mom_source_3d', mom_source_x, mom_source_y, 0)
print_expr('temp_source_3d', temp_source_3d)
