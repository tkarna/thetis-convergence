"""
MMS test for baroclinic 3D model
"""
from thetis import *
from scipy import stats

# TODO step 1: run model with non-trivial density field (lin eos) and add correct forcing term in momentum eq.


def setup1(xy, xyz, lx, ly, depth, salt_const, nu0, rho_0, g_grav, eos_params):
    """
    Constant bathymetry and zero velocty, non-trivial tracer
    """
    eos_alpha = eos_params['alpha']
    eos_beta = eos_params['beta']
    eos_t0 = eos_params['th_ref']
    eos_s0 = eos_params['s_ref']

    out = {}
    out['elevation_2d'] = Constant(0.0)
    out['uv_3d'] = as_vector((0, 0, 0))
    out['w_3d'] = as_vector((0, 0, 0))
    out['viscosity_3d'] = Constant(nu0)
    out['temperature_3d'] = 5*cos(xyz[2]/depth)*cos((2*xyz[0] + xyz[1])/lx) + 15
    out['density_3d'] = -eos_alpha*(-eos_t0 + 5*cos(xyz[2]/depth)*cos((2*xyz[0] + xyz[1])/lx) + 15) + eos_beta*(-eos_s0 + salt_const)
    out['baroc_head_3d'] = (-eos_alpha*(5*depth*sin(xyz[2]/depth)*cos((2*xyz[0] + xyz[1])/lx) - eos_t0*xyz[2] + 15*xyz[2]) + eos_beta*xyz[2]*(-eos_s0 + salt_const))/rho_0
    out['int_pg_3d'] = as_vector((-10*depth*eos_alpha*g_grav*sin(xyz[2]/depth)*sin((2*xyz[0] + xyz[1])/lx)/(lx*rho_0), -5*depth*eos_alpha*g_grav*sin(xyz[2]/depth)*sin((2*xyz[0] + xyz[1])/lx)/(lx*rho_0), 0))
    out['mom_source_3d'] = as_vector((-10*depth*eos_alpha*g_grav*sin(xyz[2]/depth)*sin((2*xyz[0] + xyz[1])/lx)/(lx*rho_0), -5*depth*eos_alpha*g_grav*sin(xyz[2]/depth)*sin((2*xyz[0] + xyz[1])/lx)/(lx*rho_0), 0))

    out['options'] = {}
    return out


def get_exact_sol(solver_obj, temp_expr, salt_expr):
    eos_params = solver_obj.options.equation_of_state_options
    rho = rho0 - eos_params['alpha']*(temp_expr - eos_params['th_ref']) + eos_params['beta']*(salt_expr - eos_params['s_ref'])
    return rho


def run(setup, refinement, polynomial_degree, do_export=True, **options):
    """Run single test and return L2 error"""
    print_output('--- running {:} refinement {:}'.format(setup.__name__, refinement))

    rho_0 = 1000.0
    physical_constants['rho0'] = rho_0
    g_grav = physical_constants['g_grav']

    # domain dimensions
    lx = 15e3
    ly = 10e3
    area = lx*ly
    depth = 40.0
    nu0 = 50.0
    t_end = 200.0

    bath_expr = Constant(depth)

    # mesh
    n_layers = 4*refinement
    nx = 4*refinement
    ny = 4*refinement
    mesh2d = RectangleMesh(nx, ny, lx, ly)

    alpha = 0.2  # thermal expansion coeff
    beta = 0.0  # haline contraction coeff
    temp_ref = 15.0
    salt_const = 10.0

    # outputs
    outputdir = 'outputs'

    # bathymetry
    p1_2d = FunctionSpace(mesh2d, 'CG', 1)
    bathymetry_2d = Function(p1_2d, name='Bathymetry')
    bathymetry_2d.project(bath_expr)

    solver_obj = solver.FlowSolver(mesh2d, bathymetry_2d, n_layers)
    options = solver_obj.options
    options.element_family = 'dg-dg'
    options.polynomial_degree = polynomial_degree
    options.use_bottom_friction = False
    options.use_baroclinic_formulation = True
    options.solve_temperature = True
    options.solve_salinity = False
    options.use_quadratic_pressure = True
    options.constant_salinity = Constant(salt_const)
    options.horizontal_velocity_scale = Constant(1.0)
    options.no_exports = not do_export
    options.output_directory = outputdir
    options.simulation_end_time = t_end
    options.fields_to_export = ['elev_2d', 'salt_3d', 'uv_3d', 'w_3d',
                                'temp_3d',
                                'density_3d', 'baroc_head_3d', 'int_pg_3d']
    options.horizontal_viscosity_scale = Constant(nu0)
    #salt = Constant(salt_const)
    eos_params = {
        'rho_ref': rho_0,
        's_ref': salt_const,
        'th_ref': temp_ref,
        'alpha': alpha,
        'beta': beta,
    }
    options.equation_of_state_type = 'linear'
    options.equation_of_state_options.update(eos_params)
    # diffusivuty
    #nu = Function(solver_obj.function_spaces.P1, name='diffusivity')
    #nu.project(sdict['nu_expr'])
    #options.horizontal_diffusivity = nu
    xyz = SpatialCoordinate(solver_obj.mesh)
    xy = SpatialCoordinate(solver_obj.mesh2d)
    sdict = setup(xy, xyz, lx, ly, depth, salt_const, nu0, rho_0, g_grav, eos_params)

    options.update(sdict['options'])
    options.update(options)

    solver_obj.create_function_spaces()

    # initial conditions
    temp_expr = sdict['temperature_3d']
    init_temp_3d = Function(solver_obj.function_spaces.H, name='initial temperature')
    init_temp_3d.project(temp_expr)

    # analytical solution
    eos_alpha = eos_params['alpha']
    eos_beta = eos_params['beta']
    eos_t0 = eos_params['th_ref']
    eos_s0 = eos_params['s_ref']

    ana_temp_3d = sdict['temperature_3d']
    f = Function(solver_obj.function_spaces.H, name='ana temp').project(ana_temp_3d)
    out = File(outputdir + '/ana_temp.pvd')
    out.write(f)
    out.write(f)

    ana_rho_3d = sdict['density_3d']
    f = Function(solver_obj.function_spaces.H, name='ana rho').project(ana_rho_3d)
    out = File(outputdir + '/ana_rho.pvd')
    out.write(f)
    out.write(f)

    ana_bhead_3d = sdict['baroc_head_3d']
    f = Function(solver_obj.function_spaces.H, name='ana bhead').project(ana_bhead_3d)
    out = File(outputdir + '/ana_bhead.pvd')
    out.write(f)
    out.write(f)

    ana_intpg_3d = sdict['int_pg_3d']
    f = Function(solver_obj.function_spaces.U, name='ana int pg').project(ana_intpg_3d)
    out = File(outputdir + '/ana_intpg.pvd')
    out.write(f)
    out.write(f)

    ana_uv_3d = sdict['uv_3d']

    options.momentum_source_3d = sdict['mom_source_3d']

    ## analytical solution in high-order space for computing L2 norms
    #h_ho = FunctionSpace(solver_obj.mesh, 'DG', order+3)
    #trac_ana_ho = Function(h_ho, name='Analytical T')
    #trac_ana_ho.project(sdict['tracer_expr'])
    ## analytical solution
    #trac_ana = Function(solver_obj.function_spaces.H, name='Analytical T')
    #trac_ana.project(sdict['tracer_expr'])

    bnd_temp = {'value': init_temp_3d}
    solver_obj.bnd_functions['temp'] = {1: bnd_temp, 2: bnd_temp,
                                        3: bnd_temp, 4: bnd_temp}
    # NOTE use symmetic uv condition to get correct w
    bnd_mom = {'symm': None}
    solver_obj.bnd_functions['momentum'] = {1: bnd_mom, 2: bnd_mom,
                                            3: bnd_mom, 4: bnd_mom}

    solver_obj.create_equations()
    solver_obj.assign_initial_conditions(temp=init_temp_3d)

    solver_obj.iterate()

    area = lx*ly
    l2_err = errornorm(ana_uv_3d, solver_obj.fields.uv_3d)/np.sqrt(area)
    print_output('L2 error {:.12f}'.format(l2_err))

    return l2_err


def run_convergence(setup, ref_list, saveplot=False, **options):
    """Runs test for a list of refinements and computes error convergence rate"""
    polynomial_degree = options.get('polynomial_degree', 1)
    space_str = options.get('element_family')
    l2_err = []
    for r in ref_list:
        l2_err.append(run(setup, r, **options))
    x_log = np.log10(np.array(ref_list, dtype=float)**-1)
    y_log_uv = np.log10(np.array(l2_err))
    setup_name = setup.__name__

    def check_convergence(x_log, y_log, expected_slope, field_str, saveplot, ax):
        slope_rtol = 0.07
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_log, y_log)
        if saveplot:
            # plot points
            ax.plot(x_log, y_log, 'k.')
            x_min = x_log.min()
            x_max = x_log.max()
            offset = 0.05*(x_max - x_min)
            npoints = 50
            xx = np.linspace(x_min - offset, x_max + offset, npoints)
            yy = intercept + slope*xx
            # plot line
            ax.plot(xx, yy, linestyle='--', linewidth=0.5, color='k')
            ax.text(xx[2*int(npoints/3)], yy[2*int(npoints/3)], '{:4.2f}'.format(slope),
                    verticalalignment='top',
                    horizontalalignment='left')
            ax.set_xlabel('log10(dx)')
            ax.set_ylabel('log10(L2 error)')
            #ax.set_title(' '.join([setup_name, field_str, 'degree={:}'.format(polynomial_degree), space_str]))
            ax.set_title(field_str)
        if expected_slope is not None:
            err_msg = '{:}: Wrong convergence rate {:.4f}, expected {:.4f}'.format(setup_name, slope, expected_slope)
            assert slope > expected_slope*(1 - slope_rtol), err_msg
            print_output('{:}: convergence rate {:.4f} PASSED'.format(setup_name, slope))
        else:
            print_output('{:}: {:} convergence rate {:.4f}'.format(setup_name, field_str, slope))
        return slope

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    try:
        check_convergence(x_log, y_log_uv, polynomial_degree+1, 'Velocity', saveplot, ax)
    except Exception as e:
        print(e)

    ref_str = 'ref-' + '-'.join([str(r) for r in ref_list])
    degree_str = 'o{:}'.format(polynomial_degree)
    imgfile = '_'.join(['convergence', setup_name, ref_str, degree_str, space_str])
    imgfile += '.png'
    imgdir = create_directory('plots')
    imgfile = os.path.join(imgdir, imgfile)
    print_output('saving figure {:}'.format(imgfile))
    plt.savefig(imgfile, dpi=200, bbox_inches='tight')


#run(setup1, refinement=1, degree=1)
run_convergence(setup1, [1, 2, 4], saveplot=True, polynomial_degree=1, element_family='dg-dg', no_exports=True)

