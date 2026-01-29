def AddEquationsAero(aero_eq_manager):
    """
    Adds the standard set of equations for aerodynamic loads and vehicle dynamics
    to the provided Equation_Manager instance.
    """


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "p{{i}} + 1/2 * rho{{i}} * v{{i}}**2 + rho{{i}} * g{{i}} * h{{i}} = p{{j}} + 1/2 * rho{{j}} * v{{j}}**2 + rho{{j}} * g{{j}} * h{{j}}"
        ],
        relevant_vars=[
            ("p",   "Static pressure at point [Pa]"),
            ("rho", "Fluid density at point [kg/m³]"),
            ("v",   "Flow velocity at point [m/s]"),
            ("g",   "Gravitational acceleration [m/s²]"),
            ("h",   "Geometric height of point [m]")
        ],
        vars_to_check=[(["p", "rho", "v"], 2)],
        name="Bernoulli Equation (two points)"
    )

    # aero_eq_manager.add_equation_template(
    #     equations_to_add=[
    #         "rho{{i}} * A{{i}} * v{{i}} = rho{{j}} * A{{j}} * v{{j}}"
    #     ],
    #     relevant_vars=[
    #         ("rho", "Fluid density [kg/m³]"),
    #         ("A",   "Cross-sectional area [m²]"),
    #         ("v",   "Flow velocity [m/s]")
    #     ],
    #     vars_to_check=[(["rho", "A", "v"], 2)],
    #     name="Continuity Equation (mass conservation)"
    # )

    #v/v=A/A continuity equation for incompressible flow
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "v{{i}} / v{{j}} = A{{j}} / A{{i}}"
        ],
        relevant_vars=[
            ("A",   "Cross-sectional area [m²]"),
            ("v",   "Flow velocity [m/s]")
        ],
        vars_to_check=[(["A", "v"], 2)],
        name="Continuity Equation for Incompressible Flow"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Re{{i}} = rho{{i}} * v{{i}} * L{{i}} / eta{{i}}"
        ],
        relevant_vars=[
            ("Re",  "Reynolds number [-]"),
            ("rho", "Fluid density [kg/m³]"),
            ("v",   "Characteristic velocity [m/s]"),
            ("L",   "Characteristic length [m]"),
            ("eta", "Dynamic viscosity [Pa·s]")
        ],
        vars_to_check=[(["Re", "rho", "v", "L"], 3)],
        name="Reynolds Number"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "nu{{i}} = eta{{i}} / rho{{i}}"
        ],
        relevant_vars=[
            ("nu",  "Kinematic viscosity [m²/s]"),
            ("eta", "Dynamic viscosity [Pa·s]"),
            ("rho", "Fluid density [kg/m³]")
        ],
        vars_to_check=[(["nu", "eta", "rho"], 2)],
        name="Kinematic Viscosity Definition"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_d{{i}} = 1/2 * rho{{i}} * C_d{{i}} * A_ref{{i}} * v{{i}}**2"
        ],
        relevant_vars=[
            ("F_d",     "Aerodynamic drag force [N]"),
            ("rho",     "Air density [kg/m³]"),
            ("C_d",     "Drag coefficient [-]"),
            ("A_ref",   "Aerodynamic reference area [m²]"),
            ("v",       "Vehicle velocity relative to air [m/s]")
        ],
        vars_to_check=[(["F_d", "rho", "C_d", "A_ref", "v"], 3)],
        name="Aerodynamic Drag Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_y_aero{{i}} = 1/2 * rho{{i}} * C_y{{i}} * A_ref{{i}} * v{{i}}**2"
        ],
        relevant_vars=[
            ("F_y_aero", "Aerodynamic side force [N]"),
            ("rho",      "Air density [kg/m³]"),
            ("C_y",      "Side force coefficient [-]"),
            ("A_ref",    "Aerodynamic reference area [m²]"),
            ("v",        "Vehicle speed relative to air [m/s]")
        ],
        vars_to_check=[(["F_y_aero","C_y"], 1)],
        name="Aerodynamic Lateral Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "M_yaw_aero{{i}} = F_y_aero{{i}} * x_rel_aero{{i}}"
        ],
        relevant_vars=[
            ("M_yaw_aero", "Aerodynamic yawing moment about CG [N·m]"),
            ("F_y_aero",   "Aerodynamic side force [N]"),
            ("x_rel_aero", "Longitudinal distance from CG to aero center [m]")
        ],
        vars_to_check=[(["M_yaw_aero", "F_y_aero", "x_rel_aero"], 2)],
        name="Aerodynamic Yawing Moment"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "v{{i}} = sqrt(v_x{{i}}**2 + v_y{{i}}**2)"
        ],
        relevant_vars=[
            ("v",   "Vehicle speed magnitude [m/s]"),
            ("v_x", "Longitudinal velocity component [m/s]"),
            ("v_y", "Lateral velocity component [m/s]")
        ],
        vars_to_check=[(["v", "v_x", "v_y"], 2)],
        name="Vehicle Speed from Velocity Components"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "beta_sideslip{{i}} = atan(v_y{{i}} / v_x{{i}})"
        ],
        relevant_vars=[
            ("beta_sideslip", "Vehicle sideslip angle at CG [rad]"),
            ("v_x",  "Longitudinal velocity [m/s]"),
            ("v_y",  "Lateral velocity [m/s]")
        ],
        vars_to_check=[(["beta_sideslip", "v_x", "v_y"], 2)],
        name="Vehicle Slip Angle at CG"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_x_drive{{i}} - F_d{{i}} - F_rr{{i}} - F_x_inertia{{i}} = 0"
        ],
        relevant_vars=[
            ("F_x_drive",   "Total driven longitudinal force at wheels [N]"),
            ("F_d",         "Aerodynamic drag force [N]"),
            ("F_rr",        "Rolling resistance force [N]"),
            ("F_x_inertia", "Longitudinal inertial force [N]")
        ],
        vars_to_check=[(["F_x_drive", "F_rr", "F_x_inertia"], 2)],
        name="Longitudinal Force Equilibrium"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_rr{{i}} = C_rr{{i}} * F_z{{i}}"
        ],
        relevant_vars=[
            ("F_rr", "Rolling resistance force [N]"),
            ("C_rr", "Rolling resistance coefficient [-]"),
            ("F_z",  "Generic normal force [N]")
        ],
        vars_to_check=[(["C_rr", "F_z"], 1)],
        name="Rolling Resistance Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "P_x{{i}} = F_x_drive{{i}} * v{{i}}"
        ],
        relevant_vars=[
            ("P_x",       "Longitudinal power at wheels [W]"),
            ("F_x_drive", "Driven longitudinal force [N]"),
            ("v",         "Vehicle speed [m/s]")
        ],
        vars_to_check=[(["P_x", "F_x_drive", "v"], 2)],
        name="Longitudinal Power from Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "M_z{{i}} = I_z{{i}} * r_dot{{i}}"
        ],
        relevant_vars=[
            ("M_z", "Sum of yaw moments about CG [N·m]"),
            ("I_z",     "Yaw moment of inertia [kg·m²]"),
            ("r_dot",   "Yaw acceleration [rad/s²]")
        ],
        vars_to_check=[(["M_z", "I_z", "r_dot"], 2)],
        name="Yaw Moment Equation of Motion"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_x_max{{i}} = mu_long{{i}} * F_z{{i}}"
        ],
        relevant_vars=[
            ("F_x_max", "Maximum longitudinal tyre force [N]"),
            ("mu_long", "Longitudinal friction coefficient [-]"),
            ("F_z",     "Generic normal force on tyre [N]")
        ],
        vars_to_check=[(["F_x_max", "mu_long", "F_z"], 2)],
        name="Maximum Longitudinal Tyre Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_f{{i}} = F_z_fl{{i}} + F_z_fr{{i}}",
            "F_z_r{{i}} = F_z_rl{{i}} + F_z_rr{{i}}"
        ],
        relevant_vars=[
            ("F_z_f",  "Front axle normal force [N]"),
            ("F_z_r",  "Rear axle normal force [N]"),
            ("F_z_fl", "Front-left wheel normal force [N]"),
            ("F_z_fr", "Front-right wheel normal force [N]"),
            ("F_z_rl", "Rear-left wheel normal force [N]"),
            ("F_z_rr", "Rear-right wheel normal force [N]")
        ],
        vars_to_check=[(["F_z_f", "F_z_fl", "F_z_fr"], 2),(["F_z_r", "F_z_rl", "F_z_rr"], 2)],
        name="Axle Normal Force from Wheel Loads"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "p{{i}} = rho{{i}} * R{{i}} * T{{i}}"
        ],
        relevant_vars=[
            ("p",   "Static pressure [Pa]"),
            ("rho", "Air density [kg/m³]"),
            ("R",   "Specific gas constant for air [J/(kg·K)]"),
            ("T",   "Absolute temperature [K]")
        ],
        vars_to_check=[(["p", "rho", "T", "R"], 2)],
        name="Ideal Gas Law"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "T{{i}} = T_0{{i}} - L_ISA{{i}} * h{{i}}"
        ],
        relevant_vars=[
            ("T",      "Ambient temperature [K]"),
            ("T_0",    "Sea level temperature [K]"),
            ("L_ISA",  "ISA temperature lapse rate [K/m]"),
            ("h",      "Altitude [m]")
        ],
        vars_to_check=[(["T", "T_0", "L_ISA", "h"], 3)],
        name="ISA Temperature Lapse Rate (Assumption)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "C_d_sphere{{i}} = 24/Re{{i}} + 4/(0.4 + (Re{{i}}/282000)**0.65) + 0.4/(1 + (Re{{i}}/2.8e6)**0.65)"
        ],
        relevant_vars=[
            ("C_d_sphere", "Drag coefficient of a sphere [-]"),
            ("Re",         "Reynolds number [-]")
        ],
        vars_to_check=[(["C_d_sphere", "Re"], 1)],
        name="Sphere Drag Coefficient (Re-dependent)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "v_terminal{{i}} = sqrt((2*m{{i}}*g{{i}})/(rho{{i}}*C_d_sphere{{i}}*A_ref{{i}}))"
        ],
        relevant_vars=[
            ("v_terminal", "Terminal velocity [m/s]"),
            ("m",          "Mass [kg]"),
            ("g",          "Gravitational acceleration [m/s²]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d_sphere", "Sphere drag coefficient [-]"),
            ("A_ref",      "Reference area [m²]")
        ],
        vars_to_check=[(["v_terminal", "m", "g", "rho", "C_d_sphere", "A_ref"], 4)],
        name="Terminal Velocity (Free Fall with Drag)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "t_fall{{i}} = (m{{i}}/(rho{{i}}*C_d_sphere{{i}}*A_ref{{i}})) * log(v_terminal{{i}}/(v_terminal{{i}} - sqrt(2*g{{i}}*h{{i}})))"
        ],
        relevant_vars=[
            ("t_fall",     "Time to fall from height h with drag [s]"),
            ("m",          "Mass [kg]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d_sphere", "Sphere drag coefficient [-]"),
            ("A_ref",      "Reference area [m²]"),
            ("v_terminal", "Terminal velocity [m/s]"),
            ("g",          "Gravitational acceleration [m/s²]"),
            ("h",          "Fall height [m]")
        ],
        vars_to_check=[(["t_fall", "m", "rho", "C_d_sphere", "A_ref", "v_terminal", "g", "h"], 6)],
        name="Time to Fall from Height with Drag"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "t_coast{{i}} = (m{{i}}/(0.5*rho{{i}}*C_d{{i}}*A_ref{{i}}*v{{i}}))"
            " * log(v{{i}}/v_final{{i}})"
        ],
        relevant_vars=[
            ("t_coast",  "Time to coast from v to v_final [s]"),
            ("m",        "Vehicle mass [kg]"),
            ("rho",      "Air density [kg/m³]"),
            ("C_d",      "Drag coefficient [-]"),
            ("A_ref",    "Reference area [m²]"),
            ("v",        "Initial speed [m/s]"),
            ("v_final",  "Final speed [m/s]")
        ],
        vars_to_check=[(["t_coast", "m", "rho", "C_d", "A_ref", "v", "v_final"], 5)],
        name="Coast-Down Time (Drag Only)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "s_coast{{i}} = (m{{i}}/(0.5*rho{{i}}*C_d{{i}}*A_ref{{i}}))"
            " * (v{{i}} - v_final{{i}})"
        ],
        relevant_vars=[
            ("s_coast", "Coast-down distance [m]"),
            ("m",       "Vehicle mass [kg]"),
            ("rho",     "Air density [kg/m³]"),
            ("C_d",     "Drag coefficient [-]"),
            ("A_ref",   "Reference area [m²]"),
            ("v",       "Initial speed [m/s]"),
            ("v_final", "Final speed [m/s]")
        ],
        vars_to_check=[(["s_coast", "m", "rho", "C_d", "A_ref", "v", "v_final"], 5)],
        name="Coast-Down Distance (Drag Only)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_thrust_req{{i}} = F_d{{i}} + F_x_inertia{{i}}"
        ],
        relevant_vars=[
            ("F_thrust_req", "Required thrust force [N]"),
            ("F_d",          "Aerodynamic drag force [N]"),
            ("F_x_inertia",  "Longitudinal inertial force [N]")
        ],
        vars_to_check=[(["F_thrust_req", "F_d", "F_x_inertia"], 2)],
        name="Required Thrust vs Speed"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "P_required{{i}} = F_thrust_req{{i}} * v{{i}}"
        ],
        relevant_vars=[
            ("P_required",  "Required propulsion power [W]"),
            ("F_thrust_req","Required thrust force [N]"),
            ("v",           "Vehicle speed [m/s]")
        ],
        vars_to_check=[(["P_required", "F_thrust_req", "v"], 2)],
        name="Required Power vs Speed"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "rho{{i}} * A{{i}} * v{{i}} = rho{{j}} * A{{j}} * v{{j}}"
        ],
        relevant_vars=[
            ("rho", "Fluid density [kg/m³]"),
            ("A",   "Flow area [m²]"),
            ("v",   "Flow velocity [m/s]")
        ],
        vars_to_check=[(["rho", "A", "v"], 2)],
        name="Continuity Equation (Mass Conservation)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_CG{{i}} = x_CG{{i}}",
            "b_CG{{i}} = l_wb{{i}} - x_CG{{i}}"
        ],
        relevant_vars=[
            ("a_CG", "CG distance to front axle [m]"),
            ("b_CG", "CG distance to rear axle [m]"),
            ("x_CG", "CG position from front axle [m]"),
            ("l_wb", "Wheelbase [m]")
        ],
        vars_to_check=[(["a_CG", "x_CG"], 1)],
        name="CG Distance Definitions (Front / Rear)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "k_tire{{i}} = F_z{{i}} / delta_z_tire{{i}}"
        ],
        relevant_vars=[
            ("k_tire",       "Tyre vertical stiffness [N/m]"),
            ("F_z",          "Generic tyre normal force [N]"),
            ("delta_z_tire", "Tyre vertical deflection [m]")
        ],
        vars_to_check=[(["k_tire", "F_z", "delta_z_tire"], 2)],
        name="Tyre Vertical Stiffness Definition"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_delta_z_tire{{i}} = Delta_F_z{{i}} / k_tire{{i}}"
        ],
        relevant_vars=[
            ("Delta_delta_z_tire", "Change in tyre deflection [m]"),
            ("Delta_F_z",          "Change in tyre normal load [N]"),
            ("k_tire",             "Tyre vertical stiffness [N/m]")
        ],
        vars_to_check=[(["Delta_delta_z_tire", "Delta_F_z", "k_tire"], 2)],
        name="Tyre Deflection Change from Load Change"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "z_wheel{{i}} = z_chassis{{i}} - delta_z_tire{{i}}"
        ],
        relevant_vars=[
            ("z_wheel",     "Wheel center vertical position [m]"),
            ("z_chassis",  "Chassis reference vertical position [m]"),
            ("delta_z_tire","Tyre vertical deflection [m]")
        ],
        vars_to_check=[(["z_wheel", "z_chassis", "delta_z_tire"], 2)],
        name="Wheel Vertical Position from Tyre Deflection"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_z_spring{{i}} = MR{{i}} * delta_z_wheel{{i}}"
        ],
        relevant_vars=[
            ("delta_z_spring", "Suspension spring compression [m]"),
            ("MR",             "Suspension motion ratio (spring / wheel) [-]"),
            ("delta_z_wheel",  "Wheel vertical displacement [m]")
        ],
        vars_to_check=[(["delta_z_spring", "MR", "delta_z_wheel"], 2)],
        name="Suspension Motion Ratio (Displacement)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_spring{{i}} = k_spring{{i}} * delta_z_spring{{i}}"
        ],
        relevant_vars=[
            ("F_spring",   "Suspension spring force [N]"),
            ("k_spring",   "Suspension spring stiffness [N/m]"),
            ("delta_z_spring", "Spring compression [m]")
        ],
        vars_to_check=[(["F_spring", "k_spring", "delta_z_spring"], 2)],
        name="Suspension Spring Force"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_susp{{i}} = F_spring{{i}} / MR{{i}}"
        ],
        relevant_vars=[
            ("F_z_susp", "Vertical force transmitted by suspension at wheel [N]"),
            ("F_spring", "Suspension spring force [N]"),
            ("MR",       "Suspension motion ratio [-]")
        ],
        vars_to_check=[(["F_z_susp", "F_spring", "MR"], 2)],
        name="Suspension Force at Wheel"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_z_wheel{{i}} = delta_z_tire{{i}} + delta_z_susp{{i}}"
        ],
        relevant_vars=[
            ("delta_z_wheel", "Total wheel vertical displacement [m]"),
            ("delta_z_tire",  "Tyre vertical deflection [m]"),
            ("delta_z_susp",  "Suspension vertical deflection [m]")
        ],
        vars_to_check=[(["delta_z_wheel", "delta_z_tire", "delta_z_susp"], 2)],
        name="Wheel Vertical Displacement Decomposition"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_z_susp{{i}} = F_z{{i}} * MR{{i}} / k_spring{{i}}"
        ],
        relevant_vars=[
            ("delta_z_susp", "Suspension deflection [m]"),
            ("F_z",          "Generic wheel normal force [N]"),
            ("MR",           "Suspension motion ratio [-]"),
            ("k_spring",     "Suspension spring stiffness [N/m]")
        ],
        vars_to_check=[(["delta_z_susp", "F_z", "MR", "k_spring"], 3)],
        name="Suspension Deflection from Wheel Load"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "k_wheel{{i}} = k_spring{{i}} / MR{{i}}**2"
        ],
        relevant_vars=[
            ("k_wheel",  "Effective wheel rate [N/m]"),
            ("k_spring", "Suspension spring stiffness [N/m]"),
            ("MR",       "Suspension motion ratio [-]")
        ],
        vars_to_check=[(["k_wheel", "k_spring", "MR"], 2)],
        name="Effective Wheel Rate"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "T_arb{{i}} = k_arb{{i}} * (delta_z_wheel_left{{i}} - delta_z_wheel_right{{i}})"
        ],
        relevant_vars=[
            ("T_arb",               "Anti-roll bar torque equivalent [N·m]"),
            ("k_arb",               "Anti-roll bar stiffness [N/m]"),
            ("delta_z_wheel_left",  "Left wheel vertical displacement [m]"),
            ("delta_z_wheel_right", "Right wheel vertical displacement [m]")
        ],
        vars_to_check=[(["T_arb", "k_arb", "delta_z_wheel_left", "delta_z_wheel_right"], 3)],
        name="Anti-Roll Bar Torque from Wheel Displacement"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_arb_left{{i}} =  T_arb{{i}} / t_track{{i}}",
            "F_z_arb_right{{i}} = -T_arb{{i}} / t_track{{i}}"
        ],
        relevant_vars=[
            ("F_z_arb_left",  "ARB vertical force at left wheel [N]"),
            ("F_z_arb_right", "ARB vertical force at right wheel [N]"),
            ("T_arb",         "Anti-roll bar torque equivalent [N·m]"),
            ("t_track",       "Track width [m]")
        ],
        vars_to_check=[(["F_z_arb_left", "T_arb", "t_track"], 2)],
        name="Anti-Roll Bar Vertical Wheel Forces"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_z_wheel_f{{i}} = delta_z_heave{{i}} + theta_pitch{{i}} * a_CG{{i}}",
            "delta_z_wheel_r{{i}} = delta_z_heave{{i}} - theta_pitch{{i}} * b_CG{{i}}"
        ],
        relevant_vars=[
            ("delta_z_wheel_f", "Front axle wheel vertical displacement [m]"),
            ("delta_z_wheel_r", "Rear axle wheel vertical displacement [m]"),
            ("delta_z_heave",   "Chassis heave displacement [m]"),
            ("theta_pitch",    "Pitch angle (small-angle) [rad]"),
            ("a_CG",            "CG distance to front axle [m]"),
            ("b_CG",            "CG distance to rear axle [m]")
        ],
        vars_to_check=[(["delta_z_wheel_f", "delta_z_heave", "theta_pitch"], 2)],
        name="Wheel Displacement from Heave and Pitch"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_z_bump{{i}} = delta_z_susp{{i}} - delta_z_clearance{{i}}"
        ],
        relevant_vars=[
            ("delta_z_bump",      "Bump stop compression [m] (≤0 = no contact)"),
            ("delta_z_susp",      "Suspension deflection [m]"),
            ("delta_z_clearance", "Bump stop clearance [m]")
        ],
        vars_to_check=[(["delta_z_bump", "delta_z_susp", "delta_z_clearance"], 2)],
        name="Bump Stop Engagement Deflection"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_bump{{i}} = k_bump{{i}} * delta_z_bump{{i}}"
        ],
        relevant_vars=[
            ("F_bump",     "Bump stop force [N] (active only if delta_z_bump > 0)"),
            ("k_bump",     "Bump stop stiffness [N/m]"),
            ("delta_z_bump","Bump stop compression [m]")
        ],
        vars_to_check=[(["F_bump", "k_bump", "delta_z_bump"], 2)],
        name="Bump Stop Force (Linear, Engaged Only)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_susp_total{{i}} = F_z_susp{{i}} + F_bump{{i}}"
        ],
        relevant_vars=[
            ("F_z_susp_total", "Total suspension vertical force at wheel [N]"),
            ("F_z_susp",       "Spring-based suspension force [N]"),
            ("F_bump",         "Bump stop force [N]")
        ],
        vars_to_check=[(["F_z_susp_total", "F_z_susp"], 1)],
        name="Total Suspension Force Including Bump Stop"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_long_geom{{i}} = F_x_tire{{i}} * tan(theta_anti{{i}})"
        ],
        relevant_vars=[
            ("Delta_Fz_long_geom", "Longitudinal load transfer via suspension geometry [N]"),
            ("F_x_tire",                "Longitudinal tyre force [N]"),
            ("theta_anti",         "Anti-dive / anti-squat angle [rad]")
        ],
        vars_to_check=[(["Delta_Fz_long_geom", "F_x_tire", "theta_anti"], 2)],
        name="Geometric Longitudinal Load Transfer (Anti-Dive / Anti-Squat)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_long_total{{i}} = Delta_Fz_long{{i}} - Delta_Fz_long_geom{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_long_total", "Net longitudinal load transfer [N]"),
            ("Delta_Fz_long",       "Load transfer from CG height [N]"),
            ("Delta_Fz_long_geom",  "Geometric load transfer component [N]")
        ],
        vars_to_check=[(["Delta_Fz_long_total", "Delta_Fz_long"], 1)],
        name="Net Longitudinal Load Transfer Including Geometry"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "1/k_vert{{i}} = 1/k_tire{{i}} + 1/k_wheel{{i}}"
        ],
        relevant_vars=[
            ("k_vert",  "Effective vertical stiffness at wheel [N/m]"),
            ("k_tire",  "Tyre vertical stiffness [N/m]"),
            ("k_wheel", "Suspension wheel rate [N/m]")
        ],
        vars_to_check=[(["k_vert", "k_tire", "k_wheel"], 2)],
        name="Effective Vertical Stiffness (Tyre + Suspension)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "phi_roll_f{{i}} = k_roll_f{{i}} / (k_roll_f{{i}} + k_roll_r{{i}})"
        ],
        relevant_vars=[
            ("phi_roll_f", "Front roll stiffness distribution [-]"),
            ("k_roll_f",   "Front axle roll stiffness [N·m/rad]"),
            ("k_roll_r",   "Rear axle roll stiffness [N·m/rad]")
        ],
        vars_to_check=[(["phi_roll_f", "k_roll_f", "k_roll_r"], 2)],
        name="Roll Stiffness Distribution (Front Fraction)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_inner{{i}} = atan(l_wb{{i}} / (R_turn{{i}} - t_track{{i}}/2))"
        ],
        relevant_vars=[
            ("delta_inner", "Inner wheel steering angle [rad]"),
            ("l_wb",        "Wheelbase [m]"),
            ("R_turn",      "Turn radius to CG [m]"),
            ("t_track",     "Track width [m]")
        ],
        vars_to_check=[(["delta_inner", "l_wb", "R_turn", "t_track"], 3)],
        name="Ackermann Inner Wheel Steering Angle"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_outer{{i}} = atan(l_wb{{i}} / (R_turn{{i}} + t_track{{i}}/2))"
        ],
        relevant_vars=[
            ("delta_outer", "Outer wheel steering angle [rad]"),
            ("l_wb",        "Wheelbase [m]"),
            ("R_turn",      "Turn radius to CG [m]"),
            ("t_track",     "Track width [m]")
        ],
        vars_to_check=[(["delta_outer", "l_wb", "R_turn", "t_track"], 3)],
        name="Ackermann Outer Wheel Steering Angle"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "delta_avg{{i}} = (delta_inner{{i}} + delta_outer{{i}}) / 2"
        ],
        relevant_vars=[
            ("delta_avg",   "Average front steering angle [rad]"),
            ("delta_inner", "Inner wheel steering angle [rad]"),
            ("delta_outer", "Outer wheel steering angle [rad]")
        ],
        vars_to_check=[(["delta_avg", "delta_inner", "delta_outer"], 2)],
        name="Average Steering Angle (Bicycle Equivalent)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "R_turn{{i}} = l_wb{{i}} / tan(delta_avg{{i}})"
        ],
        relevant_vars=[
            ("R_turn",   "Turn radius to CG [m]"),
            ("l_wb",     "Wheelbase [m]"),
            ("delta_avg","Average steering angle [rad]")
        ],
        vars_to_check=[(["R_turn", "l_wb", "delta_avg"], 2)],
        name="Turn Radius from Steering Angle"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_y_available{{i}} = F_z{{i}} * (sin(theta_bank{{i}}) + mu{{i}} * cos(theta_bank{{i}}))"
        ],
        relevant_vars=[
            ("F_y_available", "Maximum available lateral force on banked surface [N]"),
            ("F_z",           "Generic normal force [N]"),
            ("mu",            "Tyre-road friction coefficient [-]"),
            ("theta_bank",    "Bank angle [rad]")
        ],
        vars_to_check=[(["F_y_available", "F_z", "mu", "theta_bank"], 2)],
        name="Available Lateral Force on Banked Surface (for flat roads, theta_bank=0)"
    )

    #note: this not generally true. however, F_y_available is only defined if we also have the bank angle etc. therefore, this is acceptable for our use case.
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_centripetal{{i}} = F_y_available{{i}}"
        ],
        relevant_vars=[
            ("F_centripetal", "Required centripetal force [N]"),
            ("F_y_available", "Available lateral force [N]")
        ],
        vars_to_check=[(["F_centripetal", "F_y_available"], 2)],
        name="Banked Turn Lateral Force Limit"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_df_f{{i}} = balance_df{{i}} * F_z_df{{i}}",
            "F_z_df_r{{i}} = (1 - balance_df{{i}}) * F_z_df{{i}}"
        ],
        relevant_vars=[
            ("F_z_df_f",   "Front axle aerodynamic downforce [N]"),
            ("F_z_df_r",   "Rear axle aerodynamic downforce [N]"),
            ("balance_df", "Aerodynamic balance (front fraction) [-]"),
            ("F_z_df",     "Total aerodynamic downforce [N]")
        ],
        vars_to_check=[(["F_z_df_f", "F_z_df", "balance_df"], 2)],
        name="Aerodynamic Downforce Axle Split"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "eta_aero{{i}} = C_L{{i}} / C_d{{i}}"
        ],
        relevant_vars=[
            ("eta_aero", "Aerodynamic efficiency (lift-to-drag ratio) [-]"),
            ("C_L",      "Lift coefficient [-]"),
            ("C_d",      "Drag coefficient [-]")
        ],
        vars_to_check=[(["eta_aero", "C_L", "C_d"], 2)],
        name="Aerodynamic Efficiency (L/D)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=["p{{i}}=p_0{{i}}*exp(-rho_0{{i}}*g{{i}}*h{{i}}/(p_0{{i}}))"],
        relevant_vars=[
            ("p", "Pressure at altitude h [Pa]"),
            ("p_0", "Sea level standard atmospheric pressure [Pa]"),
            ("rho_0", "Sea level standard density [kg/m³]"),
            ("g", "Gravitational acceleration [m/s²]"),
            ("h", "Altitude [m]")
        ],
        vars_to_check=[(["p", "p_0", "rho_0", "g", "h"], 3)],
        name="Barometric Formula"
    )
    aero_eq_manager.add_equation_template(
        equations_to_add=["rho{{i}}=rho_0{{i}}*exp(-rho_0{{i}}*g{{i}}*h{{i}}/(p_0{{i}}))"],
        relevant_vars=[
            ("rho", "Density at altitude h [kg/m³]"),
            ("rho_0", "Sea level standard density [kg/m³]"),
            ("g", "Gravitational acceleration [m/s²]"),
            ("h", "Altitude [m]"),
            ("p_0", "Sea level standard atmospheric pressure [Pa]")
        ],
        vars_to_check=[(["rho", "rho_0", "g", "h", "p_0"], 3)],
        name="Density Variation with Altitude"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "x_rel_aero{{i}} = x_COP{{i}} - x_CG{{i}}",
            "z_rel_aero{{i}} = z_COP{{i}} - z_CG{{i}}"
        ],
        relevant_vars=[
            ("x_rel_aero", "Longitudinal lever arm from CG to aerodynamic center [m]"),
            ("z_rel_aero", "Vertical lever arm from CG to aerodynamic center [m]"),
            ("x_COP",      "Aerodynamic center / center of pressure position [m]"),
            ("z_COP",      "Aerodynamic center height [m]"),
            ("x_CG",       "Center of gravity longitudinal position [m]"),
            ("z_CG",       "Center of gravity height [m]")
        ],
        vars_to_check=[(["x_rel_aero", "x_COP", "x_CG"], 2),(["z_rel_aero", "z_COP", "z_CG"], 2)],
        name="Aerodynamic Lever Arms relative to CG"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "M_pitch_aero{{i}} = -F_z_df{{i}} * x_rel_aero{{i}} + F_d{{i}} * z_rel_aero{{i}}"
        ],
        relevant_vars=[
            ("M_pitch_aero", "Total aerodynamic pitching moment about CG [N·m]"),
            ("F_z_df",       "Total aerodynamic normal force (downforce) [N]"),
            ("F_d",          "Aerodynamic drag force [N]"),
            ("x_rel_aero",   "Longitudinal lever arm from CG to COP [m]"),
            ("z_rel_aero",   "Vertical lever arm from CG to COP [m]")
        ],
        vars_to_check=[(["M_pitch_aero", "F_z_df", "x_rel_aero"], 2)],
        name="Aerodynamic Pitching Moment about CG"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "balance_df{{i}} = (l_wb{{i}} - x_COP{{i}}) / l_wb{{i}}"
        ],
        relevant_vars=[
            ("balance_df", "Aerodynamic balance (front fraction of downforce) [-]"),
            ("x_COP",      "Center of pressure position from front axle [m]"),
            ("l_wb",       "Wheelbase [m]")
        ],
        vars_to_check=[(["balance_df", "x_COP", "l_wb"], 2)],
        name="Aerodynamic Balance from COP Location"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "M_pitch_df{{i}} = F_z_df{{i}} * (x_CG{{i}} - x_COP{{i}})"
        ],
        relevant_vars=[
            ("M_pitch_df", "Pitching moment caused by downforce distribution [N·m]"),
            ("F_z_df",     "Total aerodynamic downforce [N]"),
            ("x_CG",       "Center of gravity position [m]"),
            ("x_COP",      "Center of pressure position [m]")
        ],
        vars_to_check=[(["M_pitch_df", "F_z_df", "x_CG", "x_COP"], 3)],
        name="Pitching Moment from Downforce Distribution"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "M_pitch_drag{{i}} = F_d{{i}} * (z_CG{{i}} - z_COP{{i}})"
        ],
        relevant_vars=[
            ("M_pitch_drag", "Pitching moment caused by drag force [N·m]"),
            ("F_d",          "Aerodynamic drag force [N]"),
            ("z_CG",         "Center of gravity height [m]"),
            ("z_COP",        "Aerodynamic center height [m]")
        ],
        vars_to_check=[(["M_pitch_drag", "F_d", "z_CG", "z_COP"], 3)],
        name="Pitching Moment from Drag Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_x_inertia{{i}} = m{{i}} * a_long{{i}}"
        ],
        relevant_vars=[
            ("F_x_inertia", "Longitudinal inertial force (vehicle frame) [N]"),
            ("m",           "Vehicle mass [kg]"),
            ("a_long",      "Longitudinal acceleration (positive forward) [m/s²]")
        ],
        vars_to_check=[(["F_x_inertia", "m", "a_long"], 2)],
        name="Longitudinal Inertial Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_long{{i}} = F_x_inertia{{i}} * z_CG{{i}} / l_wb{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_long", "Total longitudinal normal load transfer [N]"),
            ("F_x_inertia",   "Longitudinal inertial force [N]"),
            ("z_CG",          "Center of gravity height [m]"),
            ("l_wb",          "Wheelbase [m]")
        ],
        vars_to_check=[(["Delta_Fz_long", "F_x_inertia", "z_CG", "l_wb"], 3)],
        name="Total Longitudinal Load Transfer (mass based)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_long_f{{i}} = -Delta_Fz_long{{i}}",
            "Delta_Fz_long_r{{i}} =  Delta_Fz_long{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_long_f", "Front axle normal load change from longitudinal effects [N]"),
            ("Delta_Fz_long_r", "Rear axle normal load change from longitudinal effects [N]"),
            ("Delta_Fz_long",   "Total longitudinal load transfer [N]")
        ],
        vars_to_check=[(["Delta_Fz_long_f", "Delta_Fz_long"], 1)],
        name="Longitudinal Load Transfer Axle Split"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_long_fl{{i}} = Delta_Fz_long_f{{i}} / 2",
            "Delta_Fz_long_fr{{i}} = Delta_Fz_long_f{{i}} / 2",
            "Delta_Fz_long_rl{{i}} = Delta_Fz_long_r{{i}} / 2",
            "Delta_Fz_long_rr{{i}} = Delta_Fz_long_r{{i}} / 2"
        ],
        relevant_vars=[
            ("Delta_Fz_long_fl", "Front-left load change from longitudinal effects [N]"),
            ("Delta_Fz_long_fr", "Front-right load change from longitudinal effects [N]"),
            ("Delta_Fz_long_rl", "Rear-left load change from longitudinal effects [N]"),
            ("Delta_Fz_long_rr", "Rear-right load change from longitudinal effects [N]")
        ],
        vars_to_check=[(["Delta_Fz_long_fl", "Delta_Fz_long_f"], 1)],
        name="Longitudinal Load Transfer Per Wheel"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_y_inertia{{i}} = m{{i}} * a_lat{{i}}"
        ],
        relevant_vars=[
            ("F_y_inertia", "Lateral inertial force (vehicle frame) [N]"),
            ("m",           "Vehicle mass [kg]"),
            ("a_lat",       "Lateral acceleration [m/s²]")
        ],
        vars_to_check=[(["F_y_inertia", "m", "a_lat"], 2)],
        name="Lateral Inertial Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_lat{{i}} = v{{i}}**2 / R_turn{{i}}"
        ],
        relevant_vars=[
            ("a_lat",   "Lateral acceleration [m/s²]"),
            ("v",       "Vehicle speed [m/s]"),
            ("R_turn",  "Turn radius to CG [m]")
        ],
        vars_to_check=[(["a_lat", "v", "R_turn"], 2)],
        name="Lateral Acceleration in a Turn"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_lat{{i}} = F_y_inertia{{i}} * z_CG{{i}} / t_track{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_lat", "Total lateral normal load transfer [N]"),
            ("F_y_inertia",  "Lateral inertial force [N]"),
            ("z_CG",         "Center of gravity height [m]"),
            ("t_track",      "Track width [m]")
        ],
        vars_to_check=[(["Delta_Fz_lat", "F_y_inertia", "z_CG", "t_track"], 3)],
        name="Total Lateral Load Transfer"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_lat_f{{i}} = balance_lat{{i}} * Delta_Fz_lat{{i}}",
            "Delta_Fz_lat_r{{i}} = (1 - balance_lat{{i}}) * Delta_Fz_lat{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_lat_f", "Front axle lateral load transfer [N]"),
            ("Delta_Fz_lat_r", "Rear axle lateral load transfer [N]"),
            ("balance_lat",    "Lateral load transfer balance (front fraction) [-]"),
            ("Delta_Fz_lat",   "Total lateral load transfer [N]")
        ],
        vars_to_check=[(["Delta_Fz_lat_f", "Delta_Fz_lat", "balance_lat"], 2)],
        name="Lateral Load Transfer Axle Distribution"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_lat_fl{{i}} =  Delta_Fz_lat_f{{i}} / 2",
            "Delta_Fz_lat_fr{{i}} = -Delta_Fz_lat_f{{i}} / 2",
            "Delta_Fz_lat_rl{{i}} =  Delta_Fz_lat_r{{i}} / 2",
            "Delta_Fz_lat_rr{{i}} = -Delta_Fz_lat_r{{i}} / 2"
        ],
        relevant_vars=[
            ("Delta_Fz_lat_fl", "Front-left load change from lateral effects [N]"),
            ("Delta_Fz_lat_fr", "Front-right load change from lateral effects [N]"),
            ("Delta_Fz_lat_rl", "Rear-left load change from lateral effects [N]"),
            ("Delta_Fz_lat_rr", "Rear-right load change from lateral effects [N]")
        ],
        vars_to_check=[(["Delta_Fz_lat_fl", "Delta_Fz_lat_f"], 1)],
        name="Lateral Load Transfer Per Wheel"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "W{{i}} = m{{i}} * g{{i}}"
        ],
        relevant_vars=[
            ("W", "Vehicle weight (gravitational force) [N]"),
            ("m", "Vehicle mass [kg]"),
            ("g", "Gravitational acceleration [m/s²]")
        ],
        vars_to_check=[(["W", "m", "g"], 2)],
        name="Vehicle Weight Definition"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static{{i}} = W{{i}}"
        ],
        relevant_vars=[
            ("F_z_static", "Total static normal force (gravity only) [N]"),
            ("W",          "Vehicle weight [N]")
        ],
        vars_to_check=[(["F_z_static", "W"], 1)],
        name="Static Normal Force (Total)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static_fl{{i}} = F_z_static_f{{i}} / 2",
            "F_z_static_fr{{i}} = F_z_static_f{{i}} / 2",
            "F_z_static_rl{{i}} = F_z_static_r{{i}} / 2",
            "F_z_static_rr{{i}} = F_z_static_r{{i}} / 2"
        ],
        relevant_vars=[
            ("F_z_static_fl", "Static front-left normal force [N]"),
            ("F_z_static_fr", "Static front-right normal force [N]"),
            ("F_z_static_rl", "Static rear-left normal force [N]"),
            ("F_z_static_rr", "Static rear-right normal force [N]")
        ],
        vars_to_check=[(["F_z_static_fl", "F_z_static_f"], 1)],
        name="Static Normal Force Per Wheel"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_drag_fl{{i}} = Delta_Fz_drag_f{{i}} / 2",
            "Delta_Fz_drag_fr{{i}} = Delta_Fz_drag_f{{i}} / 2",
            "Delta_Fz_drag_rl{{i}} = Delta_Fz_drag_r{{i}} / 2",
            "Delta_Fz_drag_rr{{i}} = Delta_Fz_drag_r{{i}} / 2"
        ],
        relevant_vars=[
            ("Delta_Fz_drag_fl", "Front-left normal load change from drag moment [N]"),
            ("Delta_Fz_drag_fr", "Front-right normal load change from drag moment [N]"),
            ("Delta_Fz_drag_rl", "Rear-left normal load change from drag moment [N]"),
            ("Delta_Fz_drag_rr", "Rear-right normal load change from drag moment [N]")
        ],
        vars_to_check=[(["Delta_Fz_drag_fl", "Delta_Fz_drag_f"], 1)],
        name="Drag-Induced Normal Load Shift Per Wheel"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static_df_f{{i}} = F_z_static_f{{i}} + F_z_df_f{{i}}",
            "F_z_static_df_r{{i}} = F_z_static_r{{i}} + F_z_df_r{{i}}"
        ],
        relevant_vars=[
            ("F_z_static_df_f", "Front axle normal force: static + downforce [N]"),
            ("F_z_static_df_r", "Rear axle normal force: static + downforce [N]")
        ],
        vars_to_check=[(["F_z_static_df_f", "F_z_static_f", "F_z_df_f"], 2)],
        name="Normal Force: Static + Downforce (Axle)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static_df_fl{{i}} = F_z_static_fl{{i}} + F_z_df_fl{{i}}",
            "F_z_static_df_fr{{i}} = F_z_static_fr{{i}} + F_z_df_fr{{i}}",
            "F_z_static_df_rl{{i}} = F_z_static_rl{{i}} + F_z_df_rl{{i}}",
            "F_z_static_df_rr{{i}} = F_z_static_rr{{i}} + F_z_df_rr{{i}}"
        ],
        relevant_vars=[
            ("F_z_static_df_fl", "Front-left normal force: static + downforce [N]"),
            ("F_z_static_df_fr", "Front-right normal force: static + downforce [N]"),
            ("F_z_static_df_rl", "Rear-left normal force: static + downforce [N]"),
            ("F_z_static_df_rr", "Rear-right normal force: static + downforce [N]")
        ],
        vars_to_check=[(["F_z_static_df_fl", "F_z_static_fl", "F_z_df_fl"], 2)],
        name="Normal Force: Static + Downforce (Per Wheel)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static_df_drag_fl{{i}} = F_z_static_df_fl{{i}} + Delta_Fz_drag_fl{{i}}",
            "F_z_static_df_drag_fr{{i}} = F_z_static_df_fr{{i}} + Delta_Fz_drag_fr{{i}}",
            "F_z_static_df_drag_rl{{i}} = F_z_static_df_rl{{i}} + Delta_Fz_drag_rl{{i}}",
            "F_z_static_df_drag_rr{{i}} = F_z_static_df_rr{{i}} + Delta_Fz_drag_rr{{i}}"
        ],
        relevant_vars=[
            ("F_z_static_df_drag_fl", "Front-left normal force: static + downforce + drag moment [N]"),
            ("F_z_static_df_drag_fr", "Front-right normal force: static + downforce + drag moment [N]"),
            ("F_z_static_df_drag_rl", "Rear-left normal force: static + downforce + drag moment [N]"),
            ("F_z_static_df_drag_rr", "Rear-right normal force: static + downforce + drag moment [N]")
        ],
        vars_to_check=[(["F_z_static_df_drag_fl", "F_z_static_df_fl", "Delta_Fz_drag_fl"], 2)],
        name="Normal Force: Static + Downforce + Drag Moment (Per Wheel)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_static_f{{i}} = W{{i}} * (l_wb{{i}} - x_CG{{i}}) / l_wb{{i}}",
            "F_z_static_r{{i}} = W{{i}} * x_CG{{i}} / l_wb{{i}}"
        ],
        relevant_vars=[
            ("F_z_static_f", "Static front axle normal force [N]"),
            ("F_z_static_r", "Static rear axle normal force [N]"),
            ("W",            "Vehicle weight [N]"),
            ("x_CG",         "CG distance from front axle [m]"),
            ("l_wb",         "Wheelbase [m]")
        ],
        vars_to_check=[(["F_z_static_f", "W", "x_CG", "l_wb"], 3)],
        name="Static Normal Force Axle Split"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_df{{i}} = 1/2 * rho{{i}} * C_L{{i}} * A_ref{{i}} * v{{i}}**2"
        ],
        relevant_vars=[
            ("F_z_df", "Total aerodynamic normal force (downforce positive) [N]"),
            ("rho",    "Air density [kg/m³]"),
            ("C_L",    "Lift coefficient [-]"),
            ("A_ref",  "Reference area [m²]"),
            ("v",      "Velocity [m/s]")
        ],
        vars_to_check=[(["F_z_df", "rho", "C_L", "A_ref", "v"], 3)],
        name="Aerodynamic Downforce (Total)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_z_df_fl{{i}} = F_z_df_f{{i}} / 2",
            "F_z_df_fr{{i}} = F_z_df_f{{i}} / 2",
            "F_z_df_rl{{i}} = F_z_df_r{{i}} / 2",
            "F_z_df_rr{{i}} = F_z_df_r{{i}} / 2"
        ],
        relevant_vars=[
            ("F_z_df_fl", "Front-left aerodynamic downforce [N]"),
            ("F_z_df_fr", "Front-right aerodynamic downforce [N]"),
            ("F_z_df_rl", "Rear-left aerodynamic downforce [N]"),
            ("F_z_df_rr", "Rear-right aerodynamic downforce [N]")
        ],
        vars_to_check=[(["F_z_df_fl", "F_z_df_f"], 1)],
        name="Aerodynamic Downforce Per Wheel"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_drag_f{{i}} =  M_pitch_drag{{i}} / l_wb{{i}}",
            "Delta_Fz_drag_r{{i}} = -M_pitch_drag{{i}} / l_wb{{i}}"
        ],
        relevant_vars=[
            ("Delta_Fz_drag_f", "Front axle normal load change from drag moment [N]"),
            ("Delta_Fz_drag_r", "Rear axle normal load change from drag moment [N]"),
            ("M_pitch_drag",    "Pitching moment from drag [N·m]"),
            ("l_wb",            "Wheelbase [m]")
        ],
        vars_to_check=[(["Delta_Fz_drag_f", "M_pitch_drag", "l_wb"], 2)],
        name="Drag-Induced Normal Load Shift (Axle)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "Delta_Fz_drag_fl{{i}} = Delta_Fz_drag_f{{i}} / 2",
            "Delta_Fz_drag_fr{{i}} = Delta_Fz_drag_f{{i}} / 2",
            "Delta_Fz_drag_rl{{i}} = Delta_Fz_drag_r{{i}} / 2",
            "Delta_Fz_drag_rr{{i}} = Delta_Fz_drag_r{{i}} / 2"
        ],
        relevant_vars=[
            ("Delta_Fz_drag_fl", "Front-left normal load change from drag moment [N]"),
            ("Delta_Fz_drag_fr", "Front-right normal load change from drag moment [N]"),
            ("Delta_Fz_drag_rl", "Rear-left normal load change from drag moment [N]"),
            ("Delta_Fz_drag_rr", "Rear-right normal load change from drag moment [N]")
        ],
        vars_to_check=[(["Delta_Fz_drag_fl", "Delta_Fz_drag_f"], 1)],
        name="Drag-Induced Normal Load Shift (Per Wheel)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_f_max{{i}} = mu{{i}} * F_z{{i}}"
        ],
        relevant_vars=[
            ("F_f_max", "Maximum available friction force magnitude [N]"),
            ("mu",      "Friction coefficient [-]"),
            ("F_z",     "Normal force [N]")
        ],
        vars_to_check=[(["F_f_max", "mu", "F_z"], 2)],
        name="Maximum Friction Force"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "mu_load_dep{{i}} = mu_0{{i}} * f_load{{i}}"
        ],
        relevant_vars=[
            ("mu_load_dep", "Load-dependent friction coefficient [-]"),
            ("mu_0",        "Reference friction coefficient [-]"),
            ("f_load",      "Dimensionless load influence factor [-]")
        ],
        vars_to_check=[(["mu_load_dep", "mu_0", "f_load"], 2)],
        name="Load-Dependent Friction Coefficient (Generic Model)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_long_max{{i}} = mu{{i}} * F_z{{i}} / m{{i}}"
        ],
        relevant_vars=[
            ("a_long_max", "Maximum longitudinal acceleration (friction-limited) [m/s²]"),
            ("mu",         "Friction coefficient [-]"),
            ("F_z",        "Normal force [N]"),
            ("m",          "Vehicle mass [kg]")
        ],
        vars_to_check=[(["a_long_max", "mu", "F_z", "m"], 3)],
        name="Maximum Longitudinal Acceleration (Friction-Limited)"
    )


    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_lat_max{{i}} = mu{{i}} * F_z{{i}} / m{{i}}"
        ],
        relevant_vars=[
            ("a_lat_max", "Maximum lateral acceleration (friction-limited) [m/s²]"),
            ("mu",        "Friction coefficient [-]"),
            ("F_z",       "Normal force [N]"),
            ("m",         "Vehicle mass [kg]")
        ],
        vars_to_check=[(["a_lat_max", "mu", "F_z", "m"], 3)],
        name="Maximum Lateral Acceleration (Friction-Limited)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "v_max_corner{{i}} = sqrt(mu{{i}} * F_z{{i}} * R_turn{{i}} / m{{i}})"
        ],
        relevant_vars=[
            ("v_max_corner", "Maximum cornering speed (friction-limited) [m/s]"),
            ("mu",           "Friction coefficient [-]"),
            ("F_z",          "Normal force [N]"),
            ("R_turn",      "Corner radius [m]"),
            ("m",            "Vehicle mass [kg]")
        ],
        vars_to_check=[(["v_max_corner", "mu", "F_z", "R_turn", "m"], 4)],
        name="Maximum Cornering Speed (Friction-Limited)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_tan{{i}} = sqrt(F_x_tire{{i}}**2 + F_y_inertia{{i}}**2)"
        ],
        relevant_vars=[
            ("F_tan", "Resultant tangential force magnitude [N]"),
            ("F_x_tire",   "Longitudinal force [N]"),
            ("F_y_inertia",   "Lateral inertial force [N]")
        ],
        vars_to_check=[(["F_tan", "F_x_tire", "F_y_inertia"], 2)],
        name="Resultant Tangential Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "friction_util{{i}} = F_tan{{i}} / (mu{{i}} * F_z{{i}})"
        ],
        relevant_vars=[
            ("friction_util", "Friction utilisation ratio (≤1 = within grip) [-]"),
            ("F_tan",         "Resultant tangential force [N]"),
            ("mu",            "Friction coefficient [-]"),
            ("F_z",           "Normal force [N]")
        ],
        vars_to_check=[(["friction_util", "F_tan", "mu", "F_z"], 3)],
        name="Friction Circle Utilisation"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "(F_x_tire{{i}} / F_x_tire_max{{i}})**2 + (F_y_inertia{{i}} / F_y_max{{i}})**2 = 1"
        ],
        relevant_vars=[
            ("F_x_tire",     "Longitudinal force [N]"),
            ("F_y_inertia",     "Lateral inertial force [N]"),
            ("F_x_tire_max", "Maximum longitudinal force [N]"),
            ("F_y_max", "Maximum lateral force [N]")
        ],
        vars_to_check=[(["F_x_tire", "F_y_inertia", "F_x_tire_max", "F_y_max"], 3)],
        name="Friction Ellipse Model"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "P_trac{{i}} = F_x_tire{{i}} * v{{i}}"
        ],
        relevant_vars=[
            ("P_trac", "Tractive power at the wheels [W]"),
            ("F_x_tire",    "Longitudinal tractive force [N]"),
            ("v",      "Vehicle speed [m/s]")
        ],
        vars_to_check=[(["P_trac", "F_x_tire", "v"], 2)],
        name="Tractive Power Definition"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_brake_aero{{i}} = mu{{i}} * (F_z_static{{i}} + F_z_df{{i}}) / m{{i}}"
        ],
        relevant_vars=[
            ("a_brake_aero", "Maximum braking deceleration with aerodynamic downforce [m/s²]"),
            ("mu",           "Friction coefficient [-]"),
            ("F_z_static",   "Static normal force [N]"),
            ("F_z_df",       "Aerodynamic downforce [N]"),
            ("m",            "Vehicle mass [kg]")
        ],
        vars_to_check=[(["a_brake_aero", "mu", "F_z_static", "F_z_df", "m"], 4)],
        name="Braking Deceleration with Aerodynamic Downforce"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_centripetal{{i}} = m{{i}} * v{{i}}**2 / R_turn{{i}}"
        ],
        relevant_vars=[
            ("F_centripetal", "Required centripetal (lateral) force [N]"),
            ("m",             "Vehicle mass [kg]"),
            ("v",             "Vehicle speed [m/s]"),
            ("R_turn",       "Curve radius [m]")
        ],
        vars_to_check=[(["F_centripetal", "m", "v", "R_turn"], 3)],
        name="Centripetal Force Requirement"
    )

    
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "F_slope{{i}} = m{{i}} * g{{i}} * sin(theta_slope{{i}})"
        ],
        relevant_vars=[
            ("F_slope",      "Longitudinal force due to road slope [N]"),
            ("m",            "Vehicle mass [kg]"),
            ("g",            "Gravitational acceleration [m/s²]"),
            ("theta_slope",  "Road slope angle (positive uphill) [rad]")
        ],
        vars_to_check=[(["F_slope", "m", "g", "theta_slope"], 3)],
        name="Road Slope Longitudinal Force"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            "a_x_max{{i}} = F_x_tire_max{{i}} / m{{i}}"
        ],
        relevant_vars=[
            ("a_x_max", "Maximum traction-limited longitudinal acceleration [m/s²]"),
            ("F_x_tire_max", "Maximum available longitudinal tyre force [N]"),
            ("m",       "Vehicle mass [kg]")
        ],
        vars_to_check=[(["a_x_max", "F_x_tire_max", "m"], 2)],
        name="Maximum Longitudinal Acceleration (Traction-Limited)"
    )

























