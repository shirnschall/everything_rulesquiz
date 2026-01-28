def AddEquationsAero(aero_eq_manager):
    """
    Adds the standard set of equations for aerodynamic loads and vehicle dynamics
    to the provided Equation_Manager instance.
    """
    # ────────────────────────────────────────────────
    # Fluid dynamics & basic relations
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("Re{{i}}", "L{{i}} * rho{{i}} * v{{i}} / eta{{i}}")],
        relevant_vars=[
            ("Re",  "Reynolds number [-]"),
            ("L",   "Characteristic length [m]"),
            ("rho", "Air density [kg/m³]"),
            ("v",   "Velocity [m/s]"),
            ("eta", "Dynamic viscosity [Pa·s]")
        ],
        vars_to_check=[(["Re", "L"], 1)],
        name="Reynolds Number"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[("nu{{i}}", "eta{{i}} / rho{{i}}")],
        relevant_vars=[
            ("nu",  "Kinematic viscosity [m²/s]"),
            ("eta", "Dynamic viscosity [Pa·s]"),
            ("rho", "Air density [kg/m³]")
        ],
        vars_to_check=[(["nu", "eta"], 1)],
        name="Kinematic Viscosity"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[("v{{i}} / v{{j}}", "A{{j}} / A{{i}}")],
        relevant_vars=[
            ("v", "Velocity [m/s]"),
            ("A", "Cross-sectional area [m²]")
        ],
        vars_to_check=[(["v", "A"], 1)],
        name="Continuity Equation"
    )


    # ────────────────────────────────────────────────
    # Static weight / gravity only
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_z_static{{i}}", "m{{i}} * g{{i}}")],
        relevant_vars=[
            ("F_z_static", "Total vehicle weight [N]"),
            ("m", "Vehicle mass [kg]"),
            ("g", "Gravitational acceleration [m/s²]")
        ],
        vars_to_check=[(["F_z_static", "m", "g"], 2)],
        name="Total Static Weight"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_static_f{{i}}", "m{{i}} * g{{i}} * (l_wy{{i}} - x_CG{{i}}) / l_wy{{i}}"),
            ("F_z_static_r{{i}}", "m{{i}} * g{{i}} * x_CG{{i}} / l_wy{{i}}"),
        ],
        relevant_vars=[
            ("F_z_static_f", "Static front axle load – gravity only [N]"),
            ("F_z_static_r", "Static rear axle load – gravity only [N]"),
            ("m", "Vehicle mass [kg]"),
            ("g", "Gravitational acceleration [m/s²]"),
            ("l_wb", "Wheelbase [m]"),
            ("x_CG", "Longitudinal CG position from front axle [m]")
        ],
        vars_to_check=[(["F_z_static_f", "m", "g", "l_wb", "x_CG"], 4)],
        name="Static Axle Loads"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_static_fl{{i}}", "F_z_static_f{{i}} / 2"),
            ("F_z_static_fr{{i}}", "F_z_static_f{{i}} / 2"),
            ("F_z_static_rl{{i}}", "F_z_static_r{{i}} / 2"),
            ("F_z_static_rr{{i}}", "F_z_static_r{{i}} / 2"),
        ],
        relevant_vars=[
            ("F_z_static_fl", "Static front-left wheel load – gravity only [N]"),
            ("F_z_static_fr", "Static front-right wheel load – gravity only [N]"),
            ("F_z_static_rl", "Static rear-left wheel load – gravity only [N]"),
            ("F_z_static_rr", "Static rear-right wheel load – gravity only [N]"),
        ],
        vars_to_check=[(["F_z_static_fl", "F_z_static_f"], 1)],
        name="Static Per-Wheel Loads (symmetric)"
    )


    # ────────────────────────────────────────────────
    # Pure aerodynamic downforce (from C_L)
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_z_df{{i}}", "1/2 * rho{{i}} * C_L{{i}} * A{{i}} * v{{i}}**2")],
        relevant_vars=[
            ("F_z_df",       "Total pure aerodynamic downforce [N]"),
            ("C_L",          "Lift coefficient (negative = downforce) [-]"),
            ("A",            "Reference area [m²]"),
            ("v",            "Velocity [m/s]"),
            ("rho",          "Air density [kg/m³]")
        ],
        vars_to_check=[(["F_z_df", "C_L"], 1)],
        name="Pure Aerodynamic Downforce"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_df_f{{i}}", "AB{{i}} * F_z_df{{i}}"),
            ("F_z_df_r{{i}}", "(1 - AB{{i}}) * F_z_df{{i}}"),
        ],
        relevant_vars=[
            ("F_z_df_f", "Pure aero downforce – front axle [N] (AB includes df moment effect)"),
            ("F_z_df_r", "Pure aero downforce – rear axle [N]"),
            ("AB",       "Aero balance – front fraction [0–1]"),
            ("F_z_df",   "Total pure aero downforce [N]")
        ],
        vars_to_check=[(["F_z_df_f", "F_z_df", "AB"], 2)],
        name="Pure Downforce – front/rear split (moment-corrected)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_df_fl{{i}}", "F_z_df_f{{i}} / 2"),
            ("F_z_df_fr{{i}}", "F_z_df_f{{i}} / 2"),
            ("F_z_df_rl{{i}}", "F_z_df_r{{i}} / 2"),
            ("F_z_df_rr{{i}}", "F_z_df_r{{i}} / 2"),
        ],
        relevant_vars=[
            ("F_z_df_fl", "Pure aero downforce – front-left wheel [N]"),
            ("F_z_df_fr", "Pure aero downforce – front-right wheel [N]"),
            ("F_z_df_rl", "Pure aero downforce – rear-left wheel [N]"),
            ("F_z_df_rr", "Pure aero downforce – rear-right wheel [N]"),
        ],
        vars_to_check=[(["F_z_df_fl", "F_z_df_f"], 1)],
        name="Pure Downforce – per wheel (symmetric)"
    )


    # ────────────────────────────────────────────────
    # Pitching moments
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("M_pitch_df{{i}}", "F_z_df{{i}} * (x_CG{{i}} - x_COP{{i}})")],
        relevant_vars=[
            ("M_pitch_df", "Pitching moment from downforce distribution [N·m]"),
            ("F_z_df",     "Total pure aero downforce [N]"),
            ("x_CG",       "CG position from front [m]"),
            ("x_COP",      "Center of pressure from front [m]")
        ],
        vars_to_check=[(["M_pitch_df", "F_z_df", "x_CG", "x_COP"], 2)],
        name="Pitching Moment – Downforce Distribution"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[("M_pitch_drag{{i}}", "F_d{{i}} * (z_CG{{i}} - z_drag{{i}})")],
        relevant_vars=[
            ("M_pitch_drag", "Pitching moment from drag height offset [N·m]"),
            ("F_d",          "Aerodynamic drag [N]"),
            ("z_CG",         "CG height [m]"),
            ("z_drag",       "Effective drag application height [m]")
        ],
        vars_to_check=[(["M_pitch_drag", "F_d", "z_CG", "z_drag"], 2)],
        name="Pitching Moment – Drag Offset"
    )


    # ────────────────────────────────────────────────
    # Load transfers from moments
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_df_f{{i}}",  "M_pitch_df{{i}} / l_wb{{i}}"),
            ("Delta_Fz_df_r{{i}}",  "-M_pitch_df{{i}} / l_wb{{i}}"),
        ],
        relevant_vars=[
            ("Delta_Fz_df_f", "Front axle load change from downforce distribution [N]"),
            ("Delta_Fz_df_r", "Rear axle load change from downforce distribution [N]"),
        ],
        vars_to_check=[(["Delta_Fz_df_f", "M_pitch_df", "l_wb"], 2)],
        name="Axle Load Shift – Downforce Distribution"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_drag_f{{i}}",  "M_pitch_drag{{i}} / l_wb{{i}}"),
            ("Delta_Fz_drag_r{{i}}",  "-M_pitch_drag{{i}} / l_wb{{i}}"),
        ],
        relevant_vars=[
            ("Delta_Fz_drag_f", "Front axle load change from drag moment [N]"),
            ("Delta_Fz_drag_r", "Rear axle load change from drag moment [N]"),
        ],
        vars_to_check=[(["Delta_Fz_drag_f", "M_pitch_drag", "l_wb"], 2)],
        name="Axle Load Shift – Drag Moment"
    )


    # ────────────────────────────────────────────────
    # Aero load incl. drag moment only (AB already includes df moment)  ← Option A
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_aero_f{{i}}", "F_z_df_f{{i}} + Delta_Fz_drag_f{{i}}"),
            ("F_z_aero_r{{i}}", "F_z_df_r{{i}} + Delta_Fz_drag_r{{i}}"),
            ("F_z_aero{{i}}",   "F_z_aero_f{{i}} + F_z_aero_r{{i}}"),
        ],
        relevant_vars=[
            ("F_z_aero",   "Total aero load incl. drag pitching moment only [N]"),
            ("F_z_aero_f", "Front axle aero load incl. drag moment [N]"),
            ("F_z_aero_r", "Rear axle aero load incl. drag moment [N]"),
        ],
        vars_to_check=[(["F_z_aero_f", "F_z_df_f"], 1)],
        name="Aero Load incl. Drag Pitching Moment Only (recommended)"
    )


    # ────────────────────────────────────────────────
    # Aero load incl. BOTH downforce distribution moment AND drag moment  ← new full-moments variant
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_aero_full_moments_f{{i}}", "F_z_df_f{{i}} + Delta_Fz_df_f{{i}} + Delta_Fz_drag_f{{i}}"),
            ("F_z_aero_full_moments_r{{i}}", "F_z_df_r{{i}} + Delta_Fz_df_r{{i}} + Delta_Fz_drag_r{{i}}"),
            ("F_z_aero_full_moments{{i}}",   "F_z_aero_full_moments_f{{i}} + F_z_aero_full_moments_r{{i}}"),
        ],
        relevant_vars=[
            ("F_z_aero_full_moments",   "Total aero load incl. both df distribution moment and drag moment [N]"),
            ("F_z_aero_full_moments_f", "Front axle aero load incl. both moments [N]"),
            ("F_z_aero_full_moments_r", "Rear axle aero load incl. both moments [N]"),
        ],
        vars_to_check=[(["F_z_aero_full_moments_f", "F_z_df_f"], 1)],
        name="Aero Load incl. Full Pitching Moments (both df and drag)"
    )

    # Per-wheel versions for both aero variants (symmetric left/right)

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_aero_fl{{i}}", "F_z_aero_f{{i}} / 2"),
            ("F_z_aero_fr{{i}}", "F_z_aero_f{{i}} / 2"),
            ("F_z_aero_rl{{i}}", "F_z_aero_r{{i}} / 2"),
            ("F_z_aero_rr{{i}}", "F_z_aero_r{{i}} / 2"),
        ],
        relevant_vars=[
            ("F_z_aero_fl", "Aero load incl. drag moment – front-left [N]"),
            ("F_z_aero_fr", "Aero load incl. drag moment – front-right [N]"),
            ("F_z_aero_rl", "Aero load incl. drag moment – rear-left [N]"),
            ("F_z_aero_rr", "Aero load incl. drag moment – rear-right [N]"),
        ],
        vars_to_check=[(["F_z_aero_fl", "F_z_aero_f"], 1)],
        name="Aero Load incl. Drag Moment – per wheel (symmetric)"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_aero_full_moments_fl{{i}}", "F_z_aero_full_moments_f{{i}} / 2"),
            ("F_z_aero_full_moments_fr{{i}}", "F_z_aero_full_moments_f{{i}} / 2"),
            ("F_z_aero_full_moments_rl{{i}}", "F_z_aero_full_moments_r{{i}} / 2"),
            ("F_z_aero_full_moments_rr{{i}}", "F_z_aero_full_moments_r{{i}} / 2"),
        ],
        relevant_vars=[
            ("F_z_aero_full_moments_fl", "Aero load incl. both moments – front-left [N]"),
            ("F_z_aero_full_moments_fr", "Aero load incl. both moments – front-right [N]"),
            ("F_z_aero_full_moments_rl", "Aero load incl. both moments – rear-left [N]"),
            ("F_z_aero_full_moments_rr", "Aero load incl. both moments – rear-right [N]"),
        ],
        vars_to_check=[(["F_z_aero_full_moments_fl", "F_z_aero_full_moments_f"], 1)],
        name="Aero Load incl. Full Moments – per wheel (symmetric)"
    )


    # ────────────────────────────────────────────────
    # Total load = static + aero incl. drag moment (recommended path)
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_total_f{{i}}", "F_z_static_f{{i}} + F_z_aero_f{{i}}"),
            ("F_z_total_r{{i}}", "F_z_static_r{{i}} + F_z_aero_r{{i}}"),
            ("F_z_total{{i}}",   "F_z_total_f{{i}} + F_z_total_r{{i}}"),
        ],
        relevant_vars=[
            ("F_z_total",   "Total normal load incl. static + aero (drag moment only) [N]"),
            ("F_z_total_f", "Front axle total load incl. aero (drag moment) [N]"),
            ("F_z_total_r", "Rear axle total load incl. aero (drag moment) [N]"),
        ],
        vars_to_check=[(["F_z_total_f", "F_z_static_f"], 1)],
        name="Total Normal Load incl. Static + Aero (recommended)"
    )


    # ────────────────────────────────────────────────
    # Lateral load transfer
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("W_lat{{i}}", "(h_CG{{i}} * m{{i}} * a_lat{{i}}) / (t_track{{i}} * g{{i}})")],
        relevant_vars=[
            ("W_lat",    "Total lateral load transfer [N]"),
            ("h_CG",     "CG height [m]"),
            ("m",        "Mass [kg]"),
            ("a_lat",    "Lateral acceleration [m/s²]"),
            ("t_track",  "Track width [m]"),
            ("g",        "Gravitational acceleration [m/s²]")
        ],
        vars_to_check=[(["W_lat", "h_CG", "m", "a_lat", "t_track"], 5)],
        name="Total Lateral Load Transfer"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_lat_f{{i}}",  "W_lat{{i}} * frac_lat_front{{i}}"),
            ("Delta_Fz_lat_r{{i}}",  "W_lat{{i}} * (1 - frac_lat_front{{i}})"),
            ("Delta_Fz_lat_fl{{i}}", "Delta_Fz_lat_f{{i}} / 2"),
            ("Delta_Fz_lat_fr{{i}}", "-Delta_Fz_lat_f{{i}} / 2"),
            ("Delta_Fz_lat_rl{{i}}", "Delta_Fz_lat_r{{i}} / 2"),
            ("Delta_Fz_lat_rr{{i}}", "-Delta_Fz_lat_r{{i}} / 2"),
        ],
        relevant_vars=[
            ("Delta_Fz_lat_f",  "Lateral load added to front axle [N]"),
            ("Delta_Fz_lat_r",  "Lateral load added to rear axle [N]"),
            ("Delta_Fz_lat_fl", "Lateral load change front-left [N]"),
            ("Delta_Fz_lat_fr", "Lateral load change front-right [N]"),
            ("Delta_Fz_lat_rl", "Lateral load change rear-left [N]"),
            ("Delta_Fz_lat_rr", "Lateral load change rear-right [N]"),
            ("frac_lat_front",  "Fraction of lateral transfer to front [0–1]"),
        ],
        vars_to_check=[(["Delta_Fz_lat_fl", "W_lat"], 1)],
        name="Lateral Load Transfer Distribution"
    )


    # ────────────────────────────────────────────────
    # Final cornering wheel loads (incl. lateral transfer)
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("F_z_total_cornering_fl{{i}}",
            "F_z_total_fl{{i}} + Delta_Fz_lat_fl{{i}}"),
            ("F_z_total_cornering_fr{{i}}",
            "F_z_total_fr{{i}} + Delta_Fz_lat_fr{{i}}"),
            ("F_z_total_cornering_rl{{i}}",
            "F_z_total_rl{{i}} + Delta_Fz_lat_rl{{i}}"),
            ("F_z_total_cornering_rr{{i}}",
            "F_z_total_rr{{i}} + Delta_Fz_lat_rr{{i}}"),
        ],
        relevant_vars=[
            ("F_z_total_cornering_fl", "Final front-left wheel normal load – cornering [N]"),
            ("F_z_total_cornering_fr", "Final front-right wheel normal load – cornering [N]"),
            ("F_z_total_cornering_rl", "Final rear-left wheel normal load – cornering [N]"),
            ("F_z_total_cornering_rr", "Final rear-right wheel normal load – cornering [N]"),
        ],
        vars_to_check=[(["F_z_total_cornering_fl", "F_z_total_fl"], 1)],
        name="Final Cornering Wheel Loads"
    )


    # ────────────────────────────────────────────────
    # Tyre model basics – using generic F_z normal load
    # ────────────────────────────────────────────────

    # 1. Tyre vertical stiffness (linear approximation)
    aero_eq_manager.add_equation_template(
        equations_to_add=[("k_tire{{i}}", "F_z{{i}} / delta_z_tire{{i}}")],
        relevant_vars=[
            ("k_tire",       "Tyre vertical stiffness [N/m]"),
            ("F_z",          "Generic normal force on the tyre [N] – set to desired level (F_z_static_*, F_z_aero_*, F_z_total_*, etc.)"),
            ("delta_z_tire", "Tyre vertical deflection / compression [m]")
        ],
        vars_to_check=[(["k_tire", "F_z"], 1)],
        name="Tyre Vertical Stiffness (linear)"
    )


    # 2. Pacejka-like peak friction coefficient scaling with load
    #    (very simplified – real Pacejka has more terms)
    aero_eq_manager.add_equation_template(
        equations_to_add=[("mu_peak{{i}}", "mu_0{{i}} * (A_load{{i}} + B_load{{i}} * (F_z{{i}} / F_z_nom{{i}})) / (1 + C_load{{i}} * (F_z{{i}} / F_z_nom{{i}}))")],
        relevant_vars=[
            ("mu_peak",    "Peak friction coefficient at current load [-]"),
            ("mu_0",       "Reference peak friction at nominal load [-]"),
            ("F_z",        "Generic normal force on the tyre [N] – set to desired level (F_z_static_*, F_z_aero_*, F_z_total_*, etc.)"),
            ("F_z_nom",    "Nominal/reference load for mu_0 [N]"),
            ("A_load", "Load sensitivity coefficients (typical: ≈1)"),
            ("B_load", "Load sensitivity coefficients (typical: ≈0)"),
            ("C_load", "Load sensitivity coefficients (typical: ≈0.1–0.3)")
            
        ],
        vars_to_check=[(["mu_peak", "mu_0", "F_z", "F_z_nom"], 3)],
        name="Load-dependent Peak Friction (simplified Pacejka style)"
    )


    # 3. Maximum lateral force (simple linear + saturation approximation)
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_y_max{{i}}", "mu_peak{{i}} * F_z{{i}}")],
        relevant_vars=[
            ("F_y_max",    "Maximum achievable lateral force [N]"),
            ("mu_peak",    "Peak friction coefficient [-]"),
            ("F_z",        "Generic normal force on the tyre [N] – set to desired level (F_z_static_*, F_z_aero_*, F_z_total_*, etc.)")
        ],
        vars_to_check=[(["F_y_max", "mu_peak", "F_z"], 2)],
        name="Maximum Lateral Force"
    )


    # 4. Linear slip angle range – cornering stiffness
    aero_eq_manager.add_equation_template(
        equations_to_add=[("C_alpha{{i}}", "C_alpha_0{{i}} * (F_z{{i}} / F_z_nom{{i}})**exponent{{i}}")],
        relevant_vars=[
            ("C_alpha",     "Tyre cornering stiffness [N/rad]"),
            ("C_alpha_0",   "Reference cornering stiffness at nominal load [N/rad]"),
            ("F_z",         "Generic normal force on the tyre [N] – set to desired level (F_z_static_*, F_z_aero_*, F_z_total_*, etc.)"),
            ("F_z_nom",     "Nominal load [N]"),
            ("exponent",    "Load exponent (typical 0.6–0.9 for passenger / race tyres)")
        ],
        vars_to_check=[(["C_alpha", "C_alpha_0", "F_z", "F_z_nom"], 3)],
        name="Load-dependent Cornering Stiffness"
    )


    # 5. Lateral force in linear region (small slip angle)
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_y_lin{{i}}", "C_alpha{{i}} * alpha{{i}}")],
        relevant_vars=[
            ("F_y_lin",    "Lateral force in linear slip range [N]"),
            ("C_alpha",    "Cornering stiffness [N/rad]"),
            ("alpha",      "Slip angle [rad]")
        ],
        vars_to_check=[(["F_y_lin", "C_alpha", "alpha"], 2)],
        name="Lateral Force – linear range"
    )


    # 6. Very simple combined slip – normalised force demand
    aero_eq_manager.add_equation_template(
        equations_to_add=[("sigma_combined{{i}}", "sqrt( (F_x{{i}} / F_x_max{{i}})**2 + (F_y{{i}} / F_y_max{{i}})**2 )")],
        relevant_vars=[
            ("sigma_combined", "Combined normalised force demand [-] (≤1 = grip available)"),
            ("F_x",            "Longitudinal force demand [N]"),
            ("F_y",            "Lateral force demand [N]"),
            ("F_x_max",        "Maximum longitudinal force [N] (often ≈ mu_peak * F_z)"),
            ("F_y_max",        "Maximum lateral force [N]")
        ],
        vars_to_check=[(["sigma_combined", "F_x", "F_y", "F_x_max", "F_y_max"], 3)],
        name="Combined Slip – normalised force utilisation"
    )


    # 7. Friction ellipse / friction circle limit check (simple)
    aero_eq_manager.add_equation_template(
        equations_to_add=[("friction_utilisation{{i}}", "sqrt( (F_x{{i}} / (mu_peak{{i}} * F_z{{i}}))**2 + (F_y{{i}} / (mu_peak{{i}} * F_z{{i}}))**2 )")],
        relevant_vars=[
            ("friction_utilisation", "Friction circle utilisation [-] (≤1 = no slip limit exceeded)"),
            ("F_x",                  "Longitudinal force [N]"),
            ("F_y",                  "Lateral force [N]"),
            ("mu_peak",              "Peak friction coefficient [-]"),
            ("F_z",                  "Generic normal force on the tyre [N] – set to desired level (F_z_static_*, F_z_aero_*, F_z_total_*, etc.)")
        ],
        vars_to_check=[(["friction_utilisation", "F_x", "F_y", "mu_peak", "F_z"], 4)],
        name="Friction Circle Utilisation"
    )

    #TODO: add:
    #cd of sphere
    #free fall with drag
    #coastdown with drag
    #thrust required vs speed
    #power required vs speed

    #cd of sphere
    aero_eq_manager.add_equation_template(
        equations_to_add=[("C_d_sphere{{i}}", "24 / Re{{i}} + (4.0 / (0.4 + (Re{{i}} / 282000)**0.65)) + (0.4 / (1 + (Re{{i}} / 2.8e6)**0.65))")],
        relevant_vars=[
            ("C_d_sphere", "Drag coefficient of a sphere [-]"),
            ("Re",         "Reynolds number [-]")
        ],
        vars_to_check=[(["C_d_sphere", "Re"], 1)],
        name="Drag Coefficient of a Sphere"
    )
    #free fall with drag
    aero_eq_manager.add_equation_template(
        equations_to_add=[("v_terminal{{i}}", "sqrt((2 * m{{i}} * g{{i}}) / (rho{{i}} * C_d_sphere{{i}} * A{{i}}))")],
        relevant_vars=[
            ("v_terminal", "Terminal velocity during free fall with drag [m/s]"),
            ("m",          "Mass [kg]"),
            ("g",          "Gravitational acceleration [m/s²]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d_sphere","Drag coefficient of a sphere [-]"),
            ("A",          "Cross-sectional area [m²]")
        ],
        vars_to_check=[(["v_terminal", "m", "g", "rho", "C_d_sphere", "A"], 5)],
        name="Terminal Velocity during Free Fall with Drag"
    )
    #free fall with drag and height
    aero_eq_manager.add_equation_template(
        equations_to_add=[("t_fall{{i}}", " (m{{i}} / (rho{{i}} * C_d_sphere{{i}} * A{{i}})) * log( (v_terminal{{i}}) / (v_terminal{{i}} - sqrt(2 * g{{i}} * h{{i}})) )")],
        relevant_vars=[
            ("t_fall",     "Time to fall from height h with drag [s]"),
            ("m",          "Mass [kg]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d_sphere","Drag coefficient of a sphere [-]"),
            ("A",          "Cross-sectional area [m²]"),
            ("v_terminal","Terminal velocity during free fall with drag [m/s]"),
            ("g",          "Gravitational acceleration [m/s²]"),
            ("h",          "Height of fall [m]")
        ],
        vars_to_check=[(["t_fall", "m", "rho", "C_d_sphere", "A", "v_terminal", "g", "h"], 7)],
        name="Time to Fall from Height with Drag"
    )
    #coastdown with drag
    aero_eq_manager.add_equation_template(
        equations_to_add=[("t_coast{{i}}", "(m{{i}} / (0.5 * rho{{i}} * C_d{{i}} * A{{i}} * v{{i}})) * log(v{{i}} / v_final{{i}})")],
        relevant_vars=[
            ("t_coast",    "Time to coast down from initial speed v to final speed v_final with drag [s]"),
            ("m",          "Mass [kg]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d",        "Drag coefficient [-]"),
            ("A",          "Cross-sectional area [m²]"),
            ("v",          "Initial velocity [m/s]"),
            ("v_final",    "Final velocity [m/s]")
        ],
        vars_to_check=[(["t_coast", "m", "rho", "C_d_sphere", "A", "v", "v_final"], 6)],
        name="Time to Coast Down with Drag"
    )
    #coastdown distance with drag
    aero_eq_manager.add_equation_template(
        equations_to_add=[("s_coast{{i}}", "(m{{i}} / (0.5 * rho{{i}} * C_d{{i}} * A{{i}})) * (v{{i}} - v_final{{i}})")],
        relevant_vars=[
            ("s_coast",    "Distance to coast down from initial speed v to final speed v_final with drag [m]"),
            ("m",          "Mass [kg]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d",        "Drag coefficient [-]"),
            ("A",          "Cross-sectional area [m²]"),
            ("v",          "Initial velocity [m/s]"),
            ("v_final",    "Final velocity [m/s]")
        ],
        vars_to_check=[(["s_coast", "m", "rho", "C_d_sphere", "A", "v", "v_final"], 6)],
        name="Distance to Coast Down with Drag"
    )
    #thrust required vs speed
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_thrust{{i}}", "F_d{{i}} + m{{i}} * a{{i}}")],
        relevant_vars=[
            ("F_thrust",   "Required thrust to maintain acceleration a at speed v [N]"),
            ("F_d",        "Aerodynamic drag [N]"),
            ("m",          "Mass [kg]"),
            ("a",          "Desired acceleration [m/s²]")
        ],
        vars_to_check=[(["F_thrust", "F_d", "m", "a"], 3)],
        name="Required Thrust vs Speed"
    )
    #F_d
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_d{{i}}", "1/2 * rho{{i}} * C_d{{i}} * A{{i}} * v{{i}}**2")],
        relevant_vars=[
            ("F_d",        "Aerodynamic drag [N]"),
            ("rho",        "Air density [kg/m³]"),
            ("C_d",        "Drag coefficient [-]"),
            ("A",          "Cross-sectional area [m²]"),
            ("v",          "Velocity [m/s]")
        ],
        vars_to_check=[(["F_d", "rho", "C_d", "A", "v"], 4)],
        name="Aerodynamic Drag"
    )
    #power required vs speed
    aero_eq_manager.add_equation_template(
        equations_to_add=[("P_required{{i}}", "F_thrust{{i}} * v{{i}}")],
        relevant_vars=[
            ("P_required", "Power required to maintain acceleration a at speed v [W]"),
            ("F_thrust",   "Required thrust to maintain acceleration a at speed v [N]"),
            ("v",          "Velocity [m/s]")
        ],
        vars_to_check=[(["P_required", "F_thrust", "v"], 2)],
        name="Power Required vs Speed"
    )
    #F_thrust, definition and name is wrong. TODO: fix this and equation above.
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_thrust{{i}}", "F_d{{i}} + m{{i}} * a{{i}}")],
        relevant_vars=[
            ("F_thrust",   "Required thrust to maintain acceleration a at speed v [N]"),
            ("F_d",        "Aerodynamic drag [N]"),
            ("m",          "Mass [kg]"),
            ("a",          "Desired acceleration [m/s²]")
        ],
        vars_to_check=[(["F_thrust", "F_d", "m", "a"], 3)],
        name="Thrust Required vs Speed"
    )

    #bernoulli
    aero_eq_manager.add_equation_template(
        equations_to_add=[("p{{i}}+1/2*rho{{i}}*v{{i}}^2+rho{{i}}*g*zi=p{{j}}+1/2*rho{{j}}*v{{j}}^2+rho{{j}}*g*h{{j}}")],
        relevant_vars=[
            ("p", "Pressure at point [Pa]"),
            ("rho", "Fluid density [kg/m³]"),
            ("v", "Velocity at point [m/s]"),
            ("g", "Gravitational acceleration [m/s²]"),
            ("h", "Height at point [m]")
        ],
        vars_to_check=[(["P", "rho", "v", "g", "h"], 3)],
        name="Bernoulli's Equation"
    )

    #barometric formula
    aero_eq_manager.add_equation_template(
        equations_to_add=[("p{{i}}=p_0{{i}}*exp(-rho_0{{i}}*g{{i}}*h{{i}}/(p_0{{i}}))")],
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
        equations_to_add=[("rho{{i}}=rho_0{{i}}*exp(-rho_0{{i}}*g{{i}}*h{{i}}/(p_0{{i}}))")],
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

    #longitudinal weight shift
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_long_f{{i}}",  "m{{i}} * a_long{{i}} * b_CG{{i}} / l_wb{{i}}"),
            ("Delta_Fz_long_r{{i}}",  "-m{{i}} * a_long{{i}} * a_CG{{i}} / l_wb{{i}}"),
        ],
        relevant_vars=[
            ("Delta_Fz_long_f", "Front axle load change from longitudinal acceleration [N]"),
            ("Delta_Fz_long_r", "Rear axle load change from longitudinal acceleration [N]"),
            ("m",               "Vehicle mass [kg]"),
            ("a_long",         "Longitudinal acceleration [m/s²]"),
            ("a_CG",           "Longitudinal distance from CG to front axle [m]"),
            ("b_CG",           "Longitudinal distance from CG to rear axle [m]"),
            ("l_wb",           "Wheelbase [m]")
        ],
        vars_to_check=[(["Delta_Fz_long_f", "m", "a_long", "b_CG", "l_wb"], 4)],
        name="Axle Load Shift – Longitudinal Acceleration"
    )

    #lateral weight shift 
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_lat_f{{i}}",  "m{{i}} * a_lat{{i}} * h_CG{{i}} / t_track{{i}}"),
            ("Delta_Fz_lat_r{{i}}",  "-m{{i}} * a_lat{{i}} * h_CG{{i}} / t_track{{i}}"),
        ],
        relevant_vars=[
            ("Delta_Fz_lat_f", "Left side load change from lateral acceleration [N]"),
            ("Delta_Fz_lat_r", "Right side load change from lateral acceleration [N]"),
            ("m",               "Vehicle mass [kg]"),
            ("a_lat",          "Lateral acceleration [m/s²]"),
            ("h_CG",           "CG height [m]"),
            ("t_track",        "Track width [m]")
        ],
        vars_to_check=[(["Delta_Fz_lat_f", "m", "a_lat", "h_CG", "t_track"], 4)],
        name="Side Load Shift – Lateral Acceleration"
    )

    #F_lat
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_lat{{i}}", "m{{i}} * a_lat{{i}}")],
        relevant_vars=[
            ("F_lat",    "Lateral force [N]"),
            ("m",        "Mass [kg]"),
            ("a_lat",    "Lateral acceleration [m/s²]")
        ],
        vars_to_check=[(["F_lat", "a_lat"], 1)],
        name="Lateral Force from Lateral Acceleration"
    )
    #F_long
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_long{{i}}", "m{{i}} * a_long{{i}}")],
        relevant_vars=[
            ("F_long",   "Longitudinal force [N]"),
            ("m",        "Mass [kg]"),
            ("a_long",   "Longitudinal acceleration [m/s²]")
        ],
        vars_to_check=[(["F_long", "a_long"], 1)],
        name="Longitudinal Force from Longitudinal Acceleration"
    )
    #friction
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_max{{i}}", "mu{{i}} * F_z{{i}}")],
        relevant_vars=[
            ("F_fmax",   "Maximum frictional force [N]"),
            ("mu",         "Coefficient of friction [-]"),
            ("F_z",        "Normal force [N]")
        ],
        vars_to_check=[(["F_fmax", "mu", "F_z"], 2)],
        name="Maximum Frictional Force from Normal Load"
    )
    #max acceleration
    aero_eq_manager.add_equation_template(
        equations_to_add=[("a_max{{i}}", "mu{{i}} * F_z{{i}} / m{{i}}")],
        relevant_vars=[
            ("a_max",    "Maximum acceleration [m/s²]"),
            ("mu",         "Coefficient of friction [-]"),
            ("F_z",        "Normal force [N]"),
            ("m",        "Mass [kg]")
        ],
        vars_to_check=[(["a_max", "mu", "F_z"], 2)],
        name="Maximum Acceleration from Friction"
    )
    #centripetal force
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_centripetal{{i}}", "m{{i}} * v{{i}}**2 / r_curve{{i}}")],
        relevant_vars=[
            ("F_centripetal", "Centripetal force [N]"),
            ("m",             "Mass [kg]"),
            ("v",             "Velocity [m/s]"),
            ("r_curve",       "Curve radius [m]")
        ],
        vars_to_check=[(["F_centripetal", "m", "v", "r_curve"], 3)],
        name="Centripetal Force"
    )
    #banked curve normal force
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_z_banked{{i}}", "m{{i}} * g{{i}} / (cos(theta_banked{{i}}) - (v{{i}}**2 / (r_curve{{i}} * g{{i}})) * sin(theta_banked{{i}}))")],
        relevant_vars=[
            ("F_z_banked",   "Normal force on banked curve [N]"),
            ("m",            "Mass [kg]"),
            ("g",            "Gravitational acceleration [m/s²]"),
            ("theta_banked", "Bank angle [rad]"),
            ("v",            "Velocity [m/s]"),
            ("r_curve",      "Curve radius [m]")
        ],
        vars_to_check=[(["F_z_banked", "m", "g", "theta_banked", "v", "r_curve"], 5)],
        name="Normal Force on Banked Curve"
    )
    #banked turn F_centripetal with friction mu
    aero_eq_manager.add_equation_template(
        equations_to_add=[("F_centripetal_banked{{i}}", "mu{{i}} * F_z_banked{{i}} + m{{i}} * g{{i}} * sin(theta_banked{{i}})")],
        relevant_vars=[
            ("F_centripetal_banked", "Centripetal force on banked curve [N]"),
            ("mu",                    "Coefficient of friction [-]"),
            ("F_z_banked",           "Normal force on banked curve [N]"),
            ("m",                    "Mass [kg]"),
            ("g",                    "Gravitational acceleration [m/s²]"),
            ("theta_banked",         "Bank angle [rad]")
        ],
        vars_to_check=[(["F_centripetal_banked", "mu", "F_z_banked", "m", "g", "theta_banked"], 5)],
        name="Centripetal Force on Banked Curve with Friction"
    )

    #weight shift banked, elastic, geometric, and total
    aero_eq_manager.add_equation_template(
        equations_to_add=[
            ("Delta_Fz_banked{{i}}",  "m{{i}} * g{{i}} * tan(theta_banked{{i}}) * h_CG{{i}} / t_track{{i}}"),
            ("Delta_Fz_elastic{{i}}",  "F_z_total{{i}} * h_CG{{i}} / (k_susp{{i}} * t_track{{i}}) * a_lat{{i}}"),
            ("Delta_Fz_geometric{{i}}",  "m{{i}} * a_lat{{i}} * h_CG{{i}} / t_track{{i}}"),
            ("Delta_Fz_total{{i}}",  "Delta_Fz_banked{{i}} + Delta_Fz_elastic{{i}} + Delta_Fz_geometric{{i}}"),
        ],
        relevant_vars=[
            ("Delta_Fz_banked",   "Lateral load change from banked curve [N]"),
            ("Delta_Fz_elastic",  "Lateral load change from suspension elasticity [N]"),
            ("Delta_Fz_geometric","Lateral load change from geometric weight shift [N]"),
            ("Delta_Fz_total",    "Total lateral load change [N]"),
            ("m",                 "Vehicle mass [kg]"),
            ("g",                 "Gravitational acceleration [m/s²]"),
            ("theta_banked",      "Bank angle [rad]"),
            ("h_CG",              "CG height [m]"),
            ("t_track",           "Track width [m]"),
            ("F_z_total",         "Total normal load [N]"),
            ("k_susp",            "Suspension roll stiffness [N/m]"),
            ("a_lat",             "Lateral acceleration [m/s²]")
        ],
        vars_to_check=[(["Delta_Fz_total", "Delta_Fz_banked", "Delta_Fz_elastic", "Delta_Fz_geometric"], 3)],
        name="Lateral Load Shift – Banked Curve, Elastic, Geometric, and Total"
    )

    #energy kinetic
    aero_eq_manager.add_equation_template(
        equations_to_add=[("E_kinetic{{i}}", "0.5 * m{{i}} * v{{i}}**2")],
        relevant_vars=[
            ("E_kinetic", "Kinetic energy [J]"),
            ("m",         "Mass [kg]"),
            ("v",         "Velocity [m/s]")
        ],
        vars_to_check=[(["E_kinetic", "m", "v"], 2)],
        name="Kinetic Energy"
    )
    #energy potential
    aero_eq_manager.add_equation_template(
        equations_to_add=[("E_potential{{i}}", "m{{i}} * g{{i}} * h{{i}}")],
        relevant_vars=[
            ("E_potential", "Potential energy [J]"),
            ("m",           "Mass [kg]"),
            ("g",           "Gravitational acceleration [m/s²]"),
            ("h",           "Height [m]")
        ],
        vars_to_check=[(["E_potential", "m", "g", "h"], 3)],
        name="Potential Energy"
    )
    #energy total
    aero_eq_manager.add_equation_template(
        equations_to_add=[("E_total{{i}}", "E_kinetic{{i}} + E_potential{{i}}")],
        relevant_vars=[
            ("E_total", "Total mechanical energy [J]"),
            ("E_kinetic", "Kinetic energy [J]"),
            ("E_potential", "Potential energy [J]")
        ],
        vars_to_check=[(["E_total", "E_kinetic", "E_potential"], 2)],
        name="Total Mechanical Energy"
    )
    