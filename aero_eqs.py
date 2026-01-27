def AddEquationsAero(aero_eq_manager):
    """
    Adds the standard set of equations for aerodynamic loads and vehicle dynamics
    to the provided Equation_Manager instance.
    """
    # ────────────────────────────────────────────────
    # Fluid dynamics & basic relations
    # ────────────────────────────────────────────────

    aero_eq_manager.add_equation_template(
        equations_to_add=[("Reii", "Lii * rhoii * vii / etaii")],
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
        equations_to_add=[("nuii", "etaii / rhoii")],
        relevant_vars=[
            ("nu",  "Kinematic viscosity [m²/s]"),
            ("eta", "Dynamic viscosity [Pa·s]"),
            ("rho", "Air density [kg/m³]")
        ],
        vars_to_check=[(["nu", "eta"], 1)],
        name="Kinematic Viscosity"
    )

    aero_eq_manager.add_equation_template(
        equations_to_add=[("vii / vjj", "Ajj / Aii")],
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
        equations_to_add=[("F_z_staticii", "mii * gii")],
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
            ("F_z_static_fii", "mii * gii * (l_wyii - x_CGii) / l_wyii"),
            ("F_z_static_rii", "mii * gii * x_CGii / l_wyii"),
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
            ("F_z_static_flii", "F_z_static_fii / 2"),
            ("F_z_static_frii", "F_z_static_fii / 2"),
            ("F_z_static_rlii", "F_z_static_rii / 2"),
            ("F_z_static_rrii", "F_z_static_rii / 2"),
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
        equations_to_add=[("F_z_dfii", "1/2 * rhoii * C_Lii * Aii * vii**2")],
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
            ("F_z_df_fii", "ABii * F_z_dfii"),
            ("F_z_df_rii", "(1 - ABii) * F_z_dfii"),
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
            ("F_z_df_flii", "F_z_df_fii / 2"),
            ("F_z_df_frii", "F_z_df_fii / 2"),
            ("F_z_df_rlii", "F_z_df_rii / 2"),
            ("F_z_df_rrii", "F_z_df_rii / 2"),
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
        equations_to_add=[("M_pitch_dfii", "F_z_dfii * (x_CGii - x_COPii)")],
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
        equations_to_add=[("M_pitch_dragii", "F_dii * (z_CGii - z_dragii)")],
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
            ("Delta_Fz_df_fii",  "M_pitch_dfii / l_wbii"),
            ("Delta_Fz_df_rii",  "-M_pitch_dfii / l_wbii"),
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
            ("Delta_Fz_drag_fii",  "M_pitch_dragii / l_wbii"),
            ("Delta_Fz_drag_rii",  "-M_pitch_dragii / l_wbii"),
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
            ("F_z_aero_fii", "F_z_df_fii + Delta_Fz_drag_fii"),
            ("F_z_aero_rii", "F_z_df_rii + Delta_Fz_drag_rii"),
            ("F_z_aeroii",   "F_z_aero_fii + F_z_aero_rii"),
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
            ("F_z_aero_full_moments_fii", "F_z_df_fii + Delta_Fz_df_fii + Delta_Fz_drag_fii"),
            ("F_z_aero_full_moments_rii", "F_z_df_rii + Delta_Fz_df_rii + Delta_Fz_drag_rii"),
            ("F_z_aero_full_momentsii",   "F_z_aero_full_moments_fii + F_z_aero_full_moments_rii"),
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
            ("F_z_aero_flii", "F_z_aero_fii / 2"),
            ("F_z_aero_frii", "F_z_aero_fii / 2"),
            ("F_z_aero_rlii", "F_z_aero_rii / 2"),
            ("F_z_aero_rrii", "F_z_aero_rii / 2"),
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
            ("F_z_aero_full_moments_flii", "F_z_aero_full_moments_fii / 2"),
            ("F_z_aero_full_moments_frii", "F_z_aero_full_moments_fii / 2"),
            ("F_z_aero_full_moments_rlii", "F_z_aero_full_moments_rii / 2"),
            ("F_z_aero_full_moments_rrii", "F_z_aero_full_moments_rii / 2"),
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
            ("F_z_total_fii", "F_z_static_fii + F_z_aero_fii"),
            ("F_z_total_rii", "F_z_static_rii + F_z_aero_rii"),
            ("F_z_totalii",   "F_z_total_fii + F_z_total_rii"),
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
        equations_to_add=[("W_latii", "(h_CGii * mii * a_latii) / (t_trackii * gii)")],
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
            ("Delta_Fz_lat_fii",  "W_latii * frac_lat_frontii"),
            ("Delta_Fz_lat_rii",  "W_latii * (1 - frac_lat_frontii)"),
            ("Delta_Fz_lat_flii", "Delta_Fz_lat_fii / 2"),
            ("Delta_Fz_lat_frii", "-Delta_Fz_lat_fii / 2"),
            ("Delta_Fz_lat_rlii", "Delta_Fz_lat_rii / 2"),
            ("Delta_Fz_lat_rrii", "-Delta_Fz_lat_rii / 2"),
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
            ("F_z_total_cornering_flii",
            "F_z_total_flii + Delta_Fz_lat_flii"),
            ("F_z_total_cornering_frii",
            "F_z_total_frii + Delta_Fz_lat_frii"),
            ("F_z_total_cornering_rlii",
            "F_z_total_rlii + Delta_Fz_lat_rlii"),
            ("F_z_total_cornering_rrii",
            "F_z_total_rrii + Delta_Fz_lat_rrii"),
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
        equations_to_add=[("k_tireii", "F_zii / delta_z_tireii")],
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
        equations_to_add=[("mu_peakii", "mu_0ii * (A_loadii + B_loadii * (F_zii / F_z_nomii)) / (1 + C_loadii * (F_zii / F_z_nomii))")],
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
        equations_to_add=[("F_y_maxii", "mu_peakii * F_zii")],
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
        equations_to_add=[("C_alphaii", "C_alpha_0ii * (F_zii / F_z_nomii)**exponentii")],
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
        equations_to_add=[("F_y_linii", "C_alphaii * alphaii")],
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
        equations_to_add=[("sigma_combinedii", "sqrt( (F_xii / F_x_maxii)**2 + (F_yii / F_y_maxii)**2 )")],
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
        equations_to_add=[("friction_utilisationii", "sqrt( (F_xii / (mu_peakii * F_zii))**2 + (F_yii / (mu_peakii * F_zii))**2 )")],
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
        equations_to_add=[("C_d_sphereii", "24 / Reii + (4.0 / (0.4 + (Reii / 282000)**0.65)) + (0.4 / (1 + (Reii / 2.8e6)**0.65))")],
        relevant_vars=[
            ("C_d_sphere", "Drag coefficient of a sphere [-]"),
            ("Re",         "Reynolds number [-]")
        ],
        vars_to_check=[(["C_d_sphere", "Re"], 1)],
        name="Drag Coefficient of a Sphere"
    )
    #free fall with drag
    aero_eq_manager.add_equation_template(
        equations_to_add=[("v_terminalii", "sqrt((2 * mii * gii) / (rhoii * C_d_sphereii * Aii))")],
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
        equations_to_add=[("t_fallii", " (mii / (rhoii * C_d_sphereii * Aii)) * log( (v_terminalii) / (v_terminalii - sqrt(2 * gii * hii)) )")],
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
        equations_to_add=[("t_coastii", "(mii / (0.5 * rhoii * C_dii * Aii * vii)) * log(vii / v_finalii)")],
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
        equations_to_add=[("s_coastii", "(mii / (0.5 * rhoii * C_dii * Aii)) * (vii - v_finalii)")],
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
        equations_to_add=[("F_thrustii", "F_dii + mii * aii")],
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
        equations_to_add=[("F_dii", "1/2 * rhoii * C_dii * Aii * vii**2")],
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
        equations_to_add=[("P_requiredii", "F_thrustii * vii")],
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
        equations_to_add=[("F_thrustii", "F_dii + mii * aii")],
        relevant_vars=[
            ("F_thrust",   "Required thrust to maintain acceleration a at speed v [N]"),
            ("F_d",        "Aerodynamic drag [N]"),
            ("m",          "Mass [kg]"),
            ("a",          "Desired acceleration [m/s²]")
        ],
        vars_to_check=[(["F_thrust", "F_d", "m", "a"], 3)],
        name="Thrust Required vs Speed"
    )
