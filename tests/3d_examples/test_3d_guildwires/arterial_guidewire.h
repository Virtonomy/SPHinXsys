// /**
//  * @file 	taylor_bar.h
//  * @brief 	This is the case setup for plastic taylor bar.
//  * @author 	Xiaojing Tang, Dong Wu and Xiangyu Hu
//  */
// #pragma once

// #include "sphinxsys.h"
// using namespace SPH;
// // UNITS

// namespace acds // transfemoral_arterial_delivery_system
// {

// static constexpr double scale_stl = 1; // in mm currently

// static constexpr double mass_scaling_phase_1 = 1e5;
// static constexpr double mass_scaling_phase_2 = 1e5;

// // UNITS
// namespace unit_system
// {
// // Units: g, mm, ms, N, MPa

// // Unit System 1
// // MASS	LENGTH	TIME	FORCE	STRESS	ENERGY	DENSITY		YOUNG's		35MPH 56.33KMPH		GRAVITY
// // g	mm		ms		N		MPa		N-mm	7.83e-03	2.07e+05	15.65				9.806e-03
// // https://www.dynasupport.com/howtos/general/consistent-units

// constexpr double mass_unit = 1e3;          // kg to g
// constexpr double length_unit = 1e3;        // m to mm
// constexpr double time_unit = 1e3;          // s to ms
// constexpr double pressure_unit = 1e-6;     // Pa to MPa
// constexpr double density_unit = 1e-6;      // kg/m^3 to g/mm^3
// constexpr double density_unit_2 = 1e6;     // tonne/mm3 to g/mm3
// constexpr double acceleration_unit = 1e-3; // m/s^2 to mm/ms^2
// constexpr double length_tolerance = 1e-9;  // TODO: decide if this is appropriate

// } // namespace unit_system

// struct delivery_system_input_data
// {
//     struct guidewire
//     {
//         // parameters collected by mongodb

//         double poisson = 0.3;
//         double full_length = 250.00;        // mm
//         double diameter = 1.066;            // mm - original: 0.889 mm
//         double resolution = diameter / 3.0; // diameter / 5.0;

//         double rho = 695.48 * unit_system::density_unit * mass_scaling_phase_1; // g/mm3
//         double youngs_modulus = 82909.0;                                        // MPa --> 83 GPa
//     };

//     // create instances of all the devices and structures in the simulation
//     delivery_system_input_data::guidewire guidewire;
// };

// } // namespace acds