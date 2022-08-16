/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.
 * Copyright (c) 2022 Virtonomy GmbH                                  *​
 * Copyright (c) 2022 John Benjamin												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file 	pkj_biventricular.h
 * @brief 	Electro-contraction of biventricle heart model.
 * @author 	Chi ZHANG and Xiangyu HU, John Benjamin
 * @version  0.3
 * @version 0.2.1
 * 			Chi Zhang
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 *@version 0.3
 *			Here, the coupling with Purkinje network will be condcuted.
 */
#pragma once
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
/** Set the file path to the stl file. */
std::string full_path_to_stl = "./input/leftventricle.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
Vec3d translation = Vec3d(0);
/** Domain bounds of the system. */
shared_ptr<TriangleMeshShapeSTL> triangle_mesh_heart_shape = make_shared<TriangleMeshShapeSTL>(full_path_to_stl, translation, length_scale);
BoundingBox system_domain_bounds = triangle_mesh_heart_shape->getBounds();
/** Paremeters and physical properties. */
Vec3d domain_lower_bound = system_domain_bounds.first;
Vec3d domain_upper_bound = system_domain_bounds.second;
Real dp_0 = 2.0;
/** Landmark coordinates */
Vec3d apex_coordinates = {111.65, -230.926, -195.884};
Vec3d av_node_coordinates = {28.2935, -215.659, -137.279};
Vec3d base_coordinates = {15.0342, -175.735, -151.837};
Vec3d center_line = base_coordinates - apex_coordinates;

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 150 * stress_scale;
Real a0[4] = {496.0 * stress_scale, 15196.0 * stress_scale, 3283.0 * stress_scale, 662.0 * stress_scale};
Real b0[4] = {7.209, 20.417, 11.176, 9.466};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
Real acceleration_factor = 27.5; /** Acceleration factor for fast dissuion on purkinje network. */
Real diffusion_coff = 0.8;
Real bias_coff = 0.0;
/** Electrophysiology parameters. */
StdVec<std::string> species_name_list{"Phi"};
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
/** Purkinje Network. */
Vec3d starting_point_RV = av_node_coordinates;
Vec3d second_point_RV = {30.0448, -216.7685, -137.8286};
Vec3d starting_point_LV = {28.6329, -211.3944, -141.0645};
Vec3d direction_vector = apex_coordinates - starting_point_LV;
Vec3d second_point_LV = starting_point_LV + 1.5 * (direction_vector / (direction_vector.norm() + TinyReal));
//----------------------------------------------------------------------
//	Define heart shape
//----------------------------------------------------------------------
class Heart : public ComplexShape
{
public:
	explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(full_path_to_stl, translation, length_scale);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class FiberDirectionDiffusion : public DiffusionReaction<LocallyOrthotropicMuscle>
{
public:
	FiberDirectionDiffusion()
		: DiffusionReaction<LocallyOrthotropicMuscle>(
			  species_name_list, rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0)
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	};
};
/** Set diffusion relaxation. */
class DiffusionRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  RelaxationOfAllDiffussionSpeciesInner<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>>
{
public:
	explicit DiffusionRelaxation(BodyRelationInner &body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
	virtual ~DiffusionRelaxation(){};
};
/** Imposing diffusion boundary condition */
class DiffusionBCs
	: public ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);

		Real angle = dot(face_norm, center_norm);
		if (angle >= 0.0)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		else
		{
			if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
				species_n_[phi_][index_i] = 0.0;
		}
	};

public:
	DiffusionBCs(SolidBody &body, BodySurface &body_part)
		: ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>(body, body_part)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
	virtual ~DiffusionBCs(){};
};
/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberandSheetDirections
	: public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	Real beta_epi_, beta_endo_;
	/** We define the centerline vector, which is parallel to the ventricular centerline and pointing  apex-to-base.*/
	Vecd center_line_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		/**
		 * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
		 * 		Present  doi.org/10.1016/j.cma.2016.05.031
		 */
		/** Probe the face norm from Levelset field. */
		Vecd dist_2_face = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
		if (dot(face_norm, center_norm) <= 0.0)
		{
			face_norm = -face_norm;
		}
		/** Compute the centerline's projection on the plane orthogonal to face norm. */
		Vecd center_line_proj = center_line_ - dot(center_line_, face_norm) * face_norm;
		Vecd center_line_proj_norm = center_line_proj / (center_line_proj.norm() + 1.0e-15);
		/** Define initial "flat" fiber field */
		Vecd circumferential_direction = SimTK::cross(center_line_, face_norm);
		Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
		/** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
		Real beta = (beta_epi_ - beta_endo_) * species_n_[phi_][index_i] + beta_endo_;
		/** Compute the rotation matrix through Rodrigues rotation formulation. */
		Vecd f_0 = cos(beta) * cd_norm + sin(beta) * SimTK::cross(face_norm, cd_norm) +
				   dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

		if (pos_n_[index_i][2] < 2.0 * body_->sph_adaptation_->ReferenceSpacing())
		{
			material_->local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
			material_->local_s0_[index_i] = face_norm;
		}
		else
		{
			material_->local_f0_[index_i] = Vecd(0);
			material_->local_s0_[index_i] = Vecd(0);
		}
	};

public:
	explicit ComputeFiberandSheetDirections(SolidBody &body)
		: DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
		center_line_ = center_line;
		beta_epi_ = -(70.0 / 180.0) * M_PI;
		beta_endo_ = (80.0 / 180.0) * M_PI;
	};
	virtual ~ComputeFiberandSheetDirections(){};
};
//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
public:
	MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
	{
		Real l = domain_upper_bound[0] - domain_lower_bound[0];
		Real w = domain_upper_bound[1] - domain_lower_bound[1];
		Real h = domain_upper_bound[2];
		halfsize_ = Vec3d(1.0 * dp_0, 0.5 * l, 0.5 * w);
		resolution_ = 20;
		translation_ = Vec3d(33.1861, -194.358, -152.545);
	}
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentToMmyocardium
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;
	Vec3d av_node_coordinates_;
	Real tolerance_length_ = 3.5;
	Real tolerance_width_ = 2.5;
	Real tolerance_depth_ = 0.5;

	void Update(size_t index_i, Real dt) override
	{
		if (std::abs(av_node_coordinates_[0] * length_scale - pos_n_[index_i][0]) <= tolerance_length_)
		{
			if (std::abs(av_node_coordinates_[1] * length_scale - pos_n_[index_i][1]) <= tolerance_width_)
			{
				if (std::abs(av_node_coordinates_[2] * length_scale - pos_n_[index_i][2]) <= tolerance_depth_)
				{
					species_n_[voltage_][index_i] = 0.92;
				}
			}
		}
	};

public:
	explicit ApplyStimulusCurrentToMmyocardium(SolidBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
		av_node_coordinates_ = av_node_coordinates;
	};
};
// Observer particle generator.
class HeartObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit HeartObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** position and volume. */
		positions_.push_back(Vecd(33.3116 * length_scale, -270.675 * length_scale, -175.672 * length_scale));
		positions_.push_back(Vecd(80.7354 * length_scale, -209.356 * length_scale, -134.835 * length_scale));
		positions_.push_back(Vecd(63.1737 * length_scale, -268.187 * length_scale, -191.3330 * length_scale));
		positions_.push_back(Vecd(100.818 * length_scale, -222.94 * length_scale, -152.37 * length_scale));
		positions_.push_back(Vecd(113.01 * length_scale, -235.75 * length_scale, 27.74 * length_scale));
	}
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentToPKJ
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		if (index_i <= 10)
		{
			species_n_[voltage_][index_i] = 1.0;
		}
	};

public:
	explicit ApplyStimulusCurrentToPKJ(RealBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};

/**
 * Derived network particle generator.
 */
class NetworkGeneratorWithExtraCheck : public ParticleGeneratorNetwork
{
protected:
	bool extraCheck(const Vecd &new_position) override
	{
		bool no_generation = false;
		if (new_position[2] > 0)
			no_generation = true;
		return no_generation;
	};

public:
	NetworkGeneratorWithExtraCheck(SPHBody &sph_body, Vecd starting_pnt, Vecd second_pnt, int iterator, Real grad_factor, Real angle, Real segments_in_branch)
		: ParticleGeneratorNetwork(sph_body, starting_pnt, second_pnt, iterator, grad_factor)
		{
			angle_ = angle;
			segments_in_branch_ = segments_in_branch;
		};
};