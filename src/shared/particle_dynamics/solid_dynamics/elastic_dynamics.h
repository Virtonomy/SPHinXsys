/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	elastic_dynamics.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef VIRTOSIM_ELASTIC_DYNAMICS_H_C0E2274F_EA9C_462F_B9FA_890BC5B12F2B
#define VIRTOSIM_ELASTIC_DYNAMICS_H_C0E2274F_EA9C_462F_B9FA_890BC5B12F2B

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "elastic_solid.h"
#include "general_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
namespace solid_dynamics
{
//----------------------------------------------------------------------
//		for elastic solid dynamics
//----------------------------------------------------------------------
typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;
typedef DataDelegateInner<ElasticSolidParticles> ElasticSolidDataInner;

/**
 * @class ElasticDynamicsInitialCondition
 * @brief  set initial condition for a solid body with different material
 * This is a abstract class to be override for case specific initial conditions.
 */
class ElasticDynamicsInitialCondition : public LocalDynamics, public ElasticSolidDataSimple
{
  public:
    explicit ElasticDynamicsInitialCondition(SPHBody &sph_body);
    virtual ~ElasticDynamicsInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
};

/**
 * @class UpdateElasticNormalDirection
 * @brief update particle normal directions for elastic solid
 */
class UpdateElasticNormalDirection : public LocalDynamics, public ElasticSolidDataSimple
{
  protected:
    StdLargeVec<Vecd> &n_, &n0_;
    StdLargeVec<Matd> &F_;

  public:
    explicit UpdateElasticNormalDirection(SPHBody &sph_body);
    virtual ~UpdateElasticNormalDirection(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class AcousticTimeStepSize
 * @brief Computing the acoustic time step size
 * computing time step size
 */
class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
                             public ElasticSolidDataSimple
{
  protected:
    Real CFL_;
    StdLargeVec<Vecd> &vel_, &acc_, &acc_prior_;
    Real smoothing_length_, c0_;

  public:
    explicit AcousticTimeStepSize(SPHBody &sph_body, Real CFL = 0.6);
    virtual ~AcousticTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0);
};

/**
 * @class DeformationGradientBySummation
 * @brief computing deformation gradient tensor by summation
 */
class DeformationGradientBySummation : public LocalDynamics, public ElasticSolidDataInner
{
  public:
    explicit DeformationGradientBySummation(BaseInnerRelation &inner_relation);
    virtual ~DeformationGradientBySummation(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd &pos_n_i = pos_[index_i];

        Matd deformation = Matd::Zero();
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];

            Vecd gradW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            deformation -= (pos_n_i - pos_[index_j]) * gradW_ijV_j.transpose();
        }

        F_[index_i] = deformation * B_[index_i];
    };

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Matd> &B_, &F_;
};

/**
 * @class BaseElasticIntegration
 * @brief base class for elastic relaxation
 */
class BaseElasticIntegration : public LocalDynamics, public ElasticSolidDataInner
{
  public:
    explicit BaseElasticIntegration(BaseInnerRelation &inner_relation);
    virtual ~BaseElasticIntegration(){};

  protected:
    StdLargeVec<Real> &rho_, &mass_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_;
    StdLargeVec<Matd> &B_, &F_, &dF_dt_;
};

/**
 * @class BaseIntegration1stHalf
 * @brief computing stress relaxation process by verlet time stepping
 * This is the first step
 */
class BaseIntegration1stHalf : public BaseElasticIntegration
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalf(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    ElasticSolid &elastic_solid_;
    Real rho0_, inv_rho0_;
    StdLargeVec<Vecd> &acc_prior_;
    Real smoothing_length_;
};

/**
 * @class Integration1stHalf
 * @brief computing stress relaxation process by verlet time stepping
 * This is the first step
 */
class Integration1stHalf : public BaseIntegration1stHalf
{
  public:
    explicit Integration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalf(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        const Vecd &vel_i = vel_[index_i];
        const Real inv_rho0 = inv_rho0_; // Assuming inv_rho0_ is a class member
        const double epsilon = std::numeric_limits<double>::epsilon();

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            const size_t index_j = inner_neighborhood.j_[n];
            const Vecd &vel_j = vel_[index_j];
            Vecd vel_jump = vel_i - vel_j;
            double vel_squared_magnitude = vel_jump.squaredNorm();

            const Vecd &e_ij = inner_neighborhood.e_ij_[n];
            Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

            if (vel_squared_magnitude > epsilon)
            {
                const Real r_ij = inner_neighborhood.r_ij_[n];
                const Real dim_r_ij_1 = Dimensions / r_ij;
                Vecd pos_jump = pos_[index_i] - pos_[index_j];
                Real strain_rate = dim_r_ij_1 * dim_r_ij_1 * pos_jump.dot(vel_jump);
                Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
                Matd numerical_stress_ij = 0.5 * (F_[index_i] + F_[index_j]) * elastic_solid_.PairNumericalDamping(strain_rate, smoothing_length_);

                acceleration += inv_rho0 * dW_ijV_j *
                                (stress_PK1_B_[index_i] + stress_PK1_B_[index_j] +
                                 numerical_dissipation_factor_ * weight * numerical_stress_ij) *
                                e_ij;
            }
            else
            {
                acceleration += inv_rho0 * dW_ijV_j *
                                (stress_PK1_B_[index_i] + stress_PK1_B_[index_j]) *
                                e_ij;
            }
        }

        acc_[index_i] = acceleration;
    };

  protected:
    StdLargeVec<Matd> stress_PK1_B_;
    Real numerical_dissipation_factor_;
    Real inv_W0_ = 1.0 / sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd);
};

/**
 * @class Integration1stHalfPK2
 * @brief Using PK2 stress constitute relation
 */
class Integration1stHalfPK2 : public Integration1stHalf
{
  public:
    explicit Integration1stHalfPK2(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalfPK2(){};
    void initialization(size_t index_i, Real dt = 0.0);
};

/** @class Integration1stHalfCauchy
 * @brief Using Cauchy stress constitute relation
 */
class Integration1stHalfCauchy : public Integration1stHalf
{
  public:
    explicit Integration1stHalfCauchy(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalfCauchy(){};
    void initialization(size_t index_i, Real dt = 0.0);
};

/**
 * @class Integration1stHalfKirchhoff
 * @brief Using Kirchhoff stress constitute relation
 */
class Integration1stHalfKirchhoff : public Integration1stHalf
{
  public:
    explicit Integration1stHalfKirchhoff(BaseInnerRelation &inner_relation);
    virtual ~Integration1stHalfKirchhoff(){};
    void initialization(size_t index_i, Real dt = 0.0);
};

/**
 * @class DecomposedIntegration1stHalf
 * @brief Decompose the stress into particle stress includes isotropic stress
 * and the stress due to non-homogeneous material properties.
 * The preliminary shear stress is introduced by particle pair to avoid
 * spurious stress and deformation.
 * Note that, for the shear stress term,
 * due to the mismatch of the divergence contribution between
 * the pair-wise second-order derivative Laplacian formulation
 * and particle-wise first-order gradient formulation,
 * a correction factor slight large than one is introduced.
 * Note that, if you see time step size goes unusually small,
 * it may be due to the determinate of deformation matrix become negative.
 * In this case, you may need decrease CFL number when computing time-step size.
 */
class DecomposedIntegration1stHalf : public BaseIntegration1stHalf
{
  public:
    explicit DecomposedIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~DecomposedIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        // including gravity and force from fluid
        Vecd acceleration = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd shear_force_ij = correction_factor_ * elastic_solid_.ShearModulus() *
                                  (J_to_minus_2_over_dimension_[index_i] + J_to_minus_2_over_dimension_[index_j]) *
                                  (pos_[index_i] - pos_[index_j]) / inner_neighborhood.r_ij_[n];
            acceleration += ((stress_on_particle_[index_i] + stress_on_particle_[index_j]) * inner_neighborhood.e_ij_[n] + shear_force_ij) *
                            inner_neighborhood.dW_ijV_j_[n] * inv_rho0_;
        }
        acc_[index_i] = acceleration;
    };

  protected:
    StdLargeVec<Real> J_to_minus_2_over_dimension_;
    StdLargeVec<Matd> stress_on_particle_, inverse_F_T_;
    const Real correction_factor_ = 1.07;
};

/**
 * @class Integration2ndHalf
 * @brief computing stress relaxation process by verlet time stepping
 * This is the second step
 */
class Integration2ndHalf : public BaseElasticIntegration
{
  public:
    explicit Integration2ndHalf(BaseInnerRelation &inner_relation)
        : BaseElasticIntegration(inner_relation){};
    virtual ~Integration2ndHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        const Vecd &vel_n_i = vel_[index_i];

        Matd deformation_gradient_change_rate = Matd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        const double epsilon = std::numeric_limits<double>::epsilon();

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            const Vecd vel_diff = vel_n_i - vel_[index_j];
            const double vel_squared_magnitude = vel_diff.squaredNorm();
            if (vel_squared_magnitude > epsilon)
            {
                Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                deformation_gradient_change_rate -= vel_diff * gradW_ij.transpose();
            }
        }

        dF_dt_[index_i] = deformation_gradient_change_rate * B_[index_i];
    };

    void update(size_t index_i, Real dt = 0.0);
};
} // namespace solid_dynamics
} // namespace SPH
#endif // VIRTOSIM_ELASTIC_DYNAMICS_H_C0E2274F_EA9C_462F_B9FA_890BC5B12F2B
