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

#ifndef VIRTOSIM_INELASTIC_SOLID_H_D1F4B911_43D9_4317_A7AE_48A4E224EEBF
#define VIRTOSIM_INELASTIC_SOLID_H_D1F4B911_43D9_4317_A7AE_48A4E224EEBF
/**
 * @file 	inelastic_solid.h
 * @brief 	These are classes for define properties of elastic solid materials.
 *			These classes are based on isotropic linear elastic solid.
 * 			Several more complex materials, including neo-hookean, FENE noe-hookean
 *			and anisotropic muscle, are derived from the basic elastic solid class.
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "elastic_solid.h"
#include <unsupported/Eigen/Splines>

namespace SPH
{
/**
 * @class PlasticSolid
 * @brief Abstract class for a generalized plastic solid
 */
class PlasticSolid : public NeoHookeanSolid
{
  protected:
    Real yield_stress_;

  public:
    /** Constructor */
    explicit PlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress)
        : NeoHookeanSolid(rho0, youngs_modulus, poisson_ratio), yield_stress_(yield_stress)
    {
        material_type_name_ = "PlasticSolid";
    };
    virtual ~PlasticSolid() {};

    Real YieldStress() { return yield_stress_; };
    /** compute the stress through deformation, and plastic relaxation. */
    virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) = 0;

    virtual PlasticSolid *ThisObjectPtr() override { return this; };
};

/**
 * @class HardeningPlasticSolid
 * @brief Class for plastic solid with hardening
 */
class HardeningPlasticSolid : public PlasticSolid
{
  protected:
    Real hardening_modulus_;
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    StdLargeVec<Matd> inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */
    StdLargeVec<Real> hardening_parameter_;    /**< hardening parameter */

  public:
    /** Constructor */
    explicit HardeningPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress, Real hardening_modulus)
        : PlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress), hardening_modulus_(hardening_modulus)
    {
        material_type_name_ = "HardeningPlasticSolid";
    };
    virtual ~HardeningPlasticSolid() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    Real HardeningModulus() { return hardening_modulus_; };
    /** compute the stress through deformation, and plastic relaxation. */
    virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) override;

    virtual HardeningPlasticSolid *ThisObjectPtr() override { return this; };
};

class MultiLinearHardeningPlasticSolid : public PlasticSolid
{
  protected:
    const Real sqrt_2_over_3_ = sqrt(2.0 / 3.0);
    StdLargeVec<Matd> inverse_plastic_strain_; /**< inverse of plastic right cauchy green strain tensor */
    StdLargeVec<Real> hardening_parameter_;    /**< hardening parameter */
    StdLargeVec<Real> hardening_modulus_;      /**< hardening modulus */

    Eigen::Spline<Real, 1, 1> hardening_stresses_;  /**< the yield stress corresponding to different hardening parameters from the spline */
    Eigen::Spline<Real, 1, 1> hardening_moduluses_; /**< the hardening modulus corresponding to different hardening parameters from the spline */

  public:
    /** Constructor */
    explicit MultiLinearHardeningPlasticSolid(Real rho0, Real youngs_modulus, Real poisson_ratio, Real yield_stress,
                                              const Eigen::Spline<Real, 1, 1> &hardening_moduluses, const Eigen::Spline<Real, 1, 1> &hardening_stresses)
        : PlasticSolid(rho0, youngs_modulus, poisson_ratio, yield_stress),
          hardening_moduluses_(hardening_moduluses), hardening_stresses_(hardening_stresses)
    {
        material_type_name_ = "MultiLinearHardeningPlasticSolid";
    };
    virtual ~MultiLinearHardeningPlasticSolid() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    /** compute the stress through deformation, and plastic relaxation. */
    virtual Matd PlasticConstitutiveRelation(const Matd &deformation, size_t index_i, Real dt = 0.0) override;

    virtual MultiLinearHardeningPlasticSolid *ThisObjectPtr() override { return this; };
};
} // namespace SPH

#endif // VIRTOSIM_INELASTIC_SOLID_H_D1F4B911_43D9_4317_A7AE_48A4E224EEBF
