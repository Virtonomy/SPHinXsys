/**
 * @file 	particle_generator_lattic_supplementary.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "particle_generator_lattice.h"

#include "complex_shape.h"
#include "base_mesh.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH {
	//=================================================================================================//
	void ParticleGeneratorLattice::createBaseParticles(BaseParticles* base_particles)
	{
		BaseMesh mesh(domain_bounds_, lattice_spacing_, 0);
		Real particle_volume = lattice_spacing_ * lattice_spacing_;
		Vecu number_of_lattices = mesh.NumberOfCellsFromNumberOfGridPoints(mesh.NumberOfGridPoints());
			for (size_t i = 0; i < number_of_lattices[0]; ++i)
				for (size_t j = 0; j < number_of_lattices[1]; ++j)
				{
				Vecd particle_position = mesh.CellPositionFromIndex(Vecu(i,j));
				if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
				{
					if (body_shape_->checkContain(particle_position))
					{
						createABaseParticle(base_particles, particle_position, particle_volume);
					}
				}
			}
	}
	//=================================================================================================//
	void ShellParticleGeneratorLattice::createBaseParticles(BaseParticles* base_particles)
	{
		// Calculate the total volume and
		// count the number of cells inside the body volume, where we might put particles.
		std::unique_ptr<BaseMesh> mesh(new BaseMesh(domain_bounds_, lattice_spacing_, 0));
		Vecu number_of_lattices = mesh->NumberOfCellsFromNumberOfGridPoints(mesh->NumberOfGridPoints());
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j) {
				Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j));
				if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
				{
					if (body_shape_->checkContain(particle_position))
					{
						number_of_cells_++;
						total_volume_ += lattice_spacing_ * lattice_spacing_;
					}
				}
			}
		Real number_of_particles = total_volume_ / avg_particle_volume_ + 0.5;
		planned_number_of_particles_ = int(number_of_particles);

		// Calculate the interval based on the number of particles.
		Real interval = (Real)planned_number_of_particles_ / (Real)number_of_cells_;
		if (interval <= 0) interval = 1;          // It has to be lager than 0.
		Real random_real = 0.0;
		// Add a particle in each interval, randomly. We will skip the last intervals if we already reach the number of particles
		for (size_t i = 0; i < number_of_lattices[0]; ++i)
			for (size_t j = 0; j < number_of_lattices[1]; ++j)
			{
				Vecd particle_position = mesh->CellPositionFromIndex(Vecu(i, j));
				if (body_shape_->checkNotFar(particle_position, lattice_spacing_))
				{
					if (body_shape_->checkContain(particle_position))
					{
						random_real = (Real)rand() / (RAND_MAX);
						// If the random_real is smaller than the interval, add a particle, only if we haven't reached the max. number of particles
						if (random_real <= interval && base_particles->total_real_particles_ < planned_number_of_particles_)
						{
							createABaseParticle(base_particles, particle_position, avg_particle_volume_);
						}
					}
				}
			}
	}
	//=================================================================================================//
}
