#include "particle_generator_lattice.h"

#include "base_body.h"
#include "base_mesh.h"
#include "base_particle_dynamics.h"
#include "base_particles.h"
#include "complex_geometry.h"
#include <tbb/mutex.h>
namespace SPH
{
//=================================================================================================//
void ParticleGenerator<BaseParticles, Lattice>::prepareGeometricData()
{
    Mesh mesh(domain_bounds_, lattice_spacing_, 0);
    Arrayi ncell = mesh.AllCells();
    Real volume = lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
    size_t total_cells = (size_t)ncell[0] * ncell[1] * ncell[2];

    // ——— 1) Count valid cells in parallel ———
    tbb::enumerable_thread_specific<size_t> tls_count(0);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, total_cells),
        [&](auto &r)
        {
            auto &local_count = tls_count.local();
            for (size_t idx = r.begin(); idx < r.end(); ++idx)
            {
                int i = idx / (ncell[1] * ncell[2]);
                int rem = idx % (ncell[1] * ncell[2]);
                int j = rem / ncell[2];
                int k = rem % ncell[2];
                Vecd p = mesh.CellPositionFromIndex({i, j, k});
                if (initial_shape_.checkNotFar(p, lattice_spacing_) &&
                    initial_shape_.checkContain(p))
                {
                    ++local_count;
                }
            }
        });

    // Sum up per‐thread counts
    size_t total_particles = std::accumulate(
        tls_count.begin(), tls_count.end(), size_t(0));

    // ——— 2) Pre‐allocate exactly needed space ———
    position_.clear();           // your internal array
    volumetric_measure_.clear(); // your internal volume array
    position_.reserve(total_particles);
    volumetric_measure_.reserve(total_particles);

    // Actually resize so we can write by index
    position_.resize(total_particles);
    volumetric_measure_.resize(total_particles);

    // ——— 3) Second pass: fill in place with atomic index ———
    std::atomic<size_t> write_idx{0};
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, total_cells),
        [&](auto &r)
        {
            for (size_t idx = r.begin(); idx < r.end(); ++idx)
            {
                int i = idx / (ncell[1] * ncell[2]);
                int rem = idx % (ncell[1] * ncell[2]);
                int j = rem / ncell[2];
                int k = rem % ncell[2];
                Vecd p = mesh.CellPositionFromIndex({i, j, k});
                if (initial_shape_.checkNotFar(p, lattice_spacing_) &&
                    initial_shape_.checkContain(p))
                {
                    size_t id = write_idx.fetch_add(1, std::memory_order_relaxed);
                    position_[id] = p;
                    volumetric_measure_[id] = volume;
                }
            }
        });

    std::cout << "finish ParticleGenerator: generated "
              << total_particles << " particles\n";
}
//=================================================================================================//
void ParticleGenerator<SurfaceParticles, Lattice>::prepareGeometricData()
{
    // Calculate the total volume and
    // count the number of cells inside the body volume, where we might put particles.
    Mesh mesh(domain_bounds_, lattice_spacing_, 0);
    Arrayi number_of_lattices = mesh.AllCells();
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
            for (int k = 0; k < number_of_lattices[2]; ++k)
            {
                Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j, k));
                if (initial_shape_.checkNotFar(particle_position, lattice_spacing_))
                {
                    if (initial_shape_.checkContain(particle_position))
                    {
                        all_cells_++;
                        total_volume_ += lattice_spacing_ * lattice_spacing_ * lattice_spacing_;
                    }
                }
            }
    Real number_of_particles = total_volume_ / avg_particle_volume_ + 0.5;
    planned_number_of_particles_ = int(number_of_particles);

    // initialize a uniform distribution between 0 (inclusive) and 1 (exclusive)
    std::mt19937_64 rng;
    std::uniform_real_distribution<Real> uniform_distr(0, 1);

    // Calculate the interval based on the number of particles.
    Real interval = planned_number_of_particles_ / (all_cells_ + TinyReal);
    if (interval <= 0)
        interval = 1; // It has to be lager than 0.

    // Add a particle in each interval, randomly. We will skip the last intervals if we already reach the number of particles.
    for (int i = 0; i < number_of_lattices[0]; ++i)
        for (int j = 0; j < number_of_lattices[1]; ++j)
            for (int k = 0; k < number_of_lattices[2]; ++k)
            {
                Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j, k));
                if (initial_shape_.checkNotFar(particle_position, lattice_spacing_))
                {
                    if (initial_shape_.checkContain(particle_position))
                    {
                        Real random_real = uniform_distr(rng);
                        // If the random_real is smaller than the interval, add a particle, only if we haven't reached the max. number of particles.
                        if (random_real <= interval && base_particles_.TotalRealParticles() < planned_number_of_particles_)
                        {
                            addPositionAndVolumetricMeasure(particle_position, avg_particle_volume_ / thickness_);
                            addSurfaceProperties(initial_shape_.findNormalDirection(particle_position), thickness_);
                        }
                    }
                }
            }
}
//=================================================================================================//
} // namespace SPH
