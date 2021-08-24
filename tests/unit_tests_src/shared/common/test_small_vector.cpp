#include <gtest/gtest.h>
#include "small_vectors.h"

using namespace SPH;

TEST(AngleBetweenTwo3DVectors, getAngleBetweenTwo3DVectors)
{
    Vec3d vector_1 = Vec3d(2, -4, 0);
    Vec3d vector_2 = Vec3d(3, 2, 5);

    Real cos_teta = getAngleBetweenTwo3DVectors(vector_1, vector_2);

    Real cos_teta_ref = -2 / (sqrt(20) * sqrt(38));
 
		EXPECT_EQ(cos_teta, cos_teta_ref);

}
//=================================================================================================//
TEST(VectorProjectionOf3DVector, getVectorProjectionOf3DVector)
{
    Vec3d vector_1 = Vec3d(1, 2, 3);
    Vec3d vector_2 = Vec3d(2, -1, 4);
    
    Vec3d proj_vector_1 = getVectorProjectionOf3DVector(vector_1, vector_2);

    Vec3d proj_vector_1_ref = Vec3d(1.142857, -0.571428, 2.285714);

	for (size_t i = 0; i < 3; i++)
	{
        EXPECT_NEAR(proj_vector_1[i], proj_vector_1_ref[i], 1e-6);
	}

}
//=================================================================================================//
int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}