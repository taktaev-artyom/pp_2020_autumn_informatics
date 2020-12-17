// Copyright 2020 Makarov Alexander
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "./gauss_filter_col.h"

TEST(GaussFilter, 7x3_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 7;
	unsigned int h = 3;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
	//std::cout << "Process rank " << rank << " done" << std::endl;
    if (rank == 0) {
		/*std::cout << "Original" << std::endl;
		std::cout << "w = " << w << std::endl;
		std::cout << "h = " << h << std::endl;
		std::cout << "size = " << size << std::endl;
		for (unsigned int i = 0; i < h; i++){
			for (unsigned int j = 0; j < w; j++)
				std::cout << image[i * w + j] << " ";
			std::cout << std::endl;
		}*/
		/*std::cout << std::endl;
		std::cout << "Transposed and expanded" << std::endl;
		std::vector<unsigned short> t_image = prepare_image(image, w, h, radius);
		for (unsigned int i = 0; i < w + diam; i++){
			for (unsigned int j = 0; j < h + diam; j++)
				std::cout << t_image[i * (h + diam) + j] << " ";
			std::cout << std::endl;
		}*/
		/*std::cout << std::endl;
		std::cout << "Parallel" << std::endl;
		for (unsigned int i = 0; i < h; i++){
			for (unsigned int j = 0; j < w; j++)
				std::cout << result_par[i * w + j] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;*/
		double parr_time = end_time - start_time;
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
		double seq_time = end_time - start_time;
		std::cout << "Sequential time = " << seq_time << " s" << std::endl;
		std::cout << "Parallel time = " << parr_time << " s" << std::endl;
        ASSERT_EQ(result_par, result_seq);
    }
}

TEST(GaussFilter, 100x200_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 100;
	unsigned int h = 200;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
    if (rank == 0) {
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
        ASSERT_EQ(result_par, result_seq);
    }
}

TEST(GaussFilter, 300x200_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 300;
	unsigned int h = 200;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
    if (rank == 0) {
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
        ASSERT_EQ(result_par, result_seq);
    }
}

TEST(GaussFilter, 1000x2000_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 1000;
	unsigned int h = 2000;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
    if (rank == 0) {
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
        ASSERT_EQ(result_par, result_seq);
    }
}

TEST(GaussFilter, 3000x2000_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 3000;
	unsigned int h = 2000;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
    if (rank == 0) {
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
        ASSERT_EQ(result_par, result_seq);
    }
}

TEST(GaussFilter, 10000x10000_generation_test) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double sigma = 1.;
	unsigned int radius = 1;
	unsigned int diam = radius * 2;
	unsigned int w = 10000;
	unsigned int h = 10000;
	unsigned int size = h * w;
	double start_time, end_time;
	std::vector<unsigned short> image;
	if (rank == 0) {
		image = generate_image(w, h);
	}
	if (rank == 0) start_time = MPI_Wtime();
	std::vector<unsigned short> result_par = gaussFilterParallel(image, w, h, sigma, radius);
	if (rank == 0) end_time = MPI_Wtime();
    if (rank == 0) {
		start_time = MPI_Wtime();
		std::vector<unsigned short> result_seq = gaussFilterSequential(image, w, h, sigma, radius);
		end_time = MPI_Wtime();
        ASSERT_EQ(result_par, result_seq);
    }
}

/*TEST(JacobiMethod, My_SLE_sequent) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::vector<double> A = {
            100., 19., -27.,
            -3., -78., 24.,
            18., 100., 274
        };
        std::vector<double> b = {
            357., -71., 44.,
        };

        double start_time = MPI_Wtime();
        std::vector<double> x = solveJacobiSequential(A, b);
        double end_time = MPI_Wtime();

        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, My_SLE_parallel) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        A = {
            100., 19., -27.,
            -3., -78., 24.,
            18., 100., 274.
        };
        b = {
            357., -71., 44.,
        };
    }
    double start_time, end_time;
    if (rank == 0) start_time = MPI_Wtime();
    std::vector<double> x = solveJacobiParallel(A, b);
    if (rank == 0) end_time = MPI_Wtime();

    if (rank == 0) {
        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 10_SLE_sequent) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int size = 10;
        std::vector<double> A = generate_A(size);
        std::vector<double> b = generate_b(size);

        double start_time = MPI_Wtime();
        std::vector<double> x = solveJacobiSequential(A, b);
        double end_time = MPI_Wtime();

        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 10_SLE_parallel) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        int size = 10;
        A = generate_A(size);
        b = generate_b(size);
    }
    double start_time, end_time;
    if (rank == 0) start_time = MPI_Wtime();
    std::vector<double> x = solveJacobiParallel(A, b);
    if (rank == 0) end_time = MPI_Wtime();

    if (rank == 0) {
        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 100_SLE_sequent) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int size = 100;
        std::vector<double> A = generate_A(size);
        std::vector<double> b = generate_b(size);

        double start_time = MPI_Wtime();
        std::vector<double> x = solveJacobiSequential(A, b);
        double end_time = MPI_Wtime();

        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 100_SLE_parallel) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        int size = 100;
        A = generate_A(size);
        b = generate_b(size);
    }
    double start_time, end_time;
    if (rank == 0) start_time = MPI_Wtime();
    std::vector<double> x = solveJacobiParallel(A, b);
    if (rank == 0) end_time = MPI_Wtime();

    if (rank == 0) {
        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 1000_SLE_sequent) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int size = 1000;
        std::vector<double> A = generate_A(size);
        std::vector<double> b = generate_b(size);

        double start_time = MPI_Wtime();
        std::vector<double> x = solveJacobiSequential(A, b);
        double end_time = MPI_Wtime();

        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << std::endl;
        std::cout << "\tTime: " << std::fixed << time << " s" << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 1000_SLE_parallel) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        int size = 1000;
        A = generate_A(size);
        b = generate_b(size);
    }
    double start_time, end_time;
    if (rank == 0) start_time = MPI_Wtime();
    std::vector<double> x = solveJacobiParallel(A, b);
    if (rank == 0) end_time = MPI_Wtime();

    if (rank == 0) {
        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 5000_SLE_sequent) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int size = 5000;
        std::vector<double> A = generate_A(size);
        std::vector<double> b = generate_b(size);

        double start_time = MPI_Wtime();
        std::vector<double> x = solveJacobiSequential(A, b);
        double end_time = MPI_Wtime();

        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << std::endl;
        std::cout << "\tTime: " << std::fixed << time << " s" << std::endl;
        ASSERT_LE(error, epsilon);
    }
}

TEST(JacobiMethod, 5000_SLE_parallel) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> A;
    std::vector<double> b;
    if (rank == 0) {
        int size = 5000;
        A = generate_A(size);
        b = generate_b(size);
    }
    double start_time, end_time;
    if (rank == 0) start_time = MPI_Wtime();
    std::vector<double> x = solveJacobiParallel(A, b);
    if (rank == 0) end_time = MPI_Wtime();

    if (rank == 0) {
        double error = calculateError(A, x, b);
        double time = end_time - start_time;
        std::cout << "\tError: " << std::fixed << error << " s" << std::endl;
        std::cout << "\tTime: " << std::fixed << time << std::endl;
        ASSERT_LE(error, epsilon);
    }
}*/

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
