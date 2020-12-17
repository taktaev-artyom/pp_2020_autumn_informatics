// Copyright 2020 Makarov Alexander
#include <mpi.h>
#include <vector>
#include <random>
#include <iostream>
#include "../../../modules/task_3/makarov_a_gauss_filter_col/gauss_filter_col.h"


std::vector<unsigned short> generate_image(unsigned int w, unsigned int h) {
	unsigned int size = w * h;
	if (size == 0) return std::vector<unsigned short>();
	std::random_device rd;
    std::mt19937 gen(rd());
	size = w * h;
	std::vector<unsigned short> image(size);
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; j++)
			image[i * w + j] = static_cast<unsigned short>(gen() % 256);
	}
	return image;
}

std::vector<double> createGaussianKernel(double sigma, unsigned int radius) {
	int m_radius = static_cast<int>(radius);
	unsigned int size = 2 * m_radius + 1;
	std::vector<double> result(size * size);
	double norm = 0;
	
	for (int i = -m_radius; i <= m_radius; i++)
		for (int j = -m_radius; j <= m_radius; j++) {
			int idx = (i + m_radius) * size + (j + m_radius);
			result[idx] = exp(-(double)(i * i + j * j) / (double)(2. * sigma * sigma));
			norm += result[idx];
		}
	for (int i = -m_radius; i <= m_radius; i++)
		for (int j = -m_radius; j <= m_radius; j++) {
			int idx = (i + m_radius) * size + (j + m_radius);
			result[idx] /= norm;
		}
	return result;
}

std::vector<unsigned short> gaussFilterSequential(const std::vector<unsigned short>& image, unsigned int w, unsigned int h,
                                       double sigma, unsigned int radius) {
	unsigned int diam = radius * 2;
	std::vector<unsigned short> expanded_img((w + diam) * (h + diam));
	for (int unsigned i = 0; i < h; i++)
		for (int unsigned j = 0; j < w; j++)
			expanded_img[(radius + i) * (w + diam) + (radius + j)] = image[i * w + j];
	//frame
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < radius; j++) {
			expanded_img[(i + radius) * (w + diam) + j] = expanded_img[(i + radius) * (w + diam) + radius];
			expanded_img[(i + radius) * (w + diam) + (w + radius) + j] = expanded_img[(i + radius) * (w + diam) + (w + radius - 1)];
		}
	}
	for (unsigned int i = 0; i < w + diam; i++) {
		for (unsigned int j = 0; j < radius; j++){
			expanded_img[(w + diam) * j + i] = expanded_img[(w + diam) * radius + i];
			expanded_img[(h + radius + j) * (w + diam) + i] = expanded_img[(h + radius - 1)* (w + diam) + i];
		}
	}
	/*std::cout << "Original img" << std::endl;
	for (unsigned int i = 0; i < h ; i++) {
		for (unsigned int j = 0; j < w; j++)
			std::cout << image[i * w + j] << " ";
		std::cout << std::endl;
	}*/
	/*std::cout << "Expanded img" << std::endl;
	for (unsigned int i = 0; i < h + diam; i++) {
		for (unsigned int j = 0; j < w + diam; j++)
			std::cout << expanded_img[i * (w + diam) + j] << " ";
		std::cout << std::endl;
	}*/
	
	//create gaussian kernel
	std::vector<double> gauss_kernel = createGaussianKernel(sigma, radius);
	std::vector<unsigned short> result(h * w);
	int kern_radius = static_cast<int>(radius);
	int kern_size = 2 * kern_radius + 1;
	/*std::cout << "Gaussian kernel" << std::endl;
	for (unsigned int i = 0; i < kern_size; i++) {
		for (unsigned int j = 0; j < kern_size; j++)
			std::cout << gauss_kernel[i * kern_size + j] << " ";
		std::cout << std::endl;
	}*/
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; j++) {
			double conv = 0.;
			for (int k = -kern_radius; k <= kern_radius; k++)
				for (int t = -kern_radius; t <= kern_radius; t++)
					conv += (static_cast<double>(expanded_img[(kern_radius + static_cast<int>(i) + k) * (w + diam) + kern_radius + static_cast<int>(j) + t]) *
					         gauss_kernel[(k + kern_radius) * kern_size + t + kern_radius]);
			//std::cout << conv << " ";
			result[i * w + j] = static_cast<unsigned short>(conv);
		}
		//std::cout << std::endl;
	}
	/*std::cout << "Gaussian img" << std::endl;
	for (unsigned int i = 0; i < h ; i++) {
		for (unsigned int j = 0; j < w; j++)
			std::cout << result[i * w + j] << " ";
		std::cout << std::endl;
	}*/
	return result;
}

std::vector<unsigned short> gaussFilterParallel(const std::vector<unsigned short>& image, unsigned int w, unsigned int h,
                                     double sigma, unsigned int radius) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int proc_count;
	MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
	if (proc_count <= 1) return gaussFilterSequential(image, w, h, sigma, radius);
	unsigned int diam = radius * 2;
	unsigned int delta = w / (proc_count - 1);
	unsigned int remain = w % (proc_count - 1);
	
	unsigned int first_count = (h + diam) * (remain + delta + diam);
	unsigned int col_count = (h + diam) * (delta + diam);
	
	std::vector<double> gauss_kernel = createGaussianKernel(sigma, radius);
	int kern_radius = static_cast<int>(radius);
	int kern_size = kern_radius * 2 + 1;
	
	if (rank == 0){
		std::vector<unsigned short> t_image = prepare_image(image, w, h, radius);
		/*std::cout << "Delta = " << delta << std::endl;
		std::cout << "Remain = " << remain << std::endl;
		std::cout << "First count = " << first_count << std::endl;
		std::cout << "Proc_count = " << proc_count << std::endl;*/
		/*std::cout << "First send = ";
		for (int i = 0; i < first_count; i++)
			std::cout << t_image[i] << " ";*/
		//std::cout << std::endl;
		/*std::cout << "Transposed and expanded" << std::endl;
		for (unsigned int i = 0; i < w + diam; i++){
			for (unsigned int j = 0; j < h + diam; j++)
				std::cout << t_image[i * (h + diam) + j] << " ";
			std::cout << std::endl;
		}*/
		MPI_Send(t_image.data(), first_count, MPI_SHORT, 1, 0, MPI_COMM_WORLD); 
		for (int i = 2; i < proc_count; i++) {
			//start_pos = (h + 2) * (1 + remain + (i - 1) * delta - 1);
			unsigned int start_pos = (h + diam) * (remain + (i - 1) * delta);
			MPI_Send(t_image.data() + start_pos, col_count, MPI_SHORT, i, 0, MPI_COMM_WORLD);
		}
		MPI_Status status;
		std::vector<unsigned short> result(h * w);
		std::vector<unsigned short> tmp((remain + delta) * h);
		MPI_Recv(tmp.data(), (remain + delta) * h, MPI_SHORT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		/*std::cout << "Recieved result" << std::endl;
		for (unsigned int i = 0; i < h; i++)
			for (unsigned int j = 0; j < remain + delta; j++)
				std::cout << tmp[i * (remain + delta) + j];*/
		for (unsigned int i = 0; i < h; i++) {
			for (unsigned int j = 0; j < remain + delta; j++)
				result[i * w + j] = tmp[i * (remain + delta) + j];
		}
		for (int i = 2; i < proc_count; i++){
			unsigned int start_col = remain + delta * (i - 1);
			MPI_Recv(tmp.data(), delta * h, MPI_SHORT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			for (unsigned int j = 0; j < h; j++) {
				for (unsigned int k = 0; k < delta; k++)
					//result[w * k + ((i - 1) * delta + remain + j)] = tmp[j * h + k];
					result[j * w + start_col + k] = tmp[j * delta + k];
			}
		}
		//std::cout << "Gaussian parr" << std::endl;
		/*for (unsigned int i = 0; i < h; i++){
			for (unsigned int j = 0; j < w; j++)
				std::cout << result[i * w + j] << " ";
			std::cout << std::endl;
		}*/
		return result;
	}
	else {
		std::vector<unsigned short> local_image;
		std::vector<unsigned short> local_result;
		MPI_Status status;
		if (rank == 1) {
			local_image.resize(first_count);
			MPI_Recv(local_image.data(), first_count, MPI_SHORT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			/*std::cout << std::endl;
			std::cout << "First send = ";
			for (int i = 0; i < first_count; i++)
				std::cout << local_image[i] << " ";
			std::cout << std::endl;*/
			//std::cout << "Process " << rank << " local image" << std::endl;
			/*for (unsigned int i = 0; i < remain + delta + diam; i++){
				for (unsigned int j = 0; j < h + diam; j++)
					std::cout << local_image[i * (h + diam) + j] << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;*/
			local_result.resize(h * (remain + delta));
			for (unsigned int i = 0; i < remain + delta; i++) {
				for (unsigned int j = 0; j < h; j++) {
					//local_result[i * h + j] = local_image[(radius + i) * h + radius + j];
					//local_result[j * (remain + delta) + i] = local_image[(radius + i) * (h + diam) + radius + j];
					double conv = 0.;
					for (int k = -kern_radius; k <= kern_radius; k++) {
						for (int t = -kern_radius; t <=kern_radius; t++) {
							conv += (static_cast<double>(local_image[(i + kern_radius + k) * (h + diam) + kern_radius + j + t]) * 
							gauss_kernel[(k + kern_radius) * kern_size + t + kern_radius]);
							/*conv += (static_cast<double>(expanded_img[(kern_radius + static_cast<int>(i) + k) * (w + diam) + kern_radius + static_cast<int>(j) + t]) *
					         gauss_kernel[(k + kern_radius) * kern_size + t + kern_radius]);*/
						}
					}
					local_result[j * (remain + delta) + i] = static_cast<unsigned short>(conv);
					//local_result[j * (remain + delta) + i] = local_image[radius * (h + diam) + radius + i * (h + diam) + j];
				}
			}
			/*std::cout << "Local result" << std::endl;
			for (unsigned int i = 0; i < h; i++){
				for (unsigned int j = 0; j < remain + delta; j++)
					std::cout << local_result[i * (remain + delta) + j] << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;*/
			MPI_Send(local_result.data(), (remain + delta) * h, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
		}
		else {
			local_image.resize(col_count);
			MPI_Recv(local_image.data(), col_count, MPI_SHORT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			local_result.resize(h * delta);
			for (unsigned int i = 0; i < delta; i++) {
				for (unsigned int j = 0; j < h; j++) {
					double conv = 0.;
					for (int k = -kern_radius; k <= kern_radius; k++) {
						for (int t = -kern_radius; t <=kern_radius; t++) {
							conv += (static_cast<double>(local_image[(i + kern_radius + k) * (h + diam) + kern_radius + j + t]) * 
							gauss_kernel[(k + kern_radius) * kern_size + t + kern_radius]);
					//local_result[j * delta + i] = local_image[(radius + i) * (h + diam) + radius + j];
					local_result[j * delta + i] = static_cast<unsigned short>(conv);
						}
					}
				}
			}
			MPI_Send(local_result.data(), delta * h, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
		}
		return local_result;
	}
}

std::vector<unsigned short> prepare_image(const std::vector<unsigned short>& image, unsigned int w, unsigned int h,
                                          unsigned int radius) {
	//expand
	unsigned int diam = radius * 2;
	unsigned int old_size = w * h;
	unsigned int new_size = (w + diam) * (h + diam);
	std::vector<unsigned short> result(new_size);
	//transpose
	for (unsigned int i = 0; i < h; i++)
		for (unsigned int j = 0; j < w; j++)
			result[(j + radius) * (h + diam) + (i + radius)] = image[i * w + j];
	//frame
	for (unsigned int i = 0; i < w; i++) {
		for (unsigned int j = 0; j < radius; j++) {
			result[(i + radius) * (h + diam) + j] = result[(i + radius) * (h + diam) + radius];
			result[(i + radius) * (h + diam) + (h + radius) + j] = result[(i + radius) * (h + diam) + (h + radius - 1)];
		}
	}
	for (unsigned int i = 0; i < h + diam; i++) {
		for (unsigned int j = 0; j < radius; j++){
			result[(h + diam) * j + i] = result[(h + diam) * radius + i];
			result[(w + radius + j) * (h + diam) + i] = result[(w + radius - 1)* (h + diam) + i];
		}
	}
	return result;
}

/*std::vector<double> solveJacobiSequential(const std::vector<double>& A,
                                          const std::vector<double>& b) {
    int iter = 0;
    int size = static_cast<int>(b.size());
    std::vector<double> x(size);
    std::vector<double> x_prev(b);
    double error = 0.;
    do {
        for (int i = 0; i < size; i++) {
            x[i] = b[i];
            for (int j = 0; j < size; j++)
                if (i != j)
                    x[i] -= A[i * size + j] * x_prev[j];
            x[i] /= A[i * size + i];
        }
        for (int i = 0; i < size; i++)
            x_prev[i] = x[i];
        iter++;
        error = calculateError(A, x, b);
    } while (error > epsilon && iter < max_iter);
    return x;
}

std::vector<double> solveJacobiParallel(const std::vector<double>& A,
                                          const std::vector<double>& b) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int proc_count;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);

    int size;
    if (rank == 0)
        size = static_cast<int>(b.size());
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // создание коммуникатора без лишних процессов
    MPI_Comm comm;
    if (rank < size) {
        MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &comm);
    } else {
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &comm);
    }

    if (comm != MPI_COMM_NULL) {
        int new_rank;
        MPI_Comm_rank(comm, &new_rank);
        int new_proc_count;
        MPI_Comm_size(comm, &new_proc_count);

        // распределение количества строк по процессам
        int* counts_b = nullptr;
        int* displs_b = nullptr;
        int* counts_A = nullptr;
        int* displs_A = nullptr;
        if (new_rank == 0) {
            int delta = size / new_proc_count;
            int remaind = size % new_proc_count;
            counts_b = new int[new_proc_count];
            displs_b = new int[new_proc_count];
            counts_A = new int[new_proc_count];
            displs_A = new int[new_proc_count];
            for (int i = 0; i < new_proc_count; i++) {
                counts_b[i] = delta;
                displs_b[i] = 0;
            }
            for (int i = 0; i < remaind; i++)
                counts_b[i]++;
            for (int i = 1; i < new_proc_count; i++)
                displs_b[i] += displs_b[i - 1] + counts_b[i - 1];
            for (int i = 0; i < new_proc_count; i++) {
                counts_A[i] = counts_b[i] * size;
                displs_A[i] = displs_b[i] * size;
            }
        }

        double error = 0.;
        int iter = 0;
        // узнаем сколько получать строк
        int local_size;
        MPI_Scatter(counts_b, 1, MPI_INT, &local_size, 1, MPI_INT, 0, comm);

        std::vector<double> x_global;
        std::vector<double> x_local(local_size);
        std::vector<double> x_prev(size);
        if (rank == 0) {
            x_global.resize(size);
            x_prev = b;
        }
        // передать/получить элементы b
        std::vector<double> b_local(local_size);
        MPI_Scatterv(b.data(), counts_b, displs_b, MPI_DOUBLE, b_local.data(),
                     local_size, MPI_DOUBLE, 0, comm);
        // передать/получить строки А
        std::vector<double> A_local(local_size * size);
        MPI_Scatterv(A.data(), counts_A, displs_A, MPI_DOUBLE, A_local.data(),
                     local_size * size, MPI_DOUBLE, 0, comm);
        // получить номер первой строки (для диагональных элементов)
        int first_row_num;
        MPI_Scatter(displs_b, 1, MPI_INT, &first_row_num, 1, MPI_INT, 0, comm);

        do {
            // передать/получить результат предыдущей итерации
            MPI_Bcast(x_prev.data(), size, MPI_DOUBLE, 0, comm);
            // обработать свои куски
            for (int i = 0; i < local_size; i++) {
                x_local[i] = b_local[i];
                for (int j = 0; j < size; j++)
                    if (j != first_row_num + i)
                        x_local[i] -= A_local[i * size + j] * x_prev[j];
                x_local[i] /= A_local[i * size + first_row_num + i];
            }
            iter++;

            // собрать результаты этой итерации
            MPI_Gatherv(x_local.data(), local_size, MPI_DOUBLE,
                        x_global.data(), counts_b, displs_b,
                        MPI_DOUBLE, 0, comm);
            if (rank == 0) {
                for (int i = 0; i < size; i++)
                    x_prev[i] = x_global[i];
                error = calculateError(A, x_global, b);
            }
            // передать/получить ошибку
            MPI_Bcast(&error, 1, MPI_DOUBLE, 0, comm);
        } while (error > epsilon && iter < max_iter);
        // принять результат от всех процессоров
        MPI_Comm_free(&comm);
        if (rank == 0) {
            delete[] counts_b;
            delete[] displs_b;
            delete[] counts_A;
            delete[] displs_A;
        }
        return x_global;
    } else {
        return std::vector<double>();
    }
}*/
