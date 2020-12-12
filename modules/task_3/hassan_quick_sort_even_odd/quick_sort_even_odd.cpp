// Copyright 2020 Hassan EzzAldeen
#include "../../../modules/task_3/hassan_quick_sort_even_odd/quick_sort_even_odd.h"
#include <vector>
#include <utility>
#include <algorithm>

std::vector<std::pair<int, int>> comparators;

std::vector<int> createRandomVector(int size_v) {
  std::mt19937 random;
  random.seed(static_cast<unsigned int>(time(0)));
  std::vector<int> vec(size_v);
  for (int i = 0; i < size_v; i++) {
    vec[i] = random() % 100;
  }
  return vec;
}

int part(std::vector<int>* vec, int first, int last) {
  int pivot = (*vec)[(first + last) / 2];
  int f = first;
  int  l = last;
  int tmp;
  while (f <= l) {
    while ((*vec)[f] < pivot) f++;
    while ((*vec)[l] > pivot) l--;
    if (f >= l)
      break;
    tmp = (*vec)[f];
    (*vec)[f] = (*vec)[l];
    (*vec)[l] = tmp;
    f++;
    l--;
  }
  return l;
}

void quickSort(std::vector<int>* vec, int first, int last) {
  if (first < last) {
    int pivot = part(vec, first, last);
    quickSort(vec, first, pivot);
    quickSort(vec, pivot + 1, last);
  }
}

void merge(std::vector<int> vec_up, std::vector<int> vec_down) {
  int vec_size = vec_up.size() + vec_down.size();
  if (vec_size == 1) {
    return;
  } else if (vec_size == 2) {
    comparators.push_back(std::pair<int, int>(vec_up[0], vec_down[0]));
    return;
  }

  std::vector<int> vec_up_odd, vec_up_even;
  std::vector<int> vec_down_odd, vec_down_even;
  std::vector<int> vec_result(vec_size);

  int vec_up_size = vec_up.size();
  for (int i = 0; i < vec_up_size; i++) {
    i % 2 ? vec_up_even.push_back(vec_up[i]) : vec_up_odd.push_back(vec_up[i]);
  }
  int vec_down_size = vec_down.size();
  for (int i = 0; i < vec_down_size; i++) {
    i % 2 ? vec_down_even.push_back(vec_down[i]) : vec_down_odd.push_back(vec_down[i]);
  }

  merge(vec_up_odd, vec_down_odd);
  merge(vec_up_even, vec_down_even);

  std::copy(vec_up.begin(), vec_up.end(), vec_result.begin());
  std::copy(vec_down.begin(), vec_down.end(), vec_result.begin() + vec_up.size());

  int res_size = vec_result.size();
  for (int i = 1; i < res_size - 1; i += 2) {
    comparators.push_back(std::pair<int, int>(vec_result[i], vec_result[i + 1]));
  }
}

void network(std::vector<int> procs) {
  int p_size = procs.size();
  if (p_size <= 1) {
    return;
  }
  std::vector<int> vec_up(p_size / 2);
  std::vector<int> vec_down(p_size / 2 + p_size % 2);

  std::copy(procs.begin(), procs.begin() + vec_up.size(), vec_up.begin());
  std::copy(procs.begin() + vec_up.size(), procs.end(), vec_down.begin());
  network(vec_up);
  network(vec_down);
  merge(vec_up, vec_down);
}

void quickSortBatcher(std::vector<int>* vec) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;

  int vec_size = vec->size();
  int fict = 0;
  while (vec_size % size) {
    vec->push_back(1000);
    vec_size++;
    fict++;
  }
  int delta = vec_size / size;
  std::vector<int> result(delta);
  std::vector<int> curr_vec(delta);
  std::vector<int> tmp_vec(delta);

  MPI_Scatter(&(*vec)[0], delta, MPI_INT, &result[0], delta, MPI_INT, 0, MPI_COMM_WORLD);
  std::sort(result.begin(), result.end());

  std::vector<int> procs(size);
  int s_procs = procs.size();
  for (int i = 0; i < s_procs; i++)
    procs[i] = i;
  network(procs);

  int comp_size = comparators.size();
  for (int i = 0; i < comp_size; i++) {
    std::pair<int, int> comparator = comparators[i];
    if (rank == comparator.first) {
      MPI_Send(result.data(), delta, MPI_INT, comparator.second, 0, MPI_COMM_WORLD);
      MPI_Recv(curr_vec.data(), delta, MPI_INT, comparator.second, 0, MPI_COMM_WORLD, &status);
      int resi = 0, curi = 0;
      for (int tmpi = 0; tmpi < delta; tmpi++) {
        int res = result[resi];
        int cur = curr_vec[curi];
        if (res < cur) {
          tmp_vec[tmpi] = res;
          resi++;
        } else {
          tmp_vec[tmpi] = cur;
          curi++;
        }
      }
      result = tmp_vec;
    } else if (rank == comparator.second) {
      MPI_Recv(curr_vec.data(), delta, MPI_INT, comparator.first, 0, MPI_COMM_WORLD, &status);
      MPI_Send(result.data(), delta, MPI_INT, comparator.first, 0, MPI_COMM_WORLD);
      int start = delta - 1;
      int resi = start, curi = start;
      for (int tmpi = start; tmpi >= 0; tmpi--) {
        int res = result[resi];
        int cur = curr_vec[curi];
        if (res > cur) {
          tmp_vec[tmpi] = res;
          resi--;
        } else {
          tmp_vec[tmpi] = cur;
          curi--;
        }
      }
      result = tmp_vec;
    }
  }
  MPI_Gather(&result[0], delta, MPI_INT, &(*vec)[0], delta, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0 && fict != 0) {
    vec->erase(vec->begin() + vec_size - fict, vec->end());
  }
  return;
}
