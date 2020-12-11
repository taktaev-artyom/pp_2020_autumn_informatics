// Copyright 2020 Grigoryan Garry
#ifndef MODULES_TASK_3_GRIGORYAN_G_DOUBLE_SORTING_SORT_H_
#define MODULES_TASK_3_GRIGORYAN_G_DOUBLE_SORTING_SORT_H_

bool is_sorted(double* source, int size);
void passing(double* source, double* dest, int size, int offset);
void last_passing(double* source, double* dest, int size, int offset);
void ordered_merge(double* source1, int size1, double* source2, int size2, double* dest);
void seq_sorting(double* source, int size);
void par_sorting(double** source, int size);

#endif  // MODULES_TASK_3_GRIGORYAN_G_DOUBLE_SORTING_SORT_H_
