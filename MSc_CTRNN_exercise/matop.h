#ifdef _MATOP_H
#define _MATOP_H

inline void random_vector (int size, double vector[size]);
inline void zero_vector (int size, double vector[size]);
inline void flatten_vector(double MAX, double MIN, int size, double vector[size]);
void normalize_vector(double ZERO, int size, double vector[size]);
inline void copy_vector(int size, double src[size], double dest[size]);
static int cmpdbl(const void *a, const void *b);
inline void sort_vector(int size, double vector[size]);
inline void print_vector (int size, double vector[size]);

void random_matrix (int size, double matrix[size][size]);
void zero_matrix (int size, double matrix[size][size]);
void flatten_matrix(double MAX, double MIN, int size, double matrix[size][size]);
void normalize_matrix(double ZERO, int size, double matrix[size][size]);
void copy_matrix(int size, double src[size][size], double dest[size][size]);
void print_matrix(int size, double matrix[size][size]);

#endif
