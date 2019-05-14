/* data_struct/low_level.h
*/
#ifndef LOW_LEVEL_H
#define LOW_LEVEL_H

int low_level_check_memleak();

typedef struct int_array{
    int n;
    int* data;
} int_array;

typedef struct double_array{
    int n;
    double* data;
} double_array;

/*******Allocate and free an array  *******/
int_array* int_array_alloc(int n);
double_array* double_array_alloc(int n);

void int_array_free(int_array* array);
void double_array_free(double_array* array);

/*******Set and get value from array*******/
int int_array_get(const int_array* array, int i);
double double_array_get(const double_array* array, int i);
void int_array_set(int_array* array, int i, int x);
void double_array_set(double_array* array, int i, double x);

/*******Initialize an array         *******/
void int_array_set_all(int_array* array, int x);
void double_array_set_all(double_array* array, double x);


/*******Get property from array     *******/
int int_array_get_size(const int_array* array);
int double_array_get_size(const double_array* array);

#endif
