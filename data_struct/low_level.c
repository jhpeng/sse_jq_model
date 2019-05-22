/* data_struct/low_level.c
*/
#include <stdlib.h>
#include <stdio.h>

#include "low_level.h"

/*The mode trying to debug and check the memory leakage.
* Please turn it off after testing the code.
*/
#if 0
#define DEBUG_MODE_LOW_LEVEL
#define CHECK_MEM_LEAK
#endif

static int int_array_count=0;
static int double_array_count=0;

int low_level_check_memleak(){
    int check=0;
    if(int_array_count!=0 || double_array_count!=0) check=1;

    return check;
}

int_array* int_array_alloc(int n){
    int_array* a = (int_array*)malloc(sizeof(int_array));
    a->data = (int*)malloc(n*sizeof(int));
    a->n = n;

#ifdef CHECK_MEM_LEAK
    int_array_count+=1;
#endif

    return a;
}

double_array* double_array_alloc(int n){
    double_array* a = (double_array*)malloc(sizeof(double_array));
    a->data = (double*)malloc(n*sizeof(double));
    a->n = n;

#ifdef CHECK_MEM_LEAK
    double_array_count+=1;
#endif

    return a;
}

void int_array_free(int_array* array){
    free(array->data);
    free(array);

#ifdef CHECK_MEM_LEAK
    int_array_count-=1;
#endif

}

void double_array_free(double_array* array){
    free(array->data);
    free(array);

#ifdef CHECK_MEM_LEAK
    double_array_count-=1;
#endif

}

int int_array_get(const int_array* array, int i){
#ifdef DEBUG_MODE_LOW_LEVEL
    if(i>=array->n){
        printf("low_level : i out of range!\n");
    }
#endif
    return array->data[i];
}

double double_array_get(const double_array* array, int i){
#ifdef DEBUG_MODE_LOW_LEVEL
    if(i>=array->n){
        printf("low_level : i out of range!\n");
    }
#endif
    return array->data[i];
}

void int_array_set(int_array* array, int i, int x){
#ifdef DEBUG_MODE_LOW_LEVEL
    if(i>=array->n){
        printf("low_level : i out of range!\n");
    }
#endif
    array->data[i] = x;
}

void double_array_set(double_array* array, int i, double x){
#ifdef DEBUG_MODE_LOW_LEVEL
    if(i>=array->n){
        printf("low_level : i out of range!\n");
    }
#endif
    array->data[i] = x;
}

void int_array_set_all(int_array* array, int x){
    int i;
    for(i=0;i<array->n;++i){
        int_array_set(array,i,x);
    }
}

void double_array_set_all(double_array* array, double x){
    int i;
    for(i=0;i<array->n;++i){
        double_array_set(array,i,x);
    }
}

int int_array_get_size(const int_array* array){
    return array->n;
}

int double_array_get_size(const double_array* array){
    return array->n;
}

#if 0
int main(){
    int i;
    for(i=0;i<1000;++i){
        int_array* a = int_array_alloc(1000000);
        double_array* b = double_array_alloc(1000000);
    
        int_array_set_all(a,-1);
        double_array_set_all(b,0.2198);

        int_array_free(a);
        double_array_free(b);

        if(low_level_check_memleak()){
            printf("memleak!\n");
        }
    }

    return 0;
}
#endif
