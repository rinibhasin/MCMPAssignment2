#define __USE_C99_MATH
#include <stdbool.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

long int product(int *, int);
float * read_array(char *, int *, int);
int *read_dims(char *);

int main(int argc, char *argv[]) {
    if(argc != 3) {
        printf("Usage: %s <filename_check> <filename_ref>\n", argv[0]);
        return -1;
    }

    char check_filename[500];
    char ref_filename[500];

    strcpy(check_filename, argv[1]);
    strcpy(ref_filename, argv[2]);

    int *c_dims_original = read_dims(check_filename);
    if(c_dims_original == NULL) {
        return -1;
    }
    int c_num_dims = c_dims_original[0];
    int *c_dims = c_dims_original+1;

    int *r_dims_original = read_dims(ref_filename);
    if(r_dims_original == NULL) {
        return -1;
    }
    int r_num_dims = r_dims_original[0];
    int *r_dims = r_dims_original+1;

    if(r_num_dims != c_num_dims) {
        printf("Number of dimensions does not match. Expected %d, but found %d. \n", r_num_dims, c_num_dims);
        return 1;
    }

    for(int i=0; i<r_num_dims; i++) {
        if(r_dims[i] != c_dims[i]) {
            printf("Dimension %d does not match. Expected %d but found %d.\n", i, r_dims[i], c_dims[i]);
            return 1;
        } 
            
    }

    float *c_data = read_array(check_filename, c_dims, c_num_dims);

    if(c_data == NULL) {
        return -1;
    }

    float *r_data = read_array(ref_filename, r_dims, r_num_dims);

    if(r_data == NULL) {
        return -1;
    }

    long output_size = product(r_dims, r_num_dims);
    int match = 1;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(&&:match)
    #endif
    for(int i=0;i<output_size; i++) {
        if(fabs(c_data[i]-r_data[i])>0.1) {
            match = 0;
            printf("At position %d, expected %f but found %f. \n", i, r_data[i], c_data[i]);
            break;
        }
    }

    free(c_dims_original);
    free(r_dims_original);
    free(c_data);
    free(r_data);
    if(match == 1) {
        return 0;
    } else{
        return 1;
    }
}
