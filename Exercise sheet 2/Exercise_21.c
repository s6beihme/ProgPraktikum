#include <stdio.h>
#include <stdlib.h>

//a)
struct _vector {
	int size;
	double* data;
};

typedef struct _vector* Vector;

//function to create a vector. first argument is the size of wanted vector, second argument is a pointer to the vector object
int vec_create(int size, Vector* vector) {
	if (size < 0 || vector == NULL) {
		printf("vec_create failed due to invalid input!\n");
		return -1;
	}
	(*vector) = malloc(sizeof((**vector))); //I do not understand, why this works and malloc(sizeof(Vector)) doesnt
	if ((*vector) == NULL) return-1;

	(*vector)->size = size;
	(*vector)->data = malloc(size * sizeof(double));
	if ((*vector)->data == NULL) {
		free((*vector));
		return -1;
	}

	if ((*vector)->data == NULL) {
		free((*vector));
		return -1;
	}
	return 0;
}

void vec_free(Vector* vector) {
	if (vector == NULL || (*vector) == NULL) return;

	free((*vector)->data);
	free((*vector));
	(*vector) = NULL;
}

//b)

// function to assign set 'data' of vector to values 'vals'. size is the size of the vals array
void vec_assemble(Vector vec, double* vals, int size) {
	if (size != vec->size) {
		printf("vec_assemble failed, because third argument didnt correspond to size of vector");
		return;
	}
	if (vec->data != NULL) {
		free(vec->data);
	}
	vec->data = vals;
}

//c)

//function to compute dot product of vectors v1 and v2
double vec_dot(Vector v1, Vector v2) {
	if (v1->size != v2->size) {
		printf("vec_dot failed, bacause vectors didnt have same size");
		exit(0);
	}
	double res = 0;
	for (int i = 0; i < v1->size; i++) {
		res += v1->data[i] * v2->data[i];
	}
	return res;
}

int main() {
	printf("H1\n");
	Vector v1, v2;
	printf("H2\n");
	vec_create(-1, &v1);
	vec_create(3, &v2);

	double a[3] = { 1,2,3 };
	double b[3] = { 2,2,2 };
	vec_assemble(v1, a, 3);
	vec_assemble(v2, b, 3);
	double d = vec_dot(v1, v2);
	printf("v1*v2 = %lf", d);
	vec_free(&v1);
	vec_free(&v2);
	return 0;
}
