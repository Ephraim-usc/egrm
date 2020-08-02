#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <threads.h>

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct matrix {
  ITYPE n;
  DTYPE *data;
} matrix;

matrix* new_matrix(ITYPE n)
{
  matrix* mat;
  mat = (matrix *)malloc(sizeof(matrix));
  mat->n = n;
  mat->data = calloc(n*n,sizeof(DTYPE));
  return mat;
}

void destry_matrix(matrix* mat)
{
  free(mat->data);
  free(mat);
}

void destry_matrix(matrix* mat)
{
  free(mat->data);
  free(mat);
}

void print_matrix(matrix* mat)
{
  ITYPE n = mat->n;
  printf("matrix of dimension %lu\n", n);
  
  ITYPE i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      printf("%lf ", mat->data[i*n + j]);
    printf("\n");
  }
}

void add_square(matrix* mat, ITYPE* idx, ITYPE len_idx, DTYPE q)
{
  ITYPE* x = idx;
  ITYPE* y = idx;
  ITYPE* end = idx + len_idx;
  DTYPE* data = mat->data;
  ITYPE n = mat->n;
  for (; x < end; x++)
    for (y=idx; y <= x; y++)
    {
      data[*x * n + *y] += q;
    }
}

static void parse_py_int_seq(PyObject *py_int_seq, ITYPE** pr, ITYPE* len)
{
  *len = (ITYPE)PySequence_Fast_GET_SIZE(py_int_seq);
  *pr = (ITYPE *) malloc(sizeof(ITYPE) * *len);
  ITYPE i;
  for (i = 0; i < *len; i++) 
  {
    PyObject *item = PySequence_Fast_GET_ITEM(py_int_seq, i);
    (*pr)[i] = (ITYPE)PyLong_AsLong(item);
  }
}

static PyObject* py_new_matrix(PyObject* self, PyObject* args)
{ 
  ITYPE n;
  PyArg_ParseTuple(args, "k", &n);
  
  matrix* mat = new_matrix((ITYPE)n);
  PyObject* mat_py = PyCapsule_New((void *)mat, "matrix._matrix_C_API", NULL);
  return mat_py;
}

static PyObject* py_destroy_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  destroy_matrix(mat);
  
  Py_RETURN_NONE;
}

static PyObject* py_add_square(PyObject* self, PyObject* args)
{
  PyObject* py_mat;
  PyObject* py_int_seq;
  PyObject* py_q;
  PyArg_UnpackTuple(args, NULL, 3, 3, &py_mat, &py_int_seq, &py_q);
  
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  ITYPE* idx; ITYPE len_idx;
  DTYPE q = (ITYPE)PyLong_AsLong(py_q);
  py_int_seq = PySequence_Fast(py_int_seq, NULL);
  
  parse_py_int_seq(py_int_seq, &idx, &len_idx);
  add_square(mat, idx, len_idx, q);
  
  Py_RETURN_NONE;
}

static PyObject* py_print_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  
  print_matrix(mat);
  Py_RETURN_NONE;
}

static PyObject* py_export_matrix(PyObject* self, PyObject* args)
{
  PyObject* py_mat; 
  PyArg_UnpackTuple(args, NULL, 1, 1, &py_mat);
  matrix* mat = (matrix *)PyCapsule_GetPointer(py_mat, "matrix._matrix_C_API");
  DTYPE* data = mat->data;
  
  PyListObject *py_list = (PyListObject *) Py_BuildValue("[]");
  int i = 0;
  for (i; i < mat->N; i++)
    PyList_Append(py_list, Py_BuildValue("i", data[i]));
  
  return (PyObject *) py_list;
}

static PyMethodDef myMethods[] = 
{
  {"new_matrix", py_new_matrix, METH_VARARGS, "new matrix"},
  {"print_matrix", py_print_matrix, METH_VARARGS, "print matrix"},
  {"add_square", py_add_square, METH_VARARGS, "add_square"},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef matrixModule =
{
  PyModuleDef_HEAD_INIT,
  "matrixModule",
  "matrix Module",
  -1,
  myMethods
};

PyMODINIT_FUNC PyInit_matrix(void)
{
  return PyModule_Create(&matrixModule);
}

/*
int main()
{
  //matrix* mat = newMatrix(10);
  //int idx[3] = {2,5,7}; 
  //addSquare(mat, idx, 3, 9);
  printf("%d", fib_(10));
  return 0;
}
*/
