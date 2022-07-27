/* File: fvframemodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Tue Jun  9 09:01:40 2020
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *fvframe_error;
static PyObject *fvframe_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define CHECKSCALAR(check,tcheck,name,show,var)\
    if (!(check)) {\
        char errstring[256];\
        sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
        PyErr_SetString(fvframe_error,errstring);\
        /*goto capi_fail;*/\
    } else 
#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
    PyObject* tmp = NULL;
    if (PyInt_Check(obj)) {
        *v = (int)PyInt_AS_LONG(obj);
        return 1;
    }
    tmp = PyNumber_Int(obj);
    if (tmp) {
        *v = PyInt_AS_LONG(tmp);
        Py_DECREF(tmp);
        return 1;
    }
    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyString_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj,0);
    if (tmp) {
        PyErr_Clear();
        if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = fvframe_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(build,BUILD)(int*,int*,double*,double*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/*********************************** build ***********************************/
static char doc_f2py_rout_fvframe_build[] = "\
pyvarea,pyngbs,pyvlenght,pydlenght,pymaxngbh = build(pyoids,pygids,pytx,pyty,pytedge,pytelmt,pyvx,pyvy,pyvedge,[pygnodes,pyonodes,pydnodes,pydedges,pydelems,pyvnodes,pyvedges])\n\nWrapper for ``build``.\
\n\nParameters\n----------\n"
"pyoids : input rank-1 array('i') with bounds (pyonodes)\n"
"pygids : input rank-1 array('i') with bounds (pygnodes)\n"
"pytx : input rank-1 array('d') with bounds (pydnodes)\n"
"pyty : input rank-1 array('d') with bounds (pydnodes)\n"
"pytedge : input rank-2 array('i') with bounds (pydedges,2)\n"
"pytelmt : input rank-2 array('i') with bounds (pydelems,3)\n"
"pyvx : input rank-1 array('d') with bounds (pyvnodes)\n"
"pyvy : input rank-1 array('d') with bounds (pyvnodes)\n"
"pyvedge : input rank-2 array('i') with bounds (pyvedges,2)\n"
"\nOther Parameters\n----------------\n"
"pygnodes : input int, optional\n    Default: len(pygids)\n"
"pyonodes : input int, optional\n    Default: len(pyoids)\n"
"pydnodes : input int, optional\n    Default: len(pytx)\n"
"pydedges : input int, optional\n    Default: shape(pytedge,0)\n"
"pydelems : input int, optional\n    Default: shape(pytelmt,0)\n"
"pyvnodes : input int, optional\n    Default: len(pyvx)\n"
"pyvedges : input int, optional\n    Default: shape(pyvedge,0)\n"
"\nReturns\n-------\n"
"pyvarea : rank-1 array('d') with bounds (pydnodes)\n"
"pyngbs : rank-2 array('i') with bounds (pydnodes,20)\n"
"pyvlenght : rank-2 array('d') with bounds (pydnodes,20)\n"
"pydlenght : rank-2 array('d') with bounds (pydnodes,20)\n"
"pymaxngbh : int";
/* extern void F_FUNC(build,BUILD)(int*,int*,double*,double*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*); */
static PyObject *f2py_rout_fvframe_build(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,double*,double*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*,int*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  int *pyoids = NULL;
  npy_intp pyoids_Dims[1] = {-1};
  const int pyoids_Rank = 1;
  PyArrayObject *capi_pyoids_tmp = NULL;
  int capi_pyoids_intent = 0;
  PyObject *pyoids_capi = Py_None;
  int *pygids = NULL;
  npy_intp pygids_Dims[1] = {-1};
  const int pygids_Rank = 1;
  PyArrayObject *capi_pygids_tmp = NULL;
  int capi_pygids_intent = 0;
  PyObject *pygids_capi = Py_None;
  double *pytx = NULL;
  npy_intp pytx_Dims[1] = {-1};
  const int pytx_Rank = 1;
  PyArrayObject *capi_pytx_tmp = NULL;
  int capi_pytx_intent = 0;
  PyObject *pytx_capi = Py_None;
  double *pyty = NULL;
  npy_intp pyty_Dims[1] = {-1};
  const int pyty_Rank = 1;
  PyArrayObject *capi_pyty_tmp = NULL;
  int capi_pyty_intent = 0;
  PyObject *pyty_capi = Py_None;
  int *pytedge = NULL;
  npy_intp pytedge_Dims[2] = {-1, -1};
  const int pytedge_Rank = 2;
  PyArrayObject *capi_pytedge_tmp = NULL;
  int capi_pytedge_intent = 0;
  PyObject *pytedge_capi = Py_None;
  int *pytelmt = NULL;
  npy_intp pytelmt_Dims[2] = {-1, -1};
  const int pytelmt_Rank = 2;
  PyArrayObject *capi_pytelmt_tmp = NULL;
  int capi_pytelmt_intent = 0;
  PyObject *pytelmt_capi = Py_None;
  double *pyvx = NULL;
  npy_intp pyvx_Dims[1] = {-1};
  const int pyvx_Rank = 1;
  PyArrayObject *capi_pyvx_tmp = NULL;
  int capi_pyvx_intent = 0;
  PyObject *pyvx_capi = Py_None;
  double *pyvy = NULL;
  npy_intp pyvy_Dims[1] = {-1};
  const int pyvy_Rank = 1;
  PyArrayObject *capi_pyvy_tmp = NULL;
  int capi_pyvy_intent = 0;
  PyObject *pyvy_capi = Py_None;
  int *pyvedge = NULL;
  npy_intp pyvedge_Dims[2] = {-1, -1};
  const int pyvedge_Rank = 2;
  PyArrayObject *capi_pyvedge_tmp = NULL;
  int capi_pyvedge_intent = 0;
  PyObject *pyvedge_capi = Py_None;
  double *pyvarea = NULL;
  npy_intp pyvarea_Dims[1] = {-1};
  const int pyvarea_Rank = 1;
  PyArrayObject *capi_pyvarea_tmp = NULL;
  int capi_pyvarea_intent = 0;
  int *pyngbs = NULL;
  npy_intp pyngbs_Dims[2] = {-1, -1};
  const int pyngbs_Rank = 2;
  PyArrayObject *capi_pyngbs_tmp = NULL;
  int capi_pyngbs_intent = 0;
  double *pyvlenght = NULL;
  npy_intp pyvlenght_Dims[2] = {-1, -1};
  const int pyvlenght_Rank = 2;
  PyArrayObject *capi_pyvlenght_tmp = NULL;
  int capi_pyvlenght_intent = 0;
  double *pydlenght = NULL;
  npy_intp pydlenght_Dims[2] = {-1, -1};
  const int pydlenght_Rank = 2;
  PyArrayObject *capi_pydlenght_tmp = NULL;
  int capi_pydlenght_intent = 0;
  int pymaxngbh = 0;
  int pygnodes = 0;
  PyObject *pygnodes_capi = Py_None;
  int pyonodes = 0;
  PyObject *pyonodes_capi = Py_None;
  int pydnodes = 0;
  PyObject *pydnodes_capi = Py_None;
  int pydedges = 0;
  PyObject *pydedges_capi = Py_None;
  int pydelems = 0;
  PyObject *pydelems_capi = Py_None;
  int pyvnodes = 0;
  PyObject *pyvnodes_capi = Py_None;
  int pyvedges = 0;
  PyObject *pyvedges_capi = Py_None;
  static char *capi_kwlist[] = {"pyoids","pygids","pytx","pyty","pytedge","pytelmt","pyvx","pyvy","pyvedge","pygnodes","pyonodes","pydnodes","pydedges","pydelems","pyvnodes","pyvedges",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOOO|OOOOOOO:fvframe.build",\
    capi_kwlist,&pyoids_capi,&pygids_capi,&pytx_capi,&pyty_capi,&pytedge_capi,&pytelmt_capi,&pyvx_capi,&pyvy_capi,&pyvedge_capi,&pygnodes_capi,&pyonodes_capi,&pydnodes_capi,&pydedges_capi,&pydelems_capi,&pyvnodes_capi,&pyvedges_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable pyoids */
  ;
  capi_pyoids_intent |= F2PY_INTENT_IN;
  capi_pyoids_tmp = array_from_pyobj(NPY_INT,pyoids_Dims,pyoids_Rank,capi_pyoids_intent,pyoids_capi);
  if (capi_pyoids_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 1st argument `pyoids' of fvframe.build to C/Fortran array" );
  } else {
    pyoids = (int *)(PyArray_DATA(capi_pyoids_tmp));

  /* Processing variable pygids */
  ;
  capi_pygids_intent |= F2PY_INTENT_IN;
  capi_pygids_tmp = array_from_pyobj(NPY_INT,pygids_Dims,pygids_Rank,capi_pygids_intent,pygids_capi);
  if (capi_pygids_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 2nd argument `pygids' of fvframe.build to C/Fortran array" );
  } else {
    pygids = (int *)(PyArray_DATA(capi_pygids_tmp));

  /* Processing variable pytx */
  ;
  capi_pytx_intent |= F2PY_INTENT_IN;
  capi_pytx_tmp = array_from_pyobj(NPY_DOUBLE,pytx_Dims,pytx_Rank,capi_pytx_intent,pytx_capi);
  if (capi_pytx_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 3rd argument `pytx' of fvframe.build to C/Fortran array" );
  } else {
    pytx = (double *)(PyArray_DATA(capi_pytx_tmp));

  /* Processing variable pytedge */
  pytedge_Dims[1]=2;
  capi_pytedge_intent |= F2PY_INTENT_IN;
  capi_pytedge_tmp = array_from_pyobj(NPY_INT,pytedge_Dims,pytedge_Rank,capi_pytedge_intent,pytedge_capi);
  if (capi_pytedge_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 5th argument `pytedge' of fvframe.build to C/Fortran array" );
  } else {
    pytedge = (int *)(PyArray_DATA(capi_pytedge_tmp));

  /* Processing variable pytelmt */
  pytelmt_Dims[1]=3;
  capi_pytelmt_intent |= F2PY_INTENT_IN;
  capi_pytelmt_tmp = array_from_pyobj(NPY_INT,pytelmt_Dims,pytelmt_Rank,capi_pytelmt_intent,pytelmt_capi);
  if (capi_pytelmt_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 6th argument `pytelmt' of fvframe.build to C/Fortran array" );
  } else {
    pytelmt = (int *)(PyArray_DATA(capi_pytelmt_tmp));

  /* Processing variable pyvx */
  ;
  capi_pyvx_intent |= F2PY_INTENT_IN;
  capi_pyvx_tmp = array_from_pyobj(NPY_DOUBLE,pyvx_Dims,pyvx_Rank,capi_pyvx_intent,pyvx_capi);
  if (capi_pyvx_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 7th argument `pyvx' of fvframe.build to C/Fortran array" );
  } else {
    pyvx = (double *)(PyArray_DATA(capi_pyvx_tmp));

  /* Processing variable pyvedge */
  pyvedge_Dims[1]=2;
  capi_pyvedge_intent |= F2PY_INTENT_IN;
  capi_pyvedge_tmp = array_from_pyobj(NPY_INT,pyvedge_Dims,pyvedge_Rank,capi_pyvedge_intent,pyvedge_capi);
  if (capi_pyvedge_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 9th argument `pyvedge' of fvframe.build to C/Fortran array" );
  } else {
    pyvedge = (int *)(PyArray_DATA(capi_pyvedge_tmp));

  /* Processing variable pymaxngbh */
  /* Processing variable pygnodes */
  if (pygnodes_capi == Py_None) pygnodes = len(pygids); else
    f2py_success = int_from_pyobj(&pygnodes,pygnodes_capi,"fvframe.build() 1st keyword (pygnodes) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(pygids)>=pygnodes,"len(pygids)>=pygnodes","1st keyword pygnodes","build:pygnodes=%d",pygnodes) {
  /* Processing variable pyonodes */
  if (pyonodes_capi == Py_None) pyonodes = len(pyoids); else
    f2py_success = int_from_pyobj(&pyonodes,pyonodes_capi,"fvframe.build() 2nd keyword (pyonodes) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(pyoids)>=pyonodes,"len(pyoids)>=pyonodes","2nd keyword pyonodes","build:pyonodes=%d",pyonodes) {
  /* Processing variable pydnodes */
  if (pydnodes_capi == Py_None) pydnodes = len(pytx); else
    f2py_success = int_from_pyobj(&pydnodes,pydnodes_capi,"fvframe.build() 3rd keyword (pydnodes) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(pytx)>=pydnodes,"len(pytx)>=pydnodes","3rd keyword pydnodes","build:pydnodes=%d",pydnodes) {
  /* Processing variable pydedges */
  if (pydedges_capi == Py_None) pydedges = shape(pytedge,0); else
    f2py_success = int_from_pyobj(&pydedges,pydedges_capi,"fvframe.build() 4th keyword (pydedges) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(pytedge,0)==pydedges,"shape(pytedge,0)==pydedges","4th keyword pydedges","build:pydedges=%d",pydedges) {
  /* Processing variable pydelems */
  if (pydelems_capi == Py_None) pydelems = shape(pytelmt,0); else
    f2py_success = int_from_pyobj(&pydelems,pydelems_capi,"fvframe.build() 5th keyword (pydelems) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(pytelmt,0)==pydelems,"shape(pytelmt,0)==pydelems","5th keyword pydelems","build:pydelems=%d",pydelems) {
  /* Processing variable pyvnodes */
  if (pyvnodes_capi == Py_None) pyvnodes = len(pyvx); else
    f2py_success = int_from_pyobj(&pyvnodes,pyvnodes_capi,"fvframe.build() 6th keyword (pyvnodes) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(pyvx)>=pyvnodes,"len(pyvx)>=pyvnodes","6th keyword pyvnodes","build:pyvnodes=%d",pyvnodes) {
  /* Processing variable pyvedges */
  if (pyvedges_capi == Py_None) pyvedges = shape(pyvedge,0); else
    f2py_success = int_from_pyobj(&pyvedges,pyvedges_capi,"fvframe.build() 7th keyword (pyvedges) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(pyvedge,0)==pyvedges,"shape(pyvedge,0)==pyvedges","7th keyword pyvedges","build:pyvedges=%d",pyvedges) {
  /* Processing variable pyty */
  pyty_Dims[0]=pydnodes;
  capi_pyty_intent |= F2PY_INTENT_IN;
  capi_pyty_tmp = array_from_pyobj(NPY_DOUBLE,pyty_Dims,pyty_Rank,capi_pyty_intent,pyty_capi);
  if (capi_pyty_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 4th argument `pyty' of fvframe.build to C/Fortran array" );
  } else {
    pyty = (double *)(PyArray_DATA(capi_pyty_tmp));

  /* Processing variable pyvy */
  pyvy_Dims[0]=pyvnodes;
  capi_pyvy_intent |= F2PY_INTENT_IN;
  capi_pyvy_tmp = array_from_pyobj(NPY_DOUBLE,pyvy_Dims,pyvy_Rank,capi_pyvy_intent,pyvy_capi);
  if (capi_pyvy_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting 8th argument `pyvy' of fvframe.build to C/Fortran array" );
  } else {
    pyvy = (double *)(PyArray_DATA(capi_pyvy_tmp));

  /* Processing variable pyvarea */
  pyvarea_Dims[0]=pydnodes;
  capi_pyvarea_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_pyvarea_tmp = array_from_pyobj(NPY_DOUBLE,pyvarea_Dims,pyvarea_Rank,capi_pyvarea_intent,Py_None);
  if (capi_pyvarea_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting hidden `pyvarea' of fvframe.build to C/Fortran array" );
  } else {
    pyvarea = (double *)(PyArray_DATA(capi_pyvarea_tmp));

  /* Processing variable pyngbs */
  pyngbs_Dims[0]=pydnodes,pyngbs_Dims[1]=20;
  capi_pyngbs_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_pyngbs_tmp = array_from_pyobj(NPY_INT,pyngbs_Dims,pyngbs_Rank,capi_pyngbs_intent,Py_None);
  if (capi_pyngbs_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting hidden `pyngbs' of fvframe.build to C/Fortran array" );
  } else {
    pyngbs = (int *)(PyArray_DATA(capi_pyngbs_tmp));

  /* Processing variable pyvlenght */
  pyvlenght_Dims[0]=pydnodes,pyvlenght_Dims[1]=20;
  capi_pyvlenght_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_pyvlenght_tmp = array_from_pyobj(NPY_DOUBLE,pyvlenght_Dims,pyvlenght_Rank,capi_pyvlenght_intent,Py_None);
  if (capi_pyvlenght_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting hidden `pyvlenght' of fvframe.build to C/Fortran array" );
  } else {
    pyvlenght = (double *)(PyArray_DATA(capi_pyvlenght_tmp));

  /* Processing variable pydlenght */
  pydlenght_Dims[0]=pydnodes,pydlenght_Dims[1]=20;
  capi_pydlenght_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_pydlenght_tmp = array_from_pyobj(NPY_DOUBLE,pydlenght_Dims,pydlenght_Rank,capi_pydlenght_intent,Py_None);
  if (capi_pydlenght_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(fvframe_error,"failed in converting hidden `pydlenght' of fvframe.build to C/Fortran array" );
  } else {
    pydlenght = (double *)(PyArray_DATA(capi_pydlenght_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(pyoids,pygids,pytx,pyty,pytedge,pytelmt,pyvx,pyvy,pyvedge,pyvarea,pyngbs,pyvlenght,pydlenght,&pymaxngbh,&pygnodes,&pyonodes,&pydnodes,&pydedges,&pydelems,&pyvnodes,&pyvedges);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("NNNNi",capi_pyvarea_tmp,capi_pyngbs_tmp,capi_pyvlenght_tmp,capi_pydlenght_tmp,pymaxngbh);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_pydlenght_tmp == NULL) ... else of pydlenght*/
  /* End of cleaning variable pydlenght */
  }  /*if (capi_pyvlenght_tmp == NULL) ... else of pyvlenght*/
  /* End of cleaning variable pyvlenght */
  }  /*if (capi_pyngbs_tmp == NULL) ... else of pyngbs*/
  /* End of cleaning variable pyngbs */
  }  /*if (capi_pyvarea_tmp == NULL) ... else of pyvarea*/
  /* End of cleaning variable pyvarea */
  if((PyObject *)capi_pyvy_tmp!=pyvy_capi) {
    Py_XDECREF(capi_pyvy_tmp); }
  }  /*if (capi_pyvy_tmp == NULL) ... else of pyvy*/
  /* End of cleaning variable pyvy */
  if((PyObject *)capi_pyty_tmp!=pyty_capi) {
    Py_XDECREF(capi_pyty_tmp); }
  }  /*if (capi_pyty_tmp == NULL) ... else of pyty*/
  /* End of cleaning variable pyty */
  } /*CHECKSCALAR(shape(pyvedge,0)==pyvedges)*/
  } /*if (f2py_success) of pyvedges*/
  /* End of cleaning variable pyvedges */
  } /*CHECKSCALAR(len(pyvx)>=pyvnodes)*/
  } /*if (f2py_success) of pyvnodes*/
  /* End of cleaning variable pyvnodes */
  } /*CHECKSCALAR(shape(pytelmt,0)==pydelems)*/
  } /*if (f2py_success) of pydelems*/
  /* End of cleaning variable pydelems */
  } /*CHECKSCALAR(shape(pytedge,0)==pydedges)*/
  } /*if (f2py_success) of pydedges*/
  /* End of cleaning variable pydedges */
  } /*CHECKSCALAR(len(pytx)>=pydnodes)*/
  } /*if (f2py_success) of pydnodes*/
  /* End of cleaning variable pydnodes */
  } /*CHECKSCALAR(len(pyoids)>=pyonodes)*/
  } /*if (f2py_success) of pyonodes*/
  /* End of cleaning variable pyonodes */
  } /*CHECKSCALAR(len(pygids)>=pygnodes)*/
  } /*if (f2py_success) of pygnodes*/
  /* End of cleaning variable pygnodes */
  /* End of cleaning variable pymaxngbh */
  if((PyObject *)capi_pyvedge_tmp!=pyvedge_capi) {
    Py_XDECREF(capi_pyvedge_tmp); }
  }  /*if (capi_pyvedge_tmp == NULL) ... else of pyvedge*/
  /* End of cleaning variable pyvedge */
  if((PyObject *)capi_pyvx_tmp!=pyvx_capi) {
    Py_XDECREF(capi_pyvx_tmp); }
  }  /*if (capi_pyvx_tmp == NULL) ... else of pyvx*/
  /* End of cleaning variable pyvx */
  if((PyObject *)capi_pytelmt_tmp!=pytelmt_capi) {
    Py_XDECREF(capi_pytelmt_tmp); }
  }  /*if (capi_pytelmt_tmp == NULL) ... else of pytelmt*/
  /* End of cleaning variable pytelmt */
  if((PyObject *)capi_pytedge_tmp!=pytedge_capi) {
    Py_XDECREF(capi_pytedge_tmp); }
  }  /*if (capi_pytedge_tmp == NULL) ... else of pytedge*/
  /* End of cleaning variable pytedge */
  if((PyObject *)capi_pytx_tmp!=pytx_capi) {
    Py_XDECREF(capi_pytx_tmp); }
  }  /*if (capi_pytx_tmp == NULL) ... else of pytx*/
  /* End of cleaning variable pytx */
  if((PyObject *)capi_pygids_tmp!=pygids_capi) {
    Py_XDECREF(capi_pygids_tmp); }
  }  /*if (capi_pygids_tmp == NULL) ... else of pygids*/
  /* End of cleaning variable pygids */
  if((PyObject *)capi_pyoids_tmp!=pyoids_capi) {
    Py_XDECREF(capi_pyoids_tmp); }
  }  /*if (capi_pyoids_tmp == NULL) ... else of pyoids*/
  /* End of cleaning variable pyoids */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/******************************** end of build ********************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"build",-1,{{-1}},0,(char *)F_FUNC(build,BUILD),(f2py_init_func)f2py_rout_fvframe_build,doc_f2py_rout_fvframe_build},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "fvframe",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyMODINIT_FUNC PyInit_fvframe(void) {
#else
#define RETVAL
PyMODINIT_FUNC initfvframe(void) {
#endif
  int i;
  PyObject *m,*d, *s, *tmp;
#if PY_VERSION_HEX >= 0x03000000
  m = fvframe_module = PyModule_Create(&moduledef);
#else
  m = fvframe_module = Py_InitModule("fvframe", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module fvframe (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'fvframe' is auto-generated with f2py (version:2).\nFunctions:\n"
"  pyvarea,pyngbs,pyvlenght,pydlenght,pymaxngbh = build(pyoids,pygids,pytx,pyty,pytedge,pytelmt,pyvx,pyvy,pyvedge,pygnodes=len(pygids),pyonodes=len(pyoids),pydnodes=len(pytx),pydedges=shape(pytedge,0),pydelems=shape(pytelmt,0),pyvnodes=len(pyvx),pyvedges=shape(pyvedge,0))\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  Py_DECREF(s);
  fvframe_error = PyErr_NewException ("fvframe.error", NULL, NULL);
  /*
   * Store the error object inside the dict, so that it could get deallocated.
   * (in practice, this is a module, so it likely will not and cannot.)
   */
  PyDict_SetItemString(d, "_fvframe_error", fvframe_error);
  Py_DECREF(fvframe_error);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
    tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
    PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
    Py_DECREF(tmp);
  }

/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"fvframe");
#endif
  return RETVAL;
}
#ifdef __cplusplus
}
#endif
