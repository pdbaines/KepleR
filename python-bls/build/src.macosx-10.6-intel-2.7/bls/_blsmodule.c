/* File: _blsmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * See http://cens.ioc.ee/projects/f2py2e/
 * Generation date: Tue Oct 29 22:08:41 2013
 * $Revision:$
 * $Date:$
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
static PyObject *_bls_error;
static PyObject *_bls_module;

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
#define old_rank(var) (((PyArrayObject *)(capi_ ## var ## _tmp))->nd)
#define old_shape(var,dim) (((PyArrayObject *)(capi_ ## var ## _tmp))->dimensions[dim])
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
    PyErr_SetString(_bls_error,errstring);\
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

static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
  PyObject* tmp = NULL;
  if (PyFloat_Check(obj)) {
#ifdef __sgi
    *v = PyFloat_AsDouble(obj);
#else
    *v = PyFloat_AS_DOUBLE(obj);
#endif
    return 1;
  }
  tmp = PyNumber_Float(obj);
  if (tmp) {
#ifdef __sgi
    *v = PyFloat_AsDouble(tmp);
#else
    *v = PyFloat_AS_DOUBLE(tmp);
#endif
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
    if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
    Py_DECREF(tmp);
  }
  {
    PyObject* err = PyErr_Occurred();
    if (err==NULL) err = _bls_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
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
    if (err==NULL) err = _bls_error;
    PyErr_SetString(err,errmess);
  }
  return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(eebls,EEBLS)(int*,double*,double*,double*,double*,int*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,int*,int*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/*********************************** eebls ***********************************/
static char doc_f2py_rout__bls_eebls[] = "\
Function signature:\n\
  p,bper,bpow,depth,qtran,in1,in2 = eebls(t,x,u,v,nf,fmin,df,nb,qmi,qma,[n])\n\
Required arguments:\n"
"  t : input rank-1 array('d') with bounds (n)\n"
"  x : input rank-1 array('d') with bounds (n)\n"
"  u : in/output rank-1 array('d') with bounds (n)\n"
"  v : in/output rank-1 array('d') with bounds (n)\n"
"  nf : input int\n"
"  fmin : input float\n"
"  df : input float\n"
"  nb : input int\n"
"  qmi : input float\n"
"  qma : input float\n"
"Optional arguments:\n"
"  n := len(t) input int\n"
"Return objects:\n"
"  p : rank-1 array('d') with bounds (nf)\n"
"  bper : float\n"
"  bpow : float\n"
"  depth : float\n"
"  qtran : float\n"
"  in1 : int\n"
"  in2 : int";
/* extern void F_FUNC(eebls,EEBLS)(int*,double*,double*,double*,double*,int*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,int*,int*); */
static PyObject *f2py_rout__bls_eebls(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,double*,double*,double*,double*,int*,double*,double*,int*,double*,double*,double*,double*,double*,double*,double*,int*,int*)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  int n = 0;
  PyObject *n_capi = Py_None;
  double *t = NULL;
  npy_intp t_Dims[1] = {-1};
  const int t_Rank = 1;
  PyArrayObject *capi_t_tmp = NULL;
  int capi_t_intent = 0;
  PyObject *t_capi = Py_None;
  double *x = NULL;
  npy_intp x_Dims[1] = {-1};
  const int x_Rank = 1;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  PyObject *x_capi = Py_None;
  double *u = NULL;
  npy_intp u_Dims[1] = {-1};
  const int u_Rank = 1;
  PyArrayObject *capi_u_tmp = NULL;
  int capi_u_intent = 0;
  PyObject *u_capi = Py_None;
  double *v = NULL;
  npy_intp v_Dims[1] = {-1};
  const int v_Rank = 1;
  PyArrayObject *capi_v_tmp = NULL;
  int capi_v_intent = 0;
  PyObject *v_capi = Py_None;
  int nf = 0;
  PyObject *nf_capi = Py_None;
  double fmin = 0;
  PyObject *fmin_capi = Py_None;
  double df = 0;
  PyObject *df_capi = Py_None;
  int nb = 0;
  PyObject *nb_capi = Py_None;
  double qmi = 0;
  PyObject *qmi_capi = Py_None;
  double qma = 0;
  PyObject *qma_capi = Py_None;
  double *p = NULL;
  npy_intp p_Dims[1] = {-1};
  const int p_Rank = 1;
  PyArrayObject *capi_p_tmp = NULL;
  int capi_p_intent = 0;
  double bper = 0;
  double bpow = 0;
  double depth = 0;
  double qtran = 0;
  int in1 = 0;
  int in2 = 0;
  static char *capi_kwlist[] = {"t","x","u","v","nf","fmin","df","nb","qmi","qma","n",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOOOOOOO|O:_bls.eebls",\
    capi_kwlist,&t_capi,&x_capi,&u_capi,&v_capi,&nf_capi,&fmin_capi,&df_capi,&nb_capi,&qmi_capi,&qma_capi,&n_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable bper */
  /* Processing variable df */
    f2py_success = double_from_pyobj(&df,df_capi,"_bls.eebls() 7th argument (df) can't be converted to double");
  if (f2py_success) {
  /* Processing variable in1 */
  /* Processing variable in2 */
  /* Processing variable qmi */
    f2py_success = double_from_pyobj(&qmi,qmi_capi,"_bls.eebls() 9th argument (qmi) can't be converted to double");
  if (f2py_success) {
  /* Processing variable qma */
    f2py_success = double_from_pyobj(&qma,qma_capi,"_bls.eebls() 10th argument (qma) can't be converted to double");
  if (f2py_success) {
  /* Processing variable nb */
    f2py_success = int_from_pyobj(&nb,nb_capi,"_bls.eebls() 8th argument (nb) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nf */
    f2py_success = int_from_pyobj(&nf,nf_capi,"_bls.eebls() 5th argument (nf) can't be converted to int");
  if (f2py_success) {
  /* Processing variable depth */
  /* Processing variable bpow */
  /* Processing variable t */
  ;
  capi_t_intent |= F2PY_INTENT_IN;
  capi_t_tmp = array_from_pyobj(NPY_DOUBLE,t_Dims,t_Rank,capi_t_intent,t_capi);
  if (capi_t_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(_bls_error,"failed in converting 1st argument `t' of _bls.eebls to C/Fortran array" );
  } else {
    t = (double *)(capi_t_tmp->data);

  /* Processing variable fmin */
    f2py_success = double_from_pyobj(&fmin,fmin_capi,"_bls.eebls() 6th argument (fmin) can't be converted to double");
  if (f2py_success) {
  /* Processing variable qtran */
  /* Processing variable n */
  if (n_capi == Py_None) n = len(t); else
    f2py_success = int_from_pyobj(&n,n_capi,"_bls.eebls() 1st keyword (n) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(len(t)>=n,"len(t)>=n","1st keyword n","eebls:n=%d",n) {
  /* Processing variable p */
  p_Dims[0]=nf;
  capi_p_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_p_tmp = array_from_pyobj(NPY_DOUBLE,p_Dims,p_Rank,capi_p_intent,Py_None);
  if (capi_p_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(_bls_error,"failed in converting hidden `p' of _bls.eebls to C/Fortran array" );
  } else {
    p = (double *)(capi_p_tmp->data);

  /* Processing variable u */
  u_Dims[0]=n;
  capi_u_intent |= F2PY_INTENT_INOUT;
  capi_u_tmp = array_from_pyobj(NPY_DOUBLE,u_Dims,u_Rank,capi_u_intent,u_capi);
  if (capi_u_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(_bls_error,"failed in converting 3rd argument `u' of _bls.eebls to C/Fortran array" );
  } else {
    u = (double *)(capi_u_tmp->data);

  /* Processing variable v */
  v_Dims[0]=n;
  capi_v_intent |= F2PY_INTENT_INOUT;
  capi_v_tmp = array_from_pyobj(NPY_DOUBLE,v_Dims,v_Rank,capi_v_intent,v_capi);
  if (capi_v_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(_bls_error,"failed in converting 4th argument `v' of _bls.eebls to C/Fortran array" );
  } else {
    v = (double *)(capi_v_tmp->data);

  /* Processing variable x */
  x_Dims[0]=n;
  capi_x_intent |= F2PY_INTENT_IN;
  capi_x_tmp = array_from_pyobj(NPY_DOUBLE,x_Dims,x_Rank,capi_x_intent,x_capi);
  if (capi_x_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(_bls_error,"failed in converting 2nd argument `x' of _bls.eebls to C/Fortran array" );
  } else {
    x = (double *)(capi_x_tmp->data);

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(&n,t,x,u,v,&nf,&fmin,&df,&nb,&qmi,&qma,p,&bper,&bpow,&depth,&qtran,&in1,&in2);
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
    capi_buildvalue = Py_BuildValue("Nddddii",capi_p_tmp,bper,bpow,depth,qtran,in1,in2);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  if((PyObject *)capi_x_tmp!=x_capi) {
    Py_XDECREF(capi_x_tmp); }
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
  if((PyObject *)capi_v_tmp!=v_capi) {
    Py_XDECREF(capi_v_tmp); }
  }  /*if (capi_v_tmp == NULL) ... else of v*/
  /* End of cleaning variable v */
  if((PyObject *)capi_u_tmp!=u_capi) {
    Py_XDECREF(capi_u_tmp); }
  }  /*if (capi_u_tmp == NULL) ... else of u*/
  /* End of cleaning variable u */
  }  /*if (capi_p_tmp == NULL) ... else of p*/
  /* End of cleaning variable p */
  } /*CHECKSCALAR(len(t)>=n)*/
  } /*if (f2py_success) of n*/
  /* End of cleaning variable n */
  /* End of cleaning variable qtran */
  } /*if (f2py_success) of fmin*/
  /* End of cleaning variable fmin */
  if((PyObject *)capi_t_tmp!=t_capi) {
    Py_XDECREF(capi_t_tmp); }
  }  /*if (capi_t_tmp == NULL) ... else of t*/
  /* End of cleaning variable t */
  /* End of cleaning variable bpow */
  /* End of cleaning variable depth */
  } /*if (f2py_success) of nf*/
  /* End of cleaning variable nf */
  } /*if (f2py_success) of nb*/
  /* End of cleaning variable nb */
  } /*if (f2py_success) of qma*/
  /* End of cleaning variable qma */
  } /*if (f2py_success) of qmi*/
  /* End of cleaning variable qmi */
  /* End of cleaning variable in2 */
  /* End of cleaning variable in1 */
  } /*if (f2py_success) of df*/
  /* End of cleaning variable df */
  /* End of cleaning variable bper */
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
/******************************** end of eebls ********************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"eebls",-1,{{-1}},0,(char *)F_FUNC(eebls,EEBLS),(f2py_init_func)f2py_rout__bls_eebls,doc_f2py_rout__bls_eebls},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "_bls",
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
PyObject *PyInit__bls(void) {
#else
#define RETVAL
PyMODINIT_FUNC init_bls(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = _bls_module = PyModule_Create(&moduledef);
#else
  m = _bls_module = Py_InitModule("_bls", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module _bls (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module '_bls' is auto-generated with f2py (version:2).\nFunctions:\n"
"  p,bper,bpow,depth,qtran,in1,in2 = eebls(t,x,u,v,nf,fmin,df,nb,qmi,qma,n=len(t))\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  _bls_error = PyErr_NewException ("_bls.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));

/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"_bls");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
