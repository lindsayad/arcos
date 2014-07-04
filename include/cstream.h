/** @file cstream.h
 *  @brief Global header file
 */

#ifndef _CSTREAM_H_

#include "tree.h"
#include "assert.h"
#include "stdlib.h"

#ifdef _OPENMP
# include "omp.h"
#endif

#ifndef TRUE
# define TRUE 1
#endif 

#ifndef FALSE
# define FALSE 0
#endif 


#ifdef DEBUG_LEVEL
# ifdef _OPENMP
#  define debug(level, ...) if (DEBUG_LEVEL >= level) {	      \
    fprintf (stderr,"[%d] %s:%d: ", omp_get_thread_num(),    \
	     __FILE__, __LINE__);			      \
    fprintf (stderr,__VA_ARGS__);			      \
  }
# else /*_OPENMP*/
#  define debug(level, ...) if (DEBUG_LEVEL >= level) {	      \
    fprintf (stderr,"%s:%d: ", __FILE__, __LINE__);	      \
    fprintf (stderr,__VA_ARGS__);			      \
  }
# endif
#else
# define debug(level, ...)
# define NDEBUG  /* To supress the asserts */
#endif

#define warning(...) do{				\
    fprintf (stderr, "%s: Warning: ", invok_name);	\
    fprintf (stderr, ## __VA_ARGS__);			\
  } while(0)

#define fatal(...) do{					\
    fprintf (stderr, "%s: Fatal error: ", invok_name);	\
    fprintf (stderr, ## __VA_ARGS__);			\
    exit(-1);						\
  } while(0)
    
/* Useful to debug. */
#define show_double(VAR_) printf (#VAR_ " = %g\n", VAR_)
#define show_int(VAR_) printf (#VAR_ " = %d\n", VAR_)

/* Beware of side-effects! */
#define MYMAX(X_, Y_)  ((X_) > (Y_)? (X_): (Y_))
#define MYMIN(X_, Y_)  ((X_) < (Y_)? (X_): (Y_))
#define MAX_AT_LEVEL(X_, Y_, L_)  MYMAX(X_, (Y_) << (L_))
#define MIN_AT_LEVEL(X_, Y_, L_)  MYMIN(X_, (Y_) << (L_))
#define MAX_AT_LEVEL_WITH_SHIFT(X_, Y_, L_, S_)  \
  MYMAX(X_, ((Y_) << (L_)) + ((L_) > 0? ((S_) << (L_ - 1)): 0))
#define MIN_AT_LEVEL_WITH_SHIFT(X_, Y_, L_, S_)  \
  MYMIN(X_, ((Y_) << (L_)) + ((L_) > 0? ((S_) << (L_ - 1)): 0))

#define XCHG(X1_, X2_) do {			\
    typeof(X1_) TMP_;				\
    TMP_ = X1_;					\
    X1_ = X2_;					\
    X2_ = TMP_;					\
  } while(0)


#define SQ(X_)  ((X_) * (X_))

///** @brief Information about each program parameter. */
//typedef struct param_t param_t;
//struct param_t {
  //char *name;
  //char *desc;
  //char *type;
  //void *value;
//};

/** @brief These are the types for the global parameters. */
typedef char* string;
typedef double* doublep;

//#ifdef ALLOC_PARAMS
//# define decl_param(TYPE, NAME, DESC, DEFAULT)		\
       //TYPE NAME = DEFAULT;				\
       //param_t NAME ## _st = {				\
	 //#NAME,						\
	 //DESC,						\
	 //#TYPE,						\
	 //(void *) &NAME};
//# define decl_deprec_param(TYPE, NAME, DESC, DEFAULT)	\
       //TYPE NAME;					\
       //param_t NAME ## _st = {				\
	 //#NAME,						\
	 //DESC,						\
	 //"deprecated",					\
	 //NULL};
//#else
//# define decl_param(TYPE, NAME, DESC, DEFAULT)	\
  //extern TYPE NAME;				\
  //extern param_t NAME ## _st;
//# define decl_deprec_param(TYPE, NAME, DESC, DEFAULT)	\
  //decl_param(TYPE, NAME, DESC, DEFAULT)
//#endif

typedef double REAL;

extern double *dr, *dz, dtheta;
/*!< The grid sizes */
extern double *w2k, *wk;
/*!< See cstream.c */

#define decl_field_comp(_DIR) \
  double (*ext_e_ ## _DIR) (double r, double z, double theta)
/*!< The three components of the external field: ext_e_r, ext_e_z and
 *   ext_e_theta. */

extern decl_field_comp(r);
extern decl_field_comp(z);
extern decl_field_comp(theta);

extern char *invok_name;

extern const double twopi;
extern const double invfourpi;
extern const double invpi32;

/**********
 * misc.c *
 **********/
void *xmalloc (size_t size);
void *xrealloc (void *ptr, size_t size);
void *xcalloc (size_t count, size_t size);

#define _CSTREAM_H_
#endif
