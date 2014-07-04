/** @file photo.h
 *  @brief Auxiliary header file for photo.c.
 */
#ifndef _PHOTO_H_

typedef struct photo_term_t photo_term_t;

struct photo_term_t
{
  double A;
  double lambda;

  /* We store the photoionization terms in a linked list. */
  photo_term_t *next;
};

/* Photoionization terms. */
extern photo_term_t *photo_terms;

#define photo_printf_str "{A = %g, lambda = %g}"
#define photo_printf_args(_T)  (_T)->A, (_T)->lambda

#define _PHOTO_H_
#endif
