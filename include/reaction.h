/** @file reaction.c
 *  @brief Header file to define reactions.
 */
#ifndef _REACTION_H_
#include "react_table.h"

typedef struct reaction_t reaction_t;

/* Maximum number of reactants that can enter into a reaction. */
#define REACTION_MAX_IN 4
#define REACTION_MAX_OUT 6

struct reaction_t {
  /* Photoionization is an special reaction where all the other parameters
   * are ignored. Be careful to put photoionization in the correct position
   * of the reaction list.
   */
  int is_photo;

  /* Number of species in and out. */
  int nin, nout;

  int input[REACTION_MAX_IN];
  int output[REACTION_MAX_OUT];

  const char *inname[REACTION_MAX_IN];
  const char *outname[REACTION_MAX_OUT];


  void (*f) (double *in, int nin, double *out, int nout, double k, double dens, react_table *rt);

  /* For reactions where k(E) is given by a reaction table. Such reactions
   * require 'f' to be f_react_table. 'k' is ignored for these reactions.
   * 'rt' is initialized with 'NULL', while 'tablefile' contains the filename
   * with the reaction table.
   *
   * This file is read in react_add().
   */
  react_table *rt;

  char *tablefile;

  double k;

  /* This allows us to define many simultaneous reactions. */
  reaction_t *next;
};

extern reaction_t *reactions_list;
extern reaction_t *reaction_index[];

#define _REACTION_H_
#endif
