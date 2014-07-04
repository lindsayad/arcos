/** @file species.h
/*  @brief include file describing the structs specie_t and seed_t
 */

#ifndef _SPECIES_H_

#include "cstream.h"

/*!< Number of virtual species. */
#define N_VIRTUAL_SPECIES 1

/*!< Information about each of the species. */
typedef struct species_t {
  double charge;
  double mass;   /* Zero or negative mass means immobile particles. */

  const char *name;
} species_t;

/*!< Initial seeds for species. Only used at start-up */
typedef struct seed_t {
  int species;
  const char *kind_species;
  double value;
  int type;
  const char *kind_type;
  double x0;
  double y0;
  double z0;
  double sigma_x;
  double sigma_y;
  double sigma_z;
} seed_t;

extern species_t spec_electrons;
extern species_t spec_ions;

/*!< species[i] points to a species_t structure for species number I */
extern species_t *spec_index[];

extern seed_t *seed_index[];

int no_seed, no_reactions, no_species;

/*!< These are the index of two  @a special species (that sounds strange, uh?)
 *
 * For example, electrons and ions are the species that are initialized
 * with a gaussian seed, and photo_ions those that are affected by photo-
 * ionization. Their values are set in the kinetic model file (usually
 * minimal.c, but more complex models are possible).
 */
int electrons, ions, photo_ions;

/*!< There is also a  @a virtual species: Just to make it easy to create
 *  reaction functions that depend on variables other than species densities.
 *  Now this trick is only used for @a eabs absolute value of the electric
 *  field.
 */
extern const int virt_eabs;

#define _SPECIES_H_
#endif
