/** @file tree.c
 *  @brief Functions to handle tree structures.
 */

#include <stdlib.h>
#include <stdio.h>

#include "species.h"

#include "tree.h"

/** @brief Allocates new memory for a leaf and gives it the corresponding
 * values.
 */
leaf_t*
leaf_new_a (leaf_t *parent, leaf_t *next, leaf_t *first_child, void *data)
{
  leaf_t *leaf;

  debug (2, "leaf_new_a\n");

  leaf = (leaf_t*) xmalloc (sizeof(leaf_t));

  leaf->parent = parent;
  leaf->next = next;
  leaf->first_child = first_child;
  leaf->data = data;

  return leaf;
}

/** @brief Creates a new child of @a parent with the given @a data. */
leaf_t*
leaf_new_child_a (leaf_t *parent, void *data)
{
  leaf_t *leaf;

  debug (2, "leaf_new_child_a\n");

  leaf = leaf_new_a (parent, parent->first_child, NULL, data);

  parent->first_child = leaf;

  leaf->level = parent->level + 1;

  return leaf;
}
