/** @file tree.h
 *  @brief Structures and function definitions to handle a tree structure.
 */

/*!< These fields are included as the first fields of any structure that
 *  we want to convert into leaves of a tree. For people familiar with object
 *  orientation, this is equivalent to defining a class  @a leaf and make
 *  other structures derive from it.
 *  The idea here is that we can share the code for the management of
 *  the tree structure of the grids both for the poisson and the cdr parts.
 */

#ifndef _TREE_H_

#define LEAF_FIELDS(X)				\
  X *parent;					\
  X *next;					\
  X *first_child;				\
  int level;					\

#define LEAF_ROOT  0

/*!< Iterate over all childs of a node. Simply write
 *
 * leaf_t leaf;
 *
 * iter_childs (node, leaf) {
 *    do_something (leaf);
 * }
 */
#define iter_childs(ROOT, LEAF) for (LEAF = ROOT->first_child; LEAF;	\
				     LEAF = LEAF->next)

/*!< iter_child may not work if in the loop body we make free(LEAF)
 * (though this will go unnoticed most of the time). So when
 * freeing grids, the following macro has to be used instead:
 */
#define free_childs(ROOT, LEAF, FREE) do {				\
    grid_t *next__;							\
    for (LEAF = ROOT->first_child; LEAF; LEAF = (typeof(LEAF)) next__) { \
      next__ = ((grid_t *) LEAF)->next;					\
      FREE (LEAF);							\
    }									\
  } while(0)

/*!< This is used to put into a leaf any set of values. */
#define set_leaf(LEAF, PARENT, NEXT, FIRST_CHILD, LEVEL) do {	\
    LEAF->parent = PARENT;					\
    LEAF->next = NEXT;						\
    LEAF->first_child = FIRST_CHILD;				\
    LEAF->level = LEVEL;					\
  } while(0)

#define init_leaf(LEAF) set_leaf(LEAF, NULL, NULL, NULL, -1)

/*!< Sometimes we have to remove all of a leaf's children. After freeing
 * them with free_childs, call this one.
 */
#define set_childless(LEAF) (LEAF)->first_child = NULL

/*!< Adds a child to a given leaf. */
#define add_child(PARENT, CHILD) do {					\
    set_leaf(CHILD, PARENT, PARENT->first_child, CHILD->first_child,	\
	     PARENT->level + 1);					\
    PARENT->first_child = CHILD;					\
  } while(0)

/*!< Given a function name, creates a function with the same name plus `_r'
 *  that (tail- in the second case) recursively calls the first.
 *  Can only be applied on functions that receive as single parameter a
 *  grid of type `type'.
 */
#define mk_recursive(func, type)		\
  void func ## _r (type *grid_)			\
  {						\
    type *child_;				\
    func (grid_);				\
    iter_childs (grid_, child_) {		\
      func ## _r (child_);			\
    }						\
  }

#define mk_tail_recursive(func, type)		\
  void func ## _r (type *grid_)			\
  {						\
    type *child_;				\
    iter_childs (grid_, child_) {		\
      func ## _r (child_);			\
    }						\
    func (grid_);				\
  }

#define _TREE_H_
#endif
