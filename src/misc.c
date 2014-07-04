/** @file misc.c
 *  @brief Miscellaneous utilities.
 */

#include <stdlib.h>
#include <stdio.h>

extern char *invok_name;

/** @brief  Allocates memory with check and (eventually) error reporting. */
void*
xmalloc (size_t size)
{
  void *value = malloc(size);
  if (value == 0) {
    fprintf (stderr, "%s: virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}

/** @brief Reallocates memory */
void *
xrealloc (void *ptr, size_t size)
{
  void *value = realloc (ptr, size);
  if (value == 0) {
    fprintf (stderr, "%s: Virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}


/** @brief Reallocates memory, but now initializes the memory to zero. */
void*
xcalloc (size_t count, size_t size)
{
  void *value;

  value = calloc (count, size);
  if (value == 0){
    fprintf (stderr, "%s: virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}
