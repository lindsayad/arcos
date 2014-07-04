/** @file react_table.c
 *  @brief  Handles the communication between the main code and the reaction tables.
 *
 * Rather than hardcoding a different function for each possible dependence
 * of the reaction-rate on the electric field, these rate/field dependencies
 * are written in a table. This allows reactions to be added/changed/removed
 * without recompiling the source.
 */

/* Structure of the reaction table file between -----

------
Optional comments, maximum line length of 100 chars!
START
is_log
e_min
e_step
points
underflow
overflow
value(0)
value(1)
...
value(points-1)
-----

Anything before 'START' is ignored. 'START' denotes the start of the actual data
and is case-sensitive. 'is_log' equals one if the electric field values are spaced
logarithmically. 'e_min' denotes the lowest value of the electric field for
which a data-point exists. 'e_step' is the difference in fieldstrength between
consecutive data-points and 'points' is the number of data-points. 'value(x)' are
the actual data-points. Note that anything after 'points' points is ignored.

Underflow and overflow are reaction rates that are returned when the supplied field
strength is out of the bounds of the table.

IMPORTANT: electric field values are 10-logarithmic. So instead of 'e_min' actually
denotes log10(e_min)!
IMPORTANT: electric field values must be (logarithmically) equi-distant! */

/* Example file

-----
Sample reaction table. k(E) = E. 4 points are defined at E = 10^i, i \in {-2,-1,0,1}.
So the 'e_min' value is -2, 'e_step' is 1 and 'points' is 4. Overflow (100) and 
underflow (0) values are added in case E-input values fall outside the bounds, in 
this case E < 10^-2 or E > 10^1.
START
-2.0
1
4
0.0
100.0
0.01
0.1
1.0
10.0
This bit is ignored by the program, so can be used for comments. Though if i had anything
useful to say, it would've probably been better placed at the start of the file!
-----
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "react_table.h"

/** @brief Reads a reaction rate table from file @a filename and stores the
 * table in the form of a 'react_table' at location @a r.
*/
void
react_table_read(char *filename, react_table *r)
{
  FILE *fp;
  char *comments, buffer[101];
  double e_min, e_step, underflow, overflow;
  
  int cnt;

  double values[MAX_TABLE_SIZE];

  fp = fopen(filename,"r");
  if (fp == NULL)
  {
    fprintf(stderr,"Unable to load reaction table in file: %s -- Shutting down\n",filename);
    exit(1);
  }
    
  while (strncmp(buffer,"START",5) != 0)
     fgets(buffer,100,fp);
 
  fgets(buffer,100,fp);
  r->e_min = atof(buffer);
  fgets(buffer,100,fp);
  r->e_step = atof(buffer);
  fgets(buffer,100,fp);
  r->steps = atoi(buffer) - 1;
  fgets(buffer,100,fp);
  r->underflow = atof(buffer);
  fgets(buffer,100,fp);
  r->overflow = atof(buffer);
 
  for (cnt = 0; cnt <= r->steps; cnt++)
  {
    fgets(buffer,100,fp);
    r->values[cnt] = atof(buffer);
  }

  fclose(fp);
}

/** @brief Computes an approximation for the reaction rate.
 *
 * Computes an approximation for the reaction rate at field-strength @a e
 * by interpolating reaction rates from lookup table @a r. Stores the result
 * in @a ra.
 */
void
react_table_lookup(react_table *r, double e, double *ra)
{
  int pos;
  double e1, e2, val1, val2, log_e, res;

  // e == 0 is possible at initialization.
    log_e = (e > 0) ? log10(e) : -1000.0;
   
    /* If the supplied fieldstrength falls outside the boundaries of the table,
       return predetermined under-/overflow values. In most cases the underflow(overflow)
       will be equal to the lowest(highest) of the specified values in the table. */
    if (log_e < r->e_min) { *ra = r->underflow; return; }
    if (log_e > r->e_min + r->e_step * r->steps) { *ra = r->overflow; return; }
  
    pos = floor((log_e - r->e_min) / r->e_step);
    val1 = r->values[pos];
    val2 = r->values[pos+1];
    e1 = pow(10, r->e_min + r->e_step * (double) pos);
    e2 = e1 * pow(10, r->e_step);

    /* Linear interpolation of the 2 table values closest to e. 
       Note that the field-values in the table are logarithmic! */

    res = val1 + ((val2 - val1) / (e2 - e1)) * (e - e1);

    /* Dirty tricks with pointers to get the return value right. The normal 'return'
       statement would give a completely different value on the "other side".

       It's not pretty, but it works, so who cares? :) */
    *ra = res; 
}  
