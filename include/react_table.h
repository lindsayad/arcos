/** @file  react_table.h
 *  @brief Maximum size of lookup tables.
 */
#ifndef _REACT_TABLE_H_
#define MAX_TABLE_SIZE 500

typedef struct react_table react_table;

struct react_table
{
  double e_min;
  double e_step;

  double underflow;
  double overflow;

  int steps;
  
  double values[MAX_TABLE_SIZE];
};

#define _REACT_TABLE_H_
#endif
