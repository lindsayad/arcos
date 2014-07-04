#include "stdlib.h"
#include "stdio.h"
#include "fishpack.h"


int main()
{
  fish_hstcyl (0.0, 1.0, 256, 
	       1, NULL, NULL,
	       0.0, 1.0, 256, 
	       1, NULL, NULL,
	       0.0, 0.0, NULL, 0);
}
