/* Some utility functions (not part of the algorithm) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to return the maximum of two variables */
float maximum (float a, float b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
float minimum (float a, float b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}
