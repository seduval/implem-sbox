#ifndef CONFIG_H
#define CONFIG_H

#define MAX 4
/* MAX represents the number of words of 32 bits we need to express a poly
For example, over 6 bits, the number of monomials is 2^6 = 64, we will then have MAX = 2 */

#define TAB_SIZE 2

//Degré 3
#define MAX_RANK_DEG3_V1 3
#define MAX_RANK_DEG3_V2 4

//Degré 4

#define MAX_RANK_DEG4_V1 4
#define MAX_RANK_DEG4_V2 7

//Degré 5

#define MAX_RANK_DEG5_V1 8
#define MAX_RANK_DEG5_V2 6

#define verbose 0

#endif
