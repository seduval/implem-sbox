#ifndef SBOX_H
#define SBOX_H

#include <cstdint>

#define ABS(x) ((x)>0?(x):(-(x)))

//Compute the mobius transform from the truth_table
void moebius_transform (uint32_t * table, uint32_t nb_elem, uint32_t nb_words);

//Compute the anf from the truth_table
void truth_table_to_anf (uint32_t * truth_table, uint32_t nb_elem, uint32_t nb_words); 

//Compute the truth_table from the anf
void anf_to_truth_table (uint32_t * anf, uint32_t nb_elem, uint32_t nb_words); 

//Compute the truth_table from the LUT
void lut_to_truth_table (uint32_t ** truth_table, uint32_t * lut, uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_words); 

//Compute the ANF from the LUT
void lut_to_anf (uint32_t ** anf, uint32_t * lut, uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_words);

//Compute the hamming weight of n
uint32_t weight (uint32_t n); 

//Compute the algebraic degree of one component of the anf
uint32_t algebraic_degree_of_anf (uint32_t * anf, uint32_t nb_elem, uint32_t nb_words); 

//Compute the DDT from the LUT
void DDT(uint32_t lut[], uint32_t **ddt, uint32_t nb_elem); 

//Compute the differential uniformity from the DDT
uint32_t Un_diff(uint32_t **tab, uint32_t nb_elem); 

//Compute the walsh transform
void walsh_transform_bool (uint32_t * f, int32_t * walsh_of_f, int size); 

//Compute the linearity from the LUT
uint32_t linearity (uint32_t * lut, int size); 

//Compute the inverse 
void inv(uint32_t lut [], uint32_t lut_inv [], uint32_t nb_elem);

//Print the anf
void print_anf (uint32_t ** anf, uint32_t size_in, uint32_t size_out);

#endif