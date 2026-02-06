#ifndef PRECOMPUTATION_H
#define PRECOMPUTATION_H

#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

typedef pair<poly_quad, poly_quad [10][2]> pair_xor;
// To store pairs associated with an non linear operation such that the xor of the pair is equal to the operation

//Return the size of the set corresponding to Q1_size
uint32_t nb_bits_to_size_set (uint32_t size);

//Create Q1_size
void parse_file_and_create_set(uint32_t size, poly_quad set_op []); 

//Find the position a specific polynomial in Q1,n
int64_t find_set (poly_quad val, poly_quad set_op [], uint64_t look_start, uint64_t look_end); 

//Return the size of the map corresponding to Q2_size
uint32_t nb_bits_to_size_map_xor (uint32_t size);

//Create Q2_size
void parse_file_and_create_map_xor(uint32_t size, pair_xor map_xor []); 

//Find the position of a specific polynomial in Q2,n
int64_t find_map (poly_quad val, pair_xor v [], uint64_t look_start, uint64_t look_end);   

//Create a vector containing the polynomials in L_size
vector<poly> create_poly_lin(uint32_t size, uint32_t nb_elem);

//Create a vector containing the polynomials in Q1,size^L
vector<poly> create_set_op_l_plus_lin(uint32_t size);

//Create a vector containing the polynomials in Q2,size^L
vector<poly> create_map_xor_l_plus_lin(uint32_t size);

struct implem_deg5_4 { //Used to factor degree-5 polynomials made up only of degree-5 and degree-4 monomials, where p = (quad * q) ^ r
    poly quad;
    poly q;
    poly r;
};

//Return a vector containing the solutions associated to poly p represented by the integer ent_p
vector<implem_deg5_4> parse_file_and_create_sol(uint32_t ent_p);

#endif