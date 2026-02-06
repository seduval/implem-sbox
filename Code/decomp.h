#ifndef DECOMP_H
#define DECOMP_H

struct implem{
    vector<poly> op_sol; /* Contain all the polynomials to implement an output bit */ 
    set<poly_quad> quad_sol; /* Contain only the quadratic polynomials to implement an output bit */
    string formula; /* Formula used to compute the implementation from the polynomials. */
};

// Remove duplicates in a vector of implem
vector<implem> remove_duplicates(vector<implem> v);

// Return the most significant bit in a quadratic polynomial
uint64_t most_significant_bit(poly_quad p);

// Compute the rank of a family of quadratic polynomials
uint32_t Rank(set<poly_quad> op_selec);

typedef vector<implem> (*fptr)(poly, vector<poly>, vector<poly>, uint32_t, uint32_t, uint32_t, poly_quad [], pair_xor [], uint32_t, uint32_t );

class decomposition {
public:

    fptr ftab[TAB_SIZE];

    decomposition(uint32_t size, uint32_t deg);

    vector<implem> add_to_op_selec(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size, uint32_t * test);
};

vector<vector<implem>> create_sets(uint32_t * a, uint32_t * test, vector<vector<uint32_t>> * T, vector<poly> y, vector<poly> l, vector<poly> set_op_l_plus_lin , uint32_t size_in, uint32_t size_out, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size);

#endif