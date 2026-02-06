#ifndef RECONSTITUTION_H
#define RECONSTITUTION_H

#include <map>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <stack>

vector<vector<poly>> Permutations(const vector<poly>& v);

void return_implem(uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_and_max, uint32_t nb_sol, vector<poly> Y, vector<poly> l, vector<poly> l2, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size, vector<poly> ANF);

#endif