#ifndef POLY_H
#define POLY_H

#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <bitset>
#include <set>
#include <omp.h>
#include <random>
#include <getopt.h>
#include <ostream>
#include <cstdint>

using namespace std;

#define POLY_ADD(a, b, apb) \
    { \
        apb.data = (a.data) ^ (b.data); \
    } \
    
//Add two poly and put it in apb

#define POLY_EQUAL(a, b) \
    ({ \
        bool equal = true; \
        for (uint32_t ii=0; ii< MAX; ii++) { \
            if ((a).get(ii) != (b).get(ii)) { \
                equal = false; \
                break; \
            } \
        } \
        equal; \
    })
//Test the egality of two poly

#define MONOM_MUL(a, b, size) \
    ({ \
        uint32_t temp_a = (a); \
        uint32_t bits_a[size]; \
        for (uint32_t ii = 0; ii < (size); ii++) { \
            bits_a[ii] = temp_a & 1ul; \
            temp_a >>= 1; \
        } \
        uint32_t temp_b = (b); \
        uint32_t bits_b[size]; \
        for (uint32_t ii = 0; ii < (size); ii++) { \
            bits_b[ii] = temp_b & 1ul; \
            temp_b >>= 1; \
        } \
        for (uint32_t jj = 0; jj < (size); jj++) { \
            if ((bits_a[jj] == bits_b[jj]) && (bits_a[jj] == 1)) { \
                (a) ^= (1u << jj); \
            } \
        } \
        (a) + (b); \
    })
//Used in POLY_MUL

uint32_t monom_mul(uint32_t a, uint32_t b , uint32_t size);

/**Multiply two poly and put it in ab**/
#define POLY_MUL(a, b, ab, nb_elem) \
    ({ \
        for (uint32_t ii = 0; ii < (nb_elem); ii++) { \
            if (((a).get(ii >> 5) >> (ii & 0x1f)) & 1u){ \
                for (uint32_t jj = 0; jj < (nb_elem); jj++) { \
                    if ((((b).get(jj >> 5) >> (jj & 0x1f)) & 1u)) { \
                        uint32_t val_mon = (ii)|(jj); \
                        uint32_t k = val_mon & 0x1f; \
                        uint32_t l = val_mon >> 5; \
                        (ab).data[l] ^= 1u << k; \
                    } \
                } \
            }\
        } \
    })

/**Copy a in b**/
#define COPY_POLY(a,b) \
    ({ \
        for (uint32_t ii=0; ii<MAX; ii++) { \
           (b).data[ii] = (a).get(ii); \
        } \
    })

/**Multiply a poly p and a monomial mon and put it in res**/
#define POLY_MUL2(mon, p, res, nb_elem) \
    ({ \
        for (uint32_t jj = 0; jj < (nb_elem); jj++) { \
            if ((((p).get(jj >> 5) >> (jj & 0x1f)) & 1u)) { \
                uint32_t val_mon = (mon)|(jj); \
                (res).data[(val_mon >> 5)] ^= 1u << (val_mon & 0x1f); \
            } \
        } \
    })


#define MONOM_DIV(a,b) (a)^(b)

typedef uint64_t poly_quad;

//Give the correspondance between the place of a quadratic monomial in poly_quad and its place in a poly 
uint32_t pos_to_power (uint32_t pos, uint32_t size);

class poly {
public:

    uint32_t data __attribute__ ((vector_size (MAX*4))); // 

    // Constructor initializing the poly to zero
    poly();

    // Destructor
    ~poly();

    // Constructor from a char*
    // cpol must look like "x41^x20^x1^x17\0" ( = x0x3x5 + x2x4 + x0 + x0x4) over 6 bits
    poly(const char * cpol);

    //Constructor from a poly_quad
    poly(poly_quad p, uint32_t size);

    bool operator<(const poly& p) const;

    bool operator==(const poly& other) const;

    bool operator!=(const poly& other) const;

    // Test if the poly is linear
    bool is_linear_monomial(uint32_t size);

    //Return the i-element of a poly
    uint32_t get(uint32_t i) const;

    //Print the poly
    void print_poly(uint32_t size) const;

    string poly_to_string(uint32_t size) const;

    //Add the current poly to an other
    void add(poly a);

    //Return the degree of the current poly
    uint32_t algebraic_degree (uint32_t nb_elem);

};

namespace std {
    template<>
    struct hash<poly> {
        std::size_t operator()(const poly& p) const noexcept {
            std::size_t h = 0;
            for (size_t i = 0; i < MAX; i++) {
                h ^= std::hash<uint32_t>()(p.data[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
}

//Convert a poly in a poly_quad
poly_quad poly_to_poly_quad (poly op, uint32_t size, uint32_t nb_elem);

//Truncate a poly from its linear part
poly truncate_lin(const poly &a, uint32_t size);

//Conserve only monomials of degree d in a poly
poly truncate_except_deg(const poly& a, uint32_t d, uint32_t size, uint32_t nb_elem);

//Create the monomial order associated to the given monomial monom_dom
void create_monomial_order(uint32_t monomial_order [] , uint32_t nb_elem, uint32_t size, uint32_t monom_dom);

//return the leading term of a poly considering the given monomial order
uint32_t leading_term(const poly& a, uint32_t monomial_order [], uint32_t nb_elem); 

//Test if two monomials are divisible
bool is_divisible(uint32_t a, uint32_t b);   

//Create a poly containing the specified monomial(7=111="x2x1x0")
poly mon_to_poly(uint32_t mon);

// Return a couple (q,r) as a = bq + r with deg(r) < deg(a) considering the given monomial order
pair<poly, poly> poly_div (poly a, poly b, uint32_t monomial_order [], uint32_t size, uint32_t nb_elem); 

pair<poly, poly> poly_div_no_truncate (const poly& a, const poly& b, uint32_t monomial_order [], uint32_t size, uint32_t nb_elem);

//Return if we need to create the monomial order associate to mon_dom
bool mon_dom_in (poly l, uint32_t mon_dom, uint32_t size, uint32_t nb_elem);

//Return a unique integer associated to a degree-5-and-4 poly
uint32_t poly_deg5_4_to_uint(poly p);

#endif