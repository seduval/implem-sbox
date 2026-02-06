#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstring>
#include <cstdint>

extern "C"
{
    #include <inttypes.h>
}

#include "sbox.h"

using namespace std;


//Compute the mobius transform from the truth_table
void moebius_transform (uint32_t * table, uint32_t nb_elem, uint32_t nb_words) {

    // For higher order bits.
    for (uint32_t k=1; k<nb_words; k<<=1) {
        for (uint32_t i=0; i<nb_words; i+=k<<1) {
            for (uint32_t j=0; j<k; j++) {
                table[i+k+j] ^= table[i+j];
            }
        }
    }

    // For lower order bits.
    for(uint32_t word_n=0;word_n<nb_words;word_n++) {
        uint32_t tmp = table[word_n];
        tmp^=(tmp&0x55555555)<<1;
        tmp^=(tmp&0x33333333)<<2;
        tmp^=(tmp&0xf0f0f0f)<<4;
        tmp^=(tmp&0xff00ff)<<8;
        tmp^=(tmp<<16);
        table[word_n]=tmp;
    }

    if (nb_elem&0x1f)
        table[nb_words-1] &= (1ul << (nb_elem&0x1f)) - 1;
}

//Compute the anf from the truth_table
void truth_table_to_anf (uint32_t * truth_table, uint32_t nb_elem, uint32_t nb_words) {
    moebius_transform(truth_table, nb_elem, nb_words);
}

//Compute the truth_table from the anf
void anf_to_truth_table (uint32_t * anf, uint32_t nb_elem, uint32_t nb_words) {
    moebius_transform(anf, nb_elem, nb_words);
}

//Compute the truth_table from the LUT
void lut_to_truth_table (uint32_t ** truth_table, uint32_t * lut, uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_words) {
    uint32_t x, fx, bit_n, fx_shift, x_word, x_ind;
    uint32_t word_n;

    for (bit_n=0; bit_n<size_out; bit_n++)
        for (word_n=0; word_n<nb_words; word_n++)
            truth_table[bit_n][word_n] = 0;

    for (x=0; x<nb_elem; x++) {
        fx = lut[x];
        x_word = x>>5; // x/32
        x_ind = x&0x1f; // x%32
        for (bit_n=0, fx_shift=fx; bit_n<size_out; bit_n++, fx_shift>>=1) {
            truth_table[bit_n][x_word] |= (fx_shift&1ul)<<x_ind;
            //cout<<" bit_n : "<<bit_n<<" ; x_word : "<<x_word<<" ; valeur : "<<truth_table[bit_n][x_word]<<endl;
        }
    }
}

//Compute the ANF from the LUT
void lut_to_anf (uint32_t ** anf, uint32_t * lut, uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_words) {
    uint32_t ** truth_table = (uint32_t**)malloc((size_out)*sizeof(uint32_t*));
    for (uint32_t comp_n=0; comp_n<size_out; comp_n++)
        truth_table[comp_n] = (uint32_t*)malloc(nb_words*sizeof(uint32_t));

    lut_to_truth_table((uint32_t**)truth_table, lut, size_in, size_out, nb_elem, nb_words);

    for (uint32_t coord_n=0; coord_n<size_out; coord_n++) {
        memcpy(anf[coord_n], truth_table[coord_n], nb_elem/8);
        truth_table_to_anf(anf[coord_n], nb_elem, nb_words);
    }

    for (uint32_t comp_n=0; comp_n<size_out; comp_n++)
        free(truth_table[comp_n]);
    free(truth_table);
}

//Compute the hamming weight of n
uint32_t weight (uint32_t n) {
    
    // each bit in n is a one-bit integer that indicates how many bits are set
    // in that bit.
    
    n = ((n & 0xAAAAAAAA) >> 1) + (n & 0x55555555);
    // Now every two bits are a two bit integer that indicate how many bits were
    // set in those two bits in the original number
    
    n = ((n & 0xCCCCCCCC) >> 2) + (n & 0x33333333);
    // Now we're at 4 bits
    
    n = ((n & 0xF0F0F0F0) >> 4) + (n & 0x0F0F0F0F);
    // 8 bits
    
    n = ((n & 0xFF00FF00) >> 8) + (n & 0x00FF00FF);
    // 16 bits
    
    n = ((n & 0xFFFF0000) >> 16) + (n & 0x0000FFFF);
    // kaboom - 32 bits
    
    return n;
}

//Compute the algebraic degree of one component of the anf
uint32_t algebraic_degree_of_anf (uint32_t * anf, uint32_t nb_elem, uint32_t nb_words) {
    
    uint32_t max_deg = 0;
    for (uint32_t u=0; u<nb_elem; u++) {
        if ((anf[u>>5]>>(u&0x1f))&1ul) {
            uint32_t w = weight(u);
            if (w > max_deg) {
                max_deg = w;
            }
        }
    }
    
    return max_deg;
}

//Compute the DDT from the LUT
void DDT(uint32_t lut[], uint32_t **ddt, uint32_t nb_elem) {          
    for (uint32_t a=0; a<nb_elem; a++) {           
        for (uint32_t x=0; x<nb_elem; x++) {
                int b = lut[x^a] ^ lut[x];
                ddt[a][b] += 1;
        }
    }
}

//Compute the differential uniformity from the DDT
uint32_t Un_diff(uint32_t **tab, uint32_t nb_elem) {                              
    for(uint32_t i=1; i <nb_elem; i++)   {
        for(uint32_t j=0; j<nb_elem; j++)   {
            if (tab[1][1] < tab [i][j])
                tab[1][1] = tab[i][j];
        }
    }
    return tab[1][1];
}

//Compute the walsh transform
void walsh_transform_bool (uint32_t * f, int32_t * walsh_of_f, int size) {
    uint32_t nb_elem = 1ul<<size;
    uint32_t nb_elem_k, nb_elem_n_k, pow_2_k, block_start, block_shifted_start, ind1, ind2, i, j;
    int memsize = nb_elem*sizeof(uint32_t);
    int32_t * buff = (int32_t*)malloc (memsize);

    for (uint32_t x=0; x<nb_elem; x++) {
    walsh_of_f[x] = (f[x]?-1:1);
    }

    // For all blocks of size 2^k
    for (int k=1; k<=size; k++) {
    nb_elem_k = 1ul<<(k-1);
    pow_2_k = nb_elem_k<<1;
    nb_elem_n_k = 1ul<<(size-k);
    for (i=0; i<nb_elem_n_k; i++) {
        block_start = i*pow_2_k;
        block_shifted_start = block_start + nb_elem_k;
        // Compute the image of the i-th 2^k-bit block.
        for (j=0; j<nb_elem_k; j++) {
            ind1 = block_start+j;
            ind2 = block_shifted_start+j;
            buff[ind1] = walsh_of_f[ind1] + walsh_of_f[ind2];
            buff[ind2] = walsh_of_f[ind1] - walsh_of_f[ind2];
        }
    }
    memcpy(walsh_of_f, buff, memsize);
    }
    free(buff);
}

//Compute the linearity from the LUT
uint32_t linearity (uint32_t * lut, int size) {

    uint32_t max = 0;

    #pragma omp parallel
    {
    uint32_t nb_elem = 1ul<<size;
    int memsize = nb_elem*sizeof(uint32_t);
    uint32_t private_max = 0;
    #pragma omp for
    for (uint32_t comp_n=1; comp_n<nb_elem; comp_n++) {
        uint32_t * component = (uint32_t*)malloc (memsize);
        int32_t * walsh_of_component = (int32_t*)malloc (memsize);

        for (uint32_t x=0; x<nb_elem; x++)
            component[x] = 0;

        uint32_t comp = comp_n;
        uint32_t coord = 0;
        while (comp > 0) {
            if (comp & 1ul) {
            for (uint32_t x=0; x<nb_elem; x++)
                component[x] ^= (lut[x]>>coord) & 1ul;
            }
            comp >>= 1;
            coord++;
        }

        walsh_transform_bool(component, walsh_of_component, size);

        free(component);

        // Updating linearity with Walsh transform of this component.
        for (uint32_t x=0; x<nb_elem; x++) {
            if (ABS(walsh_of_component[x]) > private_max) {
                private_max = ABS(walsh_of_component[x]);
            }
        }
        free(walsh_of_component);
    }

    #pragma omp critical
    {
        if (private_max > max)
            max = private_max;
        }
    }

    return max;
}

//Compute the inverse 
void inv(uint32_t lut [], uint32_t lut_inv [], uint32_t nb_elem){
    for (uint32_t i=0; i<nb_elem; i++)   {
        lut_inv[lut[i]] = i;
    }
}

//Print the ANF
void print_anf (uint32_t ** anf, uint32_t size_in, uint32_t size_out) {
    cout<<"ANF is given by : "<<endl;
    uint32_t nb_elem = 0x1<<size_in;

    bool xor_up = false;

    for (uint32_t comp_n=0; comp_n<size_out; comp_n++) {
        cout << "y" + to_string(comp_n) + " = ";
        xor_up = false;
        for (uint32_t u=0; u<nb_elem; u++) {
            if (((anf[comp_n][u>>5]>>(u&0x1f))&1ul)) {
                if (xor_up)
                    cout << " ^ ";
                if (u==0)
                    cout << "1";
                else {
                    for (uint32_t bit_n=0; bit_n<size_in; bit_n++) {
                        if ((u>>bit_n)&1ul)
                            cout << "x" + to_string(bit_n);
                    }
                }
                xor_up = true;
            }
        }
        cout << "\n";
    }
}