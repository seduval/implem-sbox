#include "config.h"
#include "sbox.h"
#include "poly.h"

using namespace std;

// Constructor by default, initializing the poly to 0
poly::poly() {
    for (uint32_t i=0; i<MAX; i++){
        data[i] = 0;
    }
}

// Destructor
poly::~poly(){}

// Constructor from a char*
// cpol must look like "x41^x20^x1^x17\0" ( = x0x3x5 + x2x4 + x0 + x0x4) over 6 bits
poly::poly(const char * cpol) {
    for (uint32_t i=0; i<MAX; i++){
        data[i] = 0;
    }
    uint32_t ind = 0;
    while (cpol[ind] != '\0')     {
        if (cpol[ind] >= '0' && cpol[ind] <= '9') {
            char c[10];
            uint32_t k=0;
            while(cpol[ind] >= '0' && cpol[ind]<= '9'){
                c[k++]=cpol[ind++];
            }
            c[k] = '\0';
            int val_mon =atoi(c);
            uint32_t i = val_mon&0x1f;
            uint32_t j = val_mon>>5;
            data[j] ^= 1<<i;
        }
        else {
            ind++;
        }
    }
}

//Give the correspondance between the place of a quadratic monomial in poly_quad and its place in a poly 
uint32_t pos_to_power (uint32_t pos, uint32_t size)    {
    if (size == 8)  {
        if (pos == 27)
            return 3;
        if (pos == 26)
            return 5;
        if (pos == 25)
            return 9;
        if (pos == 24)
            return 17;
        if (pos == 23)
            return 33;
        if (pos == 22)
            return 65;
        if (pos == 21)
            return 129;
        if (pos == 20)
            return 6;
        if (pos == 19)
            return 10;
        if (pos == 18)
            return 18;
        if (pos == 17)
            return 34;
        if (pos == 16)
            return 66;
        if (pos == 15)
            return 130;
        if (pos == 14)
            return 12;
        if (pos == 13)
            return 20;
        if (pos == 12)
            return 36;
        if (pos == 11)
            return 68;
        if (pos == 10)
            return 132;
        if (pos == 9)
            return 24;
        if (pos == 8)
            return 40;
        if (pos == 7)
            return 72;
        if (pos == 6)
            return 136;
        if (pos == 5)
            return 48;
        if (pos == 4)
            return 80;
        if (pos == 3)
            return 144;
        if (pos == 2)
            return 96;
        if (pos == 1)
            return 160;
        if (pos == 0)
            return 192;
        return 0;
    }
    if (size == 7)  {
        if (pos == 20)
            return 3;
        if (pos == 19)
            return 5;
        if (pos == 18)
            return 9;
        if (pos == 17)
            return 17;
        if (pos == 16)
            return 33;
        if (pos == 15)
            return 65;
        if (pos == 14)
            return 6;
        if (pos == 13)
            return 10;
        if (pos == 12)
            return 18;
        if (pos == 11)
            return 34;
        if (pos == 10)
            return 66;
        if (pos == 9)
            return 12;
        if (pos == 8)
            return 20;
        if (pos == 7)
            return 36;
        if (pos == 6)
            return 68;
        if (pos == 5)
            return 24;
        if (pos == 4)
            return 40;
        if (pos == 3)
            return 72;
        if (pos == 2)
            return 48;
        if (pos == 1)
            return 80;
        if (pos == 0)
            return 96;
    }
    if (size == 6)  {
        if (pos == 14)
            return 3;
        if (pos == 13)
            return 5;
        if (pos == 12)
            return 9;
        if (pos == 11)
            return 17;
        if (pos == 10)
            return 33;
        if (pos == 9)
            return 6;
        if (pos == 8)
            return 10;
        if (pos == 7)
            return 18;
        if (pos == 6)
            return 34;
        if (pos == 5)
            return 12;
        if (pos == 4)
            return 20;
        if (pos == 3)
            return 36;
        if (pos == 2)
            return 24;
        if (pos == 1)
            return 40;
        if (pos == 0)
            return 48;
    }
    if (size == 5)  {
        if (pos == 9)
            return 3;
        if (pos == 8)
            return 5;
        if (pos == 7)
            return 9;
        if (pos == 6)
            return 17;
        if (pos == 5)
            return 6;
        if (pos == 4)
            return 10;
        if (pos == 3)
            return 18;
        if (pos == 2)
            return 12;
        if (pos == 1)
            return 20;
        if (pos == 0)
            return 24;
    }
    if (size == 4)  {
        if (pos == 5)
            return 3;
        if (pos == 4)
            return 5;
        if (pos == 3)
            return 9;
        if (pos == 2)
            return 6;
        if (pos == 1)
            return 10;
        if (pos == 0)
            return 12;
    }
    return 0;
}

//Constructor from a poly_quad
poly::poly(poly_quad p, uint32_t size)   {
    for (uint32_t i=0; i<MAX; i++){
        data[i] = 0;
    }
    poly_quad temp_p = p;
    uint32_t nb_bits = (size*(size + 1))/2;
    for (uint32_t b = 0; b<nb_bits; b++) {
        if (temp_p&1ul) {  
            uint32_t val_mon = pos_to_power(b, size);
            uint32_t i = val_mon&0x1f;
            uint32_t j = val_mon>>5;
            data[j] ^= 1<<i;
        }
        temp_p >>= 1;
    }
}

bool poly::operator<(const poly& p) const {
    for (uint32_t i=0; i<MAX; i++){
        if (this->data[i] == p.data[i]){
            continue;
        }
        else {
            return this->data[i] < p.data[i];
        }
    }
    return this->data[0] < p.data[0];
}

bool poly::operator==(const poly& other) const {
    for (uint32_t i = 0; i < MAX; i++) {
        if (data[i] != other.data[i]) {
            return false;
        }
    }
    return true;
}

bool poly::operator!=(const poly& other) const {
    return !(*this == other);
}

bool poly::is_linear_monomial(uint32_t size){
    for (uint32_t i=0; i<size; i++){
        char * cvar = new char[5];
        sprintf(cvar, "x%u", 1<<i);
        poly var(cvar);
        if (*this == var){
            return true;
        }
    }
    return false;
}

//return the i-element of a poly
uint32_t poly::get(uint32_t i) const{
    return data[i];
}

void poly::print_poly(uint32_t size) const{
    // size represents the number of variables of the polynomial
    poly zero;
    if (!POLY_EQUAL(*this, zero)){
        bool xor_up = false;
        for (uint32_t u=0; u<MAX*32; u++) {
            if ((data[u>>5]>>(u&0x1f))&1ul) {
                if (xor_up){
                    cout << " ^ ";
                }
                if (u==0){
                    cout << "1";
                }
                else {
                    for (uint32_t bit_n=0; bit_n<size; bit_n++) {
                        if ((u>>bit_n)&1ul){
                            cout << "x" + to_string(bit_n);  
                        }
                    }
                }
                xor_up = true;
            }
        }
    }
    else {
        cout<<"0";
    }
}

string poly::poly_to_string(uint32_t size) const{
    string str;
    poly zero;
    uint32_t nb_elem = 1<<size;
    if (!POLY_EQUAL(*this, zero)){
        bool xor_up = false;
        for (uint32_t u=0; u<nb_elem; u++) {
            if ((data[u>>5]>>(u&0x1f))&1ul) {
                if (xor_up){
                    str += "^";
                }
                if (u==0){
                    str += "1";
                }
                else {
                    for (uint32_t bit_n=0; bit_n<size; bit_n++) {
                        if ((u>>bit_n)&1ul){
                            str += "x" + to_string(bit_n);
                        }
                    }
                }
                xor_up = true;
            }
        }
    }
    else {
        str += "0";
    }
    return str;
}

void poly::add(poly a) {
    data ^= a.data;
}

uint32_t poly::algebraic_degree (uint32_t nb_elem) {
    uint32_t max_deg = 0;
    for (uint32_t u=0; u<nb_elem; u++) {
        if ((data[u>>5]>>(u&0x1f))&1ul) {
            uint32_t w = weight(u);
            if (w > max_deg) {
                max_deg = w;
            }
        }
    }
    return max_deg;
}

uint64_t corres(uint32_t b1, uint32_t b2, int n) {
    uint32_t i;
    uint32_t j;
    if (b1 < b2)    {
        i = b1;
        j = b2;
    }
    else {
        i = b2;
        j = b1;
    }
    uint64_t res = (n *(n -1))/2 -1 -i*(n-1) + (i*(i-1)/2) - (j-i) +1;
    return res;
}

poly_quad poly_to_poly_quad (poly op, uint32_t size, uint32_t nb_elem)  { //To convert a poly into a poly_quad
    poly_quad p = 0;
    for (uint32_t u=0; u<nb_elem; u++) {
        if (((op.get(u>>5)>>(u&0x1f))&1ul)) {
            uint32_t temp_u = u;
            uint32_t bits_u[size];
            for (uint32_t b = 0 ; b<size; b++) {
                bits_u[b] = temp_u&1ul;
                temp_u >>= 1;
            }
            for (uint32_t b1 = 0; b1<size ; b1++) {
                for (uint32_t b2 = b1+1; b2<size; b2++) {
                    if ((bits_u[b1] == bits_u[b2]) && (bits_u[b1] == 1))    {
                        p ^= 1ul << (corres(b1, b2, size));
                    }
                }
            }
        }
    }
    return p;
}

poly truncate_lin(const poly& a, uint32_t size)  {
    poly res;
    COPY_POLY(a,res);
    for (uint32_t i=0; i<size; i++) {
        uint32_t u = 1<<i;
        if (a.get(u>>5)>>u&1ul) {
            res.data[u>>5] ^= 1<<u;
        }
    }
    if ((a.get(0))&1ul){
        res.data[0] ^= 1;
    }
    return res;
}

poly truncate_except_deg(const poly& a, uint32_t d, uint32_t size, uint32_t nb_elem)   {
    // return only the monomials of degree d in p
    poly res;
    for (uint32_t u=0; u<nb_elem; u++) {
        if ((a.get(u>>5)>>(u&0x1f))&1ul) {
            if (weight(u) == d)    {
                uint32_t val_mon = u;
                uint32_t k = val_mon&0x1f;
                uint32_t l = val_mon>>5;
                res.data[l] ^= 1<<k;
            }

        }
    }
    return res;
}

void create_monomial_order(uint32_t monomial_order [] ,uint32_t nb_elem,uint32_t size,uint32_t monom_dom){
    vector<uint32_t> mon_order;
    if ((size == 4) && (monom_dom == 0)){
        // Correspond to the lexical order induced by x0 > 1 > x1 > x2 > x3
        mon_order = {8,7,12,4,14,5,10,1,15,6,11,2,13,3,9,0};
    }
    if ((size == 4) && (monom_dom == 1)){
        mon_order = {8,12,7,4,14,10,5,1,15,11,6,2,13,9,3,0};
    }
    if ((size == 4) && (monom_dom == 2)){
        mon_order = {8,12,14,10,7,4,5,1,15,11,13,9,6,2,3,0};
    }
    if ((size == 4) && (monom_dom == 3)){
        mon_order = {8,12,14,10,15,11,13,9,7,4,5,1,6,2,3,0};
    }
    if ((size == 5) && (monom_dom == 0)){
        mon_order = {31,15,23,7,27,11,19,3,29,13,21,5,25,9,17,1,30,14,22,6,26,10,18,2,28,12,20,4,24,8,16,0};
    }
    if ((size == 5) && (monom_dom == 1)){
        mon_order = {31, 23, 15, 7, 27, 19, 11, 3, 29, 21, 13, 5, 25, 17, 9, 1, 30, 22, 14, 6, 26, 18, 10, 2, 28, 20, 12, 4, 24, 16, 8, 0};
    }
    if ((size == 5) && (monom_dom == 2)){
        mon_order = {31, 23, 27, 19, 15, 7, 11, 3, 29, 21, 25, 17, 13, 5, 9, 1, 30, 22, 26, 18, 14, 6, 10, 2, 28, 20, 24, 16, 12, 4, 8, 0};
    }
    if ((size == 5) && (monom_dom == 3)){
        mon_order = {31, 23, 27, 19, 29, 21, 25, 17, 15, 7, 11, 3, 13, 5, 9, 1, 30, 22, 26, 18, 28, 20, 24, 16, 14, 6, 10, 2, 12, 4, 8, 0};
    }
    if ((size == 5) && (monom_dom == 4)){
        mon_order = {31, 23, 27, 19, 29, 21, 25, 17, 30, 22, 26, 18, 28, 20, 24, 16, 15, 7, 11, 3, 13, 5, 9, 1, 14, 6, 10, 2, 12, 4, 8, 0};
    }
    if ((size == 6) && (monom_dom == 0)){
        mon_order = {63, 31, 47, 15, 55, 23, 39, 7, 59, 27, 43, 11, 51, 19, 35, 3, 61, 29, 45, 13, 53, 21, 37, 5, 57, 25, 41, 9, 49, 17, 33, 1, 62, 30, 46, 14, 54, 22, 38, 6, 58, 26, 42, 10, 50, 18, 34, 2, 60, 28, 44, 12, 52, 20, 36, 4, 56, 24, 40, 8, 48, 16, 32, 0};
    }
    if ((size == 6) && (monom_dom == 1)){
        mon_order = {63, 47, 31, 15, 55, 39, 23, 7, 59, 43, 27, 11, 51, 35, 19, 3, 61, 45, 29, 13, 53, 37, 21, 5, 57, 41, 25, 9, 49, 33, 17, 1, 62, 46, 30, 14, 54, 38, 22, 6, 58, 42, 26, 10, 50, 34, 18, 2, 60, 44, 28, 12, 52, 36, 20, 4, 56, 40, 24, 8, 48, 32, 16, 0};
    }
    if ((size == 6) && (monom_dom == 2)){
        mon_order = {63, 47, 55, 39, 31, 15, 23, 7, 59, 43, 51, 35, 27, 11, 19, 3, 61, 45, 53, 37, 29, 13, 21, 5, 57, 41, 49, 33, 25, 9, 17, 1, 62, 46, 54, 38, 30, 14, 22, 6, 58, 42, 50, 34, 26, 10, 18, 2, 60, 44, 52, 36, 28, 12, 20, 4, 56, 40, 48, 32, 24, 8, 16, 0};
    }
    if ((size == 6) && (monom_dom == 3)){
        mon_order = {63, 47, 55, 39, 59, 43, 51, 35, 31, 15, 23, 7, 27, 11, 19, 3, 61, 45, 53, 37, 57, 41, 49, 33, 29, 13, 21, 5, 25, 9, 17, 1, 62, 46, 54, 38, 58, 42, 50, 34, 30, 14, 22, 6, 26, 10, 18, 2, 60, 44, 52, 36, 56, 40, 48, 32, 28, 12, 20, 4, 24, 8, 16, 0};
    }
    if ((size == 6) && (monom_dom == 4)){
        mon_order = {63, 47, 55, 39, 59, 43, 51, 35, 61, 45, 53, 37, 57, 41, 49, 33, 31, 15, 23, 7, 27, 11, 19, 3, 29, 13, 21, 5, 25, 9, 17, 1, 62, 46, 54, 38, 58, 42, 50, 34, 60, 44, 52, 36, 56, 40, 48, 32, 30, 14, 22, 6, 26, 10, 18, 2, 28, 12, 20, 4, 24, 8, 16, 0};
    }
    if ((size == 6) && (monom_dom == 5)){
        mon_order = {63, 62, 47, 46, 55, 54, 39, 38, 59, 58, 43, 42, 51, 50, 35, 34, 61, 60, 45, 44, 53, 52, 37, 36, 57, 56, 41, 40, 49, 48, 33, 32, 31, 30, 15, 14, 23, 22, 7, 6, 27, 26, 11, 10, 19, 18, 3, 2, 29, 28, 13, 12, 21, 20, 5, 4, 25, 24, 9, 8, 17, 16, 1, 0};
    }
    if ((size == 7) && (monom_dom == 0)){
        mon_order = {127, 63, 95, 31, 111, 47, 79, 15, 119, 55, 87, 23, 103, 39, 71, 7, 123, 59, 91, 27, 107, 43, 75, 11, 115, 51, 83, 19, 99, 35, 67, 3, 125, 61, 93, 29, 109, 45, 77, 13, 117, 53, 85, 21, 101, 37, 69, 5, 121, 57, 89, 25, 105, 41, 73, 9, 113, 49, 81, 17, 97, 33, 65, 1, 126, 62, 94, 30, 110, 46, 78, 14, 118, 54, 86, 22, 102, 38, 70, 6, 122, 58, 90, 26, 106, 42, 74, 10, 114, 50, 82, 18, 98, 34, 66, 2, 124, 60, 92, 28, 108, 44, 76, 12, 116, 52, 84, 20, 100, 36, 68, 4, 120, 56, 88, 24, 104, 40, 72, 8, 112, 48, 80, 16, 96, 32, 64, 0};
    }
    if ((size == 7) && (monom_dom == 1)){
        mon_order = {127,63,95,31,111,47,79,15,119,55,87,23,103,39,71,7,123,59,91,27,107,43,75,11,115,51,83,19,99,35,67,3,126,62,94,30,110,46,78,14,118,54,86,22,102,38,70,6,122,58,90,26,106,42,74,10,114,50,82,18,98,34,66,2,125,61,93,29,109,45,77,13,117,53,85,21,101,37,69,5,121,57,89,25,105,41,73,9,113,49,81,17,97,33,65,1,124,60,92,28,108,44,76,12,116,52,84,20,100,36,68,4,120,56,88,24,104,40,72,8,112,48,80,16,96,32,64,0};
    }
    if ((size == 7) && (monom_dom == 2)){
        mon_order = {127, 95, 111, 79, 63, 31, 47, 15, 119, 87, 103, 71, 55, 23, 39, 7, 123, 91, 107, 75, 59, 27, 43, 11, 115, 83, 99, 67, 51, 19, 35, 3, 125, 93, 109, 77, 61, 29, 45, 13, 117, 85, 101, 69, 53, 21, 37, 5, 121, 89, 105, 73, 57, 25, 41, 9, 113, 81, 97, 65, 49, 17, 33, 1, 126, 94, 110, 78, 62, 30, 46, 14, 118, 86, 102, 70, 54, 22, 38, 6, 122, 90, 106, 74, 58, 26, 42, 10, 114, 82, 98, 66, 50, 18, 34, 2, 124, 92, 108, 76, 60, 28, 44, 12, 116, 84, 100, 68, 52, 20, 36, 4, 120, 88, 104, 72, 56, 24, 40, 8, 112, 80, 96, 64, 48, 16, 32, 0};
    }
    if ((size == 7) && (monom_dom == 3)){
        mon_order = {127, 95, 111, 79, 119, 87, 103, 71, 63, 31, 47, 15, 55, 23, 39, 7, 123, 91, 107, 75, 115, 83, 99, 67, 59, 27, 43, 11, 51, 19, 35, 3, 125, 93, 109, 77, 117, 85, 101, 69, 61, 29, 45, 13, 53, 21, 37, 5, 121, 89, 105, 73, 113, 81, 97, 65, 57, 25, 41, 9, 49, 17, 33, 1, 126, 94, 110, 78, 118, 86, 102, 70, 62, 30, 46, 14, 54, 22, 38, 6, 122, 90, 106, 74, 114, 82, 98, 66, 58, 26, 42, 10, 50, 18, 34, 2, 124, 92, 108, 76, 116, 84, 100, 68, 60, 28, 44, 12, 52, 20, 36, 4, 120, 88, 104, 72, 112, 80, 96, 64, 56, 24, 40, 8, 48, 16, 32, 0};
    }
    if ((size == 7) && (monom_dom == 4)){
        mon_order = {127, 95, 111, 79, 119, 87, 103, 71, 123, 91, 107, 75, 115, 83, 99, 67, 63, 31, 47, 15, 55, 23, 39, 7, 59, 27, 43, 11, 51, 19, 35, 3, 125, 93, 109, 77, 117, 85, 101, 69, 121, 89, 105, 73, 113, 81, 97, 65, 61, 29, 45, 13, 53, 21, 37, 5, 57, 25, 41, 9, 49, 17, 33, 1, 126, 94, 110, 78, 118, 86, 102, 70, 122, 90, 106, 74, 114, 82, 98, 66, 62, 30, 46, 14, 54, 22, 38, 6, 58, 26, 42, 10, 50, 18, 34, 2, 124, 92, 108, 76, 116, 84, 100, 68, 120, 88, 104, 72, 112, 80, 96, 64, 60, 28, 44, 12, 52, 20, 36, 4, 56, 24, 40, 8, 48, 16, 32, 0};
    }
    if ((size == 7) && (monom_dom == 5)){
        mon_order = {127, 95, 111, 79, 119, 87, 103, 71, 123, 91, 107, 75, 115, 83, 99, 67, 125, 93, 109, 77, 117, 85, 101, 69, 121, 89, 105, 73, 113, 81, 97, 65, 126, 94, 110, 78, 118, 86, 102, 70, 122, 90, 106, 74, 114, 82, 98, 66, 124, 92, 108, 76, 116, 84, 100, 68, 120, 88, 104, 72, 112, 80, 96, 64, 63, 31, 47, 15, 55, 23, 39, 7, 59, 27, 43, 11, 51, 19, 35, 3, 61, 29, 45, 13, 53, 21, 37, 5, 57, 25, 41, 9, 49, 17, 33, 1, 62, 30, 46, 14, 54, 22, 38, 6, 58, 26, 42, 10, 50, 18, 34, 2, 60, 28, 44, 12, 52, 20, 36, 4, 56, 24, 40, 8, 48, 16, 32, 0};
    }
    if ((size == 7) && (monom_dom == 6)){
        mon_order = {127, 126, 95, 94, 111, 110, 79, 78, 119, 118, 87, 86, 103, 102, 71, 70, 123, 122, 91, 90, 107, 106, 75, 74, 115, 114, 83, 82, 99, 98, 67, 66, 63, 62, 31, 30, 47, 46, 15, 14, 55, 54, 23, 22, 39, 38, 7, 6, 59, 58, 27, 26, 43, 42, 11, 10, 51, 50, 19, 18, 35, 34, 3, 2, 125, 124, 93, 92, 109, 108, 77, 76, 117, 116, 85, 84, 101, 100, 69, 68, 121, 120, 89, 88, 105, 104, 73, 72, 113, 112, 81, 80, 97, 96, 65, 64, 61, 60, 29, 28, 45, 44, 13, 12, 53, 52, 21, 20, 37, 36, 5, 4, 57, 56, 25, 24, 41, 40, 9, 8, 49, 48, 17, 16, 33, 32, 1, 0};
    }
    for (uint32_t i=0; i<nb_elem; i++){
        monomial_order[i] = mon_order[i];
    }
}

uint32_t leading_term(const poly& a, uint32_t monomial_order [], uint32_t nb_elem) {
    bool test = false ;
    uint32_t b=nb_elem;
    uint32_t val_res = 1;
    for (uint32_t j=0; j<MAX && !test; j++){
        if ((32 * j) > b){
            test = true;
            continue;
        }
        uint32_t temp = a.get(j);
        if (temp != 0){
            for (uint32_t i=0 ; i<32; i++,temp>>=1){
                if (temp == 0){
                    i=31;
                    continue;
                }
                else {
                    if (temp&1u){
                        uint32_t val_mon = (32 * j) + i;
                        uint32_t b2 = monomial_order[val_mon];
                        if (b2 < b){
                            b = b2;
                            val_res = val_mon;
                        }
                    }
                }
            }
        }
    }
    return val_res ;
}

bool is_divisible(uint32_t a, uint32_t b) {
    return ((a|b) == a);
}

poly mon_to_poly(uint32_t mon){
    poly res;
    res.data[mon>>5] ^= 1u<<mon;
    return res;
}

pair<poly, poly> poly_div_no_truncate (const poly& a, const poly& b, uint32_t monomial_order [], uint32_t size, uint32_t nb_elem) {

    poly q; 
    poly r; 
    COPY_POLY(a, r);

    poly zero; 

    if (b == zero)    {
        return {q,r};
    }

    uint32_t lt = leading_term(r, monomial_order, nb_elem);
    uint32_t qt = leading_term(b, monomial_order, nb_elem);

    uint32_t m = 0;
    uint32_t previous_m = nb_elem + 1;
    uint32_t compteur = 0;

    while (is_divisible(lt, qt) && !POLY_EQUAL(r,zero) && (compteur < 30)) {

        m = MONOM_DIV(qt, lt);

        if (m==previous_m){
            return {q,r};
        }

        poly p;
        POLY_MUL2(m,b,p,nb_elem);

        POLY_ADD(r, p, r);

        poly poly_m = mon_to_poly(m);
        POLY_ADD(q, poly_m, q);

        lt = leading_term(r, monomial_order, nb_elem);
        previous_m = m; 
        compteur ++;

    }
    
    return {q, r};
}

pair<poly, poly> poly_div (poly a, poly b, uint32_t monomial_order [], uint32_t size, uint32_t nb_elem) {

    poly q; 
    poly r; 
    COPY_POLY(a, r);

    poly zero; 

    if (POLY_EQUAL(b,zero))    {
        return {q,r};
    }

    uint32_t lt = leading_term(r, monomial_order, nb_elem);
    uint32_t qt = leading_term(b, monomial_order, nb_elem);

    uint32_t mon = 0;
    uint32_t previous_mon = nb_elem + 1;  
    uint32_t compteur = 0;

    // While the leading term of the dividend is divisble by the leading term of the divisor we continue
    while (is_divisible(lt, qt) && !POLY_EQUAL(r,zero) && (compteur < 30)) {

        mon = MONOM_DIV(qt, lt);

        if (mon==previous_mon){
            return {q,r};
        }

        poly p;
        POLY_MUL2(mon,b,p,nb_elem);

        POLY_ADD(r, p, r);
        r = truncate_lin(r, size);

        poly poly_m = mon_to_poly(mon);
        POLY_ADD(q, poly_m, q);

        lt = leading_term(r, monomial_order, nb_elem);
        previous_mon = mon; 
        compteur ++;

    }
    
    return {q, r};
}

bool mon_dom_in (poly l, uint32_t mon_dom, uint32_t size, uint32_t nb_elem){
    poly test;

    poly m;
    uint32_t val_mon = 1<<mon_dom;
    uint32_t i = val_mon&0x1f;
    uint32_t j = val_mon>>5;
    m.data[j] ^= 1<<i;

    test.data[j] = l.data[j]&m.data[j];

    return(test.data[j]!=0);
}

uint32_t poly_deg5_4_to_uint(poly p){
    vector<uint32_t> mon = {60,58,54,46,30,57,53,45,29,51,43,27,39,23,15,62,61,59,55,47,31};
    uint32_t res=0;
    for (uint32_t i=0; i<mon.size(); i++){
        uint32_t val_mon = mon[i];
        if ((p.get(val_mon>>5)>>(val_mon&0x1f))&1u){
            res += 1<<i;
        }
    }
    return res - 32768;
}

