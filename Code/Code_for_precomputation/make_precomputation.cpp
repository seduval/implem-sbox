#include "config.h"
#include "sbox.h"
#include "poly.h"
#include "precalcul.h"

// Here is an example of how we build the set map_xor_l_plus_lin for size 4.

using namespace std;

string print_poly2(poly p,uint32_t size) {
    // size represents the number of variables of the polynom
    string str;
    uint32_t nb_elem = 1<<size;
    poly zero;
    if (!POLY_EQUAL(p, zero)){
        bool xor_up = false;
        for (uint32_t u=0; u<nb_elem; u++) {
            if ((p.get(u>>5)>>(u&0x1f))&1ul) {
                if (xor_up){
                    //cout << " ^ ";
                    str += "^";
                }
                if (u==0){
                    //cout << "1";
                    str += "1";
                }
                else {
                    str += "x" + to_string(u);
                }
                xor_up = true;
            }
        }
    }
    else {
        //cout<<"0";
        str += "0";
    }
    return str;
}

int main(int argc, char **argv)  {

    uint32_t size_in = 4;
    uint32_t nb_elem = 16;

    uint32_t set_op_size = nb_bits_to_size_set(size_in);
    poly_quad set_op[set_op_size];
    parse_file_and_create_set(size_in, set_op);

    uint32_t map_xor_size = nb_bits_to_size_map_xor(size_in);
    pair_xor * map_xor = new pair_xor[map_xor_size];
    parse_file_and_create_map_xor(size_in, map_xor);

    vector<poly> l = create_poly_lin(size_in, nb_elem);

    set<poly> map_xor_l_plus_lin;

    for (uint32_t i=0; i<map_xor_size; i++){
        poly p(map_xor[i].first, size_in);
        map_xor_l_plus_lin.insert(p);
        for (uint32_t j=0; j<l.size(); j++){
            poly p2;
            POLY_ADD(p, l[j], p2);
            map_xor_l_plus_lin.insert(p2);
        }
    }

    for (auto it=map_xor_l_plus_lin.begin(); it != map_xor_l_plus_lin.end(); it++){
        cout<<print_poly2(*it, size_in)<<endl;
    }

    exit(0);
}