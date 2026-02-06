#include "config.h"
#include "sbox.h"
#include "poly.h"
#include "precomputation.h"
#include "decomp.h"
#include "reconstitution.h"

using namespace std;

int main(int argc, char **argv)  {

    uint32_t size_in = 0;
    uint32_t size_out = 0;
    uint32_t nb_and_max = 0;
    uint32_t nb_sol = 1;
    uint32_t nb_elem = 0;
    vector<uint32_t> v_lut;

    option longopts[] = {
        {"lut", required_argument, NULL, 'l'},
        {"solmax", optional_argument, NULL, 's'},
        {"andmax", required_argument, NULL, 'a'},
        {"help", no_argument, NULL, 'h'}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "l:s:a:h", longopts, 0)) != -1) {
        switch (opt) {
             case 'h': {
                cout<<" The option '--lut' or '-l' requires the values of the LUT separated by ','"<<endl
                <<" For example: '-l 0,1,5,4,8,3,7,2'"<<endl
                <<" The number of values must be a power of 2 beetween 16 and 64"<<endl
                <<" The option ' --andmax' or '-a' allows you to choose a maximum for the number of AND gates wanted in the implementation"<<endl
                <<" The other options are optionnal : "<<endl
                <<" The option '--solmax' or '-s' allows you to choose the number of different solutions wanted, by default the programm will try to find all of them"<<endl;
                exit(0);
            }
            case 'l':   {
                string lut_values (optarg);
                istringstream iss(lut_values);
                string value;
                while (getline(iss, value, ',')) {
                    v_lut.push_back(stoi(value));
                }
                nb_elem = v_lut.size();
                if ((nb_elem > 128) || (nb_elem < 16) || (weight(nb_elem) != 1)) {
                    cout<<"The number of values in the LUT must be a power of 2 beetween 16 and 128 "<<endl;
                    exit(0);
                }
                size_in = log2(nb_elem);
                size_out = size_in;
                break;
            }
            case 's': {
                nb_sol = atoi(optarg);
                break;
            }
            case 'a': {
                nb_and_max = atoi(optarg);
                break;
            }
            default:
                cerr << "Unknown program argument " << opt << endl;
        }
    }
    uint32_t lut [nb_elem];
    for (uint32_t i=0; i<nb_elem; i++)  {
        lut[i] = v_lut[i];
    }

    uint32_t nb_words = nb_elem/32 + (nb_elem%32?1:0);

    uint32_t ** anf = (uint32_t**)malloc(size_out * sizeof(uint32_t*));
    for(uint32_t i=0; i<size_out; i++)   {
        anf[i] = (uint32_t*)calloc(nb_words , sizeof(uint32_t));
    }

    lut_to_anf(anf, lut, size_in, size_out, nb_elem, nb_words);

    print_anf(anf,size_in, size_out);

    vector<poly> Y;
    vector<poly> ANF;
    for (uint32_t i = 0; i < size_out; i++){
        poly y;
        poly y_non_lin;
        for (uint32_t j=0; j<nb_words; j++){
            y.data[j] = anf[i][j];
        }
        for (uint32_t j=nb_words + 1; j<MAX; j++){
            y.data[j] = 0;
        }
        ANF.push_back(y);
        y_non_lin = truncate_lin(y, size_in);
        Y.push_back(y_non_lin);
    }

    uint32_t set_op_size = nb_bits_to_size_set(size_in);
    poly_quad set_op[set_op_size];
    parse_file_and_create_set(size_in, set_op);

    uint32_t map_xor_size = nb_bits_to_size_map_xor(size_in);
    pair_xor * map_xor = new pair_xor[map_xor_size];
    parse_file_and_create_map_xor(size_in, map_xor);

    vector<poly> l = create_poly_lin(size_in, nb_elem);
    vector<poly> set_op_l_plus_lin = create_set_op_l_plus_lin(size_in);

    return_implem(size_in, size_out, nb_elem, nb_and_max, nb_sol, Y, l, set_op_l_plus_lin, set_op, map_xor, set_op_size, map_xor_size, ANF);

    
    exit(0);
}

