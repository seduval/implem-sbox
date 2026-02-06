#include "config.h"
#include "sbox.h"
#include "poly.h"
#include "precomputation.h"

using namespace std;

uint32_t nb_bits_to_size_set (uint32_t size){
    uint32_t res;
    if (size == 4)  {
        res = 35;
    }
    if (size == 5)  {
        res = 155;
    }
    if (size == 6)  {
        res = 651;
    }
    if (size == 7)  {
        res = 2667;
    }
    if (size == 8)  {
        res = 10795;
    }
    if (size == 9)  {
        res = 43435;
    }
    if (size == 10)  {
        res = 174251;
    }
    return res;
}

void parse_file_and_create_set(uint32_t size, poly_quad set_op []) { // To create the tab that contains the quadratic operations
    
    string filename = string("precomputation_files/precomputed_set_") + to_string(size) + string(".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    uint32_t it = 0;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        size_t pos = line.find("{");
        while (pos != string::npos) {
            size_t endPos = line.find("}", pos); //Position du caractère '}'
            if (endPos == string::npos) 
                break;

            string set_str = line.substr(pos + 1, endPos - pos - 1); //On extrait la chaîne de caractère se trouvant entre la position de '{' et de '}'

            size_t pos_valeur = set_str.find_first_of("0123456789"); //Position de la première valeur
            while (pos_valeur != string::npos) { 
                size_t pos_virgule = set_str.find(",", pos_valeur); //position du caractère ','
                if (pos_virgule == string::npos)
                    pos_virgule = set_str.size();
                
                poly_quad value = stoul(set_str.substr(pos_valeur, pos_virgule - pos_valeur)); //On convertit la chaîne de caractère situé entre chaque virgule en valeur  

                set_op[it] = value; // on insère la valeur au set
                it++;

                pos_valeur = set_str.find_first_of("0123456789", pos_virgule + 1);
            }
            pos = line.find("{", endPos);
        }
    }

    file.close();
}

int64_t find_set (poly_quad val, poly_quad set_op [], uint64_t look_start, uint64_t look_end) {
    //To find a given element in a poly[]
    if (look_start==look_end){
        return -1;
    }
    int64_t mid = (look_start + look_end)/2;
    int64_t test = val - set_op[mid];
    return (test?(test<0?find_set(val,set_op,look_start,mid):find_set(val,set_op,mid + 1,look_end)):mid);
}

uint32_t nb_bits_to_size_map_xor (uint32_t size){
    uint32_t res;
    if (size == 4)  {
        res = 28;
    }
    if (size == 5)  {
        res = 868;
    }
    if (size == 6)  {
        res = 18228;
    }
    if (size == 7)  {
        res = 330708;
    }
    if (size == 8)  {
        res = 5622036;
    }
    if (size == 9)  {
        res = 92672916;
    }
    return res;
}

void parse_file_and_create_map_xor(uint32_t size, pair_xor map_xor []) {
    
    string filename = string("precomputation_files/precomputed_map_xor_") + to_string(size) + string(".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    uint32_t i = 0; 

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        size_t pos = line.find("[");
        if (pos == string::npos)
            continue;

        poly_quad key = stoul(line.substr(pos + 1)); 
        map_xor[i].first = key;
        
        uint32_t j=0;

        pos = line.find("{", pos); 
        while (pos != string::npos) {
            size_t endPos = line.find("}", pos);
            if (endPos == string::npos) 
                break;

            string set_str = line.substr(pos + 1, endPos - pos - 1); 
           
            uint32_t k=0; 

            size_t pos_valeur = set_str.find_first_of("0123456789"); 
            while (pos_valeur != string::npos) { 
                size_t pos_virgule = set_str.find(",", pos_valeur); 
                if (pos_virgule == string::npos)
                    pos_virgule = set_str.size();
                
                poly_quad value = stoul(set_str.substr(pos_valeur, pos_virgule - pos_valeur)); 
                map_xor[i].second[j][k] = value;

                pos_valeur = set_str.find_first_of("0123456789", pos_virgule + 1);
                k++;
            }

            pos = line.find("{", endPos);
            j++;
        }
        i++;
    }

    file.close();
}

int64_t find_map (poly_quad val, pair_xor v [], uint64_t look_start, uint64_t look_end)   {
    //To find a given element in a pair_xor []
    if (look_start==look_end){
        return -1;
    }
    int64_t mid = (look_start + look_end)/2;
    int64_t test = val - v[mid].first;
    return (test?(test<0?find_map(val,v,look_start,mid):find_map(val,v,mid + 1,look_end)):mid);
}

vector<poly> create_poly_lin(uint32_t size, uint32_t nb_elem)    {

    string filename = string("precomputation_files/precomputed_lin_") + to_string(size) + string(".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    vector<poly> tab_lin;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        const char* str = line.c_str();
        poly p(str);
        tab_lin.push_back(p);
    }
    file.close();
    return tab_lin;
}

vector<poly> create_set_op_l_plus_lin(uint32_t size)    {

    string filename = string("precomputation_files/precomputed_set_op_l_plus_lin_") + to_string(size) + string(".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    vector<poly> tab_lin;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        const char* str = line.c_str();
        poly p(str);
        tab_lin.push_back(p);
    }
    file.close();
    return tab_lin;
}

vector<poly> create_map_xor_l_plus_lin(uint32_t size)    {

    string filename = string("precomputation_files/precomputed_map_xor_l_plus_lin_") + to_string(size) + string(".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    vector<poly> tab_lin;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        const char* str = line.c_str();
        poly p(str);
        tab_lin.push_back(p);
    }
    file.close();
    return tab_lin;
}

vector<implem_deg5_4> parse_file_and_create_sol(uint32_t ent_p) { 
    
    uint32_t m = ent_p/100000; //To search in the right file
    string filename = string("precomputation_files/precomputed_sol_deg5_4/precomputed_sol_deg5_4_" + to_string(m) + ".txt");

    FILE* fp = fopen(filename.c_str(), "r");
    if (fp == NULL) {
        cerr << "Couldn't open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    vector<implem_deg5_4> vect_res;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {

        size_t pos_ent = line.find('{');
        uint32_t k = stoul(line.substr(0, pos_ent));

        if (k == ent_p){

            size_t pos = 0;
            while ((pos = line.find('{', pos)) != string::npos) {
                size_t debut = pos + 1; 
                size_t fin = line.find('}', debut); 
                if (fin != string::npos) {
                    string s_imp = line.substr(debut, fin - debut);

                    implem_deg5_4 imp;

                    size_t pos_virgule = s_imp.find(","); 
                    size_t pos_virgule2 = s_imp.find(",",pos_virgule+1);

                    string set_str = s_imp.substr(0, pos_virgule); 

                    const char* str = set_str.c_str();
                    poly pq(str);
                    imp.quad = pq;
    
                    string set_str2 = s_imp.substr(pos_virgule+1, pos_virgule2 - (pos_virgule+1)); 

                    const char* str2 = set_str2.c_str();
                    poly quotient(str2);
                    imp.q = quotient;

                    string set_str3 = s_imp.substr(pos_virgule2+1); 

                    const char* str3 = set_str3.c_str();
                    poly reste(str3);
                    imp.r = reste;

                    vect_res.push_back(imp);

                    pos = fin + 1; 
                } else {
                    break; // if no '}', we stop
                }
            }
        }
    }
    file.close();
    return vect_res;
}