#include "../../config.h"
#include "../../sbox.h"
#include "../../poly.h"
#include "../../precalcul.h"

// Here is an example of how we generate precomputation files for degree 5 functions on 6 bits.

vector<poly> create_vect_poly_deg5_4()    {

    string filename = string("../precomputed_poly_deg5_4.txt");

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

vector<poly> create_set_op_l_plus_lin_6()    {

    string filename = string("../precomputed_set_op_l_plus_lin.txt");

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

int main(){

    uint32_t size_in = 6;
    uint32_t nb_elem = 1<<size_in;

    vector<poly> set_op_l_plus_lin = create_set_op_l_plus_lin_6();
    vector<poly> vect_deg5_4 = create_vect_poly_deg5_4();

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic) 
        for (uint32_t k=1807913; k<1900000; k++){
            set<set<poly>> set_sol;
            string imp = to_string(k);
            for (uint32_t i=0; i<set_op_l_plus_lin.size(); i++){
                for (uint32_t j=0; j<size_in; j++){

                    uint32_t monomial_order [nb_elem];
                    create_monomial_order(monomial_order, nb_elem, size_in, j);                                                   
                    pair<poly,poly> p;
                    p = poly_div_no_truncate(vect_deg5_4[k], set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                    uint32_t m = p.second.algebraic_degree(nb_elem);
                    uint32_t n = p.first.algebraic_degree(nb_elem);

                    if (m<=3 && n<=3){
                        set<poly> solution;
                        solution.insert(p.first);
                        solution.insert(p.second);
                        uint32_t s = set_sol.size();
                        set_sol.insert(solution);
                        if (set_sol.size() != s){
                            imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(p.first,size_in) + "," + print_poly2(p.second,size_in) + "};";
                        }
                    }
                    else {
                        if (n!=0 && n<=3){
                            for (uint32_t j2=0; j2<size_in; j2++){

                                if (j2==j){
                                    continue;
                                }

                                create_monomial_order(monomial_order, nb_elem, size_in, j2);                                                   
                                pair<poly,poly> p2;
                                p2 = poly_div(p.second, set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                                uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                                if (m2<=3 && n2<=3){
                                    set<poly> solution;
                                    poly q;
                                    POLY_ADD(p.first, p2.first, q);
                                    solution.insert(q);
                                    solution.insert(p2.second);
                                    uint32_t s = set_sol.size();
                                    set_sol.insert(solution);
                                    if (set_sol.size() != s){
                                        imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(q,size_in) + "," + print_poly2(p2.second,size_in) + "};";
                                    }
                                }

                                create_monomial_order2(monomial_order, nb_elem, size_in, j2);                                                   
                                pair<poly,poly> pp2;
                                pp2 = poly_div(p.second, set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                                uint32_t mm2 = pp2.second.algebraic_degree(nb_elem);
                                uint32_t nn2 = pp2.first.algebraic_degree(nb_elem);

                                if (mm2<=3 && nn2<=3){
                                    set<poly> solution;
                                    poly q;
                                    POLY_ADD(p.first, pp2.first, q);
                                    solution.insert(q);
                                    solution.insert(pp2.second);
                                    uint32_t s = set_sol.size();
                                    set_sol.insert(solution);
                                    if (set_sol.size() != s){
                                        imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(q,size_in) + "," + print_poly2(pp2.second,size_in) + "};";
                                    }      
                                }
                            }
                        }
                    }

                    create_monomial_order2(monomial_order, nb_elem, size_in, j);                                                   
                    pair<poly,poly> pp;
                    pp = poly_div_no_truncate(vect_deg5_4[k], set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                    uint32_t mm = pp.second.algebraic_degree(nb_elem);
                    uint32_t nn = pp.first.algebraic_degree(nb_elem);

                    if (mm<=3 && nn<=3){
                        set<poly> solution;
                        solution.insert(pp.first);
                        solution.insert(pp.second);
                        uint32_t s = set_sol.size();
                        set_sol.insert(solution);
                        if (set_sol.size() != s){
                            imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(pp.first,size_in) + "," + print_poly2(pp.second,size_in) + "};";
                        }
                    }

                    else {
                        if (nn!=0 && nn<=3){
                            for (uint32_t j2=0; j2<size_in; j2++){
                                if (j2 == j){
                                    continue;
                                }

                                create_monomial_order(monomial_order, nb_elem, size_in, j2);                                                   
                                pair<poly,poly> p2;
                                p2 = poly_div(pp.second, set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                                uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                                if (m2<=3 && n2<=3){
                                    set<poly> solution;
                                    poly q;
                                    POLY_ADD(pp.first, p2.first, q);
                                    solution.insert(q);
                                    solution.insert(p2.second);
                                    uint32_t s = set_sol.size();
                                    set_sol.insert(solution);
                                    if (set_sol.size() != s){
                                        imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(q,size_in) + "," + print_poly2(p2.second,size_in) + "};";
                                    }
                                }

                                create_monomial_order2(monomial_order, nb_elem, size_in, j2);                                                   
                                pair<poly,poly> pp2;
                                pp2 = poly_div(pp.second, set_op_l_plus_lin[i], monomial_order, size_in, nb_elem);

                                uint32_t mm2 = pp2.second.algebraic_degree(nb_elem);
                                uint32_t nn2 = pp2.first.algebraic_degree(nb_elem);

                                if (mm2<=3 && nn2<=3){
                                    set<poly> solution;
                                    poly q;
                                    POLY_ADD(pp.first, pp2.first, q);
                                    solution.insert(q);
                                    solution.insert(pp2.second);
                                    uint32_t s = set_sol.size();
                                    set_sol.insert(solution);
                                    if (set_sol.size() != s){
                                        imp += "{" + print_poly2(set_op_l_plus_lin[i],size_in) + "," + print_poly2(q,size_in) + "," + print_poly2(pp2.second,size_in) + "};";
                                    }
                                }
                            }
                        }
                    }
                }
            }
            #pragma omp critical 
            {
                cout<<imp<<endl;
            }
        }
    }
    
    /*string imp;

        const char * cpol = "x15^x29^x43^x46^x47^x51^x54^x55^x58^x60^x61\0";
        poly p_test(cpol);
        imp += "p = " + p_test.print_poly(size_in);
        uint32_t ent_p = poly_deg5_4_to_uint(p_test);
        imp += "\nent_p = " + to_string(ent_p) + "\n";
        cout<<imp<<endl;
        
        string imp2;
        imp2 += print_poly2(p_test, size_in) + "\n";
        cout<<imp2;*/

        /*const char * cpol2 = "x23^x29^x30^x31^x39^x43^x45^x46^x47^x51^x53^x57^x62\0";
        poly p_test2(cpol2);
        imp +="p2 = ";
        p_test2.print_poly(size_in);
        imp +=endl;
        uint32_t ent_p2 = poly_deg5_4_to_uint(p_test2);
        imp +=ent_p2<<endl;

        const char * cpol = "x15^x23^x27^x29^x30^x31^x46^x47^x54^x55^x58^x59^x60^x61^x62\0";
        poly p_test(cpol);
        imp += "p = " + p_test.print_poly(size_in);
        uint32_t ent_p = poly_deg5_4_to_uint(p_test);
        imp += "\nent_p = " + to_string(ent_p) + "\n";
        cout<<imp<<endl;*/
    exit(0);
}


