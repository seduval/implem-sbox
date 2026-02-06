#include "config.h"
#include "sbox.h"
#include "poly.h"
#include "precomputation.h"
#include "decomp.h"

vector<implem> remove_duplicates(vector<implem> v){
    vector<implem> vect_res;
    set<set<poly_quad>> set_comp;

    uint32_t s_size = 0;
    for (uint32_t i=0; i<v.size(); i++){
        set<poly_quad> s = v[i].quad_sol;
        set_comp.insert(s);
        if (set_comp.size() != s_size){
            s_size = set_comp.size();
            vect_res.push_back(v[i]);
        }
    }
    return vect_res;
}

uint64_t most_significant_bit(poly_quad p){
    if (p == 0) {
        return 0;
    }
    return 1<<(63-__builtin_clzl(p));
}

uint32_t Rank(set<poly_quad> op_selec){

    uint32_t rank = 0;
    while ( !op_selec.empty() ){

        auto itf = op_selec.end();
        itf--;
        uint64_t b = most_significant_bit(*itf);
        
        set<poly_quad> new_set;
        
        for (auto it = op_selec.begin(); it != itf; it++){
            if ( (*it) >= b){
                poly_quad q = *it ^ *itf;
                new_set.insert(q);
            } 
            else {
                new_set.insert(*it); 
            }
        }

        op_selec = move(new_set);
        rank++;
    }
    return rank;
}

// Default Function

vector<implem> add_to_op_selec_dft(poly y, vector<poly> l, vector<poly> l2, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    return res;
}

// DEG 2

vector<implem> add_to_op_selec_deg_2_size_4_5_6_v1(poly y, vector<poly> l, vector<poly> l2, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    
    poly_quad pq1 = poly_to_poly_quad(y,size, nb_elem);
                    
    int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

    if (indice1 != -1)   {
        implem imp;
        imp.op_sol.push_back(y);
        imp.quad_sol.insert(pq1);
        imp.formula = "p1";
        res.push_back(imp);
    }
    else {
        indice1 = find_map(pq1, map_xor, 0, map_xor_size);
        if (indice1 == -1)   {
            return res;
        }
        else {
            for (uint32_t m=0; m<10; m++)   {
                implem imp;
                poly_quad op_q1 = map_xor[indice1].second[m][0];
                poly op1(op_q1,size);
                imp.op_sol.push_back(op1);
                imp.quad_sol.insert(op_q1);
                poly_quad op_q2 = map_xor[indice1].second[m][1];
                poly op2(op_q2,size);
                imp.op_sol.push_back(op2);
                imp.quad_sol.insert(op_q2);
                imp.formula = "p1 ^ p2";
                res.push_back(imp);
            }
        }
    }

    return res;
}

//DEG 3

vector<implem> add_to_op_selec_deg_3_v1(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    for (uint32_t i=0; i<l.size(); i++) {
        for (uint32_t j=0; j<size; j++){
            if (mon_dom_in(l[i],j, size, nb_elem)){

                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p1;
                p1 = poly_div(y, l[i], monomial_order, size, nb_elem);
                uint32_t m = p1.second.algebraic_degree(nb_elem);

                if (m<=2)   {
                // y = p1.first * l[i] + p1.second
                        
                    poly P1;
                    P1 = truncate_lin(p1.first, size);
                    poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);
                    poly P1_lin;
                    POLY_ADD(P1, p1.first, P1_lin);

                    vector<implem> v1;
                        
                    int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

                    if (indice1 != -1)   {
                        implem imp;
                        imp.op_sol.push_back(P1);
                        imp.op_sol.push_back(P1_lin);
                        imp.op_sol.push_back(l[i]);
                        imp.quad_sol.insert(pq1);
                        imp.formula = "(p1 ^ p2) * p3";
                        v1.push_back(imp);
                    }
                    else {
                        indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                        if (indice1 == -1)  {
                            //cout<<"On a des solutions pour v1 avec p1.first quadratique"<<endl;
                            continue;
                        }
                        else {
                            for (uint32_t m=0; m<10; m++)   {
                                implem imp;
                                poly_quad op_q1 = map_xor[indice1].second[m][0];
                                poly op1(op_q1,size);
                                imp.op_sol.push_back(op1);
                                imp.quad_sol.insert(op_q1);
                                poly_quad op_q2 = map_xor[indice1].second[m][1];
                                poly op2(op_q2,size);
                                imp.op_sol.push_back(op2);
                                imp.quad_sol.insert(op_q2);
                                imp.op_sol.push_back(P1_lin);
                                imp.op_sol.push_back(l[i]);
                                imp.formula = "";
                                v1.push_back(imp);
                            }
                        }
                    }

                    vector<implem> v2;

                    if (m==2)   {
                    // si p1.second est de deg 2

                        poly_quad pq2 = poly_to_poly_quad(p1.second ,size, nb_elem);
                            
                        int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);

                        if (indice2 != -1)   {
                            implem imp;
                            imp.op_sol.push_back(p1.second);
                            imp.quad_sol.insert(pq2);
                            imp.formula = "";
                            v2.push_back(imp);
                        }
                        else {
                            indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                            if (indice2 == -1)  {
                                //cout<<"On a des solutions pour v1 avec p1.second quadratique"<<endl;
                                continue;
                            }
                            else {
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                            }
                        }

                        for (uint32_t i1 = 0; i1<v1.size(); i1++) {

                            set<poly_quad> set_quad_1 = v1[i1].quad_sol;

                            for (uint32_t i2 = 0; i2<v2.size(); i2++) {

                                set<poly_quad> set_2 = v2[i2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);
                                set_quad_2.erase(0);

                                set<poly_quad> s_quad_temp = set_quad_2;
                                uint32_t r = Rank(s_quad_temp);

                                if (r<=MAX_RANK_DEG3_V1)   {
                                    implem temp;
                                    temp.quad_sol = set_quad_2;
                                    temp.op_sol.insert(temp.op_sol.end(),v1[i1].op_sol.begin(), v1[i1].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v2[i2].op_sol.begin(), v2[i2].op_sol.end());

                                    string f = "";
                                    if (v1.size()==1){
                                        if (v2.size()==1){
                                            f = "(p1 ^ p2) * (p3) ^ p4";
                                        }
                                        else {
                                            f = "(p1 ^ p2) * (p3) ^ p4 ^ p5";
                                        }
                                    }
                                    else {
                                        if (v2.size()==1){
                                            f = "(p1 ^ p2 ^ p3) * (p4) ^ p5";
                                        }
                                        else {
                                            f = "(p1 ^ p2 ^ p3) * (p4) ^ p5 ^ p6";
                                        }
                                    }
                                    temp.formula = f;
                                    res.push_back(temp);
                                }              
                            }
                        }
                    }
                    else {
                    // si p1.second est linéaire (ie = 0)
                        for (uint32_t i1 = 0; i1<v1.size(); i1++){
                            res.push_back(v1[i1]);
                        }
                    }
                }

                else {

                    uint32_t m1 = p1.first.algebraic_degree(nb_elem);
                    if (m1 == 0){
                        continue;
                    }
                    //y = p1.first * (p2.first + l[i]) + p2.second

                    for (uint32_t k=0; k<size; k++){

                        create_monomial_order(monomial_order, nb_elem, size, k);

                        pair<poly,poly> p2;
                        p2 = poly_div(p1.second, p1.first, monomial_order, size, nb_elem);

                        uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                        if (m2 <= 2){

                            poly P1;
                            P1 = truncate_lin(p1.first, size);
                            poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);
                            poly P1_lin;
                            POLY_ADD(P1, p1.first, P1_lin);

                            vector<implem> v1;
                        
                            int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

                            if (indice1 != -1)   {
                                implem imp;
                                imp.op_sol.push_back(P1);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);
                            }
                            else {
                                indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                if (indice1 == -1)  {
                                    //cout<<"On a des solutions pour v1 avec p1.first quadratique"<<endl;
                                    continue;
                                }
                                else {
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.formula = "";
                                        v1.push_back(imp);
                                    }
                                }
                            }

                            poly P;
                            POLY_ADD(p2.first, l[i], P);
                            poly P2 = truncate_lin(P, size);
                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
                            poly P2_lin;
                            POLY_ADD(P2, P, P2_lin);

                            vector<implem> v2;
                            
                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);

                            if (indice2 != -1)   {
                                implem imp;
                                imp.op_sol.push_back(P2);
                                imp.op_sol.push_back(P2_lin);
                                imp.quad_sol.insert(pq2);
                                imp.formula = "";
                                v2.push_back(imp);
                            }
                            else {
                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                if (indice2 == -1)  {
                                    continue;
                                }
                                else {
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                }
                            }

                            if (m2==2)   {
                            // si p1.second est de deg 2
                                poly P3 = truncate_lin(p2.second, size);
                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                vector<implem> v3;
                                    
                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);

                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(p2.second);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
                                }

                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;
            
                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                        set<poly_quad> set_quad_2 = set_quad_1;
                                        set_quad_2.merge(set_2);
            
                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                            set<poly_quad> set_quad_3 = set_quad_2;
                                            set_quad_3.merge(set_3);
                                            set_quad_3.erase(0);
            
                                            set<poly_quad> s_quad_temp = set_quad_3;
                                            uint32_t r = Rank(s_quad_temp);
            
                                            if (r<=MAX_RANK_DEG3_V1)   {
                                                implem temp; 
                                                temp.quad_sol = set_quad_3;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());

                                                string f = "";

                                                if (v1.size()==1){
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6 ^ p7";
                                                        }
                                                    }
                                                }
                                                else {
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7 ^ p8";
                                                        }
                                                    }
                                                }
                                                
                                                temp.formula = f;
                                                res.push_back(temp);  
                                            }            
                                        }
                                    }
                                }    
                            }

                            else {
                            // si p2.second est linéaire (ie = 0)
                            
                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                        set<poly_quad> set_quad_2 = set_quad_1;
                                        set_quad_2.merge(set_2);
                                        set_quad_2.erase(0);
    
                                        set<poly_quad> s_quad_temp = set_quad_2;
                                        uint32_t r = Rank(s_quad_temp);
    
                                        if (r<=MAX_RANK_DEG3_V1)   {
                                            implem temp; 
                                            temp.quad_sol = set_quad_2;
                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());

                                            string f = "";

                                            if (v1.size()==1){
                                                if (v2.size()==1){
                                                    f = "(p1 ^ p2) * (p3 ^ p4)";
                                                }
                                                else {
                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5)";           
                                                }
                                            }
                                            else {
                                                if (v2.size()==1){
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5)";
                                                }
                                                else {
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6)";
                                                }
                                            }

                                            temp.formula = f;
                                            res.push_back(temp); 
                                        }              
                                    }
                                }       
                            }
                        }
                    }
                }
            }
        }
    }

    //if (res.size() == 0){
        for (uint32_t i=0; i<set_op_l_plus_lin.size(); i++) {

            for (uint32_t j=0; j<size; j++){
    
                create_monomial_order(monomial_order, nb_elem, size, j);
    
                pair<poly,poly> p;
                p = poly_div(y, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
                uint32_t m = p.second.algebraic_degree(nb_elem);
                uint32_t n = p.first.algebraic_degree(nb_elem);
    
                if (m<=2 && n<=2)   {

                // y = p.first * s[i] + p.second
    
                    vector<implem> v1;
    
                    poly P1 = truncate_lin(set_op_l_plus_lin[i], size);
                    poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem);           
                    implem imp;
                    imp.op_sol.push_back(P1);
                    poly P1_lin;
                    POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);
                    imp.op_sol.push_back(P1_lin);
                    imp.quad_sol.insert(pq1);
                    imp.formula = "";
                    v1.push_back(imp);                
    
                    vector<implem> v2;
    
                    poly P2;
                    P2 = truncate_lin(p.first, size);
                    poly P2_lin;
                    POLY_ADD(P2, p.first, P2_lin);
                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
    
                    int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                    if (indice2 != -1)   {
                        implem imp;
                        imp.op_sol.push_back(P2);
                        imp.op_sol.push_back(P2_lin);
                        imp.quad_sol.insert(pq2);
                        imp.formula = "";
                        v2.push_back(imp);
                    }
                    else {
                        indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                        if (indice2 == -1)  {
                            continue;
                        }
                        else {
                            for (uint32_t m=0; m<10; m++)   {
                                implem imp;
                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                poly op1(op_q1,size);
                                imp.op_sol.push_back(op1);
                                imp.quad_sol.insert(op_q1);
                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                poly op2(op_q2,size);
                                imp.op_sol.push_back(op2);
                                imp.quad_sol.insert(op_q2);
                                imp.op_sol.push_back(P2_lin);
                                imp.formula = "";
                                v2.push_back(imp);
                            }
                        }
                    }
    
                    if (m == 0){
                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);
                                set_quad_2.erase(0);
    
                                set<poly_quad> s_quad_temp = set_quad_2;
                                uint32_t r = Rank(s_quad_temp);
    
                                if (r<=MAX_RANK_DEG3_V1)   {
                                    implem temp; 
                                    temp.quad_sol = set_quad_2;
                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                    string f = "";
                                    if (v2.size()==1){
                                        f = "(p1 ^ p2) * (p3 ^ p4)";
                                    }
                                    else {
                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5)";
                                    }
                                    temp.formula = f;
                                    res.push_back(temp);
                                }
                            }
                        }
                    }
                    else {
    
                        vector<implem> v3;
    
                        poly P3 = truncate_lin(p.second, size);
                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);
    
                        int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                        if (indice3 != -1)   {
                            implem imp;
                            imp.op_sol.push_back(p.second);
                            imp.quad_sol.insert(pq3);
                            imp.formula = "";
                            v3.push_back(imp);
                        }
                        else {
                            indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                            if (indice3 == -1)  {
                                continue;
                            }
                            else {
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                            }
                        }                           
    
                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);
    
                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                    set<poly_quad> set_quad_3 = set_quad_2;
                                    set_quad_3.merge(set_3);
                                    set_quad_3.erase(0);
    
                                    set<poly_quad> s_quad_temp = set_quad_3;
                                    uint32_t r = Rank(s_quad_temp);
    
                                    if (r<=MAX_RANK_DEG3_V1)   {
                                        implem temp; 
                                        temp.quad_sol = set_quad_3;
                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                        
                                        string f = "";
                                        if (v2.size()==1){
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ p5";
                                            }
                                            else {
                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6";
                                            }
                                        }
                                        else {
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6";
                                            }
                                            else {
                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6 ^ p7";
                                            }
                                        }
                                        temp.formula = f;
                                        res.push_back(temp);
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    if (n!=0 && n<=2){
                        for (uint32_t jj=0; jj<size; jj++){
    
                            if (jj == j){
                                continue;
                            }
                            
                            create_monomial_order(monomial_order, nb_elem, size, jj);
    
                            pair<poly,poly> pp;
                            pp = poly_div(p.second, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
                            uint32_t mm = pp.second.algebraic_degree(nb_elem);
                            uint32_t nn = pp.first.algebraic_degree(nb_elem);
    
                            if (mm <=2 && nn<=2)   {
                            // y = (p.first + pp.first) * s[i] + pp.second
    
                                vector<implem> v1;
    
                                poly P1 = truncate_lin(set_op_l_plus_lin[i], size);
                                poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem); 
                                implem imp;
                                imp.op_sol.push_back(P1);
                                poly P1_lin;
                                POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);                
    
                                vector<implem> v2;
    
                                poly PP2;
                                POLY_ADD(p.first, pp.first, PP2);
                                poly P2 = truncate_lin(PP2, size);
                                poly P2_lin;
                                POLY_ADD(P2, PP2, P2_lin);
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
    
                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                    if (indice2 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
                                }
    
                                if (mm == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);
                                            set_quad_2.erase(0);
    
                                            set<poly_quad> s_quad_temp = set_quad_2;
                                            uint32_t r = Rank(s_quad_temp);
    
                                            if (r<=MAX_RANK_DEG3_V1)   {
                                                implem temp; 
                                                temp.quad_sol = set_quad_2;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                string f = "";
                                                if (v2.size()==1){
                                                    f = "(p1 ^ p2) * (p3 ^ p4)";
                                                }
                                                else {
                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5)";
                                                }
                                                temp.formula = f;
                                                res.push_back(temp); 
                                            }
                                        }
                                    }
                                }
                                else {
    
                                    vector<implem> v3;
    
                                    poly P3;
                                    P3 = truncate_lin(pp.second, size);
                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);
    
                                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                    if (indice3 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(pp.second);
                                        imp.quad_sol.insert(pq3);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                    else {
                                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                        if (indice3 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                        }
                                    }                           
    
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);
    
                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);
    
                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);
    
                                                if (r<=MAX_RANK_DEG3_V1)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6 ^ p7";
                                                        }
                                                    }
                                                    temp.formula = f;
                                                    res.push_back(temp); 
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    //}

    if (res.size() == 0){
        vector<poly> map_xor_l_plus_lin = create_map_xor_l_plus_lin(size);

        for (uint32_t i=0; i<map_xor_l_plus_lin.size(); i++) {
            for (uint32_t j=0; j<size; j++){

                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p;
                p = poly_div(y, map_xor_l_plus_lin[i], monomial_order, size, nb_elem);
                uint32_t m = p.second.algebraic_degree(nb_elem);
                uint32_t n = p.first.algebraic_degree(nb_elem);

                if (m<=2 && n<=2)   {
                // y = p.first * mx[i] + p.second

                    vector<implem> v1;

                    poly P1 = truncate_lin(map_xor_l_plus_lin[i], size);
                    poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem);  
                    int32_t indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                    for (uint32_t m=0; m<10; m++)   {
                        implem imp;
                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                        poly op1(op_q1,size);
                        imp.op_sol.push_back(op1);
                        imp.quad_sol.insert(op_q1);
                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                        poly op2(op_q2,size);
                        imp.op_sol.push_back(op2);
                        imp.quad_sol.insert(op_q2);
                        poly D; 
                        POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                        imp.op_sol.push_back(D);
                        imp.formula = "";
                        v1.push_back(imp);
                    }              

                    vector<implem> v2;

                    poly P2;
                    P2 = truncate_lin(p.first, size);
                    poly P2_lin;
                    POLY_ADD(P2, p.first, P2_lin);
                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                    int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                    if (indice2 != -1)   {
                        implem imp;
                        imp.op_sol.push_back(P2);
                        imp.op_sol.push_back(P2_lin);
                        imp.quad_sol.insert(pq2);
                        imp.formula = "";
                        v2.push_back(imp);
                    }
                    else {
                        indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                        if (indice2 == -1)  {
                            continue;
                        }
                        else {
                            for (uint32_t m=0; m<10; m++)   {
                                implem imp;
                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                poly op1(op_q1,size);
                                imp.op_sol.push_back(op1);
                                imp.quad_sol.insert(op_q1);
                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                poly op2(op_q2,size);
                                imp.op_sol.push_back(op2);
                                imp.quad_sol.insert(op_q2);
                                imp.op_sol.push_back(P2_lin);
                                imp.formula = "";
                                v2.push_back(imp);
                            }
                        }
                    }

                    if (m == 0){
                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);
                                set_quad_2.erase(0);

                                set<poly_quad> s_quad_temp = set_quad_2;
                                uint32_t r = Rank(s_quad_temp);

                                if (r<=MAX_RANK_DEG3_V1)   {
                                    implem temp; 
                                    temp.quad_sol = set_quad_2;
                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                    string f = "";
                                    if (v2.size()==1){
                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5)";
                                    }
                                    else {
                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6)";
                                    }
                                    temp.formula = f;
                                    res.push_back(temp); 
                                }
                            }
                        }
                    }
                    else {

                        vector<implem> v3;

                        poly P3 = truncate_lin(p.second, size);
                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                        int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                        if (indice3 != -1)   {
                            implem imp;
                            imp.op_sol.push_back(p.second);
                            imp.quad_sol.insert(pq3);
                            imp.formula = "";
                            v3.push_back(imp);
                        }
                        else {
                            indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                            if (indice3 == -1)  {
                                continue;
                            }
                            else {
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                            }
                        }                           

                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);

                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                    set<poly_quad> set_quad_3 = set_quad_2;
                                    set_quad_3.merge(set_3);
                                    set_quad_3.erase(0);

                                    set<poly_quad> s_quad_temp = set_quad_3;
                                    uint32_t r = Rank(s_quad_temp);

                                    if (r<=MAX_RANK_DEG3_V1)   {
                                        implem temp; 
                                        temp.quad_sol = set_quad_3;
                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                        string f = "";
                                        if (v2.size()==1){
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6";
                                            }
                                            else {
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7";
                                            }
                                        }
                                        else {
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7";
                                            }
                                            else {
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7 ^ p8";
                                            }
                                        }
                                        temp.formula = f;
                                        res.push_back(temp); 
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    if (n!=0 && n<=2){
                        for (uint32_t jj=0; jj<size; jj++){

                            if (jj == j){
                                continue;
                            }
                            
                            create_monomial_order(monomial_order, nb_elem, size, jj);

                            pair<poly,poly> pp;
                            pp = poly_div(p.second, map_xor_l_plus_lin[i], monomial_order, size, nb_elem);
                            uint32_t mm = pp.second.algebraic_degree(nb_elem);
                            uint32_t nn = pp.first.algebraic_degree(nb_elem);

                            if (mm <=2 && nn<=2)   {
                            // y = (p.first + pp.first) * s[i] + pp.second

                                vector<implem> v1;

                                poly P1 = truncate_lin(map_xor_l_plus_lin[i], size);
                                poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem); 
                                int32_t indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice1].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice1].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    poly D; 
                                    POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                                    imp.op_sol.push_back(D);
                                    imp.formula = "";
                                    v1.push_back(imp);
                                }                 

                                vector<implem> v2;

                                poly PP2;
                                POLY_ADD(p.first, pp.first, PP2);
                                poly P2 = truncate_lin(PP2, size);
                                poly P2_lin;
                                POLY_ADD(PP2, P2, P2_lin);
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                    if (indice2 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
                                }

                                if (mm == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);
                                            set_quad_2.erase(0);

                                            set<poly_quad> s_quad_temp = set_quad_2;
                                            uint32_t r = Rank(s_quad_temp);

                                            if (r<=MAX_RANK_DEG3_V1)   {
                                                implem temp; 
                                                temp.quad_sol = set_quad_2;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                string f = "";
                                                if (v2.size()==1){
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5)";
                                                }
                                                else {
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6)";
                                                }
                                                temp.formula = f;
                                                res.push_back(temp);  
                                            }
                                        }
                                    }
                                }
                                else {

                                    vector<implem> v3;

                                    poly P3;
                                    P3 = truncate_lin(pp.second, size);
                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                    if (indice3 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(pp.second);
                                        imp.quad_sol.insert(pq3);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                    else {
                                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                        if (indice3 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                        }
                                    }                           

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);

                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG3_V1)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7 ^ p8";
                                                        }
                                                    }
                                                    temp.formula = f;
                                                    res.push_back(temp);  
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    /*cout<<"Fin deg 3 v1 pour thread : "<<omp_get_thread_num()<<endl;
    res = remove_duplicates(res);
    cout<<" Thread : "<<omp_get_thread_num()<<" res_size v1 = "<<res.size()<<endl;*/
    res = remove_duplicates(res);

    return res;
}

vector<implem> add_to_op_selec_deg_3_v2(poly y, vector<poly> l, vector<poly> l2, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    for (uint32_t i=0; i<l.size(); i++) {
        for (uint32_t j=0; j<size; j++){
            if (mon_dom_in(l[i],j, size, nb_elem)){
                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p1;
                p1 = poly_div(y, l[i], monomial_order, size, nb_elem);
                uint32_t m = p1.second.algebraic_degree(nb_elem);
                uint32_t n = p1.first.algebraic_degree(nb_elem);

                if (m<=2 || n<2)    { // Cas traité dans la v1
                    continue;
                }
                else {
                // On traite le cas où y = p1.first * l[i] + p2.first * l[i2] + p2.second

                    for (uint32_t i2=0; i2<l.size(); i2++) {
                        for (uint32_t j2=0; j2<size; j2++){

                            if (mon_dom_in(l[i2],j2, size, nb_elem)){
                                create_monomial_order(monomial_order, nb_elem, size, j2);

                                pair<poly,poly> p2;
                                p2 = poly_div(p1.second, l[i2], monomial_order, size, nb_elem);

                                uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                
                                if (m2 <= 2){

                                    poly P1;
                                    P1 = truncate_lin(p1.first, size);
                                    poly P1_lin;
                                    POLY_ADD(P1, p1.first, P1_lin);
                                    poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                                    vector<implem> v1;
                                
                                    int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

                                    if (indice1 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(P1);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.op_sol.push_back(l[i]);
                                        imp.quad_sol.insert(pq1);
                                        imp.formula = "";
                                        v1.push_back(imp);
                                    }
                                    else {
                                        indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                        if (indice1 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice1].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice1].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.op_sol.push_back(P1_lin);
                                                imp.op_sol.push_back(l[i]);
                                                imp.formula = "";
                                                v1.push_back(imp);
                                            }
                                        }
                                    }

                                    poly P2 = truncate_lin(p2.first, size);
                                    poly P2_lin;
                                    POLY_ADD(P2, p2.first, P2_lin);
                                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                    vector<implem> v2;
                                    
                                    int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);

                                    if (indice2 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(P2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.op_sol.push_back(l[i2]);
                                        imp.quad_sol.insert(pq2);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                    else {
                                        indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                        if (indice2 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.formula = "";
                                                
                                                v2.push_back(imp);
                                            }
                                        }
                                    }

                                    if (m2==2)   {
                                    // si p1.second est de deg 2
                                        poly P3 = truncate_lin(p2.second, size);
                                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                        vector<implem> v3;
                                            
                                        int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);

                                        if (indice3 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(p2.second);
                                            imp.quad_sol.insert(pq3);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                        else {
                                            indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                            if (indice3 == -1)  {
                                                continue;
                                            }
                                            else {
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.formula = ""; 
                                                    v3.push_back(imp);
                                                }
                                            }
                                        }

                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                    implem temp;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    set<poly_quad> set_temp = v1[it1].quad_sol;
                                                    temp.quad_sol.merge(set_temp);

                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    set_temp = v2[it2].quad_sol;
                                                    temp.quad_sol.merge(set_temp);

                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    set_temp = v3[it3].quad_sol;
                                                    temp.quad_sol.merge(set_temp);

                                                    set<poly_quad> s_quad_temp = temp.quad_sol;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG3_V2)   {
                                                        string f = "";

                                                        if (v1.size()==1){
                                                            if (v2.size()==1){
                                                                if (v3.size()==1){
                                                                    f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5) * (p6) ^ p7";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5) * (p6) ^ p7 ^ p8";
                                                                }
                                                            }
                                                            else {
                                                                if (v3.size()==1){
                                                                    f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5 ^ p6) * (p7) ^ p8";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5 ^ p6) * (p7) ^ p8 ^ p9";
                                                                }
                                                            }
                                                        }
                                                        else {
                                                            if (v2.size()==1){
                                                                if (v3.size()==1){
                                                                    f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6) * (p7) ^ p8";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6) * (p7) ^ p8 ^ p9";
                                                                }
                                                            }
                                                            else {
                                                                if (v3.size()==1){
                                                                    f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                }
                                                            }
                                                        }
                                                        
                                                        temp.formula = f;
                                                        res.push_back(temp);  
                                                    }            
                                                }
                                            }
                                        }    
                                    }

                                    else {
                                    // si p2.second est linéaire (ie = 0)
                                    
                                        set<set<poly_quad>> set_quad;
                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                implem temp;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                set<poly_quad> set_temp = v1[it1].quad_sol;
                                                temp.quad_sol.merge(set_temp);

                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                set_temp = v2[it2].quad_sol;
                                                temp.quad_sol.merge(set_temp);

                                                set<poly_quad> s_quad_temp = temp.quad_sol;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG3_V2)   {
                                                    string f = "";

                                                    if (v1.size()==1){
                                                        if (v2.size()==1){ 
                                                            f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5) * (p6)";
                                                        }
                                                        else {
                                                            
                                                            f = "(p1 ^ p2) * (p3) ^ (p4 ^ p5 ^ p6) * (p7)";
                                                        }
                                                    }
                                                    else {
                                                        if (v2.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6) * (p7)";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4) ^ (p5 ^ p6 ^ p7) * (p8)";
                                                        }
                                                    }
                                                        
                                                    temp.formula = f;
                                                    res.push_back(temp);
                                                }             
                                            }
                                        }       
                                    }
                                }

                                /*else {// p2.second est de deg 3, on essaye de le diviser par p1.first ou p2.first

                                    for (uint32_t k=0; k<size; k++){

                                        create_monomial_order(monomial_order, nb_elem, size, k);

                                        pair<poly,poly> p3;
                                        p3 = poly_div(p2.second, p1.first, monomial_order, size, nb_elem);

                                        uint32_t m3 = p3.second.algebraic_degree(nb_elem);

                                        if (m3 > 2){

                                            p3 = poly_div(p2.second, p2.first, monomial_order, size, nb_elem);

                                            uint32_t m4 = p3.second.algebraic_degree(nb_elem);

                                            if (m4 > 2){
                                                continue;
                                            }
                                            else {
                                            // On traite le cas où y = p1.first * l[i] + p2.first * (l[i2] + p3.first) + p3.second
                                                poly P1;
                                                P1 = truncate_lin(p1.first, size);
                                                poly P1_lin;
                                                POLY_ADD(P1, p1.first, P1_lin);
                                                poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                                                vector<implem> v1;
                            
                                                int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

                                                if (indice1 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(P1);
                                                    imp.op_sol.push_back(P1_lin);
                                                    imp.op_sol.push_back(l[i]);
                                                    imp.quad_sol.insert(pq1);
                                                    imp.formula = "";
                                                    v1.push_back(imp);
                                                }
                                                else {
                                                    indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                                    if (indice1 == -1)  {
                                                        continue;
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice1].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice1].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P1_lin);
                                                            imp.op_sol.push_back(l[i]);
                                                            imp.formula = "";
                                                            v1.push_back(imp);
                                                        }
                                                    }
                                                }

                                                poly P2;
                                                P2 = truncate_lin(p2.first, size);
                                                poly P2_lin;
                                                POLY_ADD(P2, p2.first, P2_lin);
                                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                                vector<implem> v2;
                            
                                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);

                                                if (indice2 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(P2);
                                                    imp.op_sol.push_back(P2_lin);
                                                    imp.quad_sol.insert(pq2);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                                else {
                                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                    if (indice2 == -1)  {
                                                        continue;
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                    }
                                                }

                                                poly p;
                                                POLY_ADD(l[i2], p3.first, p);
                                                poly P3;
                                                P3 = truncate_lin(p, size);
                                                poly P3_lin;
                                                POLY_ADD(P3, p, P3_lin);
                                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                                vector<implem> v3;
                            
                                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);

                                                if (indice3 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(P3);
                                                    imp.op_sol.push_back(P3_lin);
                                                    imp.quad_sol.insert(pq3);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                                else {
                                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                    if (indice3 == -1)  {
                                                        continue;
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P3_lin);
                                                            imp.formula = "";
                                                            v3.push_back(imp);
                                                        }
                                                    }
                                                }

                                                if (m4 == 0){

                                                }

                                                else {

                                                    poly P4;
                                                    P4 = truncate_lin(p3.second, size,);
                                                    poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                                    vector<implem> v4;
                                
                                                    int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);

                                                    if (indice4 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(p3.second);
                                                        imp.quad_sol.insert(pq4);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                    else {
                                                        indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                        if (indice4 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                poly D; 
                                                                POLY_ADD(P4, p3.second, D); //Pour récupérer le delta linéaire
                                                                imp.op_sol.push_back(D);
                                                                imp.formula = "";
                                                                v4.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                    implem temp;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    set<poly_quad> set_temp = v1[it1].quad_sol;
                                                                    temp.quad_sol.merge(set_temp);

                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    set_temp = v2[it2].quad_sol;
                                                                    temp.quad_sol.merge(set_temp);

                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    set_temp = v3[it3].quad_sol;
                                                                    temp.quad_sol.merge(set_temp);

                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    set_temp = v4[it4].quad_sol;
                                                                    temp.quad_sol.merge(set_temp);

                                                                    set<poly_quad> s_quad_temp = temp.quad_sol;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG3_V2)   {
                                                                        temp.formula = "323";
                                                                        res.push_back(temp);  
                                                                    }  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {
                                            
                                            // On traite le cas où y = p1.first * (l[i] + p3.first) + p2.first * l[i2] + p3.second
                                            poly P1;
                                            P1 = truncate_lin(p1.first, size);
                                            poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                                            vector<implem> v1;
                            
                                            int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);

                                            if (indice1 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p1.first);
                                                imp.quad_sol.insert(pq1);
                                                imp.formula = "";
                                                v1.push_back(imp);
                                            }
                                            else {
                                                indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                                if (indice1 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        poly D; 
                                                        POLY_ADD(P1, p1.first, D); //Pour récupérer le delta linéaire
                                                        imp.op_sol.push_back(D);
                                                        imp.formula = "";
                                                        v1.push_back(imp);
                                                    }
                                                }
                                            }

                                            poly p;
                                            POLY_ADD(l[i], p3.first, p);
                                            poly P2;
                                            P2 = truncate_lin(p, size);
                                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                            vector<implem> v2;
                            
                                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);

                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                if (indice2 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        poly D; 
                                                        POLY_ADD(P2, p, D); //Pour récupérer le delta linéaire
                                                        imp.op_sol.push_back(D);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                }
                                            }

                                            poly P3;
                                            P3 = truncate_lin(p2.first, size);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            vector<implem> v3;
                        
                                            int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);

                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p2.first);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                if (indice3 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.op_sol.push_back(l[i2]);
                                                        imp.quad_sol.insert(op_q2);
                                                        poly D; 
                                                        POLY_ADD(P3, p2.first, D); //Pour récupérer le delta linéaire
                                                        imp.op_sol.push_back(D);
                                                        imp.formula = "";
                                                        v3.push_back(imp);
                                                    }
                                                }
                                            }

                                            poly P4;
                                            P4 = truncate_lin(p3.second, size);
                                            poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                            vector<implem> v4;
                        
                                            int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);

                                            if (indice4 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p3.second);
                                                imp.quad_sol.insert(pq4);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                            else {
                                                indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                if (indice4 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                }
                                            }

                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            implem temp;
                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                            set<poly_quad> set_temp = v1[it1].quad_sol;
                                                            temp.quad_sol.merge(set_temp);

                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                            set_temp = v2[it2].quad_sol;
                                                            temp.quad_sol.merge(set_temp);

                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                            set_temp = v3[it3].quad_sol;
                                                            temp.quad_sol.merge(set_temp);

                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                            set_temp = v4[it4].quad_sol;
                                                            temp.quad_sol.merge(set_temp);

                                                            set<poly_quad> s_quad_temp = temp.quad_sol;
                                                            uint32_t r = Rank(s_quad_temp);

                                                            if (r<=MAX_RANK_DEG3_V2)   {
                                                                temp.formula = "324";
                                                                res.push_back(temp);
                                                            }    
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }*/
                            }
                        }
                    }
                }
            }
        }
    }

    /*cout<<"Fin deg 3 v2 pour thread : "<<omp_get_thread_num()<<endl;
    res = remove_duplicates(res);
    cout<<" Thread : "<<omp_get_thread_num()<<" res_size v2 = "<<res.size()<<endl;*/

    return remove_duplicates(res);
}

//DEG 4

vector<implem> add_to_op_selec_deg_4_v1(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    poly y4 = truncate_except_deg(y, 4, size, nb_elem);
    poly y3 = truncate_except_deg(y, 3, size, nb_elem);
    poly y2 = truncate_except_deg(y, 2, size, nb_elem);

    poly y4_3;
    POLY_ADD(y4, y3, y4_3);

    for (uint32_t i=0; i<set_op_l_plus_lin.size(); i++) {
        for (uint32_t j=0; j<size; j++){

            create_monomial_order(monomial_order, nb_elem, size, j);

            pair<poly,poly> p;
            p = poly_div(y4_3, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
            uint32_t m = p.second.algebraic_degree(nb_elem);
            uint32_t n = p.first.algebraic_degree(nb_elem);

            if (m<=2 && n<=2)   {
            // y = p.first * s[i] + p.second + y2

                poly P3;
                POLY_ADD(p.second, y2, P3);
                m = P3.algebraic_degree(nb_elem);

                vector<implem> v1;

                poly P1 = truncate_lin(set_op_l_plus_lin[i], size);
                poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem);
                poly P1_lin;
                POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);           
                implem imp;
                imp.op_sol.push_back(P1);
                imp.op_sol.push_back(P1_lin);
                imp.quad_sol.insert(pq1);
                imp.formula = "";
                v1.push_back(imp);                

                vector<implem> v2;

                poly P2;
                P2 = truncate_lin(p.first, size);
                poly P2_lin;
                POLY_ADD(P2, p.first, P2_lin);
                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                if (indice2 != -1)   {
                    implem imp;
                    imp.op_sol.push_back(P2);
                    imp.op_sol.push_back(P2_lin);
                    imp.quad_sol.insert(pq2);
                    imp.formula = "";
                    v2.push_back(imp);
                }
                else {
                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                    if (indice2 == -1)  {
                        continue;
                    }
                    else {
                        for (uint32_t m=0; m<10; m++)   {
                            implem imp;
                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                            poly op1(op_q1,size);
                            imp.op_sol.push_back(op1);
                            imp.quad_sol.insert(op_q1);
                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                            poly op2(op_q2,size);
                            imp.op_sol.push_back(op2);
                            imp.quad_sol.insert(op_q2);
                            imp.op_sol.push_back(P2_lin);
                            imp.formula = "";
                            v2.push_back(imp);
                        }
                    }
                }

                if (m == 0){
                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                            set<poly_quad> set_2 = v2[it2].quad_sol;
                            set<poly_quad> set_quad_2 = set_quad_1;
                            set_quad_2.merge(set_2);
                            set_quad_2.erase(0);

                            set<poly_quad> s_quad_temp = set_quad_2;
                            uint32_t r = Rank(s_quad_temp);

                            if (r<=MAX_RANK_DEG4_V1)   {
                                implem temp; 
                                temp.quad_sol = set_quad_2;
                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                string f = "";
                                if (v2.size()==1){
                                    f = "(p1 ^ p2) * (p3 ^ p4)";
                                }
                                else {
                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5)";
                                }
                                temp.formula = f;
                                res.push_back(temp);  
                            }
                        }
                    }
                }
                else {

                    vector<implem> v3;

                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                    if (indice3 != -1)   {
                        implem imp;
                        imp.op_sol.push_back(P3);
                        imp.quad_sol.insert(pq3);
                        imp.formula = "";
                        v3.push_back(imp);
                    }
                    else {
                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                        if (indice3 == -1)  {
                            continue;
                        }
                        else {
                            for (uint32_t m=0; m<10; m++)   {
                                implem imp;
                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                poly op1(op_q1,size);
                                imp.op_sol.push_back(op1);
                                imp.quad_sol.insert(op_q1);
                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                poly op2(op_q2,size);
                                imp.op_sol.push_back(op2);
                                imp.quad_sol.insert(op_q2);
                                imp.formula = "";
                                v3.push_back(imp);
                            }
                        }
                    }                           

                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                            set<poly_quad> set_2 = v2[it2].quad_sol;
                            set<poly_quad> set_quad_2 = set_quad_1;
                            set_quad_2.merge(set_2);

                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                set<poly_quad> set_quad_3 = set_quad_2;
                                set_quad_3.merge(set_3);
                                set_quad_3.erase(0);

                                set<poly_quad> s_quad_temp = set_quad_3;
                                uint32_t r = Rank(s_quad_temp);

                                if (r<=MAX_RANK_DEG4_V1)   {
                                    implem temp; 
                                    temp.quad_sol = set_quad_3;
                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                    string f = "";
                                    if (v2.size()==1){
                                        if (v3.size()==1){
                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5";
                                        }
                                        else {
                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6";
                                        }
                                    }
                                    else {
                                        if (v3.size()==1){
                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6";
                                        }
                                        else {
                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6 ^ p7";
                                        }
                                    }
                                    temp.formula = f; 
                                    res.push_back(temp);
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (n!=0 && n<=2){
                    for (uint32_t jj=0; jj<size; jj++){

                        if (jj == j){
                            continue;
                        }
                        
                        create_monomial_order(monomial_order, nb_elem, size, jj);

                        pair<poly,poly> pp;
                        pp = poly_div(p.second, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
                        uint32_t mm = pp.second.algebraic_degree(nb_elem);
                        uint32_t nn = pp.first.algebraic_degree(nb_elem);

                        if (mm <=2 && nn<=2)   {
                        // y = (p.first + pp.first) * s[i] + pp.second
                            poly P3;
                            POLY_ADD(pp.second, y2, P3);
                            mm = P3.algebraic_degree(nb_elem);

                            vector<implem> v1;

                            poly P1 = truncate_lin(set_op_l_plus_lin[i], size);
                            poly P1_lin;
                            POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);
                            poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem); 
                            implem imp;
                            imp.op_sol.push_back(P1);
                            imp.op_sol.push_back(P1_lin);
                            imp.quad_sol.insert(pq1);
                            imp.formula = "";
                            v1.push_back(imp);                

                            vector<implem> v2;

                            poly PP2;
                            POLY_ADD(p.first, pp.first, PP2);
                            poly P2 = truncate_lin(PP2, size);
                            poly P2_lin;
                            POLY_ADD(P2, PP2, P2_lin);
                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                            if (indice2 != -1)   {
                                implem imp;
                                imp.op_sol.push_back(P2);
                                imp.op_sol.push_back(P2_lin);
                                imp.quad_sol.insert(pq2);
                                imp.formula = "";
                                v2.push_back(imp);
                            }
                            else {
                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                if (indice2 == -1)  {
                                    continue;
                                }
                                else {
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                }
                            }

                            if (mm == 0){
                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                        set<poly_quad> set_quad_2 = set_quad_1;
                                        set_quad_2.merge(set_2);
                                        set_quad_2.erase(0);

                                        set<poly_quad> s_quad_temp = set_quad_2;
                                        uint32_t r = Rank(s_quad_temp);

                                        if (r<=MAX_RANK_DEG4_V1)   {
                                            implem temp; 
                                            temp.quad_sol = set_quad_2;
                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                            string f = "";
                                            if (v2.size()==1){
                                                f = "(p1 ^ p2) * (p3 ^ p4)";
                                            }
                                            else {
                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5)";
                                            }
                                            temp.formula = f;
                                            res.push_back(temp);   
                                        }
                                    }
                                }
                            }
                            else {

                                vector<implem> v3;

                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P3);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
                                }                           

                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                        set<poly_quad> set_quad_2 = set_quad_1;
                                        set_quad_2.merge(set_2);

                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                            set<poly_quad> set_quad_3 = set_quad_2;
                                            set_quad_3.merge(set_3);
                                            set_quad_3.erase(0);

                                            set<poly_quad> s_quad_temp = set_quad_3;
                                            uint32_t r = Rank(s_quad_temp);

                                            if (r<=MAX_RANK_DEG4_V1)   {
                                                implem temp; 
                                                temp.quad_sol = set_quad_3;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                string f = "";
                                                if (v2.size()==1){
                                                    if (v3.size()==1){
                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ p5";
                                                    }
                                                    else {
                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6";
                                                    }
                                                }
                                                else {
                                                    if (v3.size()==1){
                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6";
                                                    }
                                                    else {
                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ p6 ^ p7";
                                                    }
                                                }
                                                temp.formula = f; 
                                                res.push_back(temp); 
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (res.size() == 0){

        vector<poly> map_xor_l_plus_lin = create_map_xor_l_plus_lin(size);

        for (uint32_t i=0; i<map_xor_l_plus_lin.size(); i++) {
            for (uint32_t j=0; j<size; j++){

                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p;
                p = poly_div(y4_3, map_xor_l_plus_lin[i], monomial_order, size, nb_elem);
                uint32_t m = p.second.algebraic_degree(nb_elem);
                uint32_t n = p.first.algebraic_degree(nb_elem);

                if (m<=2 && n<=2)   {
                // y = p.first * mx[i] + p.second + y2

                    poly P3;
                    POLY_ADD(p.second, y2, P3);
                    m = P3.algebraic_degree(nb_elem);

                    vector<implem> v1;

                    poly P1 = truncate_lin(map_xor_l_plus_lin[i], size);
                    poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem);  
                    int32_t indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                    for (uint32_t m=0; m<10; m++)   {
                        implem imp;
                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                        poly op1(op_q1,size);
                        imp.op_sol.push_back(op1);
                        imp.quad_sol.insert(op_q1);
                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                        poly op2(op_q2,size);
                        imp.op_sol.push_back(op2);
                        imp.quad_sol.insert(op_q2);
                        poly D; 
                        POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                        imp.op_sol.push_back(D);
                        imp.formula = "";
                        v1.push_back(imp);
                    }              

                    vector<implem> v2;

                    poly P2;
                    P2 = truncate_lin(p.first, size);
                    poly P2_lin;
                    POLY_ADD(P2, p.first, P2_lin);
                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                    int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                    if (indice2 != -1)   {
                        implem imp;
                        imp.op_sol.push_back(P2);
                        imp.op_sol.push_back(P2_lin);
                        imp.quad_sol.insert(pq2);
                        imp.formula = "";
                        v2.push_back(imp);
                    }
                    else {
                        indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                        if (indice2 == -1)  {
                            continue;
                        }
                        else {
                            for (uint32_t m=0; m<10; m++)   {
                                implem imp;
                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                poly op1(op_q1,size);
                                imp.op_sol.push_back(op1);
                                imp.quad_sol.insert(op_q1);
                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                poly op2(op_q2,size);
                                imp.op_sol.push_back(op2);
                                imp.quad_sol.insert(op_q2);
                                imp.op_sol.push_back(P2_lin);
                                imp.formula = "";
                                v2.push_back(imp);
                            }
                        }
                    }

                    if (m == 0){
                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);
                                set_quad_2.erase(0);

                                set<poly_quad> s_quad_temp = set_quad_2;
                                uint32_t r = Rank(s_quad_temp);

                                if (r<=MAX_RANK_DEG4_V1)   {
                                    implem temp; 
                                    temp.quad_sol = set_quad_2;
                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                    string f = "";
                                    if (v2.size()==1){
                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5)";
                                    }
                                    else {
                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6)";
                                    }
                                    temp.formula = f;
                                    res.push_back(temp);
                                }
                            }
                        }
                    }
                    else {

                        vector<implem> v3;

                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                        int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                        if (indice3 != -1)   {
                            implem imp;
                            imp.op_sol.push_back(P3);
                            imp.quad_sol.insert(pq3);
                            imp.formula = "";
                            v3.push_back(imp);
                        }
                        else {
                            indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                            if (indice3 == -1)  {
                                continue;
                            }
                            else {
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                            }
                        }                           

                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                set<poly_quad> set_quad_2 = set_quad_1;
                                set_quad_2.merge(set_2);

                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                    set<poly_quad> set_quad_3 = set_quad_2;
                                    set_quad_3.merge(set_3);
                                    set_quad_3.erase(0);

                                    set<poly_quad> s_quad_temp = set_quad_3;
                                    uint32_t r = Rank(s_quad_temp);

                                    if (r<=MAX_RANK_DEG4_V1)   {
                                        implem temp; 
                                        temp.quad_sol = set_quad_3;
                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                        string f = "";
                                        if (v2.size()==1){
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6";
                                            }
                                            else {
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7";
                                            }
                                        }
                                        else {
                                            if (v3.size()==1){
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7";
                                            }
                                            else {
                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7 ^ p8";
                                            }
                                        }
                                        temp.formula = f;
                                        res.push_back(temp);
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    if (n!=0 && n<=2){
                        for (uint32_t jj=0; jj<size; jj++){

                            if (jj == j){
                                continue;
                            }
                            
                            create_monomial_order(monomial_order, nb_elem, size, jj);

                            pair<poly,poly> pp;
                            pp = poly_div(p.second, map_xor_l_plus_lin[i], monomial_order, size, nb_elem);
                            uint32_t mm = pp.second.algebraic_degree(nb_elem);
                            uint32_t nn = pp.first.algebraic_degree(nb_elem);

                            if (mm <=2 && nn<=2)   {
                            // y = (p.first + pp.first) * s[i] + pp.second
                                poly P3;
                                POLY_ADD(pp.second, y2, P3);
                                mm = P3.algebraic_degree(nb_elem);

                                vector<implem> v1;

                                poly P1 = truncate_lin(map_xor_l_plus_lin[i], size);
                                poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem); 
                                int32_t indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                for (uint32_t m=0; m<10; m++)   {
                                    implem imp;
                                    poly_quad op_q1 = map_xor[indice1].second[m][0];
                                    poly op1(op_q1,size);
                                    imp.op_sol.push_back(op1);
                                    imp.quad_sol.insert(op_q1);
                                    poly_quad op_q2 = map_xor[indice1].second[m][1];
                                    poly op2(op_q2,size);
                                    imp.op_sol.push_back(op2);
                                    imp.quad_sol.insert(op_q2);
                                    poly D; 
                                    POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                                    imp.op_sol.push_back(D);
                                    imp.formula = "";
                                    v1.push_back(imp);
                                }                 

                                vector<implem> v2;

                                poly PP2;
                                POLY_ADD(p.first, pp.first, PP2);
                                poly P2 = truncate_lin(PP2, size);
                                poly P2_lin;
                                POLY_ADD(P2, PP2, P2_lin);
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                    if (indice2 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
                                }

                                if (mm == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);
                                            set_quad_2.erase(0);

                                            set<poly_quad> s_quad_temp = set_quad_2;
                                            uint32_t r = Rank(s_quad_temp);

                                            if (r<=MAX_RANK_DEG4_V1)   {
                                                implem temp; 
                                                temp.quad_sol = set_quad_2;
                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                string f = "";
                                                if (v2.size()==1){
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5)";
                                                }
                                                else {
                                                    f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6)";
                                                }
                                                temp.formula = f;
                                                res.push_back(temp);  
                                            }
                                        }
                                    }
                                }
                                else {

                                    vector<implem> v3;

                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                    if (indice3 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(P3);
                                        imp.quad_sol.insert(pq3);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                    else {
                                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                        if (indice3 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                        }
                                    }                           

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);

                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG4_V1)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ p7 ^ p8";
                                                        }
                                                    }
                                                    temp.formula = f;
                                                    res.push_back(temp); 
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    if (verbose){
        cout<<"Fin deg 4 v1 pour thread : "<<omp_get_thread_num()<<endl;
        res = remove_duplicates(res);
        cout<<" Thread : "<<omp_get_thread_num()<<" res_size v1 = "<<res.size()<<endl;
    }

    return res;
}

vector<implem> add_to_op_selec_deg_4_v2(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)    {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    poly y4 = truncate_except_deg(y, 4, size, nb_elem);
    poly y3 = truncate_except_deg(y, 3, size, nb_elem);
    poly y2 = truncate_except_deg(y, 2, size, nb_elem);

    poly y3_2;
    POLY_ADD(y3, y2, y3_2);
    
    for (uint32_t i=0; i<set_op_l_plus_lin.size(); i++) {
        for (uint32_t j=0; j<size; j++){

            create_monomial_order(monomial_order, nb_elem, size, j);

            pair<poly,poly> p;
            p = poly_div(y4, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
            uint32_t m = p.second.algebraic_degree(nb_elem);
            uint32_t n = p.first.algebraic_degree(nb_elem);

            if (m<=3 && n<=2)   {
            //y = s[i] * p.first + p.second + y3_2

                poly P1;
                P1 = truncate_lin(set_op_l_plus_lin[i], size);
                poly P1_lin;
                POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);
                poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);
                
                poly P2 = truncate_lin(p.first, size);
                poly P2_lin;
                POLY_ADD(P2, p.first, P2_lin);
                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                int32_t indice_mP2 = -1;
                if (indiceP2 == -1){
                    indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                    if (indice_mP2 == -1){
                        continue;
                    }
                }

                bool test_ps = false;

                poly reste;
                POLY_ADD(p.second, y3_2, reste);

                for (uint32_t i2=0; i2<l.size(); i2++){
                    for (uint32_t j2=0; j2<size; j2++){
                        if (mon_dom_in(l[i2],j2, size, nb_elem)){

                            create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                            pair<poly,poly> p2;
                            p2 = poly_div(reste, l[i2], monomial_order, size, nb_elem);

                            uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                            uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                            if (m2 <= 2 && n2 <= 2){
                            //y = s[i] * p.first + p2.first * l[i2] + p2.second
                            bool test_ps = true;

                                vector<implem> v1;
                                                
                                implem imp;
                                imp.op_sol.push_back(P1);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);

                                vector<implem> v2;

                                int32_t indice2 = indiceP2;
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = indice_mP2;
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                }

                                vector<implem> v3;
                                poly P3;
                                P3 = truncate_lin(p2.first, size);
                                poly P3_lin;
                                POLY_ADD(P3, p2.first, P3_lin);
                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P3);
                                    imp.op_sol.push_back(P3_lin);
                                    imp.op_sol.push_back(l[i2]);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1){
                                        continue;
                                    }
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice3].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice3].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P3_lin);
                                        imp.op_sol.push_back(l[i2]);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                }

                                if (m2 == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);

                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG4_V2)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7)";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8)";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8)";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9)";
                                                        }
                                                    }
                                                      
                                                    temp.formula = f;
                                                    res.push_back(temp);  
                                                }
                                            }
                                        }
                                    }
                                }
                                else {

                                    vector<implem> v4;

                                    poly_quad pq4 = poly_to_poly_quad(p2.second ,size, nb_elem);
                                    int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                    if (indice4 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p2.second);
                                        imp.quad_sol.insert(pq4);
                                        imp.formula = "";
                                        v4.push_back(imp);
                                    }
                                    else {
                                        indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                        if (indice4 == -1){
                                            continue;
                                        }
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice4].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice4].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                    }

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);
                                                    set_quad_4.erase(0);

                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG4_V2)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_4;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                        string f = "";

                                                        if (v2.size()==1){
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7) ^ p8";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7) ^ p8 ^ p9";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                }
                                                            }
                                                        }
                                                        else {
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10 ^ p11";
                                                                }
                                                            }
                                                        }
                                                        temp.formula = f;
                                                        res.push_back(temp); 
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (!test_ps){
                    for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                        for (uint32_t j3=0; j3<size; j3++){

                            create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                            pair<poly,poly> p3;
                            p3 = poly_div(reste, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                            uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                            uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                            if (m3 <= 2 && n3 <= 2){
                            // y = s[i] * p.first + p3.first * s[i3] + p3.second

                                vector<implem> v1;
                                                
                                implem imp;
                                imp.op_sol.push_back(P1);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);

                                vector<implem> v2;

                                int32_t indice2 = indiceP2;
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = indice_mP2;
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                }

                                vector<implem> v3;

                                poly P3;
                                P3 = truncate_lin(p3.first, size);
                                poly P3_lin;
                                POLY_ADD(P3, p3.first, P3_lin);
                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P3);
                                    imp.op_sol.push_back(P3_lin);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1){
                                        continue;
                                    }
                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice3].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice3].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P3_lin);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                }

                                vector<implem> v4;

                                poly P4;
                                P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                poly P4_lin;
                                POLY_ADD(P4, set_op_l_plus_lin[i3], P4_lin);
                                poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                implem imp4;
                                imp4.op_sol.push_back(P4);
                                imp4.op_sol.push_back(P4_lin);
                                imp4.quad_sol.insert(pq4);
                                imp4.formula = "";
                                v4.push_back(imp4);

                                if (m3 == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);
                                                    set_quad_4.erase(0);

                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG4_V2)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_4;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                        string f = "";
                                                        if (v2.size()==1){
                                                            if (v3.size()==1){
                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8)";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9)";
                                                            }
                                                        }
                                                        else {
                                                            if (v3.size()==1){
                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9)";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10)";
                                                            }
                                                        }
                                                        temp.formula = f;  
                                                        res.push_back(temp);
                                                    }                                                    
                                                }                                                   
                                            }
                                        }
                                    }
                                }
                                else {

                                    vector<implem> v5;

                                    poly P5;
                                    P5 = truncate_lin(p3.second, size);
                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);

                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                    if (indice5 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p3.second);
                                        imp.quad_sol.insert(pq5);
                                        imp.formula = "";
                                        v5.push_back(imp);
                                    }
                                    else {
                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                        if (indice5 == -1){
                                            continue;
                                        }
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice5].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice5].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.formula = "";
                                            v5.push_back(imp);
                                        }
                                    }

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;
                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);

                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                        set_quad_5.merge(set_5);
                                                        set_quad_5.erase(0);

                                                        set<poly_quad> s_quad_temp = set_quad_5;
                                                        uint32_t r = Rank(s_quad_temp);

                                                        if (r<=MAX_RANK_DEG4_V2)   {
                                                            implem temp; 
                                                            temp.quad_sol = set_quad_5;
                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());

                                                            string f = "";

                                                            if (v2.size()==1){
                                                                if (v3.size()==1){
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8) ^ p9";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8) ^ p9 ^ p10";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                            }
                                                            else {
                                                                if (v3.size()==1){
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10) ^ p11";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10) ^ p11 ^ p12";
                                                                    }
                                                                }
                                                            }
                                                            temp.formula = f;
                                                            res.push_back(temp);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                if (n3 != 0 && n3 <= 2){
                                    for (uint32_t jj3=0; jj3<size; jj3++){

                                        if (jj3 == j3){
                                            continue;
                                        }

                                        create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                        pair<poly,poly> pp3;
                                        pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                        uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                        uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);

                                        if (mm3 <= 2 && nn3 <= 2){
                                        // y = p.first * s[i] + (p3.first + pp3.first) * s[i3] + pp3.second
                                            vector<implem> v1;
                                
                                            implem imp;
                                            imp.op_sol.push_back(P1);
                                            imp.op_sol.push_back(P1_lin);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            int32_t indice2 = indiceP2;
                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = indice_mP2;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P2_lin);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                            }

                                            vector<implem> v3;

                                            poly PP3;
                                            POLY_ADD(p3.first, pp3.first, PP3);
                                            poly P3 = truncate_lin(PP3, size);
                                            poly P3_lin;
                                            POLY_ADD(PP3, P3, P3_lin);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P3);
                                                imp.op_sol.push_back(P3_lin);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P3_lin);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                            }

                                            vector<implem> v4;

                                            poly P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                            poly P4_lin;
                                            POLY_ADD(P4, set_op_l_plus_lin[i3], P4_lin);
                                            poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                            implem imp4;
                                            imp4.op_sol.push_back(P4);
                                            imp4.op_sol.push_back(P4_lin);
                                            imp4.quad_sol.insert(pq4);
                                            imp4.formula = "";
                                            v4.push_back(imp4);

                                            if (mm3 == 0){
                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                            
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);
                                                                set_quad_4.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_4;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG4_V2)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_4;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    string f = "";
                                                                    
                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9)";
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10)";
                                                                        }
                                                                    }
                                                                        
                                                                    temp.formula = f;  
                                                                    res.push_back(temp); 
                                                                }
                                                            }
                                                        }                                                    
                                                    }
                                                }
                                            }
                                            else {

                                                vector<implem> v5;

                                                poly P5;
                                                P5 = truncate_lin(pp3.second, size);
                                                poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);

                                                int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                if (indice5 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(pp3.second);
                                                    imp.quad_sol.insert(pq5);
                                                    imp.formula = "";
                                                    v5.push_back(imp);
                                                }
                                                else {
                                                    indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                    if (indice5 == -1){
                                                        continue;
                                                    }
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                }

                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                                
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;
                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);

                                                                for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                    set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                    set<poly_quad> set_quad_5 = set_quad_4;
                                                                    set_quad_5.merge(set_5);
                                                                    set_quad_5.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_5;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG4_V2)   {
                                                                        implem temp; 
                                                                        temp.quad_sol = set_quad_5;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        string f = "";
                                                                        
                                                                        if (v2.size()==1){
                                                                            if (v3.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8) ^ p9";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7 ^ p8) ^ p9 ^ p10";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9) ^ p10";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8 ^ p9) ^ p10 ^ p11";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v3.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9) ^ p10";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8 ^ p9) ^ p10 ^ p11";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10) ^ p11";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9 ^ p10) ^ p11 ^ p12";
                                                                                }
                                                                            }
                                                                        }
                                                                        
                                                                        temp.formula = f;
                                                                        res.push_back(temp); 
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (n!=0 && n<=2){
                    for (uint32_t jj=0; jj<size; jj++){

                        if (jj == j){
                            continue;
                        }
                        
                        create_monomial_order(monomial_order, nb_elem, size, jj);

                        pair<poly,poly> pp;
                        pp = poly_div(p.second, set_op_l_plus_lin[i], monomial_order, size, nb_elem);
                        uint32_t mm = pp.second.algebraic_degree(nb_elem);
                        uint32_t nn = pp.first.algebraic_degree(nb_elem);

                        if (mm <=3 && nn<=2)   {
                        // y = (p.first + pp.first) * s[i] + pp.second

                            poly P1;
                            P1 = truncate_lin(set_op_l_plus_lin[i], size);
                            poly P1_lin;
                            POLY_ADD(P1, set_op_l_plus_lin[i], P1_lin);
                            poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                            poly PP2;
                            POLY_ADD(p.first, pp.first, PP2);
                            poly P2 = truncate_lin(PP2, size);
                            poly P2_lin;
                            POLY_ADD(PP2, P2, P2_lin);

                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
                            int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                            int32_t indice_mP2 = -1;
                            if (indiceP2 == -1){
                                indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                                if (indice_mP2 == -1){
                                    continue;
                                }
                            }

                            poly reste;
                            POLY_ADD(pp.second, y3_2, reste);

                            bool test_pps = false;

                            for (uint32_t i2=0; i2<l.size(); i2++){
                                for (uint32_t j2=0; j2<size; j2++){
                                    if (mon_dom_in(l[i2],j2, size, nb_elem)){

                                        create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                                        pair<poly,poly> p2;
                                        p2 = poly_div(reste, l[i2], monomial_order, size, nb_elem);

                                        uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                        uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                                        if (m2 <= 2 && n2 <= 2){
                                        // y = (p.first + pp.first) * s[i] + p2.first * l[i2] + p2.second
                                        bool test_pps = true;

                                            vector<implem> v1;
                                                
                                            implem imp;
                                            imp.op_sol.push_back(P1);
                                            imp.op_sol.push_back(P1_lin);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            int32_t indice2 = indiceP2;
                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = indice_mP2;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P2_lin);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                            }

                                            vector<implem> v3;
                                            poly P3;
                                            P3 = truncate_lin(p2.first, size);
                                            poly P3_lin;
                                            POLY_ADD(P3, p2.first, P3_lin);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P3);
                                                imp.op_sol.push_back(P3_lin);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                if (indice3 == -1){
                                                    continue;
                                                }
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P3_lin);
                                                    imp.op_sol.push_back(l[i2]);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                            }

                                            if (m2 == 0){
                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);
                                                            set_quad_3.erase(0);

                                                            set<poly_quad> s_quad_temp = set_quad_3;
                                                            uint32_t r = Rank(s_quad_temp);

                                                            if (r<=MAX_RANK_DEG4_V2)   {
                                                                implem temp; 
                                                                temp.quad_sol = set_quad_3;
                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                string f = "";
                                                                if (v2.size()==1){
                                                                    if (v3.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7)";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8)";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v3.size()==1){
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8)";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9)";
                                                                    }
                                                                }
                                                                temp.formula = f;
                                                                res.push_back(temp);  
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else {

                                                vector<implem> v4;

                                                poly P4;
                                                P4 = truncate_lin(p2.second, size);
                                                poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                                int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                                if (indice4 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(p2.second);
                                                    imp.quad_sol.insert(pq4);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                                else {
                                                    indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                    if (indice4 == -1){
                                                        continue;
                                                    }
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                }

                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);
                                                                set_quad_4.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_4;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG4_V2)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_4;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    string f = "";

                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7) ^ p8";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6) * (p7) ^ p8 ^ p9";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * (p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                            }
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * (p3 ^ p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10 ^ p11";
                                                                            }
                                                                        }
                                                                    }
                                                                    temp.formula = f;
                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            if (!test_pps){
                                for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                                    for (uint32_t j3=0; j3<size; j3++){

                                        create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                                        pair<poly,poly> p3;
                                        p3 = poly_div(reste, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                        uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                                        uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                                        if (m3 <= 2 && n3 <= 2){
                                        // y = (p.first + pp.first) * s[i] + p3.first * s[i3] + p3.second
                                            vector<implem> v1;
                                                
                                            implem imp;
                                            imp.op_sol.push_back(set_op_l_plus_lin[i]);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            int32_t indice2 = indiceP2;
                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(PP2);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = indice_mP2;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    poly D; 
                                                    POLY_ADD(PP2, P2, D); //Pour récupérer le delta linéaire
                                                    imp.op_sol.push_back(D);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                            }

                                            vector<implem> v3;

                                            poly P3;
                                            P3 = truncate_lin(p3.first, size);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p3.first);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                if (indice3 == -1){
                                                    continue;
                                                }
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    poly D; 
                                                    POLY_ADD(P3, p3.first, D); //Pour récupérer le delta linéaire
                                                    imp.op_sol.push_back(D);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                            }

                                            vector<implem> v4;

                                            poly P4;
                                            P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                            poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                            implem imp4;
                                            imp4.op_sol.push_back(set_op_l_plus_lin[i3]);
                                            imp4.quad_sol.insert(pq4);
                                            imp4.formula = "";
                                            v4.push_back(imp4);

                                            if (m3 == 0){
                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);
                                                                set_quad_4.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_4;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG4_V2)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_4;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    string f = "";
                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            f = "(p1 * p2) ^ (p3 * p4)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 * p2) ^ (p3 ^ p4 ^ p5) * (p6)";
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6)";
                                                                        }
                                                                        else {
                                                                            f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * p8";
                                                                        }
                                                                    }
                                                                    
                                                                    temp.formula = f;  
                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else {

                                                vector<implem> v5;

                                                poly P5;
                                                P5 = truncate_lin(p3.second, size);
                                                poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);

                                                int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                if (indice5 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(p3.second);
                                                    imp.quad_sol.insert(pq5);
                                                    imp.formula = "";
                                                    v5.push_back(imp);
                                                }
                                                else {
                                                    indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                    if (indice5 == -1){
                                                        continue;
                                                    }
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                }

                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);

                                                                for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                    set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                    set<poly_quad> set_quad_5 = set_quad_4;
                                                                    set_quad_5.merge(set_5);
                                                                    set_quad_5.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_5;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG4_V2)   {
                                                                        implem temp; 
                                                                        temp.quad_sol = set_quad_5;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        string f = "";
                                                                        if (v2.size()==1){
                                                                            if (v3.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 * p2) ^ (p3 * p4) ^ p5";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 * p2) ^ (p3 * p4) ^ p5 ^ p6";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 * p2) ^ (p2 ^ p3 ^ p4) * (p5) ^ p6";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 * p2) ^ (p2 ^ p3 ^ p4) * (p5) ^ p6 ^ p7";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v3.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6) ^ p7";
                                                                                }
                                                                                else {
                                                                                    f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6) ^ p7 ^ p8";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                                }
                                                                                else {
                                                                                    f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                                }
                                                                            }
                                                                        }
                                                                        
                                                                        temp.formula = f;
                                                                        res.push_back(temp); 
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {
                                            if (n3 != 0 && n3 <= 2){
                                                for (uint32_t jj3=0; jj3<size; jj3++){

                                                    if (jj3 == j3){
                                                        continue;
                                                    }

                                                    create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                                    pair<poly,poly> pp3;
                                                    pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                                    uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                                    uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);

                                                    if (mm3 <= 2 && nn3 <= 2){
                                                    // y = (p.first + pp.first) * s[i] + (p3.first + pp3.first) * s[i3] + pp3.second
                                                        vector<implem> v1;
                                                
                                                        implem imp;
                                                        imp.op_sol.push_back(set_op_l_plus_lin[i]);
                                                        imp.quad_sol.insert(pq1);
                                                        imp.formula = "";
                                                        v1.push_back(imp);

                                                        vector<implem> v2;

                                                        int32_t indice2 = indiceP2;
                                                        if (indice2 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(PP2);
                                                            imp.quad_sol.insert(pq2);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                        else {
                                                            indice2 = indice_mP2;
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                poly D; 
                                                                POLY_ADD(PP2, P2, D); //Pour récupérer le delta linéaire
                                                                imp.op_sol.push_back(D);
                                                                imp.formula = "";
                                                                v2.push_back(imp);
                                                            }
                                                        }

                                                        vector<implem> v3;
                                                        poly PP3;
                                                        POLY_ADD(p3.first, pp3.first, PP3)
                                                        poly P3 = truncate_lin(PP3, size);
                                                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                                        int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                                        if (indice3 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(PP3);
                                                            imp.quad_sol.insert(pq3);
                                                            imp.formula = "";
                                                            v3.push_back(imp);
                                                        }
                                                        else {
                                                            indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                            if (indice3 == -1){
                                                                continue;
                                                            }
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                poly D; 
                                                                POLY_ADD(PP3, P3, D); //Pour récupérer le delta linéaire
                                                                imp.op_sol.push_back(D);
                                                                imp.formula = "";
                                                                v3.push_back(imp);
                                                            }
                                                        }

                                                        vector<implem> v4;

                                                        poly P4;
                                                        P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                                        poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);

                                                        implem imp4;
                                                        imp4.op_sol.push_back(set_op_l_plus_lin[i3]);
                                                        imp4.quad_sol.insert(pq4);
                                                        imp4.formula = "";
                                                        v4.push_back(imp4);

                                                        if (mm3 == 0){
                                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                            
                                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                                    set_quad_2.merge(set_2);

                                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                                        set_quad_3.merge(set_3);

                                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                                            set_quad_4.merge(set_4);
                                                                            set_quad_4.erase(0);

                                                                            set<poly_quad> s_quad_temp = set_quad_4;
                                                                            uint32_t r = Rank(s_quad_temp);

                                                                            if (r<=MAX_RANK_DEG4_V2)   {
                                                                                implem temp; 
                                                                                temp.quad_sol = set_quad_4;
                                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                string f = "";
                                                                                if (v2.size()==1){
                                                                                    if (v3.size()==1){
                                                                                        f = "(p1 * p2) ^ (p3 * p4)";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 * p2) ^ (p3 ^ p4 ^ p5) * (p6)";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v3.size()==1){
                                                                                        f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6)";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * p8";
                                                                                    }
                                                                                }
                                                                                
                                                                                temp.formula = f;  
                                                                                res.push_back(temp);  
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        else {

                                                            vector<implem> v5;

                                                            poly P5;
                                                            P5 = truncate_lin(pp3.second, size);
                                                            poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);

                                                            int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                            if (indice5 != -1)   {
                                                                implem imp;
                                                                imp.op_sol.push_back(pp3.second);
                                                                imp.quad_sol.insert(pq5);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                            else {
                                                                indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                                if (indice5 == -1){
                                                                    continue;
                                                                }
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.formula = "";
                                                                    v5.push_back(imp);
                                                                }
                                                            }

                                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                                
                                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                                    set_quad_2.merge(set_2);

                                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                                        set_quad_3.merge(set_3);

                                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                                            set_quad_4.merge(set_4);

                                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                                set_quad_5.merge(set_5);
                                                                                set_quad_5.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_5;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG4_V2)   {
                                                                                    implem temp; 
                                                                                    temp.quad_sol = set_quad_5;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    string f = "";
                                                                                    if (v2.size()==1){
                                                                                        if (v3.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 * p2) ^ (p3 * p4) ^ p5";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 * p2) ^ (p3 * p4) ^ p5 ^ p6";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 * p2) ^ (p2 ^ p3 ^ p4) * (p5) ^ p6";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 * p2) ^ (p2 ^ p3 ^ p4) * (p5) ^ p6 ^ p7";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v3.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6) ^ p7";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 * p6) ^ p7 ^ p8";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1) * (p2 ^ p3 ^ p4) ^ (p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    
                                                                                    temp.formula = f;
                                                                                    res.push_back(temp);
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    if (res.size() == 0){

        vector<poly> map_xor_l_plus_lin = create_map_xor_l_plus_lin(size);

        for (uint32_t i=0; i<map_xor_l_plus_lin.size(); i++) {
            for (uint32_t j=0; j<size; j++){

                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p;
                p = poly_div(y4, map_xor_l_plus_lin[i], monomial_order, size, nb_elem);
                uint32_t m = p.second.algebraic_degree(nb_elem);
                uint32_t n = p.first.algebraic_degree(nb_elem);

                if (m<=3 && n<=2)   {
                    //y = p.first * m[i] + p.second + y3_2

                    poly P1 = truncate_lin(map_xor_l_plus_lin[i], size);
                    poly P1_lin;
                    POLY_ADD(P1, map_xor_l_plus_lin[i], P1_lin);
                    poly_quad pq1 = poly_to_poly_quad(P1, size, nb_elem);  
                    int32_t indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                    
                    poly P2 = truncate_lin(p.first, size);
                    poly P2_lin;
                    POLY_ADD(P2, p.first, P2_lin);
                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                    int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                    int32_t indice_mP2 = -1;
                    if (indiceP2 == -1){
                        indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                        if (indice_mP2 == -1){
                            continue;
                        }
                    }

                    bool test_ps = false;

                    poly reste;
                    POLY_ADD(p.second, y3_2, reste);

                    for (uint32_t i2=0; i2<l.size(); i2++){
                        for (uint32_t j2=0; j2<size; j2++){
                            if (mon_dom_in(l[i2],j2, size, nb_elem)){
    
                                create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                                pair<poly,poly> p2;
                                p2 = poly_div(reste, l[i2], monomial_order, size, nb_elem);
    
                                uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                uint32_t n2 = p2.first.algebraic_degree(nb_elem);
    
                                if (m2 <= 2 && n2 <= 2){
                                //y = p.first * m[i] + p2.first * l[i2] + p2.second
                                bool test_ps = true;
    
                                    vector<implem> v1;

                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.formula = "";
                                        v1.push_back(imp);
                                    }  
    
                                    vector<implem> v2;
    
                                    int32_t indice2 = indiceP2;
                                    if (indice2 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(P2);
                                        imp.op_sol.push_back(P2_lin);
                                        imp.quad_sol.insert(pq2);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                    else {
                                        indice2 = indice_mP2;
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
    
                                    vector<implem> v3;
                                    poly P3;
                                    P3 = truncate_lin(p2.first, size);
                                    poly P3_lin;
                                    POLY_ADD(P3, p2.first, P3_lin);
                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);
    
                                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                    if (indice3 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(P3);
                                        imp.op_sol.push_back(P3_lin);
                                        imp.op_sol.push_back(l[i2]);
                                        imp.quad_sol.insert(pq3);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                    else {
                                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                        if (indice3 == -1){
                                            continue;
                                        }
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P3_lin);
                                            imp.op_sol.push_back(l[i2]);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
    
                                    if (m2 == 0){
                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                set_quad_2.merge(set_2);
    
                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                    set_quad_3.merge(set_3);
                                                    set_quad_3.erase(0);
    
                                                    set<poly_quad> s_quad_temp = set_quad_3;
                                                    uint32_t r = Rank(s_quad_temp);
    
                                                    if (r<=MAX_RANK_DEG4_V2)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_3;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        
                                                        string f = "";
                                                        if (v2.size()==1){
                                                            if (v3.size()==1){
                                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7) * (p8)";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9)";
                                                            }
                                                        }
                                                        else {
                                                            if (v3.size()==1){
                                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8) * (p9)";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8 ^ p9) * (p10)";
                                                            }
                                                        }
                                                        
                                                        temp.formula = f;

                                                        res.push_back(temp);  
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else {
    
                                        vector<implem> v4;
    
                                        
                                        poly_quad pq4 = poly_to_poly_quad(p2.second ,size, nb_elem);
                                        int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                        if (indice4 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(p2.second);
                                            imp.quad_sol.insert(pq4);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                        else {
                                            indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                            if (indice4 == -1){
                                                continue;
                                            }
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                        }
    
                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                            
                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                set_quad_2.merge(set_2);
    
                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                    set_quad_3.merge(set_3);
    
                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                        set_quad_4.merge(set_4);
                                                        set_quad_4.erase(0);
    
                                                        set<poly_quad> s_quad_temp = set_quad_4;
                                                        uint32_t r = Rank(s_quad_temp);
    
                                                        if (r<=MAX_RANK_DEG4_V2)   {
                                                            implem temp; 
                                                            temp.quad_sol = set_quad_4;
                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                            
                                                            string f = "";

                                                            if (v2.size()==1){
                                                                if (v3.size()==1){
                                                                    if (v4.size()==1){
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v4.size()==1){
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5) ^ (p6 ^ p7 ^ p8) * (p9) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                            }
                                                            else {
                                                                if (v3.size()==1){
                                                                    if (v4.size()==1){
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8) * (p9) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8) * (p9) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v4.size()==1){
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8 ^ p9) * (p10) ^ p11";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2 ^ p3) * (p4 ^ p5 ^ p6) ^ (p7 ^ p8 ^ p9) * (p10) ^ p11 ^ p12";
                                                                    }
                                                                }
                                                            }
                                                            temp.formula = f;
                                                            res.push_back(temp);  
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
    
                    if (!test_ps){
                        for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                            for (uint32_t j3=0; j3<size; j3++){
    
                                create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                                pair<poly,poly> p3;
                                p3 = poly_div(reste, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);
    
                                uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                                uint32_t n3 = p3.first.algebraic_degree(nb_elem);
    
                                if (m3 <= 2 && n3 <= 2){
                                // y = p.first * s[i] + p3.first * s[i3] + p3.second
    
                                    vector<implem> v1;

                                    for (uint32_t m=0; m<10; m++)   {
                                        implem imp;
                                        poly_quad op_q1 = map_xor[indice1].second[m][0];
                                        poly op1(op_q1,size);
                                        imp.op_sol.push_back(op1);
                                        imp.quad_sol.insert(op_q1);
                                        poly_quad op_q2 = map_xor[indice1].second[m][1];
                                        poly op2(op_q2,size);
                                        imp.op_sol.push_back(op2);
                                        imp.quad_sol.insert(op_q2);
                                        poly D; 
                                        POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                                        imp.op_sol.push_back(D);
                                        imp.formula = "";
                                        v1.push_back(imp);
                                    }  

                                    vector<implem> v2;
    
                                    int32_t indice2 = indiceP2;
                                    if (indice2 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p.first);
                                        imp.quad_sol.insert(pq2);
                                        imp.formula = "";
                                        v2.push_back(imp);
                                    }
                                    else {
                                        indice2 = indice_mP2;
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            poly D; 
                                            POLY_ADD(p.first, P2, D); //Pour récupérer le delta linéaire
                                            imp.op_sol.push_back(D);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
    
                                    vector<implem> v3;
    
                                    poly P3;
                                    P3 = truncate_lin(p3.first, size);
                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);
    
                                    int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                    if (indice3 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p3.first);
                                        imp.quad_sol.insert(pq3);
                                        imp.formula = "";
                                        v3.push_back(imp);
                                    }
                                    else {
                                        indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                        if (indice3 == -1){
                                            continue;
                                        }
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            poly D; 
                                            POLY_ADD(P3, p3.first, D); //Pour récupérer le delta linéaire
                                            imp.op_sol.push_back(D);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
    
                                    vector<implem> v4;
    
                                    poly P4;
                                    P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                    poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);
    
                                    implem imp4;
                                    imp4.op_sol.push_back(set_op_l_plus_lin[i3]);
                                    imp4.quad_sol.insert(pq4);
                                    imp4.formula = "";
                                    v4.push_back(imp4);
    
                                    if (m3 == 0){
                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);
    
                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);
    
                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);
                                                            set_quad_4.erase(0);
    
                                                            set<poly_quad> s_quad_temp = set_quad_4;
                                                            uint32_t r = Rank(s_quad_temp);
    
                                                            if (r<=MAX_RANK_DEG4_V2)   {
                                                                implem temp; 
                                                                temp.quad_sol = set_quad_4;
                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                temp.formula = "4215";
                                                                res.push_back(temp);  
                                                            }
                                                        }
                                                    }                                                    
                                                }
                                        }
                                    }
                                    else {
    
                                        vector<implem> v5;
    
                                        poly P5;
                                        P5 = truncate_lin(p3.second, size);
                                        poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
    
                                        int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                        if (indice5 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(p3.second);
                                            imp.quad_sol.insert(pq5);
                                            imp.formula = "";
                                            v5.push_back(imp);
                                        }
                                        else {
                                            indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                            if (indice5 == -1){
                                                continue;
                                            }
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v5.push_back(imp);
                                            }
                                        }
    
                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                        
                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;
                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                set_quad_2.merge(set_2);
    
                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                    set_quad_3.merge(set_3);
    
                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                        set_quad_4.merge(set_4);
    
                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                            set_quad_5.merge(set_5);
                                                            set_quad_5.erase(0);
    
                                                            set<poly_quad> s_quad_temp = set_quad_5;
                                                            uint32_t r = Rank(s_quad_temp);
    
                                                            if (r<=MAX_RANK_DEG4_V2)   {
                                                                implem temp; 
                                                                temp.quad_sol = set_quad_5;
                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                temp.formula = "4216";
                                                                res.push_back(temp);  
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else {
                                    if (n3 != 0 && n3 <= 2){
                                        for (uint32_t jj3=0; jj3<size; jj3++){
    
                                            if (jj3 == j3){
                                                continue;
                                            }
    
                                            create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                            pair<poly,poly> pp3;
                                            pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);
    
                                            uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                            uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);
    
                                            if (mm3 <= 2 && nn3 <= 2){
                                            // y = p.first * m[i] + (p3.first + pp3.first) * s[i3] + pp3.second
                                                vector<implem> v1;

                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice1].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice1].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    poly D; 
                                                    POLY_ADD(map_xor_l_plus_lin[i], P1, D); //Pour récupérer le delta linéaire
                                                    imp.op_sol.push_back(D);
                                                    imp.formula = "";
                                                    v1.push_back(imp);
                                                } 

                                                vector<implem> v2;
    
                                                int32_t indice2 = indiceP2;
                                                if (indice2 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(p.first);
                                                    imp.quad_sol.insert(pq2);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                                else {
                                                    indice2 = indice_mP2;
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        poly D; 
                                                        POLY_ADD(p.first, P2, D); //Pour récupérer le delta linéaire
                                                        imp.op_sol.push_back(D);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                }
    
                                                vector<implem> v3;
    
                                                poly PP3;
                                                POLY_ADD(p3.first, pp3.first, PP3);
                                                poly P3 = truncate_lin(PP3, size);
                                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);
    
                                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                                if (indice3 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(PP3);
                                                    imp.quad_sol.insert(pq3);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                                else {
                                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        poly D; 
                                                        POLY_ADD(P3, PP3, D); //Pour récupérer le delta linéaire
                                                        imp.op_sol.push_back(D);
                                                        imp.formula = "";
                                                        v3.push_back(imp);
                                                    }
                                                }
    
                                                vector<implem> v4;
    
                                                poly P4 = truncate_lin(set_op_l_plus_lin[i3], size);
                                                poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);
    
                                                implem imp4;
                                                imp4.op_sol.push_back(set_op_l_plus_lin[i3]);
                                                imp4.quad_sol.insert(pq4);
                                                imp4.formula = "";
                                                v4.push_back(imp4);
    
                                                if (mm3 == 0){
                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                                
                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;
    
                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);
    
                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);
    
                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);
                                                                        set_quad_4.erase(0);
    
                                                                        set<poly_quad> s_quad_temp = set_quad_4;
                                                                        uint32_t r = Rank(s_quad_temp);
    
                                                                        if (r<=MAX_RANK_DEG4_V2)   {
                                                                            implem temp; 
                                                                            temp.quad_sol = set_quad_4;
                                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                            temp.formula = "4217";
                                                                            res.push_back(temp);  
                                                                        }
                                                                    }
                                                                }                                                    
                                                            }
                                                    }
                                                }
                                                else {
    
                                                    vector<implem> v5;
    
                                                    poly P5;
                                                    P5 = truncate_lin(pp3.second, size);
                                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
    
                                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                    if (indice5 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(pp3.second);
                                                        imp.quad_sol.insert(pq5);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                    else {
                                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                        if (indice5 == -1){
                                                            continue;
                                                        }
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.formula = "";
                                                            v5.push_back(imp);
                                                        }
                                                    }
    
                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                                    
                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;
                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                                            set<poly_quad> set_quad_2 = set_quad_1;
                                                            set_quad_2.merge(set_2);
    
                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                                set_quad_3.merge(set_3);
    
                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                                    set_quad_4.merge(set_4);
    
                                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                                        set_quad_5.merge(set_5);
                                                                        set_quad_5.erase(0);
    
                                                                        set<poly_quad> s_quad_temp = set_quad_5;
                                                                        uint32_t r = Rank(s_quad_temp);
    
                                                                        if (r<=MAX_RANK_DEG4_V2)   {
                                                                            implem temp; 
                                                                            temp.quad_sol = set_quad_5;
                                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                            temp.formula = "4218";
                                                                            res.push_back(temp);  
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    res = remove_duplicates(res);
    cout<<" Thread : "<<omp_get_thread_num()<<" res_size = "<<res.size()<<endl;
   
    return res;
}

//DEG 5

vector<implem> add_to_op_selec_deg_5_v1(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)   {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    if (verbose)  {
        #pragma omp critical
        {
        cout<<"Début nv add_to_op_selec v1 avec précalcul pour thread "<<omp_get_thread_num()<<" et poly :"<<endl;
        y.print_poly(size);
        cout<<endl;
        }
    }

    poly y5 = truncate_except_deg(y, 5, size, nb_elem);
    poly y4 = truncate_except_deg(y, 4, size, nb_elem);
    poly y3 = truncate_except_deg(y, 3, size, nb_elem);
    poly y2 = truncate_except_deg(y, 2, size, nb_elem);

    poly y5_4;
    POLY_ADD(y5,y4,y5_4);

    poly y3_2;
    POLY_ADD(y3,y2,y3_2);

    uint32_t ent_y5_4 = poly_deg5_4_to_uint(y5_4);

    vector<implem_deg5_4> vect_sol_deg5_4 = parse_file_and_create_sol(ent_y5_4);

    for (uint32_t i=0; i<vect_sol_deg5_4.size(); i++){

        poly quotient = vect_sol_deg5_4[i].q;

        poly reste;
        POLY_ADD(vect_sol_deg5_4[i].r, y3_2, reste);
        poly test_reste = truncate_except_deg(reste, 3, size, nb_elem);

        for (uint32_t j=0; j<size; j++){
            create_monomial_order(monomial_order, nb_elem, size, j);                                                   
            pair<poly,poly> p;
            p = poly_div(reste, vect_sol_deg5_4[i].quad, monomial_order, size, nb_elem);

            uint32_t m = p.second.algebraic_degree(nb_elem);
            uint32_t n = p.first.algebraic_degree(nb_elem);

            if (m<=2 && n<=3){
            //y = (quotient + p.first) * quad + p.second

                poly P1;
                P1 = truncate_lin(vect_sol_deg5_4[i].quad, size);
                poly P1_lin;
                POLY_ADD(P1,vect_sol_deg5_4[i].quad, P1_lin );
                poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                poly q;
                POLY_ADD(quotient, p.first, q);
                
                bool test_quotient_div_l = false;

                for (uint32_t i2=0; i2<l.size(); i2++){
                    for (uint32_t j2=0; j2<size; j2++){

                        if (mon_dom_in(l[i2],j2, size, nb_elem)){

                            create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                            pair<poly,poly> p2;
                            p2 = poly_div(q, l[i2], monomial_order, size, nb_elem);

                            uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                            uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                            if (m2 <= 2 && n2 <= 2){
                            // y = q * quad  + p.second
                            // y = (p2.first * l[i2] + p2.second) * quad  + p.second
                            test_quotient_div_l = true;

                                vector<implem> v1;
                                                
                                implem imp;
                                imp.op_sol.push_back(P1);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);

                                vector<implem> v2;

                                poly P2;
                                P2 = truncate_lin(p2.first, size);
                                poly P2_lin;
                                POLY_ADD(P2, p2.first, P2_lin );
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.op_sol.push_back(l[i2]);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                    if (indice2 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.op_sol.push_back(l[i2]);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
                                }

                                vector<implem> v3;

                                poly P3;
                                poly m1;
                                POLY_MUL(p2.first, l[i2], m1, nb_elem);
                                POLY_ADD(m1, q, P3);
                                poly P3_lin;
                                POLY_ADD(P3, p2.second, P3_lin);
                                poly_quad pq3 = poly_to_poly_quad(p2.second ,size, nb_elem);

                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(p2.second);
                                    imp.op_sol.push_back(P3_lin);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P3_lin);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
                                }

                                if (m == 0){
                                // If p.second is linear (ie == 0)
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);

                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG5_V1)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v2.size()==1){
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 )";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 )";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 )";
                                                        }
                                                        else {
                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 )";
                                                        }
                                                    }
                                                    temp.formula = f;
                                                    res.push_back(temp);  
                                                }
                                            }
                                        }
                                    }
                                }

                                else {
                                    vector<implem> v4;

                                    poly_quad pq4 = poly_to_poly_quad(p.second ,size, nb_elem);
                                    int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                    if (indice4 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p.second);
                                        imp.quad_sol.insert(pq4);
                                        imp.formula = "";
                                        v4.push_back(imp);
                                    }
                                    else {
                                        indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                        if (indice4 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                        }
                                    }

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);
                                                    set_quad_4.erase(0);

                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG5_V1)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_4;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                        string f = "";

                                                        if (v2.size()==1){
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ p8";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ p8 ^ p9";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ p9";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ p9 ^ p10";
                                                                }
                                                            }
                                                        }
                                                        else {
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ) ^ p9";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ) ^ p9 ^ p10";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                }
                                                                else {
                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                }
                                                            }
                                                        }
                                                        
                                                        temp.formula = f;         
                                                        res.push_back(temp);  
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (!test_quotient_div_l){
                    for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                        for (uint32_t j3=0; j3<size; j3++){

                            create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                            pair<poly,poly> p3;
                            p3 = poly_div(q, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                            uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                            uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                            if (m3 <= 2 && n3 <= 2){
                            // y = (p3.first * s[i3] + p3.second) * quad  + p.second

                                vector<implem> v1;
                                             
                                implem imp;
                               
                                imp.op_sol.push_back(P1);
                                imp.op_sol.push_back(P1_lin);
                                imp.quad_sol.insert(pq1);
                                imp.formula = "";
                                v1.push_back(imp);

                                vector<implem> v2;

                                poly P2;
                                P2 = truncate_lin(p3.first, size);
                                poly P2_lin;
                                POLY_ADD(P2, p3.first, P2_lin);
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                if (indice2 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P2);
                                    imp.op_sol.push_back(P2_lin);
                                    imp.quad_sol.insert(pq2);
                                    imp.formula = "";
                                    v2.push_back(imp);
                                }
                                else {
                                    indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                    if (indice2 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                    }
                                }

                                vector<implem> v3;

                                poly P3;
                                P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                                poly P3_lin;
                                POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin)
                                poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                implem imp3;
                                imp3.op_sol.push_back(P3);
                                imp3.op_sol.push_back(P3_lin);
                                imp3.quad_sol.insert(pq3);
                                imp3.formula = "";
                                v3.push_back(imp3);

                                vector<implem> v4;

                                poly P4;
                                poly m1;
                                POLY_MUL(p3.first, set_op_l_plus_lin[i3], m1, nb_elem);
                                POLY_ADD(m1, q, P4);
                                poly P4_lin;
                                POLY_ADD(P4, p3.second, P4_lin);
                                poly_quad pq4 = poly_to_poly_quad(p3.second ,size, nb_elem);

                                int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                if (indice4 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(p3.second);
                                    imp.op_sol.push_back(P4_lin);
                                    imp.quad_sol.insert(pq4);
                                    imp.formula = "";
                                    v4.push_back(imp);
                                }
                                else {
                                    indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                    if (indice4 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice4].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice4].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P4_lin);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                    }
                                }

                                if (m == 0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);
                                                    set_quad_4.erase(0);

                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG5_V1)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_4;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                        string f = "";
                                                        if (v2.size()==1){
                                                            if (v4.size()==1){
                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 )";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 )";
                                                            }
                                                        }
                                                        else {
                                                            if (v4.size()==1){
                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 )";
                                                            }
                                                            else {
                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 )";
                                                            }
                                                        }
                                                        temp.formula = f;
                                                        res.push_back(temp);  
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else {
                                    vector<implem> v5;

                                    poly_quad pq5 = poly_to_poly_quad(p.second ,size, nb_elem);
                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                    if (indice5 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p.second);
                                        imp.quad_sol.insert(pq5);
                                        imp.formula = "";
                                        v5.push_back(imp);
                                    }
                                    else {
                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                        if (indice5 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v5.push_back(imp);
                                            }
                                        }
                                    }

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);

                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                        set_quad_5.merge(set_5);
                                                        set_quad_5.erase(0);

                                                        set<poly_quad> s_quad_temp = set_quad_5;
                                                        uint32_t r = Rank(s_quad_temp);

                                                        if (r<=MAX_RANK_DEG5_V1)   {
                                                            implem temp; 
                                                            temp.quad_sol = set_quad_5;
                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                            temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                            string f = "";

                                                            if (v2.size()==1){
                                                                if (v4.size()==1){
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9 ^ p10";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v5.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11 ^ p12";
                                                                    }
                                                                }
                                                            }
                                                            
                                                            temp.formula = f;
                                                            res.push_back(temp);  
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                if (n3 != 0 && n3 <= 2){
                                    for (uint32_t jj3=0; jj3<size; jj3++){

                                        if (j3 == jj3){
                                            continue;
                                        }

                                        create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                        pair<poly,poly> pp3;
                                        pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                        uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                        uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);

                                        if (mm3 <= 2 && nn3 <= 2){
                                        // y = ((p3.first + pp3.first) * s[i3] + pp3.second) * quad  + pp.second
                                            vector<implem> v1;
                                                
                                            implem imp;
                                            imp.op_sol.push_back(vect_sol_deg5_4[i].quad);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            poly PP2;
                                            POLY_ADD(p3.first, pp3.first, PP2);
                                            poly P2 = truncate_lin(PP2, size);
                                            poly P2_lin;
                                            POLY_ADD(P2, PP2, P2_lin);
                                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                if (indice2 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                }
                                            }
                                                
                                            vector<implem> v3;

                                            poly P3;
                                            P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                                            poly P3_lin;
                                            POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            implem imp3;
                                            imp3.op_sol.push_back(P3);
                                            imp3.op_sol.push_back(P3_lin);
                                            imp3.quad_sol.insert(pq3);
                                            imp3.formula = "";
                                            v3.push_back(imp3);

                                            vector<implem> v4;

                                            if (mm3 != 2){
                                                cout<<"pp3.second est linéaire"<<endl;
                                            }
                                            else {

                                                poly P4;
                                                poly m1;
                                                POLY_MUL(PP2, set_op_l_plus_lin[i3], m1, nb_elem);
                                                POLY_ADD(m1, q, P4);
                                                poly P4_lin;
                                                POLY_ADD(P4, pp3.second ,P4_lin)
                                                poly_quad pq4 = poly_to_poly_quad(pp3.second ,size, nb_elem);

                                                int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                                if (indice4 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(pp3.second);
                                                    imp.op_sol.push_back(P4_lin);
                                                    imp.quad_sol.insert(pq4);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                                else {
                                                    indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                    if (indice4 == -1)  {
                                                        continue;
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P4_lin);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                    }
                                                }

                                                if (m == 0){
                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                
                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                                            set<poly_quad> set_quad_2 = set_quad_1;
                                                            set_quad_2.merge(set_2);

                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                                set_quad_3.merge(set_3);

                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                                    set_quad_4.merge(set_4);
                                                                    set_quad_4.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG5_V1)   {
                                                                        implem temp; 
                                                                        temp.quad_sol = set_quad_4;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                
                                                                        string f = "";
                                                                        if (v2.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 )";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 )";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 )";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 )";
                                                                            }
                                                                        }
                                                                        temp.formula = f;
                                                                        res.push_back(temp); 
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                else {
                                                    vector<implem> v5;

                                                    poly_quad pq5 = poly_to_poly_quad(p.second ,size, nb_elem);
                                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                    if (indice5 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(p.second);
                                                        imp.quad_sol.insert(pq5);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                    else {
                                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                        if (indice5 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                                            set<poly_quad> set_quad_2 = set_quad_1;
                                                            set_quad_2.merge(set_2);

                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                                set_quad_3.merge(set_3);

                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                                    set_quad_4.merge(set_4);

                                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                                        set_quad_5.merge(set_5);
                                                                        set_quad_5.erase(0);

                                                                        set<poly_quad> s_quad_temp = set_quad_5;
                                                                        uint32_t r = Rank(s_quad_temp);

                                                                        if (r<=MAX_RANK_DEG5_V1)   {
                                                                            implem temp; 
                                                                            temp.quad_sol = set_quad_5;
                                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());

                                                                            string f = "";
                                                                            if (v2.size()==1){
                                                                                if (v4.size()==1){
                                                                                    if (v5.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9 ^ p10";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v5.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v4.size()==1){
                                                                                    if (v5.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v5.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11 ^ p12";
                                                                                    }
                                                                                }
                                                                            }
                                                            
                                                                            temp.formula = f;
                                                                            res.push_back(temp);
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                if (n != 0 && n<=3){
                    for (uint32_t jj=0; jj<size; jj++){
                        if (jj == j){
                            continue;
                        }

                        create_monomial_order(monomial_order, nb_elem, size, jj);                                                   
                        pair<poly,poly> pp;
                        pp = poly_div(p.second, vect_sol_deg5_4[i].quad, monomial_order, size, nb_elem);

                        uint32_t mm = pp.second.algebraic_degree(nb_elem);
                        uint32_t nn = pp.first.algebraic_degree(nb_elem);

                        if (mm <=2 && nn <= 3){
                            // y = quad * (quotient + p.first + pp.first) + pp.second

                            poly P1;
                            P1 = truncate_lin(vect_sol_deg5_4[i].quad, size);
                            poly P1_lin;
                            POLY_ADD(P1, vect_sol_deg5_4[i].quad, P1_lin);
                            poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                            poly q;
                            POLY_ADD(quotient, p.first, q);
                            POLY_ADD(q, pp.first, q);

                            bool test_quotient_div_l = false;

                            for (uint32_t i2=0; i2<l.size(); i2++){
                                for (uint32_t j2=0; j2<size; j2++){

                                    if (mon_dom_in(l[i2],j2, size, nb_elem)){

                                        create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                                        pair<poly,poly> p2;
                                        p2 = poly_div(q, l[i2], monomial_order, size, nb_elem);

                                        uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                                        uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                                        if (m2 <= 2 && n2 <= 2){
                                        // y = q * quad  + pp.second
                                        // y = (p2.first * l[i2] + p2.second) * quad  + pp.second
                                        //cout<<"On a une solution du type : y = (p2.first * l[i2] + p2.second) * quad  + pp.second"<<endl;
                                        test_quotient_div_l = true;

                                            vector<implem> v1;
                                                
                                            implem imp;
                                            imp.op_sol.push_back(P1);
                                            imp.op_sol.push_back(P1_lin);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            poly P2;
                                            P2 = truncate_lin(p2.first, size);
                                            poly P2_lin;
                                            POLY_ADD(P2, p2.first, P2_lin);
                                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                            if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                if (indice2 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.op_sol.push_back(l[i2]);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                }
                                            }

                                            vector<implem> v3;

                                            poly P3;
                                            poly m1;
                                            POLY_MUL(p2.first, l[i2], m1, nb_elem);
                                            POLY_ADD(m1, q, P3);
                                            poly P3_lin;
                                            POLY_ADD(p2.second, P3, P3_lin);
                                            poly_quad pq3 = poly_to_poly_quad(p2.second ,size, nb_elem);

                                            int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p2.second);
                                                imp.op_sol.push_back(P3_lin);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                                if (indice3 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.op_sol.push_back(P3_lin);
                                                        imp.formula = "";
                                                        v3.push_back(imp);
                                                    }
                                                }
                                            }

                                            if (mm == 0){
                                            // If pp.second is linear (ie == 0)
                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);
                                                            set_quad_3.erase(0);

                                                            set<poly_quad> s_quad_temp = set_quad_3;
                                                            uint32_t r = Rank(s_quad_temp);

                                                            if (r<=MAX_RANK_DEG5_V1)   {
                                                                implem temp; 
                                                                temp.quad_sol = set_quad_3;
                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                
                                                                string f = "";
                                                                if (v2.size()==1){
                                                                    if (v3.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 )";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 )";
                                                                    }
                                                                }
                                                                else {
                                                                    if (v3.size()==1){
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 )";
                                                                    }
                                                                    else {
                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 )";
                                                                    }
                                                                }
                                                                temp.formula = f;

                                                                res.push_back(temp);  
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            else {
                                                vector<implem> v4;

                                                poly_quad pq4 = poly_to_poly_quad(pp.second ,size, nb_elem);
                                                int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                                if (indice4 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(pp.second);
                                                    imp.quad_sol.insert(pq4);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                                else {
                                                    indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                    if (indice4 == -1)  {
                                                        //cout<<"pp.second n'est pas dans map_xor"<<endl;
                                                        continue;
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                    }
                                                }

                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);
                                                                set_quad_4.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_4;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG5_V1)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_4;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    
                                                                    string f = "";

                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ p8";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ p8 ^ p9";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ p9";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ p9 ^ p10";
                                                                            }
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ) ^ p9";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ) ^ p9 ^ p10";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                            }
                                                                        }
                                                                    }
                                                                    
                                                                    temp.formula = f;  

                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            if (!test_quotient_div_l){
                                for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                                    for (uint32_t j3=0; j3<size; j3++){

                                        create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                                        pair<poly,poly> p3;
                                        p3 = poly_div(q, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                        uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                                        uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                                        if (m3 <= 2 && n3 <= 2){
                                        // y = quad * (p3.first * s[i3] + p3.second) * quad  + pp.second

                                            set<poly> quad_sol_imp;

                                            quad_sol_imp.insert(P1);

                                            vector<implem> v1;
                                                
                                            implem imp;
                                            imp.op_sol.push_back(P1);
                                            imp.op_sol.push_back(P1_lin);
                                            imp.quad_sol.insert(pq1);
                                            imp.formula = "";
                                            v1.push_back(imp);

                                            vector<implem> v2;

                                            poly P2;
                                            P2 = truncate_lin(p3.first, size);
                                            poly P2_lin;
                                            POLY_ADD(P2, p3.first, P2_lin);
                                            poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                            int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                            if (indice2 != -1)   {
                                                quad_sol_imp.insert(P2);

                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                            else {
                                                indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                if (indice2 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);

                                                        quad_sol_imp.insert(op1);
                                                        quad_sol_imp.insert(op2);

                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                }
                                            }

                                            vector<implem> v3;

                                            poly P3;
                                            P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                                            poly P3_lin;
                                            POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin);
                                            poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                            implem imp3;
                                            imp3.op_sol.push_back(P3);
                                            imp3.op_sol.push_back(P3_lin);
                                            imp3.quad_sol.insert(pq3);
                                            imp3.formula = "";
                                            v3.push_back(imp3);

                                            quad_sol_imp.insert(P3);

                                            vector<implem> v4;

                                            poly P4;
                                            poly m1;
                                            POLY_MUL(p3.first, set_op_l_plus_lin[i3], m1, nb_elem);
                                            POLY_ADD(m1, q, P4);
                                            poly P4_lin;
                                            POLY_ADD(P4, p3.second, P4_lin);
                                            poly_quad pq4 = poly_to_poly_quad(p3.second ,size, nb_elem);

                                            int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                            if (indice4 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P4);
                                                imp.op_sol.push_back(P4_lin);
                                                imp.quad_sol.insert(pq4);
                                                imp.formula = "";
                                                v4.push_back(imp);

                                                quad_sol_imp.insert(p3.second);
                                            }
                                            else {
                                                indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                if (indice4 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);

                                                        quad_sol_imp.insert(op1);
                                                        quad_sol_imp.insert(op2);

                                                        imp.op_sol.push_back(P4_lin);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                }
                                            }

                                            if (mm == 0){
                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);
                                                                set_quad_4.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_4;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG5_V1)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_4;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    
                                                                    string f = "";
                                                                    if (v2.size()==1){
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 )";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 )";
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 )";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 )";
                                                                        }
                                                                    }
                                                                    temp.formula = f;

                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            else {
                                                vector<implem> v5;

                                                poly_quad pq5 = poly_to_poly_quad(pp.second ,size, nb_elem);
                                                int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                if (indice5 != -1)   {
                                                    implem imp;
                                                    imp.op_sol.push_back(pp.second);
                                                    imp.quad_sol.insert(pq5);
                                                    imp.formula = "";
                                                    v5.push_back(imp);
                                                }
                                                else {
                                                    indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                    if (indice5 == -1)  {
                                                        for (auto it = quad_sol_imp.begin(); it != quad_sol_imp.end(); it++){
                                                            poly new_poly = *it;
                                                            poly test;
                                                            POLY_ADD(new_poly, pp.second,test);
                                                            poly_quad p_test = poly_to_poly_quad(test, size, nb_elem);
                                                            indice5 = find_map(p_test, map_xor, 0, map_xor_size);
                                                            if (indice5 != -1){
                                                                implem imp;
                                                                imp.op_sol.push_back(new_poly);
                                                                poly_quad p_new_poly = poly_to_poly_quad(new_poly, size, nb_elem);
                                                                imp.quad_sol.insert(p_new_poly);
                                                                
                                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                        }
                                                    }
                                                    else {
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.formula = "";
                                                            v5.push_back(imp);
                                                        }
                                                    }
                                                }

                                                for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                                    set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                    for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                        set<poly_quad> set_2 = v2[it2].quad_sol;
                                                        set<poly_quad> set_quad_2 = set_quad_1;
                                                        set_quad_2.merge(set_2);

                                                        for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                            set<poly_quad> set_3 = v3[it3].quad_sol;
                                                            set<poly_quad> set_quad_3 = set_quad_2;
                                                            set_quad_3.merge(set_3);

                                                            for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                set<poly_quad> set_quad_4 = set_quad_3;
                                                                set_quad_4.merge(set_4);

                                                                for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                    set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                    set<poly_quad> set_quad_5 = set_quad_4;
                                                                    set_quad_5.merge(set_5);
                                                                    set_quad_5.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_5;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG5_V1)   {
                                                                        implem temp; 
                                                                        temp.quad_sol = set_quad_5;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        
                                                                        string f = "";

                                                                        if (v2.size()==1){
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9 ^ p10";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11 ^ p12";
                                                                                }
                                                                            }
                                                                        }
                                                                        
                                                                        temp.formula = f;

                                                                        res.push_back(temp);  
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {
                                            if (n3 != 0 && n3 <= 2){
                                                for (uint32_t jj3=0; jj3<size; jj3++){

                                                    if (j3 == jj3){
                                                        continue;
                                                    }

                                                    create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                                    pair<poly,poly> pp3;
                                                    pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                                    uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                                    uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);

                                                    if (mm3 <= 2 && nn3 <= 2){
                                                    // y = quad * ((p3.first + pp3.first) * s[i3] + pp3.second) + pp.second
                                                        vector<implem> v1;
                                                
                                                        implem imp;
                                                        imp.op_sol.push_back(P1);
                                                        imp.op_sol.push_back(P1_lin);
                                                        imp.quad_sol.insert(pq1);
                                                        imp.formula = "";
                                                        v1.push_back(imp);

                                                        vector<implem> v2;

                                                        poly PP2;
                                                        POLY_ADD(p3.first, pp3.first, PP2);
                                                        poly P2 = truncate_lin(PP2, size);
                                                        poly P2_lin;
                                                        POLY_ADD(PP2, P2, P2_lin);
                                                        poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                                        int32_t indice2 = find_set(pq2, set_op, 0, set_op_size);
                                                        if (indice2 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.quad_sol.insert(pq2);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                        else {
                                                            indice2 = find_map(pq2, map_xor, 0, map_xor_size);
                                                            if (indice2 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.op_sol.push_back(P2_lin);
                                                                    imp.formula = "";
                                                                    v2.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        vector<implem> v3;

                                                        poly P3;
                                                        P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                                                        poly P3_lin;
                                                        POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin);
                                                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                                        implem imp3;
                                                        imp3.op_sol.push_back(P3);
                                                        imp3.op_sol.push_back(P3_lin);
                                                        imp3.quad_sol.insert(pq3);
                                                        imp3.formula = "";
                                                        v3.push_back(imp3);

                                                        vector<implem> v4;

                                                        poly P4;
                                                        poly m1;
                                                        POLY_MUL(PP2, set_op_l_plus_lin[i3], m1, nb_elem);
                                                        POLY_ADD(m1, q, P4);
                                                        poly P4_lin;
                                                        POLY_ADD(P4, pp3.second, P4_lin);
                                                        poly_quad pq4 = poly_to_poly_quad(pp3.second ,size, nb_elem);

                                                        int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                                        if (indice4 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P4);
                                                            imp.op_sol.push_back(P4_lin);
                                                            imp.quad_sol.insert(pq4);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                        else {
                                                            indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                            if (indice4 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.op_sol.push_back(P4_lin);
                                                                    imp.formula = "";
                                                                    v4.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        if (mm == 0){
                                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                
                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                                    set_quad_2.merge(set_2);

                                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                                        set_quad_3.merge(set_3);

                                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                                            set_quad_4.merge(set_4);
                                                                            set_quad_4.erase(0);

                                                                            set<poly_quad> s_quad_temp = set_quad_4;
                                                                            uint32_t r = Rank(s_quad_temp);

                                                                            if (r<=MAX_RANK_DEG5_V1)   {
                                                                                implem temp; 
                                                                                temp.quad_sol = set_quad_4;
                                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                
                                                                                string f = "";

                                                                                if (v2.size()==1){
                                                                                    if (v4.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 )";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 )";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v4.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 )";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 )";
                                                                                    }
                                                                                }
                                                                                temp.formula = f;

                                                                                res.push_back(temp);  
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        else {
                                                            vector<implem> v5;

                                                            poly_quad pq5 = poly_to_poly_quad(pp.second ,size, nb_elem);
                                                            int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                            if (indice5 != -1)   {
                                                                implem imp;
                                                                imp.op_sol.push_back(pp.second);
                                                                imp.quad_sol.insert(pq5);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                            else {
                                                                indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                                if (indice5 == -1)  {
                                                                    continue;
                                                                }
                                                                else {
                                                                    for (uint32_t m=0; m<10; m++)   {
                                                                        implem imp;
                                                                        poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                        poly op1(op_q1,size);
                                                                        imp.op_sol.push_back(op1);
                                                                        imp.quad_sol.insert(op_q1);
                                                                        poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                        poly op2(op_q2,size);
                                                                        imp.op_sol.push_back(op2);
                                                                        imp.quad_sol.insert(op_q2);
                                                                        imp.formula = "";
                                                                        v5.push_back(imp);
                                                                    }
                                                                }
                                                            }

                                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                                    
                                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                                    set_quad_2.merge(set_2);

                                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                                        set_quad_3.merge(set_3);

                                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                                            set_quad_4.merge(set_4);

                                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {
                                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                                set_quad_5.merge(set_5);
                                                                                set_quad_5.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_5;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG5_V1)   {
                                                                                    implem temp; 
                                                                                    temp.quad_sol = set_quad_5;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    
                                                                                    string f = "";

                                                                                    if (v2.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ p9 ^ p10";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ) ^ p10 ^ p11";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10 ) ^ p11 ^ p12";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    
                                                                                    temp.formula = f;

                                                                                    res.push_back(temp);  
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (uint32_t i=0; i<l.size(); i++){
        for (uint32_t j=0; j<size; j++){

            if (mon_dom_in(l[i],j, size, nb_elem)){
                create_monomial_order(monomial_order, nb_elem, size, j);

                pair<poly,poly> p1;
                p1 = poly_div(y, l[i], monomial_order, size, nb_elem);

                uint32_t m = p1.second.algebraic_degree(nb_elem);

                if (m<=2)   {
                    for (uint32_t i2=0; i2<set_op_l_plus_lin.size(); i2++){
                        for (uint32_t j2=0; j2<size; j2++){

                            create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                            pair<poly,poly> p2;
                            p2 = poly_div(p1.first, set_op_l_plus_lin[i2], monomial_order, size, nb_elem);

                            uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                            uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                            if (m2 <= 2 && n2 <= 2){
                            // y = (p2.first * s[i2] + p2.second) * l[i] + p.second
                                vector<implem> v1;
                                                
                                poly P1;
                                P1 = truncate_lin(p2.first, size);
                                poly P1_lin;
                                POLY_ADD(P1, p2.first, P1_lin );
                                poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

                                int32_t indice1 = find_set(pq1, set_op, 0, set_op_size);
                                if (indice1 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(P1);
                                    imp.op_sol.push_back(P1_lin);
                                    imp.quad_sol.insert(pq1);
                                    imp.formula = "";
                                    v1.push_back(imp);
                                }
                                else {
                                    indice1 = find_map(pq1, map_xor, 0, map_xor_size);
                                    if (indice1 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice1].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice1].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P1_lin);
                                            imp.formula = "";
                                            v1.push_back(imp);
                                        }
                                    }
                                }

                                vector<implem> v2;

                                poly P2;
                                P2 = truncate_lin(set_op_l_plus_lin[i2], size);
                                poly P2_lin;
                                POLY_ADD(P2, set_op_l_plus_lin[i2] , P2_lin);
                                poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);

                                implem imp;
                                imp.op_sol.push_back(P2);
                                imp.op_sol.push_back(P2_lin);
                                imp.op_sol.push_back(l[i2]);
                                imp.quad_sol.insert(pq2);
                                imp.formula = "";
                                v2.push_back(imp);

                                vector<implem> v3;

                                poly P3;
                                poly m1;
                                POLY_MUL(p2.first, set_op_l_plus_lin[i2], m1, nb_elem);
                                POLY_ADD(m1, p1.first, P3);
                                poly P3_lin;
                                POLY_ADD(P3, p2.second, P3_lin);
                                poly_quad pq3 = poly_to_poly_quad(p2.second ,size, nb_elem);

                                int32_t indice3 = find_set(pq3, set_op, 0, set_op_size);
                                if (indice3 != -1)   {
                                    implem imp;
                                    imp.op_sol.push_back(p2.second);
                                    imp.op_sol.push_back(P3_lin);
                                    imp.op_sol.push_back(l[i]);
                                    imp.quad_sol.insert(pq3);
                                    imp.formula = "";
                                    v3.push_back(imp);
                                }
                                else {
                                    indice3 = find_map(pq3, map_xor, 0, map_xor_size);
                                    if (indice3 == -1)  {
                                        continue;
                                    }
                                    else {
                                        for (uint32_t m=0; m<10; m++)   {
                                            implem imp;
                                            poly_quad op_q1 = map_xor[indice3].second[m][0];
                                            poly op1(op_q1,size);
                                            imp.op_sol.push_back(op1);
                                            imp.quad_sol.insert(op_q1);
                                            poly_quad op_q2 = map_xor[indice3].second[m][1];
                                            poly op2(op_q2,size);
                                            imp.op_sol.push_back(op2);
                                            imp.quad_sol.insert(op_q2);
                                            imp.op_sol.push_back(P3_lin);
                                            imp.op_sol.push_back(l[i]);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                    }
                                }

                                if (m==0){
                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);
                                                set_quad_3.erase(0);

                                                set<poly_quad> s_quad_temp = set_quad_3;
                                                uint32_t r = Rank(s_quad_temp);

                                                if (r<=MAX_RANK_DEG5_V1)   {
                                                    implem temp; 
                                                    temp.quad_sol = set_quad_3;
                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                    string f = "";
                                                    if (v1.size()==1){
                                                        if (v3.size()==1){
                                                            f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6) * (p7)";
                                                        }
                                                        else {
                                                            f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6 ^ p7) * (p8)";
                                                        }
                                                    }
                                                    else {
                                                        if (v3.size()==1){
                                                            f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7) * (p8)";
                                                        }
                                                        else {
                                                            f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7 ^ p8) * (p9)";
                                                        }
                                                    }
                                                    temp.formula = f;
                                                    res.push_back(temp);  
                                                }
                                            }
                                        }
                                    }
                                }
                                else {
                                    vector<implem> v4;

                                    poly_quad pq4 = poly_to_poly_quad(p1.second ,size, nb_elem);
                                    int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                    if (indice4 != -1)   {
                                        implem imp;
                                        imp.op_sol.push_back(p1.second);
                                        imp.quad_sol.insert(pq4);
                                        imp.formula = "";
                                        v4.push_back(imp);
                                    }
                                    else {
                                        indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                        if (indice4 == -1)  {
                                            continue;
                                        }
                                        else {
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                        }
                                    }

                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                            set<poly_quad> set_quad_2 = set_quad_1;
                                            set_quad_2.merge(set_2);

                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                set_quad_3.merge(set_3);

                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                    set_quad_4.merge(set_4);
                                                    set_quad_4.erase(0);

                                                    set<poly_quad> s_quad_temp = set_quad_4;
                                                    uint32_t r = Rank(s_quad_temp);

                                                    if (r<=MAX_RANK_DEG5_V1)   {
                                                        implem temp; 
                                                        temp.quad_sol = set_quad_4;
                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                        string f = "";

                                                        if (v1.size()==1){
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6) * (p7) ^ p8";
                                                                }
                                                                else {
                                                                    f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6) * (p7) ^ p8 ^ p9";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6 ^ p7) * (p8) ^ p9";
                                                                }
                                                                else {
                                                                    f = "( (p1 ^ p2) * (p3 ^ p4) ^ p5 ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                }
                                                            }
                                                        }
                                                        else {
                                                            if (v3.size()==1){
                                                                if (v4.size()==1){
                                                                    f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7) * (p8) ^ p9";
                                                                }
                                                                else {
                                                                    f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7) * (p8) ^ p9 ^ p10";
                                                                }
                                                            }
                                                            else {
                                                                if (v4.size()==1){
                                                                    f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7 ^ p8) * (p9) ^ p10";
                                                                }
                                                                else {
                                                                    f = "( (p1 ^ p2 ^ p3) * (p4 ^ p5) ^ p6 ^ p7 ^ p8) * (p9) ^ p10 ^ p11";
                                                                }
                                                            }
                                                        }
                                                        
                                                        temp.formula = f;         
                                                        res.push_back(temp);  
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    res = remove_duplicates(res);

    if (verbose){
        cout<<" Thread : "<<omp_get_thread_num()<<" res_size v1 = "<<res.size()<<endl;
    }

    return res;
}

vector<implem> add_to_op_selec_deg_5_v2(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size)   {
    vector<implem> res;
    uint32_t monomial_order [nb_elem];

    poly y5 = truncate_except_deg(y, 5, size, nb_elem);
    poly y4 = truncate_except_deg(y, 4, size, nb_elem);
    poly y3 = truncate_except_deg(y, 3, size, nb_elem);
    poly y2 = truncate_except_deg(y, 2, size, nb_elem);

    poly y5_4;
    POLY_ADD(y5,y4,y5_4);

    poly y3_2;
    POLY_ADD(y3,y2,y3_2);

    uint32_t ent_y5_4 = poly_deg5_4_to_uint(y5_4);

    vector<implem_deg5_4> vect_sol_deg5_4 = parse_file_and_create_sol(ent_y5_4);

    for (uint32_t i=0; i<vect_sol_deg5_4.size(); i++){

        poly quotient = vect_sol_deg5_4[i].q;
        poly test_quotient = truncate_except_deg(quotient, 3, size, nb_elem);

        poly reste;
        POLY_ADD(vect_sol_deg5_4[i].r, y3_2, reste);
        poly test_reste = truncate_except_deg(reste, 3, size, nb_elem);

        poly P1;
        P1 = truncate_lin(vect_sol_deg5_4[i].quad, size);
        poly P1_lin;
        POLY_ADD(P1, vect_sol_deg5_4[i].quad, P1_lin);
        poly_quad pq1 = poly_to_poly_quad(P1 ,size, nb_elem);

        bool test_quotient_div_l = false;

        for (uint32_t i2=0; i2<l.size(); i2++){
            for (uint32_t j2=0; j2<size; j2++){

                if (mon_dom_in(l[i2],j2, size, nb_elem)){
                    create_monomial_order(monomial_order, nb_elem, size, j2);                                                   
                    pair<poly,poly> p2;
                    p2 = poly_div(quotient, l[i2], monomial_order, size, nb_elem);

                    uint32_t m2 = p2.second.algebraic_degree(nb_elem);
                    uint32_t n2 = p2.first.algebraic_degree(nb_elem);

                    if (m2 <= 2 && n2 <= 2){
                    //y5_4 = (p2.first * l[i2] + p2.second) * .quad + reste

                        poly P2;
                        P2 = truncate_lin(p2.first, size);
                        poly P2_lin;
                        POLY_ADD(P2, p2.first, P2_lin);
                        poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
                        int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                        int32_t indice_mP2 = -1;
                        if (indiceP2 == -1){
                            indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                            if (indice_mP2 == -1){
                                continue;
                            }
                        }

                        poly_quad pq3; 
                        int32_t indiceP3 = -1;
                        int32_t indice_mP3 = -1;
                        if (m2 == 2){
                            pq3 = poly_to_poly_quad(p2.second ,size, nb_elem);
                            indiceP3 = find_set(pq3, set_op, 0, set_op_size);
                            if (indiceP3 == -1){
                                indice_mP3 = find_map(pq3, map_xor, 0, map_xor_size);
                                if (indice_mP3 == -1){
                                        continue;
                                }
                            }
                        }

                        poly P3;
                        poly m1;
                        POLY_MUL(p2.first, l[i2], m1, nb_elem);
                        POLY_ADD(m1, quotient, P3);
                        poly P3_lin;
                        POLY_ADD(P3, p2.second, P3_lin);

                        bool test_reste_div_par_l = false;
                        //bool test_reste_div_par_l = true;

                        for (uint32_t i3=0; i3<l.size(); i3++){
                            for (uint32_t j3=0; j3<size; j3++){

                                if (mon_dom_in(l[i3],j3, size, nb_elem)){

                                    create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                                    pair<poly,poly> p3;
                                    p3 = poly_div(reste, l[i3], monomial_order, size, nb_elem);

                                    uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                                    uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                                    if (m3<=2 && n3<=2){
                                    // y = (P2 * l[i2] + P3) * P1 + (P4 * l[i3] + P5)
                                    // y = .quad * (p2.first * l[i2] + p2.second) + p3.first * l[i3] + p3.second
                                        test_quotient_div_l = true;
                                        test_reste_div_par_l = true;

                                        vector<implem> v1;
                                                
                                        implem imp;
                                        imp.op_sol.push_back(P1);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.quad_sol.insert(pq1);
                                        imp.formula = "";
                                        v1.push_back(imp);

                                        vector<implem> v2;

                                        int32_t indice2 = indiceP2;
                                        if (indice2 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.quad_sol.insert(pq2);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                        }
                                        else {
                                                indice2 = indice_mP2;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P2_lin);
                                                    imp.op_sol.push_back(l[i2]);
                                                    imp.formula = "";
                                                    v2.push_back(imp);
                                                }
                                        }

                                        vector<implem> v3;

                                        if (m2 == 0){
                                            implem imp;
                                            imp.op_sol.push_back(P3);
                                            imp.quad_sol.insert(0);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                        else {
                                            int32_t indice3 = indiceP3;
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P3);
                                                imp.op_sol.push_back(P3_lin);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = indice_mP3;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P3_lin);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v4;

                                        poly P4;
                                        P4 = truncate_lin(p3.first, size);
                                        poly P4_lin;
                                        POLY_ADD(P4, p3.first, P4_lin);
                                        poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);
                                        int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                        if (indice4 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P4);
                                                imp.op_sol.push_back(P4_lin);
                                                imp.op_sol.push_back(l[i3]);
                                                imp.quad_sol.insert(pq4);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                        }
                                        else {
                                                indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                if (indice4 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.op_sol.push_back(P4_lin);
                                                        imp.op_sol.push_back(l[i3]);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                }
                                        }

                                        if (m3 == 0){
                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                    
                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);
                                                            set_quad_4.erase(0);

                                                            set<poly_quad> s_quad_temp = set_quad_4;
                                                            uint32_t r = Rank(s_quad_temp);

                                                            if (r<=MAX_RANK_DEG5_V2)   {
                                                                implem temp; 
                                                                temp.quad_sol = set_quad_4;
                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    
                                                                string f = "";

                                                                if (v2.size()==1){
                                                                    if (v3.size()==1){
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9) * (p10)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9 ^ p10) * (p11)";
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p8 ^ p9) * (p10)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p9 ^ p10 ^ p11) * (p12)";
                                                                        }
                                                                    }
                                                                }
                                                                else {
                                                                    if (v3.size()==1){
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12)";
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v4.size()==1){
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11) * (p12)";
                                                                        }
                                                                        else {
                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13)";
                                                                        }
                                                                    }
                                                                }
                                                                    
                                                                temp.formula = f;
                                                                res.push_back(temp);  
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        } 
                                            
                                        else {

                                            vector<implem> v5;

                                            poly_quad pq5 = poly_to_poly_quad(p3.second ,size, nb_elem);
                                            int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                            if (indice5 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p3.second);
                                                imp.quad_sol.insert(pq5);
                                                imp.formula = "";
                                                v5.push_back(imp);
                                            }
                                            else {
                                                indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                if (indice5 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                }
                                            }

                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {
                                                        
                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {
                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {
                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {
                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);
                                                                set_quad_5.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_5;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                    implem temp; 
                                                                    temp.quad_sol = set_quad_5;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());

                                                                    string f = "";
                                                                        
                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10) ) ^ p11";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10) ) ^ p11 ^ p12";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11) ) ^ p12";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11) ) ^ p12 ^ p13";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10) * (p11) ) ^ p12";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10) * (p11) ) ^ p12 ^ p13";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12 ^ p13";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                }
                                                                            }
                                                                        }
                                                                    }                                                                        
                                                                       
                                                                    temp.formula = f;
                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }                                                              
                                        }
                                    }
                                }
                            }
                        }

                        if (!test_reste_div_par_l){
                            for (uint32_t i4=0; i4<set_op_l_plus_lin.size(); i4++){
                                for (uint32_t j4=0; j4<size; j4++){

                                    create_monomial_order(monomial_order, nb_elem, size, j4);                                                   
                                    pair<poly,poly> p4;
                                    p4 = poly_div(reste, set_op_l_plus_lin[i4], monomial_order, size, nb_elem);

                                    uint32_t m4 = p4.second.algebraic_degree(nb_elem);
                                    uint32_t n4 = p4.first.algebraic_degree(nb_elem);

                                    if (m4 <= 2 && n4 <= 2){
                                    // y = (P2 * l[i2] + P3) * P1 + (P4 * P5 + P6)
                                    //y = quad * (p2.first * l[i2] + p2.second) + p4.first * s[i4] + p4.second
                                    test_quotient_div_l = true;

                                        vector<implem> v1;
                                            
                                        implem imp;
                                        imp.op_sol.push_back(P1);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.quad_sol.insert(pq1);
                                        imp.formula = "";
                                        v1.push_back(imp);

                                        vector<implem> v2;

                                        int32_t indice2 = indiceP2;
                                        if (indice2 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.op_sol.push_back(l[i2]);
                                            imp.quad_sol.insert(pq2);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                        else {
                                            indice2 = indice_mP2;
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.op_sol.push_back(l[i2]);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                        }

                                        vector<implem> v3;

                                        if (m2 == 0){
                                            implem imp;
                                            imp.op_sol.push_back(P3);
                                            imp.quad_sol.insert(0);
                                            imp.formula = "";
                                            v3.push_back(imp);
                                        }
                                         else {
                                            int32_t indice3 = indiceP3;
                                            if (indice3 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P3);
                                                imp.op_sol.push_back(P3_lin);
                                                imp.quad_sol.insert(pq3);
                                                imp.formula = "";
                                                v3.push_back(imp);
                                            }
                                            else {
                                                indice3 = indice_mP3;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P3_lin);
                                                    imp.formula = "";
                                                    v3.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v4;

                                        poly P4;
                                        P4 = truncate_lin(p4.first, size);
                                        poly P4_lin;
                                        POLY_ADD(P4, p4.first, P4_lin);
                                        poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);
                                        int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                        if (indice4 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P4);
                                            imp.op_sol.push_back(P4_lin);
                                            imp.quad_sol.insert(pq4);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                        else {
                                            indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                            if (indice4 == -1)  {
                                                continue;
                                            }
                                            else {
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P4_lin);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v5;

                                        poly P5;
                                        P5 = truncate_lin(set_op_l_plus_lin[i4], size);
                                        poly P5_lin;
                                        POLY_ADD(P5, set_op_l_plus_lin[i4], P5_lin)
                                        poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                        implem imp5;
                                        imp5.op_sol.push_back(P5);
                                        imp5.op_sol.push_back(P5_lin);
                                        imp5.quad_sol.insert(pq5);
                                        imp5.formula = "";
                                        v5.push_back(imp5);

                                        if (m4 == 0){
                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);
                                                                set_quad_5.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_5;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                    implem temp;
                                                                    temp.quad_sol = set_quad_5;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                    
                                                                    string f = "";

                                                                    if (v2.size()==1){
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9) * (p10 ^ p11)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9 ^ p10) * (p11 ^ p12)";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p8 ^ p9) * (p10 ^ p11)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                            }
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v3.size()==1){
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                            }
                                                                        }
                                                                    }
                                                                        
                                                                    temp.formula = f;
                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {

                                            vector<implem> v6;

                                            poly_quad pq6 = poly_to_poly_quad(p4.second ,size, nb_elem);
                                            int32_t indice6 = find_set(pq6, set_op, 0, set_op_size);
                                            if (indice6 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p4.second);
                                                imp.quad_sol.insert(pq6);
                                                imp.formula = "";
                                                v6.push_back(imp);
                                            }
                                            else {
                                                indice6 = find_map(pq6, map_xor, 0, map_xor_size);
                                                if (indice6 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice6].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice6].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v6.push_back(imp);
                                                    }
                                                }
                                            }

                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);

                                                                for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                    set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                    set<poly_quad> set_quad_6 = set_quad_5;
                                                                    set_quad_6.merge(set_6);
                                                                    set_quad_6.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_6;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG5_V2)   {
                                                                        implem temp;
                                                                        temp.quad_sol = set_quad_6;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());

                                                                        string f = "";
                                                                        
                                                                        if (v2.size()==1){
                                                                            if (v3.size()==1){
                                                                                if (v4.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12 ^ p13";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v4.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12 ^ p13";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v3.size()==1){
                                                                                if (v4.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v4.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                    }
                                                                                }
                                                                            }
                                                                        }                                                                        
                                                                        
                                                                        temp.formula = f; 
                                                                        res.push_back(temp);  
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else {
                                    //On fait la double division sur la partie degré 3
                                        if (n4 != 0 && n4 <= 2){
                                            for (uint32_t jj4=0; jj4<size; jj4++){
                                                if (jj4 == j4){
                                                    continue;
                                                }

                                                create_monomial_order(monomial_order, nb_elem, size, jj4);                                                   
                                                pair<poly,poly> pp4;
                                                pp4 = poly_div(p4.second, set_op_l_plus_lin[i4], monomial_order, size, nb_elem);

                                                uint32_t mm4 = pp4.second.algebraic_degree(nb_elem);
                                                uint32_t nn4 = pp4.first.algebraic_degree(nb_elem);

                                                if (mm4 <= 2 && nn4 <= 2){
                                                test_quotient_div_l = true;
                                                //y = (P2 * l[i2] + P3) * P1 + (PP4 * P5 + P6)
                                                //y = (p2.first * l[i2] + p2.second) * .quad + (p4.first + pp4.first) * s[i4] + pp4.second
                                                    vector<implem> v1;
                                            
                                                    implem imp;
                                                    imp.op_sol.push_back(P1);
                                                    imp.op_sol.push_back(P1_lin);
                                                    imp.quad_sol.insert(pq1);
                                                    imp.formula = "";
                                                    v1.push_back(imp);

                                                    vector<implem> v2;

                                                    int32_t indice2 = indiceP2;
                                                    if (indice2 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.op_sol.push_back(l[i2]);
                                                        imp.quad_sol.insert(pq2);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                    else {
                                                        indice2 = indice_mP2;
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.op_sol.push_back(l[i2]);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                    }

                                                    vector<implem> v3;

                                                    if (m2 == 0){
                                                        implem imp;
                                                        imp.op_sol.push_back(P3);
                                                        imp.quad_sol.insert(0);
                                                        imp.formula = "";
                                                        v3.push_back(imp);
                                                    }
                                                    else {
                                                        int32_t indice3 = indiceP3;
                                                        if (indice3 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P3);
                                                            imp.op_sol.push_back(P3_lin);
                                                            imp.quad_sol.insert(pq3);
                                                            imp.formula = "";
                                                            v3.push_back(imp);
                                                        }
                                                        else {
                                                            indice3 = indice_mP3;
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice3].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice3].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P3_lin);
                                                                imp.formula = "";
                                                                v3.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v4;

                                                    poly PP4;
                                                    POLY_ADD(p4.first, pp4.first, PP4);
                                                    poly P4;
                                                    P4 = truncate_lin(PP4, size);
                                                    poly P4_lin;
                                                    POLY_ADD(PP4, P4, P4_lin);
                                                    poly_quad pq4 = poly_to_poly_quad(P4 ,size, nb_elem);
                                                    int32_t indice4 = find_set(pq4, set_op, 0, set_op_size);
                                                    if (indice4 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P4);
                                                        imp.op_sol.push_back(P4_lin);
                                                        imp.quad_sol.insert(pq4);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                    else {
                                                        indice4 = find_map(pq4, map_xor, 0, map_xor_size);
                                                        if (indice4 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P4_lin);
                                                                imp.formula = "";
                                                                v4.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v5;

                                                    poly P5;
                                                    P5 = truncate_lin(set_op_l_plus_lin[i4], size);
                                                    poly P5_lin;
                                                    POLY_ADD(P5, set_op_l_plus_lin[i4], P5_lin);
                                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                                    implem imp5;
                                                    imp5.op_sol.push_back(P5);
                                                    imp5.op_sol.push_back(P5_lin);
                                                    imp5.quad_sol.insert(pq5);
                                                    imp5.formula = "";
                                                    v5.push_back(imp5);

                                                    if (mm4 == 0){
                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);
                                                                            set_quad_5.erase(0);

                                                                            set<poly_quad> s_quad_temp = set_quad_5;
                                                                            uint32_t r = Rank(s_quad_temp);

                                                                            if (r<=MAX_RANK_DEG5_V2)   {
                                                                                implem temp;
                                                                                temp.quad_sol = set_quad_5;
                                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                
                                                                                string f = "";

                                                                                if (v2.size()==1){
                                                                                    if (v3.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9) * (p10 ^ p11)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ (p8 ^ p9 ^ p10) * (p11 ^ p12)";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p8 ^ p9) * (p10 ^ p11)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v3.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                    
                                                                                temp.formula = f;

                                                                                res.push_back(temp);  
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }

                                                    else {

                                                        vector<implem> v6;

                                                        poly_quad pq6 = poly_to_poly_quad(pp4.second ,size, nb_elem);
                                                        int32_t indice6 = find_set(pq6, set_op, 0, set_op_size);
                                                        if (indice6 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(pp4.second);
                                                            imp.quad_sol.insert(pq6);
                                                            imp.formula = "";
                                                            v6.push_back(imp);
                                                        }
                                                        else {
                                                            indice6 = find_map(pq6, map_xor, 0, map_xor_size);
                                                            if (indice6 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice6].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice6].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.formula = "";
                                                                    v6.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);
                                                                                set_quad_6.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_6;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                                    implem temp;
                                                                                    temp.quad_sol = set_quad_6;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                    
                                                                                    string f = "";
                                                                        
                                                                                    if (v2.size()==1){
                                                                                        if (v3.size()==1){
                                                                                            if (v4.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12 ^ p13";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7) ^ ( (p8 ^ p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v4.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p8 ^ p9) * (p10 ^ p11) ) ^ p12 ^ p13";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5) ^ p6 ^ p7 ^ p8 ) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v3.size()==1){
                                                                                            if (v4.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v4.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6) ^ p7 ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }                                                                        
                                                                                    
                                                                                    temp.formula = f; 

                                                                                    res.push_back(temp);  
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!test_quotient_div_l){
            for (uint32_t i3=0; i3<set_op_l_plus_lin.size(); i3++){
                for (uint32_t j3=0; j3<size; j3++){

                    create_monomial_order(monomial_order, nb_elem, size, j3);                                                   
                    pair<poly,poly> p3;
                    p3 = poly_div(quotient, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                    uint32_t m3 = p3.second.algebraic_degree(nb_elem);
                    uint32_t n3 = p3.first.algebraic_degree(nb_elem);

                    if (m3 <= 2 && n3 <= 2){
                    //y5_4 = (p3.first * s[i3] + p3.second) * .quad + reste
                    //y5_4 = (P2 * P3 + P4) * P1 + reste

                        poly P2;
                        P2 = truncate_lin(p3.first, size);
                        poly P2_lin;
                        POLY_ADD(P2, p3.first, P2_lin);
                        poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
                        int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                        int32_t indice_mP2 = -1;
                        if (indiceP2 == -1){
                            indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                            if (indice_mP2 == -1){
                                continue;
                            }
                        }

                        poly P3;
                        P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                        poly P3_lin;
                        POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin);
                        poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                        poly_quad pq4;
                        int32_t indiceP4 = -1;
                        int32_t indice_mP4 = -1;
                        if (m3 != 0){   
                            pq4 = poly_to_poly_quad(p3.second ,size, nb_elem);
                            indiceP4 = find_set(pq4, set_op, 0, set_op_size);
                            if (indiceP4 == -1){
                                indice_mP4 = find_map(pq4, map_xor, 0, map_xor_size);
                                if (indice_mP4 == -1){
                                        continue;
                                }
                            }
                        }

                        poly P4;
                        poly m1;
                        POLY_MUL(p3.first, set_op_l_plus_lin[i3], m1, nb_elem);
                        POLY_ADD(m1, quotient, P4);
                        poly P4_lin;
                        POLY_ADD(P4, p3.second, P4_lin);

                        bool test_reste_div_par_l = false;
                        //bool test_reste_div_par_l = true;

                        for (uint32_t i4=0; i4<l.size(); i4++){
                            for (uint32_t j4=0; j4<size; j4++){

                                if (mon_dom_in(l[i4],j4, size, nb_elem)){

                                    create_monomial_order(monomial_order, nb_elem, size, j4);                                                   
                                    pair<poly,poly> p4;
                                    p4 = poly_div(reste, l[i4], monomial_order, size, nb_elem);

                                    uint32_t m4 = p4.second.algebraic_degree(nb_elem);
                                    uint32_t n4 = p4.first.algebraic_degree(nb_elem);

                                    if (m4 <= 2 && n4 <= 2){
                                    // y = ( (P2 * P3) + P4) * P1 + (P5 * l[i4]) + P6
                                    test_reste_div_par_l = true;

                                        vector<implem> v1;
                                        
                                        implem imp1;
                                        imp1.op_sol.push_back(P1);
                                        imp1.op_sol.push_back(P1_lin);
                                        imp1.quad_sol.insert(pq1);
                                        imp1.formula = "";
                                        v1.push_back(imp1);

                                        vector<implem> v2;

                                        int32_t indice2 = indiceP2;
                                        if (indice2 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.quad_sol.insert(pq2);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                        else {
                                            indice2 = indice_mP2;
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                        }

                                        vector<implem> v3;

                                        implem imp3;
                                        imp3.op_sol.push_back(P3);
                                        imp3.op_sol.push_back(P3_lin);
                                        imp3.quad_sol.insert(pq3);
                                        imp3.formula = "";
                                        v3.push_back(imp3);

                                        vector<implem> v4;

                                        if (m3 == 0){
                                            implem imp;
                                            imp.op_sol.push_back(P4);
                                            imp.quad_sol.insert(0);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                        else {
                                            int32_t indice4 = indiceP4;
                                            if (indice4 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P4);
                                                imp.op_sol.push_back(P4_lin);
                                                imp.quad_sol.insert(pq4);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                            else {
                                                indice4 = indice_mP4;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P4_lin);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v5;

                                        poly P5;
                                        P5 = truncate_lin(p4.first, size);
                                        poly P5_lin;
                                        POLY_ADD(p4.first, P5, P5_lin);
                                        poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                        int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                        if (indice5 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P5);
                                            imp.op_sol.push_back(P5_lin);
                                            imp.op_sol.push_back(l[i4]);
                                            imp.quad_sol.insert(pq5);
                                            imp.formula = "";
                                            v5.push_back(imp);
                                        }
                                        else {
                                            indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                            if (indice5 == -1)  {
                                                continue;
                                            }
                                            else {
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P5_lin);
                                                    imp.op_sol.push_back(l[i4]);
                                                    imp.formula = "";
                                                    v5.push_back(imp);
                                                }
                                            }
                                        }

                                        if (m4 == 0){
                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);
                                                                set_quad_5.erase(0);

                                                                set<poly_quad> s_quad_temp = set_quad_5;
                                                                uint32_t r = Rank(s_quad_temp);

                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                    implem temp;
                                                                    temp.quad_sol = set_quad_5;
                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                    
                                                                    string f = "";

                                                                    if (v2.size()==1){
                                                                        if (v4.size()==1){
                                                                            if (v5.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12)";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v5.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13)";
                                                                            }
                                                                        }
                                                                    }
                                                                    else {
                                                                        if (v4.size()==1){
                                                                            if (v5.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13)";
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v5.size()==1){
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13)";
                                                                            }
                                                                            else {
                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14)";
                                                                            }
                                                                        }
                                                                    }
                                                                        
                                                                    temp.formula = f;

                                                                    res.push_back(temp);  
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {

                                            vector<implem> v6;

                                            poly_quad pq6 = poly_to_poly_quad(p4.second ,size, nb_elem);
                                            int32_t indice6 = find_set(pq6, set_op, 0, set_op_size);
                                            if (indice6 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p4.second);
                                                imp.quad_sol.insert(pq6);
                                                imp.formula = "";
                                                v6.push_back(imp);
                                            }
                                            else {
                                                indice6 = find_map(pq6, map_xor, 0, map_xor_size);
                                                if (indice6 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice6].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice6].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v6.push_back(imp);
                                                    }
                                                }
                                            }

                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);

                                                                for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                    set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                    set<poly_quad> set_quad_6 = set_quad_5;
                                                                    set_quad_6.merge(set_6);
                                                                    set_quad_6.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_6;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG5_V2)   {
                                                                        implem temp;
                                                                        temp.quad_sol = set_quad_6;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                        
                                                                        string f = "";
                                                                        
                                                                        if (v2.size()==1){
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12 ^ p13";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12) ) ^ p13";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13) ) ^ p14";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v6.size()==1){
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14) ) ^ p15";
                                                                                    }
                                                                                    else {
                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14) ) ^ p15 ^ p16";
                                                                                    }
                                                                                }
                                                                            }
                                                                        }                                                                        
                                                                        
                                                                        temp.formula = f;

                                                                        res.push_back(temp);  
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (!test_reste_div_par_l){

                            for (uint32_t i5=0; i5<set_op_l_plus_lin.size(); i5++){
                                for (uint32_t j5=0; j5<size; j5++){

                                    create_monomial_order(monomial_order, nb_elem, size, j5);                                                   
                                    pair<poly,poly> p5;
                                    p5 = poly_div(reste, set_op_l_plus_lin[i5], monomial_order, size, nb_elem);

                                    uint32_t m5 = p5.second.algebraic_degree(nb_elem);
                                    uint32_t n5 = p5.first.algebraic_degree(nb_elem);

                                    if (m5 <= 2 && n5 <= 2){
                                    // y = ( (P2 * P3) + P4) * P1 + (P5 * P6) + P7
                                        vector<implem> v1;
                                        
                                        implem imp;
                                        imp.op_sol.push_back(P1);
                                        imp.op_sol.push_back(P1_lin);
                                        imp.quad_sol.insert(pq1);
                                        imp.formula = "";
                                        v1.push_back(imp);

                                        vector<implem> v2;

                                        int32_t indice2 = indiceP2;
                                        if (indice2 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P2);
                                            imp.op_sol.push_back(P2_lin);
                                            imp.quad_sol.insert(pq2);
                                            imp.formula = "";
                                            v2.push_back(imp);
                                        }
                                        else {
                                            indice2 = indice_mP2;
                                            for (uint32_t m=0; m<10; m++)   {
                                                implem imp;
                                                poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                poly op1(op_q1,size);
                                                imp.op_sol.push_back(op1);
                                                imp.quad_sol.insert(op_q1);
                                                poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                poly op2(op_q2,size);
                                                imp.op_sol.push_back(op2);
                                                imp.quad_sol.insert(op_q2);
                                                imp.op_sol.push_back(P2_lin);
                                                imp.formula = "";
                                                v2.push_back(imp);
                                            }
                                        }

                                        vector<implem> v3;

                                        implem imp3;
                                        imp3.op_sol.push_back(P3);
                                        imp3.op_sol.push_back(P3_lin);
                                        imp3.quad_sol.insert(pq3);
                                        imp3.formula = "";
                                        v3.push_back(imp3);

                                        vector<implem> v4;

                                        if (m3 == 0){
                                            implem imp;
                                            imp.op_sol.push_back(P4);
                                            imp.quad_sol.insert(0);
                                            imp.formula = "";
                                            v4.push_back(imp);
                                        }
                                        else {
                                            int32_t indice4 = indiceP4;
                                            if (indice4 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(P4);
                                                imp.op_sol.push_back(P4_lin);
                                                imp.quad_sol.insert(pq4);
                                                imp.formula = "";
                                                v4.push_back(imp);
                                            }
                                            else {
                                                indice4 = indice_mP4;
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P4_lin);
                                                    imp.formula = "";
                                                    v4.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v5;

                                        poly P5;
                                        P5 = truncate_lin(p5.first, size);
                                        poly P5_lin;
                                        POLY_ADD(P5, p5.first, P5_lin);
                                        poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                        int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                        if (indice5 != -1)   {
                                            implem imp;
                                            imp.op_sol.push_back(P5);
                                            imp.op_sol.push_back(P5_lin);
                                            imp.quad_sol.insert(pq5);
                                            imp.formula = "";
                                            v5.push_back(imp);
                                        }
                                        else {
                                            indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                            if (indice5 == -1)  {
                                                continue;
                                            }
                                            else {
                                                for (uint32_t m=0; m<10; m++)   {
                                                    implem imp;
                                                    poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                    poly op1(op_q1,size);
                                                    imp.op_sol.push_back(op1);
                                                    imp.quad_sol.insert(op_q1);
                                                    poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                    poly op2(op_q2,size);
                                                    imp.op_sol.push_back(op2);
                                                    imp.quad_sol.insert(op_q2);
                                                    imp.op_sol.push_back(P5_lin);
                                                    imp.formula = "";
                                                    v5.push_back(imp);
                                                }
                                            }
                                        }

                                        vector<implem> v6;

                                        poly P6;
                                        P6 = truncate_lin(set_op_l_plus_lin[i5], size);
                                        poly P6_lin;
                                        POLY_ADD(P6, set_op_l_plus_lin[i5], P6_lin);
                                        poly_quad pq6 = poly_to_poly_quad(P6 ,size, nb_elem);
                                        implem imp6;
                                        imp6.op_sol.push_back(P6);
                                        imp6.op_sol.push_back(P6_lin);
                                        imp6.quad_sol.insert(pq6);
                                        imp6.formula = "";
                                        v6.push_back(imp6);

                                        if (m5 == 0){
                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);

                                                                for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                    set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                    set<poly_quad> set_quad_6 = set_quad_5;
                                                                    set_quad_6.merge(set_6);
                                                                    set_quad_6.erase(0);

                                                                    set<poly_quad> s_quad_temp = set_quad_6;
                                                                    uint32_t r = Rank(s_quad_temp);

                                                                    if (r<=MAX_RANK_DEG5_V2)   {
                                                                        implem temp;
                                                                        temp.quad_sol = set_quad_6;
                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                        temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                        
                                                                        string f = "";

                                                                        if (v2.size()==1){
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                }
                                                                            }
                                                                        }
                                                                        else {
                                                                            if (v4.size()==1){
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v5.size()==1){
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13 ^ p14)";
                                                                                }
                                                                                else {
                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14 ^ p15)";
                                                                                }
                                                                            }
                                                                        }
                                                                            
                                                                        temp.formula = f;

                                                                        res.push_back(temp);  
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        else {

                                            vector<implem> v7;

                                            poly_quad pq7 = poly_to_poly_quad(p5.second ,size, nb_elem);
                                            int32_t indice7 = find_set(pq7, set_op, 0, set_op_size);
                                            if (indice7 != -1)   {
                                                implem imp;
                                                imp.op_sol.push_back(p5.second);
                                                imp.quad_sol.insert(pq7);
                                                imp.formula = "";
                                                v7.push_back(imp);
                                            }
                                            else {
                                                indice7 = find_map(pq7, map_xor, 0, map_xor_size);
                                                if (indice7 == -1)  {
                                                    continue;
                                                }
                                                else {
                                                    for (uint32_t m=0; m<10; m++)   {
                                                        implem imp;
                                                        poly_quad op_q1 = map_xor[indice7].second[m][0];
                                                        poly op1(op_q1,size);
                                                        imp.op_sol.push_back(op1);
                                                        imp.quad_sol.insert(op_q1);
                                                        poly_quad op_q2 = map_xor[indice7].second[m][1];
                                                        poly op2(op_q2,size);
                                                        imp.op_sol.push_back(op2);
                                                        imp.quad_sol.insert(op_q2);
                                                        imp.formula = "";
                                                        v7.push_back(imp);
                                                    }
                                                }
                                            }

                                            for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                    set<poly_quad> set_2 = v2[it2].quad_sol;
                                                    set<poly_quad> set_quad_2 = set_quad_1;
                                                    set_quad_2.merge(set_2);

                                                    for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                        set<poly_quad> set_3 = v3[it3].quad_sol;
                                                        set<poly_quad> set_quad_3 = set_quad_2;
                                                        set_quad_3.merge(set_3);

                                                        for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                            set<poly_quad> set_4 = v4[it4].quad_sol;
                                                            set<poly_quad> set_quad_4 = set_quad_3;
                                                            set_quad_4.merge(set_4);

                                                            for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                set<poly_quad> set_quad_5 = set_quad_4;
                                                                set_quad_5.merge(set_5);

                                                                for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                    set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                    set<poly_quad> set_quad_6 = set_quad_5;
                                                                    set_quad_6.merge(set_6);

                                                                    for (uint32_t it7 = 0; it7<v7.size(); it7++) {
                                                                                    
                                                                        set<poly_quad> set_7 = v7[it7].quad_sol;
                                                                        set<poly_quad> set_quad_7 = set_quad_6;
                                                                        set_quad_7.merge(set_7);
                                                                        set_quad_7.erase(0);

                                                                        set<poly_quad> s_quad_temp = set_quad_7;
                                                                        uint32_t r = Rank(s_quad_temp);

                                                                        if (r<=MAX_RANK_DEG5_V2)   {
                                                                            implem temp;
                                                                            temp.quad_sol = set_quad_7;
                                                                            temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                            temp.op_sol.insert(temp.op_sol.end(),v7[it7].op_sol.begin(), v7[it7].op_sol.end());
                                                                            
                                                                            string f = "";
                                                                        
                                                                            if (v2.size()==1){
                                                                                if (v4.size()==1){
                                                                                    if (v5.size()==1){
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v5.size()==1){
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                            else {
                                                                                if (v4.size()==1){
                                                                                    if (v5.size()==1){
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v5.size()==1){
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v7.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16 ^ p17";
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }                                                                        
                                                                            
                                                                            temp.formula = f;
                                                                            res.push_back(temp);  
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else {
                                        if (n5!=0 && n5<=2){
                                            for (uint32_t jj5=0; jj5<size; jj5++){

                                                create_monomial_order(monomial_order, nb_elem, size, jj5);                                                   
                                                pair<poly,poly> pp5;
                                                pp5 = poly_div(p5.second, set_op_l_plus_lin[i5], monomial_order, size, nb_elem);

                                                uint32_t mm5 = pp5.second.algebraic_degree(nb_elem);
                                                uint32_t nn5 = pp5.first.algebraic_degree(nb_elem);

                                                if (mm5<=2 && nn5<=2){
                                                // y = ( (P2 * P3) + P4) * P1 + (P5 * P6) + P7
                                                    vector<implem> v1;
                    
                                                    implem imp;
                                                    imp.op_sol.push_back(P1);
                                                    imp.op_sol.push_back(P1_lin);
                                                    imp.quad_sol.insert(pq1);
                                                    imp.formula = "";
                                                    v1.push_back(imp);

                                                    vector<implem> v2;

                                                    int32_t indice2 = indiceP2;
                                                    if (indice2 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.quad_sol.insert(pq2);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                    else {
                                                        indice2 = indice_mP2;
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                    }

                                                    vector<implem> v3;

                                                    implem imp3;
                                                    imp3.op_sol.push_back(P3);
                                                    imp3.op_sol.push_back(P3_lin);
                                                    imp3.quad_sol.insert(pq3);
                                                    imp3.formula = "";
                                                    v3.push_back(imp3);

                                                    vector<implem> v4;

                                                    if (m3 == 0){
                                                        implem imp;
                                                        imp.op_sol.push_back(P4);
                                                        imp.quad_sol.insert(0);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                    else {
                                                        int32_t indice4 = indiceP4;
                                                        if (indice4 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P4);
                                                            imp.op_sol.push_back(P4_lin);
                                                            imp.quad_sol.insert(pq4);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                        else {
                                                            indice4 = indice_mP4;
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P4_lin);
                                                                imp.formula = "";
                                                                v4.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v5;

                                                    poly PP5;
                                                    POLY_ADD(p5.first, pp5.first, PP5);
                                                    poly P5 = truncate_lin(PP5, size);
                                                    poly P5_lin;
                                                    POLY_ADD(PP5, P5, P5_lin);
                                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                    if (indice5 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P5);
                                                        imp.op_sol.push_back(P5_lin);
                                                        imp.quad_sol.insert(pq5);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                    else {
                                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                        if (indice5 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P5_lin);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v6;

                                                    poly P6;
                                                    P6 = truncate_lin(set_op_l_plus_lin[i5], size);
                                                    poly P6_lin;
                                                    POLY_ADD(P6, set_op_l_plus_lin[i5], P6_lin);
                                                    poly_quad pq6 = poly_to_poly_quad(P6 ,size, nb_elem);
                                                    implem imp6;
                                                    imp6.op_sol.push_back(P6);
                                                    imp6.op_sol.push_back(P6_lin);
                                                    imp6.quad_sol.insert(pq6);
                                                    imp6.formula = "";
                                                    v6.push_back(imp6);

                                                    if (mm5 == 0){
                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);
                                                                                set_quad_6.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_6;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                                    implem temp;
                                                                                    temp.quad_sol = set_quad_6;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                    
                                                                                    string f = "";

                                                                                    if (v2.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14 ^ p15)";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                        
                                                                                    temp.formula = f;

                                                                                    res.push_back(temp);  
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    else {
                                                        vector<implem> v7;

                                                        poly_quad pq7 = poly_to_poly_quad(pp5.second ,size, nb_elem);
                                                        int32_t indice7 = find_set(pq7, set_op, 0, set_op_size);
                                                        if (indice7 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(pp5.second);
                                                            imp.quad_sol.insert(pq7);
                                                            imp.formula = "";
                                                            v7.push_back(imp);
                                                        }
                                                        else {
                                                            indice7 = find_map(pq7, map_xor, 0, map_xor_size);
                                                            if (indice7 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice7].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice7].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.formula = "";
                                                                    v7.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);

                                                                                for (uint32_t it7 = 0; it7<v7.size(); it7++) {
                                                                                    
                                                                                    set<poly_quad> set_7 = v7[it7].quad_sol;
                                                                                    set<poly_quad> set_quad_7 = set_quad_6;
                                                                                    set_quad_7.merge(set_7);
                                                                                    set_quad_7.erase(0);

                                                                                    set<poly_quad> s_quad_temp = set_quad_7;
                                                                                    uint32_t r = Rank(s_quad_temp);

                                                                                    if (r<=MAX_RANK_DEG5_V2)   {
                                                                                        implem temp;
                                                                                        temp.quad_sol = set_quad_7;
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v7[it7].op_sol.begin(), v7[it7].op_sol.end());
                                                                                        
                                                                                        string f = "";
                                                                        
                                                                                        if (v2.size()==1){
                                                                                            if (v4.size()==1){
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v4.size()==1){
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16 ^ p17";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }                                                                        
                                                                                        
                                                                                        temp.formula = f;
                                                                                        res.push_back(temp);  
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }        
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else {
                        if (n3 != 0 && n3 <= 2){
                            for (uint32_t jj3=0; jj3<size; jj3++){
                                if (jj3 == j3){
                                    continue;
                                }

                                create_monomial_order(monomial_order, nb_elem, size, jj3);                                                   
                                pair<poly,poly> pp3;
                                pp3 = poly_div(p3.second, set_op_l_plus_lin[i3], monomial_order, size, nb_elem);

                                uint32_t mm3 = pp3.second.algebraic_degree(nb_elem);
                                uint32_t nn3 = pp3.first.algebraic_degree(nb_elem);

                                if (mm3 <= 2 && nn3 <= 2){
                                //y5_4 = .quad * ((p3.first + pp3.first)* s[i3] + pp3.second) + reste
                                    poly PP2;
                                    POLY_ADD(p3.first, pp3.first, PP2);
                                    poly P2 = truncate_lin(PP2, size);
                                    poly P2_lin;
                                    POLY_ADD(PP2, P2, P2_lin);
                                    poly_quad pq2 = poly_to_poly_quad(P2 ,size, nb_elem);
                                    int32_t indiceP2 = find_set(pq2, set_op, 0, set_op_size);
                                    int32_t indice_mP2 = -1;
                                    if (indiceP2 == -1){
                                        indice_mP2 = find_map(pq2, map_xor, 0, map_xor_size);
                                        if (indice_mP2 == -1){
                                            continue;
                                        }
                                    }

                                    poly P3;
                                    P3 = truncate_lin(set_op_l_plus_lin[i3], size);
                                    poly P3_lin;
                                    POLY_ADD(P3, set_op_l_plus_lin[i3], P3_lin);
                                    poly_quad pq3 = poly_to_poly_quad(P3 ,size, nb_elem);

                                    poly_quad pq4;
                                    int32_t indiceP4 = -1; 
                                    int32_t indice_mP4 = -1;
                                    if (mm3 != 0){
                                        pq4 = poly_to_poly_quad(pp3.second ,size, nb_elem);
                                        indiceP4 = find_set(pq4, set_op, 0, set_op_size); 
                                        if (indiceP4 == -1){
                                            indice_mP4 = find_map(pq4, map_xor, 0, map_xor_size);
                                            if (indice_mP4 == -1){
                                                continue;
                                            }
                                        }
                                    }
                                    poly P4;
                                    poly m1;
                                    POLY_MUL(PP2, set_op_l_plus_lin[i3], m1, nb_elem);
                                    POLY_ADD(m1, quotient, P4);
                                    poly P4_lin;
                                    POLY_ADD(P4, pp3.second, P4_lin);

                                    bool test_reste_div_par_l = false ;
                                    //bool test_reste_div_par_l = true ;

                                    for (uint32_t i4=0; i4<l.size(); i4++){
                                        for (uint32_t j4=0; j4<size; j4++){

                                            if (mon_dom_in(l[i4],j4, size, nb_elem)){

                                                create_monomial_order(monomial_order, nb_elem, size, j4);                                                   
                                                pair<poly,poly> p4;
                                                p4 = poly_div(reste, l[i4], monomial_order, size, nb_elem);

                                                uint32_t m4 = p4.second.algebraic_degree(nb_elem);
                                                uint32_t n4 = p4.first.algebraic_degree(nb_elem);

                                                if (m4 <= 2 && n4 <= 2){
                                                test_reste_div_par_l = true;
                                                // y = ( (PP2 * P3) + P4) * P1 + (P5 * l[i4]) + P6

                                                    vector<implem> v1;
                                        
                                                    implem imp1;
                                                    imp1.op_sol.push_back(P1);
                                                    imp1.op_sol.push_back(P1_lin);
                                                    imp1.quad_sol.insert(pq1);
                                                    imp1.formula = "";
                                                    v1.push_back(imp1);

                                                    vector<implem> v2;

                                                    int32_t indice2 = indiceP2;
                                                    if (indice2 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.quad_sol.insert(pq2);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                    else {
                                                        indice2 = indice_mP2;
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                    }

                                                    vector<implem> v3;

                                                    implem imp3;
                                                    imp3.op_sol.push_back(P3);
                                                    imp3.op_sol.push_back(P3_lin);
                                                    imp3.quad_sol.insert(pq3);
                                                    imp3.formula = "";
                                                    v3.push_back(imp3);

                                                    vector<implem> v4;

                                                    if (mm3 == 0){
                                                        implem imp;
                                                        imp.op_sol.push_back(P4);
                                                        imp.op_sol.push_back(P4_lin);
                                                        imp.quad_sol.insert(0);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                    else {
                                                        int32_t indice4 = indiceP4;
                                                        if (indice4 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P4);
                                                            imp.op_sol.push_back(P4_lin);
                                                            imp.quad_sol.insert(pq4);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                        else {
                                                            indice4 = indice_mP4;
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P4_lin);
                                                                imp.formula = "";
                                                                v4.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v5;

                                                    poly P5;
                                                    P5 = truncate_lin(p4.first, size);
                                                    poly P5_lin;
                                                    POLY_ADD(P5, p4.first, P5_lin);
                                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                    if (indice5 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P5);
                                                        imp.op_sol.push_back(P5_lin);
                                                        imp.op_sol.push_back(l[i4]);
                                                        imp.quad_sol.insert(pq5);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                    else {
                                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                        if (indice5 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P5_lin);
                                                                imp.op_sol.push_back(l[i4]);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    if (m4 == 0){
                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);
                                                                            set_quad_5.erase(0);

                                                                            set<poly_quad> s_quad_temp = set_quad_5;
                                                                            uint32_t r = Rank(s_quad_temp);

                                                                            if (r<=MAX_RANK_DEG5_V2)   {
                                                                                implem temp;
                                                                                temp.quad_sol = set_quad_5;
                                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                
                                                                                string f = "";

                                                                                if (v2.size()==1){
                                                                                    if (v4.size()==1){
                                                                                        if (v5.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12)";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v5.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13)";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                else {
                                                                                    if (v4.size()==1){
                                                                                        if (v5.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13)";
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v5.size()==1){
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13)";
                                                                                        }
                                                                                        else {
                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14)";
                                                                                        }
                                                                                    }
                                                                                }
                                                                                    
                                                                                temp.formula = f;

                                                                                res.push_back(temp);  
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    else {

                                                        vector<implem> v6;

                                                        poly_quad pq6 = poly_to_poly_quad(p4.second ,size, nb_elem);
                                                        int32_t indice6 = find_set(pq6, set_op, 0, set_op_size);
                                                        if (indice6 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(p4.second);
                                                            imp.quad_sol.insert(pq6);
                                                            imp.formula = "";
                                                            v6.push_back(imp);
                                                        }
                                                        else {
                                                            indice6 = find_map(pq6, map_xor, 0, map_xor_size);
                                                            if (indice6 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice6].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice6].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.formula = "";
                                                                    v6.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);
                                                                                set_quad_6.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_6;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                                    implem temp;
                                                                                    temp.quad_sol = set_quad_6;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                    
                                                                                    string f = "";
                                                                        
                                                                                    if (v2.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11) ) ^ p12 ^ p13";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12) ) ^ p13";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12) ) ^ p13 ^ p14";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13) ) ^ p14";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13) ) ^ p14 ^ p15";
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v6.size()==1){
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14) ) ^ p15";
                                                                                                }
                                                                                                else {
                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14) ) ^ p15 ^ p16";
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }                                                                        
                                                                                    
                                                                                    temp.formula = f;

                                                                                    res.push_back(temp);  
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    if (!test_reste_div_par_l){
                                        for (uint32_t i5=0; i5<set_op_l_plus_lin.size(); i5++){
                                            for (uint32_t j5=0; j5<size; j5++){

                                                create_monomial_order(monomial_order, nb_elem, size, j5);                                                   
                                                pair<poly,poly> p5;
                                                p5 = poly_div(reste, set_op_l_plus_lin[i5], monomial_order, size, nb_elem);

                                                uint32_t m5 = p5.second.algebraic_degree(nb_elem);
                                                uint32_t n5 = p5.first.algebraic_degree(nb_elem);

                                                if (m5 <= 2 && n5 <= 2){
                                                // y = ( (PP2 * P3) + P4) * P1 + (P5 * P6) + P7
                                                    vector<implem> v1;
                                      
                                                    implem imp1;
                                                    imp1.op_sol.push_back(P1);
                                                    imp1.op_sol.push_back(P1_lin);
                                                    imp1.quad_sol.insert(pq1);
                                                    imp1.formula = "";
                                                    v1.push_back(imp1);

                                                    vector<implem> v2;

                                                    int32_t indice2 = indiceP2;
                                                    if (indice2 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P2);
                                                        imp.op_sol.push_back(P2_lin);
                                                        imp.quad_sol.insert(pq2);
                                                        imp.formula = "";
                                                        v2.push_back(imp);
                                                    }
                                                    else {
                                                        indice2 = indice_mP2;
                                                        for (uint32_t m=0; m<10; m++)   {
                                                            implem imp;
                                                            poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                            poly op1(op_q1,size);
                                                            imp.op_sol.push_back(op1);
                                                            imp.quad_sol.insert(op_q1);
                                                            poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                            poly op2(op_q2,size);
                                                            imp.op_sol.push_back(op2);
                                                            imp.quad_sol.insert(op_q2);
                                                            imp.op_sol.push_back(P2_lin);
                                                            imp.formula = "";
                                                            v2.push_back(imp);
                                                        }
                                                    }

                                                    vector<implem> v3;

                                                    implem imp3;
                                                    imp3.op_sol.push_back(P3);
                                                    imp3.op_sol.push_back(P3_lin);
                                                    imp3.quad_sol.insert(pq3);
                                                    imp3.formula = "";
                                                    v3.push_back(imp3);

                                                    vector<implem> v4;

                                                    if (mm3 == 0){
                                                        implem imp;
                                                        imp.op_sol.push_back(P4);
                                                        imp.op_sol.push_back(P4_lin);
                                                        imp.quad_sol.insert(0);
                                                        imp.formula = "";
                                                        v4.push_back(imp);
                                                    }
                                                    else {
                                                        int32_t indice4 = indiceP4;
                                                        if (indice4 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(P4);
                                                            imp.op_sol.push_back(P4_lin);
                                                            imp.quad_sol.insert(pq4);
                                                            imp.formula = "";
                                                            v4.push_back(imp);
                                                        }
                                                        else {
                                                            indice4 = indice_mP4;
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P4_lin);
                                                                imp.formula = "";
                                                                v4.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v5;

                                                    poly P5;
                                                    P5 = truncate_lin(p5.first, size);
                                                    poly P5_lin;
                                                    POLY_ADD(P5, p5.first, P5_lin);
                                                    poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                                    int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                    if (indice5 != -1)   {
                                                        implem imp;
                                                        imp.op_sol.push_back(P5);
                                                        imp.op_sol.push_back(P5_lin);
                                                        imp.quad_sol.insert(pq5);
                                                        imp.formula = "";
                                                        v5.push_back(imp);
                                                    }
                                                    else {
                                                        indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                        if (indice5 == -1)  {
                                                            continue;
                                                        }
                                                        else {
                                                            for (uint32_t m=0; m<10; m++)   {
                                                                implem imp;
                                                                poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                poly op1(op_q1,size);
                                                                imp.op_sol.push_back(op1);
                                                                imp.quad_sol.insert(op_q1);
                                                                poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                poly op2(op_q2,size);
                                                                imp.op_sol.push_back(op2);
                                                                imp.quad_sol.insert(op_q2);
                                                                imp.op_sol.push_back(P5_lin);
                                                                imp.formula = "";
                                                                v5.push_back(imp);
                                                            }
                                                        }
                                                    }

                                                    vector<implem> v6;

                                                    poly P6;
                                                    P6 = truncate_lin(set_op_l_plus_lin[i5], size);
                                                    poly P6_lin;
                                                    POLY_ADD(P6, set_op_l_plus_lin[i5], P6_lin);
                                                    poly_quad pq6 = poly_to_poly_quad(P6 ,size, nb_elem);
                                                    implem imp6;
                                                    imp6.op_sol.push_back(P6);
                                                    imp6.op_sol.push_back(P6_lin);
                                                    imp6.quad_sol.insert(pq6);
                                                    imp6.formula = "";
                                                    v6.push_back(imp6);

                                                    if (m5 == 0){
                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);
                                                                                set_quad_6.erase(0);

                                                                                set<poly_quad> s_quad_temp = set_quad_6;
                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                                    implem temp;
                                                                                    temp.quad_sol = set_quad_6;
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                    temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                    
                                                                                    string f = "";

                                                                                    if (v2.size()==1){
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                    else {
                                                                                        if (v4.size()==1){
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v5.size()==1){
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13 ^ p14)";
                                                                                            }
                                                                                            else {
                                                                                                f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14 ^ p15)";
                                                                                            }
                                                                                        }
                                                                                    }
                                                                            
                                                                                    temp.formula = f;
                                                                                    res.push_back(temp);  
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                    else {

                                                        vector<implem> v7;

                                                        poly_quad pq7 = poly_to_poly_quad(p5.second ,size, nb_elem);
                                                        int32_t indice7 = find_set(pq7, set_op, 0, set_op_size);
                                                        if (indice7 != -1)   {
                                                            implem imp;
                                                            imp.op_sol.push_back(p5.second);
                                                            imp.quad_sol.insert(pq7);
                                                            imp.formula = "";
                                                            v7.push_back(imp);
                                                        }
                                                        else {
                                                            indice7 = find_map(pq7, map_xor, 0, map_xor_size);
                                                            if (indice7 == -1)  {
                                                                continue;
                                                            }
                                                            else {
                                                                for (uint32_t m=0; m<10; m++)   {
                                                                    implem imp;
                                                                    poly_quad op_q1 = map_xor[indice7].second[m][0];
                                                                    poly op1(op_q1,size);
                                                                    imp.op_sol.push_back(op1);
                                                                    imp.quad_sol.insert(op_q1);
                                                                    poly_quad op_q2 = map_xor[indice7].second[m][1];
                                                                    poly op2(op_q2,size);
                                                                    imp.op_sol.push_back(op2);
                                                                    imp.quad_sol.insert(op_q2);
                                                                    imp.formula = "";
                                                                    v7.push_back(imp);
                                                                }
                                                            }
                                                        }

                                                        for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                            set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                            for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                set<poly_quad> set_quad_2 = set_quad_1;
                                                                set_quad_2.merge(set_2);

                                                                for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                    set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                    set<poly_quad> set_quad_3 = set_quad_2;
                                                                    set_quad_3.merge(set_3);

                                                                    for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                        set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                        set<poly_quad> set_quad_4 = set_quad_3;
                                                                        set_quad_4.merge(set_4);

                                                                        for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                            set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                            set<poly_quad> set_quad_5 = set_quad_4;
                                                                            set_quad_5.merge(set_5);

                                                                            for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                set<poly_quad> set_quad_6 = set_quad_5;
                                                                                set_quad_6.merge(set_6);

                                                                                for (uint32_t it7 = 0; it7<v7.size(); it7++) {
                                                                                    
                                                                                    set<poly_quad> set_7 = v7[it7].quad_sol;
                                                                                    set<poly_quad> set_quad_7 = set_quad_6;
                                                                                    set_quad_7.merge(set_7);
                                                                                    set_quad_7.erase(0);

                                                                                    set<poly_quad> s_quad_temp = set_quad_7;
                                                                                    uint32_t r = Rank(s_quad_temp);

                                                                                    if (r<=MAX_RANK_DEG5_V2)   {
                                                                                        implem temp;
                                                                                        temp.quad_sol = set_quad_7;
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                        temp.op_sol.insert(temp.op_sol.end(),v7[it7].op_sol.begin(), v7[it7].op_sol.end());
                                                                                        
                                                                                        string f = "";
                                                                        
                                                                                        if (v2.size()==1){
                                                                                            if (v4.size()==1){
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                        else {
                                                                                            if (v4.size()==1){
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                            else {
                                                                                                if (v5.size()==1){
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v7.size()==1){
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16";
                                                                                                    }
                                                                                                    else {
                                                                                                        f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16 ^ p17";
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }                                                                        
                                                                                        
                                                                                        temp.formula = f;
                                                                                        res.push_back(temp);  
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                else {
                                                    if (n5 != 0 && n5 <= 2){
                                                        for (uint32_t jj5=0; jj5<size; jj5++){

                                                            if (jj5 == j5){
                                                                continue;
                                                            }

                                                            create_monomial_order(monomial_order, nb_elem, size, jj5);                                                   
                                                            pair<poly,poly> pp5;
                                                            pp5 = poly_div(p5.second, set_op_l_plus_lin[i5], monomial_order, size, nb_elem);

                                                            uint32_t mm5 = pp5.second.algebraic_degree(nb_elem);
                                                            uint32_t nn5 = pp5.first.algebraic_degree(nb_elem);

                                                            if (mm5<=2 && nn5<=2){
                                                            // y = ( (PP2 * P3) + P4) * P1 + (PP5 * P6) + P7

                                                                vector<implem> v1;
                                      
                                                                implem imp1;
                                                                imp1.op_sol.push_back(P1);
                                                                imp1.op_sol.push_back(P1_lin);
                                                                imp1.quad_sol.insert(pq1);
                                                                imp1.formula = "";
                                                                v1.push_back(imp1);

                                                                vector<implem> v2;

                                                                int32_t indice2 = indiceP2;
                                                                if (indice2 != -1)   {
                                                                    implem imp;
                                                                    imp.op_sol.push_back(P2);
                                                                    imp.op_sol.push_back(P2_lin);
                                                                    imp.quad_sol.insert(pq2);
                                                                    imp.formula = "";
                                                                    v2.push_back(imp);
                                                                }
                                                                else {
                                                                    indice2 = indice_mP2;
                                                                    for (uint32_t m=0; m<10; m++)   {
                                                                        implem imp;
                                                                        poly_quad op_q1 = map_xor[indice2].second[m][0];
                                                                        poly op1(op_q1,size);
                                                                        imp.op_sol.push_back(op1);
                                                                        imp.quad_sol.insert(op_q1);
                                                                        poly_quad op_q2 = map_xor[indice2].second[m][1];
                                                                        poly op2(op_q2,size);
                                                                        imp.op_sol.push_back(op2);
                                                                        imp.quad_sol.insert(op_q2);
                                                                        imp.op_sol.push_back(P2_lin);
                                                                        imp.formula = "";
                                                                        v2.push_back(imp);
                                                                    }
                                                                }

                                                                vector<implem> v3;

                                                                implem imp3;
                                                                imp3.op_sol.push_back(P3);
                                                                imp3.op_sol.push_back(P3_lin);
                                                                imp3.quad_sol.insert(pq3);
                                                                imp3.formula = "";
                                                                v3.push_back(imp3);

                                                                vector<implem> v4;

                                                                if (mm3 == 0){
                                                                    implem imp;
                                                                    imp.op_sol.push_back(P4);
                                                                    imp.op_sol.push_back(P4_lin);
                                                                    imp.quad_sol.insert(0);
                                                                    imp.formula = "";
                                                                    v4.push_back(imp);
                                                                }
                                                                else {

                                                                    int32_t indice4 = indiceP4;
                                                                    if (indice4 != -1)   {
                                                                        implem imp;
                                                                        imp.op_sol.push_back(P4);
                                                                        imp.op_sol.push_back(P4_lin);
                                                                        imp.quad_sol.insert(pq4);
                                                                        imp.formula = "";
                                                                        v4.push_back(imp);
                                                                    }
                                                                    else {
                                                                        indice4 = indice_mP4;
                                                                        for (uint32_t m=0; m<10; m++)   {
                                                                            implem imp;
                                                                            poly_quad op_q1 = map_xor[indice4].second[m][0];
                                                                            poly op1(op_q1,size);
                                                                            imp.op_sol.push_back(op1);
                                                                            imp.quad_sol.insert(op_q1);
                                                                            poly_quad op_q2 = map_xor[indice4].second[m][1];
                                                                            poly op2(op_q2,size);
                                                                            imp.op_sol.push_back(op2);
                                                                            imp.quad_sol.insert(op_q2);
                                                                            imp.op_sol.push_back(P4_lin);
                                                                            imp.formula = "";
                                                                            v4.push_back(imp);
                                                                        }
                                                                    }
                                                                }

                                                                vector<implem> v5;

                                                                poly PP5;
                                                                POLY_ADD(p5.first, pp5.first ,PP5);
                                                                poly P5 = truncate_lin(PP5, size);
                                                                poly P5_lin;
                                                                POLY_ADD(PP5, P5, P5_lin);
                                                                poly_quad pq5 = poly_to_poly_quad(P5 ,size, nb_elem);
                                                                int32_t indice5 = find_set(pq5, set_op, 0, set_op_size);
                                                                if (indice5 != -1)   {
                                                                    implem imp;
                                                                    imp.op_sol.push_back(P5);
                                                                    imp.op_sol.push_back(P5_lin);
                                                                    imp.quad_sol.insert(pq5);
                                                                    imp.formula = "";
                                                                    v5.push_back(imp);
                                                                }
                                                                else {
                                                                    indice5 = find_map(pq5, map_xor, 0, map_xor_size);
                                                                    if (indice5 == -1)  {
                                                                        continue;
                                                                    }
                                                                    else {
                                                                        for (uint32_t m=0; m<10; m++)   {
                                                                            implem imp;
                                                                            poly_quad op_q1 = map_xor[indice5].second[m][0];
                                                                            poly op1(op_q1,size);
                                                                            imp.op_sol.push_back(op1);
                                                                            imp.quad_sol.insert(op_q1);
                                                                            poly_quad op_q2 = map_xor[indice5].second[m][1];
                                                                            poly op2(op_q2,size);
                                                                            imp.op_sol.push_back(op2);
                                                                            imp.quad_sol.insert(op_q2);
                                                                            imp.op_sol.push_back(P5_lin);
                                                                            imp.formula = "";
                                                                            v5.push_back(imp);
                                                                        }
                                                                    }
                                                                }

                                                                vector<implem> v6;

                                                                poly P6;
                                                                P6 = truncate_lin(set_op_l_plus_lin[i5], size);
                                                                poly P6_lin;
                                                                POLY_ADD(P6, set_op_l_plus_lin[i5], P6_lin);
                                                                poly_quad pq6 = poly_to_poly_quad(P6 ,size, nb_elem);
                                                                implem imp6;
                                                                imp6.op_sol.push_back(P6);
                                                                imp6.op_sol.push_back(P6_lin);
                                                                imp6.quad_sol.insert(pq6);
                                                                imp6.formula = "";
                                                                v6.push_back(imp6);

                                                                if (mm5 == 0){
                                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                            set<poly_quad> set_quad_2 = set_quad_1;
                                                                            set_quad_2.merge(set_2);

                                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                                                set_quad_3.merge(set_3);

                                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                                                    set_quad_4.merge(set_4);

                                                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                                                        set_quad_5.merge(set_5);

                                                                                        for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                            set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                            set<poly_quad> set_quad_6 = set_quad_5;
                                                                                            set_quad_6.merge(set_6);
                                                                                            set_quad_6.erase(0);

                                                                                            set<poly_quad> s_quad_temp = set_quad_6;
                                                                                            uint32_t r = Rank(s_quad_temp);

                                                                                            if (r<=MAX_RANK_DEG5_V2)   {
                                                                                                implem temp;
                                                                                                temp.quad_sol = set_quad_6;
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                                temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                                string f = "";

                                                                                                if (v2.size()==1){
                                                                                                    if (v4.size()==1){
                                                                                                        if (v5.size()==1){
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10) * (p11 ^ p12)";
                                                                                                        }
                                                                                                        else {
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ (p9 ^ p10 ^ p11) * (p12 ^ p13)";
                                                                                                        }
                                                                                                    }
                                                                                                    else {
                                                                                                        if (v5.size()==1){
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                                        }
                                                                                                        else {
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                                else {
                                                                                                    if (v4.size()==1){
                                                                                                        if (v5.size()==1){
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11) * (p12 ^ p13)";
                                                                                                        }
                                                                                                        else {
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ (p10 ^ p11 ^ p12) * (p13 ^ p14)";
                                                                                                        }
                                                                                                    }
                                                                                                    else {
                                                                                                        if (v5.size()==1){
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12) * (p13 ^ p14)";
                                                                                                        }
                                                                                                        else {
                                                                                                            f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ (p11 ^ p12 ^ p13) * (p14 ^ p15)";
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                                    
                                                                                                temp.formula = f;

                                                                                                res.push_back(temp);  
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                                else {
                                                                    vector<implem> v7;

                                                                    poly_quad pq7 = poly_to_poly_quad(pp5.second ,size, nb_elem);
                                                                    int32_t indice7 = find_set(pq7, set_op, 0, set_op_size);
                                                                    if (indice7 != -1)   {
                                                                        implem imp;
                                                                        imp.op_sol.push_back(pp5.second);
                                                                        imp.quad_sol.insert(pq7);
                                                                        imp.formula = "";
                                                                        v7.push_back(imp);
                                                                    }
                                                                    else {
                                                                        indice7 = find_map(pq7, map_xor, 0, map_xor_size);
                                                                        if (indice7 == -1)  {
                                                                            continue;
                                                                        }
                                                                        else {
                                                                            for (uint32_t m=0; m<10; m++)   {
                                                                                implem imp;
                                                                                poly_quad op_q1 = map_xor[indice7].second[m][0];
                                                                                poly op1(op_q1,size);
                                                                                imp.op_sol.push_back(op1);
                                                                                imp.quad_sol.insert(op_q1);
                                                                                poly_quad op_q2 = map_xor[indice7].second[m][1];
                                                                                poly op2(op_q2,size);
                                                                                imp.op_sol.push_back(op2);
                                                                                imp.quad_sol.insert(op_q2);
                                                                                imp.formula = "";
                                                                                v7.push_back(imp);
                                                                            }
                                                                        }
                                                                    }

                                                                    for (uint32_t it1 = 0; it1<v1.size(); it1++) {

                                                                        set<poly_quad> set_quad_1 = v1[it1].quad_sol;

                                                                        for (uint32_t it2 = 0; it2<v2.size(); it2++) {

                                                                            set<poly_quad> set_2 = v2[it2].quad_sol;
                                                                            set<poly_quad> set_quad_2 = set_quad_1;
                                                                            set_quad_2.merge(set_2);

                                                                            for (uint32_t it3 = 0; it3<v3.size(); it3++) {

                                                                                set<poly_quad> set_3 = v3[it3].quad_sol;
                                                                                set<poly_quad> set_quad_3 = set_quad_2;
                                                                                set_quad_3.merge(set_3);

                                                                                for (uint32_t it4 = 0; it4<v4.size(); it4++) {

                                                                                    set<poly_quad> set_4 = v4[it4].quad_sol;
                                                                                    set<poly_quad> set_quad_4 = set_quad_3;
                                                                                    set_quad_4.merge(set_4);

                                                                                    for (uint32_t it5 = 0; it5<v5.size(); it5++) {

                                                                                        set<poly_quad> set_5 = v5[it5].quad_sol;
                                                                                        set<poly_quad> set_quad_5 = set_quad_4;
                                                                                        set_quad_5.merge(set_5);

                                                                                        for (uint32_t it6 = 0; it6<v6.size(); it6++) {

                                                                                            set<poly_quad> set_6 = v6[it6].quad_sol;
                                                                                            set<poly_quad> set_quad_6 = set_quad_5;
                                                                                            set_quad_6.merge(set_6);

                                                                                            for (uint32_t it7 = 0; it7<v7.size(); it7++) {
                                                                                                
                                                                                                set<poly_quad> set_7 = v7[it7].quad_sol;
                                                                                                set<poly_quad> set_quad_7 = set_quad_6;
                                                                                                set_quad_7.merge(set_7);
                                                                                                set_quad_7.erase(0);

                                                                                                set<poly_quad> s_quad_temp = set_quad_7;
                                                                                                uint32_t r = Rank(s_quad_temp);

                                                                                                if (r<=MAX_RANK_DEG5_V2)   {
                                                                                                    implem temp;
                                                                                                    temp.quad_sol = set_quad_7;
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v1[it1].op_sol.begin(), v1[it1].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v2[it2].op_sol.begin(), v2[it2].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v3[it3].op_sol.begin(), v3[it3].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v4[it4].op_sol.begin(), v4[it4].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v5[it5].op_sol.begin(), v5[it5].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v6[it6].op_sol.begin(), v6[it6].op_sol.end());
                                                                                                    temp.op_sol.insert(temp.op_sol.end(),v7[it7].op_sol.begin(), v7[it7].op_sol.end());
                                                                                                    
                                                                                                    string f = "";
                                                                        
                                                                                                    if (v2.size()==1){
                                                                                                        if (v4.size()==1){
                                                                                                            if (v5.size()==1){
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10) * (p11 ^ p12) ) ^ p13 ^ p14";
                                                                                                                }
                                                                                                            }
                                                                                                            else {
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8) ^ ( (p9 ^ p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                        else {
                                                                                                            if (v5.size()==1){
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                                }
                                                                                                            }
                                                                                                            else {
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4) * (p5 ^ p6) ^ p7 ^ p8 ^ p9 ) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                    else {
                                                                                                        if (v4.size()==1){
                                                                                                            if (v5.size()==1){
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11) * (p12 ^ p13) ) ^ p14 ^ p15";
                                                                                                                }
                                                                                                            }
                                                                                                            else {
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9) ^ ( (p10 ^ p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                        else {
                                                                                                            if (v5.size()==1){
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12) * (p13 ^ p14) ) ^ p15 ^ p16";
                                                                                                                }
                                                                                                            }
                                                                                                            else {
                                                                                                                if (v7.size()==1){
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16";
                                                                                                                }
                                                                                                                else {
                                                                                                                    f = "(p1 ^ p2) * ( (p3 ^ p4 ^ p5) * (p6 ^ p7) ^ p8 ^ p9 ^ p10) ^ ( (p11 ^ p12 ^ p13) * (p14 ^ p15) ) ^ p16 ^ p17";
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }                                                                        
                                                                                                    
                                                                                                    temp.formula = f;
                                                                                                    res.push_back(temp);  
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }  
        }
    }

    res = remove_duplicates(res);
    if (verbose){
        cout<<" Thread : "<<omp_get_thread_num()<<" res_size v2 = "<<res.size()<<endl;
    }

    return res;
}

decomposition::decomposition(uint32_t size, uint32_t deg) {
    if (deg == 2){
        ftab[0] = &add_to_op_selec_deg_2_size_4_5_6_v1;
        ftab[1] = &add_to_op_selec_dft;
    }
    if (deg == 3){
        ftab[0] = &add_to_op_selec_deg_3_v1; 
        ftab[1] = &add_to_op_selec_deg_3_v2; 
    }
    if (deg == 4){
        ftab[0] = &add_to_op_selec_deg_4_v1;
        ftab[1] = &add_to_op_selec_deg_4_v2;
    }
    if ((deg == 5) && (size == 6)){
        ftab[0] = &add_to_op_selec_deg_5_v1;
        ftab[1] = &add_to_op_selec_deg_5_v2;
    }
}

vector<implem> decomposition::add_to_op_selec(poly y, vector<poly> l, vector<poly> set_op_l_plus_lin, uint32_t size, uint32_t deg, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size, uint32_t * test)   {
    vector<implem> vect_res;
    for (uint32_t i=0; i<*test+1; i++){
        vect_res = (*(ftab[i]))(y, l, set_op_l_plus_lin, size, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size);
        if (!vect_res.empty()){
            *test = i;
            return vect_res;
        }
    }
    return vect_res;
}

uint32_t nb_mult(uint32_t deg, uint32_t test){
    if (deg == 2){
        return 0;
    }
    if (deg ==3) {
        if(test == 0){
            return 1;
        }
        if(test == 1){
            return 2;
        }
    }
    if (deg == 4) {
        if(test == 0){
            return 1;
        }
        if(test == 1){
            return 2;
        }
    }
    if (deg == 5) {
        if(test == 0){
            return 2;
        }
        if (test == 1){
            return 3;
        }
    }
    return 100;
}

uint32_t find_size(vector<vector<implem>> v, uint32_t size){
    uint32_t iterator = v.size();
    for (uint32_t i=0; i<v.size(); i++){
        if (v[i].size() >= size){
            return i;
        }
    }
    return iterator;
}

void insert_to_size_vect (vector<implem> vect_to_add, vector<vector<implem>> * V, vector<uint32_t> T_to_add, vector<vector<uint32_t>> * T ){
    if ( (*V).empty()){
        (*V).push_back(vect_to_add);
        (*T).push_back(T_to_add);
    }
    else {
        uint32_t it = find_size(*V, vect_to_add.size());
        (*V).insert((*V).begin() + it, vect_to_add);
        (*T).insert((*T).begin() + it, T_to_add);
    }
}

vector<vector<implem>> create_sets(uint32_t * a, uint32_t * test, vector<vector<uint32_t>> * T, vector<poly> y, vector<poly> l, vector<poly> l2, uint32_t size_in, uint32_t size_out, uint32_t nb_elem, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size ){
    vector<vector<implem>> imp_p;
    vector<vector<implem>> * imp;
    imp = &imp_p;

    #pragma omp parallel sections num_threads(3)
    {
    
    #pragma omp section
    {

    //y0 :
        vector<uint32_t> T0;
        *test=TAB_SIZE -1;
        uint32_t deg = y[0].algebraic_degree(nb_elem);
        decomposition dec(size_in,deg);
        vector<implem> imp_0 = dec.add_to_op_selec(y[0], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
        *a += nb_mult(deg, *test);
        T0.push_back(imp_0.size());
    //

    cout<<"y0 done : "<<imp_0.size()<<endl;
    
    //y1 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T1;
        deg = y[1].algebraic_degree(nb_elem);
        decomposition dec1(size_in,deg);

        vector<implem> imp_1 = dec1.add_to_op_selec(y[1], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T1.push_back(imp_1.size());

        poly y0py1;
        POLY_ADD(y[0],y[1],y0py1);
        uint32_t deg2 = y0py1.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {
            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }

            int t1 = *test;
            vector<implem> imp_0_1 = dec1.add_to_op_selec(y0py1, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                // Concatène les vecteurs pour obtenir y1 et y0 + y1
                imp_1.insert(imp_1.end(),imp_0_1.begin(), imp_0_1.end());
                T1.push_back(imp_0_1.size());
            }
            else {
                deg = deg2;
                imp_1 = imp_0_1;
                T1[1] = 0;
            }
        }
        *a += nb_mult(deg, *test);
    //
    
    cout<<"y1 done : "<<imp_1.size()<<endl;

    //y2 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T2;
        deg = y[2].algebraic_degree(nb_elem);
        decomposition dec2(size_in,deg);

        vector<implem> imp_2 = dec2.add_to_op_selec(y[2], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T2.push_back(imp_2.size());

        poly y0py2 ;
        POLY_ADD(y[0],y[2],y0py2);
        deg2 = y0py2.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2);

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }

            int t1 = *test;
            vector<implem> imp_0_2 = dec.add_to_op_selec(y0py2, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_2.insert(imp_2.end(),imp_0_2.begin(), imp_0_2.end());
                T2.push_back(imp_2.size());
            }
            else { //Si on a trouvé un Xor de bits de sortie de degré inférieur au bit initial:
                deg = deg2; //On met à jour la borne pour les suivants
                imp_2 = imp_0_2; //On supprime les ensembles précédents
                for (uint32_t i=0; i<T2.size(); i++){ //On remet toutes les valeurs de T2 à 0.
                    T2[i] = 0;
                }
                T2.push_back(imp_2.size());
            }
        }
        else {
            T2.push_back(0); 
        }        

        poly y1py2 ;
        POLY_ADD(y[1],y[2],y1py2);
        deg2 = y1py2.algebraic_degree(nb_elem);

        if(deg2 <= deg){
        
            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2 = dec.add_to_op_selec(y1py2, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_2.insert(imp_2.end(),imp_1_2.begin(), imp_1_2.end());

                T2.push_back(imp_2.size());
            }
            else { 
                deg = deg2; 
                imp_2 = imp_1_2; 
                for (uint32_t i=0; i<T2.size(); i++){
                    T2[i] = 0;
                }
                T2.push_back(imp_2.size());
            }
        }
        else {
            T2.push_back(0);
        }

        poly y0py1py2 ;
        POLY_ADD(y[0],y1py2,y0py1py2);
        deg2 = y0py1py2.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2 = dec.add_to_op_selec(y0py1py2, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test); 
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_2.insert(imp_2.end(),imp_0_1_2.begin(), imp_0_1_2.end());
                T2.push_back(imp_2.size());
            }
            else { 
                deg = deg2; 
                imp_2 = imp_0_1_2; 
                for (uint32_t i=0; i<T2.size(); i++){
                    T2[i] = 0;
                }
                T2.push_back(imp_2.size());
            }
        }
        else {
            T2.push_back(0);
        }

        *a += nb_mult(deg, *test);
    //

    cout<<"y2 done : "<<imp_2.size()<<endl;
    
    //y3 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T3;
        deg = y[3].algebraic_degree(nb_elem);
        decomposition dec3(size_in,deg);

        vector<implem> imp_3 = dec3.add_to_op_selec(y[3], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T3.push_back(imp_3.size());

        poly y0py3 ;
        POLY_ADD(y[0],y[3],y0py3);
        deg2 = y0py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_3 = dec.add_to_op_selec(y0py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_0_3.begin(), imp_0_3.end());    
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_0_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y1py3 ;
        POLY_ADD(y[1],y[3],y1py3);
        deg2 = y1py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){
            
            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_3 = dec.add_to_op_selec(y1py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;
            
            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_1_3.begin(), imp_1_3.end());
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_1_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y2py3 ;
        POLY_ADD(y[2],y[3],y2py3);
        deg2 = y2py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){
            
            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_3 = dec.add_to_op_selec(y2py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_2_3.begin(), imp_2_3.end());
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_2_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y0py1py3 ;
        POLY_ADD(y[0],y1py3,y0py1py3);
        deg2 = y0py1py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_3 = dec.add_to_op_selec(y0py1py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_0_1_3.begin(), imp_0_1_3.end());
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_0_1_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y0py2py3 ;
        POLY_ADD(y[0],y2py3,y0py2py3);
        deg2 = y0py2py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_3 = dec.add_to_op_selec(y0py2py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_0_2_3.begin(), imp_0_2_3.end());
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_0_2_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y1py2py3 ;
        POLY_ADD(y[1],y2py3,y1py2py3);
        deg2 = y1py2py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_3 = dec.add_to_op_selec(y1py2py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_1_2_3.begin(), imp_1_2_3.end());
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_1_2_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }

        poly y0py1py2py3 ;
        POLY_ADD(y[0],y1py2py3,y0py1py2py3);
        deg2 = y0py1py2py3.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }

            int t1 = *test;
            vector<implem> imp_0_1_2_3 = dec.add_to_op_selec(y0py1py2py3, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_3.insert(imp_3.end(),imp_0_1_2_3.begin(), imp_0_1_2_3.end()); 
                T3.push_back(imp_3.size());
            }
            else { 
                deg = deg2; 
                imp_3 = imp_0_1_2_3; 
                for (uint32_t i=0; i<T3.size(); i++){
                    T3[i] = 0;
                }
                T3.push_back(imp_3.size());
            }
        }
        else {
            T3.push_back(0);
        }
        *a += nb_mult(deg, *test);
    //

    cout<<"y3 done : "<<imp_3.size()<<endl;

    #pragma omp critical
    {
    insert_to_size_vect(imp_0, imp, T0, T);
    insert_to_size_vect(imp_1, imp, T1, T);
    insert_to_size_vect(imp_2, imp, T2, T);
    insert_to_size_vect(imp_3, imp, T3, T);
    }

    }

    #pragma omp section
    {
   
    if (size_out > 4){
    //y4 :   
        *test=TAB_SIZE -1;
        vector<uint32_t> T4;
        uint32_t deg = y[4].algebraic_degree(nb_elem);
        decomposition dec4(size_in,deg);

        vector<implem> imp_4 = dec4.add_to_op_selec(y[4], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T4.push_back(imp_4.size());

        poly y0py4 ;
        POLY_ADD(y[0],y[4],y0py4);
        uint32_t deg2 = y0py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_4 = dec.add_to_op_selec(y0py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_4.begin(), imp_0_4.end()); 
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }

        }
        else {
            T4.push_back(0);
        }         

        poly y1py4 ;
        POLY_ADD(y[1],y[4],y1py4);
        deg2 = y1py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_4 = dec.add_to_op_selec(y1py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_1_4.begin(), imp_1_4.end()); 
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_1_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        } 

        poly y2py4 ;
        POLY_ADD(y[2],y[4],y2py4);
        deg2 = y2py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_4 = dec.add_to_op_selec(y2py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_2_4.begin(), imp_2_4.end());  
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_2_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        } 

        poly y3py4 ;
        POLY_ADD(y[3],y[4],y3py4);
        deg2 = y3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_3_4 = dec.add_to_op_selec(y3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_3_4.begin(), imp_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py1py4 ;
        POLY_ADD(y[0],y1py4,y0py1py4);
        deg2 = y0py1py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_4 = dec.add_to_op_selec(y0py1py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_1_4.begin(), imp_0_1_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_1_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py2py4 ;
        POLY_ADD(y[0],y2py4,y0py2py4);
        deg2 = y0py2py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_4 = dec.add_to_op_selec(y0py2py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_2_4.begin(), imp_0_2_4.end()); 
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_2_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py3py4 ;
        POLY_ADD(y[0],y3py4,y0py3py4);
        deg2 = y0py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_3_4 = dec.add_to_op_selec(y0py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_3_4.begin(), imp_0_3_4.end()); 
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y1py2py4 ;
        POLY_ADD(y[1],y2py4,y1py2py4);
        deg2 = y1py2py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){
        
            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_4 = dec.add_to_op_selec(y1py2py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_1_2_4.begin(), imp_1_2_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_1_2_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y1py3py4 ;
        POLY_ADD(y[1],y3py4,y1py3py4);
        deg2 = y1py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_3_4 = dec.add_to_op_selec(y1py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_1_3_4.begin(), imp_1_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_1_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y2py3py4 ;
        POLY_ADD(y[2],y3py4,y2py3py4);
        deg2 = y2py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_3_4 = dec.add_to_op_selec(y2py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_2_3_4.begin(), imp_2_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_2_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py1py2py4 ;
        POLY_ADD(y[0],y1py2py4,y0py1py2py4);
        deg2 = y0py1py2py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2_4 = dec.add_to_op_selec(y0py1py2py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_1_2_4.begin(), imp_0_1_2_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_1_2_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py1py3py4 ;
        POLY_ADD(y[0],y1py3py4,y0py1py3py4);
        deg2 = y0py1py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_3_4 = dec.add_to_op_selec(y0py1py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_1_3_4.begin(), imp_0_1_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_1_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py2py3py4 ;
        POLY_ADD(y[0],y2py3py4,y0py2py3py4);
        deg2 = y0py2py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_3_4 = dec.add_to_op_selec(y0py2py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_2_3_4.begin(), imp_0_2_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_2_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y1py2py3py4 ;
        POLY_ADD(y[1],y2py3py4,y1py2py3py4);
        deg2 = y1py2py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_3_4 = dec.add_to_op_selec(y1py2py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_1_2_3_4.begin(), imp_1_2_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_1_2_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }

        poly y0py1py2py3py4 ;
        POLY_ADD(y[0],y1py2py3py4,y0py1py2py3py4);
        deg2 = y0py1py2py3py4.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_4 = dec.add_to_op_selec(y0py1py2py3py4, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_4.insert(imp_4.end(),imp_0_1_2_3_4.begin(), imp_0_1_2_3_4.end());
                T4.push_back(imp_4.size());
            }
            else { 
                deg = deg2; 
                imp_4 = imp_0_1_2_3_4; 
                for (uint32_t i=0; i<T4.size(); i++){
                    T4[i] = 0;
                }
                T4.push_back(imp_4.size());
            }
        }
        else {
            T4.push_back(0);
        }
        *a += nb_mult(deg, *test);

        cout<<"y4 done : "<<imp_4.size()<<endl;
        #pragma omp critical
        {
            insert_to_size_vect(imp_4, imp, T4, T);
        }
    }

    }

    #pragma omp section
    {
    
    if (size_out > 5){
    //y5 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T5;
        uint32_t deg = y[5].algebraic_degree(nb_elem);
        decomposition dec5(size_in,deg);

        vector<implem> imp_5 = dec5.add_to_op_selec(y[5], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T5.push_back(imp_5.size());

        poly y0py5 ;
        POLY_ADD(y[0],y[5],y0py5);
        uint32_t deg2 = y0py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_5 = dec.add_to_op_selec(y0py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_5.begin(), imp_0_5.end());       
                T5.push_back(imp_5.size());
            } 
            else { 
                deg = deg2; 
                imp_5 = imp_0_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }       

        poly y1py5 ;
        POLY_ADD(y[1],y[5],y1py5);
        deg2 = y1py5.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_5 = dec.add_to_op_selec(y1py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_5.begin(), imp_1_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y2py5 ;
        POLY_ADD(y[2],y[5],y2py5);
        deg2 = y2py5.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_5 = dec.add_to_op_selec(y2py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_2_5.begin(), imp_2_5.end());  
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_2_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y3py5 ;
        POLY_ADD(y[3],y[5],y3py5);
        deg2 = y3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_3_5 = dec.add_to_op_selec(y3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_3_5.begin(), imp_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y4py5 ;
        POLY_ADD(y[4],y[5],y4py5);
        deg2 = y4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_4_5 = dec.add_to_op_selec(y4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_4_5.begin(), imp_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py5 ;
        POLY_ADD(y[0],y1py5,y0py1py5);
        deg2 = y0py1py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_5 = dec.add_to_op_selec(y0py1py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_5.begin(), imp_0_1_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py2py5 ;
        POLY_ADD(y[0],y2py5,y0py2py5);
        deg2 = y0py2py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_5 = dec.add_to_op_selec(y0py2py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_2_5.begin(), imp_0_2_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_2_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py3py5 ;
        POLY_ADD(y[0],y3py5,y0py3py5);
        deg2 = y0py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_3_5 = dec.add_to_op_selec(y0py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_3_5.begin(), imp_0_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py4py5 ;
        POLY_ADD(y[0],y4py5,y0py4py5);
        deg2 = y0py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_4_5 = dec.add_to_op_selec(y0py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_4_5.begin(), imp_0_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py2py5 ;
        POLY_ADD(y[1],y2py5,y1py2py5);
        deg2 = y1py2py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_5 = dec.add_to_op_selec(y1py2py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_2_5.begin(), imp_1_2_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_2_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py3py5 ;
        POLY_ADD(y[1],y3py5,y1py3py5);
        deg2 = y1py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_3_5 = dec.add_to_op_selec(y1py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_3_5.begin(), imp_1_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py4py5 ;
        POLY_ADD(y[1],y4py5,y1py4py5);
        deg2 = y1py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_4_5 = dec.add_to_op_selec(y1py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_4_5.begin(), imp_1_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y2py3py5 ;
        POLY_ADD(y[2],y3py5,y2py3py5);
        deg2 = y2py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_3_5 = dec.add_to_op_selec(y2py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_2_3_5.begin(), imp_2_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_2_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }
        poly y2py4py5 ;
        POLY_ADD(y[2],y4py5,y2py4py5);
        deg2 = y2py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_4_5 = dec.add_to_op_selec(y2py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_2_4_5.begin(), imp_2_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_2_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y3py4py5 ;
        POLY_ADD(y[3],y4py5,y3py4py5);
        deg2 = y3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_3_4_5 = dec.add_to_op_selec(y3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_3_4_5.begin(), imp_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py2py5 ;
        POLY_ADD(y[0],y1py2py5,y0py1py2py5);
        deg2 = y0py1py2py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2_5 = dec.add_to_op_selec(y0py1py2py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_2_5.begin(), imp_0_1_2_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_2_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py3py5 ;
        POLY_ADD(y[0],y1py3py5,y0py1py3py5);
        deg2 = y0py1py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_3_5 = dec.add_to_op_selec(y0py1py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_3_5.begin(), imp_0_1_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py4py5 ;
        POLY_ADD(y[0],y1py4py5,y0py1py4py5);
        deg2 = y0py1py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_4_5 = dec.add_to_op_selec(y0py1py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_4_5.begin(), imp_0_1_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py2py3py5 ;
        POLY_ADD(y[0],y2py3py5,y0py2py3py5);
        deg2 = y0py2py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_3_5 = dec.add_to_op_selec(y0py2py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_2_3_5.begin(), imp_0_2_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_2_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py2py4py5 ;
        POLY_ADD(y[0],y2py4py5,y0py2py4py5);
        deg2 = y0py2py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_4_5 = dec.add_to_op_selec(y0py2py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_2_4_5.begin(), imp_0_2_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_2_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py3py4py5 ;
        POLY_ADD(y[0],y3py4py5,y0py3py4py5);
        deg2 = y0py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_3_4_5 = dec.add_to_op_selec(y0py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_3_4_5.begin(), imp_0_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py2py3py5 ;
        POLY_ADD(y[1],y2py3py5,y1py2py3py5);
        deg2 = y1py2py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_3_5 = dec.add_to_op_selec(y1py2py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_2_3_5.begin(), imp_1_2_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2;                
                imp_5 = imp_1_2_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py2py4py5 ;
        POLY_ADD(y[1],y2py4py5,y1py2py4py5);
        deg2 = y1py2py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_4_5 = dec.add_to_op_selec(y1py2py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_2_4_5.begin(), imp_1_2_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_2_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py3py4py5 ;
        POLY_ADD(y[1],y3py4py5,y1py3py4py5);
        deg2 = y1py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_3_4_5 = dec.add_to_op_selec(y1py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_3_4_5.begin(), imp_1_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y2py3py4py5 ;
        POLY_ADD(y[2],y3py4py5,y2py3py4py5);
        deg2 = y2py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_2_3_4_5 = dec.add_to_op_selec(y2py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_2_3_4_5.begin(), imp_2_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_2_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py2py3py5 ;
        POLY_ADD(y[0],y1py2py3py5,y0py1py2py3py5);
        deg2 = y0py1py2py3py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2);

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            } 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_5 = dec.add_to_op_selec(y0py1py2py3py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_2_3_5.begin(), imp_0_1_2_3_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_2_3_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py2py4py5 ;
        POLY_ADD(y[0],y1py2py4py5,y0py1py2py4py5);
        deg2 = y0py1py2py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2_4_5 = dec.add_to_op_selec(y0py1py2py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_2_4_5.begin(), imp_0_1_2_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_2_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py3py4py5 ;
        POLY_ADD(y[0],y1py3py4py5,y0py1py3py4py5);
        deg2 = y0py1py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_3_4_5 = dec.add_to_op_selec(y0py1py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_3_4_5.begin(), imp_0_1_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py2py3py4py5 ;
        POLY_ADD(y[0],y2py3py4py5,y0py2py3py4py5);
        deg2 = y0py2py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_2_3_4_5 = dec.add_to_op_selec(y0py2py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_2_3_4_5.begin(), imp_0_2_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_2_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y1py2py3py4py5 ;
        POLY_ADD(y[1],y2py3py4py5,y1py2py3py4py5);
        deg2 = y1py2py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_1_2_3_4_5 = dec.add_to_op_selec(y1py2py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_1_2_3_4_5.begin(), imp_1_2_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_1_2_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }

        poly y0py1py2py3py4py5 ;
        POLY_ADD(y[0],y1py2py3py4py5,y0py1py2py3py4py5);
        deg2 = y0py1py2py3py4py5.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 

            if (deg2 < deg){
                *test = TAB_SIZE -1;
            }
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_4_5 = dec.add_to_op_selec(y0py1py2py3py4py5, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_5.insert(imp_5.end(),imp_0_1_2_3_4_5.begin(), imp_0_1_2_3_4_5.end());
                T5.push_back(imp_5.size());
            }
            else { 
                deg = deg2; 
                imp_5 = imp_0_1_2_3_4_5; 
                for (uint32_t i=0; i<T5.size(); i++){
                    T5[i] = 0;
                }
                T5.push_back(imp_5.size());
            }
        }
        else {
            T5.push_back(0);
        }
        
        *a += nb_mult(deg, *test);

        cout<<"y5 done : "<<imp_5.size()<<endl;

        #pragma omp critical
        {
            insert_to_size_vect(imp_5, imp, T5, T);
        }
    }
    
    if (size_out > 6){
    //y6 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T6;
        uint32_t deg = y[6].algebraic_degree(nb_elem);
        decomposition dec6(size_in,deg);

        vector<implem> imp_6 = dec6.add_to_op_selec(y[6], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T6.push_back(imp_6.size());

        poly y0py6 ;
        POLY_ADD(y[0],y[6],y0py6);
        uint32_t deg2 = y0py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_6 = dec.add_to_op_selec(y0py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_6.begin(), imp_0_6.end());       
                T6.push_back(imp_6.size());
            } 
            else { 
                deg = deg2; 
                imp_6 = imp_0_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }       

        poly y1py6 ;
        POLY_ADD(y[1],y[6],y1py6);
        deg2 = y1py6.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_6 = dec.add_to_op_selec(y1py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_6.begin(), imp_1_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py6 ;
        POLY_ADD(y[2],y[6],y2py6);
        deg2 = y2py6.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_6 = dec.add_to_op_selec(y2py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_6.begin(), imp_2_6.end());  
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y3py6 ;
        POLY_ADD(y[3],y[6],y3py6);
        deg2 = y3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_3_6 = dec.add_to_op_selec(y3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_3_6.begin(), imp_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y4py6 ;
        POLY_ADD(y[4],y[6],y4py6);
        deg2 = y4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_4_6 = dec.add_to_op_selec(y4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_4_6.begin(), imp_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y5py6 ;
        POLY_ADD(y[5],y[6],y5py6);
        deg2 = y5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg){

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_5_6 = dec.add_to_op_selec(y5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_5_6.begin(), imp_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py6 ;
        POLY_ADD(y[0],y1py6,y0py1py6);
        deg2 = y0py1py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_6 = dec.add_to_op_selec(y0py1py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_6.begin(), imp_0_1_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py6 ;
        POLY_ADD(y[0],y2py6,y0py2py6);
        deg2 = y0py2py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_6 = dec.add_to_op_selec(y0py2py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_6.begin(), imp_0_2_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py3py6 ;
        POLY_ADD(y[0],y3py6,y0py3py6);
        deg2 = y0py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_3_6 = dec.add_to_op_selec(y0py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_3_6.begin(), imp_0_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py4py6 ;
        POLY_ADD(y[0],y4py6,y0py4py6);
        deg2 = y0py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_4_6 = dec.add_to_op_selec(y0py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_4_6.begin(), imp_0_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py5py6 ;
        POLY_ADD(y[0],y5py6,y0py5py6);
        deg2 = y0py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_5_6 = dec.add_to_op_selec(y0py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_5_6.begin(), imp_0_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py6 ;
        POLY_ADD(y[1],y2py6,y1py2py6);
        deg2 = y1py2py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_6 = dec.add_to_op_selec(y1py2py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_6.begin(), imp_1_2_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py3py6 ;
        POLY_ADD(y[1],y3py6,y1py3py6);
        deg2 = y1py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_3_6 = dec.add_to_op_selec(y1py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_3_6.begin(), imp_1_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py4py6 ;
        POLY_ADD(y[1],y4py6,y1py4py6);
        deg2 = y1py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_4_6 = dec.add_to_op_selec(y1py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_4_6.begin(), imp_1_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py5py6 ;
        POLY_ADD(y[1],y5py6,y1py5py6);
        deg2 = y1py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_5_6 = dec.add_to_op_selec(y1py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_5_6.begin(), imp_1_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py3py6 ;
        POLY_ADD(y[2],y3py6,y2py3py6);
        deg2 = y2py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_3_6 = dec.add_to_op_selec(y2py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_3_6.begin(), imp_2_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py4py6 ;
        POLY_ADD(y[2],y4py6,y2py4py6);
        deg2 = y2py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_4_6 = dec.add_to_op_selec(y2py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_4_6.begin(), imp_2_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py5py6 ;
        POLY_ADD(y[2],y5py6,y2py5py6);
        deg2 = y2py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_5_6 = dec.add_to_op_selec(y2py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_5_6.begin(), imp_2_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y3py4py6 ;
        POLY_ADD(y[3],y4py6,y3py4py6);
        deg2 = y3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_3_4_6 = dec.add_to_op_selec(y3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_3_4_6.begin(), imp_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y3py5py6 ;
        POLY_ADD(y[3],y5py6,y3py5py6);
        deg2 = y3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_3_5_6 = dec.add_to_op_selec(y3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_3_5_6.begin(), imp_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y4py5py6 ;
        POLY_ADD(y[4],y5py6,y4py5py6);
        deg2 = y4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_4_5_6 = dec.add_to_op_selec(y4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_4_5_6.begin(), imp_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py6 ;
        POLY_ADD(y[0],y1py2py6,y0py1py2py6);
        deg2 = y0py1py2py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_6 = dec.add_to_op_selec(y0py1py2py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_6.begin(), imp_0_1_2_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py3py6 ;
        POLY_ADD(y[0],y1py3py6,y0py1py3py6);
        deg2 = y0py1py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_3_6 = dec.add_to_op_selec(y0py1py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_3_6.begin(), imp_0_1_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py4py6 ;
        POLY_ADD(y[0],y1py4py6,y0py1py4py6);
        deg2 = y0py1py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_4_6 = dec.add_to_op_selec(y0py1py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_4_6.begin(), imp_0_1_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py5py6 ;
        POLY_ADD(y[0],y1py5py6,y0py1py5py6);
        deg2 = y0py1py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_5_6 = dec.add_to_op_selec(y0py1py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_5_6.begin(), imp_0_1_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py3py6 ;
        POLY_ADD(y[0],y2py3py6,y0py2py3py6);
        deg2 = y0py2py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_3_6 = dec.add_to_op_selec(y0py2py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_3_6.begin(), imp_0_2_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py4py6 ;
        POLY_ADD(y[0],y2py4py6,y0py2py4py6);
        deg2 = y0py2py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_4_6 = dec.add_to_op_selec(y0py2py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_4_6.begin(), imp_0_2_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py5py6 ;
        POLY_ADD(y[0],y2py5py6,y0py2py5py6);
        deg2 = y0py2py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_5_6 = dec.add_to_op_selec(y0py2py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_5_6.begin(), imp_0_2_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py3py4py6 ;
        POLY_ADD(y[0],y3py4py6,y0py3py4py6);
        deg2 = y0py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_3_4_6 = dec.add_to_op_selec(y0py2py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_3_4_6.begin(), imp_0_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py3py5py6 ;
        POLY_ADD(y[0],y3py5py6,y0py3py5py6);
        deg2 = y0py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_3_5_6 = dec.add_to_op_selec(y0py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_3_5_6.begin(), imp_0_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py4py5py6 ;
        POLY_ADD(y[0],y4py5py6,y0py4py5py6);
        deg2 = y0py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_4_5_6 = dec.add_to_op_selec(y0py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_4_5_6.begin(), imp_0_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py3py6 ;
        POLY_ADD(y[1],y2py3py6,y1py2py3py6);
        deg2 = y1py2py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_3_6 = dec.add_to_op_selec(y1py2py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_3_6.begin(), imp_1_2_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py4py6 ;
        POLY_ADD(y[1],y2py4py6,y1py2py4py6);
        deg2 = y1py2py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_4_6 = dec.add_to_op_selec(y1py2py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_4_6.begin(), imp_1_2_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py5py6 ;
        POLY_ADD(y[1],y2py5py6,y1py2py5py6);
        deg2 = y1py2py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_5_6 = dec.add_to_op_selec(y1py2py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_5_6.begin(), imp_1_2_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py3py4py6 ;
        POLY_ADD(y[1],y3py4py6,y1py3py4py6);
        deg2 = y1py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_3_4_6 = dec.add_to_op_selec(y1py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_3_4_6.begin(), imp_1_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py3py5py6 ;
        POLY_ADD(y[1],y3py5py6,y1py3py5py6);
        deg2 = y1py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_3_5_6 = dec.add_to_op_selec(y1py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_3_5_6.begin(), imp_1_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py4py5py6 ;
        POLY_ADD(y[1],y4py5py6,y1py4py5py6);
        deg2 = y1py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_4_5_6 = dec.add_to_op_selec(y1py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_4_5_6.begin(), imp_1_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py3py4py6 ;
        POLY_ADD(y[2],y3py4py6,y2py3py4py6);
        deg2 = y2py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_3_4_6 = dec.add_to_op_selec(y2py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_3_4_6.begin(), imp_2_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py3py5py6 ;
        POLY_ADD(y[2],y3py5py6,y2py3py5py6);
        deg2 = y2py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_3_5_6 = dec.add_to_op_selec(y2py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_3_5_6.begin(), imp_2_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py4py5py6 ;
        POLY_ADD(y[2],y4py5py6,y2py4py5py6);
        deg2 = y2py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_4_5_6 = dec.add_to_op_selec(y2py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_4_5_6.begin(), imp_2_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y3py4py5py6 ;
        POLY_ADD(y[3],y4py5py6,y3py4py5py6);
        deg2 = y3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_3_4_5_6 = dec.add_to_op_selec(y3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_3_4_5_6.begin(), imp_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py3py6 ;
        POLY_ADD(y[0],y1py2py3py6,y0py1py2py3py6);
        deg2 = y0py1py2py3py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_6 = dec.add_to_op_selec(y0py1py2py3py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_3_6.begin(), imp_0_1_2_3_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_3_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py4py6 ;
        POLY_ADD(y[0],y1py2py4py6,y0py1py2py4py6);
        deg2 = y0py1py2py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_4_6 = dec.add_to_op_selec(y0py1py2py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_4_6.begin(), imp_0_1_2_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py5py6 ;
        POLY_ADD(y[0],y1py2py5py6,y0py1py2py5py6);
        deg2 = y0py1py2py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_5_6 = dec.add_to_op_selec(y0py1py2py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_5_6.begin(), imp_0_1_2_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py3py4py6 ;
        POLY_ADD(y[0],y1py3py4py6,y0py1py3py4py6);
        deg2 = y0py1py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_3_4_6 = dec.add_to_op_selec(y0py1py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_3_4_6.begin(), imp_0_1_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py3py5py6 ;
        POLY_ADD(y[0],y1py3py5py6,y0py1py3py5py6);
        deg2 = y0py1py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_3_5_6 = dec.add_to_op_selec(y0py1py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_3_5_6.begin(), imp_0_1_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py4py5py6 ;
        POLY_ADD(y[0],y1py4py5py6,y0py1py4py5py6);
        deg2 = y0py1py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_4_5_6 = dec.add_to_op_selec(y0py1py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_4_5_6.begin(), imp_0_1_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py3py4py6 ;
        POLY_ADD(y[0],y2py3py4py6,y0py2py3py4py6);
        deg2 = y0py2py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_3_4_6 = dec.add_to_op_selec(y0py2py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_3_4_6.begin(), imp_0_2_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py3py5py6 ;
        POLY_ADD(y[0],y2py3py5py6,y0py2py3py5py6);
        deg2 = y0py2py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_3_5_6 = dec.add_to_op_selec(y0py2py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_3_5_6.begin(), imp_0_2_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py4py5py6 ;
        POLY_ADD(y[0],y2py4py5py6,y0py2py4py5py6);
        deg2 = y0py2py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_4_5_6 = dec.add_to_op_selec(y0py2py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_4_5_6.begin(), imp_0_2_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py3py4py5py6 ;
        POLY_ADD(y[0],y3py4py5py6,y0py3py4py5py6);
        deg2 = y0py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_3_4_5_6 = dec.add_to_op_selec(y0py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_3_4_5_6.begin(), imp_0_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py3py4py6 ;
        POLY_ADD(y[1],y2py3py4py6,y1py2py3py4py6);
        deg2 = y1py2py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_3_4_6 = dec.add_to_op_selec(y1py2py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_3_4_6.begin(), imp_1_2_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py3py5py6 ;
        POLY_ADD(y[1],y2py3py5py6,y1py2py3py5py6);
        deg2 = y1py2py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_3_5_6 = dec.add_to_op_selec(y1py2py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_3_5_6.begin(), imp_1_2_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py4py5py6 ;
        POLY_ADD(y[1],y2py4py5py6,y1py2py4py5py6);
        deg2 = y1py2py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_4_5_6 = dec.add_to_op_selec(y1py2py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_4_5_6.begin(), imp_1_2_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py3py4py5py6 ;
        POLY_ADD(y[1],y3py4py5py6,y1py3py4py5py6);
        deg2 = y1py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_3_4_5_6 = dec.add_to_op_selec(y1py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_3_4_5_6.begin(), imp_1_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y2py3py4py5py6 ;
        POLY_ADD(y[2],y3py4py5py6,y2py3py4py5py6);
        deg2 = y2py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_2_3_4_5_6 = dec.add_to_op_selec(y2py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_2_3_4_5_6.begin(), imp_2_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_2_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py3py4py6 ;
        POLY_ADD(y[0],y1py2py3py4py6,y0py1py2py3py4py6);
        deg2 = y0py1py2py3py4py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_4_6 = dec.add_to_op_selec(y0py1py2py3py4py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_3_4_6.begin(), imp_0_1_2_3_4_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_3_4_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py3py5py6 ;
        POLY_ADD(y[0],y1py2py3py5py6,y0py1py2py3py5py6);
        deg2 = y0py1py2py3py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_5_6 = dec.add_to_op_selec(y0py1py2py3py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_3_5_6.begin(), imp_0_1_2_3_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_3_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py4py5py6 ;
        POLY_ADD(y[0],y1py2py4py5py6,y0py1py2py4py5py6);
        deg2 = y0py1py2py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_4_5_6 = dec.add_to_op_selec(y0py1py2py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_4_5_6.begin(), imp_0_1_2_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py3py4py5py6 ;
        POLY_ADD(y[0],y1py3py4py5py6,y0py1py3py4py5py6);
        deg2 = y0py1py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_3_4_5_6 = dec.add_to_op_selec(y0py1py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_3_4_5_6.begin(), imp_0_1_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py2py3py4py5py6 ;
        POLY_ADD(y[0],y2py3py4py5py6,y0py2py3py4py5py6);
        deg2 = y0py2py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_2_3_4_5_6 = dec.add_to_op_selec(y0py2py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_2_3_4_5_6.begin(), imp_0_2_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_2_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y1py2py3py4py5py6 ;
        POLY_ADD(y[1],y2py3py4py5py6,y1py2py3py4py5py6);
        deg2 = y1py2py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_1_2_3_4_5_6 = dec.add_to_op_selec(y1py2py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_1_2_3_4_5_6.begin(), imp_1_2_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_1_2_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }

        poly y0py1py2py3py4py5py6 ;
        POLY_ADD(y[0],y1py2py3py4py5py6,y0py1py2py3py4py5py6);
        deg2 = y0py1py2py3py4py5py6.algebraic_degree(nb_elem);

        if(deg2 <= deg)   {

            decomposition dec(size_in,deg2); 
            
            int t1 = *test;
            vector<implem> imp_0_1_2_3_4_5_6 = dec.add_to_op_selec(y0py1py2py3py4py5py6, l, l2, size_in, deg2, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);
            int t2 = *test;

            if (deg2 == deg && t1 == t2){
                imp_6.insert(imp_6.end(),imp_0_1_2_3_4_5_6.begin(), imp_0_1_2_3_4_5_6.end());
                T6.push_back(imp_6.size());
            }
            else { 
                deg = deg2; 
                imp_6 = imp_0_1_2_3_4_5_6; 
                for (uint32_t i=0; i<T6.size(); i++){
                    T6[i] = 0;
                }
                T6.push_back(imp_6.size());
            }
        }
        else {
            T6.push_back(0);
        }
        
        *a += nb_mult(deg, *test);

        #pragma omp critical
        {
            insert_to_size_vect(imp_6, imp, T6, T);
        }
    }
    }
    }

    /*if (size_out > 7){
    //Pour y7 :
        *test=TAB_SIZE -1;
        vector<uint32_t> T7;
        deg = y[7].algebraic_degree(nb_elem);
        decomposition dec7(size_in,deg);

        vector<implem> imp_7 = dec7.add_to_op_selec(y[7], l, l2, size_in, deg, nb_elem, set_op, map_xor, set_op_size, map_xor_size, test);

        T7.push_back(imp_7.size());
        *a += nb_mult(deg, *test);
        insert_to_size_vect(imp_7, imp, T7, T);
    }*/

    cout<<"Create_sets done "<<endl;
    return *imp;
}