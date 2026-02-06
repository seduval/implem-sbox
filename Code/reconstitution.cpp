#include "config.h"
#include "sbox.h"
#include "poly.h"
#include "precomputation.h"
#include "decomp.h"
#include "reconstitution.h"

static int stop = 0;

poly evaluate_implem(const implem &imp, uint32_t nb_elem)   {
    vector<string> tokens;
    string f = imp.formula;
    
    for (size_t i = 0; i < f.size(); )
    {
        if (isspace(f[i])) { i++; continue; }

        if (f[i] == '(' || f[i] == ')' || f[i] == '*' || f[i] == '^')
        {
            tokens.push_back(string(1, f[i]));
            i++;
        }
        else if (f[i] == 'p') // token pN
        {
            size_t j = i + 1;
            while (j < f.size() && isdigit(f[j])) j++;
            tokens.push_back(f.substr(i, j - i));
            i = j;
        }
        else
        {
            throw runtime_error("Invalid character in formula");
        }
    }

    vector<string> rpn;
    stack<string> opstack;

    auto priority = [](const string &op){
        if (op == "*") return 2;
        if (op == "^") return 1;
        return 0;
    };

    for (auto &tok : tokens)
    {
        if (tok[0] == 'p')  // operand
        {
            rpn.push_back(tok);
        }
        else if (tok == "*" || tok == "^")
        {
            while (!opstack.empty() &&
                   ((opstack.top() == "*" || opstack.top() == "^") &&
                    priority(opstack.top()) >= priority(tok)))
            {
                rpn.push_back(opstack.top());
                opstack.pop();
            }
            opstack.push(tok);
        }
        else if (tok == "(")
        {
            opstack.push(tok);
        }
        else if (tok == ")")
        {
            while (!opstack.empty() && opstack.top() != "(")
            {
                rpn.push_back(opstack.top());
                opstack.pop();
            }
            if (opstack.empty())
                throw runtime_error("Mismatched parentheses");
            opstack.pop(); // pop "("
        }
    }

    while (!opstack.empty())
    {
        if (opstack.top() == "(")
            throw runtime_error("Mismatched parentheses");
        rpn.push_back(opstack.top());
        opstack.pop();
    }

    stack<poly> st;

    for (auto &tok : rpn)
    {
        if (tok[0] == 'p')  // Operand pX
        {
            int idx = stoi(tok.substr(1)) - 1; // p1 → index 0
            st.push(imp.op_sol[idx]);
        }
        else
        {
            poly a = st.top(); st.pop();
            poly b = st.top(); st.pop();
            poly res;

            if (tok == "^")
            {
                POLY_ADD(a, b, res);
            }
            else if (tok == "*")
            {
                POLY_MUL(a, b, res, nb_elem);
            }
            else
            {
                throw runtime_error("Unknown operator");
            }

            st.push(res);
        }
    }

    return st.top();
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
                    str += " ^ ";
                }
                if (u==0){
                    str += "1";
                }
                else {
                    for (uint32_t bit_n=0; bit_n<size; bit_n++) {
                        if ((u>>bit_n)&1ul){ 
                            str += "x[" + to_string(bit_n) + "]";
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

string print_mon_quad_poly(poly p, uint32_t size) {
    uint32_t nb_elem = 1 << size;
    string str;

    for (uint32_t u = 0; u < nb_elem; u++) {
        if ((p.get(u >> 5) >> (u & 0x1f)) & 1ul) {

            int i = -1, j = -1;
            for (uint32_t bit = 0; bit < size; bit++) {
                if ((u >> bit) & 1ul) {
                    if (i == -1) i = bit;
                    else j = bit;
                }
            }

            str = "x[" + to_string(i) + "] & x[" + to_string(j) + "]";
            return str;  
        }
    }

    return "0";
}

string print_imp(const implem& imp, const unordered_map<poly,string>& dict) {
    string out;
    const string& f = imp.formula;

    for (size_t i = 0; i < f.size(); ++i)
    {
        if (f[i] == 'p') {
            size_t j = i+1;
            while (j < f.size() && isdigit(f[j])) j++;

            uint64_t idx = stoi(f.substr(i+1, j-(i+1))) - 1;

            if (idx >= imp.op_sol.size())
                throw runtime_error("Bad pX index");

            const poly& p = imp.op_sol[idx];

            auto it = dict.find(p);
            if (it == dict.end()){
                out += "(?)";
                i = j-1;
            }

            out += "(" + it->second + ")";
            i = j-1;
        }
        else {
            if (f[i] == '*'){
                out += "&";
            }
            else {
                out += f[i];
            }
            
        }
    }
    return out;
}

vector<vector<poly>> Permutations(const vector<poly>& v) {
    vector<vector<poly>> result;
    vector<poly> temp = v;  // Copie de v pour éviter de modifier l'original
    sort(temp.begin(), temp.end());

    do {
        result.push_back(temp);
    } while (next_permutation(temp.begin(), temp.end()));

    return result;
}

void genComb(int start, int k, int remaining, vector<int>& current, vector<vector<int>>& result) {
    if (remaining == 0) {
        result.push_back(current);
        return;
    }

    for (int i = start; i <= k - remaining; i++) {
        current.push_back(i);
        genComb(i + 1, k, remaining - 1, current, result);
        current.pop_back();
    }
}

// Generation of integer vectors representing the sequence of linear combinations of output bits that can be implemented. If k is 2, the possible combinations are {0}, {1} and {0,1}.
vector<vector<int>> generate_combinations(int k) {
    vector<vector<int>> result;
    vector<int> cur;

    for (int i = 1; i <= k; i++) {
        genComb(0, k, i, cur, result);
    }

    return result;
}

/** Return the output bit of the permutation that is implemented **/
uint32_t bit_num_to_real_order(const vector<vector<uint32_t>>* T, uint32_t bit_num){
    return log2((*T)[bit_num].size());
}

vector<uint32_t> bit_num_to_xor_sum(const vector<vector<uint32_t>>* T, uint32_t bit_num, uint32_t ind){
    vector<uint32_t> vect_res;

    uint32_t k = bit_num_to_real_order(T, bit_num);

    vect_res.push_back(k);

    if (k == 0){
        return vect_res;
    }

    vector<vector<int>> combos = generate_combinations(k);

    int32_t i=0;
    for (uint32_t j=0; j<(*T)[bit_num].size(); j++){
        if (ind < (*T)[bit_num][j]){
            i = j-1;
            break;
        }
    }

    if (i == -1){
        return vect_res;
    }
    else {
        for (uint32_t j=0; j<combos[i].size(); j++){
            vect_res.push_back(combos[i][j]);
        }
    }
    

    return vect_res;
}

string print_details_implem(const vector<vector<uint32_t>>* T, uint32_t bit_num, uint32_t ind, unordered_map<poly, string> polyToNames, vector<poly> y) {

    string to_print;

    vector<uint32_t> vect_real_order = bit_num_to_xor_sum(T, bit_num, ind);
    uint32_t k = vect_real_order[0];

    poly current_poly = y[k];
    to_print = "\tuint32_t ty" + polyToNames[current_poly] + " = ";

    for (uint32_t j=1; j<vect_real_order.size(); j++){
        poly poly_to_print = y[vect_real_order[j]];
        to_print += "ty" + polyToNames[poly_to_print] + " ^ ";
    }

    return to_print;
}

pair<vector<poly_quad>, map<poly_quad, vector<poly_quad>>> build_xor_basis(const set<poly_quad>& nums, uint32_t r) {
    size_t R = static_cast<size_t>(r);

    vector<poly_quad> basis_orig;               // éléments choisis parmi nums (taille <= R)
    vector<poly_quad> basis_reduced;            // formes réduites (pivots) correspondantes
    vector<vector<int>> basis_mask;             // basis_mask[i] length R :
                                                //   basis_reduced[i] = XOR_{j | basis_mask[i][j]==1} basis_orig[j]

    map<poly_quad, vector<poly_quad>> representations;

    for (const auto& num : nums) {
        poly_quad x = num;
        vector<int> comb(R, 0); // représentation courante de x en termes des basis_orig (index < basis_orig.size())

        // réduire x avec les pivots existants
        for (size_t i = 0; i < basis_reduced.size(); ++i) {
            poly_quad attempt = x ^ basis_reduced[i];
            // la comparaison (attempt < x) est la technique classique pour détecter si pivot réduit x.
            // elle suppose que operator< sur poly_quad correspond à un ordre compatible avec le "pivot" binaire.
            if (attempt < x) {
                x = attempt;
                // comb ^= basis_mask[i]
                for (size_t j = 0; j < R; ++j) comb[j] ^= basis_mask[i][j];
            }
        }

        if (x != poly_quad(0)) {
            // num est indépendant => on le prend comme élément de la base (obligatoire : base doit venir du set)
            if (basis_orig.size() >= R) {
                // théoriquement impossible 
                throw runtime_error("Error in rank computation");
            }

            size_t newIdx = basis_orig.size();
            basis_orig.push_back(num);
            basis_reduced.push_back(x);

            // basis_mask pour ce pivot = comb avec le bit newIdx mis à 1
            vector<int> mask = comb;
            mask[newIdx] ^= 1;
            basis_mask.push_back(mask);

            // représentation d'un vecteur de base = lui-même
            representations[num] = vector<poly_quad>{ num };
        } else {
            // dépendant : comb contient la représentation en termes des basis_orig déjà présents
            vector<poly_quad> repr;
            for (size_t j = 0; j < basis_orig.size(); ++j) {
                if (comb[j]) repr.push_back(basis_orig[j]);
            }
            representations[num] = repr;
        }
    }

    // ici on suppose basis_orig.size() == R (garanti)
    return { basis_orig, representations };
}

void rewrite_imp(implem& imp, map<poly_quad, vector<poly_quad>>& reps, const vector<poly_quad>& basis, uint32_t size, uint32_t nb_elem) {
    string f = imp.formula;
    bool changed = true;

    while (changed)
    {
        changed = false;

        // We search for pX ^ pY
        for (size_t i=0; i<f.size(); i++)   {
            if (f[i] != 'p') continue;

            // X
            size_t j = i+1;
            while (j < f.size() && isdigit(f[j])) j++;
            int idxA = stoi(f.substr(i+1, j-(i+1))) - 1;

            // Vérify if we have " ^ pY"
            size_t k = j;
            while (k < f.size() && isspace(f[k])) k++;
            if (k >= f.size() || f[k] != '^') continue;
            k++;

            while (k < f.size() && isspace(f[k])) k++;
            if (k >= f.size() || f[k] != 'p') continue;
            k++;

            // Y
            size_t m = k;
            while (m < f.size() && isdigit(f[m])) m++;
            // Les pi sont toujours croissant dans la formule donc pas besoin de lire la valeur de idxB.

            int idxB = idxA + 1;

            // Récupérer les deux polynômes
            poly pA = imp.op_sol[idxA];
            poly pB = imp.op_sol[idxB];

            uint32_t mA = pA.algebraic_degree(nb_elem);
            uint32_t mB = pB.algebraic_degree(nb_elem);

            if ((mA != 2) || (mB != 2)){
                continue;
            }

            poly_quad qA = poly_to_poly_quad(pA, size, nb_elem);
            poly_quad qB = poly_to_poly_quad(pB, size, nb_elem);

            // Construct q = qA XOR qB 
            poly_quad qnew = qA;
            qnew ^= qB; 

            bool update_reps = true;

            auto it = reps.find(qnew);
            if (it != reps.end()){
                update_reps = false;
            }

            if (update_reps){
                // Ajouter la décomposition de qnew dans reps si elle n'y ait pas déjà en enlevant les évenuels répet, Si qA = q1 ^ q2 et qB = q2 ^ q3, il faut que qnew = q1 ^ q3
                set<poly_quad> new_reps;

                for (uint32_t s=0; s<reps[qA].size(); s++){
                    new_reps.insert(reps[qA][s]);
                }
                for (uint32_t s=0; s<reps[qB].size(); s++){
                    auto it_s = new_reps.find(reps[qB][s]);
                    if (it_s != new_reps.end()){
                        new_reps.erase(it_s);
                    }
                    else {
                        new_reps.insert(reps[qB][s]);
                    }
                }
                
                vector<poly_quad> reps_qnew;
                reps_qnew.insert(reps_qnew.end(), new_reps.begin(), new_reps.end());
                reps[qnew] = reps_qnew;

            }

            // Convert to poly
            poly pnew(qnew, size);

            // Remplacer dans op_sol
            imp.op_sol[idxA] = pnew;
            imp.op_sol.erase(imp.op_sol.begin() + idxB);
           
            // remplacer "pX ^ pY" par "p(idxB+1)"
            string newVar = "p" + to_string(idxA+1);

            f = f.substr(0, i) + newVar + f.substr(m);

            // Mettre à jour tous les pZ suivants (index décrémenté après idxB)
            for (uint64_t z = idxB+1; z <= imp.op_sol.size()+1; z++)
            {
                string oldpz = "p" + to_string(z);
                string newpz = "p" + to_string(z-1);
                size_t pos = 0;
                while ((pos = f.find(oldpz, pos)) != string::npos)
                {
                    f.replace(pos, oldpz.size(), newpz);
                    pos += newpz.size();
                }
            }

            changed = true;
            break;
        }
    }

    imp.formula = f;
}

uint32_t retrieve_linear_from_quad_basis(const vector<poly_quad>& basis, vector<poly>& linear_op, unordered_map<poly, string>& polyToNames, vector<poly> l, uint32_t size_in, uint32_t nb_elem, uint32_t counter){
    uint32_t counter_l = counter;
    uint32_t monomial_order [nb_elem];

    for (size_t i = 0; i < basis.size(); i++) {
        poly p(basis[i], size_in);
        if (weight(basis[i]) != 1){ 
            bool stop2 = false;
            for (uint32_t li = 0; li < l.size() && !stop2; li++){
                for (uint32_t jl=0; jl<size_in && !stop2; jl++){
                    if (mon_dom_in(l[li],jl, size_in, nb_elem)){
                        create_monomial_order(monomial_order, nb_elem, size_in, jl);

                        pair<poly,poly> L;
                        L = poly_div_no_truncate(p, l[li], monomial_order, size_in, nb_elem);

                        uint32_t m1 = L.first.algebraic_degree(nb_elem);
                        uint32_t m2 = L.second.algebraic_degree(nb_elem);

                        if ((m1 == 1) && (m2 <= 1)){

                            auto it2 = polyToNames.find(l[li]);
                            if (it2 == polyToNames.end()){
                                linear_op.push_back(l[li]);
                                polyToNames[l[li]] = "l" + to_string(counter_l);
                                counter_l++;
                            }

                            it2 = polyToNames.find(L.first);
                            if (it2 == polyToNames.end()){
                                linear_op.push_back(L.first);
                                polyToNames[L.first] = "l" + to_string(counter_l);
                                counter_l++;
                            }

                            if (m2 == 1){

                                it2 = polyToNames.find(L.second);
                                if (it2 == polyToNames.end()){
                                    linear_op.push_back(L.second);
                                    polyToNames[L.second] = "l" + to_string(counter_l);
                                    counter_l++;
                                } 

                                string to_print_for_qi = polyToNames[l[li]] + " & " + polyToNames[L.first] + " ^ " + polyToNames[L.second];
                                polyToNames[p] = to_print_for_qi;
                            }
                            else {
                                string to_print_for_qi = polyToNames[l[li]] + " & " + polyToNames[L.first];
                                polyToNames[p] = to_print_for_qi;
                            }
                            stop2 = true;
                            continue;
                        }
                    }
                }
                if (li == l.size() -1){
                    cerr<<"No solutions found to decompose a quadratic polynomial"<<endl;
                }
            }
        }
    }
    return counter_l;
}

void write_binary_matrix(const vector<poly>& vect_op, uint32_t size_in, const string& filename) {
    ofstream out(filename);
    if (!out) {
        throw runtime_error("Unable to open file" + filename);
    }

    string init_line = "1";
    out << init_line << '\n';

    init_line = to_string(vect_op.size()) + " " + to_string(size_in + 1);
    out << init_line << '\n';

    for (const poly& p : vect_op) {
        
        string line;
        line.reserve((size_in + 1) * 2);  // bit + space

        for (uint32_t j=0; j<size_in; j++) {

            uint32_t u = (1u << j);
            uint32_t ind_word  = u >> 5;      
            uint32_t pos = u & 31;     

            uint32_t bit = (p.data[ind_word] >> pos) & 1u;

            if (bit){
                line.push_back(char('1'));
            }
            else {
                line.push_back(char('0'));
            }
            line.push_back(char(' '));
        }

        uint32_t u = 0;
        uint32_t ind_word  = u >> 5;     
        uint32_t pos = u & 31;     
        uint32_t bit = (p.data[ind_word] >> pos) & 1u;

        if (bit){
            line.push_back(char('1'));
        }
        else {
            line.push_back(char('0'));
        }

        out << line << '\n';
    }


}

uint32_t read_lin_file(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Unable to open file " << filename << "\n";
        return 0; 
    }

    string line;
    uint32_t xorCount = 0;

    if (getline(file, line)) {
        const string prefix = "XorCount:";
        if (line.rfind(prefix, 0) == 0) { 
            try {
                xorCount = stoul(line.substr(prefix.size())); 
            } catch (...) {
                cerr << "Error in XorCount conversion\n";
                xorCount = 0;
            }
        }
    }

    
    while (getline(file, line)) {
        cout << line << "\n";
    }

    return xorCount;
}

void print_basis_and_reps(const vector<poly_quad>& basis, const map<poly_quad, vector<poly_quad>>& reps, size_t size_in, unordered_map<poly, string>& polyToNames, const string& filename) {
    // Print the basis
    for (size_t i = 0; i < basis.size(); i++) {
        poly p(basis[i], size_in);
        if (weight(basis[i]) == 1){
            cout<<"\tuint32_t q"<<i<<" = "<<print_mon_quad_poly(p, size_in)<< ";"<<endl;
        }
        else {
            cout<<"\tuint32_t q"<<i<<" = "<<polyToNames[p]<< ";"<<endl;
        }
        polyToNames[p] = "q" + to_string(i);
    }

    // Dictionnary basis -> index
    unordered_map<poly_quad, uint32_t> baseIndex;
    for (uint32_t i = 0; i < basis.size(); i++) {
        baseIndex[basis[i]] = i; // q1, q2, ...
    }

    uint32_t counter = basis.size();

    ofstream out(filename);
    if (!out) {
        throw runtime_error("Impossible d'ouvrir le fichier");
    }

    string init_line = "1";
    out << init_line << '\n';

    init_line = to_string(reps.size() - basis.size()) + " " + to_string(basis.size());
    out << init_line << '\n';

    for (auto& [num, comb] : reps) {

        if (baseIndex.count(num)) {
            continue;
        }

        poly p2(num, size_in);
        polyToNames[p2] = "q" + to_string(counter);
        counter++;
        
        string line;
        line.reserve(basis.size() * 2);  

        set<uint32_t> sorted_comb;

        for (size_t j = 0; j < comb.size(); j++) { 
            sorted_comb.insert(baseIndex[comb[j]]);
        }
        
        uint32_t last_pos = 0;
        for (auto j = sorted_comb.begin(); j != sorted_comb.end(); j++) {
            uint32_t pos = *j;
            for (uint32_t i=last_pos; i<pos; i++){
                line.push_back(char('0'));
                line.push_back(char(' '));
            }
            last_pos = pos + 1;
            
            line.push_back(char('1'));
            line.push_back(char(' '));
        }

        for (uint32_t i = last_pos; i< basis.size(); i++){
            line.push_back(char('0'));
            line.push_back(char(' '));
        }

        out << line << '\n';
    }
}

void return_implem(uint32_t size_in, uint32_t size_out, uint32_t nb_elem, uint32_t nb_and_max, uint32_t nb_sol, vector<poly> Y, vector<poly> l, vector<poly> l2, poly_quad set_op [], pair_xor map_xor [], uint32_t set_op_size, uint32_t map_xor_size, vector<poly> ANF){
    
    vector<vector<poly>> vect_perm = Permutations(ANF);

    vector<poly> perm; // Non truncated permutations

    vector<vector<uint32_t>> real_order;

    random_device rd;  // a seed source for the random number engine
    mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    shuffle (vect_perm.begin(), vect_perm.end(), gen);
    uint32_t compteur = 0;

    #pragma omp parallel shared(stop) num_threads(1)
    {
        bool stop_local = false;
        #pragma omp for schedule(dynamic) 
        for (uint32_t jj=0; jj<vect_perm.size(); jj++)   {
            if (stop_local) {
                continue;
            }
            if (stop) {
                stop_local = true;
                continue;
            }

            vector<poly> y; //Permutation of Y
            for (uint32_t s=0; s<vect_perm[jj].size(); s++){
                perm.push_back(vect_perm[jj][s]);
                y.push_back(truncate_lin(vect_perm[jj][s], size_in));
            }
        
            unordered_map<poly, string> polyToNames; 
            vector<poly> linear_op;

            poly zero;
            polyToNames[zero] = "0";

            uint32_t ind_perm_to_ind_anf[size_out] ;

            for (uint32_t i=0; i<size_out; i++){
                polyToNames[Y[i]] = to_string(i);
                for (uint32_t j=0; j<size_out; j++){
                    if (perm[i] == ANF[j]){
                        ind_perm_to_ind_anf[i] = j;
                    }
                }
            }

            for (uint32_t i=0; i<size_out; i++){
                if (perm[i] != ANF[ind_perm_to_ind_anf[i]]){
                    cerr<<"Error"<<endl;
                }
            }

            vector<vector<uint32_t>> Tp;
            vector<vector<uint32_t>> * T;
            T = &Tp;

            uint32_t ap = 0;
            uint32_t * a;
            a = &ap;

            uint32_t testp = 1;
            uint32_t * test;
            test = &testp;
            
            vector<vector<implem>> imp = create_sets(a, test, T, y, l, l2, size_in, size_out, nb_elem, set_op, map_xor, set_op_size, map_xor_size);

            #pragma omp critical
            {
                /*cout<<"Sets créés pour : "<<omp_get_thread_num()<<endl;
                cout<<imp[0].size()<<endl;
                cout<<imp[1].size()<<endl;
                cout<<imp[2].size()<<endl;
                cout<<imp[3].size()<<endl;
                if (size_out > 4){
                    cout<<imp[4].size()<<endl;
                }
                if (size_out > 5){
                    cout<<imp[5].size()<<endl;
                }
                if (size_out > 6){
                    cout<<imp[6].size()<<endl;
                }
                cout<<*a<<endl;   */
            }

            if (imp[0].size() == 0){
                #pragma omp critical
                {
                    cout<<"Next permutation for thread "<<omp_get_thread_num()<<" : an output bit has no solutions"<<endl;
                    for (uint32_t s=0; s<imp.size(); s++){
                        imp[s].clear();
                    }
                }
                continue;
            }

            uint32_t nb_and_min = *a;
            uint32_t nb_and = nb_and_max - nb_and_min ;
            uint32_t nb_xor = 0;
            uint32_t nb_and_final = nb_and_min;

            cout<<"The given s-box can be implemented by :"<<endl;
            cout<<endl;
            cout<<"uint32_t Sbox(uint32_t X, uint8_t size) {"<<endl;
            cout<<"\tuint32_t x[size];"<<endl;
            cout<<"\tfor(uint8_t b=0; b<size; b++)  {"<<endl;
            cout<<"\t\tx[b] = X&1ul;"<<endl;
            cout<<"\t\tX >>= 1;"<<endl;
            cout<<"\t}"<<endl;
            cout<<"\tuint32_t y[size];"<<endl;
            
            for (uint32_t i=0; i<imp[0].size() && !stop; i++)    {

                set<poly_quad> op_0 = imp[0][i].quad_sol;                                  

                set<poly_quad> temp_0 = op_0;
                uint32_t r0 = Rank(temp_0);

                if (r0 > nb_and ){
                    continue;
                }

                for (uint32_t j=0; j<imp[1].size() && !stop; j++)    {

                    set<poly_quad> op_1 = imp[1][j].quad_sol;

                    set<poly_quad> to_merge_0 = op_0;
                    op_1.merge(to_merge_0);

                    set<poly_quad> temp_1 = op_1;
                    uint32_t r1 = Rank(temp_1);

                    if (r1 > nb_and ){
                        continue;
                    }

                    for (uint32_t k=0; k<imp[2].size() && !stop; k++)    {

                        set<poly_quad> op_2 = imp[2][k].quad_sol;

                        set<poly_quad> to_merge_1 = op_1;
                        op_2.merge(to_merge_1);
                        
                        set<poly_quad> temp_2 = op_2;
                        uint32_t r2 = Rank(temp_2);

                        if (r2 > nb_and ){
                            continue;
                        }

                        for (uint32_t m=0; m<imp[3].size() && !stop; m++)    {

                            set<poly_quad> op_3 = imp[3][m].quad_sol;

                            set<poly_quad> to_merge_2 = op_2;
                            op_3.merge(to_merge_2);

                            set<poly_quad> temp_3 = op_3;
                            uint32_t r3 = Rank(temp_3); 

                            if (r3 > nb_and ){
                                continue;
                            }

                            if (size_out == 4){
                                if (r3 <= nb_and)  {
                                    #pragma omp critical 
                                    {
                                        compteur++;

                                        vector<pair<size_t, size_t>> index = {{0, i}, {1, j}, {2, k}, {3, m}};

                                                    /* Contain the actual implementations, may be a XOR sum */
                                                    vector<poly> implem_evaluated;
                                                    for (const auto& [id, num] : index) {
                                                        implem_evaluated.push_back(evaluate_implem(imp[id][num], nb_elem));
                                                        real_order.push_back(bit_num_to_xor_sum(T, id, num));
                                                    }

                                                    /* counter to print the linear operations */
                                                    uint32_t counter = 0;

                                                    /* To retrieve the linear delta between the ANF and what is actually implemented */
                                                    
                                                    poly linear_parts_of_ob [size_out] ;

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        poly p;
                                                        uint32_t i = real_order[j][0];

                                                        p.add(implem_evaluated[j]);
                                                        for (uint32_t s=1; s<real_order[j].size(); s++){
                                                           p.add(y[real_order[j][s]]);
                                                        }

                                                        poly poly_in_ANF = ANF[ind_perm_to_ind_anf[i]];
                                                        p.add(poly_in_ANF);

                                                        if (p.algebraic_degree(nb_elem) > 1){
                                                            cerr<<"Error in linear parts"<<endl;
                                                            /*p.print_poly(size_in);
                                                            cout<<endl;*/
                                                        }
                                                        else {
                                                            linear_parts_of_ob[ind_perm_to_ind_anf[i]] = p;
                                                            if ( (p != zero) && (!p.is_linear_monomial(size_in)) ){
                                                                linear_op.push_back(p);
                                                                polyToNames[p] = "l" + to_string(counter);
                                                                counter++;
                                                            }
                                                            else{
                                                                polyToNames[p] = print_poly2(p, size_in);
                                                            }
                                                        }
                                                    }

                                                    nb_and_final += r3;

                                                    auto [basis, reps] = build_xor_basis(op_3, r3); 

                                                    counter = retrieve_linear_from_quad_basis(basis, linear_op, polyToNames, l, size_in, nb_elem, counter);

                                                    for (const auto& [id, num] : index) {
                                                        for (uint32_t idx = 0; idx < imp[id][num].op_sol.size(); idx++) {
                                                            if (imp[id][num].op_sol[idx].algebraic_degree(nb_elem) == 1) {
                                                                nb_xor++;
                                                                if (!imp[id][num].op_sol[idx].is_linear_monomial(size_in)) {
                                                                    auto it = polyToNames.find(imp[id][num].op_sol[idx]);
                                                                    if (it == polyToNames.end()) {
                                                                        linear_op.push_back(imp[id][num].op_sol[idx]);
                                                                        polyToNames[imp[id][num].op_sol[idx]] = "l" + to_string(counter);
                                                                        counter++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }


                                                    string linear_filename = "mat_lin.txt";
                                                    write_binary_matrix(linear_op, size_in, linear_filename);

                                                    string cmd = "./my_slp_heuristic_lin < mat_lin.txt > res_lin.txt";
                                                    int ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command: : " + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_lin.txt");

                                                    for (const auto& [id, num] : index) {
                                                        rewrite_imp(imp[id][num], reps, basis, size_in, nb_elem);
                                                    }

                                                    string quadratic_filename = "mat_quad.txt";
                                                    print_basis_and_reps(basis, reps, size_in, polyToNames, quadratic_filename);

                                                    cmd = "./my_slp_heuristic_quad < mat_quad.txt > res_quad.txt";
                                                    ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command:" + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_quad.txt");

                                                    for (const auto& [id, num] : index) {
                                                        string detail_implems = print_details_implem(T, id, num, polyToNames, y);
                                                        cout<<detail_implems;
                                                        string expr_y = print_imp(imp[id][num], polyToNames);
                                                        cout<<expr_y + ";"<<endl;
                                                    }

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        uint32_t i = real_order[j][0];
                                                        if (linear_parts_of_ob[ind_perm_to_ind_anf[i]] != zero){
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]];
                                                            for (uint32_t s=0; s<real_order[j].size(); s++){
                                                                cout<<" ^ " <<polyToNames[linear_parts_of_ob[ind_perm_to_ind_anf[real_order[j][s]]]];
                                                                nb_xor++;
                                                            }
                                                            cout<<";"<<endl;
                                                            
                                                        }
                                                        else {
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]]<<";"<<endl;
                                                        }
                                                        
                                                    } 
                                    }

                                    if (compteur >= nb_sol){
                                        #pragma omp atomic write
                                        stop = 1;
                                        //cout<<"Number of wanted solutions reached"<<endl;
                                    }
                                }
                            }
                            else {
                                if (size_out == 5){

                                    for (uint32_t p=0; p<imp[4].size() && !stop; p++)    {

                                        set<poly_quad> op_4 = imp[4][p].quad_sol;
                                        set<poly_quad> to_merge_3 = op_3;
                                        op_4.merge(to_merge_3);

                                        set<poly_quad> temp_4 = op_4;
                                        uint32_t r4 = Rank(temp_4);

                                        if (r4 <= nb_and)  {
                                            #pragma omp critical 
                                            {
                                                compteur++;

                                                vector<pair<size_t, size_t>> index = {{0, i}, {1, j}, {2, k}, {3, m}, {4, p}};

                                                    /* Contain the actual implementations, may be a XOR sum */
                                                    vector<poly> implem_evaluated;
                                                    for (const auto& [id, num] : index) {
                                                        implem_evaluated.push_back(evaluate_implem(imp[id][num], nb_elem));
                                                        real_order.push_back(bit_num_to_xor_sum(T, id, num));
                                                    }

                                                    /* counter to print the linear operations */
                                                    uint32_t counter = 0;

                                                    /* To retrieve the linear delta between the ANF and what is actually implemented */
                                                    
                                                    poly linear_parts_of_ob [size_out] ;

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        poly p;
                                                        uint32_t i = real_order[j][0];

                                                        p.add(implem_evaluated[j]);
                                                        for (uint32_t s=1; s<real_order[j].size(); s++){
                                                           p.add(y[real_order[j][s]]);
                                                        }

                                                        poly poly_in_ANF = ANF[ind_perm_to_ind_anf[i]];
                                                        p.add(poly_in_ANF);

                                                        if (p.algebraic_degree(nb_elem) > 1){
                                                            cerr<<"Error in linear parts"<<endl;
                                                            /*p.print_poly(size_in);
                                                            cout<<endl;*/
                                                        }
                                                        else {
                                                            linear_parts_of_ob[ind_perm_to_ind_anf[i]] = p;
                                                            if ( (p != zero) && (!p.is_linear_monomial(size_in)) ){
                                                                linear_op.push_back(p);
                                                                polyToNames[p] = "l" + to_string(counter);
                                                                counter++;
                                                            }
                                                            else{
                                                                polyToNames[p] = print_poly2(p, size_in);
                                                            }
                                                        }
                                                    }

                                                    nb_and_final += r4;

                                                    auto [basis, reps] = build_xor_basis(op_4, r4); 

                                                    counter = retrieve_linear_from_quad_basis(basis, linear_op, polyToNames, l, size_in, nb_elem, counter);

                                                    for (const auto& [id, num] : index) {
                                                        for (uint32_t idx = 0; idx < imp[id][num].op_sol.size(); idx++) {
                                                            if (imp[id][num].op_sol[idx].algebraic_degree(nb_elem) == 1) {
                                                                nb_xor++;
                                                                if (!imp[id][num].op_sol[idx].is_linear_monomial(size_in)) {
                                                                    auto it = polyToNames.find(imp[id][num].op_sol[idx]);
                                                                    if (it == polyToNames.end()) {
                                                                        linear_op.push_back(imp[id][num].op_sol[idx]);
                                                                        polyToNames[imp[id][num].op_sol[idx]] = "l" + to_string(counter);
                                                                        counter++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }


                                                    string linear_filename = "mat_lin.txt";
                                                    write_binary_matrix(linear_op, size_in, linear_filename);

                                                    string cmd = "./my_slp_heuristic_lin < mat_lin.txt > res_lin.txt";
                                                    int ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command: : " + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_lin.txt");

                                                    for (const auto& [id, num] : index) {
                                                        rewrite_imp(imp[id][num], reps, basis, size_in, nb_elem);
                                                    }

                                                    string quadratic_filename = "mat_quad.txt";
                                                    print_basis_and_reps(basis, reps, size_in, polyToNames, quadratic_filename);

                                                    cmd = "./my_slp_heuristic_quad < mat_quad.txt > res_quad.txt";
                                                    ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command:" + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_quad.txt");

                                                    for (const auto& [id, num] : index) {
                                                        string detail_implems = print_details_implem(T, id, num, polyToNames, y);
                                                        cout<<detail_implems;
                                                        string expr_y = print_imp(imp[id][num], polyToNames);
                                                        cout<<expr_y + ";"<<endl;
                                                    }

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        uint32_t i = real_order[j][0];
                                                        if (linear_parts_of_ob[ind_perm_to_ind_anf[i]] != zero){
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]];
                                                            for (uint32_t s=0; s<real_order[j].size(); s++){
                                                                cout<<" ^ " <<polyToNames[linear_parts_of_ob[ind_perm_to_ind_anf[real_order[j][s]]]];
                                                                nb_xor++;
                                                            }
                                                            cout<<";"<<endl;
                                                            
                                                        }
                                                        else {
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]]<<";"<<endl;
                                                        }
                                                        
                                                    } 

                                                if (compteur >= nb_sol){
                                                    stop = 1;
                                                    //cout<<"Number of wanted solutions reached"<<endl;
                                                }
                                            }
                                        }
                                    }
                                }
                                else {
                                    if (size_out == 6){

                                        for (uint32_t p=0; p<imp[4].size() && !stop; p++)    {

                                            set<poly_quad> op_4 = imp[4][p].quad_sol;

                                            set<poly_quad> to_merge_3 = op_3;
                                            op_4.merge(to_merge_3);

                                            set<poly_quad> temp_4 = op_4;
                                            uint32_t r4 = Rank(temp_4); 

                                            if (r4 > nb_and ){
                                                continue;
                                            }

                                            for (uint32_t q=0; q<imp[5].size() && !stop; q++)    {

                                                set<poly_quad> op_5 = imp[5][q].quad_sol;

                                                set<poly_quad> to_merge_4 = op_4;
                                                op_5.merge(to_merge_4);

                                                set<poly_quad> temp_5 = op_5;
                                                uint32_t r = Rank(temp_5);

                                                if (r <= nb_and)  {
                                                    #pragma omp critical 
                                                    {
                                                    compteur++;
   
                                                    vector<pair<size_t, size_t>> index = {{0, i}, {1, j}, {2, k}, {3, m}, {4, p}, {5, q}};

                                                    /* Contain the actual implementations, may be a XOR sum */
                                                    vector<poly> implem_evaluated;
                                                    for (const auto& [id, num] : index) {
                                                        implem_evaluated.push_back(evaluate_implem(imp[id][num], nb_elem));
                                                        real_order.push_back(bit_num_to_xor_sum(T, id, num));
                                                    }

                                                    /* counter to print the linear operations */
                                                    uint32_t counter = 0;

                                                    /* To retrieve the linear delta between the ANF and what is actually implemented */
                                                    
                                                    poly linear_parts_of_ob [size_out] ;

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        poly p;
                                                        uint32_t i = real_order[j][0];

                                                        p.add(implem_evaluated[j]);
                                                        for (uint32_t s=1; s<real_order[j].size(); s++){
                                                           p.add(y[real_order[j][s]]);
                                                        }

                                                        poly poly_in_ANF = ANF[ind_perm_to_ind_anf[i]];
                                                        p.add(poly_in_ANF);

                                                        if (p.algebraic_degree(nb_elem) > 1){
                                                            cerr<<"Error in linear parts"<<endl;
                                                            /*p.print_poly(size_in);
                                                            cout<<endl;*/
                                                        }
                                                        else {
                                                            linear_parts_of_ob[ind_perm_to_ind_anf[i]] = p;
                                                            if ( (p != zero) && (!p.is_linear_monomial(size_in)) ){
                                                                linear_op.push_back(p);
                                                                polyToNames[p] = "l" + to_string(counter);
                                                                counter++;
                                                            }
                                                            else{
                                                                polyToNames[p] = print_poly2(p, size_in);
                                                            }
                                                        }
                                                    }

                                                    nb_and_final += r;

                                                    auto [basis, reps] = build_xor_basis(op_5, r); 

                                                    counter = retrieve_linear_from_quad_basis(basis, linear_op, polyToNames, l, size_in, nb_elem, counter);

                                                    for (const auto& [id, num] : index) {
                                                        for (uint32_t idx = 0; idx < imp[id][num].op_sol.size(); idx++) {
                                                            if (imp[id][num].op_sol[idx].algebraic_degree(nb_elem) == 1) {
                                                                nb_xor++;
                                                                if (!imp[id][num].op_sol[idx].is_linear_monomial(size_in)) {
                                                                    auto it = polyToNames.find(imp[id][num].op_sol[idx]);
                                                                    if (it == polyToNames.end()) {
                                                                        linear_op.push_back(imp[id][num].op_sol[idx]);
                                                                        polyToNames[imp[id][num].op_sol[idx]] = "l" + to_string(counter);
                                                                        counter++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }


                                                    string linear_filename = "mat_lin.txt";
                                                    write_binary_matrix(linear_op, size_in, linear_filename);

                                                    string cmd = "./my_slp_heuristic_lin < mat_lin.txt > res_lin.txt";
                                                    int ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command: : " + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_lin.txt");

                                                    for (const auto& [id, num] : index) {
                                                        rewrite_imp(imp[id][num], reps, basis, size_in, nb_elem);
                                                    }

                                                    string quadratic_filename = "mat_quad.txt";
                                                    print_basis_and_reps(basis, reps, size_in, polyToNames, quadratic_filename);

                                                    cmd = "./my_slp_heuristic_quad < mat_quad.txt > res_quad.txt";
                                                    ret = system(cmd.c_str());

                                                    if ( ret != 0 ) {
                                                        cerr<<"Unable to run the command:" + cmd<<endl;
                                                    }

                                                    nb_xor += read_lin_file("res_quad.txt");

                                                    for (const auto& [id, num] : index) {
                                                        string detail_implems = print_details_implem(T, id, num, polyToNames, y);
                                                        cout<<detail_implems;
                                                        string expr_y = print_imp(imp[id][num], polyToNames);
                                                        cout<<expr_y + ";"<<endl;
                                                    }

                                                    for (uint32_t j=0; j<size_out; j++){
                                                        uint32_t i = real_order[j][0];
                                                        if (linear_parts_of_ob[ind_perm_to_ind_anf[i]] != zero){
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]];
                                                            for (uint32_t s=0; s<real_order[j].size(); s++){
                                                                cout<<" ^ " <<polyToNames[linear_parts_of_ob[ind_perm_to_ind_anf[real_order[j][s]]]];
                                                                nb_xor++;
                                                            }
                                                            cout<<";"<<endl;
                                                            
                                                        }
                                                        else {
                                                            cout<<"\ty["<<polyToNames[y[i]]<<"] = ty"<<polyToNames[y[i]]<<";"<<endl;
                                                        }
                                                        
                                                    } 

                                                    if (compteur >= nb_sol){
                                                        stop = 1;
                                                        //cout<<"Number of wanted solutions reached"<<endl;
                                                    }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else {
                                        if (size_out ==7){
                                            for (uint32_t p=0; p<imp[4].size() && !stop; p++)    {

                                                set<poly_quad> op_4 = imp[4][p].quad_sol;

                                                set<poly_quad> to_merge_3 = op_3;
                                                op_4.merge(to_merge_3);

                                                set<poly_quad> temp_4 = op_4;
                                                uint32_t r4 = Rank(temp_4); 

                                                if (r4 > nb_and-1){
                                                    continue;
                                                }

                                                for (uint32_t q=0; q<imp[5].size() && !stop; q++)    {

                                                    set<poly_quad> op_5 = imp[5][q].quad_sol;

                                                    set<poly_quad> to_merge_4 = op_4;
                                                    op_5.merge(to_merge_4);

                                                    set<poly_quad> temp_5 = op_5;
                                                    uint32_t r5 = Rank(temp_5);

                                                    if (r5 > nb_and-1){
                                                        continue;
                                                    }

                                                    for (uint32_t s=0; s<imp[6].size() && !stop; s++){

                                                        set<poly_quad> op_6 = imp[6][s].quad_sol;

                                                        set<poly_quad> to_merge_5 = op_5;
                                                        op_6.merge(to_merge_5);

                                                        set<poly_quad> temp_6 = op_6;
                                                        uint32_t r = Rank(temp_6);

                                                        if (r <= nb_and)  {

                                                            #pragma omp critical 
                                                            {
                                                            compteur++;
                                                           

                                                            if (compteur >= nb_sol){
                                                                stop = 1;
                                                                //cout<<"Number of wanted solutions reached"<<endl;
                                                            }
                                                            }
                                                        }
                                                    }
                                                }
                                            }      
                                        }
                                        else {
                                            if (size_out == 8){
                                                for (uint32_t p=0; p<imp[4].size() && !stop; p++)    {

                                                    set<poly_quad> op_4 = imp[4][p].quad_sol;

                                                    set<poly_quad> to_merge_3 = op_3;
                                                    op_4.merge(to_merge_3);

                                                    set<poly_quad> temp_4 = op_4;
                                                    uint32_t r4 = Rank(temp_4); 

                                                    if (r4 > nb_and){
                                                        continue;
                                                    }

                                                    for (uint32_t q=0; q<imp[5].size() && !stop; q++)    {

                                                        set<poly_quad> op_5 = imp[5][q].quad_sol;

                                                        set<poly_quad> to_merge_4 = op_4;
                                                        op_5.merge(to_merge_4);

                                                        set<poly_quad> temp_5 = op_5;
                                                        uint32_t r5 = Rank(temp_5);

                                                        if (r5 > nb_and){
                                                            continue;
                                                        }

                                                        for (uint32_t s=0; s<imp[6].size() && !stop; s++){

                                                            set<poly_quad> op_6 = imp[6][s].quad_sol;

                                                            set<poly_quad> to_merge_5 = op_5;
                                                            op_6.merge(to_merge_5);

                                                            set<poly_quad> temp_6 = op_6;
                                                            uint32_t r6 = Rank(temp_6);

                                                            if (r6 > nb_and){
                                                                continue;
                                                            }

                                                            for (uint32_t u=0; u<imp[7].size() && !stop; u++){

                                                                set<poly_quad> op_7 = imp[7][u].quad_sol;

                                                                set<poly_quad> to_merge_6 = op_6;
                                                                op_7.merge(to_merge_6);

                                                                set<poly_quad> temp_7 = op_7;
                                                                uint32_t r = Rank(temp_7);

                                                                if (r <= nb_and)  {

                                                                    #pragma omp critical 
                                                                    {
                                                                    compteur++;
                                                                    

                                                                    if (compteur >= nb_sol){
                                                                        stop = 1;
                                                                        //cout<<"Number of wanted solutions reached"<<endl;
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

            cout<<"\tuint32_t Y = 0;"<<endl;
            cout<<"\tfor(uint8_t b=0; b<size; b++) {"<<endl;
            cout<<"\t\tY ^= ((Y>>b) ^ y[b]) << b;"<<endl;
            cout<<"\t}"<<endl;
            cout<<"\treturn Y;"<<endl;
            cout<<"}"<<endl;

            cout<<"Nb_and = "<<nb_and_final<<endl;
            cout<<"Nb_xor = "<<nb_xor<<endl;

            
            /*#pragma omp critical
            {
                cout<<"Next permutation for thread "<<omp_get_thread_num()<<endl;
            }*/
        }
    }
}