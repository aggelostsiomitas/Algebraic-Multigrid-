#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <mkl.h>


class CSR{
public:
    int n;
    int m;
    std::vector<int>col_idx;
    std::vector<int>row_ptr;
    std::vector<double>vals;
    sparse_matrix_t mkl_handle=nullptr;
public:
    void update_mkl_handle(){
        if(mkl_handle){
             mkl_sparse_destroy(mkl_handle);
             mkl_handle=nullptr;
        }
         mkl_sparse_d_create_csr(&mkl_handle, SPARSE_INDEX_BASE_ZERO,
                        n, m, row_ptr.data(), row_ptr.data() + 1,
                        col_idx.data(), vals.data());
        }
    

    CSR()=default;
    CSR(const CSR& other): row_ptr(other.row_ptr),col_idx(other.col_idx),
    vals(other.vals),n(other.n),m(other.m){
        if(other.mkl_handle){
             update_mkl_handle();
        }
    }
    CSR& operator=(const CSR& other){
        if(this!=&other){
            n=other.n;
            m=other.m;
            row_ptr=other.row_ptr;
            col_idx=other.col_idx;
            vals=other.vals;
            update_mkl_handle();
        }
        return *this;
    }

    ~CSR(){
        if(mkl_handle){
            mkl_sparse_destroy(mkl_handle);
        }
    }
};


class AMGSOLVER{
    public:

    //level matrices 
        struct Level{
            CSR A;
            CSR P;
        };

        //parameters for the method 
        std::vector<Level> levels;//a vector that will hold the matrices of its level
        int max_levels=5; //maximum levels created
        double theta=0.5; //theta varible for connections

        //get fucion to aquire the elements of A
        double get(const CSR&A,int row,int col){
            double value=0.0;
            for(size_t k=A.row_ptr[row];k<A.row_ptr[row+1];++k){
                if(A.col_idx[k]==col){
                    return A.vals[k];
                }
            }

            return 0.0;
        }

        //function that will compute the intersection between two vectors a,b
        static  std::vector<int> intersect(const std::vector<int>&a ,const std::vector<int>&b){
            std::vector<int> res; //vector that wioll sotre the result of the intersection 

            //copies of the two vectors
            std::vector<int>aa=a;
            std::vector<int>bb=b;

            //sort the vectors to avoid  wrong intersection 
            std::sort(aa.begin(),aa.end());
            std::sort(bb.begin(),bb.end());

            //compute the intersection 
            std::set_intersection(aa.begin(),aa.end(),bb.begin(),bb.end(),std::back_inserter(res));
            return res;
        }

        //function that will compute the difference between two vectors a,b
        static std::vector<int> set_difference(const std::vector<int>&a,const std::vector<int>&b){
            std::vector<int> res; //vector that will store the result 
            for(int x:a){
                //if x does not exist on b add it 
                if(std::find(b.begin(),b.end(),x)==b.end()){
                    res.push_back(x);
                }
            }
            return res;
        }

        //function that searches the i rw of the matrix and find the non zero elements
        //it returns two vectors with the cols and vals of these elements 
        std::pair<std::vector<int>,std::vector<double>> find(const CSR&A,int row){
            std::vector<int>cols;
            std::vector<double>vals;
            for(int k=A.row_ptr[row];k<A.row_ptr[row+1];++k){
                int j=A.col_idx[k]; //column
                if(j!=row){
                    cols.push_back(j);
                    vals.push_back(A.vals[k]);
                }
            }
            return {cols,vals};
        }

        std::pair<std::vector<std::vector<int>> ,std::vector<std::vector<int>>> connections(const CSR&A){
            //construct the S,ST matrices
            int n=A.n;
            std::vector<std::vector<int>>S(n);
            std::vector<std::vector<int>>ST(n);

            for(size_t i=0;i<n;i++){
                //search the i row 
                auto [cols,vals]=find(A,i);

                //strength=-vals
                std::vector<double>strength;
                strength.reserve(vals.size());
                for(double val:vals){
                    strength.push_back(-val);
                }

                //alpha = max(strength);
                double alpha = std::numeric_limits<double>::lowest();
                for(double v:strength){
                    if(v>alpha){
                        alpha=v;
                    }
                }

                std::vector<int> strong_idx;
                strong_idx.reserve(cols.size());
                for(size_t k=0;k<cols.size();++k ){
                    if(strength[k]>=theta*alpha){
                        strong_idx.push_back(cols[k]);
                    }
                }

                S[i]=strong_idx;

                for(int node:S[i]){
                    ST[node].push_back(i);
                }
            }
            return {S,ST};
        }

        //function that will erform the C/F split 
        //based on Runge-Stuben
        std::pair<std::vector<int>,std::vector<int>> coarsening(const std::vector<std::vector<int>>&S,
        const std::vector<std::vector<int>>&ST){

            int n=S.size();
            std::vector<int>state(n,0); //vectors that will store the state 
            std::vector<int>lambda(n,0);//vector that ill sotre the mabda values 
            //that will be use to perform the C/F splitting 

            for(size_t i=0;i<n;++i){
                lambda[i]=ST[i].size();
            }

            while(std::any_of(state.begin(),state.end(),[](int val){return val==0;})){
                std::vector<int>score=lambda;
                for(size_t i=0;i<score.size();++i){
                    if(state[i]!=0){
                        score[i]=-1;
                    }
                }

                int idx=std::distance(score.begin(),std::max_element(score.begin(),score.end()));
                state[idx]=1;

                for(int j:S[idx]){
                    if(state[j]==0){
                        state[j]=-1;
                    }
                    for(int k:ST[j]){
                        if(state[k]==0){
                            lambda[k]+=1;
                        }
                    }
                }
            }


            std::vector<int> Fpoints; //state==-1
            std::vector<int> Cpoints; //state==1
            for(size_t i=0;i<n;++i){
                if (state[i]==-1){Fpoints.push_back(i);}
                if (state[i]==1){Cpoints.push_back(i);}
            }

            //i=find(state==-1)
            for(int i:Fpoints){
                //j = intersect(S{i}, find(state == -1))'
                std::vector<int>Si_F=intersect(S[i],Fpoints);
                for(int j:Si_F){
                    std::vector<int>Sj_C=intersect(S[j],Cpoints);
                    std::vector<int>Si_Sj_C=intersect(S[i],Sj_C);
                    if(Si_Sj_C.empty()){
                        state[i]=1;
                    }
                }
            }

            std::vector<int>C;
            std::vector<int>F;
            for (int i = 0; i < n; i++) {
                if (state[i] == 1) C.push_back(i);
                if (state[i] == -1) F.push_back(i);
            }
            return {C,F};
        }

        //function that will perform the interpolation 
        CSR interpolate(
            const CSR&A,
            const std::vector<std::vector<int>>&S,
            const std::vector<int>C,
            const std::vector<int>F){

                int n=A.n;
                int nc=C.size();
                
                std::vector<int>map(n,-1);
                for(size_t k=0;k<nc;++k){
                    map[C[k]]=k;
                }

                std::vector<std::vector<int>> cols(n);
                std::vector<std::vector<double>> vals(n);

                   for (int k = 0; k < nc; ++k) {
                        int i = C[k];
                        cols[i].push_back(k);
                        vals[i].push_back(1.0);
                    }

                for(int i:F){
                    std::vector<int>Fi=intersect(S[i],F);
                    std::vector<int>Ci=intersect(S[i],C);

                    if(Ci.empty()){continue;}

                    //off_diag=fnd(A(i,:))
                    auto [col_off,val_off] =find(A,i);
                    std::vector<int> weak_nodes=set_difference(col_off,S[i]);

                    double eps=1e-14;
                    double down=0.0;
                    down=get(A,i,i);
                    for(int w:weak_nodes){
                        down+=get(A,i,w);
                    }
                     if (std::abs(down) < eps){continue;}
                    

                    
                    for(int j:Ci){
                        double up=get(A,i,j);
                        for(int m:Fi){
                            double denom2=0.0;
                            for(int d:Ci){
                                 denom2+=get(A,m,d);
                            }
                            if (std::abs(denom2)>eps){
                                up+=get(A,i,m)*get(A,m,j)/denom2;
                            }
                        }
                        double wij=-up/down;
                        cols[i].push_back(map[j]);
                        vals[i].push_back(wij);
                    }
                }

            CSR P;
            P.n = n;
            P.m = nc;
            P.row_ptr.assign(n + 1, 0);
            int nnz = 0;
            for (int i = 0; i < n; ++i) {
                nnz += (int)cols[i].size();
                P.row_ptr[i + 1] = nnz;
            }
            P.col_idx.resize(nnz);
            P.vals.resize(nnz);
            
            int idx = 0;
            for(size_t i=0;i<n;++i){
                for(size_t k=0;k<cols[i].size();++k){
                    P.col_idx[idx] = cols[i][k];
                    P.vals[idx] = vals[i][k];
                    idx++;
                }
            }
            cols.clear();
            vals.clear();

            return P;
    };

    //function that will perform the galerkin opearation Ac=P'AP
    sparse_matrix_t Galerkin(const sparse_matrix_t & A,const sparse_matrix_t & P){
       
        sparse_matrix_t Ac_mkl,AP=nullptr;

        //AP=A*P
        mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE,A,P,&AP);

        //Ac=P'*AP
        mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE,P,AP,&Ac_mkl);

        mkl_sparse_optimize(Ac_mkl);
        mkl_sparse_destroy(AP);
        return Ac_mkl;
    }

    void Gauss_seidel(const CSR&A,std::vector<double>& x,const std::vector<double>&b,int iter){
        int n=A.n;
        for(int it=0;it<iter;++it){
            for(int i=0;i<n;++i){
                //b[i]-sum)/A[i][i]
                double sum=0.0;
                double diag=0.0;
                for(int k=A.row_ptr[i];k<A.row_ptr[i+1];++k){
                    int j=A.col_idx[k];
                    if(j==i){diag=A.vals[k];}
                    else 
                    {
                        sum+=A.vals[k]*x[j];
                    }
                }
                x[i]=(b[i]-sum)/diag;
            }
        }
    }

    void damped_jacobi(
         sparse_matrix_t A,
         std::vector<double>&x,
         const std::vector<double>&b,
         double omega,
         int iter
    ){
        int n=x.size();

        //extrect CSR data 
        sparse_index_base_t indexing;
        MKL_INT rows,cols;
        MKL_INT *row_start,*row_end;
        MKL_INT *col_indx;
        double *values;

        mkl_sparse_d_export_csr(
            A,
            &indexing,
            &rows,
            &cols,
            &row_start,
            &row_end,
            &col_indx,
            &values
        );

        std::vector<double>D(n);
        for(size_t i=0;i<n;++i){
            for(size_t j=row_start[i];j<row_end[i];++j){
                if (col_indx[j]==i){
                    D[i]=values[j];
                    break;
                }
            }
        }


        //x=x+omega*(b-A*x)./D
        matrix_descr descrA;
        descrA.type=SPARSE_MATRIX_TYPE_GENERAL;

        //vector for the residual
        std::vector<double>r(n);

        for(int it=0;it<iter;++it){

            //r=b
            cblas_dcopy(n,b.data(),1,r.data(),1);

            //r=r-A*x
            mkl_sparse_d_mv(
                SPARSE_OPERATION_NON_TRANSPOSE,
                -1.0,
                A,
                descrA,
                x.data(),
                1.0,
                r.data()
            );

            //x=x+omega*r./D
            for(size_t i=0;i<n;++i){
                if(D[i]!=0.0){
                    x[i]+=omega*r[i]/D[i];
                }
            }

        }
    }


    CSR mkl_to_csr(sparse_matrix_t handle) {
    CSR result;
    sparse_index_base_t indexing;
    MKL_INT rows, cols;
    MKL_INT *row_start = nullptr, *row_end = nullptr, *col_idx = nullptr;
    double *vals = nullptr;

    mkl_sparse_d_export_csr(
        handle, &indexing, &rows, &cols, 
        &row_start, &row_end, &col_idx, &vals
    );

    result.n = rows;
    result.m = cols;

    int base = (indexing == SPARSE_INDEX_BASE_ONE) ? 1 : 0;
    
    // Copy row_ptr
    result.row_ptr.resize(rows + 1);
    for (int i = 0; i < rows; ++i) {
        result.row_ptr[i] = row_start[i] - base;
    }
    result.row_ptr[rows] = row_end[rows - 1] - base;

    // Copy col_idx and values
    int nnz = result.row_ptr[rows];
    result.col_idx.resize(nnz);
    result.vals.resize(nnz);

    for (int k = 0; k < nnz; ++k) {
        result.col_idx[k] = col_idx[k] - base;
        result.vals[k] = vals[k];
    }

    // Keep handle reference synchronized
    result.mkl_handle = handle;
    return result;
}


    void amg_setup(CSR& A0){
        levels.clear();
        Level lvl;
        lvl.A=A0;
        lvl.A.update_mkl_handle();
        levels.push_back(lvl);
        
        int lvl_idx = 0;
        while(levels[lvl_idx].A.n>50 &&lvl_idx<max_levels-1){
                      
            auto [S,ST]=connections(levels[lvl_idx].A);
            auto [C,F]=coarsening(S,ST);

            CSR P=interpolate(levels[lvl_idx].A,S,C,F);
            P.update_mkl_handle();

            levels[lvl_idx].P = P;
            sparse_matrix_t Ac_mkl = Galerkin(levels[lvl_idx].A.mkl_handle, P.mkl_handle);

            Level coarse_lvl;
            coarse_lvl.A = mkl_to_csr(Ac_mkl);

            
            levels.push_back(coarse_lvl);
            lvl_idx++;
        }
    }


    void v_cycle(int lvl_idx,std::vector<double>&x,const std::vector<double>&b){
         Level& lvl = levels[lvl_idx];
        int n = lvl.A.n;

        if(lvl_idx==levels.size()-1){
            Gauss_seidel(lvl.A,x,b,5);
            return;
        }
         //smooth
        Gauss_seidel(lvl.A,x,b,4);

        std::vector<double>r(n);
        //r=b
        cblas_dcopy(n,b.data(),1,r.data(),1);

        //r=r-A*x
        matrix_descr descr;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        mkl_sparse_d_mv(
            SPARSE_OPERATION_NON_TRANSPOSE,
            -1.0,
            lvl.A.mkl_handle,
            descr,
            x.data(),
            1.0,
            r.data()
        );

        //restrict r_c = P' * r;
        int nc = lvl.P.m;
        std::vector<double> r_c(nc, 0.0);
        mkl_sparse_d_mv(
            SPARSE_OPERATION_TRANSPOSE,
            1.0,
            lvl.P.mkl_handle,
            descr,
            r.data(),
            0.0,
            r_c.data()
        );

    //coarse correction
    // e_c = zeros(size(r_c));
    // e_c = V(AMG, level+1, e_c, r_c)
    std::vector<double>e_c(nc,0.0);
    v_cycle(lvl_idx+1,e_c,r_c);
    

    //prologate x=x+P*e_c
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, lvl.P.mkl_handle, descr, e_c.data(), 1.0, x.data());

    //smmoth again 
    Gauss_seidel(lvl.A,x,b,3);
    }

    
};


CSR poisson(int n){
        int N = n * n;
    CSR A;
    A.n = N;
    A.m = N;
    A.row_ptr.reserve(N + 1);

    int nnz = 0;
    A.row_ptr.push_back(nnz);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int row = i * n + j;

            if (i > 0) {
                A.col_idx.push_back((i - 1) * n + j);
                A.vals.push_back(-1.0);
                nnz++;
            }
            if (j > 0) {
                A.col_idx.push_back(i * n + (j - 1));
                A.vals.push_back(-1.0);
                nnz++;
            }
            // Diagonal entry
            A.col_idx.push_back(row);
            A.vals.push_back(4.0);
            nnz++;

            if (j < n - 1) {
                A.col_idx.push_back(i * n + (j + 1));
                A.vals.push_back(-1.0);
                nnz++;
            }
            if (i < n - 1) {
                A.col_idx.push_back((i + 1) * n + j);
                A.vals.push_back(-1.0);
                nnz++;
            }
            A.row_ptr.push_back(nnz);
        }
    }
    A.update_mkl_handle();
    return A;
}

void amg_solve(AMGSOLVER&A,const std::vector<double>&b,std::vector<double>&x,int max_cycles,double tol){
    int n = b.size();
    x.assign(n, 0.0);

    //norm of b 
    double normb=cblas_dnrm2(n,b.data(),1);

    matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    for(int k=0;k<max_cycles;++k){
        A.v_cycle(0, x, b);

        // Compute residual r = b - A1 * x
        std::vector<double> r = b;
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, A.levels[0].A.mkl_handle, descr, x.data(), 1.0, r.data());
        //norm of r
        double normr=cblas_dnrm2(n,r.data(),1);

         if (k % 10 == 0 || normr < tol *normb) {
            std::cout << "Iteration : " << k << " with Residual: " << normr << std::endl;
            
        }

         // Stopping criterion
        if (normr < tol * normb) {
            std::cout << "Converged at cycle " << k << std::endl;
            std::cout << "Residual: " << normr << std::endl;
            break;
        }
    }
}
int main(){

    int n= 60;
    double tol=1e-6;

    std::cout << "Generating Poisson matrix of grid size " << n << "x" << n << "..." << std::endl;
    CSR A=poisson(n);

    AMGSOLVER amg;
    std::cout<<"setting up levels\n";
    amg.amg_setup(A);

    std::vector<double>x(n*n,0.0);
    std::vector<double>b(n*n,1.0);

    amg_solve(amg,b,x,200,tol);
    return 0;
}