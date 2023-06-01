#include <iostream>
#include <string>
#include <vector> 
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdint.h>
#include <thread>
#include <stack>
#include <ctime>


using namespace Eigen;
using namespace std;

std::stack<clock_t> tictoc_stack;

void tic() {
    tictoc_stack.push(clock());
}

void toc() {
    std::cout << "Time elapsed: "
              << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    tictoc_stack.pop();
}


typedef SparseMatrix<int8_t> SparseMatrix8i;
typedef Triplet<double> T;


void usage(){
    cout << "\nUsage:  lpmd_calc <TOE_lv3.txt> <lpmd_dist_map> <output_prifix>\n\n"
        << "TOE_lv3.txt\t:  Methylation array level3 file\n"
        << "lpmd_dist_map\t:  LPMD distance mapping file consisting of {CHR, CGID1, CGID2, DISTANCE}\n"
        << "output_prefix\t:  Prefix for <prefix>_whole_list.txt, <prefix>_stat.txt\n" << endl;
}

void lpmd_header(string file_name, vector<string> header_list){
    ofstream output;
    output.open(file_name);

    output << "chr\tcgid1\tcgid2\tdistance";
    for(int i=1; i < header_list.size();i++){
        output << '\t' << header_list[i];
    }
    output << endl;
}



vector<string> get_columnNames(ifstream& file_stream){
    vector<string> return_vector;
    string line; 
    getline(file_stream, line);
    istringstream iss(line);
    string column; 
    while(getline(iss, column, '\t')){
        return_vector.push_back(column);
    }
    return return_vector;
}





int main(int argc, char* argv[]) {
    if(argc<4) {usage(); return 1;}

    tic();
    /* Reads Methylation file into {met_colnames, met_rownames, met_fields} */; cerr << "Read:\t" << argv[1];
    vector<string> met_colnames, met_rownames;
    unordered_map<string, int> met_rownames_index_M;
    vector<vector<double>> met_fields;

    ifstream file(argv[1]);                                             //open Methylation file
    met_colnames = get_columnNames(file);                               //get column names {met_colnames}
    string line;     
    int row_counter = 0;
    while(getline(file, line)){
        int col_counter = 0;
        vector<double> fields_; 
        istringstream iss1(line);
        string line_;
        while(getline(iss1, line_, '\t')){                          // allocate first column into {met_rownames_index_M, met_rownames}
            if(col_counter % met_colnames.size() == 0){
                met_rownames_index_M.emplace(line_, row_counter);
                met_rownames.push_back(line_);
            }
            else{
                fields_.push_back(stod(line_));
            }
            col_counter++;
        }
        met_fields.push_back(fields_);
        row_counter++;
    }
    file.close();
    cout << "\t\t\t\t\t";
    toc();



    tic();
    /*  <ilmnMap_mat> *ptr allocation */; cerr << "Build:\t" << "Matrix of " << argv[1] << " : " << met_fields.size() << " x " << met_colnames.size()-1;
    MatrixXd ilmnMap_mat(met_fields.size(), met_colnames.size()-1);
    for(int i=0; i < met_fields.size(); i++){
        for(int j=1; j < met_colnames.size(); j++){
            double *ptr_field = &met_fields[i][j-1];
            ilmnMap_mat(i, j-1) = *ptr_field;
        }
    }
    cout << "\t\t";
    toc();


    tic(); 
    /* Reading distMap into map_dist */; cerr << "Build:\tdistMap into distVec" ;
    vector<vector<string>> distVec, map_dist;
    vector<string> chr_vec, cgid1_vec, cgid2_vec;
    unordered_map<string, int> map_chr_M, map_cgid1_M, map_cgid2_M, map_dist_M;

    ifstream distFile(argv[2]);
    int row_index = 0;
    while(getline(distFile, line)){
        istringstream iss(line);
        string chr, distance, cgid1, cgid2;
        getline(iss, chr, '\t');
        getline(iss, cgid1, '\t');
        getline(iss, cgid2, '\t');
        getline(iss, distance, '\t');

        distVec.push_back({chr, cgid1, cgid2, distance});
        map_cgid1_M.emplace(cgid1, row_index);
        map_cgid2_M.emplace(cgid2, row_index);
        row_index++;
    }
    distFile.close();
    cout << "\t\t\t\t\t\t\t\t\t";
    toc();
    
    

    tic();
    /* Parsing bi-existance pairs from mat_rownames with map_dist in {pairIdx_dist_M}*/; cerr << "Build:\tSparse substraction matrix" ;
    vector<vector<int>> met_pairIdxV;
    map<int, vector<vector<string>>> met_pairIdx_dist_M;
    map<int, vector<T>> triplets;
    for (int i=0; i < distVec.size(); i++){
        string& chr = distVec[i][0];
        string& cgid_this = distVec[i][1];
        string& cgid_that = distVec[i][2];
        int cgid_dist = stoi(distVec[i][3]);
        int index_cgid1 = met_rownames_index_M[cgid_this];
        int index_cgid2 = met_rownames_index_M[cgid_that];

        if(index_cgid1 != 0 && index_cgid2 != 0){
            met_pairIdxV.push_back({index_cgid1, index_cgid2, cgid_dist});
            vector<string> temp_met_pairIdx_dist_M = {chr, to_string(index_cgid1), to_string(index_cgid2)};
            met_pairIdx_dist_M[cgid_dist].push_back(temp_met_pairIdx_dist_M);
        }
    }

    for(int dist = 2; dist < 41; dist++){
        for(int row = 0; row < met_pairIdx_dist_M[dist].size(); row++){
            triplets[dist].push_back(T(row, stoi(met_pairIdx_dist_M[dist][row][1]), 1.0));
            triplets[dist].push_back(T(row, stoi(met_pairIdx_dist_M[dist][row][2]), -1.0));
        }
    }
    SparseMatrix8i pair_matrix(triplets[2].size()/2, met_rownames.size());
    pair_matrix.setFromTriplets(triplets[2].begin(), triplets[2].end());
    cout << "\t\t\t\t\t\t\t\t";
    toc();



    /* substraction matrix build */; cerr << "Calc:\tMatrix multiplication by distance" << endl;
    unordered_map<int, MatrixXd> lpmds;
    unordered_map<int, SparseMatrix<double>> pair_matrixM;
    for(int dist = 2; dist < 41; dist++){
        //lpmds.emplace_back(dist, met_colnames.size());
        //lpmds.back().setConstant(dist);
        SparseMatrix<double> pair_matrix(triplets[dist].size()/2, met_rownames.size());
        pair_matrix.setFromTriplets(triplets[dist].begin(), triplets[dist].end());
        pair_matrixM.emplace(dist, pair_matrix);
        cout << "pair_matrixM[" << dist << "]: [" << pair_matrixM[dist].rows()  << 'x' << pair_matrixM[dist].cols() << ']' << endl;        
        lpmds[dist] = pair_matrixM[dist] * ilmnMap_mat;
    }



    cerr << "Write:\twhole LPMD result" << endl;
    // writing whole lpmd output
    string lpmd_res = argv[3];
    lpmd_res += "_whole_list.txt";

    ofstream(output);
    output.open(lpmd_res);
    output << "chr\tcgids\tdist";
    for(int i=1; i < met_colnames.size(); i++){   // writing first column names
        output << '\t' << met_colnames[i];
    } output << endl; 

    for(int dist = 2; dist < 41; dist++){
        for(int row = 0; row < lpmds[dist].rows(); row++){
            output << met_pairIdx_dist_M[dist][row][0] 
                << '\t' << met_rownames[stoi(met_pairIdx_dist_M[dist][row][1])]
                << ':' << met_rownames[stoi(met_pairIdx_dist_M[dist][row][2])] 
                << '\t' << dist;


            for(int col = 0; col < lpmds[dist].cols(); col++){
                output << '\t' << lpmds[dist](row,col);
            }
            output << endl;
        }
    }
    output.close();





    cerr << "Write:\tCumulative lpmd" << endl;
    string cum_lpmd_res = argv[3];
    cum_lpmd_res += "_cumul.txt";
    ofstream(output_cum_lpmd_res);
    output_cum_lpmd_res.open(cum_lpmd_res); 
    output_cum_lpmd_res << "Dist";
    for(int i=1; i < met_colnames.size(); i++){
        output_cum_lpmd_res << '\t' << met_colnames[i];
    }
    output_cum_lpmd_res << endl;
    for(int i=2; i < 41; i++){
        output_cum_lpmd_res << i;
        for(int j = 0; j < met_colnames.size()-1;j++){
            double mean_val = lpmds[i].col(j).cwiseAbs().mean();
            output_cum_lpmd_res << '\t' << mean_val;
        }
        output_cum_lpmd_res << endl;
    }
    output_cum_lpmd_res.close();



// 
//     cerr << "Write:\tLPMD stats result" << endl;
//     string lpmd_stats = argv[3];
//     lpmd_stats += "_stats.txt";
//     ofstream(output_stats);
//     output_stats.open(lpmd_stats);
//     output_stats << "SampleID\tMean\tMedian\tIQR\tsd" << endl;
//     auto stats = calculateColumnStats(lpmd.cwiseAbs(), met_colnames);
//     for(const auto& var : stats){
//         output_stats << var.first << '\t' << var.second.mean << '\t' << var.second.median << '\t' <<  var.second.iqr << '\t' << var.second.std_dev << endl;
//     }
//     output_stats.close();
    return 0;
}


