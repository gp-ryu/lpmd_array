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

using namespace Eigen;
using namespace std;

typedef SparseMatrix<int8_t> SparseMatrix8i;
typedef Triplet<double> T;

struct ColumnStats {
    double mean;
    double median;
    double iqr;
    double std_dev;
};


int binarySearch(vector<int> arr, int target) {
    int left = 0;
    int right = arr.size() - 1;
    int mid;

    while (left <= right) {
        mid = (left + right) / 2;

        if (arr[mid] == target) {
            return mid;
        } else if (arr[mid] > target) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }

    // If target is not found in arr
    return -1;
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

map<string, ColumnStats> calculateColumnStats(const MatrixXd& matrix, vector<string> columnNames){
    map<string, ColumnStats> res;
    int numCols = matrix.cols(); // get number of columns
    ColumnStats stats; // create struct to hold stats

    for (int i = 0; i < numCols; i++) {
        VectorXd col = matrix.col(i); // get current column
        std::sort(col.data(), col.data() + col.size()); // sort column data
        int mid = col.size() / 2; // calculate middle index

        // calculate mean
        stats.mean = col.mean();

        // calculate median
        stats.median = col.size() % 2 == 0 ? (col(mid - 1) + col(mid)) / 2 : col(mid);

        // calculate interquartile range
        int q1_idx = col.size() / 4;
        int q3_idx = 3 * col.size() / 4;
        double q1 = col(q1_idx);
        double q3 = col(q3_idx);
        stats.iqr = q3 - q1;

        // calculate standard deviation
        double variance = ((col.array() - stats.mean).square().sum()) / (col.size() - 1);
        stats.std_dev = std::sqrt(variance);
        res.insert(make_pair(columnNames[i+1], stats));
    }


    return res; // return struct of column stats
}

void usage(){
    cout << "\nUsage:  lpmd_calc <TOE_lv3.txt> <lpmd_dist_map> <output_prifix>\n\n"
        << "TOE_lv3.txt\t:  Methylation array level3 file\n"
        << "lpmd_dist_map\t:  LPMD distance mapping file consisting of {CHR, CGID1, CGID2, DISTANCE}\n"
        << "output_prefix\t:  Prefix for <prefix>_whole_list.txt, <prefix>_stat.txt\n" << endl;
}

int main(int argc, char* argv[]) {
    if(argc<4) {usage(); return 1;}

    ifstream file(argv[1]);
    string line; 

    // <columnNames> allocation
    vector<string> columnNames;
    getline(file, line);
    istringstream iss1(line);
    string column;
    while (getline(iss1, column, '\t')){
        columnNames.push_back(column);
    }
    

    // <rownames , fields> read and allocation
    unordered_map<string, int> rownames;
    vector<vector<double>> fields;
    int jj = 0;
    cerr << "Read:\t" << argv[1] << endl;
    while(getline(file, line)){
        string line_;
        istringstream iss1(line);
        getline(iss1, line_, '\t');
        rownames.emplace(line_,jj);

        string record;
        vector<double> lines;
        int j=0;
        while(getline(iss1, record, '\t')){
            lines.push_back(stod(record));
            j++;
        }
        fields.push_back(lines);
        jj++;
    }


    // <ilmnMap_mat> *ptr allocation
    cerr << "Build:\t" << "Matrix of " << argv[1] << " : " << fields.size() << " x " << columnNames.size()-1 << endl;
    MatrixXd ilmnMap_mat(fields.size(), columnNames.size()-1);
    for(int i=0; i < fields.size(); i++){
        for(int j=1; j < columnNames.size(); j++){
            double *ptr_field = &fields[i][j-1];
            ilmnMap_mat(i, j-1) = *ptr_field;
        }
    }
    file.close();


    cerr << "Build:\tdistMap into distVec" << endl;
    // Reading distMap into distVec
    vector<vector<string>> distVec;
    vector<string> chr_vec, cgid1_vec, cgid2_vec, distance_vec;
    ifstream distFile(argv[2]);
    while(getline(distFile, line)){
        istringstream iss(line);
        string chr, distance, cgid1, cgid2;
        getline(iss, chr, '\t');
        getline(iss, cgid1, '\t');
        getline(iss, cgid2, '\t');
        getline(iss, distance, '\t');

        distVec.push_back({chr, cgid1, cgid2, distance});
        chr_vec.push_back(chr);
        cgid1_vec.push_back(cgid1);
        cgid2_vec.push_back(cgid2);
        distance_vec.push_back(distance);
    }
    distFile.close();

    cerr << "Calc:\tSparse substraction matrix" << endl;

    // Making sparse substraction matrix

    vector<vector<int>> A_index;
    for (int i=0; i < distVec.size(); i++){
        string& cgid_this = distVec[i][1];
        string& cgid_that = distVec[i][2];
        string& cgid_dist = distVec[i][3];

        auto index_this_it = rownames.find(cgid_this);
        auto index_that_it = rownames.find(cgid_that);
        if(index_that_it != rownames.end() && index_this_it != rownames.end()){
            int index_this = rownames[cgid_this];
            int index_that = rownames[cgid_that];
            A_index.push_back({i, index_this, index_that, stoi(cgid_dist)});
        }
    }

    vector<Triplet<double>> triplets;
    for (int i = 0; i < A_index.size(); i++) {
        triplets.push_back(Triplet<double>(i, A_index[i][1], 1.0));
        triplets.push_back(Triplet<double>(i, A_index[i][2], -1.0));
    }

    SparseMatrix<double> A(A_index.size(), rownames.size());
    A.setFromTriplets(triplets.begin(), triplets.end());

    SparseMatrix<double, ColMajor> A_ccs = A;

    A.makeCompressed();


    // substraction matrix build
    MatrixXd lpmd(A.rows(), columnNames.size()-1);
    lpmd = A * ilmnMap_mat;


    cerr << "Write:\twhole LPMD result" << endl;
    // writing whole lpmd output
    string lpmd_res = argv[3];
    lpmd_res += "_whole_list.txt";
    lpmd_header(lpmd_res, columnNames);

    ofstream(output);
    output.open(lpmd_res);
    output << "chrom\tcgids\tpos_dist";
    for(int i=1; i < columnNames.size(); i++){
        output << '\t' << columnNames[i];
    }
    output << endl; 
    for(int i=0; i < A_index.size(); i++){
        output << distVec[i][0] << '\t' << distVec[i][1] + ":" + distVec[i][2] << '\t' << distVec[i][3];
        for(int j=0; j < lpmd.cols();j++){
            output << '\t' << lpmd(i,j);
        }
        output << endl;
    }
    output.close();



    cerr << "Calc:\tCumulative lpmd output" << endl;
    // cumulative lpmd score Matrix build
    vector<Triplet<double>> triplets_cum_lpmd;
    vector<int> t_cum_lpmd;
    for(int i =0; i<100; i++){
        t_cum_lpmd.push_back(0);
    }
    for(int i=0; i < A_index.size(); i++){
        triplets_cum_lpmd.push_back(Triplet<double>(A_index[i][3], i, 1.0));
        t_cum_lpmd[A_index[i][3]-1]++;
    }
    SparseMatrix<double> B(100,A_index.size());
    B.setFromTriplets(triplets_cum_lpmd.begin(), triplets_cum_lpmd.end());
    B.makeCompressed();

    // Matrix multiplication
    MatrixXd cum_lpmd(100, columnNames.size()-1);
    MatrixXd cum_lpmd_mean(100, columnNames.size()-1);
    cum_lpmd = B * lpmd.cwiseAbs();

    for(int i = 0; i < 100; i++){
        cum_lpmd_mean.row(i) = cum_lpmd.row(i) / t_cum_lpmd[i];

    }

    cerr << "Write:\tCumulative lpmd" << endl;
    string cum_lpmd_res = argv[3];
    cum_lpmd_res += "_cumul.txt";
    ofstream(output_cum_lpmd_res);
    output_cum_lpmd_res.open(cum_lpmd_res); 
    output_cum_lpmd_res << "Dist";
    for(int i=1; i < columnNames.size(); i++){
        output_cum_lpmd_res << '\t' << columnNames[i];
    }
    output_cum_lpmd_res << endl;
    for(int i=0; i < 100; i++){
        output_cum_lpmd_res << i+1;
        for(int j=0; j < cum_lpmd_mean.cols(); j++){
            output_cum_lpmd_res << '\t' << cum_lpmd_mean(i,j);
        }
        output_cum_lpmd_res << endl;
    }
    output_cum_lpmd_res.close();




    cerr << "Write:\tLPMD stats result" << endl;
    string lpmd_stats = argv[3];
    lpmd_stats += "_stats.txt";
    ofstream(output_stats);
    output_stats.open(lpmd_stats);
    output_stats << "SampleID\tMean\tMedian\tIQR\tsd" << endl;
    auto stats = calculateColumnStats(lpmd.cwiseAbs(), columnNames);
    for(const auto& var : stats){
        output_stats << var.first << '\t' << var.second.mean << '\t' << var.second.median << '\t' <<  var.second.iqr << '\t' << var.second.std_dev << endl;
    }
    output_stats.close();

    return 0;
}

