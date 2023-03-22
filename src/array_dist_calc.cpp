#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <omp.h>

using namespace std;

// 
unordered_map <string, vector<pair<string, int>>> read_array_pos_map_file(string file_name){
   ifstream file(file_name);

   string line;
   unordered_map <string, vector<pair<string, int>>> cgID_info;

   while(getline(file, line)){
       istringstream iss(line);
       string chr, cgID;
       double pos; 

       getline(iss, chr, '\t');
       getline(iss, cgID, '\t');
       iss >> pos;
       cgID_info[chr].push_back(make_pair(cgID, pos));

   }
   return cgID_info;
}

// Function to calculate the distance between two coordinates
double distance(double x1, double x2) {
    return abs(x1 - x2);
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        cerr << "Please provide the input file name as a command-line argumant." << endl;
        return 1;
    }
    ifstream file(argv[1]);
    if ( !file.is_open()) {
        cerr << "Could not open file " << argv[1] << endl;
        return 1;
    }
    string line;
    unordered_map<string, vector<pair<string, int>>> groups;
    groups = read_array_pos_map_file(argv[1]);

    // Read the data from the TSV file and group by the first column

    // Loop through each group and calculate the distance between each pair of records
    #pragma omp parallel for schedule(dynamic)
    for (auto& group : groups) {
        vector<pair<string, int>>& records = group.second;
        for (int i = 0; i < records.size(); i++) {
            for (int j = i + 1; j < records.size(); j++) {
                double dist = distance(records[i].second, records[j].second);
                if (dist < 100) {
                    cout << group.first << '\t' << records[i].first << '\t' << records[j].first << '\t' << dist << endl;
                }
            }
        }
    }
    return 0;
}

