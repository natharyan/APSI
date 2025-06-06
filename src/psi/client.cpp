#include "../utils/util.hpp"
#include "../hashing/simple_hashing.hpp"
#include "main.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

struct Fp12Hash {
    std::size_t operator()(const mcl::Fp12& x) const {
        std::stringstream ss;
        x.save(ss, mcl::IoSerialize);
        return std::hash<std::string>{}(ss.str());
    }
};

struct Fp12Eq {
    bool operator()(const mcl::Fp12& a, const mcl::Fp12& b) const {
        return a == b;
    }
};

uint32_t random_value = []() {
    srand(time(nullptr));  
    return static_cast<uint32_t>(rand());
}();

BinPosition binPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + random_value) % nbins,
                      (pos.bin_offset + random_value) % (max_j + 1));
}

BinPosition invertBinPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + nbins - (random_value % nbins)) % nbins, 
                      (pos.bin_offset + max_j + 1 - (random_value % (max_j + 1))) % (max_j + 1));
}

int main() {
    struct stat st = {0};
    if (stat("output", &st) == -1) {
        mkdir("output", 0700);
    }

    const int rows = 1024;
    const int cols = 256;

    vector<vector<int>> database(rows, vector<int>(cols, 0));

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            database[i][j] = i * cols + j;

    vector<vector<int>> M1_matrix(rows, vector<int>(cols));
    auto start_m1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M1_matrix[i][j] = database[i][j] * 2;
    auto end_m1 = chrono::high_resolution_clock::now();
    long long m1_time = chrono::duration_cast<chrono::microseconds>(end_m1 - start_m1).count();
    cout << "Time to create M1: " << m1_time << " μs" << endl;

    vector<vector<int>> M2(rows, vector<int>(cols));
    auto start_m2 = chrono::high_resolution_clock::now();
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            M2[i][j] = M1_matrix[i][j] + 1; 
    auto end_m2 = chrono::high_resolution_clock::now();
    long long m2_time = chrono::duration_cast<chrono::microseconds>(end_m2 - start_m2).count();
    cout << "Time to create M2: " << m2_time << " μs" << endl;

    ofstream matrix_outfile("output/matrix_db.txt");
    if (matrix_outfile.is_open()) {
        matrix_outfile << "M1[0][0]: " << M1_matrix[0][0] << endl;
        matrix_outfile << "M2[0][0]: " << M2[0][0] << endl;
        matrix_outfile.close();
    }

    initPairing(); 
    
    string filename = "sample_sets/emails_client.txt";
    vector<string> emails = readEmailsFromFile(filename);
    cout << "Successfully read " << emails.size() << " emails from " << filename << endl;
    
    // get the simple_hashing table
    uint32_t nbins = 1000;
    auto hs = make_shared<hs_t>(hs_t{123456789, 987654321, 1122334455});
    vector<Element> elements;
    elements.reserve(emails.size());  
    
    for (const auto& email : emails) {
        elements.push_back(emailToUint8Array(email));
    }
    SimpleHashingContext T1 = simple_hashing(elements, nbins, hs);

    mcl::bn::Fr m = randomFr(); 

    mcl::bn::G2 g2;
    mcl::bn::mapToG2(g2, 1);
    mcl::bn::G2 m_g_2;
    mcl::bn::G2::mul(m_g_2, g2, m);

    const int repeat_count = 1;
    unordered_map<Element, BinPosition, ElementHash> G1; 

    long long total_g1_time = 0;
    for (int rep = 0; rep < repeat_count; ++rep) {
        unordered_map<Element, BinPosition, ElementHash> G1_tmp;
        G1_tmp.reserve(T1.table.size());  
        
        auto start_g1 = chrono::high_resolution_clock::now();
        for (const auto& elem : T1.table) {
            const Element& element = elem.first;
            const vector<BinPosition>& positions = elem.second;
            if (!positions.empty()) {
                G1_tmp[element] = positions[0];
            }
        }
        auto end_g1 = chrono::high_resolution_clock::now();
        total_g1_time += chrono::duration_cast<chrono::microseconds>(end_g1 - start_g1).count();
        if (rep == 0) {
            G1 = move(G1_tmp);  
        }
    }
    cout << "Average time to create G1: " << (total_g1_time / repeat_count) << " μs" << endl;

    uint32_t max_j = 0;
    for (const auto& [_, pos] : G1) {
        max_j = max(max_j, pos.bin_offset);
    }
    unordered_map<Element, mcl::bn::G1, ElementHash> curve_cache;
    curve_cache.reserve(G1.size());
    
    auto start_hash = chrono::high_resolution_clock::now();
    for (const auto& [elem, _] : G1) {
        curve_cache[elem] = hashToCurve(elem);
    }
    auto end_hash = chrono::high_resolution_clock::now();
    long long hash_time = chrono::duration_cast<chrono::microseconds>(end_hash - start_hash).count();
    cout << "Time to pre-compute curve hashes: " << hash_time << " μs" << endl;

    unordered_map<mcl::Fp12, BinPosition, Fp12Hash, Fp12Eq> M1_map; 
    long long total_m1_time = 0;
    
    for (int rep = 0; rep < repeat_count; ++rep) {
        unordered_map<mcl::Fp12, BinPosition, Fp12Hash, Fp12Eq> M1_tmp;
        M1_tmp.reserve(G1.size());  
        
        auto start_m1 = chrono::high_resolution_clock::now();
        for (const auto& [elem, pos] : G1) {
            const mcl::bn::G1& point_elem = curve_cache[elem];
            
            mcl::Fp12 pairing_res;
            pairing(pairing_res, point_elem, m_g_2);
            M1_tmp[pairing_res] = binPermutation(pos, nbins, max_j);
        }
        auto end_m1 = chrono::high_resolution_clock::now();
        total_m1_time += chrono::duration_cast<chrono::microseconds>(end_m1 - start_m1).count();
        
        if (rep == 0) {
            M1_map = move(M1_tmp);  
        }
    }
    cout << "Average time to create M1 (pairings only): " << (total_m1_time / repeat_count) << " μs" << endl;
    cout << "Total M1 time (including pre-computation): " << (total_m1_time / repeat_count + hash_time) << " μs" << endl;

    ofstream outfile("output/client_data.txt");
    if (!outfile.is_open()) {
        cerr << "Failed to open output file!" << endl;
        return 1;
    }

    stringstream ss_g2;
    m_g_2.save(ss_g2, mcl::IoSerialize);
    outfile << "m_g_2: " << ss_g2.str() << endl;

    outfile << "M1:" << endl;
    for (const auto& entry : M1_map) {
        const auto& pair_val = entry.first;
        const auto& pos = entry.second;
        stringstream ss_fp12;
        pair_val.save(ss_fp12, mcl::IoSerialize);
        outfile << ss_fp12.str() << " " << pos.bin_index << " " << pos.bin_offset << endl;
    }

    outfile.close();
    cout << "Saved m_g_2 and M1 to output/client_data.txt" << endl;

    return 0;
}