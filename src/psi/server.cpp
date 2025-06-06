#include "../utils/util.hpp"
#include "../hashing/simple_hashing.hpp"
#include "main.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono> 

using namespace std;

uint32_t random_value = static_cast<uint16_t>(rand());

BinPosition binPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + random_value) % nbins,(pos.bin_offset + random_value) % (max_j + 1));
}

BinPosition invertBinPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + nbins - random_value) % nbins, (pos.bin_offset + max_j + 1 - random_value) % (max_j + 1));
}

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

int main() {
    initPairing();

    string filename = "sample_sets/emails_server.txt";
    vector<string> emails = readEmailsFromFile(filename);
    cout << "Successfully read " << emails.size() << " emails from " << filename << endl;
    
    // get the simple_hashing table
    uint32_t nbins = 1000;
    auto hs = make_shared<hs_t>(hs_t{123456789, 987654321, 1122334455});
    vector<Element> elements;
    for (const auto& email : emails) {
        elements.push_back(emailToUint8Array(email));
    }
    SimpleHashingContext T2 = simple_hashing(elements, nbins, hs);

    const int repeat_count = 10;
    unordered_map<Element, BinPosition, ElementHash> G2;
    long long total_g2_time = 0;
    for (int rep = 0; rep < repeat_count; ++rep) {
        unordered_map<Element, BinPosition, ElementHash> G2_tmp;
        auto start_g2 = chrono::high_resolution_clock::now();
        for (const auto& elem : T2.table) {
            const Element& element = elem.first;
            const vector<BinPosition>& positions = elem.second;
            if (!positions.empty()) {
                G2_tmp[element] = positions[0]; // take the first position as representative
            }
        }
        auto end_g2 = chrono::high_resolution_clock::now();
        total_g2_time += chrono::duration_cast<chrono::microseconds>(end_g2 - start_g2).count();
        if (rep == 0) {
            G2 = G2_tmp;
        }
    }
    cout << "Average time to create G2: " << (total_g2_time / repeat_count) << " μs" << endl;

    // in max_j store the max bin_offset from G2
    uint32_t max_j = 0;
    for (const auto& [_, pos] : G2) {
        max_j = max(max_j, pos.bin_offset);
    }

    vector<BinPosition> V;
    unordered_map<mcl::Fp12, BinPosition, Fp12Hash, Fp12Eq> M1;
    mcl::bn::G2 m_g_2;

    ifstream infile("output/client_data.txt");
    if (!infile.is_open()) {
        cerr << "Failed to open Client's data file!" << endl;
        return 1;
    }

    string line;
    getline(infile, line);
    if (line.rfind("m_g_2: ", 0) == 0) {
        string g2_str = line.substr(8);
        stringstream ss_g2(g2_str);
        m_g_2.load(ss_g2, mcl::IoSerialize);
    } else {
        cerr << "Invalid file format: m_g_2 not found" << endl;
        return 1;
    }

    // Read the M1 header
    getline(infile, line);
    if (line != "M1:") {
        cerr << "Invalid file format: M1 header not found" << endl;
        return 1;
    }

    // Read all entries of M1
    while (getline(infile, line)) {
        size_t tab1 = line.find('\t');
        size_t tab2 = line.find('\t', tab1 + 1);
        if (tab1 == string::npos || tab2 == string::npos) continue;
        string fp12_str = line.substr(0, tab1);
        uint32_t bin_index = stoi(line.substr(tab1 + 1, tab2 - tab1 - 1));
        uint32_t bin_offset = stoi(line.substr(tab2 + 1));

        mcl::Fp12 val;
        stringstream ss_fp12(fp12_str);
        val.load(ss_fp12, mcl::IoSerialize);
        M1[val] = BinPosition(bin_index, bin_offset);
    }
    infile.close();

    unordered_map<mcl::Fp12, BinPosition, Fp12Hash, Fp12Eq> M2;
    long long total_m2_time = 0;
    for (int rep = 0; rep < repeat_count; ++rep) {
        unordered_map<mcl::Fp12, BinPosition, Fp12Hash, Fp12Eq> M2_tmp;
        auto start_m2 = chrono::high_resolution_clock::now();
        for (const auto& dict_elem : G2) {
            Element elem = dict_elem.first;
            BinPosition pos = dict_elem.second;
            mcl::bn::G1 point_elem = hashToCurve(elem);
            mcl::Fp12 pairing_res;
            pairing(pairing_res, point_elem, m_g_2);
            M2_tmp[pairing_res] = binPermutation(pos, nbins, max_j);
        }
        auto end_m2 = chrono::high_resolution_clock::now();
        total_m2_time += chrono::duration_cast<chrono::microseconds>(end_m2 - start_m2).count();
        if (rep == 0) {
            M2 = M2_tmp;
        }
    }
    cout << "Average time to create M2 (pairings): " << (total_m2_time / repeat_count) << " μs" << endl;

    for (const auto& [pair_val, pos] : M2) {
        if (M1.find(pair_val) != M1.end()) {
            V.push_back(pos);
        }
    }

    cout << "Intersection size: " << V.size() << endl;

    ofstream outfile("output/intersection_result.txt");
    if (!outfile.is_open()) {
        cerr << "Failed to open output file!" << endl;
        return 1;
    }

    for (const auto& [pair_val, pos] : M1) {
        stringstream ss_fp12;
        pair_val.save(ss_fp12, mcl::IoSerialize);
        outfile << ss_fp12.str() << "\t" << pos.bin_index << "\t" << pos.bin_offset << endl;
    }

    outfile.close();

    return 0;
}