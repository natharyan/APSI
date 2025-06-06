#include "../utils/util.hpp"
#include "../hashing/simple_hashing.hpp"
#include "main.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream> // Added for ofstream

using namespace std;


uint32_t random_value = static_cast<uint16_t>(rand());

BinPosition binPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + random_value) % nbins,(pos.bin_offset + random_value) % (max_j + 1));
}

BinPosition invertBinPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + nbins - random_value) % nbins, (pos.bin_offset + max_j + 1 - random_value) % (max_j + 1));
}

int main() {
    
    string filename = "sample_sets/emails_text.txt";
    vector<string> emails = readEmailsFromFile(filename);
    cout << "Successfully read " << emails.size() << " emails from " << filename << endl;
    
    // get the simple_hashing table
    uint32_t nbins = 1000;
    auto hs = make_shared<hs_t>(hs_t{123456789, 987654321, 1122334455});
    vector<Element> elements;
    for (const auto& email : emails) {
        elements.push_back(emailToUint8Array(email));
    }
    SimpleHashingContext T1 = simple_hashing(elements, nbins, hs);

    mcl::bn::Fr m = randomFr(); 

    unordered_map<Element, BinPosition> G1;

    for (const auto& elem : T1.table) {
        const Element& element = elem.first;
        const vector<BinPosition>& positions = elem.second;
        if (!positions.empty()) {
            G1[element] = positions[0]; // take the first position as representative
        }
    }

    // in max_j store the max bin_offset from G1
    uint32_t max_j = 0;
    for (const auto& [_, pos] : G1) {
        max_j = max(max_j, pos.bin_offset);
    }

    unordered_map<mcl::Fp12, BinPosition> M1;
    mcl::bn::G2 g2;
    mcl::bn::mapToG2(g2, 1);
    mcl::bn::G2 m_g_2;
    mcl::bn::G2::mul(m_g_2, g2, m);
    for (const auto& dict_elem : G1) {
        Element elem = dict_elem.first;
        BinPosition pos = dict_elem.second;
        mcl::bn::G1 point_elem = hashToCurve(elem);
        mcl::Fp12 pairing_res;
        pairing(pairing_res, point_elem, m_g_2);
        M1[pairing_res] = binPermutation(pos, nbins, max_j);
    }

    //  sending m_g_2, M1 to Bob.
    ofstream outfile("output/alice_data.txt");
    if (!outfile.is_open()) {
        cerr << "Failed to open output file!" << endl;
        return 1;
    }

    stringstream ss_g2;
    m_g_2.save(ss_g2, mcl::IoSerialize);
    outfile << "m_g_2: " << ss_g2.str() << endl;

    outfile << "M1:" << endl;
    for (const auto& [pair_val, pos] : M1) {
        stringstream ss_fp12;
        pair_val.save(ss_fp12, mcl::IoSerialize);
        outfile << ss_fp12.str() << " " << pos.bin_index << " " << pos.bin_offset << endl;
    }

    outfile.close();
    cout << "Saved m_g_2 and M1 to output/alice_data.txt" << endl;

    return 0;
}