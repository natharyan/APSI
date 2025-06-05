#include "../utils/util.hpp"
#include "../hashing/simple_hashing.hpp"
#include "main.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;


uint32_t random_value = static_cast<uint16_t>(rand());

BinPosition binPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + random_value) % nbins,(pos.bin_offset + random_value) % (max_j + 1));
}

BinPosition invertBinPermutation(const BinPosition& pos, uint32_t nbins, uint32_t max_j) {
    return BinPosition((pos.bin_index + nbins - random_value) % nbins, (pos.bin_offset + max_j + 1 - random_value) % (max_j + 1));
}

int main() {
    
    string filename = "sample_sets/emails_bob.txt";
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

    unordered_map<Element, BinPosition> G2;

    for (const auto& elem : T2.table) {
        const Element& element = elem.first;
        const vector<BinPosition>& positions = elem.second;
        if (!positions.empty()) {
            G2[element] = positions[0]; // take the first position as representative
        }
    }

    vector<BinPosition> V;
    // TODO: add retrieval of M1 and m_g_2 from Alice
    unordered_map<mcl::Fp12, BinPosition> M1; // TODO: change
    unordered_map<mcl::Fp12, BinPosition> M2;
    mcl::bn::G2 m_g_2;
    for (const auto& dict_elem : G2) {
        Element elem = dict_elem.first;
        BinPosition pos = dict_elem.second;
        mcl::bn::G1 point_elem = hashToCurve(elem);
        mcl::Fp12 pairing_res;
        pairing(pairing_res, point_elem, m_g_2);
        // check if pairing_res is in M1
        if (M1.find(pairing_res) != M1.end()) {
            V.push_back(pos);
        }
    }

    // output V

    return 0;
}