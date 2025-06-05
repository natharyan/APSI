#include "util.hpp"
#include <random>
#include <fstream>
#include <vector>

using namespace std; 

unordered_map<string, mcl::bn::G1> hashOracle;

vector<unsigned char> readFileToBytes(const string& filename) {
    ifstream file(filename, ios::binary);
    return vector<unsigned char>(istreambuf_iterator<char>(file), {});
}

// initializing the MCL pairing library with the BLS12-381 curve
void initPairing() {
    mcl::bn::initPairing(mcl::BLS12_381);
}

mcl::bn::Fr randomFr() {
    mcl::bn::Fr x;
    x.setByCSPRNG();
    return x;
}

string bytesToString(const vector<unsigned char>& bytes) {
    return string(bytes.begin(), bytes.end());
}

mcl::bn::G1 hashToCurve(const vector<unsigned char>& input) {
    
    string key = bytesToString(input);
    auto it = hashOracle.find(key);
    if (it != hashOracle.end()) {
        return it->second;
    }
    
    mcl::bn::Fr s = randomFr(); 
    mcl::bn::G1 generator;
    mcl::bn::mapToG1(generator, 1);
    mcl::bn::G1 point;
    mcl::bn::G1::mul(point, generator, s);  // s * G
    
    hashOracle[key] = point;
    
    return point;
}

vector<string> readEmailsFromFile(const string& filename) {
    vector<string> emails;
    
    vector<unsigned char> fileBytes = readFileToBytes(filename);
    if (fileBytes.empty()) {
        cerr << "Error: Could not read file or file is empty: " << filename << endl;
        return emails;
    }
    
    string fileContent(fileBytes.begin(), fileBytes.end());
    stringstream ss(fileContent);
    string email;
    while (getline(ss, email)) {
        if (!email.empty()) {
            emails.push_back(email);
        }
    }
    
    return emails;
}

vector<uint8_t> emailToUint8Array(const string& email) {
    return vector<uint8_t>(email.begin(), email.end());
}