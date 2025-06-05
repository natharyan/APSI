#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <mcl/bn.hpp>

std::vector<unsigned char> readFileToBytes(const std::string& filename);
void initPairing();
mcl::bn::Fr randomFr();
mcl::bn::G1 hashToCurve(const std::vector<unsigned char>& input);
std::string bytesToString(const std::vector<unsigned char>& bytes);
std::vector<std::string> readEmailsFromFile(const std::string& filename);
std::vector<uint8_t> emailToUint8Array(const std::string& email);