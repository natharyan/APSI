#ifndef SIMPLE_HASHING_H_
#define SIMPLE_HASHING_H_
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>

using Element = std::vector<uint8_t>;

struct ElementHash {
    std::size_t operator()(const Element& elem) const {
        std::size_t seed = elem.size();
        for (auto& b : elem) {
            seed ^= b + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct hs_t {
    uint64_t seed1, seed2, seed3;
    uint32_t hash(const Element& elem, uint32_t nbins, int which) const;
};

// (i, j) pair
struct BinPosition {
    uint32_t bin_index;
    uint32_t bin_offset;
    
    BinPosition(uint32_t bin_idx, uint32_t bin_off) 
        : bin_index(bin_idx), bin_offset(bin_off) {}
    BinPosition() : bin_index(0), bin_offset(0) {}
};

// dictionary: element -> (i, j)
using SimpleHashDict = std::unordered_map<Element, std::vector<BinPosition>, ElementHash>;

struct SimpleHashingContext {
    SimpleHashDict table;
    uint32_t nbins;
    std::shared_ptr<hs_t> hs;
    
    SimpleHashingContext(const SimpleHashDict& dict, uint32_t nb, const std::shared_ptr<hs_t>& hash_struct): table(dict), nbins(nb), hs(hash_struct) {}
};

SimpleHashingContext simple_hashing(
    const std::vector<Element>& elements,
    uint32_t nbins,
    const std::shared_ptr<hs_t>& hs
);

#endif // SIMPLE_HASHING_H_