#include "simple_hashing.hpp"
#include <algorithm>
#include <iostream>

using namespace std;

uint32_t hs_t::hash(const Element& elem, uint32_t nbins, int which) const {
    uint64_t seed = (which == 0) ? seed1 : (which == 1) ? seed2 : seed3;
    uint32_t h = static_cast<uint32_t>(seed);
    for (vector<uint8_t>::const_iterator it = elem.begin(); it != elem.end(); ++it) {
        h = h * 31 + (*it);
    }
    return h % nbins;
}

SimpleHashingContext simple_hashing(
    const vector<Element>& elements,
    uint32_t nbins,
    const shared_ptr<hs_t>& hs
) {
    // bin elements
    vector<vector<Element> > bins(nbins);
    for (vector<Element>::const_iterator elem_it = elements.begin(); elem_it != elements.end(); ++elem_it) {
        for (int which = 0; which < 3; ++which) {
            uint32_t i = hs->hash(*elem_it, nbins, which);
            bins[i].push_back(*elem_it);
        }
    }
    
    // order each bin
    SimpleHashDict dict;
    for (uint32_t i = 0; i < nbins; ++i) {
        std::vector<Element>& bin = bins[i];
        std::sort(bin.begin(), bin.end());
        for (uint32_t j = 0; j < bin.size(); ++j) {
            dict[bin[j]].push_back(BinPosition(i, j)); // dict[element] -> 3 BinPosition(i, j) of element
        }
    }
    
    return SimpleHashingContext(dict, nbins, hs);
}