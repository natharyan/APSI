#include "util.hpp"
#include <random>

//initializing the MCL pairing library with the BLS12-381 curve
void initPairing() {
    mcl::bn::initPairing(mcl::BLS12_381);
}

mcl::bn::Fr randomFr() {
    mcl::bn::Fr x;
    x.setByCSPRNG();
    return x;
}