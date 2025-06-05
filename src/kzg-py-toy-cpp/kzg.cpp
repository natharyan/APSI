#include "kzg.hpp"
#include "util.hpp"
#include "polynomial.hpp" 

//generates trusted setup for kzg commitments
//auto [g1_points, s_g2] = trusted_setup();
//use g1_points for commitment generation
//use s_g2 for verification

std::pair<std::vector<G1>, G2> trusted_setup() {
    initPairing(mcl::BLS12_381); 
    Fr s; //secret scalar for the setup
    s.setByCSPRNG(); //generating random scalar 

    std::vector<G1> g1_points(1024);
    G1 g1;
    mapToG1(g1, 1); 
    g1_points[0] = g1;

    for (size_t i = 1; i < g1_points.size(); ++i) {
        G1::mul(g1_points[i], g1_points[i - 1], s);
    }

    G2 g2;
    mapToG2(g2, 1);
    G2 s_g2;
    G2::mul(s_g2, g2, s);

    return {g1_points, s_g2};
}

G1 commit(const Polynomial& poly, const std::vector<G1>& setup_g1) {
    G1 result;
    result.clear();
    for (size_t i = 0; i < poly.size(); ++i) {
        G1 tmp;
        G1::mul(tmp, setup_g1[i], poly[i]);
        G1::add(result, result, tmp);
    }
    return result;
}

G1 create_proof(const Polynomial& poly, const Point& point, const std::vector<G1>& setup_g1) {
    Polynomial divisor = {Fr(), Fr()};
    divisor[0] = Fr(-point.first);  
    divisor[1] = Fr(1);

    Polynomial numerator = poly;
    Fr::sub(numerator[0], numerator[0], point.second);

    Polynomial quotient = poly_div(numerator, divisor);
    return commit(quotient, setup_g1);
}

bool verify(const G1& C, const G1& pi, const Point& point, const G2& s_g2) {
    G2 g2;
    mapToG2(g2, 1);

    G1 lhs;
    G1 y_g1;
    G1 g1;
    mapToG1(g1, 1); 
    G1::mul(y_g1, g1, point.second);  
    G1::sub(lhs, C, y_g1);

    G2 x_g2;
    Fr x_fr;
    x_fr.setStr(std::to_string(point.first));
    G2::mul(x_g2, g2, x_fr);

    G2 rhs_g2;
    G2::sub(rhs_g2, s_g2, x_g2);

    Fp12 e1, e2;
    pairing(e1, lhs, g2);
    pairing(e2, pi, rhs_g2);
    return e1 == e2;
}