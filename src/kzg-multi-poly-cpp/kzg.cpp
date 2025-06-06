#include "kzg.hpp"
#include "util.hpp"
#include "polynomial.hpp"
#include <random> 

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
    for (size_t i = 0; i < poly.size() && i < setup_g1.size(); ++i) {
        G1 tmp;
        G1::mul(tmp, setup_g1[i], poly[i]);
        G1::add(result, result, tmp);
    }
    return result;
}

//multi-commitment for a vector of polynomials
std::vector<G1> multi_commit(const MultiPolynomial& polys, const std::vector<G1>& setup_g1) {
    std::vector<G1> commitments;
    commitments.reserve(polys.size());
    
    for (const auto& poly : polys) {
        commitments.push_back(commit(poly, setup_g1));
    }
    
    return commitments;
}

MultiProof create_multi_proof(const MultiPolynomial& polys, int eval_point, const std::vector<G1>& setup_g1) {
    MultiProof result;
    result.evaluations.reserve(polys.size());
    
    for (const auto& poly : polys) {
        Fr eval = evaluate_polynomial(poly, eval_point);
        result.evaluations.push_back(eval);
    }
    
    Fr gamma;
    gamma.setStr("12345"); //fixed challenge for testing, in practice Fiat-Shamir is used 
    Fr eval_point_fr;
    eval_point_fr.setStr(std::to_string(eval_point));
    //create combined polynomial: Σ(γ^i * p_i(x))
    Polynomial combined_poly;
    Fr gamma_power;
    gamma_power = Fr(1);
    
    for (size_t i = 0; i < polys.size(); ++i) {
        if (i == 0) {
            combined_poly = polys[i];
        } else {
            size_t max_size = std::max(combined_poly.size(), polys[i].size());
            if (combined_poly.size() < max_size) {
                combined_poly.resize(max_size);
            }
            
            for (size_t j = 0; j < polys[i].size(); ++j) {
                Fr term;
                Fr::mul(term, polys[i][j], gamma_power);
                Fr::add(combined_poly[j], combined_poly[j], term);
            }
        }
        Fr::mul(gamma_power, gamma_power, gamma);
    }
    
    Fr combined_eval;
    combined_eval.setStr("0");
    gamma_power.setStr("1");
    
    for (size_t i = 0; i < result.evaluations.size(); ++i) {
        Fr term;
        Fr::mul(term, result.evaluations[i], gamma_power);
        Fr::add(combined_eval, combined_eval, term);
        Fr::mul(gamma_power, gamma_power, gamma);
    }
    
    Polynomial divisor(2);
    divisor[0].setStr(std::to_string(-eval_point));
    divisor[1].setStr("1");
    
    Polynomial numerator = combined_poly;
    if (!numerator.empty()) {
        Fr::sub(numerator[0], numerator[0], combined_eval);
    }
    
    Polynomial quotient = poly_div(numerator, divisor);
    
    result.proof = commit(quotient, setup_g1);
    
    return result;
}

bool verify_multi_proof(const std::vector<G1>& commitments, const MultiProof& proof, int eval_point, const G2& s_g2) {
    if (commitments.size() != proof.evaluations.size()) {
        return false;
    }
    
    Fr gamma;
    gamma.setStr("12345"); 
    
    // Compute combined commitment: Σ(γ^i * C_i)
    G1 combined_commitment;
    combined_commitment.clear();
    Fr gamma_power;
    gamma_power.setStr("1");
    
    for (size_t i = 0; i < commitments.size(); ++i) {
        G1 term;
        G1::mul(term, commitments[i], gamma_power);
        G1::add(combined_commitment, combined_commitment, term);
        Fr::mul(gamma_power, gamma_power, gamma);
    }
    
    // Compute combined evaluation: Σ(γ^i * y_i)
    Fr combined_eval;
    combined_eval.setStr("0");
    gamma_power.setStr("1");
    
    for (size_t i = 0; i < proof.evaluations.size(); ++i) {
        Fr term;
        Fr::mul(term, proof.evaluations[i], gamma_power);
        Fr::add(combined_eval, combined_eval, term);
        Fr::mul(gamma_power, gamma_power, gamma);
    }
    
    // Verification: e(C - y*G1, G2) = e(π, (s - x)*G2)
    G2 g2;
    mapToG2(g2, 1);
    
    // Left side: combined_commitment - combined_eval * G1
    G1 lhs;
    G1 g1_gen;
    mapToG1(g1_gen, 1);
    G1 eval_g1;
    G1::mul(eval_g1, g1_gen, combined_eval);
    G1::sub(lhs, combined_commitment, eval_g1);
    
    // Right side: (s - eval_point) * G2
    Fr eval_point_fr;
    eval_point_fr.setStr(std::to_string(eval_point));
    G2 eval_point_g2;
    G2::mul(eval_point_g2, g2, eval_point_fr);
    G2 rhs_g2;
    G2::sub(rhs_g2, s_g2, eval_point_g2);
    
    //verify pairing equation
    mcl::Fp12 e1, e2;
    pairing(e1, lhs, g2);
    pairing(e2, proof.proof, rhs_g2);
    return e1 == e2;
}