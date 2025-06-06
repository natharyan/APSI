#pragma once
#include <vector>
#include <tuple>
#include <utility>
#include <mcl/bn.hpp>

using namespace mcl::bn;

typedef std::vector<Fr> Polynomial;
typedef std::pair<int, Fr> Point;
typedef std::vector<Polynomial> MultiPolynomial;
typedef std::vector<Fr> Evaluations;

struct MultiProof {
    G1 proof;
    std::vector<Fr> evaluations;
};

std::pair<std::vector<G1>, G2> trusted_setup();
G1 commit(const Polynomial& poly, const std::vector<G1>& setup_g1);
std::vector<G1> multi_commit(const MultiPolynomial& polys, const std::vector<G1>& setup_g1);
MultiProof create_multi_proof(const MultiPolynomial& polys, int eval_point, const std::vector<G1>& setup_g1);
bool verify_multi_proof(const std::vector<G1>& commitments, const MultiProof& proof, int eval_point, const G2& s_g2);