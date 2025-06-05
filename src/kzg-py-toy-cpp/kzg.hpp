#pragma once
#include <vector>
#include <tuple>
#include <utility>
#include <mcl/bn.hpp>

using namespace mcl::bn;

typedef std::vector<Fr> Polynomial;
typedef std::pair<int, Fr> Point;

std::pair<std::vector<G1>, G2> trusted_setup();
G1 commit(const Polynomial& poly, const std::vector<G1>& setup_g1);
G1 create_proof(const Polynomial& poly, const Point& point, const std::vector<G1>& setup_g1);
bool verify(const G1& commitment, const G1& proof, const Point& point, const G2& s_g2);