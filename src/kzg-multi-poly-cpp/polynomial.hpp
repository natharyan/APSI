#pragma once
#include <vector>
#include <cstdint>
#include <utility>
#include <mcl/bn.hpp>

using namespace mcl::bn;

typedef std::vector<Fr> Polynomial;
typedef std::pair<int, Fr> Point;

std::pair<std::vector<Point>, Polynomial> encode_as_polynomial(const std::vector<uint8_t>& data);
Fr evaluate_polynomial(const Polynomial& poly, int x);
Polynomial poly_div(const Polynomial& a, const Polynomial& b);