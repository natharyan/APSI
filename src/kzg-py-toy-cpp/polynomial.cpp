#include "polynomial.hpp"

//Evaluates a polynomial at a given integer x using Horner's method
Fr evaluate_polynomial(const Polynomial& poly, int x) {
    Fr result = 0;
    Fr power = 1;
    Fr x_fr = x;
    
    for (const auto& coeff : poly) {
        Fr term;
        Fr::mul(term, coeff, power);
        Fr::add(result, result, term);
        Fr::mul(power, power, x_fr);
    }
    return result;
}

//encodes a byte vector as a polynomial and returns evaluation points and the polynomial
std::pair<std::vector<Point>, Polynomial> encode_as_polynomial(const std::vector<uint8_t>& data) {
    Polynomial poly(data.size());
    std::vector<Point> points;
    
    for (size_t i = 0; i < data.size(); ++i) {
        poly[i] = data[i];
        Fr eval = evaluate_polynomial(poly, i);  
        points.emplace_back(i, eval); 
    }
    return {points, poly}; 
}

//divides polynomial a by polynomial b and returns the quotient polynomial
Polynomial poly_div(const Polynomial& a, const Polynomial& b) {
    if (b.empty() || (b.size() == 1 && b[0].isZero())) {
        throw std::invalid_argument("Div by zero polynomial"); 
    }
    
    Polynomial q;
    Polynomial r = a;
    
    Polynomial divisor = b;
    while (!divisor.empty() && divisor.back().isZero()) {
        divisor.pop_back(); 
    }
    
    if (divisor.empty()) {
        throw std::invalid_argument("Div by zero polynomial");
    }
    
    while (!r.empty() && r.size() >= divisor.size()) {
        while (!r.empty() && r.back().isZero()) {
            r.pop_back(); 
        }
        
        if (r.empty() || r.size() < divisor.size()) break;
        
        Fr coef;
        Fr::div(coef, r.back(), divisor.back()); 
        size_t deg = r.size() - divisor.size();
        
        if (q.size() <= deg) {
            q.resize(deg + 1); 
        }
        q[deg] = coef; 
        
        for (size_t i = 0; i < divisor.size(); ++i) {
            Fr temp;
            Fr::mul(temp, coef, divisor[i]);
            Fr::sub(r[deg + i], r[deg + i], temp); 
        }
    }
    
    while (!q.empty() && q.back().isZero()) {
        q.pop_back(); 
    }
    
    return q; 
}
