#include "kzg.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <cassert>

//change paths in makefile if needed 
//run with: make && ./kzg_test
int main() {
    std::cout << "Run KZG proof test:\n";

    auto [setup_g1, s_g2] = trusted_setup(); //run trusted setup to get public parameters

    std::vector<uint8_t> data = {1, 2, 3, 4, 5}; //input data to encode
    auto [points, poly] = encode_as_polynomial(data); //encode data as polynomial

    auto C = commit(poly, setup_g1); //commit to the polynomial

    int eval_x = 2; //choose evaluation point
    Fr eval_y = evaluate_polynomial(poly, eval_x); //evaluate polynomial at eval_x
    auto eval_point = std::make_pair(eval_x, eval_y); //pair (x, y) for proof

    auto pi = create_proof(poly, eval_point, setup_g1); //create KZG proof for evaluation

    bool verification_result = verify(C, pi, eval_point, s_g2); //verify proof
    assert(verification_result); //should pass for correct proof
    std::cout << "Correct proof verified.\n";

    //test with wrong evaluation
    Fr wrong_y = eval_y; //copy correct y
    Fr one = Fr(1); //Fr element 1
    Fr::add(wrong_y, wrong_y, one); //add 1 to y to make it wrong
    auto wrong_point = std::make_pair(eval_x, wrong_y); //pair with wrong y
    
    assert(!verify(C, pi, wrong_point, s_g2)); //should fail for wrong proof
    std::cout << "Incorrect proof rejected.\n";

    std::cout << "SUCCESS!\n";
    return 0;
}