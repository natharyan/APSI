#include "kzg.hpp"
#include "polynomial.hpp"
#include <iostream>
#include <cassert>

//change paths in makefile if needed 
//run with: make && ./kzg_test
int main() {
    std::cout << "Run Multi-Polynomial KZG proof test:\n";

    auto [setup_g1, s_g2] = trusted_setup();

    //create multiple polynomials 
    MultiPolynomial polys;
    
    Polynomial poly1(3);
    poly1[0] = 1;
    poly1[1] = 2;
    poly1[2] = 3;
    polys.push_back(poly1);
    
    Polynomial poly2(2);
    poly2[0] = 4;
    poly2[1] = 5;
    polys.push_back(poly2);
    
    Polynomial poly3(4);
    poly3[0] = 6;
    poly3[1] = 7;
    poly3[2] = 8;
    poly3[3] = 9;
    polys.push_back(poly3);

    auto commitments = multi_commit(polys, setup_g1); //run trusted setup to get public parameters
    
    std::cout << "Created " << commitments.size() << " commitments\n";

    //creating a multi-proof at evaluation point x = 2
    int eval_point = 2;
    auto multi_proof = create_multi_proof(polys, eval_point, setup_g1);
    
    std::cout << "Created multi-proof with " << multi_proof.evaluations.size() << " evaluations\n";

    //verifying the multi-proof
    bool verification_result = verify_multi_proof(commitments, multi_proof, eval_point, s_g2);
    
    std::cout << "Evaluation point: " << eval_point << std::endl;
    std::cout << "Number of evaluations: " << multi_proof.evaluations.size() << std::endl;
    std::cout << "Verification result: " << (verification_result ? "PASS" : "FAIL") << std::endl;
    
    if (!verification_result) {
        std::cout << "Testing individual polynomial verifications:\n";
        for (size_t i = 0; i < polys.size(); ++i) {
            G1 single_commitment = commit(polys[i], setup_g1);
            Fr eval = evaluate_polynomial(polys[i], eval_point);
            
            Polynomial divisor(2);
            divisor[0] = Fr(-eval_point);
            divisor[1] = Fr(1);
            
            Polynomial numerator = polys[i];
            Fr::sub(numerator[0], numerator[0], eval);
            
            Polynomial quotient = poly_div(numerator, divisor);
            G1 single_proof = commit(quotient, setup_g1);
            
            G2 g2;
            mapToG2(g2, 1);
            
            G1 lhs;
            G1 g1_gen;
            mapToG1(g1_gen, 1);
            G1 eval_g1;
            G1::mul(eval_g1, g1_gen, eval);
            G1::sub(lhs, single_commitment, eval_g1);
            
            Fr eval_point_fr;
            eval_point_fr.setStr(std::to_string(eval_point));
            G2 eval_point_g2;
            G2::mul(eval_point_g2, g2, eval_point_fr);
            G2 rhs_g2;
            G2::sub(rhs_g2, s_g2, eval_point_g2);
            
            mcl::Fp12 e1, e2;
            pairing(e1, lhs, g2);
            pairing(e2, single_proof, rhs_g2);
            bool single_verify = (e1 == e2);
            std::cout << "Polynomial " << i << " individual verification: " << (single_verify ? "PASS" : "FAIL") << std::endl;
        }
    }
    
    assert(verification_result);
    std::cout << "Multi-proof verified successfully!\n";
    
    MultiProof tampered_proof = multi_proof;
    Fr one;
    one.setStr("1");
    Fr::add(tampered_proof.evaluations[0], tampered_proof.evaluations[0], one);
    
    assert(!verify_multi_proof(commitments, tampered_proof, eval_point, s_g2));
    std::cout << "Tampered proof correctly rejected!\n";

    std::cout << "SUCCESS!\n";
    return 0;
}