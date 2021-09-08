// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include <chrono>
#include "examples.h"
#include <seal/util/polyarithsmallmod.h>

using namespace std;
using namespace seal;

void example_ntt()
{
    print_example_banner("Example of NTT core");

    uint64_t q = 1152921504598720513;
    Modulus modulus(q);

    size_t poly_modulus_degree = 65536;

    seal::util::NTTTables ntt_table(16, modulus);

    cout << "This is NTT example." << endl;

    vector<uint64_t> test_vector(poly_modulus_degree);
    vector<uint64_t> test_vector1(poly_modulus_degree);
    vector<uint64_t> test_vector2(poly_modulus_degree);
    vector<uint64_t> test_vector3(poly_modulus_degree);

    uint32_t count_psi = 0;
    ifstream idt, idtt;
    idt.open("D:/SEAL/tv.csv", ios::in | ios::binary);
    idtt.open("D:/SEAL/inttt.csv", ios::in | ios::binary);
    uint64_t cur;
    while (idt >> cur)
    {
        test_vector1[count_psi] = cur;
        test_vector2[count_psi] = cur;
        count_psi++;
    }
    idt.close();

    //count_psi = 0;
    //uint64_t count = 0;
    //uint64_t curt;
    //while (idtt >> curt)
    //{
    //    test_vector3[count_psi] = curt;
    //    if (test_vector1[count_psi] != test_vector3[count_psi])
    //    {
    //        count++;
    //        cout << "test_vector1[" << count_psi << "] = " << test_vector1[count_psi] << endl;

    //    }
    //    count_psi++;
    //}
    //idtt.close();

    //cout << "Count = " << count << endl;

    for (uint64_t i = 0; i < poly_modulus_degree; i++)
    {
        test_vector[i] = i;
    }

    //fstream tv;
    //tv.open("D:/SEAL/op1.csv", ios::out);

    //for (int i = 0; i < poly_modulus_degree; i++)
    //{
    //    tv << test_vector1[i] << endl;
    //}
    //tv.close();

    fstream psi;
    psi.open("D:/SEAL/op2.csv", ios::out);

    for (int j = 0; j < poly_modulus_degree; j++)
    {
        psi << ntt_table.get_from_root_powers(j).operand << endl;
    }
    psi.close();

    fstream psiinv;
    psiinv.open("D:/SEAL/op3.csv", ios::out);

    for (int j = 0; j < poly_modulus_degree; j++)
    {
        psiinv << ntt_table.get_from_inv_root_powers(j).operand << endl;
    }
    psiinv.close();

    //for (int i = 0; i < poly_modulus_degree; i++)
    //{
    //    if (test_vector1[i] >= 1152921504598720513)
    //    {
    //        test_vector1[i] = test_vector1[i] - 1152921504598720513;
    //    }
    //     else
    //    {
    //         test_vector1[i] = test_vector1[i];
    //    }
    //}

    util::ntt_negacyclic_harvey(test_vector1, ntt_table);

    fstream ntt;
    ntt.open("D:/SEAL/ntt.csv", ios::out);

    for (int i = 0; i < poly_modulus_degree; i++)
    {
        ntt << test_vector1[i] << endl;
    }
    ntt.close();

    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$:" << endl;
    auto start = std::chrono::steady_clock::now();
    util::inverse_ntt_negacyclic_harvey(test_vector, ntt_table);
    auto diff = std::chrono::steady_clock::now() - start;
    std::cout << "NTT time = " << chrono::duration<double, milli>(diff).count() << " ms" << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$:" << endl;


    fstream intt;
    intt.open("D:/SEAL/intt.csv", ios::out);

    for (int i = 0; i < poly_modulus_degree; i++)
    {
        intt << test_vector[i] << endl;
    }
    intt.close();

    uint64_t root = ntt_table.get_root();
    uint64_t inverse_root;
    util::try_invert_uint_mod(root, modulus, inverse_root);
    cout << "modulus value: " << modulus.value() << endl;
    cout << "ratio 0: " << modulus.const_ratio()[0] << endl;
    cout << "ratio 1: " << modulus.const_ratio()[1] << endl;
    cout << "psi value: " << root << endl;
    cout << "psi inverse value: " << inverse_root << endl;
    cout << "inverse n value: " << ntt_table.inv_degree_modulo().operand << endl << endl;

    ///////////////////
    print_example_banner("Example of DYADIC core");

    vector<uint64_t> result(poly_modulus_degree);

    for (uint64_t i = 0; i < poly_modulus_degree; i++)
    {
        test_vector1[i] = i * i * i;
        test_vector2[i] = i * i + i + i;
    }

    util::dyadic_product_coeffmod(test_vector1, test_vector2, poly_modulus_degree, modulus, result);

    cout << "modulus value: " << modulus.value() << endl;
    cout << "ratio 0: " << modulus.const_ratio()[0] << endl;
    cout << "ratio 1: " << modulus.const_ratio()[1] << endl;

    fstream tv1;
    tv1.open("op11.csv", ios::out);
    fstream tv2;
    tv2.open("op21.csv", ios::out);
    fstream res;
    res.open("result.csv", ios::out);

    for (int i = 0; i < poly_modulus_degree; i++)
    {
        tv1 << test_vector1[i] << endl;
        tv2 << test_vector2[i] << endl;
        res << result[i] << endl;
    }

    tv1.close();
    tv2.close();
    res.close();

}

