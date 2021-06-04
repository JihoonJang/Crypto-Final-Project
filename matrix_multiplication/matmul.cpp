/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#include <helib/binaryArith.h>
#include <helib/helib.h>
#include <helib/intraSlot.h>

#include <ctime>
#include <iostream>
#include <random>

long bitSize = 2;
long outSize = 2 * bitSize;

class PlainMatrix {
 public:
  explicit PlainMatrix(int n) {
    mat_.assign(n, std::vector<long>());
    for (auto& vec : mat_) {
      vec.resize(n);
    }
    n_ = n;
  }

  void init_randomly() {
    for (int i = 0; i < n_; ++i) {
      for (int j = 0; j < n_; ++j) {
        mat_[i][j] = rand() % 2;
      }
    }
  }

  std::vector<long>& operator[](size_t i) { return mat_[i]; }

  void print() {
    std::cout << "Matrix : " << std::endl;
    for (size_t i = 0; i < mat_.size(); ++i) {
      for (size_t j = 0; j < mat_[i].size(); ++j) {
        std::cout << mat_[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

 private:
  std::vector<std::vector<long>> mat_;
  size_t n_;
};

template <typename t = std::vector<helib::Ctxt>>
class CipherMatrix {
 public:
  explicit CipherMatrix(int n) {
    mat_.resize(n);
    for (auto& vec : mat_) {
      vec.resize(n);
    }
    n_ = n;
  }

  void encrypt(PlainMatrix& other, const helib::PubKey& public_key,
               const helib::EncryptedArray& ea) {
    helib::Ctxt scratch(public_key);

    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        mat_[i][j].resize(bitSize, scratch);
        for (long k = 0; k < bitSize; ++k) {
          std::vector<long> a_vec(ea.size());
          // Extract the i'th bit of a,b,c.
          for (auto& slot : a_vec) slot = (other[i][j] >> k) & 1;
          ea.encrypt(mat_[i][j][k], public_key, a_vec);
        }
      }
    }
  }

  std::vector<t>& operator[](size_t i) { return mat_[i]; }

  CipherMatrix multiply(CipherMatrix& result,
                        CipherMatrix& other,  // Build the unpack slot encoding.
                        std::vector<helib::zzX>& unpackSlotEncoding) {
    CipherMatrix ret(n_);
    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        std::vector<std::vector<helib::Ctxt>> summands;
        for (size_t k = 0; k < n_; ++k) {
          std::vector<helib::Ctxt> encrypted_product;
          helib::CtPtrs_vectorCt product_wrapper(encrypted_product);
          helib::multTwoNumbers(
              product_wrapper, helib::CtPtrs_vectorCt((*this)[i][k]),
              helib::CtPtrs_vectorCt(other[k][j]),
              /*rhsTwosComplement=*/false,  // This means the rhs is unsigned
                                            // rather than 2's complement.
              outSize,  // Outsize is the limit on the number of bits in the
                        // output.
              &unpackSlotEncoding);  // Information needed for bootstrapping.
          summands.push_back(encrypted_product);
        }
        std::vector<helib::Ctxt> encrypted_result;
        helib::CtPtrs_vectorCt result_wrapper(encrypted_result);
        helib::CtPtrMat_vectorCt summands_wrapper(summands);
        helib::addManyNumbers(
            result_wrapper, summands_wrapper,
            0,  // sizeLimit=0 means use as many bits as needed.
            &unpackSlotEncoding);  // Information needed for bootstrapping.
        result[i][j] = std::move(encrypted_result);
      }
    }
    return ret;
  }

  PlainMatrix decrypt(helib::SecKey& secret_key,
                      const helib::EncryptedArray& ea) {
    std::vector<long> decrypted_result;

    PlainMatrix ret(n_);

    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        helib::decryptBinaryNums(decrypted_result,
                                 helib::CtPtrs_vectorCt(mat_[i][j]), secret_key,
                                 ea);
        ret[i][j] = decrypted_result.back();
      }
    }

    return ret;
  }

 private:
  std::vector<std::vector<t>> mat_;
  size_t n_;
};

int main(int argc, char* argv[]) {
  /*  Example of binary arithmetic using the BGV scheme  */

  // First set up parameters.
  //
  // NOTE: The parameters used in this example code are for demonstration only.
  // They were chosen to provide the best performance of execution while
  // providing the context to demonstrate how to use the "Binary Arithmetic
  // APIs". The parameters do not provide the security level that might be
  // required by real use/application scenarios.

  // Plaintext prime modulus.
  long p = 2;
  // Cyclotomic polynomial - defines phi(m).
  long m = 4095;
  // Hensel lifting (default = 1).
  long r = 1;
  // Number of bits of the modulus chain.
  long bits = 500;
  // Number of columns of Key-Switching matrix (typically 2 or 3).
  long c = 2;
  // Factorisation of m required for bootstrapping.
  std::vector<long> mvec = {7, 5, 9, 13};
  // Generating set of Zm* group.
  std::vector<long> gens = {2341, 3277, 911};
  // Orders of the previous generators.
  std::vector<long> ords = {6, 4, 6};

  std::cout << "Initialising context object..." << std::endl;
  // Initialize the context.
  // This object will hold information about the algebra created from the
  // previously set parameters.
  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .gens(gens)
                               .ords(ords)
                               .bits(bits)
                               .c(c)
                               .bootstrappable(true)
                               .mvec(mvec)
                               .build();

  // Print the context.
  context.printout();
  std::cout << std::endl;

  // Print the security level.
  std::cout << "Security: " << context.securityLevel() << std::endl;

  // Secret key management.
  std::cout << "Creating secret key..." << std::endl;
  // Create a secret key associated with the context.
  helib::SecKey secret_key(context);
  // Generate the secret key.
  secret_key.GenSecKey();

  // Generate bootstrapping data.
  secret_key.genRecryptData();

  // Public key management.
  // Set the secret key (upcast: SecKey is a subclass of PubKey).
  const helib::PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context.
  const helib::EncryptedArray& ea = context.getEA();

  // Build the unpack slot encoding.
  std::vector<helib::zzX> unpackSlotEncoding;
  buildUnpackSlotEncoding(unpackSlotEncoding, ea);

  // Get the number of slot (phi(m)).
  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;

  helib::Ctxt scratch(public_key);

  int matrix_size = 4;
  PlainMatrix pm1(matrix_size);

  pm1.init_randomly();
  pm1.print();

  PlainMatrix pm2(matrix_size);

  pm2.init_randomly();
  pm2.print();

  CipherMatrix cm1(matrix_size);

  cm1.encrypt(pm1, public_key, ea);

  CipherMatrix cm2(matrix_size);

  cm2.encrypt(pm2, public_key, ea);

  CipherMatrix cm3(matrix_size);

  clock_t start_time = clock();
  cm1.multiply(cm3, cm2, unpackSlotEncoding);
  clock_t end_time = clock();

  std::cout << "Matrix multiplication time : "
            << (double)(end_time - start_time) / CLOCKS_PER_SEC << " sec"
            << std::endl;

  start_time = clock();
  PlainMatrix pm3 = cm3.decrypt(secret_key, ea);
  end_time = clock();

  std::cout << "Matrix decryption time : "
            << (double)(end_time - start_time) / CLOCKS_PER_SEC << " sec"
            << std::endl;
  pm3.print();

  return 0;
}
