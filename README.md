# RSA-homework

## A fast implementation of RSA.

1.It use Newton Method to optimize modular and division (during modular exponentiation operation).

2.Use CRT & Newton Method to accelerate the encoding and decoding process.

3.One unsigned int record 0 to 65536 (16 bits data).

4.Use 8 threads to generate primes.

*5.Use NTT to accelerate multiple but slower than simple method (It may be better to limit unsigned int in 32768).


## Points to optimize which not done:

1.Montgomery algorithm.

2.Variable length array for big integer.
