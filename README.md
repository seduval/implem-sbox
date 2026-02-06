# implem-sbox
A tool for finding good side-channel-protected circuits implementing small functions (cryptographic S-box).

This tool and its purpose are described in more detailed in the TCHES article (link will be provided soon).
A file containing the look-up-tables of many interesting S-boxes is given in `list_of_test_luts.txt`.

## Install instructions

- clone the repository
- build the program:
`cd Code`
`make`

- download the precomputation files for all S-box sizes from the archive (takes about 1h):
`wget -r -np -nH --cut-dirs=2 -A txt  https://caramba.loria.fr/sbox/precomputation_files/`
- or download the files individually for the S-box sizes you intend to study

## Program manual

 - The option '--lut' or '-l' requires the values of the LUT separated by ','
 - For example: '-l 0,1,5,4,8,3,7,2'
 - The number of values must be a power of 2 beetween 16 and 64
 - The option ' --andmax' or '-a' allows you to choose a maximum for the number of AND gates wanted in the implementation
 - The option '--solmax' or '-s' (optionnal) allows you to choose the number of different solutions wanted, by default the programm will try to find all of them

 Note: We have identified a small bug which only occurs with very simple Boolean functions. It is very rare on S-boxes on more than 5 bits. We are working on a fix.

## Example execution

This execution gives an implementation of the S-box corresponding to the inverse of the Q2256 representative (6 bits, degree 3), with at most 16 AND gates. This should run within a few seconds on a laptop.

```
./implem --lut 0,1,2,3,4,7,5,6,8,16,32,58,9,19,61,37,10,17,50,43,11,18,45,54,41,52,63,34,47,48,56,39,12,29,62,49,53,57,31,13,20,27,55,36,51,33,23,25,59,42,24,21,26,22,60,44,30,14,40,38,35,46,15,28 --andmax 16 -s 1
ANF is given by :
y0 = x0 ^ x1x2 ^ x0x3 ^ x2x3 ^ x1x2x3 ^ x2x4 ^ x1x2x4 ^ x3x4 ^ x0x3x4 ^ x2x5 ^ x0x2x5 ^ x1x2x5 ^ x0x3x5 ^ x1x3x5 ^ x2x3x5 ^ x4x5 ^ x1x4x5 ^ x2x4x5
y1 = x1 ^ x0x2 ^ x1x2 ^ x1x3 ^ x0x1x3 ^ x1x2x3 ^ x4 ^ x0x4 ^ x1x4 ^ x0x1x4 ^ x3x4 ^ x0x3x4 ^ x2x3x4 ^ x0x1x5 ^ x0x2x5 ^ x1x2x5 ^ x0x3x5 ^ x1x3x5 ^ x2x3x5 ^ x0x4x5 ^ x1x4x5 ^ x3x4x5
y2 = x2 ^ x2x3 ^ x1x2x3 ^ x2x4 ^ x1x2x4 ^ x0x3x4 ^ x1x3x4 ^ x5 ^ x0x1x5 ^ x2x5 ^ x0x2x5 ^ x0x3x5 ^ x4x5 ^ x2x4x5 ^ x3x4x5
y3 = x3 ^ x0x3 ^ x1x3 ^ x1x2x3 ^ x4 ^ x0x4 ^ x1x4 ^ x1x2x4 ^ x3x4 ^ x0x3x4 ^ x5 ^ x0x1x5 ^ x2x5 ^ x0x2x5 ^ x1x2x5 ^ x1x3x5 ^ x2x3x5 ^ x4x5 ^ x0x4x5 ^ x1x4x5 ^ x2x4x5
y4 = x0x3 ^ x1x2x3 ^ x0x4 ^ x1x4 ^ x1x2x4 ^ x0x3x4 ^ x0x5 ^ x1x5 ^ x0x1x5 ^ x2x5 ^ x0x2x5 ^ x1x2x5 ^ x3x5 ^ x1x3x5 ^ x2x3x5 ^ x4x5 ^ x0x4x5 ^ x2x4x5 ^ x3x4x5
y5 = x1x3 ^ x1x4 ^ x3x4 ^ x1x5 ^ x2x5 ^ x1x3x5 ^ x4x5 ^ x1x4x5
The given s-box can be implemented by :

uint32_t Sbox(uint32_t X, uint8_t size) {
        uint32_t x[size];
        for(uint8_t b=0; b<size; b++)  {
                x[b] = X&1ul;
                X >>= 1;
        }
        uint32_t y[size];
        uint32_t l1 = x[3] ^ x[5];
        uint32_t l0 = x[4] ^ l1;
        uint32_t l2 = x[5];
        uint32_t l3 = x[1] ^ x[3];
        uint32_t l4 = x[4] ^ x[5];
        uint32_t l5 = x[5] ^ l3;
        uint32_t l6 = x[3];
        uint32_t l7 = l1 ^ l3;
        uint32_t l8 = l4 ^ l7;
        uint32_t l9 = x[4];
        uint32_t l10 = x[0] ^ x[1];
        uint32_t l11 = x[5] ^ l10;
        uint32_t l13 = x[2] ^ l7;
        uint32_t l14 = l4 ^ l13;
        uint32_t l16 = l7 ^ l10;
        uint32_t l17 = l8 ^ l11;
        uint32_t l18 = x[2] ^ x[5];
        uint32_t l19 = l4 ^ l18;
        uint32_t l20 = l3 ^ l17;
        uint32_t l21 = l14 ^ l20;
        uint32_t l12 = x[5] ^ l21;
        uint32_t l15 = l7 ^ l21;
        uint32_t q0 = l2 & l3;
        uint32_t q1 = l4 & l5 ^ l2;
        uint32_t q2 = l6 & l7;
        uint32_t q3 = l6 & l8;
        uint32_t q4 = x[2] & l7;
        uint32_t q5 = x[0] & x[5];
        uint32_t q6 = l9 & l10;
        uint32_t q7 = l6 & l11;
        uint32_t q8 = l12 & l13 ^ x[2];
        uint32_t q9 = l10 & l14 ^ x[1];
        uint32_t q11 = q1 ^ q3;
        uint32_t q10 = q0 ^ q11;
        uint32_t q13 = q0 ^ q5;
        uint32_t q14 = q2 ^ q5;
        uint32_t q15 = q2 ^ q13;
        uint32_t q17 = q4 ^ q7;
        uint32_t tq0 = q2 ^ q4;
        uint32_t q12 = q11 ^ tq0;
        uint32_t q16 = q6 ^ q12;
        uint32_t tq1 = q7 ^ q8;
        uint32_t q18 = q13 ^ tq1;
        uint32_t ty4 = ((q9) ^ (l15)) & ((q16) ^ (l16)) ^ (q12);
        uint32_t ty2 = ((q8) ^ (l17)) & ((q17) ^ (l18)) ^ (q7);
        uint32_t ty1 = ty2 ^ ((q10) ^ (l19)) & ((q1) ^ (l20)) ^ (q3);
        uint32_t ty3 = ty4 ^ ((q0) ^ (l6)) & ((l4)) ^ (q14);
        uint32_t ty0 = ty4 ^ ((q13) ^ (l15)) & ((l17)) ^ (q18);
        uint32_t ty5 = ty4 ^ ty2 ^ ((q0) ^ (l0)) & ((q15) ^ (l21)) ^ (q13);
        y[4] = ty4 ^ x[0];
        y[2] = ty2 ^ x[2];
        y[1] = ty1 ^ x[1] ^ x[2];
        y[3] = ty3 ^ l0 ^ x[0];
        y[0] = ty0;
        y[5] = ty5 ^ l1 ^ x[0] ^ x[2];
        uint32_t Y = 0;
        for(uint8_t b=0; b<size; b++) {
                Y ^= ((Y>>b) ^ y[b]) << b;
        }
        return Y;
}
Nb_and = 16
Nb_xor = 51
```

## Expected timings

Degree 4 on 6 bits takes half an hour, degree 5 on 6 bits takes around a day, degree 3 on 7 bits takes few hours, all other results are obtained within seconds on a laptop.
