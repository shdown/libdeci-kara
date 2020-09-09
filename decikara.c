/*
 * Copyright (C) 2020  libdeci-kara developers
 *
 * This file is part of libdeci-kara.
 *
 * libdeci-kara is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libdeci-kara is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libdeci-kara.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "decikara.h"

#include <stdbool.h>

// See https://gmplib.org/manual/Karatsuba-Multiplication -- it talks about splitting at a power of
// two, but the formulas hold for arbitrary base as well.
//
// >  high              low
// > +----------+----------+
// > |    x1    |    x0    |
// > +----------+----------+
// >
// > +----------+----------+
// > |    y1    |    y0    |
// > +----------+----------+
// >
// > Let b be the power of 2 where the split occurs, i.e. if x0 is k limbs (y0 the same) then
// > b=2^(k*mp_bits_per_limb). With that x=x1*b+x0 and y=y1*b+y0, and the following holds,
// >
// >     x*y = (b^2+b)*x1*y1 - b*(x1-x0)*(y1-y0) + (b+1)*x0*y0
// >
// > This formula means doing only three multiplies of (N/2)x(N/2) limbs, whereas a basecase
// > multiply of NxN limbs is equivalent to four multiplies of (N/2)x(N/2). The factors (b^2+b) etc
// > represent the positions where the three products must be added.
// >
// >  high                              low
// > +--------+--------+ +--------+--------+
// > |      x1*y1      | |      x0*y0      |
// > +--------+--------+ +--------+--------+
// >           +--------+--------+
// >       add |      x1*y1      |
// >           +--------+--------+
// >           +--------+--------+
// >       add |      x0*y0      |
// >           +--------+--------+
// >           +--------+--------+
// >       sub | (x1-x0)*(y1-y0) |
// >           +--------+--------+
// >
// > The term (x1-x0)*(y1-y0) is best calculated as an absolute value, and the sign used to choose
// > to add or subtract. Notice the sum high(x0*y0)+low(x1*y1) occurs twice, so it’s possible to do
// > 5*k limb additions, rather than 6*k, but in GMP extra function call overheads outweigh the
// > saving.
// [...]
// > A similar formula for both multiplying and squaring can be constructed with a middle term
// > (x1+x0)*(y1+y0). But those sums can exceed k limbs, leading to more carry handling and
// > additions than the form above.

// So we use the formula above; and also introduce use the following notation:
//
//  * n_1 and n_2 are the sizes of the multiplicands in words;
//
//  * m = floor(min(n_1, n_2) / 2);
//
//  * the split of the multiplicands, x->(x1,x0) and y->(y1,y0), is done at b=DECI_BASE^m;
//
//  * K(n_1, n_2) is the extra space in words needed to multiply n_1 words by n_2 words with our
//    algorithm.
//
// So let's calculate the amout of extra space we need. We are going to perform the following
// low-level steps:
//
// 1. Multiply x0*y0 -> low part of output.
//
//    Extra space needed: K(m, m).
//
// 2. Multiply x1*y1 -> high part of output.
//
//    Extra space needed: K(n_1 - m, n_2 - m).
//
// 3. Copying some stuff around and perform the two additions to the middle (x1*y1 and x0*y0).
//
//    Extra space needed: m.
//
// 4. Multiply (x1-x0)*(y1-y0) -> scratch.
//
//    Extra space needed: 2*(n_1 + n_2 - 2*m) + K(n_1 - m, n_2 - m).
//
// 5. Subtract the result of the multiplication above from the middle.
//
//    Extra space needed: no.
//
// So, if both n_1 and n_2 are not less than the cutoff, then
//   K(n_1, n_2) = 2*(h_1 + h_2) + K(h_1, h_2),
// where:
//   * h_1 = n_1 - m;
//   * h_2 = n_2 - m;
//   * m = floor(min(n_1, n_2) / 2).

// Saturating addition for 'size_t'.
static inline size_t add_zu(size_t a, size_t b)
{
    size_t r = a + b;
    return r < a ? ((size_t) -1) : r;
}

size_t decikara_nscratch(
    size_t nwa,
    size_t nwb,
    size_t cutoff)
{
    size_t min_size = nwa < nwb ? nwa : nwb;

    if (min_size < cutoff)
        return 0;

    size_t m = min_size / 2;
    size_t nhi_a = nwa - m;
    size_t nhi_b = nwb - m;

    size_t extra = decikara_nscratch(nhi_a, nhi_b, cutoff);

    size_t nhi = add_zu(nhi_a, nhi_b);
    size_t twice_nhi = add_zu(nhi, nhi);
    return add_zu(extra, twice_nhi);
}

static bool copy_hi_sub_lo(deci_UWORD *src, size_t nlo, size_t nhi, deci_UWORD *out)
{
    deci_memcpy(out, src + nlo, nhi);
    return deci_sub(out, out + nhi, src, src + nlo);
}

static void add_or_sub_shifted(
        deci_UWORD *wa, size_t nwa,
        deci_UWORD *wb, size_t nwb,
        size_t shift, bool add)
{
    if (add) {
        (void) deci_add(
            wa + shift,     wa + nwa,
            wb,             wb + nwb);
    } else {
        bool borrow = deci_sub_raw(
            wa + shift,     wa + nwa,
            wb,             wb + nwb);
        if (borrow)
            deci_uncomplement(wa, wa + nwa);
    }
}

void decikara_mul(
    deci_UWORD *wa, size_t nwa,
    deci_UWORD *wb, size_t nwb,
    deci_UWORD *scratch,
    deci_UWORD *out,
    size_t cutoff)
{
    size_t min_size = nwa < nwb ? nwa : nwb;

    if (min_size < cutoff) {
        deci_zero_out_n(out, nwa + nwb);
        deci_mul(wa, wa + nwa, wb, wb + nwb, out);
        return;
    }

    size_t m = min_size / 2;
    size_t nhi_a = nwa - m;
    size_t nhi_b = nwb - m;
    size_t nhi = nhi_a + nhi_b;

    // step 1

    decikara_mul(
        /*wa=*/wa, /*nwa=*/m,
        /*wb=*/wb, /*nwb=*/m,
        /*scratch=*/scratch,
        /*out=*/out,
        /*cutoff=*/cutoff);

    // step 2

    decikara_mul(
        /*wa=*/wa + m, /*nwa=*/nhi_a,
        /*wb=*/wb + m, /*nwb=*/nhi_b,
        /*scratch=*/scratch,
        /*out=*/out + m * 2,
        /*cutoff=*/cutoff);

    // step 3

    deci_memcpy(scratch, out, m * 2);

    add_or_sub_shifted(
        /*wa=*/out,         /*nwa=*/nwa + nwb,
        /*wb=*/out + m * 2, /*nwb=*/nhi,
        /*shift=*/m,
        /*add=*/true);

    add_or_sub_shifted(
        /*wa=*/out,     /*nwa=*/nwa + nwb,
        /*wb=*/scratch, /*nwb=*/m * 2,
        /*shift=*/m,
        /*add=*/true);

    // step 4

    bool sign = copy_hi_sub_lo(
        /*src=*/wa, /*nlo=*/m, /*nhi=*/nhi_a, /*out=*/scratch);
    sign ^= copy_hi_sub_lo(
        /*src=*/wb, /*nlo=*/m, /*nhi=*/nhi_b, /*out=*/scratch + nhi_a);

    decikara_mul(
        /*wa=*/scratch,         /*nwa=*/nhi_a,
        /*wb=*/scratch + nhi_a, /*nwb=*/nhi_b,
        /*scratch=*/scratch + 2 * nhi,
        /*out=*/scratch + nhi,
        /*cutoff=*/cutoff);

    // step 5

    add_or_sub_shifted(
        /*wa=*/out,           /*nwa=*/nwa + nwb,
        /*wb=*/scratch + nhi, /*nwb=*/nhi,
        /*shift=*/m,
        /*add=*/sign);
}
