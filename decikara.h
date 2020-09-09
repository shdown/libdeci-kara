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

#pragma once

#include <stddef.h>
#include "deci.h"

// The 'cutoff' argument is the threshold below which the multiplication should be done using the
// "basecase"/"dumb"/quadratic algorithm.

// Returns '(size_t) -1' on overflow.
size_t decikara_nscratch(
    size_t nwa,
    size_t nwb,
    size_t cutoff);

// Does not modify either (wa ... wa+nwa) or (wb... wb+nwb), writes result into
// (out ... out+nwa+nwb).
//
// 'scratch' must have capacity of 'decikara_nscratch(nwa, nwb, cutoff)'.
void decikara_mul(
    deci_UWORD *wa, size_t nwa,
    deci_UWORD *wb, size_t nwb,
    deci_UWORD *scratch,
    deci_UWORD *out,
    size_t cutoff);
