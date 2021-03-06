/* This file is part of nfftls.
 *
 * Nfftls is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Nfftls is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with nfftls.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) 2012, 2016 by B. Leroy
 */
#ifndef NFFT_H
#define NFFT_H

#include <complex.h>
#include <stdbool.h>

void nfft(const double* t, const double* y, int* n, int* m, complex double* d);

#endif /* !NFFT_H */
