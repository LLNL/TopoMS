/*
 * Copyright (c) 2017 University of Utah 
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#define INDEX_TYPE long long	// type to index elements of a mesh - need to support big meshes! - keep this signed to make arithmetic consistent
//#define FLOATTYPE float				// type used for computing numerical integration
#define FLOATTYPE double				// type used for computing numerical integration
#define INT_TYPE int					// regular old ints
#define DIM_TYPE unsigned char			// used to query the dimension of a cell - we usually have values between 0-3
#define BYTE_TYPE unsigned char			// just a regular old byte
#define BOUNDARY_TYPE unsigned char		// used to query the boundary classification of a cell - usually a small integer
#define ASSIGNED_TYPE unsigned char		// used to query if a cell has been assigned - usually 0, 1 values

enum DestType {
	BACKGROUND,
	UNASSIGNED,
	ASSIGNED,
	CERTAIN_TERMINAL,
	CERTAIN_NONTERMINAL
};
#endif
