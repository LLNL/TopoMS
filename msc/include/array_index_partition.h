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

#ifndef ARRAY_INDEX_PARTITION_H
#define ARRAY_INDEX_PARTITION_H


#include <math.h>
#include <stdio.h>

namespace MSC {

	class ArrayIndexPartitioner {
	protected:

	public:

		// fills the partition array with numbers such that the 0th element is 0, and last element
		// is the number of elements total. A task with index i should do work between
		// the interval from index >= partition[i] to index < partition[i+1] 
		// equalizes the amount of work between the numthreads indices 
		static void EvenChunkSplit(INDEX_TYPE num_elements, int numThreads, std::vector<INDEX_TYPE>& partition)  {
			partition.clear();
			partition.reserve(numThreads + 1);
			INDEX_TYPE chunksize = num_elements / numThreads;
			INDEX_TYPE remainder = num_elements % numThreads;
			INDEX_TYPE start = 0;
			partition.push_back(start);
			for (int i = 0; i < numThreads; i++) {
				start += chunksize;
				if (i < remainder) start += 1;
				partition.push_back(start);
			}
		}
	};


}

#endif