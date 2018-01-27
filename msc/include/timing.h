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

#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>
#include <vector>
#include <chrono>
#include "basic_types.h"


namespace MSC {

    typedef std::chrono::time_point < std::chrono::high_resolution_clock> TIMEP;
    class ThreadedTimer {
    protected:
        struct Timing {
            int activity;
            TIMEP start_time;
        };
        struct ThreadTimings {
            int thread_id;
            std::vector < Timing > timings;
        };

        ThreadTimings* m_thread_timings;
        int m_num_threads;
        TIMEP m_global_start;
        TIMEP m_global_end;
    public:



        ThreadedTimer(int num_threads) : m_num_threads(num_threads) {
            m_thread_timings = new ThreadTimings[m_num_threads];
        }

        void StartGlobal() {
            m_global_start = std::chrono::high_resolution_clock::now();
        }

        void EndGlobal() {
            m_global_end = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < m_num_threads; i++) this->RecordThreadActivity(i, 0);
        }

        void PrintAll() {
            printf(" took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(m_global_end - m_global_start).count());
        }

        void PrintAllToFile(char* fname, const char** strings) {

            FILE* fout = fopen(fname, "w");
            for (int i = 0; i < m_num_threads; i++) {
                ThreadTimings& tt = m_thread_timings[i];
                for (int a = 0; a < tt.timings.size() - 1; a++) {
                    fprintf(fout, "%d %d %d %s\n", i,
                        std::chrono::duration_cast<std::chrono::milliseconds>(tt.timings[a].start_time - m_global_start),
                        std::chrono::duration_cast<std::chrono::milliseconds>(tt.timings[a + 1].start_time - m_global_start), strings[tt.timings[a].activity]);


                }
            }

            fclose(fout);
        }
        void RecordThreadActivity(int thread_num, int activity) {
            Timing t = { activity, std::chrono::high_resolution_clock::now() };
            m_thread_timings[thread_num].timings.push_back(t);
        }

    };

}

#endif
