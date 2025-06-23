#include <iostream>
#include <sstream>
#include <vector>
#include <experimental/filesystem>
#include <algorithm>
#include <fstream>
#include <sys/stat.h>
#include <cmath>
#include <cstring>
#include <cassert>
#include <bitset>
#include <array>
#include <map>
#include <chrono>
#include <string.h>
#include <omp.h>
#include <bits/stdc++.h>
#include <cctype>

typedef uint32_t converted_row_t;

// Useful for computing frequency of occurrence
struct UInt64DefaultedToZero
{
    uint64_t val = 0;
};

#define BASE64_TO_INT(X) std::bitset<6>(static_cast<int>((X)) - 65)

#define UPPER_BOUND_CSV_LINE_LENGTH_EXCLUDING_TIMESTAMP 50

#if DEBUG
    #define NUM_WINDOWS_PER_THREAD 500000
    #define NUM_SAMPLES_PER_THREAD_FORMATTING 500000
    #define NUM_SAMPLES_PER_THREAD_TRIMMING 500000
    #define NUM_PARALLEL_THREADS 1
    #define DEBUG_PRINTF(f,...) printf(f,##__VA_ARGS__)
#else
    // For parallelization
    #define NUM_WINDOWS_PER_THREAD 100000
    #define NUM_SAMPLES_PER_THREAD_FORMATTING 100000
    #define NUM_SAMPLES_PER_THREAD_TRIMMING 100000
    #define NUM_PARALLEL_THREADS 32
    #define DEBUG_PRINTF(f,...)
#endif

template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>> &param)
{
    assert(param.size() > 0);

    std::vector<std::vector<T>> ret(param[0].size(), std::vector<T>());

    for (uint64_t i = 0; i < param.size(); i++)
    {
        for (uint64_t j = 0; j < param[i].size(); j++)
        {
            ret[j].push_back(param[i][j]);
        }
    }
    return ret;
}

inline bool is_row_full_of_x(const std::vector<size_t> &considered_col_ids, const std::vector<size_t> &samples_line, size_t x)
{
    for (size_t considered_col_id : considered_col_ids)
        if (samples_line[considered_col_id] != x)
            return false;
    return true;
}

void print_row(const std::vector<size_t> &samples_line)
{
    for (auto elem : samples_line)
    {
        printf("%ld ", elem);
    }
    printf("\n");
}

size_t get_col_id_by_label(const std::vector<std::string> req_labels, const std::vector<std::string> &labels)
{
    std::stringstream in;
    for (size_t label_id = 0; label_id < labels.size(); label_id++)
    {
        std::string label = labels[label_id];
        std::transform(label.begin(), label.end(),
                       label.begin(),
                       [](unsigned char c)
                       { return std::tolower(c); });

        in << labels[label_id] << ",";
        if (std::find(req_labels.begin(), req_labels.end(), label) != req_labels.end())
        {
            return label_id;
        }
    }

    std::stringstream out;
    for (size_t req_labels_id = 0; req_labels_id < req_labels.size(); req_labels_id++)
    {
        out << req_labels[req_labels_id] << ",";
    }
    std::cerr << "[-] Could not find required column(s) " << out.str() << " in the input file columns: " << in.str() << std::endl;
    exit(EXIT_FAILURE);
}

// Useful for computing frequency of occurrence
uint64_t row_to_int(const std::vector<size_t> &considered_col_ids, const std::vector<size_t> &samples_line)
{
    size_t ret = 0;
    for (size_t i = 0; i < considered_col_ids.size(); i++)
    {
        ret |= (samples_line[considered_col_ids[i]] << i);
    }
    return ret;
}

std::vector<size_t> get_considered_column_ids(const std::vector<std::string> &labels)
{
    std::vector<size_t> ret;

    // size_t ignored_col_id_0 = get_col_id_by_label("Time", labels);
    size_t ignored_col_id_0 = get_col_id_by_label({"ck0", "clk", "ck", "clk0"}, labels) - 1;
    size_t ignored_col_id_1 = get_col_id_by_label({"cs", "cs0", "cs0n", "cs0_n"}, labels) - 1;

    // labels.size()-1 because of the "Time" label
    for (size_t col_id = 0; col_id < labels.size() - 1; col_id++)
    {
        if (col_id != ignored_col_id_0 && col_id != ignored_col_id_1)
            ret.push_back(col_id);
    }
    return ret;
}

// @param lines are actually columns
void write_csv_fast(const std::string &filename_out, double hor_start, double hor_step_size,
                    std::vector<std::string> &labels, std::vector<std::vector<size_t>> &all_samples_per_line_transpose, bool trim)
{
    std::chrono::steady_clock::time_point chrono_begin;
    std::chrono::steady_clock::time_point chrono_end;

    chrono_begin = std::chrono::steady_clock::now();

    uint64_t num_samples = all_samples_per_line_transpose.size();
    printf("[+] Trimming output data. Found %lu samples.\n", num_samples);

    size_t clk_col_id = get_col_id_by_label({"ck0", "clk", "ck", "clk0"}, labels) - 1;
    size_t cs_col_id = get_col_id_by_label({"cs", "cs0", "cs0n", "cs0_n"}, labels) - 1;

    std::vector<size_t> considered_col_ids = get_considered_column_ids(labels);

#if DEBUG
    int num_threads_trimming = 1;
#else
    int num_threads_trimming = (num_samples - 1 + NUM_SAMPLES_PER_THREAD_TRIMMING - 1) / NUM_SAMPLES_PER_THREAD_TRIMMING;
#endif

    std::vector<std::vector<uint64_t>> rowstokeep_line_numbers(num_threads_trimming);         // One per thread
    std::vector<std::vector<uint64_t>> rowstokeep_rownumber2cyclecount(num_threads_trimming); // One per thread
    std::vector<uint64_t> cycle_cnt_per_thread(num_threads_trimming);                         // One per thread

    if (trim)
    {
#pragma omp parallel for num_threads(NUM_PARALLEL_THREADS)
        for (int thread_id = 0; thread_id < num_threads_trimming; ++thread_id)
        {
            int64_t sample_clk_rising_row = -1;
            int64_t sample_cs_changed_row = -1;

            cycle_cnt_per_thread[thread_id] = 0;

            std::vector<uint64_t> samples_same_01_cycle;

            uint64_t max_sample_id = (thread_id + 1) * NUM_SAMPLES_PER_THREAD_TRIMMING;
            if (max_sample_id > num_samples - 1)
                max_sample_id = num_samples - 1;

            // Go back to the left until we find the previous 1 to 0 transition or hit zero
            uint64_t sample_id = thread_id * NUM_SAMPLES_PER_THREAD_TRIMMING;
            bool should_start_with_zero = true;

            
            for (uint64_t offset_left = 0; offset_left < sample_id; offset_left++)
            {
                if (all_samples_per_line_transpose[sample_id - offset_left - 1][clk_col_id] == 1 && all_samples_per_line_transpose[sample_id - offset_left][clk_col_id] == 0)
                {
                    should_start_with_zero = false;
                    sample_id = sample_id - offset_left;
                    break;
                }
            }
            if (should_start_with_zero)
            sample_id = 0;
            
            for (; sample_id < max_sample_id; sample_id++)
            {
                DEBUG_PRINTF("ts=%.8e sample_id=%lu\n", hor_start + hor_step_size * sample_id, sample_id);
                // Check clock rising edge
                if (all_samples_per_line_transpose[sample_id][clk_col_id] == 0 && all_samples_per_line_transpose[sample_id + 1][clk_col_id] == 1)
                {
                    sample_clk_rising_row = sample_id + 1;
                }

                samples_same_01_cycle.push_back(sample_id);

                // Check changes in CS
                if (all_samples_per_line_transpose[sample_id][cs_col_id] != all_samples_per_line_transpose[sample_id + 1][cs_col_id]) {
                    // sample_cs_changed_row = sample_id + 1;
                    sample_cs_changed_row = sample_id;
                }
                
                DEBUG_PRINTF("sample_cs_changed_row=%ld sample_clk_rising_row=%ld\n", sample_cs_changed_row, sample_clk_rising_row);
                // Do processing on the falling clock edge
                if (all_samples_per_line_transpose[sample_id][clk_col_id] == 1 && all_samples_per_line_transpose[sample_id + 1][clk_col_id] == 0)
                {
                    DEBUG_PRINTF("IN all_samples_per_line_transpose[sample_id][clk_col_id] == 1 && all_samples_per_line_transpose[sample_id + 1][clk_col_id] == 0\n");
                    if (sample_cs_changed_row >= 0 && !is_row_full_of_x(considered_col_ids, all_samples_per_line_transpose[sample_cs_changed_row], 1))
                    {
                        rowstokeep_line_numbers[thread_id].push_back(sample_cs_changed_row);
                        rowstokeep_rownumber2cyclecount[thread_id].push_back(cycle_cnt_per_thread[thread_id]);
                        DEBUG_PRINTF("1| sample %lu keeping\n", sample_cs_changed_row);
                    }
                    else
                    {
                        DEBUG_PRINTF("ELSE\n");
                        std::map<uint64_t, UInt64DefaultedToZero> occurrences; // occurrences[row_to_int] = num_occurrences
                        std::map<converted_row_t, uint64_t> converted_row_to_row_id;
                        uint64_t total = 0; // num of non-full-one rows in this block
                        bool valid = false;
                        // Capture the occurrences
                        // We exclude full-ones rows from the block
                        for (auto ridx : samples_same_01_cycle)
                        {
                            total++;
                            converted_row_t converted_row = row_to_int(considered_col_ids, all_samples_per_line_transpose[ridx]);
                            if (!is_row_full_of_x(considered_col_ids, all_samples_per_line_transpose[ridx], 1))
                            { // FUTURE: to be even faster, we could compare the converted row against a very specific integer full of a given number of ones
                                occurrences[converted_row].val++;
                                DEBUG_PRINTF("converted_row_to_row_id[converted_row] = converted_row_to_row_id[%d] = %lu\n", converted_row, ridx);
                                converted_row_to_row_id[converted_row] = ridx;
                                valid = true;
                                DEBUG_PRINTF("valid = true\n");
                            }
                        }
                        // If there was at least one non-full-ones row in the block
                        uint64_t max_occurrence_id = 0;
                        uint64_t max_occurrence_amount = 0;
                        if (valid)
                        {
                            // Find the index of max occurrence
                            auto max_occurrence = std::max_element(occurrences.begin(), occurrences.end(),
                                                                   [](const auto &a, const auto &b)
                                                                   {
                                                                       return a.second.val < b.second.val;
                                                                   });
                            max_occurrence_id = max_occurrence->first;
                            max_occurrence_amount = max_occurrence->second.val;
                            DEBUG_PRINTF("max_occurrence_id=%lu max_occurrence_amount=%lu total=%lu\n", max_occurrence_id, max_occurrence_amount, total);
                        }
                        // The rule "0.4 * total" is problematic because if the value changed before the current 01 clock-cycle, then we will count
                        // fewer occurences than there actually are. So we are more conservative here and just require at least 2 occurences.
                        // In a future version, we should consider the number of occurences in the previous 01 clock-cycle if the signal stabilized before.
                        // if (valid && max_occurrence_amount >= 0.4 * total)
                        if (valid && max_occurrence_amount >= 2)
                        {
                            assert(!is_row_full_of_x(considered_col_ids, all_samples_per_line_transpose[converted_row_to_row_id[max_occurrence_id]], 1));
                            rowstokeep_line_numbers[thread_id].push_back(converted_row_to_row_id[max_occurrence_id]);
                            DEBUG_PRINTF("2| sample %lu keeping\n", converted_row_to_row_id[max_occurrence_id]);
                            rowstokeep_rownumber2cyclecount[thread_id].push_back(cycle_cnt_per_thread[thread_id]);
                        }
                        else if (sample_clk_rising_row >= 0 && !is_row_full_of_x(considered_col_ids, all_samples_per_line_transpose[sample_clk_rising_row], 1))
                        {
                            rowstokeep_line_numbers[thread_id].push_back(sample_clk_rising_row);
                            DEBUG_PRINTF("3| sample %lu keeping\n", sample_clk_rising_row);
                            rowstokeep_rownumber2cyclecount[thread_id].push_back(cycle_cnt_per_thread[thread_id]);
                        }
                    }
                    sample_cs_changed_row = -1;
                    sample_clk_rising_row = -1;
                    samples_same_01_cycle.clear();
                    cycle_cnt_per_thread[thread_id]++;
                }
            }
        }
    }
    else
    {
        int thread_id = 0;
        // Keep all rows
        for (size_t sample_id = 0; sample_id < num_samples; sample_id++)
        {

            if (!is_row_full_of_x(considered_col_ids, all_samples_per_line_transpose[sample_id], 0))
            {
                rowstokeep_line_numbers[thread_id].push_back(sample_id);
                rowstokeep_rownumber2cyclecount[thread_id].push_back(cycle_cnt_per_thread[thread_id]);
            }

            // Do processing on the falling clock edge
            if (all_samples_per_line_transpose[sample_id][clk_col_id] == 1 && all_samples_per_line_transpose[sample_id + 1][clk_col_id] == 0)
            {
                cycle_cnt_per_thread[thread_id]++;
            }
        }
    }

    // Fill the trimmed lines and add the cycle count
    std::vector<std::vector<size_t>> all_samples_per_line_transpose_keep;
    for (int thread_id = 0; thread_id < num_threads_trimming; ++thread_id)
    {
        // Compute the cycle count offset
        uint64_t cycle_cnt_offset = 0;
        for (int prev_thread_id = 0; prev_thread_id < thread_id; ++prev_thread_id)
        {
            cycle_cnt_offset += cycle_cnt_per_thread[prev_thread_id];
        }

        uint64_t id_in_rows_to_keep = 0;
        for (auto id_to_keep : rowstokeep_line_numbers[thread_id])
        {
            auto new_row = all_samples_per_line_transpose[id_to_keep];
            new_row.push_back(cycle_cnt_offset + rowstokeep_rownumber2cyclecount[thread_id][id_in_rows_to_keep]);
            // std::cout << "cycle_cnt_cumul " << cycle_cnt_offset + rowstokeep_rownumber2cyclecount[thread_id][id_in_rows_to_keep] << std::endl;

            // print_row(new_row);
            all_samples_per_line_transpose_keep.push_back(new_row);
            id_in_rows_to_keep++;
        }
    }

    // Merge rowstokeep_line_numbers since it is used for formatting
    std::vector<size_t> rowstokeep_line_numbers_flat;
    for (int thread_id = 0; thread_id < num_threads_trimming; ++thread_id)
    {
        for (auto id_to_keep : rowstokeep_line_numbers[thread_id])
        {
            rowstokeep_line_numbers_flat.push_back(id_to_keep);
        }
    }

    // Finally, re-transpose. That's not super efficient but it allows to re-use the legacy code below, and anyways at this stage the data is relatively small.
    std::vector<std::vector<size_t>> all_samples_per_line;
    if (all_samples_per_line_transpose_keep.size() > 0)
        all_samples_per_line = transpose(all_samples_per_line_transpose_keep);

    chrono_end = std::chrono::steady_clock::now();
    std::cout << "[+] Trimming done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

    FILE *f = fopen(filename_out.c_str(), "w");
    // Add the cycle_cnt label
    labels.push_back("cycle_cnt");
    // now print data as CSV to file
    // -- header

    bool first = true;
    for (const auto &l : labels)
    {
        if (!first)
            fprintf(f, ",");
        first = false;
        fprintf(f, "%s", l.c_str());
    }
    fprintf(f, "\n");
    // -- data

    uint64_t num_lines = labels.size();
    uint64_t num_samples_per_line;
    if (all_samples_per_line_transpose_keep.size() > 0)
        num_samples_per_line = all_samples_per_line[0].size();
    else
        num_samples_per_line = 0;

    printf("[+] Formatting the timestamps.\n");
    chrono_begin = std::chrono::steady_clock::now();
    std::vector<std::array<char, 18>> curTimes(num_samples_per_line);

    int num_threads_formatting = (num_samples_per_line + NUM_SAMPLES_PER_THREAD_FORMATTING - 1) / NUM_SAMPLES_PER_THREAD_FORMATTING;
#pragma omp parallel for
    for (int thread_id = 0; thread_id < num_threads_formatting; ++thread_id)
    {
        uint64_t max_sample_index = (thread_id + 1) * NUM_SAMPLES_PER_THREAD_FORMATTING;
        if (max_sample_index > num_samples_per_line)
            max_sample_index = num_samples_per_line;
        for (uint64_t sample_no = thread_id * NUM_SAMPLES_PER_THREAD_FORMATTING; sample_no < max_sample_index; ++sample_no)
        {
            // all_samples_per_line[num_lines-2][sample_no] is the cycle identifier
            sprintf(curTimes[sample_no].data(), "%.8e,", hor_start + hor_step_size * rowstokeep_line_numbers_flat[sample_no]);
        }
    }

    chrono_end = std::chrono::steady_clock::now();
    std::cout << "[+] Formatting the timestamps done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

    printf("[+] Generating the data lines.\n");
    chrono_begin = std::chrono::steady_clock::now();
    std::vector<std::array<char, 50>> datalines(num_samples_per_line);

#pragma omp parallel for
    for (int thread_id = 0; thread_id < num_threads_formatting; ++thread_id)
    {
        uint64_t max_row_index = (thread_id + 1) * NUM_SAMPLES_PER_THREAD_FORMATTING;
        if (max_row_index > num_samples_per_line)
            max_row_index = num_samples_per_line;
        for (uint64_t cur_row = thread_id * NUM_SAMPLES_PER_THREAD_FORMATTING; cur_row < max_row_index; ++cur_row)
        {
            for (uint64_t cur_index_in_row = 0; cur_index_in_row < num_lines - 2; ++cur_index_in_row)
            {
                datalines[cur_row][2 * cur_index_in_row] = '0' + all_samples_per_line.at(cur_index_in_row).at(cur_row);
                datalines[cur_row][2 * cur_index_in_row + 1] = ',';
            }
            // Print the cycle count
            size_t cur_index_in_row = num_lines - 2;
            std::string curr_index_in_row_str = std::to_string(all_samples_per_line.at(cur_index_in_row).at(cur_row)) + "\n\0";
            strncpy(datalines[cur_row].data() + 2 * cur_index_in_row, curr_index_in_row_str.c_str(), 50 - 2 * cur_index_in_row);
        }
    }

    chrono_end = std::chrono::steady_clock::now();
    std::cout << "[+] Generating the data lines done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

    printf("[+] Printing the data.\n");
    chrono_begin = std::chrono::steady_clock::now();

    for (size_t sample_no = 0; sample_no < num_samples_per_line; ++sample_no)
    {
        fputs(curTimes[sample_no].data(), f);
        fputs(datalines[sample_no].data(), f);
    }
    chrono_end = std::chrono::steady_clock::now();
    std::cout << "[+] Printing the data done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;
    fclose(f);
}
