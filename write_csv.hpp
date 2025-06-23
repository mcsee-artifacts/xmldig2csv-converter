#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
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
#include <chrono>
#include <omp.h>

#define BASE64_TO_INT(X) std::bitset<6>(static_cast<int>((X))-65)

#define UPPER_BOUND_CSV_LINE_LENGTH_EXCLUDING_TIMESTAMP 50
// For parallelization
#define NUM_WINDOWS_PER_THREAD 100000
#define NUM_SAMPLES_PER_THREAD_FORMATTING 100000


// @param lines are actually columns
void write_csv_fast(const std::string &filename_out, long num_lines, long num_samples_per_line, double hor_start, double hor_step_size,
               std::vector<std::string> &labels, std::vector<std::vector<size_t>> &all_samples_per_line) {
  FILE *f = fopen(filename_out.c_str(), "w");
  // now print data as CSV to file
  // -- header
  bool first = true;
  for (const auto &l: labels) {
    if (!first) fprintf(f, ",");
    first = false;
    fprintf(f, "%s", l.c_str());
  }
  fprintf(f, "\n");
  // -- data

  printf("Formatting the timestamps.\n");
  std::chrono::steady_clock::time_point chrono_begin = std::chrono::steady_clock::now();
  std::vector<std::array<char, 18>> curTimes(num_samples_per_line);

  int num_threads_formatting = (num_samples_per_line + NUM_SAMPLES_PER_THREAD_FORMATTING - 1) / NUM_SAMPLES_PER_THREAD_FORMATTING;
  #pragma omp parallel for
  for (size_t thread_id = 0; thread_id < num_threads_formatting; ++thread_id) {
    size_t max_sample_index = (thread_id+1)*NUM_SAMPLES_PER_THREAD_FORMATTING;
    if (max_sample_index > num_samples_per_line)
      max_sample_index = num_samples_per_line;
    for (size_t sample_no = thread_id*NUM_SAMPLES_PER_THREAD_FORMATTING; sample_no < max_sample_index; ++sample_no) {
      sprintf(curTimes[sample_no].data(), "%.8e,", hor_start + hor_step_size*sample_no);
    }
  }

  std::chrono::steady_clock::time_point chrono_end = std::chrono::steady_clock::now();
  std::cout << "Formatting the timestamps done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  printf("Generating the data lines.\n");
  chrono_begin = std::chrono::steady_clock::now();
  std::vector<std::array<char,40>> datalines(num_samples_per_line);

  #pragma omp parallel for
  for (size_t thread_id = 0; thread_id < num_threads_formatting; ++thread_id) {
    size_t cur_line;
    size_t cur_offset_in_line;
    size_t max_row_index = (thread_id+1)*NUM_SAMPLES_PER_THREAD_FORMATTING;
    if (max_row_index > num_samples_per_line)
      max_row_index = num_samples_per_line;
    for (size_t cur_row = thread_id*NUM_SAMPLES_PER_THREAD_FORMATTING; cur_row < max_row_index; ++cur_row) {
      for (size_t cur_sample_in_row = 0; cur_sample_in_row < num_lines-1; ++cur_sample_in_row) {
        datalines[cur_row][2*cur_sample_in_row] = '0'+all_samples_per_line.at(cur_sample_in_row).at(cur_row);      
        datalines[cur_row][2*cur_sample_in_row+1] = ',';
      }
      size_t cur_sample_in_row = num_lines-1;
      datalines[cur_row][2*cur_sample_in_row] = '0'+all_samples_per_line.at(cur_sample_in_row).at(cur_row);      
      datalines[cur_row][2*cur_sample_in_row+1] = '\n';      
      datalines[cur_row][2*cur_sample_in_row+2] = '\0';      
    }
  }

  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Generating the data lines done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  printf("Printing the data.\n");
  chrono_begin = std::chrono::steady_clock::now();

  for (size_t sample_no = 0; sample_no < num_samples_per_line; ++sample_no) {
    fputs(curTimes[sample_no].data(), f);
    fputs(datalines[sample_no].data(), f);
  }
  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Printing the data done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;
  fclose(f);
}
