#include "third_party/argparse.hpp"
#include "write_csv.hpp"
#include "xmlfunctions.hpp"

#include <fcntl.h>
#include <sys/mman.h>

#define XMLDIG2CSV_HUGE_SUPPORT 1

// This function returns two things (using out parameters):
// * an in-memory XMLdig file with an empty "BinaryData" tag (used for extracting metadata using an XML parser), and
// * a pointer (+ length) to the "BinaryData".
void extract_binary_data(char const* filename,
                         char *& out_xmldig_nocontent_ptr,size_t& out_xmldig_nocontent_len,
                         char const*& out_binary_data_ptr, size_t& out_binary_data_len) {
  int fd = open(filename, O_RDONLY);
  auto full_len = std::experimental::filesystem::file_size(filename);
  auto* full_ptr = (char const*)mmap(nullptr, full_len, PROT_READ, MAP_SHARED, fd, 0);
  char const* DATA_START_TAG = "<BinaryData>";
  char const* DATA_END_TAG = "</BinaryData>";
  // Find start.
  auto* start_pos = (char const*)memmem(full_ptr, full_len, DATA_START_TAG, strlen(DATA_START_TAG));
  out_binary_data_ptr = start_pos + strlen(DATA_START_TAG);
  // Find end.
  auto* end_pos = (char const*)memmem(full_ptr, full_len, DATA_END_TAG, strlen(DATA_END_TAG));
  out_binary_data_len = end_pos - out_binary_data_ptr;

  printf("Found binary data in XMLdig file: %zu bytes.\n", out_binary_data_len);

  // Create new in-memory copy of the XMLdig file, without the binary data.
  out_xmldig_nocontent_len = full_len - out_binary_data_len;
  out_xmldig_nocontent_ptr = (char*)malloc(out_xmldig_nocontent_len);
  auto before_content_len = (size_t)(out_binary_data_ptr - full_ptr);
  auto after_content_len = (size_t)(out_xmldig_nocontent_len - before_content_len);
  assert(before_content_len + after_content_len == out_xmldig_nocontent_len);
  memcpy(out_xmldig_nocontent_ptr, full_ptr, before_content_len);
  memcpy(out_xmldig_nocontent_ptr + before_content_len, full_ptr + before_content_len + out_binary_data_len, after_content_len);

  auto binary_data_mib = (double)out_binary_data_len / 1024 / 1024;
  printf("[+] Extracted binary data from XML tree (content-less XML size is %zu B, binary data size is %zu B (%.1f MiB)).\n",
         out_xmldig_nocontent_len, out_binary_data_len, binary_data_mib);
}

int main(int argc, char **argv) {
  argparse::ArgumentParser program("xmldig2csv");

  program.add_argument("path")
    .help("display the square of a given integer")
    .required();

  program.add_argument("--no-trimming")
    .help("trim the output to the first N samples")
    .implicit_value(true)
    .default_value(false);

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  std::string filename_in = program.get<std::string>("path");
  bool trim = program.get<bool>("--no-trimming");

  // Check that filename_in is a file that exists.
  if (!file_exists(filename_in)) {
    std::cerr << "[-] Input file " << filename_in << " does not exist!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Define the output filename.
  std::string filename_out(filename_in);
  if (filename_out.find(".XMLdig") == std::string::npos) {
    std::cerr << "[-] input file must end with .XMLdig (same capitalization!)" << std::endl;
  }
  replace_first(filename_out, ".XMLdig", ".csv");

  // do not do anything if the output file already exists
  if (file_exists(filename_out)) {
    std::cout << "[+] destination file (" << filename_out << ") already exists, skipping conversion." << std::endl;
    exit(EXIT_SUCCESS);
  }

  printf("[+] Parsing the XMLdig file: %s\n", filename_in.c_str());
  std::chrono::steady_clock::time_point chrono_begin = std::chrono::steady_clock::now();

  // === PARSING ===
  xmlInitParser();

#if XMLDIG2CSV_HUGE_SUPPORT
  char* xmldig_nocontent_ptr;
  size_t xmldig_nocontent_len;
  char const* binary_data_ptr;
  size_t binary_data_len;
  extract_binary_data(filename_in.c_str(), xmldig_nocontent_ptr, xmldig_nocontent_len, binary_data_ptr, binary_data_len);
  auto doc = xmlReadMemory(xmldig_nocontent_ptr, (int)xmldig_nocontent_len, nullptr, nullptr, XML_PARSE_HUGE);
#else
  auto doc = xmlReadFile(filename_in, nullptr, XML_PARSE_HUGE);
#endif
  if (doc == nullptr) {
    fprintf(stderr, "[-] Unable to parse file \"%s\"\n", filename_in.c_str());
    xmlFreeDoc(doc);
    exit(EXIT_FAILURE);
  }
  std::chrono::steady_clock::time_point chrono_end = std::chrono::steady_clock::now();
  std::cout << "[+] Parsing done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  auto xpathCtx = xmlXPathNewContext(doc);
  if (xpathCtx==nullptr) {
    fprintf(stderr, "Error: unable to create new XPath context\n");
    xmlFreeDoc(doc);
    exit(EXIT_FAILURE);
  }

  auto *xpeNumLines = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/BusInfo/NumLines";
  auto num_lines = static_cast<int>(evalXpathGetLong(doc, xpathCtx, xpeNumLines));

  auto *xpNumSamples = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/NumSamples";
  auto num_samples_per_line = evalXpathGetLong(doc, xpathCtx, xpNumSamples); // WARNING: A line is a column!

  auto *xpHorStart = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/HorScale/HorStart";
  auto hor_start = evalXpathGetFloat(doc, xpathCtx, xpHorStart);

  auto *xpHorStep = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/HorScale/HorPerStep";
  auto hor_step = evalXpathGetFloat(doc, xpathCtx, xpHorStep);

  std::stringstream ss;
  std::vector<std::string> labels = {"Time"};
  for (int i = 1; i < num_lines + 1; ++i) {
    std::stringstream sstream;
    sstream << "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/LineInfo/Line[" << i << "]/Name";
    auto label = evalXpathGetString(doc, xpathCtx, (xmlChar *) sstream.str().c_str(), 256);
    labels.push_back(label);
    ss << label << ",";
  }
  std::cout << "[+] Found columns: " << ss.str() << "\n";

#if XMLDIG2CSV_HUGE_SUPPORT
  // FIXME: Due to the buffer overrun described below, we need to provide a larger buffer to avoid a segfault.
  char* binaryData = (char*)malloc(binary_data_len + 4096);
  memcpy(binaryData, binary_data_ptr, binary_data_len);
  size_t num_ascii_chars_input = binary_data_len;
#else
  auto *xpBinaryData = (xmlChar *) "//LeCroyXStreamDSOdata/Data/Segment/BinaryData";
  std::uintmax_t size = std::experimental::filesystem::file_size(filename_in)*64;

  printf("Zeroing binarydata.\n");
  chrono_begin = std::chrono::steady_clock::now();
  char *binaryData = (char*)malloc(size);
  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Zeroing binarydata in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  evalXpathGetString(doc, xpathCtx, xpBinaryData, binaryData);
  const auto num_ascii_chars_input = strlen(binaryData);
#endif

  xmlXPathFreeContext(xpathCtx);
  xmlFreeDoc(doc);
  xmlCleanupParser();

  const auto MAX_WINDOW_SIZE = 8; // ASCII characters
  const auto BITS_PER_ASCII_CHAR = 6;
  const auto BITS_PER_DIGITAL_VAL = 8;
  const auto NUM_SAMPLES_PER_MAX_WINDOW = MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR / BITS_PER_DIGITAL_VAL;

  assert(MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR >= 8 && MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR <= 64);
  assert(MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR % BITS_PER_DIGITAL_VAL == 0);
  // std::vector<size_t> all_samples_per_line_flat(num_lines*num_samples_per_line, 0);
  std::vector<std::vector<size_t>> all_samples_per_line_transpose(num_samples_per_line, std::vector<size_t>(num_lines, 0));

  // For parallelization.
  // We use a separate thread to treat the trailing outputs (mayube it's not extremely efficient but should be enough for now)
  size_t num_windows = (num_ascii_chars_input+MAX_WINDOW_SIZE-1) / MAX_WINDOW_SIZE;
  size_t num_threads = (num_windows+NUM_WINDOWS_PER_THREAD-1) / NUM_WINDOWS_PER_THREAD;

  omp_set_num_threads(8); // Use 16 threads for all consecutive parallel regions

  printf("[+] Starting the pre-decoding.\n");
  chrono_begin = std::chrono::steady_clock::now();

  // FIXME: This does currently not handle the end of binaryData correctly, in two ways:
  //        (1) It overruns the buffer at the end while reading.
  //        (2) Depending on trailing whitespace (e.g., CRLF vs LF), it produces very slightly different timestamp
  //            values for the first few (~3) samples.
  #pragma omp parallel for default(none) shared(num_threads,num_samples_per_line,binaryData,num_lines,all_samples_per_line_transpose)
  for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
    for (size_t window_id_in_thread = 0; window_id_in_thread < NUM_WINDOWS_PER_THREAD; window_id_in_thread++) {
        // int cur_window_size;
        // printf("Thread: %ld, Window: %ld\n", thread_id, window_id_in_thread); 

        uint64_t decoded_values = 0;
        size_t window_id_absolute = thread_id*NUM_WINDOWS_PER_THREAD+window_id_in_thread;

        // // Check whether the window is complete, or if it is the last one.
        // if ((window_id_absolute+1)*MAX_WINDOW_SIZE > num_ascii_chars_input) {
        //   int remaining_samples = (num_ascii_chars_input - (window_id_absolute+1)*MAX_WINDOW_SIZE) * BITS_PER_ASCII_CHAR / BITS_PER_DIGITAL_VAL;
        //   cur_window_size = get_next_bigger(remaining_samples/BITS_PER_ASCII_CHAR, BITS_PER_ASCII_CHAR);
        // } else {
        //   cur_window_size = MAX_WINDOW_SIZE;
        // }

        for (size_t char_id_in_window = 0; char_id_in_window < MAX_WINDOW_SIZE; char_id_in_window++) {
          size_t char_id_absolute = window_id_absolute*MAX_WINDOW_SIZE + char_id_in_window;
          auto base64_value = binaryData[char_id_absolute];
          auto decoded_value = BASE64_TO_INT(base64_value).to_ullong();
          auto shift_offt = BITS_PER_ASCII_CHAR*(MAX_WINDOW_SIZE - (char_id_in_window + 1));
          auto decoded_value_shifted = (decoded_value << shift_offt);
          decoded_values |= decoded_value_shifted;
        }

        // interpret decoded values: each 8 bits represent one digital value (i.e., 0 or 1)
        for (size_t sample_id_in_window = 0; sample_id_in_window < NUM_SAMPLES_PER_MAX_WINDOW; ++sample_id_in_window) {
          // Util
          int shift_offt_down = (MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR)-(BITS_PER_DIGITAL_VAL * (sample_id_in_window+1));
          uint64_t mask = ((uint64_t)0xffu) << ((MAX_WINDOW_SIZE*BITS_PER_ASCII_CHAR)-(BITS_PER_DIGITAL_VAL * (sample_id_in_window+1)));

          uint8_t extracted_num = (((decoded_values & mask) >> shift_offt_down));
          size_t sample_id_absolute = window_id_absolute*NUM_SAMPLES_PER_MAX_WINDOW + sample_id_in_window;

          auto cur_line = static_cast<int>(sample_id_absolute / num_samples_per_line);
          if (cur_line >= num_lines)
            break;
          size_t cur_offset_in_line = sample_id_absolute % num_samples_per_line;

          all_samples_per_line_transpose[cur_offset_in_line][cur_line] = extracted_num;
          // all_samples_per_line_flat[cur_line*num_samples_per_line + cur_offset_in_line] = extracted_num;
      }
    }
  }

  chrono_end = std::chrono::steady_clock::now();
  std::cout << "[+] Pre-decoding done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  // std::vector<std::vector<size_t>> all_samples_per_line(num_lines, std::vector<size_t>(num_samples_per_line, 0));

  write_csv_fast(filename_out, hor_start, hor_step, labels, all_samples_per_line_transpose, !trim);

#if XMLDIG2CSV_HUGE_SUPPORT
  // Do nothing.
#else
  printf("Freeing binary data.\n");
  chrono_begin = std::chrono::steady_clock::now();
  free(binaryData);
  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Freeing binary data done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;
#endif
  return EXIT_SUCCESS;
}
