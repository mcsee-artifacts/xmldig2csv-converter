#include "write_csv.hpp"
#include "xmlfunctions.hpp"

int get_next_bigger(int number, int base) {
  return (number+base-1)/base;
}

int main(int argc, char **argv) {
  // make sure all mandatory arguments are present
  if (argc < 2) {
    std::cerr << "[-] positional argument with path to .XMLdig file is missing!" << std::endl;
    exit(EXIT_FAILURE);
  }
  auto filename_in = argv[1];

  // define the output filename
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

  printf("Parsing the XML file.\n");
  std::chrono::steady_clock::time_point chrono_begin = std::chrono::steady_clock::now();

  // parse the XML file
  xmlInitParser();
  auto doc = xmlReadFile(filename_in, nullptr, XML_PARSE_HUGE);
  if (doc == nullptr) {
    fprintf(stderr, "Error: unable to parse file \"%s\"\n", filename_in);
    xmlFreeDoc(doc);
    exit(EXIT_FAILURE);
  }
  std::chrono::steady_clock::time_point chrono_end = std::chrono::steady_clock::now();
  std::cout << "Parsing done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  auto xpathCtx = xmlXPathNewContext(doc);
  if (xpathCtx==nullptr) {
    fprintf(stderr, "Error: unable to create new XPath context\n");
    xmlFreeDoc(doc);
    exit(EXIT_FAILURE);
  }

  auto *xpeNumLines = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/BusInfo/NumLines";
  auto num_lines = evalXpathGetLong(doc, xpathCtx, xpeNumLines);

  auto *xpNumSamples = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/NumSamples";
  auto num_samples_per_line = evalXpathGetLong(doc, xpathCtx, xpNumSamples); // WARNING: A line is a column!

  auto *xpHorStart = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/HorScale/HorStart";
  auto hor_start = evalXpathGetFloat(doc, xpathCtx, xpHorStart);

  auto *xpHorStep = (xmlChar *) "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/HorScale/HorPerStep";
  auto hor_step = evalXpathGetFloat(doc, xpathCtx, xpHorStep);

  std::vector<std::string> labels = {"Time"};
  for (size_t i = 1; i < num_lines + 1; ++i) {
    std::stringstream sstream;
    sstream << "//LeCroyXStreamDSOdata/Header/ILecDigitalResult/LineInfo/Line[" << i << "]/Name";
    auto label = evalXpathGetString(doc, xpathCtx, (xmlChar *) sstream.str().c_str(), 256);
    labels.push_back(label);
  }

  auto *xpBinaryData = (xmlChar *) "//LeCroyXStreamDSOdata/Data/Segment/BinaryData";
  std::uintmax_t size = std::experimental::filesystem::file_size(filename_in)*64;

  printf("Zeroing binarydata.\n");
  chrono_begin = std::chrono::steady_clock::now();
  char *binaryData = (char*)malloc(size);
  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Zeroing binarydata in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  evalXpathGetString(doc, xpathCtx, xpBinaryData, binaryData);
  const auto num_ascii_chars_input = strlen(binaryData);

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
  std::vector<std::vector<size_t>> all_samples_per_line(num_lines, std::vector<size_t>(num_samples_per_line, 0));

// #ifdef DEBUG
//   std::cout << "num_samples: " << num_samples << std::endl;
//   std::cout << "num_lines: " << num_lines << std::endl;
// #endif

  // For parallelization.
  // We use a separate thread to treat the trailing outputs (mayube it's not extremely efficient but should be enough for now)
  size_t num_windows = (num_ascii_chars_input+MAX_WINDOW_SIZE-1) / MAX_WINDOW_SIZE;
  size_t num_threads = (num_windows+NUM_WINDOWS_PER_THREAD-1) / NUM_WINDOWS_PER_THREAD;

  omp_set_num_threads(16); // Use 16 threads for all consecutive parallel regions

  printf("Starting the decoding.\n");
  chrono_begin = std::chrono::steady_clock::now();

    #pragma omp parallel for
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

          size_t cur_line = sample_id_absolute / num_samples_per_line;
          if (cur_line >= num_lines)
            break;
          size_t cur_offset_in_line = sample_id_absolute % num_samples_per_line;

          all_samples_per_line[cur_line][cur_offset_in_line] = extracted_num;
          // all_samples_per_line_flat[cur_line*num_samples_per_line + cur_offset_in_line] = extracted_num;
      }
    }
  }

  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Decoding done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;

  // std::vector<std::vector<size_t>> all_samples_per_line(num_lines, std::vector<size_t>(num_samples_per_line, 0));

  write_csv_fast(filename_out, num_lines, num_samples_per_line, hor_start, hor_step, labels, all_samples_per_line);

  printf("Freeing binary data.\n");
  chrono_begin = std::chrono::steady_clock::now();
  free(binaryData);
  chrono_end = std::chrono::steady_clock::now();
  std::cout << "Freeing binary data done in " << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_end - chrono_begin).count() << "ms" << std::endl;
  return EXIT_SUCCESS;
}
