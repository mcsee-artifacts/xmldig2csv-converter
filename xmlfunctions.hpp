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


void get_string(xmlDocPtr doc, xmlXPathContextPtr xpathCtx, const xmlChar *xpathExpr, unsigned char *outstr) {
  auto xpathObj = xmlXPathEvalExpression(BAD_CAST xpathExpr, xpathCtx);
  if (xpathObj!=nullptr) {
    auto nodeset = xpathObj->nodesetval;
    for (auto i = 0; i < nodeset->nodeNr; i++) {
      auto keyword = xmlNodeListGetString(doc, nodeset->nodeTab[i]->xmlChildrenNode, 1);
      sprintf(reinterpret_cast<char *>(outstr), "%s", keyword);
      xmlFree(keyword);
    }
  } else {
    fprintf(stderr, "Error: unable to evaluate xpath expression \"%s\"\n", xpathExpr);
    xmlXPathFreeObject(xpathObj);
    exit(EXIT_FAILURE);
  }
  xmlXPathFreeObject(xpathObj);
}

long evalXpathGetLong(xmlDocPtr doc, xmlXPathContextPtr xpathCtx, xmlChar *xpePathExpr) {
  const auto outstr_numChars = 256;
  char outstr[outstr_numChars];
  char *endptr;
  get_string(doc, xpathCtx, xpePathExpr, reinterpret_cast<unsigned char *>(outstr));
  auto result = strtol(outstr, &endptr, 10);
  if (endptr==outstr || ((result==LONG_MAX || result==LONG_MIN) && errno==ERANGE)) {
    std::cerr << "[-] could not eval XPath " << xpePathExpr << std::endl;
    exit(EXIT_FAILURE);
  }
  return result;
}

float evalXpathGetFloat(xmlDocPtr doc, xmlXPathContextPtr xpathCtx, xmlChar *xpePathExpr) {
  const auto outstr_numChars = 256;
  char outstr[outstr_numChars];
  char *endptr;
  get_string(doc, xpathCtx, xpePathExpr, reinterpret_cast<unsigned char *>(outstr));
  auto result = strtof(outstr, &endptr);
  if (endptr==outstr
      || ((result==std::numeric_limits<float>::max() || result==std::numeric_limits<float>::min()) && errno==ERANGE)) {
    std::cerr << "[-] could not parse HorStart from XML" << std::endl;
    exit(EXIT_FAILURE);
  }
  return result;
}

inline bool file_exists(const std::string& name) {
  struct stat buffer{};
  return (stat (name.c_str(), &buffer) == 0);
}

std::string evalXpathGetString(xmlDocPtr doc, xmlXPathContextPtr xpathCtx, xmlChar *xpePathExpr, size_t expLength) {
  char *out_buf = (char*)malloc(expLength);
  get_string(doc, xpathCtx, xpePathExpr, reinterpret_cast<unsigned char *>(out_buf));
  std::string out(out_buf);
  free(out_buf);
  return out;
}

void evalXpathGetString(xmlDocPtr doc, xmlXPathContextPtr xpathCtx, xmlChar *xpePathExpr, char* buf) {
  get_string(doc, xpathCtx, xpePathExpr, reinterpret_cast<unsigned char *>(buf));
}

void replace_first(std::string &s, std::string const &toReplace, std::string const &replaceWith) {
  std::size_t pos = s.find(toReplace);
  if (pos==std::string::npos)
    return;
  s.replace(pos, toReplace.length(), replaceWith);
}
