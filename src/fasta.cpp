// FASTA parsing and simple sequence utilities.
#include <Rcpp.h>
#include <string>
#include <vector>
#include <cctype>

using namespace Rcpp;

// Parse the raw bytes of a (possibly multi-line) FASTA file.
// Returns a list with the header lines (without '>') and the sequences
// (upper case, whitespace removed). Lines starting with ';' are ignored.
// [[Rcpp::export(rng = false)]]
List parse_fasta_raw(RawVector x) {
  const unsigned char *p = RAW(x);
  const size_t n = Rf_xlength(x);
  std::vector<std::string> headers;
  std::vector<std::string> seqs;
  std::string cur;
  bool in_record = false;
  size_t i = 0;
  // skip UTF-8 byte order mark
  if (n >= 3 && p[0] == 0xEF && p[1] == 0xBB && p[2] == 0xBF) i = 3;
  while (i < n) {
    size_t j = i;
    while (j < n && p[j] != '\n') j++;
    size_t end = j;
    if (end > i && p[end - 1] == '\r') end--;
    if (end > i && p[i] == '>') {
      if (in_record) seqs.push_back(cur);
      cur.clear();
      headers.emplace_back(reinterpret_cast<const char *>(p + i + 1), end - i - 1);
      in_record = true;
    } else if (end > i && p[i] == ';') {
      // comment line: ignore
    } else {
      for (size_t k = i; k < end; k++) {
        unsigned char c = p[k];
        if (std::isspace(c)) continue;
        if (!in_record) {
          stop("This does not look like a FASTA file: sequence data found before the first '>' header line.");
        }
        cur.push_back(static_cast<char>(std::toupper(c)));
      }
    }
    i = j + 1;
  }
  if (in_record) seqs.push_back(cur);
  return List::create(_["header"] = wrap(headers), _["sequence"] = wrap(seqs));
}

static inline char complement(char c) {
  switch (c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';
  case 'U': return 'A';
  case 'R': return 'Y';
  case 'Y': return 'R';
  case 'S': return 'S';
  case 'W': return 'W';
  case 'K': return 'M';
  case 'M': return 'K';
  case 'B': return 'V';
  case 'V': return 'B';
  case 'D': return 'H';
  case 'H': return 'D';
  case 'N': return 'N';
  case 'a': return 't';
  case 'c': return 'g';
  case 'g': return 'c';
  case 't': return 'a';
  case 'u': return 'a';
  default: return c;
  }
}

// Reverse complement of each element (IUPAC aware).
// [[Rcpp::export(rng = false)]]
CharacterVector revcomp_cpp(CharacterVector x) {
  const R_xlen_t n = x.size();
  CharacterVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    if (CharacterVector::is_na(x[i])) {
      out[i] = NA_STRING;
      continue;
    }
    std::string s = Rcpp::as<std::string>(x[i]);
    std::string r(s.size(), 'N');
    const size_t m = s.size();
    for (size_t k = 0; k < m; k++) r[m - 1 - k] = complement(s[k]);
    out[i] = r;
  }
  return out;
}
