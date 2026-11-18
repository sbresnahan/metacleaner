// Native nucleotide search engine for metacleaner.
//
// Reference sequences are indexed by k-mers (2-bit encoded). For each query
// (both strands) candidate references are found with a lossless q-gram
// filter: given the identity and coverage thresholds, any reference that can
// produce a qualifying alignment must share at least T k-mers with the query.
// The most frequent query k-mers are skipped (at most T - 1 of them), which
// keeps the filter lossless while avoiding very common (conserved) k-mers.
// Candidates are verified with a banded local alignment using
// megablast-like scores (match +1, mismatch -2, linear gap 2.5; scaled by 2).
// When both thresholds are 100 an exact containment test is used instead.

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <thread>
#include <atomic>

using namespace Rcpp;

namespace {

const int SC_MATCH = 2;
const int SC_MISMATCH = -4;
const int SC_GAP = -5;

inline uint8_t enc(char c) {
  switch (c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': case 'U': case 'u': return 3;
  default: return 4;
  }
}

inline char upc(char c) { return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c; }

std::string revcomp(const std::string &s) {
  std::string r(s.size(), 'N');
  const size_t m = s.size();
  for (size_t k = 0; k < m; k++) {
    char c = upc(s[k]), o;
    switch (c) {
    case 'A': o = 'T'; break;
    case 'C': o = 'G'; break;
    case 'G': o = 'C'; break;
    case 'T': case 'U': o = 'A'; break;
    case 'R': o = 'Y'; break;
    case 'Y': o = 'R'; break;
    case 'K': o = 'M'; break;
    case 'M': o = 'K'; break;
    case 'B': o = 'V'; break;
    case 'V': o = 'B'; break;
    case 'D': o = 'H'; break;
    case 'H': o = 'D'; break;
    default: o = c;
    }
    r[m - 1 - k] = o;
  }
  return r;
}

struct KmerIndex {
  int k;
  uint32_t mask;
  std::vector<uint32_t> offsets;   // 4^k + 1
  std::vector<uint32_t> positions; // global coordinates
  std::vector<uint32_t> starts;    // subject start coordinates (n + 1)

  inline uint32_t occ(uint32_t code) const { return offsets[code + 1] - offsets[code]; }
  inline uint32_t subject_of(uint32_t g) const {
    return static_cast<uint32_t>(std::upper_bound(starts.begin(), starts.end(), g) - starts.begin() - 1);
  }
};

void build_index(const std::vector<std::string> &subj, int k, KmerIndex &idx) {
  idx.k = k;
  const uint64_t nb = 1ULL << (2 * k);
  idx.mask = static_cast<uint32_t>(nb - 1);
  idx.starts.assign(subj.size() + 1, 0);
  uint64_t total = 0;
  for (size_t s = 0; s < subj.size(); s++) {
    idx.starts[s] = static_cast<uint32_t>(total);
    total += subj[s].size();
    if (total >= std::numeric_limits<uint32_t>::max())
      stop("The reference database is too large for the native engine (> 4.29 billion bases). Split it into parts or use engine = \"blast\".");
  }
  idx.starts[subj.size()] = static_cast<uint32_t>(total);
  std::vector<uint32_t> counts(nb + 1, 0);
  for (size_t s = 0; s < subj.size(); s++) {
    const std::string &x = subj[s];
    uint32_t code = 0;
    int run = 0;
    for (size_t i = 0; i < x.size(); i++) {
      uint8_t e = enc(x[i]);
      if (e > 3) { run = 0; code = 0; continue; }
      code = ((code << 2) | e) & idx.mask;
      if (++run >= k) counts[code + 1]++;
    }
  }
  for (uint64_t c = 0; c < nb; c++) counts[c + 1] += counts[c];
  idx.offsets = counts;
  idx.positions.assign(idx.offsets[nb], 0);
  std::vector<uint32_t> cursor(idx.offsets.begin(), idx.offsets.end() - 1);
  for (size_t s = 0; s < subj.size(); s++) {
    const std::string &x = subj[s];
    uint32_t code = 0;
    int run = 0;
    for (size_t i = 0; i < x.size(); i++) {
      uint8_t e = enc(x[i]);
      if (e > 3) { run = 0; code = 0; continue; }
      code = ((code << 2) | e) & idx.mask;
      if (++run >= k) idx.positions[cursor[code]++] = idx.starts[s] + static_cast<uint32_t>(i + 1 - k);
    }
  }
}

struct Hit {
  int subject = -1;
  double identity = 0;
  double coverage = 0;
  int score = 0;
  int qstart = 0, qend = 0, sstart = 0, send = 0; // 0-based inclusive, query on + strand
  bool minus = false;
};

inline bool hit_better(const Hit &a, const Hit &b) {
  if (a.score != b.score) return a.score > b.score;
  if (a.identity != b.identity) return a.identity > b.identity;
  return a.subject < b.subject;
}

// Per-thread scratch space, reused across queries.
struct Workspace {
  std::vector<uint32_t> cnt;      // shared k-mer count per subject
  std::vector<uint32_t> touched;  // subjects with cnt > 0
  std::vector<int> prev, curr;
  std::vector<uint8_t> tb;
  explicit Workspace(size_t nsub) : cnt(nsub, 0) {}
};

struct QueryResult {
  std::vector<Hit> pass;
  int n_pass = 0;
  bool has_best_fail = false;
  Hit best_fail;
};

// Banded local alignment of query q against subject s around diagonal d
// (subject position = query position + d), half band width w.
bool banded_local(const std::string &q, const std::string &s, long d, int w,
                  std::vector<int> &prev, std::vector<int> &curr,
                  std::vector<uint8_t> &tb, Hit &out) {
  const long L = static_cast<long>(q.size());
  const long S = static_cast<long>(s.size());
  const int W = 2 * w + 1;
  prev.assign(W + 2, 0);
  curr.assign(W + 2, 0);
  tb.assign(static_cast<size_t>(L) * W, 0);
  int best = 0;
  long bi = -1, bo = -1;
  // offset o in [0, W): j = i + d - w + o ; arrays are shifted by 1 so that
  // index o+1 holds offset o and indices 0 and W+1 act as zero borders.
  // Sequences are upper case with U converted to T (done in R), so two
  // positions are identical exactly when their characters are equal.
  int *pv = prev.data();
  int *cv = curr.data();
  const char *sp = s.data();
  for (long i = 0; i < L; i++) {
    std::fill(cv, cv + W + 2, 0);
    const long jbase = i + d - w;
    const int olo = static_cast<int>(std::max<long>(0, -jbase));
    const int ohi = static_cast<int>(std::min<long>(W - 1, S - 1 - jbase));
    const char qc = q[i];
    uint8_t *tbr = tb.data() + static_cast<size_t>(i) * W;
    for (int o = olo; o <= ohi; o++) {
      const bool same = (qc == sp[jbase + o]);
      const int hd = pv[o + 1] + (same ? SC_MATCH : SC_MISMATCH); // diagonal (i-1, j-1)
      const int hu = pv[o + 2] + SC_GAP;                          // up (i-1, j): gap in subject
      const int hl = cv[o] + SC_GAP;                              // left (i, j-1): gap in query
      int h = 0;
      uint8_t dir = 0;
      if (hd > h) { h = hd; dir = same ? 1 : 4; }
      if (hu > h) { h = hu; dir = 2; }
      if (hl > h) { h = hl; dir = 3; }
      cv[o + 1] = h;
      tbr[o] = dir;
      if (h > best) { best = h; bi = i; bo = o; }
    }
    std::swap(pv, cv);
  }
  if (best <= 0) return false;
  // traceback
  long i = bi;
  int o = static_cast<int>(bo);
  long cols = 0, matches = 0;
  long qe = bi, se = bi + d - w + bo, qs = bi, ss = se;
  while (i >= 0 && o >= 0 && o < W) {
    uint8_t dir = tb[static_cast<size_t>(i) * W + o];
    if (dir == 0) break;
    long j = i + d - w + o;
    qs = i;
    ss = j;
    cols++;
    if (dir == 1) { matches++; i--; }
    else if (dir == 4) { i--; }
    else if (dir == 2) { i--; o++; }
    else { o--; }
  }
  out.score = best;
  out.identity = cols > 0 ? 100.0 * static_cast<double>(matches) / static_cast<double>(cols) : 0.0;
  out.coverage = 100.0 * static_cast<double>(qe - qs + 1) / static_cast<double>(L);
  out.qstart = static_cast<int>(qs);
  out.qend = static_cast<int>(qe);
  out.sstart = static_cast<int>(ss);
  out.send = static_cast<int>(se);
  return true;
}

struct Params {
  double identity, coverage;
  int max_hits;
  bool best_fail;
  int occ_skip; // k-mers with more occurrences than this may be skipped
  int heuristic_candidates;
  int min_score;
};

inline bool passes(const Hit &h, const Params &p) {
  return h.identity + 1e-9 >= p.identity && h.coverage + 1e-9 >= p.coverage;
}

// Search one strand of one query. Adds hits (best per subject) into 'best'.
void search_strand(const std::string &q, bool minus, const KmerIndex &idx,
                   const std::vector<std::string> &subj, const Params &p,
                   int exclude, std::vector<Hit> &hits, Workspace &ws,
                   bool heuristic_pass) {
  std::vector<int> &prev = ws.prev;
  std::vector<int> &curr = ws.curr;
  std::vector<uint8_t> &tb = ws.tb;
  const int k = idx.k;
  const long L = static_cast<long>(q.size());
  if (L < k) return;
  // collect valid k-mers
  std::vector<std::pair<uint32_t, int>> km; // (code, qpos)
  km.reserve(L);
  uint32_t code = 0;
  int run = 0;
  for (long i = 0; i < L; i++) {
    uint8_t e = enc(q[i]);
    if (e > 3) { run = 0; code = 0; continue; }
    code = ((code << 2) | e) & idx.mask;
    if (++run >= k) km.emplace_back(code, static_cast<int>(i + 1 - k));
  }
  const long n = static_cast<long>(km.size());
  if (n == 0) return;
  const long windows = L - k + 1;
  const long invalid = windows - n;

  const bool exact = (p.identity >= 100.0 && p.coverage >= 100.0 && !heuristic_pass);

  // order k-mers by frequency (rarest first)
  std::sort(km.begin(), km.end(), [&idx](const std::pair<uint32_t, int> &a, const std::pair<uint32_t, int> &b) {
    uint32_t oa = idx.occ(a.first), ob = idx.occ(b.first);
    if (oa != ob) return oa < ob;
    return a.second < b.second;
  });

  if (exact) {
    const uint32_t c = km[0].first;
    const int qp = km[0].second;
    for (uint32_t t = idx.offsets[c]; t < idx.offsets[c + 1]; t++) {
      uint32_t g = idx.positions[t];
      uint32_t s = idx.subject_of(g);
      if (static_cast<int>(s) == exclude) continue;
      long spos = static_cast<long>(g - idx.starts[s]);
      long start = spos - qp;
      const std::string &x = subj[s];
      if (start < 0 || start + L > static_cast<long>(x.size())) continue;
      bool ok = true;
      for (long i = 0; i < L; i++) {
        if (q[i] != x[start + i]) { ok = false; break; }
      }
      if (!ok) continue;
      Hit h;
      h.subject = static_cast<int>(s);
      h.identity = 100.0;
      h.coverage = 100.0;
      h.score = static_cast<int>(SC_MATCH * L);
      h.qstart = 0; h.qend = static_cast<int>(L - 1);
      h.sstart = static_cast<int>(start); h.send = static_cast<int>(start + L - 1);
      h.minus = minus;
      hits.push_back(h);
    }
    return;
  }

  // minimum number of shared k-mers (q-gram lemma), minimised over the
  // possible aligned query lengths
  long need = 1;
  long skip = 0;
  int w;
  const double t = p.identity / 100.0;
  if (!heuristic_pass && t > 0) {
    long minLc = static_cast<long>(std::ceil(p.coverage / 100.0 * L - 1e-9));
    if (minLc < 1) minLc = 1;
    long Tmin = std::numeric_limits<long>::max();
    for (long Lc = minLc; Lc <= L; Lc++) {
      long e = static_cast<long>(std::floor(Lc * (1.0 - t) / t + 1e-9));
      long T = (Lc - k + 1) - static_cast<long>(k) * e - invalid;
      if (T < Tmin) Tmin = T;
    }
    if (Tmin < 1) Tmin = 1;
    // Lossless filter: at least n - (Tmin - 1) of the rarest k-mers must be
    // examined. Beyond that, all k-mers that are not very frequent are also
    // examined, which raises the required count and removes candidates.
    long processed = std::max<long>(1, n - (Tmin - 1));
    while (processed < n && static_cast<long>(idx.occ(km[processed].first)) <= p.occ_skip) processed++;
    skip = n - processed;
    need = Tmin - skip;
    long emax = static_cast<long>(std::floor(L * (1.0 - t) / t + 1e-9));
    w = static_cast<int>(std::min<long>(emax + 8, L));
  } else {
    // heuristic search for the best available hit: ignore very frequent k-mers
    long nfreq = 0;
    for (long i = n - 1; i >= 0; i--) {
      if (static_cast<long>(idx.occ(km[i].first)) > p.occ_skip * 8L) nfreq++; else break;
    }
    skip = std::min(nfreq, n - 1);
    need = 1;
    w = static_cast<int>(std::min<long>(std::max<long>(16, L / 6), L));
  }
  if (w > 1000) w = 1000; // keep memory bounded for very long queries

  // pass 1: count shared k-mers per subject with a dense counter
  const long nproc = n - skip;
  for (long r = 0; r < nproc; r++) {
    const uint32_t c = km[r].first;
    for (uint32_t tt = idx.offsets[c]; tt < idx.offsets[c + 1]; tt++) {
      uint32_t sbj = idx.subject_of(idx.positions[tt]);
      if (static_cast<int>(sbj) == exclude) continue;
      if (ws.cnt[sbj] == 0) ws.touched.push_back(sbj);
      ws.cnt[sbj]++;
    }
  }
  if (ws.touched.empty()) return;
  std::vector<std::pair<uint32_t, uint32_t>> cc; // (count, subject) of candidates
  for (uint32_t sbj : ws.touched) {
    if (static_cast<long>(ws.cnt[sbj]) >= need) cc.emplace_back(ws.cnt[sbj], sbj);
    ws.cnt[sbj] = 0;
  }
  ws.touched.clear();
  if (cc.empty()) return;
  if (heuristic_pass && static_cast<int>(cc.size()) > p.heuristic_candidates) {
    std::partial_sort(cc.begin(), cc.begin() + p.heuristic_candidates, cc.end(),
                      [](const std::pair<uint32_t, uint32_t> &x, const std::pair<uint32_t, uint32_t> &y) {
                        if (x.first != y.first) return x.first > y.first;
                        return x.second < y.second;
                      });
    cc.resize(p.heuristic_candidates);
  }
  // pass 2: diagonals of the candidates only (cnt is used as a marker)
  for (const auto &x : cc) ws.cnt[x.second] = 1;
  std::vector<std::pair<uint32_t, int32_t>> sd; // (subject, diagonal)
  for (long r = 0; r < nproc; r++) {
    const uint32_t c = km[r].first;
    const int qp = km[r].second;
    for (uint32_t tt = idx.offsets[c]; tt < idx.offsets[c + 1]; tt++) {
      uint32_t g = idx.positions[tt];
      uint32_t sbj = idx.subject_of(g);
      if (!ws.cnt[sbj]) continue;
      sd.emplace_back(sbj, static_cast<int32_t>(static_cast<long>(g - idx.starts[sbj]) - qp));
    }
  }
  for (const auto &x : cc) ws.cnt[x.second] = 0;
  std::sort(sd.begin(), sd.end());
  // most populated diagonal per candidate
  std::vector<std::pair<long, std::pair<uint32_t, int32_t>>> cand;
  size_t a = 0;
  while (a < sd.size()) {
    size_t b = a;
    while (b < sd.size() && sd[b].first == sd[a].first) b++;
    size_t bestc = 0, cs = a;
    int32_t bd = sd[a].second;
    while (cs < b) {
      size_t ce = cs;
      while (ce < b && sd[ce].second == sd[cs].second) ce++;
      if (ce - cs > bestc) { bestc = ce - cs; bd = sd[cs].second; }
      cs = ce;
    }
    cand.push_back({static_cast<long>(b - a), {sd[a].first, bd}});
    a = b;
  }
  for (size_t ci = 0; ci < cand.size(); ci++) {
    const uint32_t s = cand[ci].second.first;
    const long d = cand[ci].second.second;
    Hit h;
    if (!banded_local(q, subj[s], d, w, prev, curr, tb, h)) continue;
    h.subject = static_cast<int>(s);
    h.minus = minus;
    if (minus) {
      int qs = static_cast<int>(L - 1 - h.qend), qe = static_cast<int>(L - 1 - h.qstart);
      h.qstart = qs;
      h.qend = qe;
    }
    hits.push_back(h);
  }
}

// keep the best hit per subject
void dedup_subjects(std::vector<Hit> &hits) {
  std::sort(hits.begin(), hits.end(), [](const Hit &x, const Hit &y) {
    if (x.subject != y.subject) return x.subject < y.subject;
    return hit_better(x, y);
  });
  std::vector<Hit> out;
  out.reserve(hits.size());
  for (size_t i = 0; i < hits.size(); i++) {
    if (i == 0 || hits[i].subject != hits[i - 1].subject) out.push_back(hits[i]);
  }
  hits.swap(out);
}

void search_query(const std::string &q, const KmerIndex &idx, const std::vector<std::string> &subj,
                  const Params &p, int exclude, QueryResult &res, Workspace &ws) {
  std::vector<Hit> hits;
  const std::string rc = revcomp(q);
  search_strand(q, false, idx, subj, p, exclude, hits, ws, false);
  search_strand(rc, true, idx, subj, p, exclude, hits, ws, false);
  dedup_subjects(hits);
  std::vector<Hit> pass, fail;
  for (const Hit &h : hits) {
    if (passes(h, p)) pass.push_back(h); else fail.push_back(h);
  }
  std::sort(pass.begin(), pass.end(), hit_better);
  res.n_pass = static_cast<int>(pass.size());
  if (static_cast<int>(pass.size()) > p.max_hits) pass.resize(p.max_hits);
  res.pass.swap(pass);
  if (p.best_fail && res.n_pass == 0) {
    if (fail.empty()) {
      std::vector<Hit> h2;
      search_strand(q, false, idx, subj, p, exclude, h2, ws, true);
      search_strand(rc, true, idx, subj, p, exclude, h2, ws, true);
      for (const Hit &h : h2) if (!passes(h, p)) fail.push_back(h);
    }
    if (!fail.empty()) {
      std::sort(fail.begin(), fail.end(), hit_better);
      if (fail[0].score >= p.min_score) {
        res.has_best_fail = true;
        res.best_fail = fail[0];
      }
    }
  }
}

} // namespace

// [[Rcpp::export(rng = false)]]
List native_search_cpp(CharacterVector query, CharacterVector subject,
                       double identity, double coverage, int k, int max_hits,
                       bool best_fail, IntegerVector exclude, int threads,
                       int occ_skip, int batch_size, bool verbose) {
  const R_xlen_t nq = query.size();
  const R_xlen_t ns = subject.size();
  std::vector<std::string> q(nq), s(ns);
  for (R_xlen_t i = 0; i < nq; i++) q[i] = Rcpp::as<std::string>(query[i]);
  for (R_xlen_t i = 0; i < ns; i++) s[i] = Rcpp::as<std::string>(subject[i]);
  std::vector<int> excl(nq, -1);
  if (exclude.size() == nq) {
    for (R_xlen_t i = 0; i < nq; i++) excl[i] = (exclude[i] == NA_INTEGER) ? -1 : exclude[i] - 1;
  }

  KmerIndex idx;
  build_index(s, k, idx);

  Params p;
  p.identity = identity;
  p.coverage = coverage;
  p.max_hits = max_hits;
  p.best_fail = best_fail;
  p.occ_skip = occ_skip;
  p.heuristic_candidates = 25;
  p.min_score = 2 * (k + 8);

  std::vector<QueryResult> res(nq);
  int nthreads = threads < 1 ? 1 : threads;
  std::atomic<bool> failed(false);
  for (R_xlen_t b0 = 0; b0 < nq; b0 += batch_size) {
    R_xlen_t b1 = std::min<R_xlen_t>(nq, b0 + batch_size);
    // Portable C++11 threads (no OpenMP needed, so this also works on macOS).
    // Queries are handed out one at a time from an atomic counter; workers
    // never touch the R API.
    std::atomic<R_xlen_t> next(b0);
    auto worker = [&]() {
      try {
        Workspace ws(s.size());
        for (;;) {
          R_xlen_t i = next.fetch_add(1);
          if (i >= b1) break;
          search_query(q[i], idx, s, p, excl[i], res[i], ws);
        }
      } catch (...) {
        failed = true;
      }
    };
    int nt = static_cast<int>(std::min<R_xlen_t>(nthreads, b1 - b0));
    if (nt <= 1) {
      worker();
    } else {
      std::vector<std::thread> pool;
      pool.reserve(nt);
      for (int t = 0; t < nt; t++) pool.emplace_back(worker);
      for (auto &th : pool) th.join();
    }
    if (failed) stop("The search ran out of memory or failed; try fewer threads or a smaller database.");
    Rcpp::checkUserInterrupt();
    if (verbose && nq > batch_size) {
      Rcpp::Rcout << "\r  searched " << b1 << " / " << nq << " sequences" << std::flush;
    }
  }
  if (verbose && nq > batch_size) Rcpp::Rcout << std::endl;

  // flatten
  size_t total = 0;
  for (R_xlen_t i = 0; i < nq; i++) total += res[i].pass.size() + (res[i].has_best_fail ? 1 : 0);
  IntegerVector oq(total), os(total), osc(total), oqs(total), oqe(total), oss(total), ose(total);
  NumericVector oid(total), ocov(total);
  CharacterVector ostr(total);
  LogicalVector opass(total);
  IntegerVector npass(nq);
  size_t r = 0;
  for (R_xlen_t i = 0; i < nq; i++) {
    npass[i] = res[i].n_pass;
    auto put = [&](const Hit &h, bool ps) {
      oq[r] = static_cast<int>(i + 1);
      os[r] = h.subject + 1;
      oid[r] = h.identity;
      ocov[r] = h.coverage;
      osc[r] = h.score;
      ostr[r] = h.minus ? "-" : "+";
      oqs[r] = h.qstart + 1;
      oqe[r] = h.qend + 1;
      oss[r] = h.sstart + 1;
      ose[r] = h.send + 1;
      opass[r] = ps;
      r++;
    };
    for (const Hit &h : res[i].pass) put(h, true);
    if (res[i].has_best_fail) put(res[i].best_fail, false);
  }
  DataFrame hits = DataFrame::create(
      _["query"] = oq, _["subject"] = os, _["identity"] = oid, _["coverage"] = ocov,
      _["score"] = osc, _["strand"] = ostr, _["qstart"] = oqs, _["qend"] = oqe,
      _["sstart"] = oss, _["send"] = ose, _["pass"] = opass,
      _["stringsAsFactors"] = false);
  return List::create(_["hits"] = hits, _["n_pass"] = npass);
}
