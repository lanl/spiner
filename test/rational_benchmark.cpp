// Copyright 2026 Triad National Security, LLC. All rights reserved.
// Generative AI was used to assist with this file.

// JMM: Over-engineered AI slop, but it'll do.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "rational_benchmark_metadata.hpp"
#include <ports-of-call/portability.hpp>
#include <spiner/databox.hpp>

namespace {
using RegularGrid = Spiner::RegularGrid1D<Real>;
using FastGrid = Spiner::FastNonUniformGrid1D<double>;
using Rational = Spiner::InterpolationStencils::Rational;
using Linear = Spiner::InterpolationStencils::Linear;
template <class Policy, class Grid>
using Box = Spiner::DataBox<Real, Grid, Policy>;
using Clock = std::chrono::steady_clock;
struct Options {
  int nx1d = 1024, nx2d = 256, ny2d = 192, queries = 1048576;
  int passes = 20, trials = 7, warmup = 3;
  std::uint64_t seed = 12345;
  std::string dim = "both", policy = "both", distribution = "interior-shuffled";
  std::string grid = "regular";
  std::string output = "rational-baseline.csv", hardware = "unspecified";
};
void usage() {
  std::cout
      << "rational_benchmark.bin [--dim 1|2|both] [--policy "
         "rational|linear|both]\n"
         "  [--grid regular|fast-nonuniform]\n"
         "  [--nx1d 1024] [--nx2d 256] [--ny2d 192] [--queries 1048576]\n"
         "  [--passes 20] [--trials 7] [--warmup 3] [--seed 12345]\n"
         "  [--distribution interior-shuffled|interior-ordered|boundary]\n"
         "  [--output rational-baseline.csv] [--hardware-label NAME]\n"
         "Times include query loads, interpolation, output stores, launch and "
         "fence.\n"
         "Checksums and validation are outside timing.\n"
         "The 'traffic' control copies query coordinates into output arrays:\n"
         "  1D writes (x, x); 2D writes (x, x, y). No grid lookup, table "
         "reads,\n"
         "  or interpolation. Its ns/call is time per processed query.\n"
         "It measures query/output traffic plus launch and fence overhead.\n"
         "Do not subtract it from interpolation timings.\n";
}
std::uint64_t number(const std::string &s) {
  if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
    throw std::runtime_error("Expected unsigned integer: " + s);
  std::size_t n;
  auto v = std::stoull(s, &n);
  if (n != s.size()) throw std::runtime_error("Invalid integer: " + s);
  return v;
}
Options parse(int argc, char **argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if (++i == argc) throw std::runtime_error("Missing value for " + key);
    const std::string v = argv[i];
    if (key == "--dim")
      o.dim = v;
    else if (key == "--grid")
      o.grid = v;
    else if (key == "--policy")
      o.policy = v;
    else if (key == "--distribution")
      o.distribution = v;
    else if (key == "--output")
      o.output = v;
    else if (key == "--hardware-label")
      o.hardware = v;
    else if (key == "--seed")
      o.seed = number(v);
    else {
      int *p = nullptr;
      if (key == "--nx1d")
        p = &o.nx1d;
      else if (key == "--nx2d")
        p = &o.nx2d;
      else if (key == "--ny2d")
        p = &o.ny2d;
      else if (key == "--queries")
        p = &o.queries;
      else if (key == "--passes")
        p = &o.passes;
      else if (key == "--trials")
        p = &o.trials;
      else if (key == "--warmup")
        p = &o.warmup;
      else
        throw std::runtime_error("Unknown option: " + key);
      auto n = number(v);
      if (n > std::uint64_t(std::numeric_limits<int>::max()))
        throw std::runtime_error("Option exceeds launch index limit: " + key);
      *p = static_cast<int>(n);
    }
  }
  if (o.nx1d < 4 || o.nx2d < 4 || o.ny2d < 4 || o.queries < 1 || o.passes < 1 ||
      o.trials < 1 || o.warmup < 0)
    throw std::runtime_error(
        "Require grids >=4, queries/passes/trials >=1, warmup >=0");
  if (o.dim != "1" && o.dim != "2" && o.dim != "both")
    throw std::runtime_error("Invalid dimension");
  if (o.policy != "rational" && o.policy != "linear" && o.policy != "both")
    throw std::runtime_error("Invalid policy");
  if (o.distribution != "interior-shuffled" &&
      o.distribution != "interior-ordered" && o.distribution != "boundary")
    throw std::runtime_error("Invalid distribution");
  if (o.grid != "regular" && o.grid != "fast-nonuniform")
    throw std::runtime_error("Invalid grid");
  if (o.output.empty()) throw std::runtime_error("Output path is empty");
  return o;
}
void fence() {
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::fence();
#elif defined(PORTABILITY_STRATEGY_CUDA)
  cudaDeviceSynchronize();
#endif
}
// SplitMix64: specified integer generator and explicit 53-bit fraction mapping.
struct Random {
  std::uint64_t state;
  std::uint64_t next() {
    auto z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
  }
  double fraction() { return (next() >> 11) * 0x1.0p-53; }
  std::uint64_t bounded(std::uint64_t n) {
    const auto threshold = -n % n;
    std::uint64_t x;
    do {
      x = next();
    } while (x < threshold);
    return x % n;
  }
};
std::size_t bytes(std::uint64_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(Real))
    throw std::runtime_error("Allocation size overflow");
  return static_cast<std::size_t>(count) * sizeof(Real);
}
// Scratch ownership only; no change to DataBox or Grid ownership semantics.
struct Buffer {
  Real *p;
  explicit Buffer(std::size_t n)
      : p(static_cast<Real *>(PORTABLE_MALLOC(bytes(n)))) {
    if (!p) throw std::bad_alloc();
  }
  ~Buffer() { PORTABLE_FREE(p); }
  Buffer(const Buffer &) = delete;
  Buffer &operator=(const Buffer &) = delete;
};
void hashBytes(std::uint64_t &h, const void *p, std::size_t n) {
  const auto *b = static_cast<const unsigned char *>(p);
  for (std::size_t i = 0; i < n; ++i)
    h = (h ^ b[i]) * UINT64_C(1099511628211);
}
struct Queries {
  std::vector<Real> x, y;
};
template <class Grid>
Queries queries(const Options &o, const Grid &gx, const Grid &gy, int dim) {
  Random rng{o.seed};
  Queries q{std::vector<Real>(o.queries),
            std::vector<Real>(dim == 2 ? o.queries : 0)};
  auto coordinate = [&](const Grid &g, int region) -> Real {
    const Real dx = (region == 1 || region == 3)
                        ? g.x(g.nPoints() - 1) - g.x(g.nPoints() - 2)
                        : g.x(1) - g.x(0);
    if (o.distribution == "boundary") {
      const Real f = Real(0.1 + 0.8 * rng.fraction());
      if (region == 0) return g.min() + f * dx;
      if (region == 1) return g.max() - f * dx;
      const Real distance = Real(0.1 + 0.9 * rng.fraction()) * dx;
      return region == 2 ? g.min() - distance : g.max() + distance;
    }
    const int cell = 1 + static_cast<int>(rng.bounded(g.nPoints() - 3));
    return Real(g.x(cell) +
                Real(0.1 + 0.8 * rng.fraction()) * (g.x(cell + 1) - g.x(cell)));
  };
  for (int i = 0; i < o.queries; ++i) {
    q.x[i] = coordinate(gx, i % 4);
    if (dim == 2) q.y[i] = coordinate(gy, (i / 4) % 4);
  }
  std::vector<int> order(o.queries);
  std::iota(order.begin(), order.end(), 0);
  if (o.distribution == "interior-ordered") {
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
      const int ja = dim == 2 ? gy.index(q.y[a]) : 0;
      const int jb = dim == 2 ? gy.index(q.y[b]) : 0;
      return ja != jb ? ja < jb : gx.index(q.x[a]) < gx.index(q.x[b]);
    });
  } else {
    for (int i = o.queries - 1; i > 0; --i)
      std::swap(order[i], order[rng.bounded(i + 1)]);
  }
  Queries result = q;
  for (int i = 0; i < o.queries; ++i) {
    result.x[i] = q.x[order[i]];
    if (dim == 2) result.y[i] = q.y[order[i]];
  }
  return result;
}

// Independent full-table host reference: endpoint slopes from Eqs. (2),(9),
// rational value differentiated via its quotient, rather than the stencil
// helper.
struct Ref {
  double value, deriv;
};
template <class Grid>
Ref reference1(double x, const Grid &g, const std::vector<double> &f,
               bool rational) {
  const int n = static_cast<int>(f.size());
  const int i = std::clamp(g.index(Real(x)), 0, n - 2);
  auto slope = [&](int k) { return (f[k + 1] - f[k]) / (g.x(k + 1) - g.x(k)); };
  const double d = g.x(i + 1) - g.x(i), q = x - g.x(i), r = d - q, s = slope(i);
  if (!rational || x < g.min() || x > g.max()) return {f[i] + q * s, s};
  if (i == 0 || i == n - 2) {
    double c = i == 0 ? (slope(1) - s) / (g.x(2) - g.x(0))
                      : (s - slope(i - 1)) / (g.x(i + 1) - g.x(i - 1));
    if (i == 0 && s * (s - d * c) <= 0) c = s / d;
    return {f[i] + q * (s - r * c), s + (q - r) * c};
  }
  double c1 = (s - slope(i - 1)) / (g.x(i + 1) - g.x(i - 1));
  const double c2 = (slope(i + 1) - s) / (g.x(i + 2) - g.x(i));
  if (i == 1 && slope(0) * (slope(0) - (g.x(1) - g.x(0)) * c1) <= 0)
    c1 = (s - 2 * slope(0)) / d;
  const double a = std::abs(c2), b = std::abs(c1), den = a * r + b * q;
  if (den == 0) return {f[i] + q * s, s};
  const double num = c1 * a * r + c2 * b * q, c = num / den;
  const double dc = ((c2 * b - c1 * a) * den - num * (b - a)) / (den * den);
  return {f[i] + q * (s - r * c), s + (q - r) * c - q * r * dc};
}
template <int Dim, class Grid>
std::array<double, 3> reference(Real x, Real y, const Grid &gx, const Grid &gy,
                                const std::vector<Real> &table, bool rational) {
  const int nx = gx.nPoints();
  if constexpr (Dim == 1) {
    const auto r = reference1(
        x, gx, std::vector<double>(table.begin(), table.end()), rational);
    return {r.value, r.deriv, 0};
  } else {
    const int i = std::clamp(gx.index(x), 0, nx - 2),
              j = std::clamp(gy.index(y), 0,
                             static_cast<int>(gy.nPoints()) - 2);
    auto row = [&](int k) {
      return reference1(x, gx,
                        std::vector<double>(table.begin() + k * nx,
                                            table.begin() + (k + 1) * nx),
                        rational);
    };
    auto col = [&](int k) {
      std::vector<double> f(gy.nPoints());
      for (int l = 0; l < gy.nPoints(); ++l)
        f[l] = table[l * nx + k];
      return reference1(y, gy, f, rational);
    };
    const auto a = row(j), b = row(j + 1), c = col(i), e = col(i + 1);
    const double dx = gx.x(i + 1) - gx.x(i), dy = gy.x(j + 1) - gy.x(j);
    const double u = (x - gx.x(i)) / dx, v = (y - gy.x(j)) / dy;
    const double f00 = table[j * nx + i], f10 = table[j * nx + i + 1],
                 f01 = table[(j + 1) * nx + i],
                 f11 = table[(j + 1) * nx + i + 1];
    return {
        (1 - v) * a.value + v * b.value + (1 - u) * c.value + u * e.value -
            ((1 - u) * (1 - v) * f00 + u * (1 - v) * f10 + (1 - u) * v * f01 +
             u * v * f11),
        (1 - v) * a.deriv + v * b.deriv +
            (e.value - c.value - (1 - v) * (f10 - f00) - v * (f11 - f01)) / dx,
        (1 - u) * c.deriv + u * e.deriv +
            (b.value - a.value - (1 - u) * (f01 - f00) - u * (f11 - f10)) / dy};
  }
}

template <int Dim, class Policy, bool Traffic = false, class Grid>
void launch(Box<Policy, Grid> db, const Real *x, const Real *y, Real *out,
            int n) {
  // Dispatch outside the lambdas: extended __host__ __device__ lambdas
  // cannot first-capture variables inside an if constexpr branch.
  if constexpr (Traffic && Dim == 1) {
    portableFor(
        "Rational benchmark traffic", 0, n, PORTABLE_LAMBDA(const int i) {
          out[i] = x[i];
          out[std::size_t(n) + i] = x[i];
        });
  } else if constexpr (Traffic) {
    portableFor(
        "Rational benchmark traffic", 0, n, PORTABLE_LAMBDA(const int i) {
          out[i] = x[i];
          out[std::size_t(n) + i] = x[i];
          out[2 * std::size_t(n) + i] = y[i];
        });
  } else if constexpr (Dim == 1) {
    portableFor(
        "Rational benchmark interpolation", 0, n,
        PORTABLE_LAMBDA(const int i) {
          const auto r = db.interpToRealWithGrads(x[i]);
          out[i] = r.value;
          out[std::size_t(n) + i] = r.d_dx1;
        });
  } else {
    portableFor(
        "Rational benchmark interpolation", 0, n,
        PORTABLE_LAMBDA(const int i) {
          const auto r = db.interpToRealWithGrads(y[i], x[i]);
          out[i] = r.value;
          out[std::size_t(n) + i] = r.d_dx1;
          out[2 * std::size_t(n) + i] = r.d_dx2;
        });
  }
}
std::array<double, 3> consume(const Real *out, int n, int components) {
  std::array<double, 3> sums{};
  for (int c = 0; c < components; ++c) {
    const Real *p = out + std::size_t(c) * n;
    portableReduce(
        "Benchmark checksum", 0, n,
        PORTABLE_LAMBDA(const int i, double &sum) {
          sum += static_cast<double>(p[i]);
        },
        sums[c]);
    if (!std::isfinite(sums[c])) throw std::runtime_error("Nonfinite checksum");
  }
  fence();
  return sums;
}
struct Trial {
  double seconds = 0;
  std::array<double, 3> sums{};
};
template <int Dim, class Policy, bool Traffic = false, class Grid>
Trial measure(const Options &o, Box<Policy, Grid> db, const Real *x,
              const Real *y, Real *out, int passes) {
  Trial t;
  for (int k = 0; k < passes; ++k) {
    fence();
    const auto start = Clock::now();
    launch<Dim, Policy, Traffic>(db, x, y, out, o.queries);
    fence();
    t.seconds += std::chrono::duration<double>(Clock::now() - start).count();
    const auto s = consume(out, o.queries, Dim + 1);
    for (int c = 0; c < Dim + 1; ++c)
      t.sums[c] += s[c];
  }
  return t;
}
double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double p = q * (v.size() - 1);
  const auto i = static_cast<std::size_t>(p);
  return v[i] + (p - i) * (v[std::min(i + 1, v.size() - 1)] - v[i]);
}

template <class Grid>
Grid makeGrid(int n, double extent) {
  if constexpr (std::is_same_v<Grid, RegularGrid>) {
    return Grid(0, extent, n);
  } else {
    // Smoothly varying intervals, with max/min spacing approximately three.
    std::vector<double> points(n);
    for (int i = 0; i < n; ++i) {
      const double t = double(i) / (n - 1);
      points[i] = extent * (t + t * t) / 2;
    }
    return Grid(points,
                typename Grid::Settings{.policy = Grid::Policy::RequireFast});
  }
}

template <int Dim, class Grid>
void run(const Options &o, std::ostream &csv, std::ostream &meta) {
  const int nx = Dim == 1 ? o.nx1d : o.nx2d, ny = Dim == 1 ? 1 : o.ny2d;
  const std::uint64_t count = std::uint64_t(nx) * ny;
  // DataBox's flattened indexing currently uses int.
  if (count > std::uint64_t(std::numeric_limits<int>::max()))
    throw std::runtime_error("Table exceeds DataBox flattened index limit");
  bytes(count);
  bytes(std::uint64_t(Dim + 1) * o.queries);
  Grid gx = makeGrid<Grid>(nx, 1), gy = makeGrid<Grid>(Dim == 1 ? 4 : ny, 2);
  const auto q = queries(o, gx, gy, Dim);
  std::vector<Real> table(count);
  const double pi = std::acos(-1.0);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const double x = gx.x(i), y = gy.x(j);
      if constexpr (Dim == 1)
        table[i] =
            2 + std::sin(6 * pi * x) + 0.2 * std::cos(14 * pi * x) + x * x;
      else
        table[std::size_t(j) * nx + i] =
            2 + std::sin(6 * pi * x) * std::cos(2 * pi * y) +
            0.2 * std::sin(4 * pi * x * y) + x * x - 0.3 * y * y;
    }
  std::uint64_t hash = UINT64_C(14695981039346656037);
  hashBytes(hash, table.data(), bytes(count));
  hashBytes(hash, q.x.data(), bytes(q.x.size()));
  hashBytes(hash, q.y.data(), bytes(q.y.size()));
  auto makeBox = [&](auto policy) {
    using P = decltype(policy);
    std::unique_ptr<Box<P, Grid>, Spiner::DBDeleter<true>> db;
    if constexpr (Dim == 1)
      db.reset(new Box<P, Grid>(Spiner::AllocationTarget::Device, nx));
    else
      db.reset(new Box<P, Grid>(Spiner::AllocationTarget::Device, ny, nx));
    db->setRange(0, gx.getOnDevice());
    if constexpr (Dim == 2) db->setRange(1, gy.getOnDevice());
    portableCopyToDevice(db->data(), table.data(), bytes(count));
    return db;
  };
  auto rational = makeBox(Rational{});
  auto linear = makeBox(Linear{});
  Buffer x(o.queries), y(Dim == 2 ? o.queries : 1),
      out(std::size_t(Dim + 1) * o.queries);
  portableCopyToDevice(x.p, q.x.data(), bytes(q.x.size()));
  if constexpr (Dim == 2)
    portableCopyToDevice(y.p, q.y.data(), bytes(q.y.size()));
  auto validate = [&](auto db, bool isRational) {
    using P = typename decltype(db)::GradStencilType;
    launch<Dim, P>(db, x.p, y.p, out.p, o.queries);
    fence();
    std::vector<Real> values(std::size_t(Dim + 1) * o.queries);
    portableCopyToHost(values.data(), out.p, bytes(values.size()));
    const int samples = std::min(o.queries, 257);
    const double tol = sizeof(Real) == sizeof(float) ? 2e-3 : 2e-9;
    for (int k = 0; k < samples; ++k) {
      const int i = samples == 1
                        ? 0
                        : static_cast<int>(std::uint64_t(k) * (o.queries - 1) /
                                           (samples - 1));
      const auto expected = reference<Dim>(q.x[i], Dim == 2 ? q.y[i] : Real(0),
                                           gx, gy, table, isRational);
      for (int c = 0; c < Dim + 1; ++c) {
        const double actual = values[std::size_t(c) * o.queries + i];
        if (!std::isfinite(actual) ||
            std::abs(actual - expected[c]) > tol * (1 + std::abs(expected[c])))
          throw std::runtime_error(
              "Reference validation failed: dimension=" + std::to_string(Dim) +
              " policy=" + (isRational ? "rational" : "linear") +
              " component=" + std::to_string(c) + " query=" +
              std::to_string(i) + " actual=" + std::to_string(actual) +
              " expected=" + std::to_string(expected[c]));
      }
    }
  };
  if (o.policy != "linear") validate(*rational, true);
  if (o.policy != "rational") validate(*linear, false);
  meta << "validated_dimension=" << Dim
       << " samples=" << std::min(o.queries, 257) << " input_hash=" << hash
       << '\n';
  std::array<std::vector<double>, 3> times;
  const std::array<std::string, 3> names = {"rational", "linear", "traffic"};
  auto enabled = [&](int p) {
    return p == 2 || o.policy == "both" || o.policy == names[p];
  };
  auto execute = [&](int p, int passes) {
    if (p == 0)
      return measure<Dim, Rational>(o, *rational, x.p, y.p, out.p, passes);
    if (p == 1)
      return measure<Dim, Linear>(o, *linear, x.p, y.p, out.p, passes);
    return measure<Dim, Linear, true>(o, *linear, x.p, y.p, out.p, passes);
  };
  const std::uint64_t calls = std::uint64_t(o.queries) * o.passes;
  auto row = [&](const std::string &kind, int p, int trial, double seconds,
                 std::array<double, 3> sums, double minimum, double iqr,
                 double ratio) {
    csv << kind << ',' << names[p] << ',' << Dim << ',' << o.grid << ',' << nx
        << ',' << ny << ',' << o.distribution << ',' << o.seed << ',' << hash
        << ',' << o.queries << ',' << o.passes << ',' << trial << ',' << calls
        << ',' << bytes(count) << ',' << bytes(std::uint64_t(Dim) * o.queries)
        << ',' << bytes(std::uint64_t(Dim + 1) * o.queries) << ',' << seconds
        << ',' << seconds * 1e9 / calls << ',' << calls / seconds / 1e6 << ','
        << sums[0] << ',' << sums[1] << ',' << sums[2] << ',' << minimum << ','
        << iqr << ',' << ratio << '\n';
  };
  for (int p = 0; p < 3; ++p)
    if (enabled(p)) execute(p, o.warmup);
  for (int trial = 0; trial < o.trials; ++trial) {
    // Alternate policy order; control runs last.
    for (int k = 0; k < 3; ++k) {
      const int p = k == 2 ? 2 : (trial % 2 ? 1 - k : k);
      if (!enabled(p)) continue;
      const auto t = execute(p, o.passes);
      if (t.seconds <= 0)
        throw std::runtime_error("Timer resolution insufficient");
      times[p].push_back(t.seconds);
      row("trial", p, trial, t.seconds, t.sums, 0, 0, 0);
      if (trial == 0 && t.seconds / o.passes < 0.001)
        std::cerr << "Warning: " << Dim << "D " << names[p]
                  << " passes <1 ms; increase --queries for stable timing.\n";
    }
  }
  const double ratio = o.policy == "both"
                           ? quantile(times[1], 0.5) / quantile(times[0], 0.5)
                           : 0;
  for (int p = 0; p < 3; ++p)
    if (enabled(p)) {
      const double med = quantile(times[p], 0.5), min = quantile(times[p], 0),
                   iqr = quantile(times[p], 0.75) - quantile(times[p], 0.25);
      row("summary", p, -1, med, {}, min * 1e9 / calls, iqr * 1e9 / calls,
          p == 0 ? ratio : 0);
      std::cout << Dim << "D " << names[p] << ": median " << med * 1e9 / calls
                << " ns/call; IQR " << iqr * 1e9 / calls << " ns/call\n";
    }
  fence();
  gx.finalize();
  gy.finalize();
}
} // namespace

int main(int argc, char **argv) {
  // Help also works without initializing an accelerator.
  for (int i = 1; i < argc; ++i)
    if (std::string(argv[i]) == "--help") {
      usage();
      return 0;
    }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  int status = 0;
  try {
    const auto o = parse(argc, argv);
    if (std::string(benchmarkBuildMetadata).find("build_type=Debug") !=
        std::string::npos)
      std::cerr
          << "Warning: Debug build; use Release for performance baselines.\n";
    std::ofstream csv(o.output), meta(o.output + ".metadata.txt");
    if (!csv || !meta) throw std::runtime_error("Cannot open output files");
    csv << std::setprecision(17);
    csv << "kind,policy,dimension,grid,nx,ny,distribution,seed,input_hash,"
           "queries,passes,trial,calls,table_bytes,query_bytes,output_bytes,"
           "seconds,ns_per_call,million_calls_per_second,checksum_value,"
           "checksum_dx,checksum_dy,min_ns_per_call,iqr_ns_per_call,rational_"
           "over_linear_throughput\n";
    const auto now = std::time(nullptr);
    meta
        << "Copyright 2026 Triad National Security, LLC. All rights reserved.\n"
           "Generative AI assisted the benchmark producing this metadata.\n"
        << "utc=" << std::put_time(std::gmtime(&now), "%Y-%m-%dT%H:%M:%SZ")
        << '\n';
    meta << benchmarkBuildMetadata << "\nhardware_label=" << o.hardware
         << "\nreal_bytes=" << sizeof(Real) << "\nwarmup=" << o.warmup
         << "\ntrials=" << o.trials << "\ncommand=";
    for (int i = 0; i < argc; ++i)
      meta << std::quoted(argv[i]) << ' ';
    meta << "\nchecksums=per-trial sum of per-pass reductions; summary "
            "checksum fields unused\n"
            "summary_seconds=median trial seconds; ratio unused fields are "
            "zero\n"
            "input_hash=FNV1a64 of native Real table then ordered x/y query "
            "bytes\n";
#ifdef PORTABILITY_STRATEGY_KOKKOS
    meta << "kokkos_version=" << KOKKOS_VERSION
         << "\nbackend=" << Kokkos::DefaultExecutionSpace::name()
         << "\nconcurrency=" << Kokkos::DefaultExecutionSpace().concurrency()
         << '\n';
    Kokkos::print_configuration(meta, true);
#else
    meta << "backend=Ports-of-Call native\n";
#endif
    for (const char *name : {"OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES",
                             "CUDA_VISIBLE_DEVICES"})
      if (const char *v = std::getenv(name)) meta << name << '=' << v << '\n';
    if (o.grid == "regular") {
      if (o.dim != "2") run<1, RegularGrid>(o, csv, meta);
      if (o.dim != "1") run<2, RegularGrid>(o, csv, meta);
    } else {
      if (o.dim != "2") run<1, FastGrid>(o, csv, meta);
      if (o.dim != "1") run<2, FastGrid>(o, csv, meta);
    }
    csv.flush();
    meta.flush();
    if (!csv || !meta) throw std::runtime_error("Failed writing results");
    std::cout << "Results: " << o.output << " (metadata: " << o.output
              << ".metadata.txt)\n";
  } catch (const std::exception &e) {
    std::cerr << "rational benchmark: " << e.what() << '\n';
    status = 1;
  }
#ifdef PORTABILITY_STRATEGY_KOKKOS
  Kokkos::finalize();
#endif
  return status;
}
