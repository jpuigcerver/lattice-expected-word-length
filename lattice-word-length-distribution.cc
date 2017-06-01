// latbin/lattice-word-length-distribution.cc

// MIT License

// Copyright (c) 2017 Joan Puigcerver

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/fstext-lib.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace kaldi {

void AddInsPenToLattice(BaseFloat penalty, CompactLattice *lat) {
  typedef typename CompactLattice::Arc Arc;
  typedef typename CompactLattice::Weight::W Weight;
  for (int32 state = 0; state < lat->NumStates(); ++state) {
    for (fst::MutableArcIterator<CompactLattice> aiter(lat, state);
         !aiter.Done(); aiter.Next()) {
      Arc arc(aiter.Value());
      if (arc.olabel != 0) {
        Weight weight = arc.weight.Weight();
        weight.SetValue1(weight.Value1() + penalty);
        arc.weight.SetWeight(weight);
        aiter.SetValue(arc);
      }
    }
  }
}

}  // namespace kaldi

namespace fst {

template <typename Float1, typename Float2, typename Int>
inline void ConvertLatticeWeight(
    const CompactLatticeWeightTpl<LatticeWeightTpl<Float1>, Int> &w_in,
    LogWeightTpl<Float2> *w_out){
  LogWeightTpl<Float2> w1(w_in.Weight().Value1());
  LogWeightTpl<Float2> w2(w_in.Weight().Value2());
  *w_out = Times(w1, w2);
}

}  // namespace fst

int main(int argc, char *argv[]) {
  try {
    using namespace kaldi;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::LogArc;
    using fst::StdVectorFst;
    typedef VectorFst<LogArc> LogVectorFst;

    const char *usage =
        "Prints the distribution of the length of the transcriptions "
        "in a lattice. This typically means the distribution over the number "
        "of words in the transcription.\n"
        "\n"
        "Usage: lattice-word-length-distribution [options] lattice-rspecifier1 "
        "[lattice-rspecifier2 ...]\n"
        " e.g.: lattice-word-length-distribution ark:1.lats ark:2.lats\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat insertion_penalty = 0.0;
    int32 nbest = std::numeric_limits<int32>::max();

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("insertion-penalty", &insertion_penalty,
                "Add this penalty to the lattice arcs with non-epsilon output "
                "label (typically, equivalent to word insertion penalty).");
    po.Register("nbest", &nbest,
                "Limit the distribution to this number of n-best lengths.");

    po.Read(argc, argv);

    if (po.NumArgs() < 1) {
      po.PrintUsage();
      exit(1);
    }

    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;

    for (int32 a = 1; a <= po.NumArgs(); ++a) {
      for (SequentialCompactLatticeReader lattice_reader(po.GetArg(a));
           !lattice_reader.Done(); lattice_reader.Next()) {
        const std::string key = lattice_reader.Key();
        LogVectorFst fst;
        {
          CompactLattice clat = lattice_reader.Value();
          lattice_reader.FreeCurrent();
          // TopSort compact lattice
          TopSortCompactLatticeIfNeeded(&clat);
          // Acoustic scale
          if (acoustic_scale != 1.0 || graph_scale != 1.0)
            fst::ScaleLattice(scale, &clat);
          // Word insertion penalty
          if (insertion_penalty != 0.0)
            AddInsPenToLattice(insertion_penalty, &clat);
          // Remove the alignments from the lattice.
          RemoveAlignmentsFromCompactLattice(&clat);
          // Convert CompactLattice to Fst in the log semiring
          ConvertLattice(clat, &fst);
        }

        // Remove epsilons from the Fst, since it is acyclic: O(V^2 + VE).
        RmEpsilon(&fst);

        // Replace all the labels with 1, O(V + E).
        for (fst::StateIterator<LogVectorFst> siter(fst); !siter.Done();
             siter.Next()) {
          for (fst::MutableArcIterator<LogVectorFst> aiter(&fst, siter.Value());
               !aiter.Done(); aiter.Next()) {
            LogArc arc = aiter.Value();
            arc.ilabel = arc.olabel = 1;
            aiter.SetValue(arc);
          }
        }

        // Compute the backward cost of each state, so that bw[fst.Start()] has
        // the total cost of the fst (i.e. likelihood = -cost), O(V + E).
        std::vector<LogArc::Weight> bw;
        fst::ShortestDistance<LogArc>(fst, &bw, true);
        const double total_cost = bw[fst.Start()].Value();

        // 1. Determinize the fst in the Log semiring, so that each word length
        // is represented by a single path whose likelihood is the sum of
        // all likelihoods with the same word length.
        // 2. Convert from LogArc to StdArc (tropical semiring).
        // 3. Find the n-best paths in the tropical semiring.
        typedef fst::WeightConvertMapper<fst::LogArc, fst::StdArc> WeightMapper;
        StdVectorFst nbest_fst;
        fst::ShortestPath(
            fst::ArcMapFst<fst::LogArc, fst::StdArc, WeightMapper>(
                fst::DeterminizeFst<fst::LogArc>(fst),
                WeightMapper()), &nbest_fst, nbest);

        // Get each path as separate fst.
        std::vector<StdVectorFst> nbest_fst_vector;
        fst::ConvertNbestToVector(nbest_fst, &nbest_fst_vector);

        // Print the length of the path and it's log-probability.
        std::vector<std::pair<size_t, double>> nbest_pairs;
        std::cout << key;
        for (const StdVectorFst& path_fst : nbest_fst_vector) {
          std::vector<int32> path_isymbs;
          std::vector<int32> path_osymbs;
          StdArc::Weight     path_cost;
          fst::GetLinearSymbolSequence<StdArc, int32>(
              path_fst, &path_isymbs, &path_osymbs, &path_cost);
          // logprobability = path_loglikelihood - total_loglikelihood
          //                = total_cost - path_cost
          std::cout << " "
                    << path_isymbs.size() << " "
                    << total_cost - path_cost.Value()
                    << " ;";
        }
        std::cout << std::endl;
      }
    }

    return 0;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return -1;
  }
}
