#pragma once

#include "experiments.h"

// Results are written to json_root
namespace exp_impl {
auto do_wim_experiment(const WIMAdjacencyListPair& graph, const CommonExperimentParams& common,
                       const WIMSketchingParams& sketching, json* json_root) -> ResultVoid;

auto do_wbim_experiment(const WBIMAdjacencyListPair& graph, const VertexSet& seeds,
                        const CommonExperimentParams& common, const WBIMSketchingParams& sketching, json* json_root)
    -> ResultVoid;
} // namespace exp_impl
