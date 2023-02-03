// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2022-2022 Ioannis Iakovidis

// Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of
// Code 2022 program.

// Licensed under GNU LGPL.3, see LICENCE file

// References
// Yunbum Kook, Yin Tat Lee, Ruoqi Shen, Santosh S. Vempala. "Sampling with
// Riemannian Hamiltonian
// Monte Carlo in a Constrained Space"
#ifndef AUTO_TUNER_HPP
#define AUTO_TUNER_HPP
#include "random_walks/crhmc/additional_units/dynamic_regularizer.hpp"
#include "random_walks/crhmc/additional_units/dynamic_step_size.hpp"
#include "random_walks/crhmc/additional_units/dynamic_weight.hpp"

// This class is responsible for calling the additional crhmc modules for:
// modifying the weights, ode step size and regularizer factor addaptively.

template <typename Sampler, typename RandomNumberGenerator>
class auto_tuner {
  using weight_tuner = dynamic_weight<Sampler, RandomNumberGenerator>;
  using regularizion_tuner =
      dynamic_regularizer<Sampler, RandomNumberGenerator>;
  using step_size_tuner = dynamic_step_size<Sampler>;

  using Opts = typename Sampler::Opts;

public:
  Opts options;
  std::unique_ptr<weight_tuner> tune_weights;
  std::unique_ptr<regularizion_tuner> tune_regularization;
  std::unique_ptr<step_size_tuner> tune_step_size;
  auto_tuner(Sampler &s) :
    options(s.params.options)
  {
    if (options.DynamicWeight) {
      tune_weights = std::unique_ptr<weight_tuner>(new weight_tuner(s));
    }
    if (options.DynamicRegularizer) {
      tune_regularization = std::unique_ptr<regularizion_tuner>(new regularizion_tuner(s));
    }
    if (options.DynamicStepSize) {
      tune_step_size = std::unique_ptr<step_size_tuner>(new step_size_tuner(s));
    }
  }
  void updateModules(Sampler &s, RandomNumberGenerator &rng) {
    if (options.DynamicWeight) {
      tune_weights->update_weights(s, rng);
    }
    if (options.DynamicRegularizer) {
      tune_regularization->update_regularization_factor(s, rng);
    }
    if (options.DynamicStepSize) {
      tune_step_size->update_step_size(s);
    }
  }
};
#endif
